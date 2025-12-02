//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/constraints/utils.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/algorithm/ipm-solver.hpp"
#include "pinocchio/algorithm/contact-jacobian.hpp"

#include "pinocchio/algorithm/admm-solver.hpp"
#include "pinocchio/algorithm/clarabel-solver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

double mu = 1e-4;

void buildStackOfCubesModel(
  std::vector<double> masses,
  ::pinocchio::Model & model,
  std::vector<PointContactConstraintModel> & constraint_models)
{
  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  const int n_cubes = (int)masses.size();

  for (int i = 0; i < n_cubes; i++)
  {
    const double box_mass = masses[(std::size_t)i];
    const Inertia box_inertia = Inertia::FromBox(box_mass, box_dims[0], box_dims[1], box_dims[2]);
    JointIndex joint_id =
      model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "free_flyer_" + std::to_string(i));
    model.appendBodyToJoint(joint_id, box_inertia);
  }

  const double friction_value = 0.4;
  for (int i = 0; i < n_cubes; i++)
  {
    const SE3 local_placement_box_1(
      SE3::Matrix3::Identity(), 0.5 * SE3::Vector3(box_dims[0], box_dims[1], box_dims[2]));
    const SE3 local_placement_box_2(
      SE3::Matrix3::Identity(), 0.5 * SE3::Vector3(box_dims[0], box_dims[1], -box_dims[2]));
    SE3::Matrix3 rot = SE3::Matrix3::Identity();
    for (int j = 0; j < 4; ++j)
    {
      const SE3 local_placement_1(
        SE3::Matrix3::Identity(), rot * local_placement_box_1.translation());
      const SE3 local_placement_2(
        SE3::Matrix3::Identity(), rot * local_placement_box_2.translation());
      PointContactConstraintModel cm(
        model, (JointIndex)i, local_placement_1, (JointIndex)i + 1, local_placement_2);
      cm.setFriction(friction_value);
      constraint_models.push_back(cm);
      rot = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix() * rot;
    }
  }
}

Eigen::Vector3d computeFtotOfFirstBoxInStackOfBoxes(const Eigen::VectorXd & contact_forces)
{
  // we make the assumption that each box has contact constraints at each corner
  PINOCCHIO_THROW_IF(
    contact_forces.size() < 3 * 4, std::logic_error, "Invalid number of contact forces");

  Eigen::Vector3d f_tot = Eigen::Vector3d::Zero();
  for (int k = 0; k < 4; ++k)
  {
    f_tot += contact_forces.segment(3 * k, 3);
  }
  return f_tot;
}

template<typename _ConstraintModel>
struct TestBoxTpl
{
  typedef _ConstraintModel ConstraintModel;

  typedef typename ConstraintModel::ConstraintData ConstraintData;
  typedef typename ConstraintModel::ConstraintSet ConstraintSet;

  TestBoxTpl(const Model & model, const std::vector<ConstraintModel> & constraint_models)
  : model(model)
  , data(model)
  , constraint_models(constraint_models)
  , v_next(Eigen::VectorXd::Zero(model.nv))
  {
    for (const auto & cm : constraint_models)
    {
      constraint_datas.push_back(cm.createData());
    }

    const Eigen::DenseIndex constraint_size = getTotalConstraintMaxResidualSize(constraint_models);
    primal_solution = dual_solution = dual_solution_sparse = Eigen::VectorXd::Zero(constraint_size);
  }

  void operator()(
    const Eigen::VectorXd & q0,
    const Eigen::VectorXd & v0,
    const Eigen::VectorXd & tau0,
    const Force & fext,
    const double dt)
  {
    std::vector<Force> external_forces(size_t(model.njoints), Force::Zero());
    external_forces[1] = fext;

    const Eigen::VectorXd v_free =
      v0 + dt * aba(model, data, q0, v0, tau0, external_forces, Convention::WORLD);
    data.q_in = q0;
    data.v_in = v0;
    data.tau_in = tau0;
    calc(model, data, constraint_models, constraint_datas);

    // Cholesky of the Delassus matrix
    crba(model, data, q0, Convention::WORLD);
    ContactCholeskyDecomposition chol(model, data, constraint_models, constraint_datas);
    chol.resize(model, constraint_models, constraint_datas);
    chol.compute(model, data, constraint_models, constraint_datas, 1e-10);
    const Eigen::MatrixXd delassus_matrix_plain = chol.getDelassusCholeskyExpression().matrix();
    auto delassus_dense = chol.getDelassusCholeskyExpression().dense();
    Eigen::VectorXd compliance =
      Eigen::VectorXd::Zero(int(residualSize(constraint_models, constraint_datas)));
    delassus_dense.updateCompliance(compliance);
    delassus_dense.updateDamping(0.);

    const double tol_abs = 1e-12;
    const double tol_rel = 1e-12;
    const bool solve_ncp = false;

    // Compute g = constraints drift (constraints free velocity)
    Eigen::MatrixXd constraint_jacobian(delassus_matrix_plain.rows(), model.nv);
    constraint_jacobian.setZero();
    getConstraintsJacobian(model, data, constraint_models, constraint_datas, constraint_jacobian);
    const Eigen::VectorXd g = constraint_jacobian * v_free;

    // Configure the IP solver
    IPMConstraintSolverTpl<double> ip_solver(
      int(residualSize(constraint_models, constraint_datas)));
    ip_solver.setMaxIterations(10000);
    ip_solver.setAbsolutePrecision(tol_abs);
    ip_solver.setRelativePrecision(tol_rel);
    has_converged = ip_solver.solve(
      delassus_dense, g, constraint_models, constraint_datas, boost::none, boost::none, solve_ncp);
    primal_solution = ip_solver.getPrimalSolution();
    dual_solution = ip_solver.getDualSolution();
    n_iter = ip_solver.getIterationCount();
    const Eigen::VectorXd tau_ext = constraint_jacobian.transpose() * primal_solution / dt;

    // TODO: test warmstart

    v_next =
      v0
      + dt * aba(model, data, q0, v0, (tau0 + tau_ext).eval(), external_forces, Convention::WORLD);

    // Run ADMM
    ADMMConstraintSolverTpl<double> admm_solver(
      int(residualSize(constraint_models, constraint_datas)));
    admm_solver.setMaxIterations(10000);
    admm_solver.setAbsolutePrecision(tol_abs);
    admm_solver.setRelativePrecision(tol_rel);
    admm_solver.solve(
      delassus_dense, g, constraint_models, constraint_datas, boost::none, boost::none, boost::none,
      solve_ncp, ADMMUpdateRule::OSQP);

#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT
    // Run CLARABEL
    ClarabelContactSolverTpl<double> clarabel_solver(
      int(residualSize(constraint_models, constraint_datas)));
    clarabel_solver.setMaxIterations(10000);
    clarabel_solver.setAbsolutePrecision(tol_abs);
    clarabel_solver.setRelativePrecision(tol_rel);
    auto G_expression = chol.getDelassusCholeskyExpression();
    clarabel_solver.solve(
      G_expression, g, constraint_models, constraint_datas, boost::none, boost::none, boost::none,
      solve_ncp, false, false /*verbose*/);

    std::cout << "NUMIT CLARABEL: " << clarabel_solver.getIterationCount()
              << " / NUMIT IPSolver: " << ip_solver.getIterationCount()
              << " / NUMIT ADMMSolver: " << admm_solver.getIterationCount()
              << " (chol updates: " << admm_solver.getDelassusDecompositionUpdateCount() << ")\n";
#else
    std::cout << " / NUMIT IPSolver: " << ip_solver.getIterationCount()
              << " / NUMIT ADMMSolver: " << admm_solver.getIterationCount()
              << " (chol updates: " << admm_solver.getDelassusDecompositionUpdateCount() << ")\n";
#endif

    std::cout << "TIMINGS IPSolver: " << ip_solver.getCPUTimes().user << "\n";
    std::cout << "TIMINGS ADMMSolver: " << admm_solver.getCPUTimes().user << "\n";
  }

  Model model;
  Data data;
  std::vector<ConstraintModel> constraint_models;
  std::vector<ConstraintData> constraint_datas;
  Eigen::VectorXd v_next;

  Eigen::VectorXd primal_solution, dual_solution, dual_solution_sparse;
  bool has_converged;
  int n_iter;
};

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(ball)
{
  Model model;
  model.addJoint(0, JointModelFreeFlyer(), SE3::Identity(), "free_flyer");

  const double ball_dim = 1;
  const double ball_mass = 10;
  const Inertia ball_inertia = Inertia::FromSphere(ball_mass, ball_dim);

  model.appendBodyToJoint(1, ball_inertia);

  BOOST_CHECK(model.check(model.createData()));

  Eigen::VectorXd q0 = neutral(model);
  q0.const_cast_derived()[2] += ball_dim / 2;
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  typedef PointContactConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;

  const double friction_value = 0.4;
  {
    const SE3 local_placement_ball(SE3::Matrix3::Identity(), SE3::Vector3(0, 0, -ball_dim));
    ConstraintModel cm(model, 0, SE3::Identity(), 1, local_placement_ball);
    cm.setFriction(friction_value);
    constraint_models.push_back(cm);
  }

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    const Force::Vector3 f_tot_ref = -ball_mass * Model::gravity981 - fext.linear();
    Force::Vector3 f_tot = test.primal_solution.head(3) / dt;
    BOOST_CHECK(f_tot.isApprox(f_tot_ref, 1e-8));
    BOOST_CHECK(test.v_next.isZero(2e-10));
  }
}

BOOST_AUTO_TEST_CASE(box)
{
  Model model;
  typedef PointContactConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;
  const double box_mass = 1e1;
  const std::vector<double> masses = {box_mass};

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();
  buildStackOfCubesModel(masses, model, constraint_models);

  const int num_tests =
#ifdef NDEBUG
    1
#else
    1
#endif
    ;

  BOOST_CHECK(model.check(model.createData()));

  Eigen::VectorXd q0 = neutral(model);
  q0[2] += box_dims[2] / 2;
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    const Force::Vector3 f_tot_ref = -box_mass * Model::gravity981 - fext.linear();
    const Force::Vector3 f_tot = computeFtotOfFirstBoxInStackOfBoxes(test.primal_solution / dt);
    BOOST_CHECK(f_tot.isApprox(f_tot_ref, 1e-8));
    BOOST_CHECK(test.v_next.isZero(2e-10));
  }

  const double friction_value = 0.4;
  const double f_sliding = friction_value * Model::gravity981.norm() * box_mass;

  // Test static motion with small external force
  for (int k = 0; k < num_tests; ++k)
  {
    const double scaling = 0.9;
    Force fext = Force::Zero();
    fext.linear().head<2>().setRandom().normalize();
    fext.linear() *= scaling * f_sliding;

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(1e-8));
    const Force::Vector3 f_tot_ref = -box_mass * Model::gravity981 - fext.linear();
    const Force::Vector3 f_tot = computeFtotOfFirstBoxInStackOfBoxes(test.primal_solution / dt);
    BOOST_CHECK(f_tot.isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(test.v_next.isZero(1e-8));
  }

  // Test sliding motion
  for (int k = 0; k < num_tests; ++k)
  {
    const double scaling = 1.1;
    Force fext = Force::Zero();
    fext.linear().head<2>().setRandom().normalize();
    fext.linear() *= scaling * f_sliding;

    std::cout << "-- SLIDING MOTION --\n";
    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    const Force::Vector3 f_tot_ref = -box_mass * Model::gravity981 - 1 / scaling * fext.linear();
    const Force::Vector3 f_tot = computeFtotOfFirstBoxInStackOfBoxes(test.primal_solution / dt);
    std::cout << "-- STATIC SMALL FORCE --\n";
    std::cout << "f_tot_ref = " << f_tot_ref.transpose() << "\n";
    std::cout << "f_tot = " << f_tot.transpose() << "\n";
    std::cout << "v_next = " << test.v_next.transpose() << "\n";
    BOOST_CHECK(f_tot.isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(
      math::fabs(Motion(test.v_next).linear().norm() - (f_sliding * 0.1 / box_mass * dt)) <= 1e-6);
    BOOST_CHECK(Motion(test.v_next).angular().isZero(1e-6));
  }
}

BOOST_AUTO_TEST_CASE(stack_of_boxes)
{
  const int n_cubes = 10;
  // const double conditionning = 1e6;
  const double conditionning = 1e3;
  const double mass_factor = std::pow(conditionning, 1. / (n_cubes - 1));
  std::vector<double> masses;
  double mass_tot = 0;
  for (int i = 0; i < n_cubes; i++)
  {
    const double box_mass = 1e-3 * std::pow(mass_factor, i);
    masses.push_back(box_mass);
    mass_tot += box_mass;
  }

  Model model;
  typedef PointContactConstraintModel ConstraintModel;
  typedef TestBoxTpl<ConstraintModel> TestBox;
  std::vector<ConstraintModel> constraint_models;

  buildStackOfCubesModel(masses, model, constraint_models);
  BOOST_CHECK(model.check(model.createData()));

  const SE3::Vector3 box_dims = SE3::Vector3::Ones();

  Eigen::VectorXd q0 = neutral(model);
  for (int i = 0; i < n_cubes; i++)
  {
    q0[7 * i + 2] = i * box_dims[2];
    q0[7 * i + 2] += box_dims[2] / 2;
  }
  const Eigen::VectorXd v0 = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd tau0 = Eigen::VectorXd::Zero(model.nv);

  const double dt = 1e-3;

  // Test static motion with zero external force
  {
    const Force fext = Force::Zero();

    TestBox test(model, constraint_models);
    test(q0, v0, tau0, fext, dt);

    BOOST_CHECK(test.has_converged == true);
    BOOST_CHECK(test.dual_solution.isZero(2e-10));
    // We check the total force applied on the bottom box of the stack
    const Force::Vector3 f_tot_ref = -mass_tot * Model::gravity981;
    const Force::Vector3 f_tot = computeFtotOfFirstBoxInStackOfBoxes(test.primal_solution / dt);
    std::cout << "f_tot_ref = " << f_tot_ref.transpose() << "\n";
    std::cout << "f_tot     = " << f_tot.transpose() << "\n";
    BOOST_CHECK(f_tot.isApprox(f_tot_ref, 1e-6));
    BOOST_CHECK(test.v_next.isZero(1e-8));
  }
}

BOOST_AUTO_TEST_SUITE_END()
