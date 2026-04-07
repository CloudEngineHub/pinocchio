#
# Copyright (c) 2024-2025 INRIA
#

import unittest

import numpy as np
import pinocchio as pin
from test_case import PinocchioTestCase as TestCase


def build_four_bar_model():
    """Build the four-bar linkage model from the simulation-lcaba example.

    Returns (model, constraint_model, q_sol) where q_sol satisfies the
    closed-loop constraint.
    """
    import coal

    height = 0.1
    width = 0.01
    radius = 0.05

    mass_link_A = 10.0
    length_link_A = 1.0

    mass_link_B = 5.0
    length_link_B = 0.6

    inertia_link_A = pin.Inertia.FromBox(mass_link_A, length_link_A, width, height)
    placement_center_link_A = pin.SE3.Identity()
    placement_center_link_A.translation = pin.XAxis * length_link_A / 2.0

    inertia_link_B = pin.Inertia.FromBox(mass_link_B, length_link_B, width, height)
    placement_center_link_B = pin.SE3.Identity()
    placement_center_link_B.translation = pin.XAxis * length_link_B / 2.0

    model = pin.Model()
    base_joint_id = 0

    joint1_placement = pin.SE3.Identity()
    joint1_placement.translation = pin.XAxis * length_link_A / 2.0
    joint1_id = model.addJoint(
        base_joint_id, pin.JointModelRY(), joint1_placement, "link_B1"
    )
    model.appendBodyToJoint(joint1_id, inertia_link_B, placement_center_link_B)

    joint2_placement = pin.SE3.Identity()
    joint2_placement.translation = pin.XAxis * length_link_B
    joint2_id = model.addJoint(
        joint1_id, pin.JointModelRY(), joint2_placement, "link_A2"
    )
    model.appendBodyToJoint(joint2_id, inertia_link_A, placement_center_link_A)

    joint3_placement = pin.SE3.Identity()
    joint3_placement.translation = pin.XAxis * length_link_A
    joint3_id = model.addJoint(
        joint2_id, pin.JointModelRY(), joint3_placement, "link_B2"
    )
    model.appendBodyToJoint(joint3_id, inertia_link_B, placement_center_link_B)

    # Closed-loop constraint: tip of link_B2 connects back to base
    constraint1_joint1_placement = pin.SE3.Identity()
    constraint1_joint1_placement.translation = pin.XAxis * length_link_B

    constraint1_joint2_placement = pin.SE3.Identity()
    constraint1_joint2_placement.translation = -pin.XAxis * length_link_A / 2.0

    constraint_model = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,
        model,
        joint3_id,
        constraint1_joint1_placement,
        base_joint_id,
        constraint1_joint2_placement,
    )

    # ── Inverse geometry: satisfy constraints at the initial configuration ──
    q = pin.neutral(model)
    constraint_data = constraint_model.createData()
    constraint_size = constraint_model.residualSize()
    constraint_models = [constraint_model]
    constraint_datas = [constraint_data]

    rho = 1e-10
    mu = 1e-4
    y = np.ones(constraint_size)
    data = model.createData()
    data.M = np.eye(model.nv) * rho
    kkt_constraint = pin.ContactCholeskyDecomposition(
        model, data, constraint_models, constraint_datas
    )

    for _ in range(100):
        pin.computeJointJacobians(model, data, q)
        data.q_in = q
        constraint_model.calc(model, data, constraint_data)
        kkt_constraint.compute(model, data, constraint_models, constraint_datas, mu)
        constraint_value = constraint_data.c1Mc2.translation
        J = pin.getFrameJacobian(
            model,
            data,
            constraint_model.joint1_id,
            constraint_model.joint1_placement,
            constraint_model.reference_frame,
        )[:3, :]
        if (
            np.linalg.norm(constraint_value, np.inf) < 1e-10
            and np.linalg.norm(J.T.dot(constraint_value + y), np.inf) < 1e-10
        ):
            break
        rhs = np.concatenate([-constraint_value - y * mu, np.zeros(model.nv)])
        dz = kkt_constraint.solve(rhs)
        dy = dz[:constraint_size]
        dq = dz[constraint_size:]
        alpha = 1.0
        q = pin.integrate(model, q, -alpha * dq)
        y -= alpha * (-dy + y)

    q_sol = (q[:] + np.pi) % np.pi - np.pi
    return model, constraint_model, q_sol


class TestLCABABindings(TestCase):
    def setUp(self):
        try:
            import coal  # noqa: F401
        except ImportError:
            self.skipTest("coal is not available")

        self.model, self.constraint_model, self.q_sol = build_four_bar_model()
        self.data = self.model.createData()

        self.constraint_data = self.constraint_model.createData()
        self.constraint_models = [self.constraint_model]
        self.constraint_datas = [self.constraint_data]

        self.v = np.zeros(self.model.nv)
        self.tau = np.zeros(self.model.nv)

        self.mu_sim = 1e-6
        self.prox_settings = pin.ProximalSettings(1e-12, self.mu_sim, 100)

    def _make_prox_settings(self):
        """Return a fresh ProximalSettings (settings are mutated by lcaba)."""
        return pin.ProximalSettings(
            self.prox_settings.absolute_accuracy,
            self.prox_settings.mu,
            self.prox_settings.max_iter,
        )

    # ── computeJointMinimalOrdering ────────────────────────────────────────

    def test_computeJointMinimalOrdering_runs(self):
        """computeJointMinimalOrdering must run without error."""
        data = self.model.createData()
        pin.computeJointMinimalOrdering(
            self.model, data, self.constraint_models
        )

    def test_computeJointMinimalOrdering_empty_constraints(self):
        """computeJointMinimalOrdering with no constraints must succeed."""
        data = self.model.createData()
        pin.computeJointMinimalOrdering(self.model, data, [])

    # ── lcaba ──────────────────────────────────────────────────────────────

    def test_lcaba_returns_correct_size(self):
        """lcaba must return an acceleration vector of size nv."""
        data = self.model.createData()
        pin.computeJointMinimalOrdering(
            self.model, data, self.constraint_models
        )
        ddq = pin.lcaba(
            self.model,
            data,
            self.q_sol,
            self.v,
            self.tau,
            self.constraint_models,
            self.constraint_datas,
            self._make_prox_settings(),
        )
        self.assertEqual(ddq.shape, (self.model.nv,))

    def test_lcaba_result_stored_in_data(self):
        """lcaba must store the result in data.ddq."""
        data = self.model.createData()
        pin.computeJointMinimalOrdering(
            self.model, data, self.constraint_models
        )
        ddq = pin.lcaba(
            self.model,
            data,
            self.q_sol,
            self.v,
            self.tau,
            self.constraint_models,
            self.constraint_datas,
            self._make_prox_settings(),
        )
        self.assertApprox(ddq, data.ddq)

    def test_lcaba_consistent_with_constraintDynamics(self):
        """lcaba and constraintDynamics must agree on ddq for a 3D constraint."""
        # Reference: constraintDynamics
        data_ref = self.model.createData()
        constraint_datas_ref = [self.constraint_model.createData()]
        prox_ref = pin.ProximalSettings(1e-14, 1e-5, 100)
        pin.initConstraintDynamics(
            self.model, data_ref, self.constraint_models, constraint_datas_ref
        )
        pin.constraintDynamics(
            self.model,
            data_ref,
            self.q_sol,
            self.v,
            self.tau,
            self.constraint_models,
            constraint_datas_ref,
            prox_ref,
        )

        # Tested: lcaba
        data_lcaba = self.model.createData()
        constraint_datas_lcaba = [self.constraint_model.createData()]
        prox_lcaba = pin.ProximalSettings(1e-14, 1e-5, 100)
        pin.computeJointMinimalOrdering(
            self.model, data_lcaba, self.constraint_models
        )
        pin.lcaba(
            self.model,
            data_lcaba,
            self.q_sol,
            self.v,
            self.tau,
            self.constraint_models,
            constraint_datas_lcaba,
            prox_lcaba,
        )

        self.assertApprox(data_ref.ddq, data_lcaba.ddq, eps=1e-6)

    def test_lcaba_no_constraints_matches_aba(self):
        """Without constraints, lcaba must match standard ABA."""
        # ABA reference
        data_aba = self.model.createData()
        ddq_aba = pin.aba(
            self.model, data_aba, self.q_sol, self.v, self.tau
        )

        # lcaba with empty constraint list
        data_lcaba = self.model.createData()
        pin.computeJointMinimalOrdering(self.model, data_lcaba, [])
        ddq_lcaba = pin.lcaba(
            self.model,
            data_lcaba,
            self.q_sol,
            self.v,
            self.tau,
            [],
            [],
            self._make_prox_settings(),
        )

        self.assertApprox(ddq_aba, ddq_lcaba, eps=1e-10)

    # ── Baumgarte stabilisation ────────────────────────────────────────────

    def test_baumgarte_parameters_accessible(self):
        """m_baumgarte_parameters must be accessible and settable."""
        cm = self.constraint_model
        self.assertTrue(hasattr(cm, "m_baumgarte_parameters"))
        cm.m_baumgarte_parameters.Kp = 10.0
        cm.m_baumgarte_parameters.Kd = 2.0 * np.sqrt(cm.m_baumgarte_parameters.Kp)
        self.assertAlmostEqual(cm.m_baumgarte_parameters.Kp, 10.0)
        self.assertAlmostEqual(
            cm.m_baumgarte_parameters.Kd, 2.0 * np.sqrt(10.0)
        )

    # ── Multi-step simulation ──────────────────────────────────────────────

    def test_simulation_steps_run(self):
        """A short simulation loop using lcaba must run without error."""
        model = self.model
        constraint_models = self.constraint_models
        constraint_datas = [self.constraint_model.createData()]

        # Baumgarte stabilisation
        constraint_models[0].m_baumgarte_parameters.Kp = 10.0
        constraint_models[0].m_baumgarte_parameters.Kd = 2.0 * np.sqrt(10.0)

        data_sim = model.createData()
        pin.computeJointMinimalOrdering(model, data_sim, constraint_models)

        q = self.q_sol.copy()
        v = np.zeros(model.nv)
        tau = np.zeros(model.nv)
        dt = 5e-3

        for _ in range(5):
            prox_iter = pin.ProximalSettings(1e-12, self.mu_sim, 100)
            a = pin.lcaba(
                model,
                data_sim,
                q,
                v,
                tau,
                constraint_models,
                constraint_datas,
                prox_iter,
            )
            self.assertEqual(a.shape, (model.nv,))
            v += a * dt
            q = pin.integrate(model, q, v * dt)

    def test_simulation_energy_bounded(self):
        """Under zero torque and Baumgarte stabilisation, kinetic energy
        must stay bounded (no unbounded growth) over a few steps."""
        model = self.model
        constraint_models = [self.constraint_model]
        constraint_datas = [self.constraint_model.createData()]

        constraint_models[0].m_baumgarte_parameters.Kp = 10.0
        constraint_models[0].m_baumgarte_parameters.Kd = 2.0 * np.sqrt(10.0)

        data_sim = model.createData()
        pin.computeJointMinimalOrdering(model, data_sim, constraint_models)

        q = self.q_sol.copy()
        v = np.zeros(model.nv)
        tau = np.zeros(model.nv)
        dt = 5e-3
        n_steps = 20

        energies = []
        for _ in range(n_steps):
            pin.crba(model, data_sim, q)
            ke = 0.5 * v @ data_sim.M @ v
            energies.append(ke)

            prox_iter = pin.ProximalSettings(1e-12, self.mu_sim, 100)
            a = pin.lcaba(
                model,
                data_sim,
                q,
                v,
                tau,
                constraint_models,
                constraint_datas,
                prox_iter,
            )
            v += a * dt
            q = pin.integrate(model, q, v * dt)

        # Starting from rest – energy should remain small (< 10 J over 20 steps)
        self.assertLess(max(energies), 10.0)


if __name__ == "__main__":
    unittest.main()
