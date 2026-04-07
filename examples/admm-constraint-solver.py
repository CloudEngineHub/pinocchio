"""
Stack of two cubes solved with the ADMM constraint solver.

This example demonstrates how to:
  - Build a kinematic model with two free-floating cubes
  - Build a geometry model with the cube shapes and a ground plane
  - Manually define 8 PointContactConstraintModel constraints:
      4 for the floor-cube 1 interaction (bottom corners of cube 1)
      4 for the cube 1-cube 2 interaction (corners at the shared face)
  - Build the Delassus operator via the Cholesky decomposition
  - Compute the constraint drift vector  g = J * v_free
  - Solve the constrained problem with the ADMM solver
  - Print convergence statistics
"""

import sys

import numpy as np
import pinocchio as pin

# ─── 1. Kinematic / dynamic model ────────────────────────────────────────────

model = pin.Model()
# Default gravity is already (0, 0, -9.81, 0, 0, 0); shown here for clarity.
model.gravity = pin.Motion(np.array([0.0, 0.0, -9.81, 0.0, 0.0, 0.0]))

box_size = 1.0  # edge length of each cube [m]
box_half = box_size / 2.0
box_mass = 1.0  # mass of each cube [kg]

box_inertia = pin.Inertia.FromBox(box_mass, box_size, box_size, box_size)

# Cube 1 - free-flyer joint branching from the universe (joint 0)
joint1_id = model.addJoint(
    0, pin.JointModelFreeFlyer(), pin.SE3.Identity(), "box1_joint"
)
model.appendBodyToJoint(joint1_id, box_inertia, pin.SE3.Identity())

# Cube 2 - free-flyer joint, also branching from the universe
joint2_id = model.addJoint(
    0, pin.JointModelFreeFlyer(), pin.SE3.Identity(), "box2_joint"
)
model.appendBodyToJoint(joint2_id, box_inertia, pin.SE3.Identity())

# ─── 2. Geometry model ────────────────────────────────────────────────────────

geom_model = pin.GeometryModel()

try:
    import coal

    cube_shape = coal.Box(box_size, box_size, box_size)

    # Cube 1 - geometry in joint1's local frame (centred at the joint origin)
    geom_box1 = pin.GeometryObject("box1", joint1_id, 0, pin.SE3.Identity(), cube_shape)
    geom_box1.meshColor = np.array([0.2, 0.6, 0.2, 1.0])
    geom_model.addGeometryObject(geom_box1)

    # Cube 2 - geometry in joint2's local frame
    geom_box2 = pin.GeometryObject("box2", joint2_id, 0, pin.SE3.Identity(), cube_shape)
    geom_box2.meshColor = np.array([0.2, 0.2, 0.8, 1.0])
    geom_model.addGeometryObject(geom_box2)

    # Floor - half-space  { p : p_z >= 0 } attached to the universe joint (id 0)
    floor_shape = coal.Halfspace(np.array([0.0, 0.0, 1.0]), 0.0)
    geom_floor = pin.GeometryObject("floor", 0, 0, pin.SE3.Identity(), floor_shape)
    geom_floor.meshColor = np.array([0.5, 0.5, 0.5, 0.5])
    geom_model.addGeometryObject(geom_floor)

except ImportError:
    print("coal not found - geometry model will be empty (solver still runs).")

# ─── 3. Initial configuration ─────────────────────────────────────────────────

# Free-flyer q = [tx, ty, tz,  qx, qy, qz, qw]
# Cube 1 centre at z = box_half  (bottom face touching z = 0)
# Cube 2 centre at z = 3*box_half (bottom face touching the top of cube 1)
q_box1 = np.array([0.0, 0.0, box_half, 0.0, 0.0, 0.0, 1.0])
q_box2 = np.array([0.0, 0.0, 3.0 * box_half, 0.0, 0.0, 0.0, 1.0])
q0 = np.concatenate([q_box1, q_box2])

v0 = np.zeros(model.nv)  # zero initial velocity
tau0 = np.zeros(model.nv)  # no external torques
dt = 1e-3  # time-step [s]

# ─── 4. Build the 8 contact constraints ──────────────────────────────────────

friction_coeff = 0.8


def make_corner_constraints(model, jid1, jid2, z1, z2):
    """
    Return 4 PointContactConstraintModel objects for the corners of a planar
    contact interface between two bodies.

    Parameters
    ----------
    jid1, jid2 : joint indices of the two contacting bodies
    z1         : z-coordinate of the contact plane in jid1's local frame
                 (e.g. +box_half for the top face, 0.0 for the world floor)
    z2         : z-coordinate of the contact plane in jid2's local frame
                 (e.g. -box_half for the bottom face)
    """
    corners_xy = [
        np.array([+box_half, +box_half]),
        np.array([-box_half, +box_half]),
        np.array([-box_half, -box_half]),
        np.array([+box_half, -box_half]),
    ]
    cms = []
    for xy in corners_xy:
        p1 = np.array([xy[0], xy[1], z1])
        p2 = np.array([xy[0], xy[1], z2])
        placement1 = pin.SE3(np.eye(3), p1)
        placement2 = pin.SE3(np.eye(3), p2)
        cm = pin.PointContactConstraintModel(model, jid1, placement1, jid2, placement2)
        cm.set = pin.CoulombFrictionCone(friction_coeff)
        cms.append(cm)
    return cms


constraint_models = pin.StdVec_ConstraintModel()

# 4 constraints: floor (universe, jid=0) ↔ cube 1 bottom face
# - in universe frame: contact points are on the floor plane  z = 0
# - in box1 local frame: contact points are at the bottom corners  z = -box_half
for cm in make_corner_constraints(model, 0, joint1_id, 0.0, -box_half):
    constraint_models.append(pin.ConstraintModel(cm))

# 4 constraints: cube 1 top face ↔ cube 2 bottom face
# - in box1 local frame: top corners at  z = +box_half
# - in box2 local frame: bottom corners at  z = -box_half
for cm in make_corner_constraints(model, joint1_id, joint2_id, +box_half, -box_half):
    constraint_models.append(pin.ConstraintModel(cm))

total_residual_size = sum(cm.residualSize() for cm in constraint_models)
print(f"Number of constraints:          {len(constraint_models)}")
print(f"Total constraint residual size: {total_residual_size}")

# ─── 5. Delassus operator and drift vector ────────────────────────────────────

data = model.createData()
fext = [pin.Force.Zero() for _ in range(model.njoints)]

# CRBA is required before building the Cholesky decomposition of the Delassus
# matrix G = J M⁻¹ Jᵀ.
pin.crba(model, data, q0, pin.Convention.WORLD)

# Free acceleration: velocity the system would reach in one time-step without contacts.
v_free = v0 + dt * pin.aba(model, data, q0, v0, tau0, fext)

# Initialise constraint data and evaluate constraint Jacobians at q0.
constraint_datas = pin.StdVec_ConstraintData()
for cmodel in constraint_models:
    cdata = cmodel.createData()
    cmodel.calc(model, data, cdata)
    constraint_datas.append(cdata)

# Cholesky decomposition of the Delassus matrix.
chol = pin.ConstraintCholeskyDecomposition(
    model, data, constraint_models, constraint_datas
)
chol.compute(model, data, constraint_models, constraint_datas, 1e-10)

# DelassusCholeskyExpression wraps the Cholesky factors for efficient solves.
delassus_expr = chol.getDelassusCholeskyExpression()

# Constraint Jacobian and drift  g = J v_free.
Jc = pin.getConstraintsJacobian(model, data, constraint_models, constraint_datas)
g = Jc @ v_free

print(f"Delassus matrix size:           {delassus_expr.matrix().shape}")
print(f"Drift vector ‖g‖:               {np.linalg.norm(g):.4e}")

# ─── 6. ADMM solver ───────────────────────────────────────────────────────────

solver = pin.ADMMConstraintSolver()

settings = pin.ADMMSolverSettings()
settings.max_iterations = 1000
settings.absolute_feasibility_tol = 1e-10
settings.relative_feasibility_tol = 1e-12
settings.absolute_complementarity_tol = 1e-10
settings.relative_complementarity_tol = 1e-12
settings.admm_update_rule = pin.ADMMUpdateRule.SPECTRAL
settings.mu_prox = 1e-6
settings.stat_record = (
    True  # per-iteration statistics. Turn off for faster solves if not needed.
)
settings.solve_ncp = True

result = pin.ADMMSolverResult()
has_converged = solver.solve(
    delassus_expr, g, constraint_models, constraint_datas, settings, result
)

# ─── 7. Results ───────────────────────────────────────────────────────────────

print()
print("── ADMM solver results ──────────────────────────────────────────────")
print(f"  Converged:            {result.converged}")
print(f"  Iterations:           {result.iterations}")
print(f"  Primal feasibility:   {result.primal_feasibility:.4e}")
print(f"  Dual feasibility:     {result.dual_feasibility:.4e}")
print(f"  Complementarity:      {result.complementarity:.4e}")
print(f"  Final rho:            {result.rho:.4e}")
print("─────────────────────────────────────────────────────────────────────")

impulses = result.retrieveConstraintImpulses()
velocities = result.retrieveConstraintVelocities()
print(f"\nConstraint impulses   ‖λ‖: {np.linalg.norm(impulses):.4e}")
print(f"Constraint velocities ‖v‖: {np.linalg.norm(velocities):.4e}")

if not has_converged:
    print("\nWarning: solver did not converge within the iteration budget.")
    sys.exit(1)
