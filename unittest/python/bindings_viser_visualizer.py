"""
Tests for pinocchio.visualize.ViserVisualizer.

Focused on the meshScale handling: URDF ``<mesh scale="..."/>`` must be
baked into the mesh vertices at load time, and must NOT be multiplied
into the frame translation in updatePlacements.
"""

import unittest
from pathlib import Path

import numpy as np
import pinocchio as pin

try:
    import viser
    from pinocchio.visualize import ViserVisualizer

    WITH_VISER = True
except ImportError:
    WITH_VISER = False


# Resolve the pinocchio repo's shipped test mesh.
MODELS_DIR = (
    Path(__file__).resolve().parent.parent.parent / "models"
)  # unittest/python/../../models -> models
BOX_STL = MODELS_DIR / "simple_humanoid_description" / "box.stl"


URDF_WITH_SCALE = """<?xml version="1.0"?>
<robot name="scaled_box">
  <link name="base_link">
    <visual>
      <origin xyz="1 2 3" rpy="0 0 0"/>
      <geometry>
        <mesh filename="package://box.stl" scale="2 3 4"/>
      </geometry>
    </visual>
  </link>
</robot>
"""


@unittest.skipUnless(WITH_VISER, "Needs viser")
@unittest.skipUnless(pin.WITH_URDFDOM, "Needs URDFDOM")
class TestViserVisualizerMeshScale(unittest.TestCase):
    def setUp(self):
        if not BOX_STL.is_file():
            self.skipTest(f"box.stl not found at {BOX_STL}")

        self.model = pin.buildModelFromXML(URDF_WITH_SCALE)
        self.visual_model = pin.buildGeomFromUrdfString(
            self.model,
            URDF_WITH_SCALE,
            pin.GeometryType.VISUAL,
            package_dirs=[str(BOX_STL.parent)],
        )
        # Reuse visual as collision; the ViserVisualizer updates both on load.
        self.collision_model = self.visual_model
        self.expected_scale = np.array([2.0, 3.0, 4.0])
        self.expected_translation = np.array([1.0, 2.0, 3.0])

        geom = self.visual_model.geometryObjects[0]
        np.testing.assert_allclose(geom.meshScale, self.expected_scale)

    def _make_visualizer(self):
        # A headless viser server on an ephemeral port.
        server = viser.ViserServer(host="localhost", server_port=0)
        viz = ViserVisualizer(
            self.model,
            collision_model=self.collision_model,
            visual_model=self.visual_model,
        )
        viz.initViewer(viewer=server)
        return server, viz

    def test_scale_is_forwarded_to_mesh_loader(self):
        """
        Check that ``_add_mesh_from_path`` receives the geometry meshScale.
        """
        captured = {}
        orig = ViserVisualizer._add_mesh_from_path

        def spy(self, name, mesh_path, color, scale=None):
            captured["scale"] = None if scale is None else np.array(scale)
            return orig(self, name, mesh_path, color, scale=scale)

        ViserVisualizer._add_mesh_from_path = spy
        try:
            server, viz = self._make_visualizer()
            try:
                viz.loadViewerModel()
            finally:
                server.stop()
        finally:
            ViserVisualizer._add_mesh_from_path = orig

        self.assertIsNotNone(captured.get("scale"))
        np.testing.assert_allclose(captured["scale"], self.expected_scale)

    def test_translation_is_not_multiplied_by_mesh_scale(self):
        """
        Check that frame.position equals the link translation.

        After ``display(q)``, the viser frame translation must match the
        URDF link origin, *not* origin * meshScale (the previous bug).
        """
        server, viz = self._make_visualizer()
        try:
            viz.loadViewerModel()
            viz.display(pin.neutral(self.model))

            geom_name = self.visual_model.geometryObjects[0].name
            frame = viz.frames[viz.viewerRootNodeName + "/visual/" + geom_name]

            np.testing.assert_allclose(
                np.asarray(frame.position), self.expected_translation, atol=1e-9
            )
            buggy = self.expected_translation * self.expected_scale
            self.assertFalse(np.allclose(np.asarray(frame.position), buggy))
        finally:
            server.stop()


if __name__ == "__main__":
    unittest.main()
