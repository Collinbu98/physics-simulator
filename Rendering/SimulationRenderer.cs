using System.Collections.Generic;
using Godot;
using PhysicsSimulator.Simulation;

namespace PhysicsSimulator.Rendering;

/// <summary>
/// Renders simulation particles using MultiMeshInstance3D.
/// Also renders a wireframe representation of the container boundaries.
/// </summary>
public partial class SimulationRenderer : MultiMeshInstance3D
{
    private MultiMesh _multiMesh = null!;
    private MeshInstance3D? _wireframeNode;

    public SimulationRenderer()
    {
        Name = "SimulationRenderer";
    }

    public override void _Ready()
    {
        _multiMesh = new MultiMesh
        {
            TransformFormat = MultiMesh.TransformFormatEnum.Transform3D,
            InstanceCount = 0,
            Mesh = new SphereMesh
            {
                Radius = 0.02f,
                Height = 0.04f,
                RadialSegments = 8,
                Rings = 4,
            },
        };

        Multimesh = _multiMesh;
    }

    /// <summary>
    /// Updates the visual representation to match the current particle positions.
    /// Call this each frame after the simulation step.
    /// </summary>
    public void UpdateParticles(IReadOnlyList<Particle> particles)
    {
        if (_multiMesh == null)
            return;

        int count = particles.Count;

        // Resize the multimesh if needed (grow but don't shrink)
        if (count > _multiMesh.InstanceCount)
        {
            _multiMesh.InstanceCount = count;
        }

        for (int i = 0; i < count; i++)
        {
            var pos = particles[i].Position;
            _multiMesh.SetInstanceTransform(i, new Transform3D(Basis.Identity, new Vector3(pos.X, pos.Y, pos.Z)));
        }
    }

    /// <summary>
    /// Creates or updates a wireframe box showing the container boundaries.
    /// The box is drawn with bottom at y = 0, extending upward to height,
    /// and centered at the origin in X and Z.
    /// </summary>
    public void UpdateContainerWireframe(float halfX, float height, float halfZ)
    {
        if (_wireframeNode == null)
        {
            _wireframeNode = new MeshInstance3D { Name = "ContainerWireframe" };
            AddChild(_wireframeNode);
        }

        var mesh = new ImmediateMesh();
        mesh.SurfaceBegin(Mesh.PrimitiveType.Lines);

        // Bottom face (y = 0)
        mesh.SurfaceAddVertex(new Vector3(-halfX, 0.0f, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, 0.0f, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, 0.0f, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, 0.0f, halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, 0.0f, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, 0.0f, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, 0.0f, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, 0.0f, -halfZ));

        // Top face (y = height)
        mesh.SurfaceAddVertex(new Vector3(-halfX, height, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, height, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, height, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, height, halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, height, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, height, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, height, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, height, -halfZ));

        // Vertical edges connecting bottom to top
        mesh.SurfaceAddVertex(new Vector3(-halfX, 0.0f, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, height, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, 0.0f, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, height, -halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, 0.0f, halfZ));
        mesh.SurfaceAddVertex(new Vector3(halfX, height, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, 0.0f, halfZ));
        mesh.SurfaceAddVertex(new Vector3(-halfX, height, halfZ));

        mesh.SurfaceEnd();

        _wireframeNode.Mesh = mesh;
        _wireframeNode.MaterialOverride = new StandardMaterial3D
        {
            AlbedoColor = new Color(0.5f, 0.8f, 1.0f),
            ShadingMode = BaseMaterial3D.ShadingModeEnum.Unshaded,
        };
    }
}
