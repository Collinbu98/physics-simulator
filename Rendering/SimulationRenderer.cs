using System.Collections.Generic;
using Godot;
using PhysicsSimulator.Simulation;

namespace PhysicsSimulator.Rendering;

/// <summary>
/// Renders simulation particles using MultiMeshInstance3D.
/// Fluid particles are rendered as opaque spheres; boundary particles as
/// smaller, semi-transparent spheres so they are visible for debugging
/// without visually dominating the fluid.
/// Also renders a wireframe representation of the container boundaries.
/// </summary>
public partial class SimulationRenderer : MultiMeshInstance3D
{
    private MultiMesh _multiMesh = null!;
    private MultiMeshInstance3D? _boundaryNode;
    private MultiMesh? _boundaryMesh;
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

        // Boundary particle multimesh — smaller, semi-transparent, distinct color
        _boundaryNode = new MultiMeshInstance3D { Name = "BoundaryParticles" };
        AddChild(_boundaryNode);

        _boundaryMesh = new MultiMesh
        {
            TransformFormat = MultiMesh.TransformFormatEnum.Transform3D,
            InstanceCount = 0,
            Mesh = new SphereMesh
            {
                Radius = 0.012f,
                Height = 0.024f,
                RadialSegments = 6,
                Rings = 3,
            },
        };

        _boundaryNode.Multimesh = _boundaryMesh;
        _boundaryNode.MaterialOverride = new StandardMaterial3D
        {
            AlbedoColor = new Color(0.3f, 0.5f, 0.8f, 0.35f),
            ShadingMode = BaseMaterial3D.ShadingModeEnum.Unshaded,
            Transparency = BaseMaterial3D.TransparencyEnum.Alpha,
            CullMode = BaseMaterial3D.CullModeEnum.Disabled,
        };
    }

    /// <summary>
    /// Updates the visual representation to match the current particle positions.
    /// Fluid and boundary particles are rendered on separate multimeshes.
    /// Call this each frame after the simulation step.
    /// </summary>
    public void UpdateParticles(IReadOnlyList<Particle> particles)
    {
        if (_multiMesh == null || _boundaryMesh == null || _boundaryNode == null)
            return;

        // Separate fluid and boundary particles
        int fluidCount = 0;
        int boundaryCount = 0;
        for (int i = 0; i < particles.Count; i++)
        {
            if (particles[i].IsFluid)
                fluidCount++;
            else
                boundaryCount++;
        }

        // Fluid multimesh
        if (fluidCount > _multiMesh.InstanceCount)
            _multiMesh.InstanceCount = fluidCount;

        int fi = 0;
        for (int i = 0; i < particles.Count; i++)
        {
            if (!particles[i].IsFluid) continue;
            var pos = particles[i].Position;
            _multiMesh.SetInstanceTransform(fi, new Transform3D(Basis.Identity, new Vector3(pos.X, pos.Y, pos.Z)));
            fi++;
        }

        // Boundary multimesh
        if (boundaryCount > _boundaryMesh.InstanceCount)
            _boundaryMesh.InstanceCount = boundaryCount;

        int bi = 0;
        for (int i = 0; i < particles.Count; i++)
        {
            if (particles[i].IsFluid) continue;
            var pos = particles[i].Position;
            _boundaryMesh.SetInstanceTransform(bi, new Transform3D(Basis.Identity, new Vector3(pos.X, pos.Y, pos.Z)));
            bi++;
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
