using System.Collections.Generic;
using Godot;
using PhysicsSimulator.Simulation;

namespace PhysicsSimulator.Rendering;

/// <summary>
/// Renders simulation particles using MultiMeshInstance3D.
/// This is a data-oriented approach: one draw call for all particles, rather than one Node3D per particle.
/// </summary>
public partial class SimulationRenderer : MultiMeshInstance3D
{
    private MultiMesh _multiMesh = null!;

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
}
