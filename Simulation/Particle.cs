using System.Numerics;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// Core data for a single simulation particle.
/// Uses System.Numerics types to keep the simulation layer independent of Godot.
/// </summary>
public struct Particle
{
    public Vector3 Position;
    public Vector3 Velocity;
    public Vector3 Acceleration;
    public float Density;
    public float Pressure;
    public float Mass;
}
