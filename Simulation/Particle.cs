using System.Numerics;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// Distinguishes fluid particles from static boundary particles.
/// </summary>
public enum ParticleType : byte
{
    Fluid,
    Boundary,
}

/// <summary>
/// Core data for a single simulation particle.
/// Uses System.Numerics types to keep the simulation layer independent of Godot.
/// </summary>
public struct Particle
{
    public ParticleType Type;
    public Vector3 Position;
    public Vector3 Velocity;
    public Vector3 Acceleration;
    public float Density;
    public float Pressure;
    public float Mass;

    /// <summary>
    /// True if this is a fluid particle (participates in dynamics).
    /// False for static boundary particles.
    /// </summary>
    public readonly bool IsFluid => Type == ParticleType.Fluid;
}
