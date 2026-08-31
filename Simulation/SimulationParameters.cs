namespace PhysicsSimulator.Simulation;

/// <summary>
/// Tunable constants for the simulation.
/// Separated from the simulation class so they can be easily adjusted or exposed in an editor.
/// </summary>
public class SimulationParameters
{
    /// <summary>
    /// Fixed timestep for the simulation in seconds.
    /// SPH simulations typically use a fixed timestep for stability.
    /// </summary>
    public float TimeStep { get; set; } = 1.0f / 60.0f;

    /// <summary>
    /// Rest density of the fluid in kg/m^3.
    /// Water is approximately 1000 kg/m^3.
    /// </summary>
    public float RestDensity { get; set; } = 1000.0f;

    /// <summary>
    /// Stiffness constant for the equation of state (pressure calculation).
    /// </summary>
    public float Stiffness { get; set; } = 200.0f;

    /// <summary>
    /// Dynamic viscosity coefficient.
    /// Water is approximately 0.001 Pa·s.
    /// </summary>
    public float Viscosity { get; set; } = 0.1f;

    /// <summary>
    /// Particle mass in kg. Used for SPH density and force calculations.
    /// </summary>
    public float ParticleMass { get; set; } = 0.02f;

    /// <summary>
    /// Smoothness radius (smoothing kernel support) in meters.
    /// Defines the neighborhood for SPH interpolation.
    /// </summary>
    public float SmoothingRadius { get; set; } = 0.1f;

    /// <summary>
    /// Gravity vector in m/s^2.
    /// </summary>
    public float Gravity { get; set; } = -9.81f;

    /// <summary>
    /// Global speed multiplier for the simulation.
    /// Set to 0 to pause the simulation.
    /// </summary>
    public float TimeScale { get; set; } = 1.0f;
}
