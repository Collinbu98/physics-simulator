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
	public float TimeStep { get; set; } = 0.0005f;

	/// <summary>
	/// Rest density of the fluid in kg/m^3.
	/// Water is approximately 1000 kg/m^3.
	/// </summary>
	public float RestDensity { get; set; } = 1000.0f;

	/// <summary>
	/// Pressure stiffness constant (k) for the equation of state: P = k * (rho - rho_0).
	/// Higher values make the fluid more incompressible but may cause instability.
	/// This is an initial value for experimentation; tune as needed.
	/// </summary>
	public float PressureStiffness { get; set; } = 50.0f;

	/// <summary>
	/// Dynamic viscosity coefficient.
	/// Water is approximately 0.001 Pa·s.
	/// </summary>
	public float Viscosity { get; set; } = 0.1f;

	/// <summary>
	/// Kinematic viscosity for SPH viscous diffusion (m²/s).
	/// Controls how strongly velocity differences between neighboring particles are smoothed.
	/// The SPH viscosity force on particle i is:
	///   F_i^visc = Σ_j m_j × ν × (v_j - v_i)/ρ_j × ∇²W_visc(r_ij, h)
	/// where ∇²W_visc is the viscosity kernel Laplacian and ν is this parameter.
	/// Higher values cause faster velocity equalization (more "syrupy" fluid).
	/// Water is approximately 1e-6 m²/s; values of 1e-5 to 1e-4 produce more
	/// visible damping in small-scale particle simulations.
	/// </summary>
	public float KinematicViscosity { get; set; } = 1e-4f;

	/// <summary>
	/// Particle mass in kg. Used for SPH density and force calculations.
	/// </summary>
	public float ParticleMass { get; set; } = 0.27f;

	/// <summary>
	/// Smoothness radius (smoothing kernel support) in meters.
	/// Defines the neighborhood for SPH interpolation.
	/// </summary>
	public float SmoothingRadius { get; set; } = 0.075f;

	/// <summary>
	/// Gravity vector in m/s^2.
	/// </summary>
	public float Gravity { get; set; } = -9.81f;

	/// <summary>
	/// Global speed multiplier for the simulation.
	/// Set to 0 to pause the simulation.
	/// </summary>
	public float TimeScale { get; set; } = 1.0f;

	// ── Container boundaries ───────────────────────────────────────────

	/// <summary>
	/// Width of the container in meters (X axis).
	/// The container extends from -Width/2 to +Width/2 in X.
	/// Default of 0.6 m comfortably contains the initial 4×4×4 particle block.
	/// </summary>
	public float ContainerWidth { get; set; } = 0.6f;

	/// <summary>
	/// Height of the container in meters (Y axis).
	/// The container extends from 0 to Height in Y (bottom at y = 0, open top).
	/// </summary>
	public float ContainerHeight { get; set; } = 0.6f;

	/// <summary>
	/// Depth of the container in meters (Z axis).
	/// The container extends from -Depth/2 to +Depth/2 in Z.
	/// </summary>
	public float ContainerDepth { get; set; } = 0.6f;

	/// <summary>
	/// Coefficient of restitution for boundary collisions.
	/// 0 = perfectly inelastic (no bounce), 1 = perfectly elastic (full bounce).
	/// A low value like 0.3 allows particles to settle without excessive bounciness.
	/// Applied only to the velocity component normal to the wall.
	/// </summary>
	public float BoundaryRestitution { get; set; } = 0.3f;
}
