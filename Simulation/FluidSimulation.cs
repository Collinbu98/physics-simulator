using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// Core simulation class. Owns particle state and advances it each step.
/// This class has no dependency on Godot or any other engine.
/// </summary>
public class FluidSimulation
{
    private readonly List<Particle> _particles = new();
    private readonly SimulationParameters _parameters;
    private readonly SpatialGrid _grid;
    private readonly List<int> _neighborBuffer = new();
    private float _timeAccumulator;
    private int _boundaryParticleCount;

    // Diagnostic storage for pressure force magnitudes (updated each step)
    private float _minPressureForceMagnitude = float.MaxValue;
    private float _maxPressureForceMagnitude = float.MinValue;
    private float _totalPressureForceMagnitude;
    private int _pressureForceCount;
    private bool _pressureForceNonFinite;
    private float _maxAccelerationMagnitude;
    private float _maxVelocityMagnitude;

    // Diagnostic storage for viscosity force magnitudes (updated each step)
    private float _minViscosityForceMagnitude = float.MaxValue;
    private float _maxViscosityForceMagnitude = float.MinValue;
    private float _totalViscosityForceMagnitude;
    private int _viscosityForceCount;

    private long _stepCount;
    private int _lastBoundaryCollisions;

    // Stage-level profiling accumulators
    private const int ProfileInterval = 2000;
    private long _profileSteps;
    private double _profileGridMs;
    private double _profileDensityMs;
    private double _profilePressureMs;
    private double _profilePressureForceMs;
    private double _profileViscosityMs;
    private double _profileIntegrationMs;
    private double _profileBoundaryMs;
    private readonly Stopwatch _sw = new();

    /// <summary>
    /// Set to a profile report string after every ProfileInterval steps, then cleared by reader.
    /// SimulationNode polls this and prints via GD.Print.
    /// </summary>
    public string? LastProfileReport { get; private set; }

    /// <summary>
    /// Clears the profile report after it has been consumed by the caller.
    /// </summary>
    public void ClearProfileReport() => LastProfileReport = null;

    /// <summary>
    /// Read-only access to particle data.
    /// </summary>
    public IReadOnlyList<Particle> Particles => _particles;

    /// <summary>
    /// The current number of particles in the simulation (fluid + boundary).
    /// </summary>
    public int ParticleCount => _particles.Count;

    /// <summary>
    /// Number of static boundary particles. Remaining particles are fluid.
    /// </summary>
    public int BoundaryParticleCount => _boundaryParticleCount;

    /// <summary>
    /// Number of fluid (dynamic) particles.
    /// </summary>
    public int FluidParticleCount => _particles.Count - _boundaryParticleCount;

    public SimulationParameters Parameters => _parameters;

    /// <summary>
    /// Total number of fixed timesteps executed since the last Reset().
    /// </summary>
    public long StepCount => _stepCount;

    /// <summary>
    /// Number of boundary collisions in the most recent fixed step.
    /// </summary>
    public int LastBoundaryCollisions => _lastBoundaryCollisions;

    /// <summary>
    /// The spatial grid used for neighbor searching. Rebuilt each fixed step.
    /// </summary>
    public SpatialGrid Grid => _grid;

    public FluidSimulation(SimulationParameters parameters)
    {
        _parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));
        _grid = new SpatialGrid(_parameters.SmoothingRadius);
    }

    /// <summary>
    /// Resets the simulation to an empty state.
    /// </summary>
    public void Reset()
    {
        _particles.Clear();
        _boundaryParticleCount = 0;
        _grid.Clear();
        _timeAccumulator = 0.0f;
        _stepCount = 0;
        _profileSteps = 0;
        _profileGridMs = 0;
        _profileDensityMs = 0;
        _profilePressureMs = 0;
        _profilePressureForceMs = 0;
        _profileViscosityMs = 0;
        _profileIntegrationMs = 0;
        _profileBoundaryMs = 0;
    }

    /// <summary>
    /// Adds a single fluid particle with the given state.
    /// </summary>
    public void AddParticle(Vector3 position, Vector3 velocity, float mass)
    {
        _particles.Add(new Particle
        {
            Type = ParticleType.Fluid,
            Position = position,
            Velocity = velocity,
            Acceleration = Vector3.Zero,
            Density = 0.0f,
            Pressure = 0.0f,
            Mass = mass,
        });
    }

    /// <summary>
    /// Moves a single particle to a new position. Used by diagnostics that need
    /// to perturb one particle after spawning a clean lattice.
    /// </summary>
    public void SetParticlePosition(int index, Vector3 newPosition)
    {
        var p = _particles[index];
        p.Position = newPosition;
        _particles[index] = p;
    }

    /// <summary>
    /// Generates static boundary particles along the five interior faces of the
    /// container (bottom, left, right, front, back). The top is open.
    /// Boundary particles participate in SPH neighbor queries so fluid particles
    /// near a wall receive density and pressure support from the solid boundary.
    /// Particles are placed at <see cref="SimulationParameters.BoundaryParticleSpacing"/>
    /// intervals, offset inward from each wall by half the spacing.
    /// </summary>
    public void GenerateBoundaryParticles()
    {
        float spacing = _parameters.BoundaryParticleSpacing;
        float halfSpacing = spacing * 0.5f;
        float halfW = _parameters.ContainerWidth * 0.5f;
        float halfD = _parameters.ContainerDepth * 0.5f;
        float mass = _parameters.ParticleMass;
        float boundaryDensity = _parameters.RestDensity;
        int countBeforeBoundary = _particles.Count;

        // Bottom wall (y = 0) — particles placed at y = -halfSpacing (behind the collision plane)
        for (float x = -halfW + halfSpacing; x < halfW; x += spacing)
        for (float z = -halfD + halfSpacing; z < halfD; z += spacing)
        {
            _particles.Add(new Particle
            {
                Type = ParticleType.Boundary,
                Position = new Vector3(x, -halfSpacing, z),
                Velocity = Vector3.Zero,
                Acceleration = Vector3.Zero,
                Density = boundaryDensity,
                Pressure = 0.0f,
                Mass = mass,
            });
        }

        // Left wall (x = -halfW)
        for (float y = halfSpacing; y < _parameters.ContainerHeight; y += spacing)
        for (float z = -halfD + halfSpacing; z < halfD; z += spacing)
        {
            _particles.Add(new Particle
            {
                Type = ParticleType.Boundary,
                Position = new Vector3(-halfW + halfSpacing, y, z),
                Velocity = Vector3.Zero,
                Acceleration = Vector3.Zero,
                Density = boundaryDensity,
                Pressure = 0.0f,
                Mass = mass,
            });
        }

        // Right wall (x = +halfW)
        for (float y = halfSpacing; y < _parameters.ContainerHeight; y += spacing)
        for (float z = -halfD + halfSpacing; z < halfD; z += spacing)
        {
            _particles.Add(new Particle
            {
                Type = ParticleType.Boundary,
                Position = new Vector3(halfW - halfSpacing, y, z),
                Velocity = Vector3.Zero,
                Acceleration = Vector3.Zero,
                Density = boundaryDensity,
                Pressure = 0.0f,
                Mass = mass,
            });
        }

        // Front wall (z = -halfD)
        for (float x = -halfW + halfSpacing; x < halfW; x += spacing)
        for (float y = halfSpacing; y < _parameters.ContainerHeight; y += spacing)
        {
            _particles.Add(new Particle
            {
                Type = ParticleType.Boundary,
                Position = new Vector3(x, y, -halfD + halfSpacing),
                Velocity = Vector3.Zero,
                Acceleration = Vector3.Zero,
                Density = boundaryDensity,
                Pressure = 0.0f,
                Mass = mass,
            });
        }

        // Back wall (z = +halfD)
        for (float x = -halfW + halfSpacing; x < halfW; x += spacing)
        for (float y = halfSpacing; y < _parameters.ContainerHeight; y += spacing)
        {
            _particles.Add(new Particle
            {
                Type = ParticleType.Boundary,
                Position = new Vector3(x, y, halfD - halfSpacing),
                Velocity = Vector3.Zero,
                Acceleration = Vector3.Zero,
                Density = boundaryDensity,
                Pressure = 0.0f,
                Mass = mass,
            });
        }

        _boundaryParticleCount = _particles.Count - countBeforeBoundary;
    }

    /// <summary>
    /// Advances the simulation by the given real-time delta (in seconds).
    /// Uses a fixed-timestep accumulator for deterministic behavior.
    /// Does nothing if TimeScale is zero.
    /// </summary>
    public void Step(float deltaTime)
    {
        if (_parameters.TimeScale <= 0.0f)
            return;

        _timeAccumulator += deltaTime * _parameters.TimeScale;

        while (_timeAccumulator >= _parameters.TimeStep)
        {
            FixedStep(_parameters.TimeStep);
            _timeAccumulator -= _parameters.TimeStep;
        }
    }

    /// <summary>
    /// Advances the simulation by exactly one fixed timestep, regardless of
    /// TimeScale or the accumulator. Intended for single-step debugging while paused.
    /// </summary>
    public void StepOnce()
    {
        FixedStep(_parameters.TimeStep);
    }

    /// <summary>
    /// One fixed-size simulation step. Computes SPH pressure forces and integrates motion.
    /// Pipeline: grid → density → pressure → pressure force → acceleration → velocity → position
    /// </summary>
    private void FixedStep(float dt)
    {
        _stepCount++;
        _lastBoundaryCollisions = 0;

        // --- Profiled stages ---
        _sw.Restart();
        RebuildGrid();
        _sw.Stop();
        _profileGridMs += _sw.Elapsed.TotalMilliseconds;

        _sw.Restart();
        ComputeAllDensities();
        _sw.Stop();
        _profileDensityMs += _sw.Elapsed.TotalMilliseconds;

        _sw.Restart();
        ComputeAllPressures();
        _sw.Stop();
        _profilePressureMs += _sw.Elapsed.TotalMilliseconds;

        _sw.Restart();
        ComputeAllPressureForces();
        _sw.Stop();
        _profilePressureForceMs += _sw.Elapsed.TotalMilliseconds;

        _sw.Restart();
        ComputeAllViscosityForces();
        _sw.Stop();
        _profileViscosityMs += _sw.Elapsed.TotalMilliseconds;

        _sw.Restart();
        var gravity = new Vector3(0.0f, _parameters.Gravity, 0.0f);
        _maxVelocityMagnitude = 0.0f;
        int fluidCount = _particles.Count - _boundaryParticleCount;
        for (int i = 0; i < fluidCount; i++)
        {
            var p = _particles[i];
            p.Acceleration += gravity;
            p.Velocity += p.Acceleration * dt;
            p.Position += p.Velocity * dt;
            _particles[i] = p;
            float velMag = p.Velocity.Length();
            if (velMag > _maxVelocityMagnitude)
                _maxVelocityMagnitude = velMag;
        }
        _sw.Stop();
        _profileIntegrationMs += _sw.Elapsed.TotalMilliseconds;

        _sw.Restart();
        HandleBoundaryCollisions();
        _sw.Stop();
        _profileBoundaryMs += _sw.Elapsed.TotalMilliseconds;

        _profileSteps++;

        if (_profileSteps >= ProfileInterval)
        {
            double grid = _profileGridMs;
            double density = _profileDensityMs;
            double pressure = _profilePressureMs;
            double pressF = _profilePressureForceMs;
            double visc = _profileViscosityMs;
            double integ = _profileIntegrationMs;
            double bndry = _profileBoundaryMs;
            double total = grid + density + pressure + pressF + visc + integ + bndry;
            double avg = total / _profileSteps;
            double pct(double v) => total > 0 ? v / total * 100.0 : 0;

            LastProfileReport =
                $"[Profile] {FluidParticleCount}f+{_boundaryParticleCount}b | {_profileSteps} steps\n" +
                $"  Grid:         {grid,8:F2} ms  ({pct(grid):F1}%)\n" +
                $"  Density:      {density,8:F2} ms  ({pct(density):F1}%)\n" +
                $"  Pressure:     {pressure,8:F2} ms  ({pct(pressure):F1}%)\n" +
                $"  PressureF:    {pressF,8:F2} ms  ({pct(pressF):F1}%)\n" +
                $"  Viscosity:    {visc,8:F2} ms  ({pct(visc):F1}%)\n" +
                $"  Integration:  {integ,8:F2} ms  ({pct(integ):F1}%)\n" +
                $"  Boundary:     {bndry,8:F2} ms  ({pct(bndry):F1}%)\n" +
                $"  ─────────────────────────────────\n" +
                $"  Total:        {total,8:F2} ms  avg={avg:F4} ms/step\n" +
                $"  Steps/sec:    {_profileSteps / (total / 1000.0),8:F0}  (target {(int)(1.0 / _parameters.TimeStep)})";

            _profileSteps = 0;
            _profileGridMs = 0;
            _profileDensityMs = 0;
            _profilePressureMs = 0;
            _profilePressureForceMs = 0;
            _profileViscosityMs = 0;
            _profileIntegrationMs = 0;
            _profileBoundaryMs = 0;
        }
    }

    /// <summary>
    /// Clears and repopulates the spatial grid from current particle positions.
    /// Must be called at the start of every fixed step before any queries.
    /// </summary>
    private void RebuildGrid()
    {
        _grid.Clear();

        for (int i = 0; i < _particles.Count; i++)
            _grid.Insert(i, _particles[i].Position);
    }

    /// <summary>
    /// Computes density for all particles (fluid and boundary) using the Poly6 kernel.
    /// Boundary particles need density computed so they are recognized as valid
    /// neighbors in force calculations (the density guard checks pj.Density > 0).
    /// Must be called after RebuildGrid so the grid is current.
    /// Public so external code (e.g. diagnostics) can trigger density computation
    /// without running a full simulation step.
    /// </summary>
    public void ComputeAllDensities()
    {
        float h = _parameters.SmoothingRadius;
        float h2 = h * h;
        float mass = _parameters.ParticleMass;
        float poly6Coeff = SPHKernels.Poly6Coefficient(h);
        float poly6Origin = SPHKernels.Poly6AtOrigin(h);

        for (int i = 0; i < _particles.Count; i++)
        {
            // Self-contribution: mass * W(0, h)
            float density = mass * poly6Origin;

            // Neighbor contributions (GetNeighbors excludes self)
            GetNeighbors(i, _neighborBuffer);

            Vector3 posI = _particles[i].Position;
            for (int n = 0; n < _neighborBuffer.Count; n++)
            {
                int j = _neighborBuffer[n];
                Vector3 diff = _particles[j].Position - posI;
                float distSq = diff.LengthSquared();
                density += mass * SPHKernels.Poly6(distSq, h2, poly6Coeff);
            }

            var p = _particles[i];
            p.Density = density;
            _particles[i] = p;
        }
    }

    /// <summary>
    /// Computes pressure for fluid particles using the equation of state: P = k * (rho - rho_0).
    /// Boundary particles keep their pre-assigned density and zero pressure.
    /// Must be called after ComputeAllDensities so that Particle.Density is populated.
    /// Does not affect velocity, acceleration, or position.
    /// </summary>
    public void ComputeAllPressures()
    {
        float k = _parameters.PressureStiffness;
        float rho0 = _parameters.RestDensity;

        int fluidCount = _particles.Count - _boundaryParticleCount;
        for (int i = 0; i < fluidCount; i++)
        {
            var p = _particles[i];
            p.Pressure = k * (p.Density - rho0);
            _particles[i] = p;
        }
    }

    /// <summary>
    /// Computes SPH pressure forces for fluid particles using the symmetric formulation:
    ///   F_i = Σ_j m_j × (P_i/ρ_i² + P_j/ρ_j²) × ∇W_spiky(r_ij, h)
    ///
    /// Only fluid particles are processed in the outer loop. Boundary particles
    /// appear as neighbors: they contribute via their pre-assigned rest density
    /// and zero pressure. Since P_j = 0 for boundary particles, a fluid particle
    /// near a wall receives a pressure push away from the solid material.
    /// Boundary particles never receive force or acceleration.
    /// Must be called after ComputeAllDensities and ComputeAllPressures.
    /// Stores the force-per-unit-mass (acceleration) in Particle.Acceleration.
    /// </summary>
    public void ComputeAllPressureForces()
    {
        float h = _parameters.SmoothingRadius;
        float h2 = h * h;
        float mass = _parameters.ParticleMass;
        float spikyCoeff = SPHKernels.SpikyGradientCoefficient(h);

        // Reset diagnostic accumulators
        _minPressureForceMagnitude = float.MaxValue;
        _maxPressureForceMagnitude = float.MinValue;
        _totalPressureForceMagnitude = 0.0f;
        _pressureForceCount = 0;
        _pressureForceNonFinite = false;
        _maxAccelerationMagnitude = 0.0f;
        _maxVelocityMagnitude = 0.0f;

        // Only compute pressure forces for fluid particles.
        // Boundary particles (indices >= fluidCount) are static and never receive forces.
        int fluidCount = _particles.Count - _boundaryParticleCount;
        for (int i = 0; i < fluidCount; i++)
        {
            var pi = _particles[i];

            // Guard: zero or negative density → skip force (will be caught by diagnostic)
            if (pi.Density <= 0.0f)
            {
                pi.Acceleration = Vector3.Zero;
                _particles[i] = pi;
                continue;
            }

            float rhoI = pi.Density;
            float pOverRhoI2 = pi.Pressure / (rhoI * rhoI);

            Vector3 totalForce = Vector3.Zero;

            GetNeighbors(i, _neighborBuffer);
            Vector3 posI = pi.Position;

            for (int n = 0; n < _neighborBuffer.Count; n++)
            {
                int j = _neighborBuffer[n];
                var pj = _particles[j];

                // Guard: zero or negative density for neighbor
                if (pj.Density <= 0.0f)
                    continue;

                float rhoJ = pj.Density;
                float pOverRhoJ2 = pj.Pressure / (rhoJ * rhoJ);

                // r_ij = r_j - r_i (points from i toward j)
                Vector3 rVec = pj.Position - posI;

                // Spiky gradient: ∇W_spiky(r_ij, h)
                Vector3 gradW = SPHKernels.SpikyGradient(rVec, h2, spikyCoeff);

                // Symmetric pressure force term:
                // F_i += m_j × (P_i/ρ_i² + P_j/ρ_j²) × ∇W_spiky
                totalForce += mass * (pOverRhoI2 + pOverRhoJ2) * gradW;
            }

            // Convert force to acceleration: a = F / m_i
            Vector3 acceleration = totalForce / pi.Mass;

            // Diagnostic: track magnitudes
            float forceMag = totalForce.Length();
            if (_pressureForceCount == 0 || forceMag < _minPressureForceMagnitude)
                _minPressureForceMagnitude = forceMag;
            if (forceMag > _maxPressureForceMagnitude)
                _maxPressureForceMagnitude = forceMag;
            _totalPressureForceMagnitude += forceMag;
            _pressureForceCount++;

            float accelMag = acceleration.Length();
            if (accelMag > _maxAccelerationMagnitude)
                _maxAccelerationMagnitude = accelMag;

            // Check for non-finite values
            if (!float.IsFinite(acceleration.X) || !float.IsFinite(acceleration.Y) || !float.IsFinite(acceleration.Z))
                _pressureForceNonFinite = true;
            if (!float.IsFinite(totalForce.Length()))
                _pressureForceNonFinite = true;

            // Track velocity magnitude (needed when F4 calls this standalone)
            float velMag = pi.Velocity.Length();
            if (velMag > _maxVelocityMagnitude)
                _maxVelocityMagnitude = velMag;

            pi.Acceleration = acceleration;
            _particles[i] = pi;
        }
    }

    /// <summary>
    /// Computes SPH viscosity forces for fluid particles using the standard formulation:
    ///   a_i^visc = (1/ρ_i) × Σ_j m_j × ν × (v_j - v_i)/ρ_j × ∇²W_visc(r_ij, h)
    ///
    /// Only fluid particles are processed in the outer loop. Boundary particles
    /// appear as neighbors with zero velocity, producing a no-slip wall condition:
    /// fluid particles near a wall are decelerated toward the wall's zero velocity.
    /// Must be called after RebuildGrid and ComputeAllDensities.
    /// Adds the viscous acceleration to the existing Particle.Acceleration (which already contains
    /// pressure acceleration from ComputeAllPressureForces).
    /// </summary>
    public void ComputeAllViscosityForces()
    {
        float h = _parameters.SmoothingRadius;
        float h2 = h * h;
        float mass = _parameters.ParticleMass;
        float nu = _parameters.KinematicViscosity;
        float viscLapCoeff = SPHKernels.ViscosityLaplacianCoefficient(h);

        // Reset viscosity diagnostic accumulators
        _minViscosityForceMagnitude = float.MaxValue;
        _maxViscosityForceMagnitude = float.MinValue;
        _totalViscosityForceMagnitude = 0.0f;
        _viscosityForceCount = 0;

        // Only compute viscosity for fluid particles. Boundary particles are static.
        int fluidCount = _particles.Count - _boundaryParticleCount;
        for (int i = 0; i < fluidCount; i++)
        {
            var pi = _particles[i];

            // Guard: zero or negative density → skip
            if (pi.Density <= 0.0f)
                continue;

            float rhoI = pi.Density;
            Vector3 posI = pi.Position;
            Vector3 velI = pi.Velocity;
            Vector3 viscAccel = Vector3.Zero;

            GetNeighbors(i, _neighborBuffer);

            for (int n = 0; n < _neighborBuffer.Count; n++)
            {
                int j = _neighborBuffer[n];
                var pj = _particles[j];

                // Guard: zero or negative density for neighbor
                if (pj.Density <= 0.0f)
                    continue;

                float rhoJ = pj.Density;
                Vector3 rVec = pj.Position - posI;
                float dist = rVec.Length();

                // Viscosity kernel Laplacian: ∇²W_visc(r, h) = 45/(π h⁶) × (h - |r|)
                float laplacian = SPHKernels.ViscosityLaplacian(dist, h, viscLapCoeff);

                // Velocity difference: v_j - v_i
                Vector3 velDiff = pj.Velocity - velI;

                // Viscosity term: m_j × ν × (v_j - v_i)/ρ_j × ∇²W
                viscAccel += (mass * nu / rhoJ) * velDiff * laplacian;
            }

            // Divide by own density to complete the SPH average
            viscAccel = viscAccel / rhoI;

            // Diagnostic: track force magnitude (convert accel to force: F = a × m)
            float forceMag = (viscAccel * pi.Mass).Length();
            if (_viscosityForceCount == 0 || forceMag < _minViscosityForceMagnitude)
                _minViscosityForceMagnitude = forceMag;
            if (forceMag > _maxViscosityForceMagnitude)
                _maxViscosityForceMagnitude = forceMag;
            _totalViscosityForceMagnitude += forceMag;
            _viscosityForceCount++;

            // Add viscous acceleration to existing acceleration (pressure + gravity will be added later)
            pi.Acceleration += viscAccel;
            _particles[i] = pi;
        }
    }

    /// <summary>
    /// Finds all actual neighbors of particle at the given index.
    /// A neighbor is a particle within smoothing radius distance.
    /// Uses the spatial grid to find candidates, then verifies with squared distance.
    ///
    /// Results are written to the provided list (cleared first). Passing a reusable
    /// list avoids per-call allocation.
    /// </summary>
    public void GetNeighbors(int particleIndex, List<int> results)
    {
        results.Clear();

        Vector3 pos = _particles[particleIndex].Position;
        float radiusSq = _parameters.SmoothingRadius * _parameters.SmoothingRadius;

        List<int> candidates = _grid.QueryCandidates(pos);

        for (int i = 0; i < candidates.Count; i++)
        {
            int candidateIndex = candidates[i];

            // Skip self
            if (candidateIndex == particleIndex)
                continue;

            Vector3 diff = _particles[candidateIndex].Position - pos;
            float distSq = diff.LengthSquared();

            if (distSq <= radiusSq)
                results.Add(candidateIndex);
        }
    }

    /// <summary>
    /// Returns the half-extents of the container for rendering.
    /// X = width/2, Y = height, Z = depth/2. Bottom is at y = 0.
    /// </summary>
    public (float HalfX, float Height, float HalfZ) GetContainerHalfExtents()
    {
        return (_parameters.ContainerWidth * 0.5f, _parameters.ContainerHeight, _parameters.ContainerDepth * 0.5f);
    }

    /// <summary>
    /// Enforces static container boundaries after position integration.
    /// For each particle that has crossed a wall, the position is clamped back
    /// inside the container and the normal velocity component is reflected
    /// using the configured coefficient of restitution.
    ///
    /// Container bounds: x ∈ [-w/2, w/2], y ∈ [0, height] (open top), z ∈ [-d/2, d/2].
    /// </summary>
    private void HandleBoundaryCollisions()
    {
        float halfW = _parameters.ContainerWidth * 0.5f;
        float halfD = _parameters.ContainerDepth * 0.5f;
        float height = _parameters.ContainerHeight;
        float e = _parameters.BoundaryRestitution;
        int collisions = 0;

        // Only enforce geometric boundaries on fluid particles.
        // Boundary particles are placed inside the container and never move.
        int fluidCount = _particles.Count - _boundaryParticleCount;
        for (int i = 0; i < fluidCount; i++)
        {
            var p = _particles[i];
            bool hit = false;

            // Bottom wall (y = 0)
            if (p.Position.Y < 0.0f)
            {
                p.Position = new Vector3(p.Position.X, 0.0f, p.Position.Z);
                if (p.Velocity.Y < 0.0f)
                {
                    p.Velocity = new Vector3(p.Velocity.X, -p.Velocity.Y * e, p.Velocity.Z);
                    hit = true;
                }
            }

            // Left wall (x = -halfW)
            if (p.Position.X < -halfW)
            {
                p.Position = new Vector3(-halfW, p.Position.Y, p.Position.Z);
                if (p.Velocity.X < 0.0f)
                {
                    p.Velocity = new Vector3(-p.Velocity.X * e, p.Velocity.Y, p.Velocity.Z);
                    hit = true;
                }
            }

            // Right wall (x = +halfW)
            if (p.Position.X > halfW)
            {
                p.Position = new Vector3(halfW, p.Position.Y, p.Position.Z);
                if (p.Velocity.X > 0.0f)
                {
                    p.Velocity = new Vector3(-p.Velocity.X * e, p.Velocity.Y, p.Velocity.Z);
                    hit = true;
                }
            }

            // Front wall (z = -halfD)
            if (p.Position.Z < -halfD)
            {
                p.Position = new Vector3(p.Position.X, p.Position.Y, -halfD);
                if (p.Velocity.Z < 0.0f)
                {
                    p.Velocity = new Vector3(p.Velocity.X, p.Velocity.Y, -p.Velocity.Z * e);
                    hit = true;
                }
            }

            // Back wall (z = +halfD)
            if (p.Position.Z > halfD)
            {
                p.Position = new Vector3(p.Position.X, p.Position.Y, halfD);
                if (p.Velocity.Z > 0.0f)
                {
                    p.Velocity = new Vector3(p.Velocity.X, p.Velocity.Y, -p.Velocity.Z * e);
                    hit = true;
                }
            }

            if (hit)
                collisions++;

            _particles[i] = p;
        }

        _lastBoundaryCollisions = collisions;
    }

    /// <summary>
    /// Diagnostic method: compares grid-based neighbor search against brute-force O(N²).
    /// Returns a string with statistics and any mismatches. Intended for debugging only.
    /// Assumes the grid has been rebuilt (call after RebuildGrid / at the start of a step).
    /// </summary>
    public string RunNeighborSearchDiagnostic()
    {
        int n = _particles.Count;
        float radiusSq = _parameters.SmoothingRadius * _parameters.SmoothingRadius;
        var gridResults = new List<int>();
        int totalGridNeighbors = 0;
        int totalBruteNeighbors = 0;
        int mismatches = 0;
        int particlesWithZeroNeighbors = 0;

        for (int i = 0; i < n; i++)
        {
            // Grid-based search
            GetNeighbors(i, gridResults);
            var gridSet = new HashSet<int>(gridResults);
            totalGridNeighbors += gridResults.Count;

            // Brute-force search
            var bruteSet = new HashSet<int>();
            Vector3 posI = _particles[i].Position;
            for (int j = 0; j < n; j++)
            {
                if (i == j) continue;
                Vector3 diff = _particles[j].Position - posI;
                if (diff.LengthSquared() <= radiusSq)
                    bruteSet.Add(j);
            }
            totalBruteNeighbors += bruteSet.Count;

            // Compare
            if (!gridSet.SetEquals(bruteSet))
            {
                mismatches++;

                // Find which indices differ
                var missing = new List<int>();
                var extra = new List<int>();
                foreach (int idx in bruteSet)
                    if (!gridSet.Contains(idx)) missing.Add(idx);
                foreach (int idx in gridSet)
                    if (!bruteSet.Contains(idx)) extra.Add(idx);

                if (mismatches <= 5) // Print first few mismatches
                {
                    System.Diagnostics.Debug.WriteLine(
                        $"  Particle {i} at ({posI.X:F4}, {posI.Y:F4}, {posI.Z:F4}): " +
                        $"grid={gridResults.Count} brute={bruteSet.Count} " +
                        $"missing=[{string.Join(",", missing)}] extra=[{string.Join(",", extra)}]");
                }
            }

            if (gridResults.Count == 0)
                particlesWithZeroNeighbors++;
        }

        return $"Particles: {n}, " +
               $"Grid total neighbors: {totalGridNeighbors}, " +
               $"Brute total neighbors: {totalBruteNeighbors}, " +
               $"Mismatches: {mismatches}, " +
               $"Zero-neighbor particles: {particlesWithZeroNeighbors}, " +
               $"Grid cells: {_grid.CellCount}";
    }

    /// <summary>
    /// Diagnostic: compares grid-based density against brute-force O(N²) density.
    /// Returns a string with density statistics and any mismatches. Intended for debugging.
    /// </summary>
    public string RunDensityDiagnostic()
    {
        int n = _particles.Count;
        float h = _parameters.SmoothingRadius;
        float h2 = h * h;
        float mass = _parameters.ParticleMass;
        float poly6Coeff = SPHKernels.Poly6Coefficient(h);
        float poly6Origin = SPHKernels.Poly6AtOrigin(h);

        float minDensity = float.MaxValue;
        float maxDensity = float.MinValue;
        float totalDensity = 0.0f;
        int mismatches = 0;

        for (int i = 0; i < n; i++)
        {
            // Grid-based density (already computed in particle.Density)
            float gridDensity = _particles[i].Density;

            // Brute-force density
            float bruteDensity = mass * poly6Origin;
            Vector3 posI = _particles[i].Position;
            for (int j = 0; j < n; j++)
            {
                if (i == j) continue;
                Vector3 diff = _particles[j].Position - posI;
                float distSq = diff.LengthSquared();
                bruteDensity += mass * SPHKernels.Poly6(distSq, h2, poly6Coeff);
            }

            // Compare
            float densityDiff = MathF.Abs(gridDensity - bruteDensity);
            if (densityDiff > 0.001f) // tolerance for floating point
                mismatches++;

            // Stats
            if (gridDensity < minDensity) minDensity = gridDensity;
            if (gridDensity > maxDensity) maxDensity = gridDensity;
            totalDensity += gridDensity;
        }

        float avgDensity = n > 0 ? totalDensity / n : 0.0f;

        return $"Particles: {n}, " +
               $"Density — min: {minDensity:F4}, max: {maxDensity:F4}, avg: {avgDensity:F4}, " +
               $"Rest: {_parameters.RestDensity:F4}, " +
               $"Mismatches: {mismatches}/{n}, " +
               $"h: {h}, mass: {mass}";
    }

    /// <summary>
    /// Diagnostic: reports pressure statistics computed from the equation of state.
    /// Assumes ComputeAllDensities and ComputeAllPressures have been called.
    /// Verifies P = k * (rho - rho_0) by recalculating from stored density values.
    /// Reports only on fluid particles (boundary particles have zero pressure by design).
    /// </summary>
    public string RunPressureDiagnostic()
    {
        int n = _particles.Count;
        int fluidCount = n - _boundaryParticleCount;
        float k = _parameters.PressureStiffness;
        float rho0 = _parameters.RestDensity;

        float minPressure = float.MaxValue;
        float maxPressure = float.MinValue;
        float totalPressure = 0.0f;
        float minDensity = float.MaxValue;
        float maxDensity = float.MinValue;
        float totalDensity = 0.0f;
        int eosMismatches = 0;

        // Only report pressure stats for fluid particles
        for (int i = 0; i < fluidCount; i++)
        {
            float density = _particles[i].Density;
            float pressure = _particles[i].Pressure;

            // Verify equation of state independently
            float expectedPressure = k * (density - rho0);
            if (MathF.Abs(pressure - expectedPressure) > 0.001f)
                eosMismatches++;

            if (pressure < minPressure) minPressure = pressure;
            if (pressure > maxPressure) maxPressure = pressure;
            totalPressure += pressure;

            if (density < minDensity) minDensity = density;
            if (density > maxDensity) maxDensity = density;
            totalDensity += density;
        }

        float avgPressure = fluidCount > 0 ? totalPressure / fluidCount : 0.0f;
        float avgDensity = fluidCount > 0 ? totalDensity / fluidCount : 0.0f;

        return $"Particles: {n} (fluid: {fluidCount}, boundary: {_boundaryParticleCount})\n" +
               $"  Pressure — min: {minPressure:F4}, max: {maxPressure:F4}, avg: {avgPressure:F4}\n" +
               $"  Density  — min: {minDensity:F4}, max: {maxDensity:F4}, avg: {avgDensity:F4}\n" +
               $"  Rest density: {rho0:F4}, Pressure stiffness (k): {k:F4}\n" +
               $"  EOS mismatches: {eosMismatches}/{fluidCount} " +
               $"(P = k * (rho - rho_0))";
    }

    /// <summary>
    /// Diagnostic (F5): pressure-restoration test.
    /// After a clean lattice is spawned and one particle is displaced inward, this method
    /// runs the full SPH pipeline (grid, density, pressure, force) and reports detailed
    /// information about the perturbed particle and its neighbors.
    /// Validates that the pressure force opposes the compression.
    /// 
    /// Call after: RebuildGrid, ComputeAllDensities, ComputeAllPressures, ComputeAllPressureForces.
    /// </summary>
    /// <param name="perturbedIndex">Index of the particle that was displaced.</param>
    /// <param name="displacement">The displacement vector that was applied (original → perturbed).</param>
    public string RunPressureRestorationDiagnostic(int perturbedIndex, Vector3 displacement)
    {
        var sb = new System.Text.StringBuilder();

        sb.AppendLine("═══════════════════════════════════════════════════════════");
        sb.AppendLine("  PRESSURE RESTORATION TEST");
        sb.AppendLine("═══════════════════════════════════════════════════════════");

        // --- Perturbed particle ---
        var pp = _particles[perturbedIndex];
        float dispMag = displacement.Length();

        sb.AppendLine();
        sb.AppendLine($"Perturbed particle index: {perturbedIndex}");
        sb.AppendLine($"  Position:          ({pp.Position.X:F6}, {pp.Position.Y:F6}, {pp.Position.Z:F6})");
        sb.AppendLine($"  Displacement:      ({displacement.X:F6}, {displacement.Y:F6}, {displacement.Z:F6})  |d| = {dispMag:F6} m");
        sb.AppendLine($"  Density:           {pp.Density:F6}  (rest = {_parameters.RestDensity:F4})");
        sb.AppendLine($"  Pressure:          {pp.Pressure:F6}  (k={_parameters.PressureStiffness}, P = k*(ρ-ρ₀))");
        sb.AppendLine($"  Pressure force:    ({pp.Acceleration.X * pp.Mass:F6}, {pp.Acceleration.Y * pp.Mass:F6}, {pp.Acceleration.Z * pp.Mass:F6})  |F| = {(pp.Acceleration * pp.Mass).Length():F6} N");
        sb.AppendLine($"  Acceleration:      ({pp.Acceleration.X:F6}, {pp.Acceleration.Y:F6}, {pp.Acceleration.Z:F6})  |a| = {pp.Acceleration.Length():F6} m/s²");

        // Direction the particle was pushed (inward) vs force direction
        Vector3 forceOnPerturbed = pp.Acceleration * pp.Mass;
        if (dispMag > 1e-10f)
        {
            Vector3 dispDir = displacement / dispMag;
            float alignment = Vector3.Dot(forceOnPerturbed, dispDir);
            sb.AppendLine($"  Force·displacement: {alignment:F6}  (negative = restoring, positive = amplifying)");
            sb.AppendLine($"  Force opposes displacement: {(alignment < 0 ? "YES" : "NO")}");
        }

        // --- Neighbors ---
        GetNeighbors(perturbedIndex, _neighborBuffer);
        sb.AppendLine();
        sb.AppendLine($"Immediate SPH neighbors: {_neighborBuffer.Count}");
        sb.AppendLine($"  h = {_parameters.SmoothingRadius:F4} m,  h² = {_parameters.SmoothingRadius * _parameters.SmoothingRadius:F6}");

        float totalNeighborDensity = 0.0f;
        float totalNeighborPressure = 0.0f;
        float minNeighborDist = float.MaxValue;
        float maxNeighborDist = float.MinValue;

        for (int n = 0; n < _neighborBuffer.Count; n++)
        {
            int j = _neighborBuffer[n];
            var pj = _particles[j];
            Vector3 diff = pj.Position - pp.Position;
            float dist = diff.Length();
            totalNeighborDensity += pj.Density;
            totalNeighborPressure += pj.Pressure;
            if (dist < minNeighborDist) minNeighborDist = dist;
            if (dist > maxNeighborDist) maxNeighborDist = dist;

            Vector3 forceJ = pj.Acceleration * pj.Mass;
            sb.AppendLine($"  [{j,3}] pos=({pj.Position.X:F6}, {pj.Position.Y:F6}, {pj.Position.Z:F6})" +
                         $"  dist={dist:F6}" +
                         $"  ρ={pj.Density:F4}  P={pj.Pressure:F4}" +
                         $"  F=({forceJ.X:E4}, {forceJ.Y:E4}, {forceJ.Z:E4})");
        }

        if (_neighborBuffer.Count > 0)
        {
            sb.AppendLine($"  Neighbor distance range: [{minNeighborDist:F6}, {maxNeighborDist:F6}]");
            sb.AppendLine($"  Avg neighbor density:   {totalNeighborDensity / _neighborBuffer.Count:F4}");
            sb.AppendLine($"  Avg neighbor pressure:  {totalNeighborPressure / _neighborBuffer.Count:F4}");
        }

        // --- Symmetry check: compare with a mirror particle if it exists ---
        // Find the particle at the mirror position relative to center (approx)
        // For a clean lattice, the center of the block is at (0, 0.5, 0).
        // We look for a particle near the displaced particle's original position.
        sb.AppendLine();
        sb.AppendLine("--- Symmetry check (mirror particle) ---");

        Vector3 originalPos = pp.Position - displacement;
        // The "opposite" particle from center would be at originalPos reflected through center
        Vector3 center = new Vector3(0.0f, 0.5f, 0.0f);
        Vector3 mirrorTarget = 2.0f * center - originalPos;

        // Search for closest particle to mirrorTarget
        int mirrorIndex = -1;
        float mirrorDistSq = float.MaxValue;
        for (int i = 0; i < _particles.Count; i++)
        {
            if (i == perturbedIndex) continue;
            float dSq = (_particles[i].Position - mirrorTarget).LengthSquared();
            if (dSq < mirrorDistSq)
            {
                mirrorDistSq = dSq;
                mirrorIndex = i;
            }
        }

        if (mirrorIndex >= 0 && mirrorDistSq < 0.01f * 0.01f)
        {
            var pm = _particles[mirrorIndex];
            Vector3 mirrorDisp = -displacement; // mirror particle moved in opposite direction
            Vector3 mirrorForceVec = pm.Acceleration * pm.Mass;
            float mirrorForceAlign = 0.0f;
            if (mirrorDisp.Length() > 1e-10f)
            {
                mirrorForceAlign = Vector3.Dot(mirrorForceVec, mirrorDisp / mirrorDisp.Length());
            }
            sb.AppendLine($"  Mirror particle [{mirrorIndex}] at ({pm.Position.X:F6}, {pm.Position.Y:F6}, {pm.Position.Z:F6})");
            sb.AppendLine($"    Density: {pm.Density:F6},  Pressure: {pm.Pressure:F6}");
            sb.AppendLine($"    Force: ({mirrorForceVec.X:E4}, {mirrorForceVec.Y:E4}, {mirrorForceVec.Z:E4})  |F| = {mirrorForceVec.Length():F6}");
            sb.AppendLine($"    Force·mirror-displacement: {mirrorForceAlign:F6}  (should be restoring)");
        }
        else
        {
            sb.AppendLine($"  No mirror particle found within tolerance (closest dist²={mirrorDistSq:F6})");
        }

        // --- Summary ---
        sb.AppendLine();
        sb.AppendLine("--- Summary ---");
        sb.AppendLine($"  Particles: {_particles.Count}");
        sb.AppendLine($"  Grid cells: {_grid.CellCount}");
        if (dispMag > 1e-10f)
        {
            Vector3 forceDir = forceOnPerturbed.Length() > 1e-10f ? forceOnPerturbed / forceOnPerturbed.Length() : Vector3.Zero;
            Vector3 dispDir = displacement / dispMag;
            float cosAngle = Vector3.Dot(forceDir, dispDir);
            sb.AppendLine($"  Force magnitude:   {forceOnPerturbed.Length():E6} N");
            sb.AppendLine($"  Displacement magnitude: {dispMag:F6} m");
            sb.AppendLine($"  cos(angle) between force and displacement: {cosAngle:F6}");
            sb.AppendLine($"  Interpretation: force {(cosAngle < 0 ? "RESTORES" : "AMPLIFIES")} the perturbation");
        }

        sb.AppendLine("═══════════════════════════════════════════════════════════");
        return sb.ToString();
    }

    /// <summary>
    /// Diagnostic (F4): reports pressure and viscosity force statistics after
    /// ComputeAllPressureForces and ComputeAllViscosityForces have been called.
    /// Includes min/max/average force magnitudes, max acceleration, max velocity,
    /// and non-finite value detection.
    /// </summary>
    public string RunPressureForceDiagnostic()
    {
        int n = _particles.Count;
        int fluidCount = n - _boundaryParticleCount;
        float avgPressureForce = _pressureForceCount > 0
            ? _totalPressureForceMagnitude / _pressureForceCount
            : 0.0f;
        float avgViscosityForce = _viscosityForceCount > 0
            ? _totalViscosityForceMagnitude / _viscosityForceCount
            : 0.0f;

        // Count zero-density fluid particles (boundary particles don't need density for diagnostics)
        int zeroDensityCount = 0;
        for (int i = 0; i < fluidCount; i++)
        {
            if (_particles[i].Density <= 0.0f)
                zeroDensityCount++;
        }

        // Verify that acceleration was actually set by pressure force
        // (check that at least one fluid particle has non-zero acceleration from pressure)
        int particlesWithPressureAccel = 0;
        for (int i = 0; i < fluidCount; i++)
        {
            if (_particles[i].Acceleration.LengthSquared() > 1e-10f)
                particlesWithPressureAccel++;
        }

        return $"Particles: {n} (fluid: {fluidCount}, boundary: {_boundaryParticleCount})\n" +
               $"  Pressure force magnitude — min: {_minPressureForceMagnitude:E4}, " +
               $"max: {_maxPressureForceMagnitude:E4}, avg: {avgPressureForce:E4}\n" +
               $"  Viscosity force magnitude — min: {_minViscosityForceMagnitude:E4}, " +
               $"max: {_maxViscosityForceMagnitude:E4}, avg: {avgViscosityForce:E4}\n" +
               $"  Max acceleration (pressure+visc): {_maxAccelerationMagnitude:E4} m/s²\n" +
               $"  Max velocity: {_maxVelocityMagnitude:E4} m/s\n" +
               $"  Zero-density particles: {zeroDensityCount}/{fluidCount}\n" +
               $"  Particles with non-zero acceleration: {particlesWithPressureAccel}/{fluidCount}\n" +
               $"  Non-finite values encountered: {_pressureForceNonFinite}\n" +
               $"  Parameters: h={_parameters.SmoothingRadius}, mass={_parameters.ParticleMass}, " +
               $"k={_parameters.PressureStiffness}, restDensity={_parameters.RestDensity}, " +
               $"ν={_parameters.KinematicViscosity:E2}";
    }
}
