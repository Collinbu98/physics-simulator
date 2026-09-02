using System;
using System.Collections.Generic;
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

    // Diagnostic storage for pressure force magnitudes (updated each step)
    private float _minPressureForceMagnitude = float.MaxValue;
    private float _maxPressureForceMagnitude = float.MinValue;
    private float _totalPressureForceMagnitude;
    private int _pressureForceCount;
    private bool _pressureForceNonFinite;
    private float _maxAccelerationMagnitude;
    private float _maxVelocityMagnitude;
    private long _stepCount;

    /// <summary>
    /// Read-only access to particle data.
    /// </summary>
    public IReadOnlyList<Particle> Particles => _particles;

    /// <summary>
    /// The current number of particles in the simulation.
    /// </summary>
    public int ParticleCount => _particles.Count;

    public SimulationParameters Parameters => _parameters;

    /// <summary>
    /// Total number of fixed timesteps executed since the last Reset().
    /// </summary>
    public long StepCount => _stepCount;

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
        _grid.Clear();
        _timeAccumulator = 0.0f;
        _stepCount = 0;
    }

    /// <summary>
    /// Adds a single particle with the given state.
    /// </summary>
    public void AddParticle(Vector3 position, Vector3 velocity, float mass)
    {
        _particles.Add(new Particle
        {
            Position = position,
            Velocity = velocity,
            Acceleration = Vector3.Zero,
            Density = 0.0f,
            Pressure = 0.0f,
            Mass = mass,
        });
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
        RebuildGrid();
        ComputeAllDensities();
        ComputeAllPressures();
        ComputeAllPressureForces();

        var gravity = new Vector3(0.0f, _parameters.Gravity, 0.0f);

        // Reset velocity diagnostic before integration
        _maxVelocityMagnitude = 0.0f;

        for (int i = 0; i < _particles.Count; i++)
        {
            var p = _particles[i];

            // Combine gravity and pressure acceleration (already stored in p.Acceleration)
            // p.Acceleration was set by ComputeAllPressureForces; add gravity here
            p.Acceleration += gravity;

            // Euler integration
            p.Velocity += p.Acceleration * dt;
            p.Position += p.Velocity * dt;

            _particles[i] = p;

            // Track max velocity magnitude for diagnostics
            float velMag = p.Velocity.Length();
            if (velMag > _maxVelocityMagnitude)
                _maxVelocityMagnitude = velMag;
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
    /// Computes density for all particles using Poly6 kernel.
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
    /// Computes pressure for all particles using the equation of state: P = k * (rho - rho_0).
    /// Must be called after ComputeAllDensities so that Particle.Density is populated.
    /// Does not affect velocity, acceleration, or position.
    /// </summary>
    public void ComputeAllPressures()
    {
        float k = _parameters.PressureStiffness;
        float rho0 = _parameters.RestDensity;

        for (int i = 0; i < _particles.Count; i++)
        {
            var p = _particles[i];
            p.Pressure = k * (p.Density - rho0);
            _particles[i] = p;
        }
    }

    /// <summary>
    /// Computes SPH pressure forces for all particles using the symmetric formulation:
    ///   F_i = Σ_j m_j × (P_i/ρ_i² + P_j/ρ_j²) × ∇W_spiky(r_ij, h)
    ///
    /// Uses the Spiky kernel gradient. Handles zero-density and coincident particles safely.
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

        for (int i = 0; i < _particles.Count; i++)
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
    /// </summary>
    public string RunPressureDiagnostic()
    {
        int n = _particles.Count;
        float k = _parameters.PressureStiffness;
        float rho0 = _parameters.RestDensity;

        float minPressure = float.MaxValue;
        float maxPressure = float.MinValue;
        float totalPressure = 0.0f;
        float minDensity = float.MaxValue;
        float maxDensity = float.MinValue;
        float totalDensity = 0.0f;
        int eosMismatches = 0;

        for (int i = 0; i < n; i++)
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

        float avgPressure = n > 0 ? totalPressure / n : 0.0f;
        float avgDensity = n > 0 ? totalDensity / n : 0.0f;

        return $"Particles: {n}\n" +
               $"  Pressure — min: {minPressure:F4}, max: {maxPressure:F4}, avg: {avgPressure:F4}\n" +
               $"  Density  — min: {minDensity:F4}, max: {maxDensity:F4}, avg: {avgDensity:F4}\n" +
               $"  Rest density: {rho0:F4}, Pressure stiffness (k): {k:F4}\n" +
               $"  EOS mismatches: {eosMismatches}/{n} " +
               $"(P = k * (rho - rho_0))";
    }

    /// <summary>
    /// Diagnostic (F4): reports pressure force statistics after ComputeAllPressureForces has been called.
    /// Includes min/max/average pressure force magnitude, max acceleration, max velocity,
    /// and non-finite value detection.
    /// </summary>
    public string RunPressureForceDiagnostic()
    {
        int n = _particles.Count;
        float avgForce = _pressureForceCount > 0
            ? _totalPressureForceMagnitude / _pressureForceCount
            : 0.0f;

        // Count zero-density particles
        int zeroDensityCount = 0;
        for (int i = 0; i < n; i++)
        {
            if (_particles[i].Density <= 0.0f)
                zeroDensityCount++;
        }

        // Verify that acceleration was actually set by pressure force
        // (check that at least one particle has non-zero acceleration from pressure)
        int particlesWithPressureAccel = 0;
        for (int i = 0; i < n; i++)
        {
            if (_particles[i].Acceleration.LengthSquared() > 1e-10f)
                particlesWithPressureAccel++;
        }

        return $"Particles: {n}\n" +
               $"  Pressure force magnitude — min: {_minPressureForceMagnitude:E4}, " +
               $"max: {_maxPressureForceMagnitude:E4}, avg: {avgForce:E4}\n" +
               $"  Max acceleration (pressure only): {_maxAccelerationMagnitude:E4} m/s²\n" +
               $"  Max velocity: {_maxVelocityMagnitude:E4} m/s\n" +
               $"  Zero-density particles: {zeroDensityCount}/{n}\n" +
               $"  Particles with non-zero acceleration: {particlesWithPressureAccel}/{n}\n" +
               $"  Non-finite values encountered: {_pressureForceNonFinite}\n" +
               $"  Parameters: h={_parameters.SmoothingRadius}, mass={_parameters.ParticleMass}, " +
               $"k={_parameters.PressureStiffness}, restDensity={_parameters.RestDensity}";
    }
}
