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
    /// One fixed-size simulation step. This is where SPH forces will be computed.
    /// Currently a placeholder: integrates gravity and applies basic Euler integration.
    /// </summary>
    private void FixedStep(float dt)
    {
        RebuildGrid();
        ComputeAllDensities();

        var gravity = new Vector3(0.0f, _parameters.Gravity, 0.0f);

        for (int i = 0; i < _particles.Count; i++)
        {
            var p = _particles[i];

            // Placeholder: apply gravity
            p.Acceleration = gravity;

            // Euler integration
            p.Velocity += p.Acceleration * dt;
            p.Position += p.Velocity * dt;

            _particles[i] = p;
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
}
