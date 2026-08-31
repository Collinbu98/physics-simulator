using System;
using System.Collections.Generic;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// Core simulation class. Owns particle state and advances it each step.
/// This class has no dependency on Godot or any other engine.
/// </summary>
public class FluidSimulation
{
    private readonly List<Particle> _particles = new();
    private readonly SimulationParameters _parameters;
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

    public FluidSimulation(SimulationParameters parameters)
    {
        _parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));
    }

    /// <summary>
    /// Resets the simulation to an empty state.
    /// </summary>
    public void Reset()
    {
        _particles.Clear();
        _timeAccumulator = 0.0f;
    }

    /// <summary>
    /// Adds a single particle with the given state.
    /// </summary>
    public void AddParticle(System.Numerics.Vector3 position, System.Numerics.Vector3 velocity, float mass)
    {
        _particles.Add(new Particle
        {
            Position = position,
            Velocity = velocity,
            Acceleration = System.Numerics.Vector3.Zero,
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
        var gravity = new System.Numerics.Vector3(0.0f, _parameters.Gravity, 0.0f);

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
}
