using Godot;
using PhysicsSimulator.Rendering;
using PhysicsSimulator.Simulation;

namespace PhysicsSimulator.Godot;

/// <summary>
/// Main simulation node. Bridges the simulation layer to the Godot scene tree.
/// Owns the FluidSimulation and SimulationRenderer, drives the simulation each frame.
/// </summary>
public partial class SimulationNode : Node3D
{
    private FluidSimulation _simulation = null!;
    private SimulationRenderer _renderer = null!;

    /// <summary>
    /// Whether the simulation is currently running.
    /// </summary>
    public bool IsRunning { get; private set; }

    public override void _Ready()
    {
        var parameters = new SimulationParameters();
        _simulation = new FluidSimulation(parameters);

        _renderer = new SimulationRenderer();
        AddChild(_renderer);

        // Create a small test block of particles so we can verify rendering
        SpawnTestParticles();

        // Start paused so the user can see the initial state
        IsRunning = false;
        GD.Print($"[SimulationNode] Ready. {_simulation.ParticleCount} particles. Simulation paused.");
    }

    public override void _Process(double delta)
    {
        if (!IsRunning)
            return;

        _simulation.Step((float)delta);
        _renderer.UpdateParticles(_simulation.Particles);
    }

    public override void _UnhandledInput(InputEvent @event)
    {
        if (@event is InputEventKey key && key.Pressed)
        {
            switch (key.Keycode)
            {
                case Key.Space:
                    ToggleSimulation();
                    GetViewport().SetInputAsHandled();
                    break;
                case Key.R:
                    ResetSimulation();
                    GetViewport().SetInputAsHandled();
                    break;
                case Key.F1:
                    RunNeighborDiagnostic();
                    GetViewport().SetInputAsHandled();
                    break;
                case Key.F2:
                    RunDensityDiagnostic();
                    GetViewport().SetInputAsHandled();
                    break;
                case Key.F3:
                    RunPressureDiagnostic();
                    GetViewport().SetInputAsHandled();
                    break;
                case Key.F4:
                    RunPressureForceDiagnostic();
                    GetViewport().SetInputAsHandled();
                    break;
                case Key.N:
                    StepOnce();
                    GetViewport().SetInputAsHandled();
                    break;
            }
        }
    }

    public void ToggleSimulation()
    {
        IsRunning = !IsRunning;
        GD.Print($"[SimulationNode] Simulation {(IsRunning ? "started" : "paused")}.");
    }

    public void StepOnce()
    {
        _simulation.StepOnce();
        _renderer.UpdateParticles(_simulation.Particles);
        GD.Print($"[Simulation] Step {_simulation.StepCount} (single step)");
    }

    public void ResetSimulation()
    {
        IsRunning = false;
        _simulation.Reset();
        SpawnTestParticles();
        _renderer.UpdateParticles(_simulation.Particles);

        if (GetNode<Camera3D>("../Camera3D") is CameraController cam)
            cam.ResetCamera();

        GD.Print($"[SimulationNode] Simulation reset. {_simulation.ParticleCount} particles.");
    }

    /// <summary>
    /// Creates a small 4x4x4 cube of particles for testing.
    /// </summary>
    private void SpawnTestParticles()
    {
        float spacing = 0.08f;
        int countPerAxis = 4;
        float offset = (countPerAxis - 1) * spacing * 0.5f;
        var rng = new System.Random();

        for (int x = 0; x < countPerAxis; x++)
        for (int y = 0; y < countPerAxis; y++)
        for (int z = 0; z < countPerAxis; z++)
        {
            float noise = 0.005f;
            var pos = new System.Numerics.Vector3(
                x * spacing - offset + (float)(rng.NextDouble() * 2 - 1) * noise,
                y * spacing + 0.5f + (float)(rng.NextDouble() * 2 - 1) * noise,
                z * spacing - offset + (float)(rng.NextDouble() * 2 - 1) * noise
            );
            _simulation.AddParticle(pos, System.Numerics.Vector3.Zero, _simulation.Parameters.ParticleMass);
        }

        _renderer.UpdateParticles(_simulation.Particles);
    }

    private void RunNeighborDiagnostic()
    {
        // Rebuild grid at current positions so we can query
        _simulation.Grid.Clear();
        for (int i = 0; i < _simulation.ParticleCount; i++)
            _simulation.Grid.Insert(i, _simulation.Particles[i].Position);

        string result = _simulation.RunNeighborSearchDiagnostic();
        GD.Print($"[NeighborDiagnostic] {result}");
    }

    private void RunDensityDiagnostic()
    {
        // Rebuild grid and compute density at current positions
        _simulation.Grid.Clear();
        for (int i = 0; i < _simulation.ParticleCount; i++)
            _simulation.Grid.Insert(i, _simulation.Particles[i].Position);

        _simulation.ComputeAllDensities();

        string result = _simulation.RunDensityDiagnostic();
        GD.Print($"[DensityDiagnostic] {result}");
    }

    private void RunPressureDiagnostic()
    {
        _simulation.Grid.Clear();
        for (int i = 0; i < _simulation.ParticleCount; i++)
            _simulation.Grid.Insert(i, _simulation.Particles[i].Position);

        _simulation.ComputeAllDensities();
        _simulation.ComputeAllPressures();

        string result = _simulation.RunPressureDiagnostic();
        GD.Print($"[PressureDiagnostic]\n{result}");
    }

    private void RunPressureForceDiagnostic()
    {
        _simulation.Grid.Clear();
        for (int i = 0; i < _simulation.ParticleCount; i++)
            _simulation.Grid.Insert(i, _simulation.Particles[i].Position);

        _simulation.ComputeAllDensities();
        _simulation.ComputeAllPressures();
        _simulation.ComputeAllPressureForces();

        string result = _simulation.RunPressureForceDiagnostic();
        GD.Print($"[PressureForceDiagnostic]\n{result}");
    }
}
