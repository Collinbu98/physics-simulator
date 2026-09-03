using System;
using System.Diagnostics;
using System.Numerics;
using PhysicsSimulator.Simulation;

namespace PhysicsSimulator.Benchmark;

class Program
{
    static void Main(string[] args)
    {
        float skinWidth = 0.02f;
        if (args.Length > 0 && float.TryParse(args[0], out float parsed))
            skinWidth = parsed;

        var parameters = new SimulationParameters
        {
            TimeStep = 0.0005f,
            RestDensity = 1000.0f,
            PressureStiffness = 50.0f,
            KinematicViscosity = 1e-4f,
            ParticleMass = 0.27f,
            SmoothingRadius = 0.075f,
            Gravity = -9.81f,
            TimeScale = 1.0f,
            ContainerWidth = 1.0f,
            ContainerHeight = 1.0f,
            ContainerDepth = 1.0f,
            BoundaryRestitution = 0.3f,
            BoundaryParticleSpacing = 0.06f,
        };

        var sim = new FluidSimulation(parameters, skinWidth);

        float fluidSpacing = 0.05f;
        float halfBlock = 4.5f * fluidSpacing;
        for (int x = 0; x < 10; x++)
        for (int y = 0; y < 10; y++)
        for (int z = 0; z < 10; z++)
        {
            float px = -halfBlock + x * fluidSpacing;
            float py = fluidSpacing + y * fluidSpacing;
            float pz = -halfBlock + z * fluidSpacing;
            sim.AddParticle(new Vector3(px, py, pz), Vector3.Zero, parameters.ParticleMass);
        }
        sim.GenerateBoundaryParticles();

        Console.WriteLine($"Skin={skinWidth:F3} m  Particles: {sim.FluidParticleCount}f + {sim.BoundaryParticleCount}b");

        const int steps = 20000;

        var totalSw = Stopwatch.StartNew();
        for (int i = 0; i < steps; i++)
        {
            sim.StepOnce();
            if (sim.LastProfileReport is { } report)
            {
                Console.WriteLine(report);
                sim.ClearProfileReport();
            }
        }
        totalSw.Stop();

        double totalTimeMs = totalSw.Elapsed.TotalMilliseconds;
        double stepsPerSec = steps / (totalTimeMs / 1000.0);

        // Correctness
        int fluidCount = sim.FluidParticleCount;
        int bndCount = sim.BoundaryParticleCount;
        float minDensity = float.MaxValue, maxDensity = float.MinValue, totalDensity = 0;
        float minPressure = float.MaxValue, maxPressure = float.MinValue, totalPressure = 0;
        int zeroDensityCount = 0, nonFiniteCount = 0;
        float maxVelocity = 0;
        for (int i = 0; i < fluidCount; i++)
        {
            var p = sim.Particles[i];
            if (p.Density < minDensity) minDensity = p.Density;
            if (p.Density > maxDensity) maxDensity = p.Density;
            totalDensity += p.Density;
            if (p.Density <= 0) zeroDensityCount++;
            if (p.Pressure < minPressure) minPressure = p.Pressure;
            if (p.Pressure > maxPressure) maxPressure = p.Pressure;
            totalPressure += p.Pressure;
            float v = p.Velocity.Length();
            if (!float.IsFinite(p.Density) || !float.IsFinite(p.Pressure) || !float.IsFinite(v))
                nonFiniteCount++;
            if (v > maxVelocity) maxVelocity = v;
        }
        float avgDensity = fluidCount > 0 ? totalDensity / fluidCount : 0;
        float avgPressure = fluidCount > 0 ? totalPressure / fluidCount : 0;

        // Output CSV-friendly line
        Console.WriteLine($"  Total: {totalTimeMs:F1} ms | Steps/sec: {stepsPerSec:F0}");
        Console.WriteLine($"  Density: min={minDensity:F2} max={maxDensity:F2} avg={avgDensity:F2}");
        Console.WriteLine($"  Pressure: min={minPressure:F2} max={maxPressure:F2} avg={avgPressure:F2}");
        Console.WriteLine($"  Max velocity: {maxVelocity:F4} | Zero-density: {zeroDensityCount}/{fluidCount} | Non-finite: {nonFiniteCount}/{fluidCount}");
    }
}
