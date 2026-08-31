# Physics Simulator

A real-time material simulator built with Godot 4.6 and C#, beginning with SPH (Smoothed Particle Hydrodynamics) for water.

## Current Status

Initial project scaffold. 64 test particles rendered via MultiMesh. No SPH physics yet.

## Architecture

```
Simulation/           # Engine-independent simulation logic
  Particle.cs         #   Particle data struct (System.Numerics)
  SimulationParameters.cs # Tunable constants
  FluidSimulation.cs  #   Core simulation (state + step)

Godot/                # Godot-specific bridge
  SimulationNode.cs   #   Node3D that owns the simulation and drives updates

Rendering/            # Visualization
  SimulationRenderer.cs #  MultiMeshInstance3D-based particle renderer
```

**Key design decisions:**

- **Simulation layer uses `System.Numerics` types**, not Godot types. This keeps the simulation portable and extractable.
- **MultiMeshInstance3D** for rendering: one draw call for all particles, not one Node3D per particle. Scales to large particle counts.
- **Fixed-timestep accumulator** in the simulation for deterministic behavior regardless of frame rate.
- **No Godot dependency** in the Simulation/ namespace. Godot is only referenced in Godot/ and Rendering/.

## How to Run

### Prerequisites

1. **.NET SDK 8.0**: `sudo dnf install -y dotnet-sdk-8.0`
2. **Godot 4.6.3 .NET build** - You MUST use the .NET build, not the standard build. The standard Fedora `godot` package does NOT include C# support and will produce `No loader found for resource` errors.

   Download from: https://godotengine.org/download/archive/4.6.3-stable/
   Choose: **Linux - .NET - x86_64** (`mono_linux_x86_64.zip`)

   Extract and run the `Godot_v4.6.3-stable_mono_linux.x86_64` binary.

### Running

```bash
# Build C# project
dotnet build

# Run (use the .NET build, NOT the standard 'godot' command)
/path/to/Godot_v4.6.3-stable_mono_linux.x86_64 --path .
```

A shell alias `godot` -> `godot-dotnet` has been added to your `~/.bashrc` for convenience.

**Controls:**
- `Space` - Toggle simulation on/off
- `R` - Reset simulation

## What Has Been Implemented

- Project structure and C# configuration
- Particle data struct
- Simulation parameters (timestep, density, viscosity, gravity, etc.)
- FluidSimulation class with fixed-timestep accumulator and gravity placeholder
- SimulationNode (Godot bridge) with play/pause/reset
- SimulationRenderer using MultiMeshInstance3D
- 4x4x4 test cube of particles
- Scene with camera, lighting, and sky

## What Has NOT Been Implemented (Intentionally)

- SPH density calculations
- SPH pressure calculations
- Viscosity forces
- Surface tension
- Neighbor search (spatial hashing, grid, etc.)
- Fluid dynamics (kernel functions, force accumulation)
- Boundary conditions / collision
- GPU compute shaders
- Surface reconstruction / meshing
- Advanced rendering (transparency, reflections, refractions)
- UI for parameter tweaking
- Profiling / performance tools

## Planned Development Stages

1. **Neighbor search** - Spatial hashing or uniform grid
2. **SPH kernel functions** - Poly6, Spiky, Viscosity kernels
3. **Density & pressure** - Compute particle densities and pressure from equation of state
4. **Force accumulation** - Pressure forces, viscosity, gravity integration
5. **Boundary conditions** - Contain fluid in a volume
6. **Surface tension** - Color field estimation and surface tension forces
7. **Rendering improvements** - Better particle appearance, potential surface meshing
8. **Performance** - Spatial hashing optimization, eventual GPU acceleration
