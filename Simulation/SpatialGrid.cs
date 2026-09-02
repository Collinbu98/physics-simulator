using System;
using System.Collections.Generic;
using System.Numerics;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// Uniform spatial grid for fast neighbor queries.
/// Maps particle positions to integer 3D cell coordinates.
/// Cell size is set to the SPH smoothing radius so that any neighbor
/// (distance <= smoothing radius) must reside in the same cell or one
/// of the 26 adjacent cells.
///
/// This structure is rebuilt every simulation step. It does NOT decide
/// whether two particles are neighbors — it only produces candidates.
/// The caller must verify actual distance.
/// </summary>
public class SpatialGrid
{
	private readonly Dictionary<(int X, int Y, int Z), List<int>> _cells = new();
	private readonly float _cellSize;
	private readonly float _inverseCellSize;

	// Reusable buffer for query results to avoid per-query allocation.
	private readonly List<int> _queryResults = new();

	public SpatialGrid(float cellSize)
	{
		if (cellSize <= 0f)
			throw new ArgumentOutOfRangeException(nameof(cellSize), "Cell size must be positive.");

		_cellSize = cellSize;
		_inverseCellSize = 1.0f / cellSize;
	}

	/// <summary>
	/// The cell size used by this grid (typically equal to the smoothing radius).
	/// </summary>
	public float CellSize => _cellSize;

	/// <summary>
	/// Clears all cells. Call this at the start of each rebuild.
	/// </summary>
	public void Clear()
	{
		// Clear each list individually to reuse the List<int> allocations.
		foreach (var kvp in _cells)
			kvp.Value.Clear();
	}

	/// <summary>
	/// Inserts a particle index at the cell corresponding to its position.
	/// </summary>
	public void Insert(int particleIndex, Vector3 position)
	{
		var cell = PositionToCell(position);

		if (!_cells.TryGetValue(cell, out var list))
		{
			list = new List<int>();
			_cells[cell] = list;
		}

		list.Add(particleIndex);
	}

	/// <summary>
	/// Converts a world-space position to integer cell coordinates.
	/// Uses MathF.Floor so that negative positions map to negative cell indices
	/// correctly (e.g., position -0.05 with cellSize 0.1 -> cell -1).
	/// </summary>
	public (int X, int Y, int Z) PositionToCell(Vector3 position)
	{
		return (
			(int)MathF.Floor(position.X * _inverseCellSize),
			(int)MathF.Floor(position.Y * _inverseCellSize),
			(int)MathF.Floor(position.Z * _inverseCellSize)
		);
	}

	/// <summary>
	/// Queries candidate neighbor indices for a particle at the given position.
	/// Returns indices of all particles in the 27 cells (own cell + 26 neighbors).
	/// The caller MUST verify actual squared distance <= smoothingRadius².
	///
	/// The returned list is a reusable internal buffer — copy it if you need to
	/// hold onto the results beyond the next query call.
	/// </summary>
	public List<int> QueryCandidates(Vector3 position)
	{
		_queryResults.Clear();

		var (cx, cy, cz) = PositionToCell(position);

		for (int dx = -1; dx <= 1; dx++)
		{
			for (int dy = -1; dy <= 1; dy++)
			{
				for (int dz = -1; dz <= 1; dz++)
				{
					var neighbor = (cx + dx, cy + dy, cz + dz);

					if (_cells.TryGetValue(neighbor, out var list))
					{
						for (int i = 0; i < list.Count; i++)
							_queryResults.Add(list[i]);
					}
				}
			}
		}

		return _queryResults;
	}

	/// <summary>
	/// Returns the total number of cell entries (non-empty cells) in the grid.
	/// Useful for debugging.
	/// </summary>
	public int CellCount => _cells.Count;
}
