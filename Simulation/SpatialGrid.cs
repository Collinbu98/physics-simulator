using System;
using System.Collections.Generic;
using System.Numerics;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// Uniform spatial grid for fast neighbor queries using a flat array.
/// Cell coordinates are shifted to non-negative indices for direct array access.
/// Cell size equals the SPH smoothing radius.
/// </summary>
public class SpatialGrid
{
	private readonly float _cellSize;
	private readonly float _inverseCellSize;

	// Flat 3D grid: index = (cx + Offset) + (cy + Offset) * Dim + (cz + Offset) * Dim * Dim
	// Offset must be large enough so that all possible cell coordinates + Offset >= 0.
	// With h=0.075 in a 1m container, cell coords span roughly -11..+13.
	// Offset=32 gives Dim=64 → 64³ = 262144 cells ≈ 2 MB (fits in L2 cache).
	private const int Offset = 32;
	private const int Dim = Offset * 2;          // 64
	private const int StrideY = Dim;             // 64
	private const int StrideZ = Dim * Dim;       // 4096
	private const int TotalCells = Dim * Dim * Dim; // 262144

	private readonly List<int>?[] _cells = new List<int>?[TotalCells];

	// Tracks which cells have been written to, so Clear() only touches those.
	// Persists across steps — never cleared, only grows as new cells are used.
	private readonly List<int> _occupiedCells = new();

	public SpatialGrid(float cellSize)
	{
		if (cellSize <= 0f)
			throw new ArgumentOutOfRangeException(nameof(cellSize), "Cell size must be positive.");
		_cellSize = cellSize;
		_inverseCellSize = 1.0f / cellSize;
	}

	public float CellSize => _cellSize;

	/// <summary>
	/// Clears all occupied cell lists for reuse. Does NOT clear the occupancy
	/// tracker — cells persist across steps so their List objects can be reused.
	/// </summary>
	public void Clear()
	{
		for (int i = 0; i < _occupiedCells.Count; i++)
			_cells[_occupiedCells[i]]!.Clear();
	}

	public void Insert(int particleIndex, Vector3 position)
	{
		int cx = (int)MathF.Floor(position.X * _inverseCellSize) + Offset;
		int cy = (int)MathF.Floor(position.Y * _inverseCellSize) + Offset;
		int cz = (int)MathF.Floor(position.Z * _inverseCellSize) + Offset;

		if ((uint)cx >= Dim || (uint)cy >= Dim || (uint)cz >= Dim)
			return;

		int idx = cx + cy * StrideY + cz * StrideZ;
		var list = _cells[idx];
		if (list == null)
		{
			list = new List<int>(8);
			_cells[idx] = list;
			_occupiedCells.Add(idx);
		}
		list.Add(particleIndex);
	}

	public void QueryCandidates(Vector3 position, List<int> results)
	{
		results.Clear();

		int cx = (int)MathF.Floor(position.X * _inverseCellSize) + Offset;
		int cy = (int)MathF.Floor(position.Y * _inverseCellSize) + Offset;
		int cz = (int)MathF.Floor(position.Z * _inverseCellSize) + Offset;

		for (int dx = -1; dx <= 1; dx++)
		{
			int nx = cx + dx;
			if ((uint)nx >= Dim) continue;

			for (int dy = -1; dy <= 1; dy++)
			{
				int ny = cy + dy;
				if ((uint)ny >= Dim) continue;

				int baseIdx = nx + ny * StrideY;

				for (int dz = -1; dz <= 1; dz++)
				{
					int nz = cz + dz;
					if ((uint)nz >= Dim) continue;

					var list = _cells[baseIdx + nz * StrideZ];
					if (list != null)
					{
						for (int i = 0; i < list.Count; i++)
							results.Add(list[i]);
					}
				}
			}
		}
	}

	public int CellCount => _occupiedCells.Count;
}
