#!/bin/bash
# Controlled benchmark: 5 runs per skin size, Release mode
set -e

SKINS="0.030 0.045 0.055 0.060"
RUNS=5

echo "=== Build Release ==="
dotnet build /home/collin/physics-simulator/Benchmark/Benchmark.csproj -c Release -v quiet 2>&1 | tail -1

for skin in $SKINS; do
    echo ""
    echo "========================================"
    echo "SKIN = $skin m  ($RUNS runs)"
    echo "========================================"
    for run in $(seq 1 $RUNS); do
        echo "--- Run $run ---"
        dotnet run --project /home/collin/physics-simulator/Benchmark/Benchmark.csproj -c Release --no-build -- $skin 2>&1
        echo ""
    done
done
