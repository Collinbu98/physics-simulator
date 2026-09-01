using System;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// SPH kernel functions. All methods are pure math with no state.
/// Currently implements the Poly6 kernel for density estimation.
/// </summary>
public static class SPHKernels
{
    /// <summary>
    /// Poly6 kernel in 3D: W(r, h) = 315/(64π h⁹) × (h² - r²)³ for r ≤ h, else 0.
    /// 
    /// This kernel is symmetric, peaks at r=0, and drops to zero at r=h.
    /// It is used for density estimation because it evaluates efficiently
    /// when only squared distances are known (no sqrt needed).
    /// </summary>
    /// <param name="distSq">Squared distance between two particles.</param>
    /// <param name="h">Smoothing radius.</param>
    /// <returns>Kernel weight at this distance.</returns>
    public static float Poly6(float distSq, float h)
    {
        if (distSq >= h * h)
            return 0.0f;

        float h2 = h * h;
        float diff = h2 - distSq;
        float diff3 = diff * diff * diff;

        return Poly6Coefficient(h) * diff3;
    }

    /// <summary>
    /// Poly6 kernel evaluated with a precomputed coefficient.
    /// Use Poly6Coefficient(h) to obtain the coefficient once,
    /// then call this in tight loops to avoid redundant h⁹ computation.
    /// </summary>
    /// <param name="distSq">Squared distance between two particles.</param>
    /// <param name="h2">Smoothing radius squared (h*h).</param>
    /// <param name="coefficient">Precomputed 315/(64π h⁹).</param>
    /// <returns>Kernel weight at this distance.</returns>
    public static float Poly6(float distSq, float h2, float coefficient)
    {
        if (distSq >= h2)
            return 0.0f;

        float diff = h2 - distSq;
        float diff3 = diff * diff * diff;
        return coefficient * diff3;
    }

    /// <summary>
    /// Poly6 kernel evaluated at r=0. This is the self-contribution weight
    /// used when a particle contributes to its own density.
    /// W(0, h) = 315/(64π h³).
    /// </summary>
    /// <param name="h">Smoothing radius.</param>
    /// <returns>Kernel weight at zero distance.</returns>
    public static float Poly6AtOrigin(float h)
    {
        return 315.0f / (64.0f * MathF.PI * h * h * h);
    }

    /// <summary>
    /// Precomputed Poly6 coefficient: 315 / (64π h⁹).
    /// Call once before a loop, then pass to the overloaded Poly6(distSq, h2, coefficient).
    /// </summary>
    /// <param name="h">Smoothing radius.</param>
    /// <returns>The coefficient 315/(64π h⁹).</returns>
    public static float Poly6Coefficient(float h)
    {
        float h3 = h * h * h;
        float h9 = h3 * h3 * h3;
        return 315.0f / (64.0f * MathF.PI * h9);
    }
}
