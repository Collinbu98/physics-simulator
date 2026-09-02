using System;
using System.Numerics;

namespace PhysicsSimulator.Simulation;

/// <summary>
/// SPH kernel functions. All methods are pure math with no state.
/// Implements the Poly6 kernel for density estimation, Spiky kernel gradient for
/// pressure forces, and the viscosity kernel Laplacian for viscous diffusion.
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

    // -----------------------------------------------------------------------
    // Spiky kernel gradient
    // -----------------------------------------------------------------------

    /// <summary>
    /// Spiky kernel gradient in 3D:
    ///   W(r, h) = 15/(π h⁶) × (h - |r|)³
    ///   ∇W(r, h) = -45/(π h⁶) × (h - |r|)² × r̂
    /// where r̂ = r/|r| is the unit direction vector.
    ///
    /// The gradient points inward (toward the particle), so a positive pressure
    /// difference pushes particles apart.
    /// </summary>
    /// <param name="rVec">Vector from particle i to particle j (r_j - r_i).</param>
    /// <param name="h">Smoothing radius.</param>
    /// <returns>Gradient of the Spiky kernel at this separation.</returns>
    public static Vector3 SpikyGradient(Vector3 rVec, float h)
    {
        float distSq = rVec.LengthSquared();
        if (distSq >= h * h || distSq < 1e-20f)
            return Vector3.Zero;

        float dist = MathF.Sqrt(distSq);
        float hMinusR = h - dist;
        float hMinusR2 = hMinusR * hMinusR;
        float coeff = -45.0f / (MathF.PI * h * h * h * h * h * h);

        // r̂ = rVec / dist (unit direction)
        Vector3 rHat = rVec / dist;

        return coeff * hMinusR2 * rHat;
    }

    /// <summary>
    /// Spiky kernel gradient evaluated with a precomputed coefficient.
    /// Precompute with SpikyGradientCoefficient(h), then call this in tight loops.
    /// </summary>
    /// <param name="rVec">Vector from particle i to particle j (r_j - r_i).</param>
    /// <param name="h2">Smoothing radius squared (h*h).</param>
    /// <param name="coefficient">Precomputed -45/(π h⁶).</param>
    /// <returns>Gradient of the Spiky kernel at this separation.</returns>
    public static Vector3 SpikyGradient(Vector3 rVec, float h2, float coefficient)
    {
        float distSq = rVec.LengthSquared();
        if (distSq >= h2 || distSq < 1e-20f)
            return Vector3.Zero;

        float dist = MathF.Sqrt(distSq);
        float h = MathF.Sqrt(h2);
        float hMinusR = h - dist;
        float hMinusR2 = hMinusR * hMinusR;
        Vector3 rHat = rVec / dist;

        return coefficient * hMinusR2 * rHat;
    }

    /// <summary>
    /// Precomputed Spiky gradient coefficient: -45 / (π h⁶).
    /// Call once before a loop, then pass to SpikyGradient(rVec, h2, coefficient).
    /// </summary>
    /// <param name="h">Smoothing radius.</param>
    /// <returns>The coefficient -45/(π h⁶).</returns>
    public static float SpikyGradientCoefficient(float h)
    {
        float h6 = h * h * h * h * h * h;
        return -45.0f / (MathF.PI * h6);
    }

    // -----------------------------------------------------------------------
    // Viscosity kernel Laplacian
    // -----------------------------------------------------------------------

    /// <summary>
    /// Laplacian of the SPH viscosity kernel in 3D:
    ///   ∇²W_visc(r, h) = 45 / (π h⁶) × (h - |r|)  for r ≤ h, else 0.
    ///
    /// This kernel is designed so its Laplacian is positive everywhere inside the
    /// support radius, making it ideal for diffusion terms like viscosity.
    /// When used in the viscosity force:
    ///   F_i^visc = Σ_j m_j × ν × (v_j - v_i)/ρ_j × ∇²W(r_ij, h)
    /// it pulls particle i's velocity toward its neighbors' velocities.
    /// </summary>
    /// <param name="rVec">Vector from particle i to particle j (r_j - r_i).</param>
    /// <param name="h">Smoothing radius.</param>
    /// <returns>Laplacian of the viscosity kernel at this separation.</returns>
    public static float ViscosityLaplacian(Vector3 rVec, float h)
    {
        float distSq = rVec.LengthSquared();
        if (distSq >= h * h)
            return 0.0f;

        float dist = MathF.Sqrt(distSq);
        float coeff = 45.0f / (MathF.PI * h * h * h * h * h * h);
        return coeff * (h - dist);
    }

    /// <summary>
    /// Laplacian of the SPH viscosity kernel evaluated with a precomputed coefficient.
    /// Precompute with ViscosityLaplacianCoefficient(h), then call this in tight loops.
    /// </summary>
    /// <param name="dist">Distance between particles (not squared).</param>
    /// <param name="h">Smoothing radius.</param>
    /// <param name="coefficient">Precomputed 45/(π h⁶).</param>
    /// <returns>Laplacian of the viscosity kernel at this distance.</returns>
    public static float ViscosityLaplacian(float dist, float h, float coefficient)
    {
        if (dist >= h)
            return 0.0f;

        return coefficient * (h - dist);
    }

    /// <summary>
    /// Precomputed viscosity Laplacian coefficient: 45 / (π h⁶).
    /// Call once before a loop, then pass to ViscosityLaplacian(dist, h, coefficient).
    /// </summary>
    /// <param name="h">Smoothing radius.</param>
    /// <returns>The coefficient 45/(π h⁶).</returns>
    public static float ViscosityLaplacianCoefficient(float h)
    {
        float h6 = h * h * h * h * h * h;
        return 45.0f / (MathF.PI * h6);
    }
}
