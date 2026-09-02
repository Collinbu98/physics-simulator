using Godot;

namespace PhysicsSimulator.Godot;

/// <summary>
/// Simple orbit/free-fly camera for inspecting particles.
/// Right-click drag: orbit around target. Scroll: zoom. WASD+QE: free-fly.
/// </summary>
public partial class CameraController : Camera3D
{
    /// <summary>
    /// World-space point the camera orbits around.
    /// </summary>
    public Vector3 OrbitTarget { get; set; } = Vector3.Zero;

    /// <summary>
    /// Current distance from the orbit target.
    /// </summary>
    private float _orbitDistance = 4.24f; // matches initial position (0,3,3) → dist ≈ 4.24

    /// <summary>
    /// Orbit angles in radians.
    /// </summary>
    private float _yaw;   // horizontal rotation
    private float _pitch; // vertical rotation (clamped)

    /// <summary>
    /// Orbit sensitivity (radians per pixel).
    /// </summary>
    private const float OrbitSensitivity = 0.005f;

    /// <summary>
    /// Zoom sensitivity (fraction of distance per scroll step).
    /// </summary>
    private const float ZoomFactor = 0.15f;

    /// <summary>
    /// Free-fly speed in meters per second.
    /// </summary>
    private const float FlySpeed = 2.0f;

    /// <summary>
    /// Fast fly speed (Shift held).
    /// </summary>
    private const float FastFlySpeed = 6.0f;

    /// <summary>
    /// Minimum orbit distance (prevents flipping through the target).
    /// </summary>
    private const float MinDistance = 0.1f;

    /// <summary>
    /// Maximum orbit distance.
    /// </summary>
    private const float MaxDistance = 100.0f;

    public override void _Ready()
    {
        // Compute initial orbit angles from current transform
        Vector3 toCamera = GlobalPosition - OrbitTarget;
        _orbitDistance = toCamera.Length();
        if (_orbitDistance < 0.001f)
        {
            _orbitDistance = 3.0f;
            toCamera = new Vector3(0, 0, 3);
        }

        _yaw = Mathf.Atan2(toCamera.X, toCamera.Z);
        _pitch = Mathf.Asin(toCamera.Y / _orbitDistance);
        ClampPitch();
        UpdateOrbitPosition();
    }

    /// <summary>
    /// Resets the camera to its initial position and orientation.
    /// </summary>
    public void ResetCamera()
    {
        OrbitTarget = Vector3.Zero;
        _orbitDistance = 4.24f;
        _yaw = Mathf.Atan2(0, 3); // matches initial (0,3,3) → yaw ≈ 0
        _pitch = Mathf.Asin(3.0f / 4.24f); // matches initial pitch
        ClampPitch();
        UpdateOrbitPosition();
    }

    public override void _UnhandledInput(InputEvent @event)
    {
        // Right-click drag to orbit
        if (@event is InputEventMouseMotion motion && motion.ButtonMask == MouseButtonMask.Right)
        {
            _yaw -= motion.Relative.X * OrbitSensitivity;
            _pitch -= motion.Relative.Y * OrbitSensitivity;
            ClampPitch();
            UpdateOrbitPosition();
            GetViewport().SetInputAsHandled();
            return;
        }

        // Scroll to zoom
        if (@event is InputEventMouseButton scroll && scroll.Pressed)
        {
            if (scroll.ButtonIndex == MouseButton.WheelUp)
            {
                _orbitDistance *= (1.0f - ZoomFactor);
                if (_orbitDistance < MinDistance) _orbitDistance = MinDistance;
                UpdateOrbitPosition();
                GetViewport().SetInputAsHandled();
                return;
            }
            if (scroll.ButtonIndex == MouseButton.WheelDown)
            {
                _orbitDistance *= (1.0f + ZoomFactor);
                if (_orbitDistance > MaxDistance) _orbitDistance = MaxDistance;
                UpdateOrbitPosition();
                GetViewport().SetInputAsHandled();
                return;
            }
        }
    }

    public override void _Process(double delta)
    {
        float speed = Input.IsKeyPressed(Key.Shift) ? FastFlySpeed : FlySpeed;
        float dt = (float)delta;

        Vector3 forward = -GlobalTransform.Basis.Z;
        Vector3 right = GlobalTransform.Basis.X;
        Vector3 up = GlobalTransform.Basis.Y;

        Vector3 move = Vector3.Zero;
        if (Input.IsKeyPressed(Key.W)) move += forward;
        if (Input.IsKeyPressed(Key.S)) move -= forward;
        if (Input.IsKeyPressed(Key.D)) move += right;
        if (Input.IsKeyPressed(Key.A)) move -= right;
        if (Input.IsKeyPressed(Key.E)) move += up;
        if (Input.IsKeyPressed(Key.Q)) move -= up;

        if (move.LengthSquared() > 0.0001f)
        {
            move = move.Normalized() * speed * dt;
            GlobalPosition += move;
            OrbitTarget = GlobalPosition + forward * _orbitDistance;
        }
    }

    private void ClampPitch()
    {
        if (_pitch > Mathf.Pi / 2.0f - 0.01f) _pitch = Mathf.Pi / 2.0f - 0.01f;
        if (_pitch < -Mathf.Pi / 2.0f + 0.01f) _pitch = -Mathf.Pi / 2.0f + 0.01f;
    }

    private void UpdateOrbitPosition()
    {
        float cosPitch = Mathf.Cos(_pitch);
        float x = _orbitDistance * cosPitch * Mathf.Sin(_yaw);
        float y = _orbitDistance * Mathf.Sin(_pitch);
        float z = _orbitDistance * cosPitch * Mathf.Cos(_yaw);

        GlobalPosition = OrbitTarget + new Vector3(x, y, z);
        LookAt(OrbitTarget, Vector3.Up);
    }
}
