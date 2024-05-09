import numpy as np
import matplotlib.pyplot as plt


def calculate_velocity_triangle_and_chord_length(H: float, Q: float, N: float, eff: float, ri: float, rh: float, rt: float, z: float) -> float:
    """
    Resolves the velocity triangle assuming free vortex at leading edge. Calculates the airfoil chord length and
    the blade stagger angle.

    Function arguments:
    H: head (m)
    Q: flow rate (m^3/s)
    N: angualar velocity (rpm)
    eff: turbine efficiency (-)
    rh: hub radius (m)
    rt: tip radius (m)
    z: number of blades
    """

    # Constants
    g = 9.81  # Gravitational acceleration (m/s^2)
    k = (g*H*eff*60) / (N*2*np.pi)  # Free vortex constant (m^2/s)
    omega = N * (2*np.pi/60)  # Angular velocity in rad/s
    cl = 0.57692  # Lift coefficient
    lift_to_drag = 59  # Drag-to-lift coefficient
    cd = cl / lift_to_drag  # Drag coefficient
    drag_to_lift = cd / cl

    # Velocities
    va = Q / (np.pi * (rt**2 - rh**2))  # Axial velocity
    vt = ri * omega  # Tangential velocity
    vc = k / ri  # Circumferential velocity

    # Blade angles
    beta1 = np.arctan(va / (vt - vc))
    beta2 = np.arctan(va / vt)

    # Mean beta angle
    beta_mean = np.arctan(0.5 * (np.tan(beta1) + np.tan(beta2)))

    # theta angle
    theta = beta_mean - np.deg2rad(2)  # AoA = 2°

    # Stagger angle
    stagger_angle = np.deg2rad(270) - theta

    # Mean relative velocity
    w_mean = va / np.cos(beta_mean)

    # Blade pitch at ri
    pitch = (2*np.pi*ri) / z

    # Chord length chord_length
    chord_length = (2 * g * pitch * eff * H) / (cl * vt * w_mean * (1-(drag_to_lift/np.tan(beta_mean))))

    # Calculate axial chord for bladegen to use on the meridional view. In bladegen, the axial chord is Ca/2.
    ca = (chord_length * np.sin(np.deg2rad(90) - theta)) * 1000  # Converted to mm

    # Blade solidity
    solidity = chord_length / pitch

    return va, vt, vc, np.rad2deg(beta1), np.rad2deg(beta2), np.rad2deg(beta_mean), np.rad2deg(stagger_angle), chord_length, solidity, ca


def generate_naca_4_digit(C: float, ri: float, stagger_angle: float) -> float:
    """
    NACA 4-digit airfoil generator. It projects the 2D profile into 3D coordinates

    Function arguments:
    C : chord length
    r : design radius
    gamma : stagger angle

    Variables subindex meaning
    c : camber
    u : upper
    chord_length : lower
    """

    # Step 1. Project the airfoil profile to a x,y 2D coordinate
    # MEL031 airfoil parameters
    M = 0.04  # Maximum camber
    P = 0.5  # Maximum camber location
    T = 0.12  # Maximum thickness

    # Thickness constants
    A0 = 0.2969
    A1 = -0.126
    A2 = -0.3516
    A3 = 0.2843
    A4 = -0.1015  # Opened trailing edge. Change to -0.1036 for a closed trailing edge.

    # Normalized x-coordinate. Cosine spacing is used to group more points at the leading and trailing edges
    beta = np.linspace(0, np.pi, 100)
    x = (1-np.cos(beta))/2

    # Camber line yc and its derivative dyc_dx
    yc_2d = np.where(x < P, (M/(P**2)) * (2*P*x - x**2), (M/(1-P)**2) * (1 - 2*P + 2*P*x - x**2))
    dyc_dx = np.where(x < P, (2*M/(P**2)) * (P - x), (2*M/(1-P)**2) * (P - x))

    camber_2d = yc_2d * C  # Variable created to plot the 2D profile

    # Thickness distribution yt
    yt = 5*T*C * (A0*np.sqrt(x) + A1*(x) + A2*(x)**2 + A3*(x)**3 + A4*(x)**4)

    # Upper surface coordinates
    xu_2d = x * C - yt * np.sin(np.arctan(dyc_dx))
    yu_2d = yc_2d * C + yt * np.cos(np.arctan(dyc_dx))

    # Lower surface coordinates
    xl_2d = x * C + yt * np.sin(np.arctan(dyc_dx))
    yl_2d = yc_2d * C - yt * np.cos(np.arctan(dyc_dx))

    # ----------------------------------------------------------------------------------------------

    # Step 2. Rotate the profiles by 270-theta. Where theta = beta_avg - AoA

    # rotating the camber
    xc_2d_rotated = xl_2d*np.cos(np.deg2rad(stagger_angle)) - yc_2d*C*np.sin(np.deg2rad(stagger_angle))
    yc_2d_rotated = xl_2d*np.sin(np.deg2rad(stagger_angle)) + yc_2d*C*np.cos(np.deg2rad(stagger_angle))

    # Rotating the upper
    xu_2d_rotated = xu_2d*np.cos(np.deg2rad(stagger_angle)) - yu_2d*np.sin(np.deg2rad(stagger_angle))
    yu_2d_rotated = xu_2d*np.sin(np.deg2rad(stagger_angle)) + yu_2d*np.cos(np.deg2rad(stagger_angle))

    # Rotating thde lower
    xl_2d_rotated = xl_2d*np.cos(np.deg2rad(stagger_angle)) - yl_2d*np.sin(np.deg2rad(stagger_angle))
    yl_2d_rotated = xl_2d*np.sin(np.deg2rad(stagger_angle)) + yl_2d*np.cos(np.deg2rad(stagger_angle))

    # ----------------------------------------------------------------------------------------------

    # Step 3. Project the profile to 3D cylindrical coordinates (r, theta, z)

    # Projecting the camber
    theta_camber = xc_2d_rotated/ri
    zc_3d = yc_2d_rotated

    # Projecting to cylindrical coordinates the upper surface
    theta_upper = xu_2d_rotated/ri
    zu_3d = yu_2d_rotated

    # Projecting to cylindrical coordinates the lower surface
    theta_lower = xl_2d_rotated/ri
    zl_3d = yl_2d_rotated

    # ----------------------------------------------------------------------------------------------

    # Step 4. Transform the cylindrical coordinates to x,y,z. The coordinate z is equal to z_cylindrical

    xc_3d = ri*np.cos(theta_camber)
    yc_3d = ri*np.sin(theta_camber)

    xu_3d = ri*np.cos(theta_upper)
    yu_3d = ri*np.sin(theta_upper)

    xl_3d = ri*np.cos(theta_lower)
    yl_3d = ri*np.sin(theta_lower)

    return (
        # 2D profile (not rotated)
        camber_2d, xu_2d, yu_2d, xl_2d, yl_2d,
        # 2D rotated profile
        xc_2d_rotated, yc_2d_rotated, xu_2d_rotated, yu_2d_rotated, xl_2d_rotated, yl_2d_rotated,
        # 3D profile projection
        xc_3d, yc_3d, zc_3d, xu_3d, yu_3d, zu_3d, xl_3d, yl_3d, zl_3d
    )


def cartesian_to_meridional(x: float, y: float, z: float) -> float:
    """Transforms 3D cartesian coordinates (x, y, z) to meridional coordinates (m, theta)"""

    m = []  # Store the m_i values
    theta = []  # Stores the theta_i values

    # Initial values for m and r (Assuming the first node as the origin)
    m_prev = 0
    r_prev = 0

    for i in range(len(x)):
        # Calculate r_i for the current node
        r_i = np.sqrt(x[i]**2 + y[i]**2)

        if i == 0:
            # Special case for the first node
            m_i = 0
        else:
            m_i = m_prev + np.sqrt((z[i]-z[i-1])**2 + (r_i-r_prev)**2) / r_i

        # Calculate theta_i for the current node
        theta_i = np.rad2deg((np.arctan2(x[i], y[i]) - np.arctan2(x[0], y[0])))  # Theta_i in degrees [°]

        m_prev = m_i
        r_prev = r_i

        m.append(m_i)
        theta.append(theta_i)

    # Calculate %m_i
    max_m = max(m)
    m_i_percentage = [i/max_m for i in m]

    return m_i_percentage, theta


def interpolate_polynomial(m_prime_percentage: float, theta: float, points_to_evaluate: float) -> float:
    """Interpolates the meridional coordinates (%m_prime, theta) and evaluates at %m_prime(0, 25, 50, 75, 100) %"""

    # 5th grade polynomial fit
    coefficients = np.polyfit(m_prime_percentage, theta, 5)

    # Evaluate the polynomial at all data points and points_to_evaluate
    interpolated_values_all = np.polyval(coefficients, m_prime_percentage)
    interpolated_values_points = np.polyval(coefficients, points_to_evaluate)

    # Calculate R^2
    mean_theta = np.mean(theta)
    total_sum_of_squares = np.sum((theta - mean_theta)**2)
    residual_sum_of_squares = np.sum((interpolated_values_all - theta)**2)
    r_squared = 1 - (residual_sum_of_squares / total_sum_of_squares)

    return interpolated_values_points, r_squared

# ----------------------------------------------------------------------------------------------


"""
User inputs:
    H: head (m)
    Q: flow rate (m^3/s)
    N: angualar velocity (rpm)
    eff: turbine efficiency (-)
    r_hub: hub radius (m)
    r_tip: tip radius (m)
    z: number of blades (-)
"""

# Hydraulic parameters: the followin parameters are required:
H = 3  # m
Q = 0.015  # m^3/s
N = 3600  # rpm
eff = 0.65
r_hub = 0.02259  # m
r_tip = 0.03765  # m
z = 5

# Do not change the following definitions
rh = r_hub
rt = r_tip
r_mid = (r_hub + r_tip) / 2

# ----------------------------------------------------------------------------------------------

# Solve the turbine velocity triangle and calculate the chord length and the stagger angle
(
    va, vt_hub, vc_hub,
    beta1_hub, beta2_hub, beta_mean_hub, stagger_angle_hub,
    chord_length_hub, solidity_hub, ca_hub
) = calculate_velocity_triangle_and_chord_length(H, Q, N, eff, r_hub, rh, rt, z)

(
    va, vt_mid, vc_mid,
    beta1_mid, beta2_mid, beta_mean_mid, stagger_angle_mid,
    chord_length_mid, solidity_mid, ca_mid
) = calculate_velocity_triangle_and_chord_length(H, Q, N, eff, r_mid, rh, rt, z)

(
    va, vt_tip, vc_tip,
    beta1_tip, beta2_tip, beta_mean_tip, stagger_angle_tip,
    chord_length_tip, solidity_tip, ca_tip
) = calculate_velocity_triangle_and_chord_length(H, Q, N, eff, r_tip, rh, rt, z)


# Generate the MEL031 4-digit profile
(
    camber_2d_hub, xu_2d_hub, yu_2d_hub, xl_2d_hub, yl_2d_hub,
    xc_2d_rotated_hub, yc_2d_rotated_hub, xu_2d_rotated_hub, yu_2d_rotated_hub, xl_2d_rotated_hub, yl_2d_rotated_hub,
    xc_3d_hub, yc_3d_hub, zc_3d_hub, xu_3d_hub, yu_3d_hub, zu_3d_hub, xl_3d_hub, yl_3d_hub, zl_3d_hub
) = generate_naca_4_digit(chord_length_hub, r_hub, stagger_angle_hub)

(
    camber_2d_mid, xu_2d_mid, yu_2d_mid, xl_2d_mid, yl_2d_mid,
    xc_2d_rotated_mid, yc_2d_rotated_mid, xu_2d_rotated_mid, yu_2d_rotated_mid, xl_2d_rotated_mid, yl_2d_rotated_mid,
    xc_3d_mid, yc_3d_mid, zc_3d_mid, xu_3d_mid, yu_3d_mid, zu_3d_mid, xl_3d_mid, yl_3d_mid, zl_3d_mid
) = generate_naca_4_digit(chord_length_mid, r_mid, stagger_angle_mid)

(
    camber_2d_tip, xu_2d_tip, yu_2d_tip, xl_2d_tip, yl_2d_tip,
    xc_2d_rotated_tip, yc_2d_rotated_tip, xu_2d_rotated_tip, yu_2d_rotated_tip, xl_2d_rotated_tip, yl_2d_rotated_tip,
    xc_3d_tip, yc_3d_tip, zc_3d_tip, xu_3d_tip, yu_3d_tip, zu_3d_tip, xl_3d_tip, yl_3d_tip, zl_3d_tip
) = generate_naca_4_digit(chord_length_tip, r_tip, stagger_angle_tip)


# Transform x,y,z to Bladegen meridional coordinates m, theta
m_c_hub, theta_c_hub = cartesian_to_meridional(xc_3d_hub, yc_3d_hub, zc_3d_hub)
m_c_mid, theta_c_mid = cartesian_to_meridional(xc_3d_mid, yc_3d_mid, zc_3d_mid)
m_c_tip, theta_c_tip = cartesian_to_meridional(xc_3d_tip, yc_3d_tip, zc_3d_tip)

# Interpolate the meridional coordinates (%m_prime, theta) with 5th grade polynomial at specified points
percentage_positions = [0.0, 0.25, 0.5, 0.75, 1.0]
interpolated_values_hub, r_squared_hub = interpolate_polynomial(m_c_hub, theta_c_hub, percentage_positions)
interpolated_values_mid, r_squared_mid = interpolate_polynomial(m_c_mid, theta_c_mid, percentage_positions)
interpolated_values_tip, r_squared_tip = interpolate_polynomial(m_c_tip, theta_c_tip, percentage_positions)

# ----------------------------------------------------------------------------------------------

# Printing the velocity triangle results and the turbine geometry parameters
print(f"Axial velocity: {va:.4f} m/s\n")
print("Tangential velocities (m/s):")
print(f"Hub: {vt_hub:.4f}\nMid: {vt_mid:.4f}\nTip: {vt_tip:.4f}\n")

print("Circumferential velocities (m/s):")
print(f"Hub: {vc_hub:.4f}\nMid: {vc_mid:.4f}\nTip: {vc_tip:.4f}\n")

print("Leading edge angle beta 1 (°):")
print(f"Hub: {beta1_hub:.4f}\nMid: {beta1_mid:.4f}\nTip: {beta1_tip:.4f}\n")
print("Trailing edge angle beta 2 (°):")
print(f"Hub: {beta2_hub:.4f}\nMid: {beta2_mid:.4f}\nTip: {beta2_tip:.4f}\n")

print("Beta mean angle (°):")
print(f"Hub: {beta_mean_hub:.4f}\nMid: {beta_mean_mid:.4f}\nTip: {beta_mean_tip:.4f}\n")

print("Stagger angles (°):")
print(f"Hub: {stagger_angle_hub:.4f}\nMid: {stagger_angle_mid:.4f}\nTip: {stagger_angle_tip:.4f}\n")

print("Solidity (-):")
print(f"Hub: {solidity_hub:.4f}\nMid: {solidity_mid:.4f}\nTip: {solidity_tip:.4f}\n")

print("Chord length (m):")
print(f"Hub: {chord_length_hub:.5f}\nMid: {chord_length_mid:.5f}\nTip: {chord_length_tip:.5f}\n")

# ----------------------------------------------------------------------------------------------

# Printing interpolated values of the function theta(%m_prime) and R^2
print("Bladegen meriodional coordinates (%m_prime vs. theta):\n\n")
print("Coefficients of determination for interpolation of coordinates:")
print(f"R^2 at hub = {r_squared_hub}")
print(f"R^2 at mid = {r_squared_mid}")
print(f"R^2 at tip = {r_squared_tip}")
print("---------------------------------")

# Printing half axial chord length "ca/2" to be used in the Bladegen Meridional view
print("Half axial chord \"ca/2\" for Bladegen (mm):")
print(f"Hub: {ca_hub/2:.4f}\nMid: {ca_mid/2:.4f}\nTip: {ca_tip/2:.4f}")
print("---------------------------------")

# Printing Bladegen coordinates (theta vs. %m_prime)
for i, pos in enumerate(percentage_positions):
    print(f"Interpolated value at %m_prime = {int(pos*100)}%")
    print(f"theta_hub = {interpolated_values_hub[i]:.4f} ")
    print(f"theta_mid = {interpolated_values_mid[i]:.4f} ")
    print(f"theta_tip = {interpolated_values_tip[i]:.4f} ")
    print("---------------------------------")

# ----------------------------------------------------------------------------------------------

# 2D plots
# Profile without rotation

plt.figure(figsize=(16, 4))
plt.plot(xl_2d_hub, camber_2d_hub, "g", label="Camber")
plt.plot(xu_2d_hub, yu_2d_hub, "r", label="Upper surface")
plt.plot(xl_2d_hub, yl_2d_hub, "b", label="Lower surface")
plt.title(fr"MEL031 profile at $r_{{hub}}$. Chord length C={chord_length_hub:.5f} m")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.axis("scaled")
plt.grid(True)
plt.show()

# # Scatter
plt.figure(figsize=(16, 4))
plt.scatter(xl_2d_hub, camber_2d_hub, c="g", s=5, label="Camber")
plt.scatter(xu_2d_hub, yu_2d_hub, c="r", s=5, label="Upper surface")
plt.scatter(xl_2d_hub, yl_2d_hub, c="b", s=5, label="Lower surface")
plt.title("Visualization of cosine spacing. It generates more points at the leading and trailing edge")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.axis("scaled")
plt.grid(True)
plt.show()

# Export profile 2D data in mm for the thesis methodology
# np.savetxt("lower_x.txt", xl_2d_hub * 1000, delimiter=",")
# np.savetxt("camber_y.txt", camber_2d_hub * 1000, delimiter=",")
# np.savetxt("lower_y.txt", yl_2d_hub * 1000, delimiter=",")
# np.savetxt("upper_x.txt", xu_2d_hub * 1000, delimiter=",")
# np.savetxt("upper_y.txt", yu_2d_hub * 1000, delimiter=",")


# ----------------------------------------------------------------------------------------------

# 2D rotated profiles

plt.figure(figsize=(9, 8))
plt.plot(xu_2d_rotated_hub, yu_2d_rotated_hub, "b", label="Hub")
plt.plot(xl_2d_rotated_hub, yl_2d_rotated_hub, "b")
plt.plot(xu_2d_rotated_mid, yu_2d_rotated_mid, "g", label="Mid")
plt.plot(xl_2d_rotated_mid, yl_2d_rotated_mid, "g")
plt.plot(xu_2d_rotated_tip, yu_2d_rotated_tip, "r", label="Tip")
plt.plot(xl_2d_rotated_tip, yl_2d_rotated_tip, "r")
plt.title(r"Rotated profiles at the stagger angle $\Gamma = 270°-\lambda$")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.axis("scaled")
plt.grid(True)
plt.show()

# Export rotated profile 2D data in mm for the thesis methodology
# Upper and lower surfaces
# np.savetxt("xu_rotated_hub.txt", xu_2d_rotated_hub * 1000, delimiter=",")
# np.savetxt("yu_rotated_hub.txt", yu_2d_rotated_hub * 1000, delimiter=",")
# np.savetxt("xl_rotated_hub.txt", xl_2d_rotated_hub * 1000, delimiter=",")
# np.savetxt("yl_rotated_hub.txt", yl_2d_rotated_hub * 1000, delimiter=",")

# np.savetxt("xu_rotated_mid.txt", xu_2d_rotated_mid * 1000, delimiter=",")
# np.savetxt("yu_rotated_mid.txt", yu_2d_rotated_mid * 1000, delimiter=",")
# np.savetxt("xl_rotated_mid.txt", xl_2d_rotated_mid * 1000, delimiter=",")
# np.savetxt("yl_rotated_mid.txt", yl_2d_rotated_mid * 1000, delimiter=",")

# np.savetxt("xu_rotated_tip.txt", xu_2d_rotated_tip * 1000, delimiter=",")
# np.savetxt("yu_rotated_tip.txt", yu_2d_rotated_tip * 1000, delimiter=",")
# np.savetxt("xl_rotated_tip.txt", xl_2d_rotated_tip * 1000, delimiter=",")
# np.savetxt("yl_rotated_tip.txt", yl_2d_rotated_tip * 1000, delimiter=",")

# # Camber line
# np.savetxt("xc_rotated_hub.txt", xc_2d_rotated_hub * 1000, delimiter=",")
# np.savetxt("yc_rotated_hub.txt", yc_2d_rotated_hub * 1000, delimiter=",")

# np.savetxt("xc_rotated_mid.txt", xc_2d_rotated_mid * 1000, delimiter=",")
# np.savetxt("yc_rotated_mid.txt", yc_2d_rotated_mid * 1000, delimiter=",")

# np.savetxt("xc_rotated_tip.txt", xc_2d_rotated_tip * 1000, delimiter=",")
# np.savetxt("yc_rotated_tip.txt", yc_2d_rotated_tip * 1000, delimiter=",")

# ----------------------------------------------------------------------------------------------

# 3D Plotting

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', proj_type="ortho")

# Plotting cambers
ax.plot(xc_3d_hub, yc_3d_hub, zc_3d_hub, "g", label="Camber")
ax.plot(xc_3d_mid, yc_3d_mid, zc_3d_mid, "g")
ax.plot(xc_3d_tip, yc_3d_tip, zc_3d_tip, "g")

# Plotting upper surfaces
ax.plot(xu_3d_hub, yu_3d_hub, zu_3d_hub, "r", label="Upper surface")
ax.plot(xu_3d_mid, yu_3d_mid, zu_3d_mid, "r")
ax.plot(xu_3d_tip, yu_3d_tip, zu_3d_tip, "r")

# Plotting lower surfaces
ax.plot(xl_3d_hub, yl_3d_hub, zl_3d_hub, "b", label="Lower surface")
ax.plot(xl_3d_mid, yl_3d_mid, zl_3d_mid, "b")
ax.plot(xl_3d_tip, yl_3d_tip, zl_3d_tip, "b")

# Plot settings
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title(fr"MEL031 3D projected profiles. $C_{{hub}}$={chord_length_hub:.5f} m, $C_{{mid}}$={chord_length_mid:.5f} m, $C_{{tip}}$={chord_length_tip:.5f} m")
ax.legend()
plt.show()

# Export rotated profile 2D data in mm for the thesis methodology
# Camber line
# np.savetxt("xc_3d_hub.txt", xc_3d_hub * 1000, delimiter=",")
# np.savetxt("yc_3d_hub.txt", yc_3d_hub * 1000, delimiter=",")
# np.savetxt("zc_3d_hub.txt", zc_3d_hub * 1000, delimiter=",")

# np.savetxt("xc_3d_mid.txt", xc_3d_mid * 1000, delimiter=",")
# np.savetxt("yc_3d_mid.txt", yc_3d_mid * 1000, delimiter=",")
# np.savetxt("zc_3d_mid.txt", zc_3d_mid * 1000, delimiter=",")

# np.savetxt("xc_3d_tip.txt", xc_3d_tip * 1000, delimiter=",")
# np.savetxt("yc_3d_tip.txt", yc_3d_tip * 1000, delimiter=",")
# np.savetxt("zc_3d_tip.txt", zc_3d_tip * 1000, delimiter=",")

# # Upper surface
# np.savetxt("xu_3d_hub.txt", xu_3d_hub * 1000, delimiter=",")
# np.savetxt("yu_3d_hub.txt", yu_3d_hub * 1000, delimiter=",")
# np.savetxt("zu_3d_hub.txt", zu_3d_hub * 1000, delimiter=",")

# np.savetxt("xu_3d_mid.txt", xu_3d_mid * 1000, delimiter=",")
# np.savetxt("yu_3d_mid.txt", yu_3d_mid * 1000, delimiter=",")
# np.savetxt("zu_3d_mid.txt", zu_3d_mid * 1000, delimiter=",")

# np.savetxt("xu_3d_tip.txt", xu_3d_tip * 1000, delimiter=",")
# np.savetxt("yu_3d_tip.txt", yu_3d_tip * 1000, delimiter=",")
# np.savetxt("zu_3d_tip.txt", zu_3d_tip * 1000, delimiter=",")

# # Lower surface
# np.savetxt("xl_3d_hub.txt", xl_3d_hub * 1000, delimiter=",")
# np.savetxt("yl_3d_hub.txt", yl_3d_hub * 1000, delimiter=",")
# np.savetxt("zl_3d_hub.txt", zl_3d_hub * 1000, delimiter=",")

# np.savetxt("xl_3d_mid.txt", xl_3d_mid * 1000, delimiter=",")
# np.savetxt("yl_3d_mid.txt", yl_3d_mid * 1000, delimiter=",")
# np.savetxt("zl_3d_mid.txt", zl_3d_mid * 1000, delimiter=",")

# np.savetxt("xl_3d_tip.txt", xl_3d_tip * 1000, delimiter=",")
# np.savetxt("yl_3d_tip.txt", yl_3d_tip * 1000, delimiter=",")
# np.savetxt("zl_3d_tip.txt", zl_3d_tip * 1000, delimiter=",")

# ----------------------------------------------------------------------------------------------

# Plotting %m_prime vs. theta
plt.plot(m_c_hub, theta_c_hub, label="Hub")
plt.plot(m_c_mid, theta_c_mid, label="Mid")
plt.plot(m_c_tip, theta_c_tip, label="Tip")
plt.xticks(percentage_positions)
plt.xlim(0, 1)
plt.ylim(0,)
plt.title("Bladegen meridional coordinates")
plt.xlabel(r"$\%m_{prime}$")
plt.ylabel(r"$\theta$")
plt.grid(True)
plt.legend()
plt.show()

# Export meridional coordinates data for the thesis methodology
# np.savetxt("m_c_hub.txt", m_c_hub, delimiter=",")
# np.savetxt("m_c_mid.txt", m_c_mid, delimiter=",")
# np.savetxt("m_c_tip.txt", m_c_tip, delimiter=",")

# np.savetxt("theta_c_hub.txt", theta_c_hub, delimiter=",")
# np.savetxt("theta_c_mid.txt", theta_c_mid, delimiter=",")
# np.savetxt("theta_c_tip.txt", theta_c_tip, delimiter=",")
