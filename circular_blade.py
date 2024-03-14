import numpy as np
import matplotlib.pyplot as plt


def generate_circular_profile(H, Q, N, eff, ri, rh, rt, z):
    """
    Circular profile generator. It projects the 2D mean line into 3D coordinates

    Function arguments:
    H: head (m)
    Q: flow rate (m^3/s)
    N: angualar velocity (rpm)
    eff: turbine efficiency (-)
    ri: radius design point (m)
    rh: hub radius (m)
    rt: tip radius (m)
    z: number of blades
    """

    # Step 1. Project the circular profile to 2D cartesian coordinates x,y

    # Velocity triangle calculations

    # Constants
    # rho = 997 # Water density (kg/m^3)
    g = 9.81  # Gravitational acceleration (m/s^2)
    k = (g*H*eff*60) / (N*2*np.pi)  # Free vortex constant (m^2/s)
    wrap_angle = np.deg2rad(360 / z)  # Wrap angle (To insert in Bladegen for the first window in angle/thickness mode)
    omega = N * (2*np.pi/60)  # Angular velocity in rad/s

    # Velocities
    va = Q / (np.pi * (rt**2 - rh**2))  # Axial velocity
    vt = omega * ri  # Tangential velocity
    vc = k / ri  # Circumferential velocity

    # Leading and trailing edge angles
    beta1 = np.arctan(va / vt)
    beta2 = np.arctan(va / (vc + vt))

    # Circualar profile x,y coodinates calculation
    L = ri * wrap_angle
    x1 = -L / 2
    x2 = L / 2

    rc = L / ((np.sin(np.pi/2 - beta1) / (np.tan(np.pi/2 - beta1))) - np.sin(beta2))  # Profile circle radius
    xc = L/2 + rc*np.sin(beta2)  # Circle center at x coordinate
    yc = rc * np.sin(np.pi/2 - beta1)  # Circle center at y coordinate

    # Calculate the half axial chord for bladegen to use on the meridional view. In bladegen, the axial chord is Ca/2.
    ca = ((rc * np.cos(beta2) - yc)) * 1000  # mm

    # 2D projection x,y

    number_of_points = 40
    # x coordinate generation
    x_2d = []
    for i in range(number_of_points, -1, -1):
        xi_iter = x2 + (i/number_of_points) * (x1 - x2)
        x_2d.append(xi_iter)

    # y coordinate generation
    y_2d = yc - np.sqrt(rc**2 - (x_2d - xc)**2)

    # Displace x to the right by L/2 to make the first values of x,y equal 0. This is mandatory to transform to bladegen coord
    # x_2d = x_2d + L/2

    # ----------------------------------------------------------------------------------------------

    # Step 2. Project the profile to 3D cylindrical coordinates (ri, theta, z)

    # Projecting the camber

    theta = [i / ri for i in x_2d]
    # theta = x_2d/ri
    z_3d = y_2d

    # ----------------------------------------------------------------------------------------------

    # Step 3. Transform the cylindrical coordinates to x,y,z. The coordinate z is equal to z_cylindrical

    x_3d = ri*np.cos(theta)
    y_3d = ri*np.sin(theta)

    return x_3d, y_3d, z_3d, x_2d, y_2d, np.rad2deg(beta1), np.rad2deg(beta2), L, x1, x2, rc, xc, yc, ca


def cartesian_to_meridional(x, y, z):
    """Transforms 3D cartesian coordinates (x, y, z) to meridional coordinates (%m_prime, theta)"""

    m = []
    theta = []

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
        theta.append(theta_i*-1)  # theta_i was multiplied by -1 to obtain positive values for theta

    # Calculate %m_i
    max_m = max(m)
    m_i_percentage = [i/max_m for i in m]

    return m_i_percentage, theta


def interpolate_polynomial(m_prime_percentage, theta, points_to_evaluate):
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
    z: number of blades
"""

# Hydraulic parameters
H = 3  # m
Q = 0.015  # m^3\s
N = 3600  # rpm
eff = 0.65
ri_hub = 0.02259  # m
ri_tip = 0.03765  # m
z = 5
# ri_mid = 0.0318

# Do not change the following definitions
rh = ri_hub
rt = ri_tip
ri_mid = (ri_hub + ri_tip) / 2

# ----------------------------------------------------------------------------------------------

# Generate the 3D circular profile with user inputs

(
    x_3d_hub, y_3d_hub, z_3d_hub,
    x_hub_2d, y_hub_2d,
    beta1_hub, beta2_hub,
    L_hub, x1_hub, x2_hub, rc_hub, xc_hub, yc_hub, ca_hub
) = generate_circular_profile(H, Q, N, eff, ri_hub, rh, rt, z)

(
    x_3d_mid, y_3d_mid, z_3d_mid,
    x_mid_2d, y_mid_2d,
    beta1_mid, beta2_mid,
    L_mid, x1_mid, x2_mid, rc_mid, xc_mid, yc_mid, ca_mid
) = generate_circular_profile(H, Q, N, eff, ri_mid, rh, rt, z)

(
    x_3d_tip, y_3d_tip, z_3d_tip,
    x_tip_2d, y_tip_2d,
    beta1_tip, beta2_tip,
    L_tip, x1_tip, x2_tip, rc_tip, xc_tip, yc_tip, ca_tip
) = generate_circular_profile(H, Q, N, eff, ri_tip, rh, rt, z)

# Transform x,y,z to Bladegen meridional coordinates %m_prime, theta
m_hub, theta_hub = cartesian_to_meridional(x_3d_hub, y_3d_hub, z_3d_hub)
m_mid, theta_mid = cartesian_to_meridional(x_3d_mid, y_3d_mid, z_3d_mid)
m_tip, theta_tip = cartesian_to_meridional(x_3d_tip, y_3d_tip, z_3d_tip)

# # Interpolate the meridional coordinates (%m_prime, theta) with 5th grade polynomial at specified points
percentage_positions = [0.0, 0.25, 0.5, 0.75, 1.0]  # Points to evaluate the fitted fucntion theta(%m_prime)
interpolated_values_hub, r_squared_hub = interpolate_polynomial(m_hub, theta_hub, percentage_positions)
interpolated_values_mid, r_squared_mid = interpolate_polynomial(m_mid, theta_mid, percentage_positions)
interpolated_values_tip, r_squared_tip = interpolate_polynomial(m_tip, theta_tip, percentage_positions)

# ----------------------------------------------------------------------------------------------

# Print 2D coordinates
# print("2D x-coordinate")
# for i in x_tip_2d:
#     print(f"{i*1000:.2f}")

# print("2D y-coordinate")
# for i in y_tip_2d:
#     print(f"{i*1000:.2f}")

# Print 3D coordinates
# print("z-coordinate")
# for i in z_3d_tip:
#     print(i)
# print("x-coordinate")
# for i in x_3d_tip:
#     print(i)
# print("y-coordinate")
# for i in y_3d_tip:
#     print(i)

# Print circular profile geometry results
print(f"""beta1_hub : {beta1_hub:.4f} °, beta2_hub : {beta2_hub:.4f} °, L_hub : {L_hub:.4f}, x1_hub : {x1_hub:.4f}, x2_hub : {x2_hub:.4f},
      rc_hub : {rc_hub:.4f}, xc_hub : {xc_hub:.4f}, yc_hub : {yc_hub:.4f}""")
print(f"""beta1_mid : {beta1_mid:.4f} °, beta2_mid : {beta2_mid:.4f} °, L_mid : {L_mid:.4f}, x1_mid : {x1_mid:.4f}, x2_mid : {x2_mid:.4f},
      rc_mid : {rc_mid:.4f}, xc_mid : {xc_mid:.4f}, yc_mid : {yc_mid:.4f}""")
print(f"""beta1_tip : {beta1_tip:.4f} °, beta2_tip : {beta2_tip:.4f} °, L_tip : {L_tip:.4f}, x1_tip : {x1_tip:.4f}, x2_tip : {x2_tip:.4f},
      rc_tip : {rc_tip:.4f}, xc_tip : {xc_tip:.4f}, yc_tip : {yc_tip:.4f}\n""")
print("---------------------------------")

# Print R^2 for the interpolated values of theta(%m_prime)
print("Coefficient of determination R^2 for the interpolated values of %m vs theta")
print(f"R^2 at hub = {r_squared_hub}")
print(f"R^2 at mid = {r_squared_mid}")
print(f"R^2 at tip = {r_squared_tip}\n")
print("---------------------------------")

# Print the blade axial chord Ca/2 for bladegen in mm
print("Axial chord Ca/2 for the Bladegen meridional view:")
print(f"Ca/2_hub = {ca_hub/2:.4f}")
print(f"Ca/2_mid = {ca_mid/2:.4f}")
print(f"Ca/2_tip = {ca_tip/2:.4f}")
print("---------------------------------")

# Printing interpolated values of the function theta(%m_prime) and R^2
for i, pos in enumerate(percentage_positions):
    print(f"Interpolated value at %m_prime = {int(pos*100)}%")
    print(f"theta_hub = {interpolated_values_hub[i]:.4f} ")
    print(f"theta_mid = {interpolated_values_mid[i]:.4f} ")
    print(f"theta_tip = {interpolated_values_tip[i]:.4f} ")
    print("---------------------------------")

# ----------------------------------------------------------------------------------------------

# 2D plotting
plt.figure(figsize=(16, 4))
plt.plot(x_hub_2d, y_hub_2d, label="Hub")
plt.plot(x_mid_2d, y_mid_2d, label="Mid")
plt.plot(x_tip_2d, y_tip_2d, label="Tip")
plt.title("2D cartesian projection of circular profile")
plt.xlabel("x")
plt.ylabel("y")
# plt.xlim(0,)
plt.grid(True)
# plt.axis("scaled")
plt.legend()
plt.show()

# 3D plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', proj_type="ortho")
ax.plot(x_3d_hub, y_3d_hub, z_3d_hub, label="Hub")
ax.plot(x_3d_mid, y_3d_mid, z_3d_mid, label="Mid")
ax.plot(x_3d_tip, y_3d_tip, z_3d_tip, label="Tip")
ax.set_title("3D projection of circular profiles")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_aspect("equal")
ax.legend()
plt.show()

# Export profile 3D data in mm for the thesis methodology
# np.savetxt("x_3D_hub.txt", x_3d_hub * 1000, delimiter=",")
# np.savetxt("y_3D_hub.txt", y_3d_hub * 1000, delimiter=",")
# np.savetxt("z_3D_hub.txt", z_3d_hub * 1000, delimiter=",")

# np.savetxt("x_3D_mid.txt", x_3d_mid * 1000, delimiter=",")
# np.savetxt("y_3D_mid.txt", y_3d_mid * 1000, delimiter=",")
# np.savetxt("z_3D_mid.txt", z_3d_mid * 1000, delimiter=",")

# np.savetxt("x_3D_tip.txt", x_3d_tip * 1000, delimiter=",")
# np.savetxt("y_3D_tip.txt", y_3d_tip * 1000, delimiter=",")
# np.savetxt("z_3D_tip.txt", z_3d_tip * 1000, delimiter=",")

# ----------------------------------------------------------------------------------------------

# # Plotting %m_prime vs. theta
plt.plot(m_hub, theta_hub, label="Hub")
plt.plot(m_mid, theta_mid, label="Mid")
plt.plot(m_tip, theta_tip, label="Tip")
plt.xticks(percentage_positions)
# plt.xlim(0, 1)
# plt.ylim(0,)
plt.xlabel(r"$\%m_{prime}$")
plt.ylabel(r"$\theta$")
plt.grid(True)
plt.legend()
plt.show()

# # Export meridional coordinates data for the thesis methodology
# np.savetxt("m_hub.txt", m_hub, delimiter=",")
# np.savetxt("m_mid.txt", m_mid, delimiter=",")
# np.savetxt("m_tip.txt", m_tip, delimiter=",")

# np.savetxt("theta_hub.txt", theta_hub, delimiter=",")
# np.savetxt("theta_mid.txt", theta_mid, delimiter=",")
# np.savetxt("theta_tip.txt", theta_tip, delimiter=",")
