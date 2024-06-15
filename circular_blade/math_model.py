import numpy as np


def generate_circular_profile(H: float, Q: float, N: float, eff: float, ri: float, rh: float, rt: float, z: float) -> float:
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
    wrap_angle = np.deg2rad(360 / z)  # Wrap angle (Input parameter for Bladegen in angle/thickness mode)
    omega = N * (2*np.pi/60)  # Angular velocity in rad/s

    # Velocities
    va = Q / (np.pi * (rt**2 - rh**2))  # Axial velocity
    vt = omega * ri  # Tangential velocity
    vc = k / ri  # Circumferential velocity

    # Leading and trailing edge angles
    beta1 = np.arctan(va / vt)
    beta2 = np.arctan(va / (vc + vt))
    beta_avg = np.arctan(0.5*(np.tan(beta1) + np.tan(beta2)))

    # Circular profile x,y coodinates calculation
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

    # Displace x to the right by L/2 to make the first values of x,y equal 0.
    # x_2d = x_2d + L/2

    # Step 2. Project the profile to 3D cylindrical coordinates (ri, theta, z)

    theta = [i / ri for i in x_2d]
    z_3d = y_2d

    # Step 3. Transform the cylindrical coordinates to x,y,z. The coordinate z is equal to z_cylindrical

    x_3d = ri*np.cos(theta)
    y_3d = ri*np.sin(theta)

    return x_3d, y_3d, z_3d, x_2d, y_2d, np.rad2deg(beta1), np.rad2deg(beta2), L, x1, x2, rc, xc, yc, ca, np.rad2deg(beta_avg)


def cartesian_to_meridional(x: float, y: float, z: float) -> float:
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
        theta_i = np.rad2deg((np.arctan2(y[i], x[i]) - np.arctan2(y[0], x[0])))  # Theta_i in degrees [°]

        m_prev = m_i
        r_prev = r_i

        m.append(m_i)
        theta.append(theta_i)  # theta_i was multiplied by -1 to obtain positive values for theta

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


def cord_length(x1: float, x2: float, y1: float, y2: float) -> float:
    """Returns the profile's cord length from the 2D coordinates of the leading and trailing edge"""
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
