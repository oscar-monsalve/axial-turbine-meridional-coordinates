import numpy as np


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
