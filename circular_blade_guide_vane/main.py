import matplotlib.pyplot as plt
import math_model

if __name__ == "__main__":

    """
    User inputs:
        H: head (m)
        Q: flow rate (m^3/s)
        N: angular velocity (rpm)
        eff: turbine efficiency (-)
        r_hub: hub radius (m)
        r_tip: tip radius (m)
        z: number of blades (-)
    """

    # Hydraulic parameters: the followin parameters are required:
    H = 3  # m
    Q = 0.015  # m^3\s
    N = 3600  # rpm
    eff = 0.65
    ri_hub = 0.02259  # m
    ri_tip = 0.03765  # m
    z = 5
    # todo: generalize the anti-clock-wise pre-rotation to set the blade leading and trailing edges a the horizontal plane.
    stagger_angle_hub = 28.96 - 69.30  # degress
    stagger_angle_mid = 23.44 - 63.77  # degress
    stagger_angle_tip = 19.53 - 59.87  # degress

    # Do not change the following definitions
    rh = ri_hub
    rt = ri_tip
    ri_mid = (ri_hub + ri_tip) / 2

    # Generate the circular profile with user inputs creating 2D and 3D coordinates
    (
        x_3d_hub, y_3d_hub, z_3d_hub,
        x_hub_2d, y_hub_2d,
        beta1_hub, beta2_hub,
        L_hub, x1_hub, x2_hub, rc_hub, xc_hub, yc_hub, ca_hub, beta_avg_hub
    ) = math_model.generate_circular_profile(H, Q, N, eff, ri_hub, rh, rt, z, stagger_angle_hub)

    (
        x_3d_mid, y_3d_mid, z_3d_mid,
        x_mid_2d, y_mid_2d,
        beta1_mid, beta2_mid,
        L_mid, x1_mid, x2_mid, rc_mid, xc_mid, yc_mid, ca_mid, beta_avg_mid
    ) = math_model.generate_circular_profile(H, Q, N, eff, ri_mid, rh, rt, z, stagger_angle_mid)

    (
        x_3d_tip, y_3d_tip, z_3d_tip,
        x_tip_2d, y_tip_2d,
        beta1_tip, beta2_tip,
        L_tip, x1_tip, x2_tip, rc_tip, xc_tip, yc_tip, ca_tip, beta_avg_tip
    ) = math_model.generate_circular_profile(H, Q, N, eff, ri_tip, rh, rt, z, stagger_angle_tip)

    # Transform x,y,z to Bladegen meridional coordinates %m_prime, theta
    m_hub, theta_hub = math_model.cartesian_to_meridional(x_3d_hub, y_3d_hub, z_3d_hub)
    m_mid, theta_mid = math_model.cartesian_to_meridional(x_3d_mid, y_3d_mid, z_3d_mid)
    m_tip, theta_tip = math_model.cartesian_to_meridional(x_3d_tip, y_3d_tip, z_3d_tip)

    # Interpolate the meridional coordinates (%m_prime, theta) with 5th grade polynomial at specified points
    percentage_positions = [0.0, 0.25, 0.5, 0.75, 1.0]  # Points to evaluate the fitted fucntion theta(%m_prime)
    interpolated_values_hub, r_squared_hub = math_model.interpolate_polynomial(m_hub, theta_hub, percentage_positions)
    interpolated_values_mid, r_squared_mid = math_model.interpolate_polynomial(m_mid, theta_mid, percentage_positions)
    interpolated_values_tip, r_squared_tip = math_model.interpolate_polynomial(m_tip, theta_tip, percentage_positions)

    # Blade geometry calculations
    # Calculate the profile's cord legnth
    cord_length_hub_rotated = math_model.cord_length(x_hub_2d[0], x_hub_2d[-1], y_hub_2d[0], y_hub_2d[-1])
    cord_length_mid_rotated = math_model.cord_length(x_mid_2d[0], x_mid_2d[-1], y_mid_2d[0], y_mid_2d[-1])
    cord_length_tip_rotated = math_model.cord_length(x_tip_2d[0], x_tip_2d[-1], y_tip_2d[0], y_tip_2d[-1])

    # Calculate the stagger angle with arctan2 of the rotated circular profile
    rotated_stagger_angle_hub_arctan2 = math_model.stagger_angle_rotated_arctan2(x_hub_2d[0], x_hub_2d[-1], y_hub_2d[0], y_hub_2d[-1])
    rotated_stagger_angle_mid_arctan2 = math_model.stagger_angle_rotated_arctan2(x_mid_2d[0], x_mid_2d[-1], y_mid_2d[0], y_mid_2d[-1])
    rotated_stagger_angle_tip_arctan2 = math_model.stagger_angle_rotated_arctan2(x_tip_2d[0], x_tip_2d[-1], y_tip_2d[0], y_tip_2d[-1])

    # Calculate the rotated blade's axial cord
    ca_hub_rotated = math_model.axial_cord_rotated_blade(cord_length_hub_rotated, rotated_stagger_angle_hub_arctan2)
    ca_mid_rotated = math_model.axial_cord_rotated_blade(cord_length_mid_rotated, rotated_stagger_angle_mid_arctan2)
    ca_tip_rotated = math_model.axial_cord_rotated_blade(cord_length_tip_rotated, rotated_stagger_angle_tip_arctan2)

    # Print data
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
    print(f"beta1_hub: {beta1_hub:.2f} °, beta2_hub: {beta2_hub:.2f} °, L_hub: {L_hub*1000:.2f} mm, x1_hub: {x1_hub*1000:.2f} mm, x2_hub: {x2_hub*1000:.2f} mm, rc_hub: {rc_hub*1000:.2f} mm, xc_hub: {xc_hub*1000:.2f} mm, yc_hub: {yc_hub*1000:.2f} mm")
    print(f"beta1_mid: {beta1_mid:.2f} °, beta2_mid: {beta2_mid:.2f} °, L_mid: {L_mid*1000:.2f} mm, x1_mid: {x1_mid*1000:.2f} mm, x2_mid: {x2_mid*1000:.2f} mm, rc_mid: {rc_mid*1000:.2f} mm, xc_mid: {xc_mid*1000:.2f} mm, yc_mid: {yc_mid*1000:.2f} mm")
    print(f"beta1_tip: {beta1_tip:.2f} °, beta2_tip: {beta2_tip:.2f} °, L_tip: {L_tip*1000:.2f} mm, x1_tip: {x1_tip*1000:.2f} mm, x2_tip: {x2_tip*1000:.2f} mm, rc_tip: {rc_tip*1000:.2f} mm, xc_tip: {xc_tip*1000:.2f} mm, yc_tip: {yc_tip*1000:.2f} mm")
    print("---------------------------------")

    # Print R^2 for the interpolated values of theta(%m_prime)
    # print("Coefficient of determination R^2 for the interpolated values of %m vs theta")
    # print(f"R^2 at hub = {r_squared_hub}")
    # print(f"R^2 at mid = {r_squared_mid}")
    # print(f"R^2 at tip = {r_squared_tip}\n")
    # print("---------------------------------")

    # Printing interpolated values of the function theta(%m_prime) and R^2
    for i, pos in enumerate(percentage_positions):
        print(f"Interpolated value at %m_prime = {int(pos*100)}%")
        print(f"    theta_hub = {interpolated_values_hub[i]:.4f} °")
        print(f"    theta_mid = {interpolated_values_mid[i]:.4f} °")
        print(f"    theta_tip = {interpolated_values_tip[i]:.4f} °")

    print("---------------------------------")

    # Print the blade axial chord Ca/2 for bladegen in mm
    print("Original blade half axial chord (Ca/2) for Bladegen meridional view:")
    print(f"    hub: {(ca_hub/2)*1000:.4f} mm")
    print(f"    mid: {(ca_mid/2)*1000:.4f} mm")
    print(f"    tip: {(ca_tip/2)*1000:.4f} mm")

    print("Rotated blade half axial chord (Ca/2) for Bladegen meridional view:")
    print(f"    hub: {(ca_hub_rotated/2)*1000:.4f} mm")
    print(f"    mid: {(ca_mid_rotated/2)*1000:.4f} mm")
    print(f"    tip: {(ca_tip_rotated/2)*1000:.4f} mm")

    print("Rotated profile cord lengths (C):")
    print(f"    hub: {cord_length_hub_rotated*1000:.2f} mm")
    print(f"    mid: {cord_length_mid_rotated*1000:.2f} mm")
    print(f"    tip: {cord_length_tip_rotated*1000:.2f} mm")
    print("---------------------------------")

    print("Original blade stagger angles:")
    print(f"    hub: {beta_avg_hub:.2f} °")
    print(f"    mid: {beta_avg_mid:.2f} °")
    print(f"    tip: {beta_avg_tip:.2f} °")

    print("Rotated blade stagger angles (arctan2):")
    print(f"    hub: {rotated_stagger_angle_hub_arctan2:.2f} °")
    print(f"    mid: {rotated_stagger_angle_mid_arctan2:.2f} °")
    print(f"    tip: {rotated_stagger_angle_tip_arctan2:.2f} °")

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

    # # Plotting %m_prime vs. theta
    plt.plot(m_hub, theta_hub, label="Hub")
    plt.plot(m_mid, theta_mid, label="Mid")
    plt.plot(m_tip, theta_tip, label="Tip")
    plt.xticks(percentage_positions)
    # plt.xlim(0, 1)
    # plt.ylim(0,)
    plt.xlabel(r"$m'\; (-)$")
    plt.ylabel(r"$\theta\; (^\circ)$")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Export profile 2D and 3D data in mm for the thesis methodology
    # np.savetxt("x_2d_rotated_hub.txt", x_hub_2d * 1000, delimiter=",")
    # np.savetxt("y_2d_rotated_hub.txt", y_hub_2d * 1000, delimiter=",")

    # np.savetxt("x_2d_rotated_mid.txt", x_mid_2d * 1000, delimiter=",")
    # np.savetxt("y_2d_rotated_mid.txt", y_mid_2d * 1000, delimiter=",")

    # np.savetxt("x_2d_rotated_tip.txt", x_tip_2d * 1000, delimiter=",")
    # np.savetxt("y_2d_rotated_tip.txt", y_tip_2d * 1000, delimiter=",")

    # np.savetxt("x_3D_hub.txt", x_3d_hub * 1000, delimiter=",")
    # np.savetxt("y_3D_hub.txt", y_3d_hub * 1000, delimiter=",")
    # np.savetxt("z_3D_hub.txt", z_3d_hub * 1000, delimiter=",")

    # np.savetxt("x_3D_mid.txt", x_3d_mid * 1000, delimiter=",")
    # np.savetxt("y_3D_mid.txt", y_3d_mid * 1000, delimiter=",")
    # np.savetxt("z_3D_mid.txt", z_3d_mid * 1000, delimiter=",")

    # np.savetxt("x_3D_tip.txt", x_3d_tip * 1000, delimiter=",")
    # np.savetxt("y_3D_tip.txt", y_3d_tip * 1000, delimiter=",")
    # np.savetxt("z_3D_tip.txt", z_3d_tip * 1000, delimiter=",")

    # # Export meridional coordinates data for the thesis methodology
    # np.savetxt("m_hub.txt", m_hub, delimiter=",")
    # np.savetxt("m_mid.txt", m_mid, delimiter=",")
    # np.savetxt("m_tip.txt", m_tip, delimiter=",")

    # np.savetxt("theta_hub.txt", theta_hub, delimiter=",")
    # np.savetxt("theta_mid.txt", theta_mid, delimiter=",")
    # np.savetxt("theta_tip.txt", theta_tip, delimiter=",")
