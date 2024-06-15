import matplotlib.pyplot as plt
import math_model


def main() -> None:

    # User inputs:
    # H: head (m)
    # Q: flow rate (m^3/s)
    # N: angular velocity (rpm)
    # eff: turbine efficiency (-)
    # r_hub: hub radius (m)
    # r_tip: tip radius (m)
    # z: number of blades (-)

    # The following hydraulic parameters can be modified to change the circular blade geometry design:
    H = 3  # m
    Q = 0.015  # m^3\s
    N = 3600  # rpm
    eff = 0.65
    ri_hub = 0.02259  # m
    ri_tip = 0.03765  # m
    z = 5

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
    ) = math_model.generate_circular_profile(H, Q, N, eff, ri_hub, rh, rt, z)

    (
        x_3d_mid, y_3d_mid, z_3d_mid,
        x_mid_2d, y_mid_2d,
        beta1_mid, beta2_mid,
        L_mid, x1_mid, x2_mid, rc_mid, xc_mid, yc_mid, ca_mid, beta_avg_mid
    ) = math_model.generate_circular_profile(H, Q, N, eff, ri_mid, rh, rt, z)

    (
        x_3d_tip, y_3d_tip, z_3d_tip,
        x_tip_2d, y_tip_2d,
        beta1_tip, beta2_tip,
        L_tip, x1_tip, x2_tip, rc_tip, xc_tip, yc_tip, ca_tip, beta_avg_tip
    ) = math_model.generate_circular_profile(H, Q, N, eff, ri_tip, rh, rt, z)

    # Transform x,y,z to Bladegen meridional coordinates %m_prime, theta
    m_hub, theta_hub = math_model.cartesian_to_meridional(x_3d_hub, y_3d_hub, z_3d_hub)
    m_mid, theta_mid = math_model.cartesian_to_meridional(x_3d_mid, y_3d_mid, z_3d_mid)
    m_tip, theta_tip = math_model.cartesian_to_meridional(x_3d_tip, y_3d_tip, z_3d_tip)

    # Interpolate the meridional coordinates (%m_prime, theta) with 5th grade polynomial at specified points
    percentage_positions = [0.0, 0.25, 0.5, 0.75, 1.0]  # Points to evaluate the fitted fucntion theta(%m_prime)
    interpolated_values_hub, r_squared_hub = math_model.interpolate_polynomial(m_hub, theta_hub, percentage_positions)
    interpolated_values_mid, r_squared_mid = math_model.interpolate_polynomial(m_mid, theta_mid, percentage_positions)
    interpolated_values_tip, r_squared_tip = math_model.interpolate_polynomial(m_tip, theta_tip, percentage_positions)

    # Calculate the profile's cord legnth
    cord_length_hub_original = math_model.cord_length(x_hub_2d[0], x_hub_2d[-1], y_hub_2d[0], y_hub_2d[-1])
    cord_length_mid_original = math_model.cord_length(x_mid_2d[0], x_mid_2d[-1], y_mid_2d[0], y_mid_2d[-1])
    cord_length_tip_original = math_model.cord_length(x_tip_2d[0], x_tip_2d[-1], y_tip_2d[0], y_tip_2d[-1])

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
    print(f"""beta1_hub : {beta1_hub:.6f} °, beta2_hub : {beta2_hub:.6f} °, L_hub : {L_hub:.6f}, x1_hub : {x1_hub:.6f}, x2_hub : {x2_hub:.6f}, rc_hub : {rc_hub:.6f}, xc_hub : {xc_hub:.6f}, yc_hub : {yc_hub:.6f}""")
    print(f"""beta1_mid : {beta1_mid:.6f} °, beta2_mid : {beta2_mid:.6f} °, L_mid : {L_mid:.6f}, x1_mid : {x1_mid:.6f}, x2_mid : {x2_mid:.6f}, rc_mid : {rc_mid:.6f}, xc_mid : {xc_mid:.6f}, yc_mid : {yc_mid:.6f}""")
    print(f"""beta1_tip : {beta1_tip:.6f} °, beta2_tip : {beta2_tip:.6f} °, L_tip : {L_tip:.6f}, x1_tip : {x1_tip:.6f}, x2_tip : {x2_tip:.6f}, rc_tip : {rc_tip:.6f}, xc_tip : {xc_tip:.6f}, yc_tip : {yc_tip:.6f}\n""")
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

    print("Original profile cord lengths:")
    print(f"    hub -> {cord_length_hub_original:.5f} m")
    print(f"    mid -> {cord_length_mid_original:.5f} m")
    print(f"    tip -> {cord_length_tip_original:.5f} m")

    print("Stagger angles:")
    print(f"    hub -> {beta_avg_hub:.2f}°")
    print(f"    mid -> {beta_avg_mid:.2f}°")
    print(f"    tip -> {beta_avg_tip:.2f}°")

    # 2D plotting of unrotated profile
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

    # Plotting %m_prime vs. theta
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

    # # Export meridional coordinates data for the thesis methodology
    # np.savetxt("m_hub.txt", m_hub, delimiter=",")
    # np.savetxt("m_mid.txt", m_mid, delimiter=",")
    # np.savetxt("m_tip.txt", m_tip, delimiter=",")

    # np.savetxt("theta_hub.txt", theta_hub, delimiter=",")
    # np.savetxt("theta_mid.txt", theta_mid, delimiter=",")
    # np.savetxt("theta_tip.txt", theta_tip, delimiter=",")


if __name__ == "__main__":
    main()
