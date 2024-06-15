# Meridional coordinates generator

The codes available in this repository generate axial turbine blades' meridional coordinates. These coordinates are $m, \\; \theta$, where $m$
is the meridional coordinate and $\theta$ is the blade wrap angle.

1. In directory "/circular_blade" is the code for a curved and constant thickness turbine blade.
2. In directory "/circular_blade_guide_vane" is the code for a curved and constant thickness guide vane. It generates the same circular profile
   for a desired stagger angle (blade inclination) at the hub, mid, and tip radii.
3. In directory "naca_4_digit" is the code for a 4-digit NACA profile turbine.

The input parameters for all the programs are:

- Head pressure $H\\; (m)$.
- Flow rate $Q\\; (m^3/s)$.
- Turbine angular velocity $N\\; (min^{-1})$.
- Turbine assumed hydraulic efficiency $\eta_t\\; (-)$.
- Turbine hub radius $r_{hub}\\; (m)$.
- Turbine tip radius $r_{tip}\\; (m)$.
- Blade number $z\\; (-)$.
- Only for the code in "/circular_blade_guide_vane", the stagger angles at the hub, mid, and tip radii have to be provided.

Note: this code is the first implementation iteration, meaning it can be improved.
