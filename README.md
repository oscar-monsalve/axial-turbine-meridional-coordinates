# Blagen coordinates generator

The codes available in this repository generate turbine blade meridional coordinates. These coordinates are $m, \\; \theta$, where $m$ is
the meridional coordinate and $\theta$ is the blade wrap angle.

1. The "circular_blade.py" file generates meridional coordinates $m,\\; \theta$ for a curved and constant thickness blade.
2. The "naca_4_digit.py" file generates meridional coordinates $m,\\; \theta$ for a 4-digit NACA profile with curvature and variable thickness.

Note: this code is the first implementation iteration, meaning it can be improved.

The input parameters for both programs are:

- Head pressure $H\\; (m)$.
- Flow rate $Q\\; (m^3/s)$.
- Turbine angular velocity $N\\; (min^{-1})$.
- Turbine assumed hydraulic efficiency $\eta_t\\; (-)$.
- Turbine hub radius $r_{hub}\\; (m)$.
- Turbine tip radius $r_{tip}\\; (m)$.
- Blade number $z\\; (-)$.
