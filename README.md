# Blagen coordinates generator

The programs available on this repo generate turbine blade coordinates needed in Bladegen of ANSYS.

The coordinates Bladegen uses are $$m, \; \theta$$, where $$m$$ is the meridional coordinate and $$\theta$$ is the blade wrap angle.

- The "circular_blade" program generates the Bladegen coordinates $$m,\; \theta$$ for a curved and constant thickness blade.
- The "naca_4_digit" program generates the Bladegen coordinates $$m,\; \theta$$ for a 4-digit NACA profile with curvature and variable thickness.
