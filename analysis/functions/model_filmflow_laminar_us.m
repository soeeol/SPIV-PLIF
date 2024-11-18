##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## surface velocity of film
##
## input:
## nu     .. kinematic viscosity in m^2 / s
## alpha  .. plate inclination vs. the horizontal in rad
## Re     .. inlet related Reynolds number (nd_re_inlet.m)
##
## Author: Sören J. Gerke
##

function u_s = model_filmflow_laminar_us (nu, alpha, Re)

  g = 9.81; # m / s^2

  u_s = 1 / 2 * (9 * nu * Re .^ 2 * g * sin (alpha)) .^ (1 / 3); # \label{eq:filmflow_surface_velocity}

endfunction
