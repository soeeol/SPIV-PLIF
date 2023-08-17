##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## parameters for Nusselt velocity profile (model_filmflow_laminar_u_profile.m)
##
## input:
## nu     .. kinematic viscosity in m^2 / s
## alpha  .. plate inclination vs. the horizontal in rad
## Re     .. Reynolds number (nd_re_inlet.m)
##
## output:
## h      .. film height in m
## u_s    .. surface velocity at y = h in m / s
## u_avg  .. average velocity in m / s
##
## Author: Sören J. Gerke

function [h, u_s, u_avg] = model_filmflow_laminar_u_profile_p (nu, alpha, Re)
  g = 9.81; # m / s^2
  h = (3 * nu ^ 2 / (g * sin(alpha)) * Re) .^ (1/3);
  u_s = g * sin(alpha) * h.^2 / (2 * nu);
  u_avg = 2 / 3 * u_s;
endfunction
