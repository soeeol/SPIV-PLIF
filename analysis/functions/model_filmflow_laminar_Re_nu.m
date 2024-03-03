##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## input:
## u_s      .. surface velocity at y = delta_u in m / s
## delta_u  .. film height in m
## alpha    .. plate inclination vs. the horizontal in rad
##
## output:
## nu     .. kinematic viscosity in m^2 / s
## re     .. film Reynolds number (nd_re_inlet.m)
##

##
## Author: Sören J. Gerke
##

function [re, nu] = model_filmflow_laminar_Re_nu (u_s, delta_u, alpha)

  g = 9.81;

  nu = g * sin (alpha) .* delta_u .^ 2 ./ (2 * u_s);

  re = 4 / 3 * u_s .^ 2 ./ (g * sin (alpha) .* delta_u);

endfunction
