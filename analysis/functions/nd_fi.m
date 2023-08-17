##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## film number
##
## input:
## rho    .. density in kg / m^3
## sigma  .. surface tension in N / m
## eta    .. dyn. viscosity in Pa * s
##
## Author: Sören J. Gerke
##

function Fi = nd_fi (rho, sigma, eta)
  g = 9.81; # m / s^2
  Fi = (rho * sigma^3) / (eta^4 * g);
endfunction

