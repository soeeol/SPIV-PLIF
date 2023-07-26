##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## local mass transfer coefficient beta in m / s
##
## input
## x    .. downstream coordinate in m
## D    .. diffusivity in m^2 / s
## u_s  .. surface velocity in m / s
##
## Author: Sören J. Gerke
##

function beta = model_filmflow_laminar_beta_x (x, D, u_s)
  beta = sqrt ( (u_s' * D) ./ (pi * x) );
endfunction
