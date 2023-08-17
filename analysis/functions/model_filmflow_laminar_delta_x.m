##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## local boundary layer thickness
##
## input:
## x    .. downstream coordinate in m
## D    .. diffusivity in m^2 / s
## u_s  .. surface velocity in m / s
##
## output:
## delta  .. concentration layer thickness in m
##
## Author: Sören J. Gerke
##

function delta = model_filmflow_laminar_delta_x (x, D, u_s)
  delta = sqrt ( (pi * x * D) ./ u_s' );
endfunction
