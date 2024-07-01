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
## delta_c  .. concentration layer thickness in m
##
## Author: Sören J. Gerke
##

function delta_c = model_filmflow_laminar_deltac (x, D, u_s)

  x = vec (x);
  u_s = (vec (u_s))';
##  D = D(1);

  delta_c = sqrt ( (pi * x * D) ./ u_s ); # \label{eq:filmflow_deltac}

endfunction
