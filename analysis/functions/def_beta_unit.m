##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## integral liquid side mass transfer coefficient of a transfer unit
##
## input:
## vfr    .. volumetric flow rate in m^3 / s
## A      .. effective mass transfer area in m^2
## c_eq   .. equilibrium concentration
## c_in   .. inlet concentration
## c_out  .. outlet concentration
##
## output:
## beta_c  .. in m / s
##
## Author: Sören J. Gerke
##

function beta_c = def_beta_unit (vfr, A, c_eq, c_in, c_out)

  beta_c = vfr ./ A .* log ( (c_eq - c_in) ./ (c_eq - c_out) ); # natural logarithm

endfunction
