##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## u (y) - velocity profile of Nusselt's solution
##
## input:
## y    .. y in m
## u_s  .. u at y = h in m / s
## h    .. h in m
## (get h and u_s from model_filmflow_laminar_u_profile_p.m)
##
## output:
## u  .. velocity profile in m / s
##
## Author: Sören J. Gerke
##

function u = model_filmflow_laminar_u_profile (y, u_s, h)
  u = u_s .* ( (2 .* y ./ h) - (y ./ h) .^2 );
endfunction
