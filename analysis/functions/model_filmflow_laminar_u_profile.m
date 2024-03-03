##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## u (y) - velocity profile of Nusselt's solution
##
## input:
## y        .. profile coordinate in m
## u_s      .. velocity at y = delta_u in m / s
## delta_u  .. delta_u in m
## (get delta_u and u_s from model_filmflow_laminar_u_profile_p.m)
##
## output:
## u  .. velocity profile in m / s
##
## Author: Sören J. Gerke
##

function u = model_filmflow_laminar_u_profile (y, u_s, delta_u)

  u = u_s .* ((2 .* y ./ delta_u) - (y ./ delta_u) .^ 2); # \label{eq:filmflow_u_profile}

endfunction
