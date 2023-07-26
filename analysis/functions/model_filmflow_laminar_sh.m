##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## avg. Sherwood number laminar film delta_c << delta_u
##
## input:
## D    .. diffusivity in m^2 / s
## u_s  .. surface velocity m / s
## L    .. length inlet to outlet downstream in m
## h    .. film height in m
##
## Author: Sören J. Gerke
##

function Sh = model_filmflow_laminar_sh (u_s, D, L, h)
  Sh = sqrt ( (4 .* u_s .* h.^2) ./ (pi .* D .* L) );
endfunction
