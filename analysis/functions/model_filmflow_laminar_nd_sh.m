##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Average Sherwood number related to film height for laminar film flow
## (delta_c << delta_u)
##
## input:
## x_nd  .. entrance length (nondimensional)
##
## Author: Sören J. Gerke
##

function Sh_h = model_filmflow_laminar_nd_sh (x_nd)
  Sh_h = sqrt ( 6 ./ (pi .* x_nd) );
endfunction
