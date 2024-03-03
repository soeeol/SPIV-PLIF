##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## analytical nondimensional concentration profile
##
## Author: Sören J. Gerke
##

function c_nd = model_filmflow_laminar_c_profile_nd (y_nd)

  c_nd = erfc (y_nd * sqrt (pi / 4));

endfunction
