##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## PDMS refractive index vs. temperatue for the matching experiment
##
## Author: Sören J. Gerke
##

function ri = ri_match_PDMS_ri_T (T_get)
  p_RI_PDMS_T_fit = [-3.2367e-04 1.5061e+00];
  ri = polyval (p_RI_PDMS_T_fit, T_get);
endfunction
