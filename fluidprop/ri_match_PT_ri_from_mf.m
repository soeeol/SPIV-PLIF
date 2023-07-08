##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## refractive index from mass fraction for the matching experiment
##
## water - glycerol
##
## Author: Sören J. Gerke
##

function ri =  ri_match_PT_ri_from_mf (mf_get, T_get)
  p_RI_PT = [
    -1.6929e-04
     1.9505e-01
    -1.0479e-04
     1.3568e+00];
  ri = mf_get .* (p_RI_PT(1).*T_get + p_RI_PT(2)) + (p_RI_PT(3).*T_get + p_RI_PT(4));
endfunction
