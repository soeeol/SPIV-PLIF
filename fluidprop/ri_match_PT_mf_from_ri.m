##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## mass fraction from refractive index matching experiment
##
## water - glycerol
##
## Author: Sören J. Gerke
##

function mf =  ri_match_PT_mf_from_ri (ri_get, T_get)
  p_RI_PT = [
    -1.6929e-04
     1.9505e-01
    -1.0479e-04
     1.3568e+00];
  mf = (ri_get - (p_RI_PT(3).*T_get + p_RI_PT(4))) ./ (p_RI_PT(1).*T_get + p_RI_PT(2));
endfunction
