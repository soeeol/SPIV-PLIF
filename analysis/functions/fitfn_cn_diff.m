##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fit function diffusivity measurements
##
## p_1, p_2 .. fit parameters
## sn       .. coordinate
## p_sc     .. scaling factor of denominator aiding the numerical optimization
##
## Author: Sören J. Gerke
##

function cout = fitfn_cn_diff (p, sn, p_sc)
   cout = p(2) * erfc (sn / (p(1) * p_sc));
endfunction
