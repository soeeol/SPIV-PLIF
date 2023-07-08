##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## molar mass binary mixture from mass fraction
##
## Author: Sören J. Gerke
##

function mm = fp_mm_mf (mf_A, mm_A, mm_B)
   mm = 1 ./ (mf_A./mm_A + (1-mf_A)./mm_B); # kg / mol
endfunction

