##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## mass fraction from molar fraction
##
## Author: Sören J. Gerke
##

function mf = fp_mf_mx (mx, mm_A, mm_B)
  mf = mx .* mm_A ./ (mx .* mm_A + (1 - mx) .* mm_B);
endfunction
