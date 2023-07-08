##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## molar fraction from mass fraction mol / mol
##
## Author: Sören J. Gerke
##

function mx = fp_mx_mf (mf, mm_A, mm_AB)
  mx = mf ./ mm_A .* mm_AB;
endfunction
