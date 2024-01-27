##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## PDMS refractive index vs. temperatue for the matching experiment
##
## Author: Sören J. Gerke
##

function ri = ri_PDMS_T (T, p)
  ri = polyval (p, T);
endfunction
