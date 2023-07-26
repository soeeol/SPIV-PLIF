##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## trivial static hold-up estimate of structure
##
## Author: Sören J. Gerke
##

function HU = hu_liq_stat (id_cell, alpha)
  switch id_cell
    case {"2d-r10"} #, "2d-r10-20", "2d-r10-40", "2d-r10-60"}
      hms = 1; # mm
      wls = tan (alpha) * hms;
      HU = 0.5 * wls * hms; # mm^2
##    case {"2d-t10"}
##      type = "t10";
##    case {"2d-s10"}
##      type = "s10";
    case {"flat"}
      HU = 0.0;
  endswitch
endfunction
