##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## interface indicator for fluorescence map with quenched interface
##
## Author: Sören J. Gerke
##

function map_out = ind_if (map_in)

  maxc = max (max (imsmooth (map_in, "Gaussian", 8)));

  map_out = (abs (map_in - maxc) / maxc);

endfunction
