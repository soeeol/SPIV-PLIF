##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## wall indicator map for wall threshold method
##
## Author: Sören J. Gerke
##

function map_out = ind_wall_c (map_in)

  map_out = map_in / max (max (map_in));

endfunction
