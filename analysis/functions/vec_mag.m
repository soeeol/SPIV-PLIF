##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## magnitude of vector field
##
## Author: Sören J. Gerke
##

function um = vec_mag (ux, uy)
  um = reshape (norm ([ux(:) uy(:)], "rows"), rows(ux), columns(ux));
endfunction
