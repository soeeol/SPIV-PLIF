##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## dimension 1 and 2 scaling factors from equidistant mesh
##
## Author: Sören J. Gerke
##

function [sf] = get_sf (msh)
  sf = [];
  if !(isempty(msh))
    sf = [abs(msh{1}(1,1) - msh{1}(1,2));
          abs(msh{2}(1,1) - msh{2}(2,1))];
    sf = 1e-3 * round (1e3 * sf); # rounding for nice mm/px value
  endif
endfunction
