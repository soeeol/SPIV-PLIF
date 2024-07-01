##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## ideal wall contour based on the cell identifier
##
## Author: Sören J. Gerke
##

function y_wall = xy_wall_cell (id_cell, x)

  y_wall = zeros (size (x));

  switch (id_cell)
    case {"2d-r10"}
      id_ms = "r10";
      x_ms_l = - 1; # mm
      x_ms_u = + 1; # mm
      xoff = 0;
  endswitch

  for i_ms = 1 : numel (xoff)
    idx_ms = (x > x_ms_l + xoff(i_ms)) & (x < x_ms_u + xoff(i_ms));
  endfor

  switch (id_ms)
    case {"r10"}
      y_wall(idx_ms) = 1;
  endswitch

endfunction

