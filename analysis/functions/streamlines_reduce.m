##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## filter stream lines generated with
##
## Author: Sören J. Gerke
##

function [x_new y_new] = streamlines_reduce (line_x, line_y, new_x, new_y, skip, minLength)

  ## Reduce resolution and omit lines with low number of points
  keep_indices = 1:skip:length(line_x);
  if length (keep_indices) < minLength
    keep_indices = [];
  endif
  ## Add reduced line to new vectors
  x_new = [new_x; line_x(keep_indices); NaN];
  y_new = [new_y; line_y(keep_indices); NaN];

endfunction
