##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## mean of all cells, used to compute average of images loaded in cell
##
## Author: Sören J. Gerke
##

function out = calc_im_avg_cells (im_cell, method)

  n_c = numel (im_cell);

  dummy = [];
  for i_c = 1:n_c
    dummy(:,:,i_c) = im_cell{i_c};
  endfor

  switch (method)
    case {"mean"}
      out = mean (dummy, 3);
    case {"median"}
      out = median (dummy, 3);
    otherwise
      error ("calc_im_avg_cells: method unknown.");
  endswitch

endfunction
