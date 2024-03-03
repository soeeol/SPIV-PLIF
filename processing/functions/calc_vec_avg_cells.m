##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## mean of all cells, used to compute average of vectors stored in cell
##
## Author: Sören J. Gerke
##

function [vec_mean, vec_std] = calc_vec_avg_cells (vec_cell, method)

  n_c = numel (vec_cell);
  n_p = numel (vec_cell{1});

  vec_mat = reshape (cell2mat (vec_cell), n_p, n_c);
  switch (method)
    case {"mean"}
      vec_mean = mean (vec_mat, 2);
    case {"median"}
      vec_mean = median (vec_mat, 2);
    otherwise
      error ("calc_vec_avg_cells: method unknown.");
  endswitch

  vec_std = std (vec_mat, [], 2);

endfunction
