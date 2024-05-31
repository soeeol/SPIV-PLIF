##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## interface statistics along x
##
## Author: Sören J. Gerke
##

function if_stats = calc_if_stats (delta_u, in_dom_x)

  n_t = numel (delta_u);

  if_vals = zeros (numel (delta_u{1}), n_t);
  for i_t = 1:n_t
    delta_u_dom = delta_u{i_t};
    for i_x = 1 : numel (delta_u_dom)
      if_vals(i_x,i_t) = delta_u_dom(i_x);
    endfor
  endfor

  if_stats.max_dev = max (abs (if_vals - median (if_vals,  2)), [], 2);
  if_stats.mad_dev = mad (if_vals, 1, 2); # 0:mean 1:median
  if_stats.std_dev = std (if_vals, [], 2);
  if_stats.mmax_dev = median (if_stats.max_dev(in_dom_x));
  if_stats.mmad_dev = median (if_stats.mad_dev(in_dom_x));
  if_stats.mstd_dev = median (if_stats.std_dev(in_dom_x));

endfunction
