##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## interface length estimate
## - alpha .. inclination of line segment to horizontal in radians
##
## Author: Sören J. Gerke
##

function [l_seg, l_total, incl_seg, t_s] = calc_if_len (x, delta_u)

  l_seg = incl_seg = zeros (1, length(x));

  x_in = vec (x);
  y_in = vec (delta_u);

  dx = abs (x_in(2:end) - x_in(1:end-1));
  dy = abs (y_in(2:end) - y_in(1:end-1));

  ## line segments min distance
  l_seg(2:end) = sqrt (dx.^2 + dy.^2);

  ## interface tangent coordinate
  t_s = [];
  for i = 1 : length (x)
    t_s(i) = sum (l_seg(1:i));
  endfor

  ## total length
  l_total = sum (l_seg);

  ## local angle of inclination
  alpha = angle2Points ([x_in(2:end) y_in(2:end)], [x_in(1:end-1) y_in(1:end-1)]);
  incl_seg(2:end) = alpha - pi; # inclination to horizontal in radians

endfunction
