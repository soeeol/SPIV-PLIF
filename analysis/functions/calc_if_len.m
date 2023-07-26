##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## interface length estimate
##
## Author: Sören J. Gerke
##

function len = calc_if_len (c_msh_x, h_if, lims_x)
  x = c_msh_x(1,:);
  idx = ((x<=lims_x(2)) & (x>=lims_x(1)));
  x = x(idx);
  ##y = c_h.gas(idx);
  if (iscell(h_if))
    y = h_if{1}(idx);
  else
    y = h_if(idx);
  endif
  ## line segments
  dx = x(2)-x(1);
  yn = y(2:end);
  ynm = y(1:end-1);
  dy = abs (yn - ynm);
  len = sqrt (dx.^2 + dy.^2);
  ##6 * sqrt (4^2 + 1^2)
  ##sum(len(!isnan(len)));
  ##sum(len);
  ##plot(x,y)
  ##hold on
  ##plot(dy)
endfunction
