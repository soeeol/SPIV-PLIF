##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function x_off = calc_xsec_if_offset_x (x_l, x_r, delta_u_l, delta_u_r, xpos_l, xpos_r)

  delta_u_l = delta_u_l;
  delta_u_r = delta_u_r;

  x_l = x_l + xpos_l;
  x_r = x_r + xpos_r;

  ol_msh = abs (min(x_r) - max(x_l));

  [~, xidxl_1] = min (abs (x_l - (max (x_l) - ol_msh)));
  [~, xidxu_1] = min (abs (x_l - (max (x_l) + ol_msh)));
  [~, xidxl_2] = min (abs (x_r - (min (x_r) - ol_msh)));
  [~, xidxu_2] = min (abs (x_r - (min (x_r) + ol_msh)));

  x_1_o = x_l(xidxl_1:xidxu_1);
  x_2_o = x_r(xidxl_2:xidxu_2);
  opt1 = delta_u_l(xidxl_1:xidxu_1);
  opt2 = delta_u_r(xidxl_2:xidxu_2);

  mopt1 = min (opt1);
  mopt2 = max (opt2);
  midopt = mean ([mopt1 mopt2]);

  pts_l = intersectPolylines ([vec(x_1_o) vec(opt1)], [[min(x_1_o) max(x_1_o)]' midopt*[1 1]']);
  pts_r = intersectPolylines ([vec(x_2_o) vec(opt2)], [[min(x_2_o) max(x_2_o)]' midopt*[1 1]']);
  x_off = max (pts_l(:,1)) - min (pts_r(:,1));

  if isempty (x_off)
    x_off = 0;
  endif
  if (abs (mad (opt1) + mad (opt2)) < 1e-3)
    x_off = 0;
  endif

  figure ();
  hold on;
  plot (x_l, delta_u_l, "k")
  plot (x_r, delta_u_r, "r")
  plot (x_1_o, opt1, "xk")
  plot (x_2_o, opt2, "or")
  plot ([min(x_l) max(x_l)], mopt1*[1 1], "k--")
  plot ([min(x_r) max(x_r)], mopt2*[1 1], "r--")
  plot ([min(x_r) max(x_r)], ((mopt1 + mopt2) / 2 ) * [1 1], "b--")
  plot (x_l-x_off, delta_u_l, "g")

endfunction
