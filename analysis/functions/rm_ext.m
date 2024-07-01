##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## remove outliers in sc_vec at the extremes of the domain defined by x_vec
##
## Author: Sören J. Gerke
##

function sc_vec_rm_ext = rm_ext (x_vec, sc_vec, wlen)

  sc_vec_rm_ext = outlier_rm (sc_vec, movmedian (sc_vec, wlen));

  sf_x = abs (x_vec(2) - x_vec(1));

  ## valid inside idx
  idx = (x_vec > min (x_vec) + sf_x * wlen / 2) & (x_vec < max (x_vec) - sf_x * wlen / 2);

  ## keep inside data unchanged
  sc_vec_rm_ext(idx) = sc_vec(idx);

  ## decide for lower and upper border
  idx_l = (x_vec < min (x_vec) + sf_x * wlen / 2);
  idx_u = (x_vec > max (x_vec) - sf_x * wlen / 2);
  for i = 1:2
    switch (i)
      case 1
        idx_border = idx_l;
      case 2
        idx_border = idx_u;
    endswitch
    sc_min = min (sc_vec(idx_border));
    sc_max = max (sc_vec(idx_border));
    sc_m = median (sc_vec(idx_border));
    sc_dmin = (sc_min - sc_m);
    sc_dmax = (sc_max - sc_m);
    sc_delta = sc_max - sc_min;
    if (! sc_dmin * sc_dmax == 0)
      if (abs (1 - abs (sc_dmin / sc_dmax)) < 0.5) # otherwise strongly asymmetric
            sc_vec_rm_ext(idx_border) = sc_vec(idx_border); # keep original data in border region
      endif
    endif
  endfor

##  figure ();
##  hold on;
##  plot (x_vec, sc_vec, "b");
##  plot(x_vec(idx_u), sc_vec(idx_u), "r");
##  plot(x_vec(idx_l), sc_vec(idx_l), "r");
##  plot (x_vec, sc_vec_rm_ext, "g");

endfunction

