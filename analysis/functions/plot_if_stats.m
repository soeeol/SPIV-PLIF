##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## plot interface statistics
##
## Author: Sören J. Gerke
##

function fh = plot_if_stats (if_stats, x, dom)

  in_dom_x = (x < dom.xmax) & (x > dom.xmin);

  fh = figure ();
  hold on;
  plot (x, if_stats.max_dev, "k;max abs deviation;");
  plot (x, if_stats.std_dev, "r;std deviation;");
  plot (x, if_stats.mad_dev, "b;mean abs deviation;");
  lh = legend ("autoupdate", "off");
  plot ([dom.xmin dom.xmax], [1 1] .* if_stats.mmax_dev, "k-.");
  plot ([dom.xmin dom.xmax], [1 1] .* if_stats.mmad_dev, "b-.");
  plot ([dom.xmin dom.xmax], [1 1] .* if_stats.mstd_dev, "r-.");
  xlabel ("x in mm");
  ylabel ("dispersion measures in mm");
  title ("deviation from mean interface position");

endfunction
