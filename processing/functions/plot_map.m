##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function fh = plot_map (map)

  fh = figure ();
  surf (map);
  view ([0 0 1]);
  shading flat;
  colormap viridis;

endfunction
