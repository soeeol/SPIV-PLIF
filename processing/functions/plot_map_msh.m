##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function fh = plot_map_msh (msh, map)
  fh = figure ();
  sf = get_sf (msh);
  surf (msh{1}-sf(1)/2, msh{2}-sf(2)/2, msh{3}, map);
  view ([0 0 1]);
  shading flat;
  colormap viridis;
  xlabel ("y in mm");
  ylabel ("x in mm");
endfunction
