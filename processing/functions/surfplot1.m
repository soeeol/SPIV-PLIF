##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function surfplot1 (msh, scmap, param)

  surf (msh{1}, msh{2}, msh{3}, scmap);
  shading flat;
  colormap viridis;
  view ([0 0 1]);
  caxis ("auto");
  axis image;

endfunction
