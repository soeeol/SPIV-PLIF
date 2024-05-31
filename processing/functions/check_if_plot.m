##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function fh = check_if_plot (xy_if, ispeak, c_msh, c_dat)

  fh = plot_map_msh (c_msh, c_dat, []);
  title ("interface detection result");
  hold on;
  ## detected if peak
  plot3 ((xy_if(ispeak==1,1)), (xy_if(ispeak==1,2)), ones(1,numel(xy_if(ispeak==1,2))), "-m");
  ## no peak detected
  plot3 ((xy_if(ispeak==0,1)), (xy_if(ispeak==0,2)), ones(1,numel(xy_if(ispeak==0,2))), "oc");
  grid off;
  hold off;
  legend (["Ic map"], "if peak detected", "no if peak found", "interpreter", "none");

endfunction
