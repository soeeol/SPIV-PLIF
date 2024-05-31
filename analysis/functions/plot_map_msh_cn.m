##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## plot cn avg with interface
##
## Author: Sören J. Gerke
##

function fh = plot_map_msh_cn (c_msh, cn_dyn_avg, cax_max, x, dX, y_wall, delta_u_avg, delta_u_fit_avg, ids_C)

  fh = plot_map_msh ({c_msh{1}+dX c_msh{2} c_msh{3}}, cn_dyn_avg, []);
  caxis ([0 cax_max])
  colorbar;
  hold on;
  draw_cell (ids_C, 0, 1);
  plot (x+dX, y_wall, "w;wall measured;");
  plot (x+dX, delta_u_avg, "m;delta_u measured;");
  plot (x+dX, delta_u_fit_avg, "r;delta_u spline fit;");
  xlabel ("x* in mm");
  ylabel ("y in mm");
  title ("cn avg");
  xlim ([min(x) max(x)] + dX);
  ylim ([0.9*min(delta_u_fit_avg) 1.1*max(delta_u_fit_avg)]);

endfunction
