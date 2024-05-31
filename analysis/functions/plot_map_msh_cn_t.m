##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## plot cn dyn with interface, dyn vs. avg
##
## Author: Sören J. Gerke
##

function plot_map_msh_cn_t (c_msh, cn_dyn, cax_max, x, delta_u_fit_avg, delta_u, delta_u_fit, i_t, fh)

  clf (fh);
  plot_map_msh (c_msh, cn_dyn, fh);
  caxis ([0 cax_max])
  colorbar;
  hold on;
  plot3 (x, delta_u_fit_avg, ones(1, numel (x)), "-k");
  plot3 (x, delta_u, ones(1, numel (x)), "-r");
  plot3 (x, delta_u_fit, ones(1, numel (x)), "-m");
  title (["cn dyn t#" num2str(i_t)]);
  xlabel ("x in mm");
  ylabel ("y in mm");
  xlim ([min(x) max(x)]);
  ylim ([0.9*min(delta_u_fit_avg) 1.1*max(delta_u_fit_avg)]);
  zlim ([-1 1]);

endfunction
