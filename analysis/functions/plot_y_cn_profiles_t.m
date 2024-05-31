##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## test plot cn (y), single frame vs. avg
##
## Author: Sören J. Gerke
##

function plot_y_cn_profiles_t (fh, y, cn_dyn_avg, cn_dyn, delta_u_fit_avg, delta_u, x_idx, i_t)

  clf (fh)
  hold on;
  plot (y, cn_dyn_avg(:,x_idx), "k;cn dyn mean;");
  plot (y, cn_dyn(:,x_idx), "b;cn dyn;");
  plot (delta_u_fit_avg(x_idx)*[1 1], [0 1], "k--;delta_u dyn mean;");
  plot (delta_u(x_idx)*[1 1], [0 1], "b--;delta_u dyn;");
  xlim ([delta_u_fit_avg(x_idx)-0.5 1.1*delta_u_fit_avg(x_idx)]);
  xlabel ("y in mm");
  ylabel ("cn dyn in -");
  title (["cn dyn t#" num2str(i_t) " at pixel row " num2str(x_idx)]);
  legend ("location", "northwest");

endfunction

