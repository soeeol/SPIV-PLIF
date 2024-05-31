##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## plot some phi_y_profiles to check alignment of phi_sat vs. phi/phi_des
##
## Author: Sören J. Gerke
##
function fh = plot_y_phi_profiles_test (y, phi_avg, phi_des, phi_sat, n_x, delta_u)

  d_nx = round (n_x / 4);

  fh = figure ();
  for i_x = d_nx:d_nx:n_x-d_nx
    hold on;
    plot (y, phi_avg(:,i_x), "k;phi;")
    plot (y, phi_des(:,i_x), "g;phi des;")
    plot (y, phi_sat(:,i_x), "r;phi sat;")
    legend ("autoupdate", "off")
    plot ([1 1]*delta_u(i_x), max(max(phi_avg))*[0 1], "--b")
    plot (y, 10 * 4 * del2 (imsmooth (phi_avg(:,i_x), 2)), "m;del2;")
  endfor
  plot ([0 0], max(max(phi_avg))*[0 1], "-b;wall;")
  xlim ([-0.6 max(y)])
  xlabel ("y in mm")
  ylabel ("intensity in a.u.")
  legend ("location", "northwest")

endfunction
