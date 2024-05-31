##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## plot avg fluorescence maps
##
## Author: Sören J. Gerke
##

function fh = plot_avg_phi_maps (c_msh, phi_avg, phi_des, phi_sat)

  sf = get_sf (c_msh);

  titles = {"phi avg\n in a.u.", "phi des\n in a.u.", "phi sat\n in a.u."};

  fh = figure ();
  for i = 1:3
    subplot (3, 1, i);
    switch (i)
      case 1
        plot_map = phi_avg;
      case 2
        plot_map = phi_des;
      case 3
        plot_map = phi_sat;
    endswitch
    plot_map_msh (c_msh, plot_map, fh);
    caxis ([0 max(max(phi_avg))]);
    grid off;
    axis image
    xlabel ("x in mm");
    ylabel ("y in mm");
    hax = colorbar ("location", "EastOutside");
    title (hax, titles{i})
  endfor

endfunction
