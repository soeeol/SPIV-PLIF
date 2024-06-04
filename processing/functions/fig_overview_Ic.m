##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## section processing result overview: fluorescence maps
##
## Author: Sören J. Gerke
##

function fh = fig_overview_Ic (pp, c_msh, c_dat, c_h)

  nmap = numel (c_dat);

  titles = {"Ic \n in a.u.", "Ic0 \n in a.u.", "Ic1 \n in a.u."};

  Icmax = 0;
  for i = 1:nmap
   Icmax = max ([Icmax, max(max(c_dat{i}))]);
  endfor

  fh = figure ();
  for i = 1:nmap
    subplot (nmap, 1, i)
    plot_map_msh (c_msh, c_dat{i}, fh);
    caxis ([0 Icmax]);
    hold on
    plot (c_msh{1}(1,:), c_h.gas, "k");
    plot (c_msh{1}(1,:), c_h.wall, "k");
    axis image
    ylabel ("y in mm");
    hax = colorbar ("location", "eastoutside");
    title (hax, titles{i})
    hold off;
  endfor
  xlabel ("x in mm");

endfunction
