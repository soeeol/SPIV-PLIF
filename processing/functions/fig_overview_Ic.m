##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## section processing result overview: fluorescence maps
##
## Author: Sören J. Gerke
##

function fh = fig_overview_Ic (pp, c_msh, c_dat, c_h)
  fp = fig_param ({"combo", "double", "elsevier"});
  domain = get_domain (pp);
  fig_size_x = fp.xsize; # cm
  fig_size_y = 8;
  fh = figure ("DefaultAxesFontSize", fp.fs_min, "DefaultAxesFontName", fp.fonts{1});
  set (gcf, "PaperPositionMode", "manual");
  set (gcf, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
  set (gcf, "PaperSize", [fig_size_x fig_size_y]);
  set (gcf, "renderer", "opengl");
  titles = {"Ic \n in a.u.", "Ic0 \n in a.u.", "Ic1 \n in a.u."};
  Icmax = 0;
  for i = 1:numel(c_dat)
   Icmax = max([Icmax, max(max(c_dat{i}))]);
  endfor
  for i = 1:numel(c_dat)
    subplot (2, 2, i)
    surf (c_msh{1}, c_msh{2}, c_msh{3}, c_dat{i}); shading flat; view ([0 0 1]);
    hold on
    plot (c_msh{1}(1,:), c_h.gas, "k");
    plot (c_msh{1}(1,:), c_h.wall, "k");
    axis image
    xlabel ("x in mm"); ylabel ("y in mm");
    yticks (0:1:max(max(c_msh{2}))), xticks (domain.xmin:2:domain.xmin)
    set (gca, "TickDir" , "out")
    set (gca, "XMinorTick" , "on" , "YMinorTick" , "on")
    set (gca, "XMinorTickValues", domain.xmin:1:domain.xmax);
    set (gca, "YMinorTickValues", 0:1:max(max(c_msh{2})));
    hax = colorbar ("location", "EastOutside"); colormap viridis
    title (hax, titles{i})
    caxis ([0 Icmax]);
    hold off
  endfor
endfunction
