##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## section processing result overview: velocity field
##
## Author: Sören J. Gerke
##

function fh = fig_overview_u (pp, u_msh, u_dat, u_masks, c_msh, c_h, c_masks)
  fp = fig_param ({"combo", "double", "elsevier"});
  domain = get_domain (pp);
  fig_size_x = fp.xsize; # cm
  fig_size_y = 8;
  fh = figure ("DefaultAxesFontSize", fp.fs_min, "DefaultAxesFontName", fp.fonts{1});
  set (gcf, "PaperPositionMode", "manual");
  set (gcf, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
  set (gcf, "PaperSize", [fig_size_x fig_size_y]);
  set (gcf, "renderer", "opengl");
  titles = {"ux \n in m/s", "uy \n in m/s", "uz \n in m/s", "um \n in m/s"};
  for i = 1:4
    subplot (2, 2, i); hold on;
    if !(i==4)
      mask = 1;
    else
      mask = u_masks.wall .* u_masks.gas;
    endif
    surf (u_msh{1}, u_msh{2}, u_msh{3}, u_dat{i}.*mask); shading flat; view ([0 0 1]);
    caxis auto
    u_max = max(max(abs(u_dat{i} .* u_masks.wall .* u_masks.gas)));
    switch (i)
      case 2 # y
        caxis ([-u_max u_max])
      case 3 # z
        caxis ([-u_max u_max])
      case 4 # mag
        caxis ([0 u_max])
    endswitch
    hold on
    if (i==4) # some streamlines
      xpos = 0; skip = 2;
      ymax = max(max(u_msh{2}));#1;
      NX = 1; NY = round( 1 / (skip) * ymax / (u_msh{1}(1,2)-u_msh{1}(1,1)));
      SY = reshape (linspace(0,ymax,NY)' .* ones (1,NX),NX*NY,1);
      SX = reshape (linspace(xpos,xpos,NX) .* ones (NY,1),NX*NY,1);
      hstream2 = streamline (u_msh{1}, u_msh{2}, u_dat{1}.*u_masks.gas, u_dat{2}.*mask, SX, SY,[1 5e2]);
      set (hstream2, "Color", "w");
      hstream2 = streamline (u_msh{1}, u_msh{2}, u_dat{1}.*u_masks.gas, u_dat{2}.*mask, SX, SY,[-1 5e2]);
      set (hstream2, "Color", "w");
    endif
    plot (c_msh{1}(1,:), c_h.gas, "k");
    plot (c_msh{1}(1,:), c_h.wall, "k");
    axis image
    xlabel ("x in mm"); ylabel ("y in mm");
    yticks (0:1:max(max(c_msh{2}))), xticks (domain.xmin:2:domain.xmax)
    set (gca, "TickDir" , "out")
    set (gca, "XMinorTick" , "on" , "YMinorTick" , "on")
    set (gca, "XMinorTickValues", domain.xmin:1:domain.xmax);
    set (gca, "YMinorTickValues", 0:1:max(max(c_msh{2})));
    hax = colorbar ("location", "EastOutside"); colormap viridis
    title (hax, titles{i})
    hold off
  endfor
endfunction
