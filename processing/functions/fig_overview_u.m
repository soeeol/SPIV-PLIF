##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## section processing result overview: velocity field
##
## Author: Sören J. Gerke
##

function fh = fig_overview_u (pp, u_msh, u_dat, u_masks, c_msh, c_h, c_masks)

  titles = {"ux \n in m/s", "uy \n in m/s", "uz \n in m/s", "um \n in m/s"};

  fh = figure ();
  for i = 1:4

    subplot (4, 1, i);

    if !(i==4)
      mask = 1;
    else
      mask = u_masks.wall .* u_masks.gas;
    endif

    plot_map_msh (u_msh, u_dat{i} .* mask, fh);

    caxis auto
    u_max = max (max (abs (u_dat{i} .* u_masks.wall .* u_masks.gas)));
    switch (i)
      case 2 # y
        caxis ([-u_max u_max])
      case 3 # z
        caxis ([-u_max u_max])
      case 4 # mag
        caxis ([0 u_max])
    endswitch

    hold on;

    if (i==4) # some streamlines
      xpos = 0;
      skip = 2;
      ymax = max (max (u_msh{2}));
      NX = 1;
      NY = round (1 / skip * ymax / (u_msh{1}(1,2) - u_msh{1}(1,1)));
      SY = reshape (linspace (0, ymax, NY)' .* ones (1, NX), NX*NY, 1);
      SX = reshape (linspace (xpos, xpos, NX) .* ones (NY,1 ), NX*NY, 1);
      hstream2 = streamline (u_msh{1}, u_msh{2}, u_dat{1}.*u_masks.gas, u_dat{2}.*mask, SX, SY, [1 5e2]);
      set (hstream2, "Color", "w");
      hstream2 = streamline (u_msh{1}, u_msh{2}, u_dat{1}.*u_masks.gas, u_dat{2}.*mask, SX, SY, [-1 5e2]);
      set (hstream2, "Color", "w");
    endif

    plot (c_msh{1}(1,:), c_h.gas, "k");
    plot (c_msh{1}(1,:), c_h.wall, "k");

    axis image
    ylabel ("y in mm");
    hax = colorbar ("location", "EastOutside");
    title (hax, titles{i})
    hold off

  endfor
  xlabel ("x in mm");

endfunction
