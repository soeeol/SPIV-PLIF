##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function fh = plot_cn_if_norm (x, h_g, msh_n, c_msh, cn, id_cell, xpos)

  c_msh{1} = c_msh{1} + xpos;

  fh = plot_map_msh (c_msh, cn, []);
  hold on;
  draw_cell (id_cell, 0, true);
  for i_p = 1 : 10 : numel (h_g)
    plot (msh_n{1}(:,i_p) + xpos, msh_n{2}(:,i_p), "m");
  endfor
  plot (x + xpos, h_g, "r");

  axis image;

  xlim ([min(x) max(x)] + xpos);
  ylim ([0 1.1*max(h_g)]);

  xlabel ("x* in mm");
  ylabel ("y in mm");

endfunction
