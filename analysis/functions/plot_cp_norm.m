##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function fh = plot_cp_norm (p_msh, cp_n)

  fh = plot_map_msh (p_msh, cp_n, []);

  caxis ([0 1]);
  colorbar;

  set (gca (), "ydir", "reverse");

  xlabel ("x in mm");
  ylabel ("s_n in mm");

endfunction
