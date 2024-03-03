##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## interface normal line mesh
##
## px, py .. coordinates of interface
## snp .. interface normal coordinates
##
## requires the matgeom package
##
## Author: Sören J. Gerke
##

function msh_n = if_norm_mesh (px, py, snp)

  px = vec (px);
  py = vec (py);
  snp = vec (snp);

  ## local angle of inclination
  ap = angle2Points ([px(1:end-1) py(1:end-1)], [px(2:end) py(2:end)]);

  ## number of lines
  n_l = numel (px);
  ## number of points along edge
  n_p = numel (snp);

  ## lines normal to interface
  norm_lines = createLine ([px py], pi / 2 + [ap; ap(end)] .* ones (n_l, 1));

  ## create the line mesh
  msh_n_x = msh_n_y = zeros (n_p, n_l);
  for i_l = 1 : n_l
    edge = createEdge (norm_lines(i_l,:) .* ones (n_p, 1), -1 * snp);
    msh_n_x(:,i_l) = edge(:,3);
    msh_n_y(:,i_l) = edge(:,4);
  endfor
  msh_n = {msh_n_x, msh_n_y, zeros (size (msh_n_x))};

endfunction
