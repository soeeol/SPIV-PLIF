##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## extract from scalar 2d field on interface normal lines
##
## Author: Sören J. Gerke
##

function [cp, msh_p, msh_n] = interp_if_norm (snp, h_g, msh, cn)

  x = msh{1}(1,:);
  z = msh{3}(1,1);

  if (isfloat (cn) && (ndims (cn) == 2))
    cn = {cn};
    h_g = {h_g};
  endif

  if (! iscell (cn))
    error ("interp_if_norm: expected type for cn: cell")
  endif

  n_t = numel (cn);

  if (! (numel (h_g) == n_t))
    error ("interp_if_norm: number of fields (cn) does not match number of interfaces (h_g)")
  endif

  cp = cell (1, n_t);

  ## interpolate normalized concentration field on line mesh
  for i_t = 1:n_t # time steps
    printf ([">>> interp_if_norm: " num2str(i_t) " of " num2str(n_t) " fields\n"]);

    ## interface normal line mesh
    msh_n = if_norm_mesh (x, h_g{i_t}, snp);

    ## interpolate scalar on interface normal line mesh
    cp{i_t} = interp2 (msh{1}, msh{2}, cn{i_t}, msh_n{1}, msh_n{2}, "pchip", 0.0);
  endfor

  ## regular mesh for the interface normal profiles
  [XX, YY, ZZ] = meshgrid (x, snp, z);
  msh_p = {XX, YY, ZZ};

endfunction

