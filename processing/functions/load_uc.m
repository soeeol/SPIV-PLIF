##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## load processing result
##
## Author: Sören J. Gerke
##

function [msh_c cdata masks_c h_c msh_u udata masks_u h_u pp] = load_uc (fullpath, workdir)

  cd (fullpath);

  c_tmp = load ("c.v7");
  msh_c = c_tmp.c_msh;
  cdata = c_tmp.c_dat;
  masks_c = c_tmp.c_masks;
  h_c = c_tmp.c_h;

  u_tmp = load ("u.v7");
  msh_u = u_tmp.u_msh;
  udata = u_tmp.u_dat;
  masks_u = u_tmp.u_masks;
  h_u = u_tmp.u_h;

  pp_tmp = load ("pp.v7");
  pp = pp_tmp.pp;

  cd (workdir);

endfunction
