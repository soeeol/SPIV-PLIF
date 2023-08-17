##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Iterative call to assemble the u data
##
## Author: Sören J. Gerke
##

function gl_u = stitch_all_2d_u (pdir, aid, pdat, it_A, it_C, it_M, it_X, ncoords, wlen, rmbk)
  for id_A = it_A
    for id_C = it_C
      for id_M = it_M
        printf([">>> u: stitching id_A = " num2str(id_A) ", id_C = " num2str(id_C) ", id_M = " num2str(id_M) " \n"]);
        [u_msh_gl, u_dat_gl, u_mask_g_gl, u_mask_w_gl, u_h_g_gl, u_h_w_gl] = ...
          stitch_x (pdir, pdat.pp(id_A, id_C, id_M, it_X), aid, pdat.u_msh(id_A, id_C, id_M, it_X), pdat.u_dat(id_A, id_C, id_M, it_X), pdat.u_masks(id_A, id_C, id_M, it_X), pdat.u_h(id_A, id_C, id_M, it_X), ncoords, wlen, rmbk);
        ## collect stitched data
        u_g_mask_gl_all{id_A, id_C, id_M} = u_mask_g_gl;
        u_w_mask_gl_all{id_A, id_C, id_M} = u_mask_w_gl;
        u_h_g_gl_all{id_A, id_C, id_M} = u_h_g_gl;
        u_h_w_gl_all{id_A, id_C, id_M} = u_h_w_gl;
        u_dat_gl_all{id_A, id_C, id_M} = u_dat_gl;
        u_msh_gl_all{id_A, id_C, id_M} = u_msh_gl;
      endfor
    endfor
  endfor
  gl_u.u_g_mask = u_g_mask_gl_all;
  gl_u.u_w_mask = u_w_mask_gl_all;
  gl_u.u_h_g = u_h_g_gl_all;
  gl_u.u_h_w = u_h_w_gl_all;
  gl_u.u_dat = u_dat_gl_all;
  gl_u.u_msh = u_msh_gl_all;
endfunction
