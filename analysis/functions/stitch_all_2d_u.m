##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Iterative call to assemble the u data
##
## Author: Sören J. Gerke
##

function gl_u = stitch_all_2d_u (pdir, aid, pdat, it_A, it_C, it_M, it_X, ncoords, wlen, rmbk)
  for i_A = it_A
    for i_C = it_C
      for i_M = it_M
        printf([">>> u: stitching i_A = " num2str(i_A) ", i_C = " num2str(i_C) ", i_M = " num2str(i_M) " \n"]);
        [u_msh_gl, u_dat_gl, u_mask_g_gl, u_mask_w_gl, u_h_g_gl, u_h_w_gl] = ...
          stitch_x (pdir, aid, pdat.u_msh(i_A, i_C, i_M, it_X), pdat.u_dat(i_A, i_C, i_M, it_X), pdat.u_masks(i_A, i_C, i_M, it_X), pdat.u_h(i_A, i_C, i_M, it_X), ncoords, wlen, rmbk);
        ## collect stitched data
        u_g_mask_gl_all{i_A,i_C,i_M} = u_mask_g_gl;
        u_w_mask_gl_all{i_A,i_C,i_M} = u_mask_w_gl;
        u_h_g_gl_all{i_A,i_C,i_M} = u_h_g_gl;
        u_h_w_gl_all{i_A,i_C,i_M} = u_h_w_gl;
        u_dat_gl_all{i_A,i_C,i_M} = u_dat_gl;
        u_msh_gl_all{i_A,i_C,i_M} = u_msh_gl;
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
