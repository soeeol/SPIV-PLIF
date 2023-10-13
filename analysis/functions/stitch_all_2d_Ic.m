##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Iterative call to x-stitch the Ic record derived data
##
## Author: Sören J. Gerke
##

function gl_Ic = stitch_all_2d_Ic (pdir, aid, pdat, it_A, it_C, it_M, it_X, ncoords, wlen, rmbk)
  for i_A = it_A
    for i_C = it_C
      for i_M = it_M
        printf([">>> c: stitching i_A = " num2str(i_A) ", i_C = " num2str(i_C) ", i_M = " num2str(i_M) " \n"]);
        [c_msh_gl, c_dat_gl, c_mask_g_gl, c_mask_w_gl, c_h_g_gl, c_h_w_gl] = ...
          stitch_x (pdir, aid, pdat.c_msh(i_A, i_C, i_M, it_X), pdat.c_dat(i_A, i_C, i_M, it_X), pdat.c_masks(i_A, i_C, i_M, it_X), pdat.c_h(i_A, i_C, i_M, it_X), ncoords, wlen, rmbk);
        ## collect stitched data
        c_mask_g_gl_all{i_A,i_C,i_M} = c_mask_g_gl;
        c_mask_w_gl_all{i_A,i_C,i_M} = c_mask_w_gl;
        c_h_g_gl_all{i_A,i_C,i_M} = c_h_g_gl;
        c_h_w_gl_all{i_A,i_C,i_M} = c_h_w_gl;
        c_dat_gl_all{i_A,i_C,i_M} = c_dat_gl;
        c_msh_gl_all{i_A,i_C,i_M} = c_msh_gl;
      endfor
    endfor
  endfor
  gl_Ic.c_g_mask = c_mask_g_gl_all;
  gl_Ic.c_w_mask = c_mask_w_gl_all;
  gl_Ic.c_h_g = c_h_g_gl_all;
  gl_Ic.c_h_w = c_h_w_gl_all;
  gl_Ic.c_dat = c_dat_gl_all;
  gl_Ic.c_msh = c_msh_gl_all;
endfunction
