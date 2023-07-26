##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Iterative call to x-stitch the Ic record derived data
##
## Author: Sören J. Gerke
##

function gl_Ic = stitch_all_2d_Ic (pdir, aid, pdat, it_A, it_C, it_M, it_X, ncoords, wlen, rmbk)
  for id_A = it_A
    for id_C = it_C
      for id_M = it_M
        printf([">>> c: stitching id_A = " num2str(id_A) ", id_C = " num2str(id_C) ", id_M = " num2str(id_M) " \n"]);
        [c_msh_gl, c_dat_gl, c_mask_g_gl, c_mask_w_gl, c_h_g_gl, c_h_w_gl] = ...
          stitch_x (pdir, pdat.pp(id_A, id_C, id_M, it_X), aid, pdat.c_msh(id_A, id_C, id_M, it_X), pdat.c_dat(id_A, id_C, id_M, it_X), pdat.c_masks(id_A, id_C, id_M, it_X), pdat.c_h(id_A, id_C, id_M, it_X), ncoords, wlen, rmbk);
        ## collect stitched data
        c_mask_g_gl_all{id_A, id_C, id_M} = c_mask_g_gl;
        c_mask_w_gl_all{id_A, id_C, id_M} = c_mask_w_gl;
        c_h_g_gl_all{id_A, id_C, id_M} = c_h_g_gl;
        c_h_w_gl_all{id_A, id_C, id_M} = c_h_w_gl;
        c_dat_gl_all{id_A, id_C, id_M} = c_dat_gl;
        c_msh_gl_all{id_A, id_C, id_M} = c_msh_gl;
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
