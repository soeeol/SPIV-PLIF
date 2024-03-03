##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## stitch sections based on their relative meshes with known absolute offsets
##
## Author: Sören J. Gerke
##

function [msh_gl dat_gl mask_g_gl mask_w_gl h_g_gl h_w_gl] = stitch_x (pdir, aid, msh, dat, masks, h_gw, ncoords, p_seam, rmbk)

  dat_gl = mask_g_gl = mask_w_gl = h_g_gl = h_w_gl = [];

  it_X = 1 : numel (aid.ids_X);

  [msh msh_gl msh_gl_sec lims_x sf] = stitch_x_msh (pdir, aid, msh, ncoords);

  ## x stitch Ic maps
  if (! isempty (dat))
    ## remove black response of camera chip
    if (rmbk)
      for i_X = it_X
        dat{i_X} = rm_blkr (pdir, dat{i_X});
      endfor
    endif

    dat_gl = stitch_x_dim (msh, dat, msh_gl_sec, it_X);

    ## smooth area around seams
    if (! isempty (p_seam))
      dat_gl = stitch_smooth_seams (p_seam, dat_gl, msh_gl, lims_x, sf);
    endif
  endif

  ## x stitch masks
  if (! isempty (masks))
    for i_X = it_X
        masks_wall_st{i_X} = {masks{i_X}.wall};
        masks_gas_st{i_X} = {masks{i_X}.gas};
    endfor
    mask_g_gl = stitch_x_dim (msh, masks_wall_st, msh_gl_sec, it_X);
    mask_w_gl = stitch_x_dim (msh, masks_gas_st, msh_gl_sec, it_X);
  endif

  ## x stitch contours
  if (! isempty (h_gw))
    for i_X = it_X
        h_wall_st{i_X} = h_gw{i_X}.wall;
        h_gas_st{i_X} = h_gw{i_X}.gas;
    endfor
    h_g_gl = stitch_x_dim (msh, h_wall_st, msh_gl_sec, it_X);
    h_w_gl = stitch_x_dim (msh, h_gas_st, msh_gl_sec, it_X);
  endif

endfunction
