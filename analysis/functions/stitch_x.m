##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## stitch sections based on their relative meshes with known absolute offsets
##
## Author: Sören J. Gerke
##

function [msh_gl dat_gl mask_g_gl mask_w_gl h_g_gl h_w_gl] = stitch_x (pdir, aid, msh, dat, masks, h_gw, ncoords, p_seam, rmbk)
  dat_gl = mask_g_gl = mask_w_gl = h_g_gl = h_w_gl = [];
  it_X = 1:numel(aid.ids_X);
  pp.optset.data = aid.ids_O{1}; # same optical setup for all sections per flow rate for now
  ## shift meshes
  for i_X = it_X
    msh{i_X}{1} = msh{i_X}{1} + aid.ids_X(i_X);
  endfor
  ## estimate scaling factors for new global mesh
  for i_X = it_X
    sfs(i_X,:) = get_sf (msh{i_X});
  endfor
  sf = min (sfs)
  ## stitching limits
  for i_X = it_X
    lims_y_min(i_X) = min (min (msh{i_X}{2}));
    lims_y_max(i_X) = max (max (msh{i_X}{2}));
  endfor
  ## limits for common grid
  domain = get_domain (pp);
  lims_y = [min(lims_y_min) max(lims_y_max)];
  lims_x_dom = [domain.xmin+0*sf(1), domain.xmax-1*sf(1)];
  x_end_plus = 0;
  for i_X = it_X
    if (i_X ==it_X(end))
      x_end_plus = sf(1);
    endif
    lims_x{i_X} = lims_x_dom + aid.ids_X(i_X) + x_end_plus;
  endfor
  lims_x{1}(1) = lims_x{1}(1,1) - domain.border;
  lims_x{end}(end) = lims_x{end}(end) + domain.border;
  ## TODO: detect overlap of data and interpolate
  ##
  ## stitching mesh
  for i_X = it_X
    msh_st{i_X} = reg_mesh (lims_x{i_X}, lims_y, sf);
  endfor
  msh_gl = stitch_x_cat (ncoords, it_X, msh_st);
  ## x stitch Ic maps
  if !isempty(dat)
    ## interpolate on part of global mesh
    for i_X = it_X
      nmap = numel(dat{i_X});
      ## remove black response of camera chip
      if (rmbk)
        dat{i_X} = rm_blkr (pdir, dat{i_X});
      endif
      for j = 1:nmap
        dat_st{i_X,j} = interp2 (msh{i_X}{1}, msh{i_X}{2}, dat{i_X}{j}, msh_st{i_X}{1}, msh_st{i_X}{2}, "pchip", NaN);
      endfor
    endfor
    ## combine maps
    dat_gl = stitch_x_cat (nmap, it_X, dat_st);
    ## smooth area around seams
    if !isempty(p_seam)
      ww = p_seam(1);
      ws = (ww - 1) / 2 + 1;
      for i = 1:nmap
        dat_gl_mm = movmedian (dat_gl{i}, 1*ww, 2, "Endpoints", 0.0);
        for j = 1:numel(aid.ids_X)-1
          xls_u = (lims_x{j+1}(1) + 1*ws*sf(1));
          xls_l = (lims_x{j}(2)   - 1*ws*sf(1));
          ## section replacing original data around stiching seam
          idx_xdom_s = (msh_gl{1} <= xls_u) & (msh_gl{1} >= xls_l);
          dat_gl{i} = dat_gl{i}.*(idx_xdom_s==0) + dat_gl_mm.*(idx_xdom_s==1);
        endfor
      endfor
    endif
  endif
  ## x stitch masks
  if !isempty(masks)
    for i_X = it_X
      ##
      mask_g_st{i_X} = interp2 (msh{i_X}{1}, msh{i_X}{2}, masks{i_X}.gas, msh_st{i_X}{1}, msh_st{i_X}{2}, "pchip", NaN);
      mask_w_st{i_X} = interp2 (msh{i_X}{1}, msh{i_X}{2}, masks{i_X}.wall, msh_st{i_X}{1}, msh_st{i_X}{2}, "pchip", NaN);
    endfor
    ## combine stitched masks
    mask_g_gl = stitch_x_cat (1, it_X, mask_g_st);
    mask_w_gl = stitch_x_cat (1, it_X, mask_w_st);
  endif
  ## x stitch contours
  if !isempty(h_gw)
    for i_X = it_X
      ##
      h_g_st{i_X} = interp1 (msh{i_X}{1}(1,:), h_gw{i_X}.gas, msh_st{i_X}{1}(1,:), "nearest", "extrap");
      h_w_st{i_X} = interp1 (msh{i_X}{1}(1,:), h_gw{i_X}.wall, msh_st{i_X}{1}(1,:), "nearest", "extrap");
    endfor
    ## combine wall and interface contour data
    h_g_gl = stitch_x_cat (1, it_X, h_g_st);
    h_w_gl = stitch_x_cat (1, it_X, h_w_st);
  endif
endfunction
