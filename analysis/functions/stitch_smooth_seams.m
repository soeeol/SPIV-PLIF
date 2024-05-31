##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## smooth stitching seams
## .. not useful if replacement section borders introduce jumps
##
## Author: Sören J. Gerke
##

function dat_gl = stitch_smooth_seams (p_seam, dat_gl, msh_gl, lims_x, sf)

  ww = p_seam(1);
  ws = (ww - 1) / 2 + 1;
  nmaps = numel (dat_gl);

  for i = 1 : nmaps
    dat_gl_mm = movmedian (dat_gl{i}, 1*ww, 2, "endpoints", 0.0);
    for j = 1 : numel (lims_x) - 1
      ## section replacing original data around stiching seam
      xls_u = (lims_x{j+1}(1) + 1*ws*sf(1));
      xls_l = (lims_x{j}(2)   - 1*ws*sf(1));
      idx_xdom_s = (msh_gl{1} <= xls_u) & (msh_gl{1} >= xls_l);
      ## replace section with moving average
      dat_gl{i} = dat_gl{i}.*(idx_xdom_s==0) + dat_gl_mm.*(idx_xdom_s==1);
    endfor
  endfor

endfunction
