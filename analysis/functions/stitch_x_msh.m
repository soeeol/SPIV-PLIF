##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## prepare meshes for stitching
##
## msh        .. relative coordinates of measured data for several x scan positions
##                 x relative to calibtration origin
##                 y relative to wall
## msh_xs     .. msh offset in x dir by x scan positions required to interpolate on msh_gl_sec
## msh_gl     .. global coordinates including all x scan positions
##                 x relative to microstructure origin
## msh_gl_sec .. section parts of global coordinates for interpolation of x scan data
##
## Author: Sören J. Gerke
##

function [msh_xs msh_gl msh_gl_sec lims_x sf] = stitch_x_msh (msh, ap, ncoords)

  it_X = 1 : numel (ap.ids_X);

  pp.optset.data = ap.ids_O{1}; # same optical setup for all sections per flow rate for now

  msh_xs = msh;

  ## x shift msh to x scan center positions relative to microstructure
  for i_X = it_X
    msh_xs{i_X}{1} = msh{i_X}{1} + ap.ids_X(i_X);
  endfor

  ## estimate scaling factors for new global mesh
  for i_X = it_X
    sfs(i_X,:) = get_sf (msh{i_X});
  endfor

  sf = min (sfs);

  ## stitching limits
   for i_X = it_X
    ## x
    lims_x_min(i_X) = min (min (msh_xs{i_X}{1}));
    lims_x_max(i_X) = max (max (msh_xs{i_X}{1}));
    ## y
    lims_y_min(i_X) = min (min (msh_xs{i_X}{2}));
    lims_y_max(i_X) = max (max (msh_xs{i_X}{2}));
  endfor

  ## to strictly use the designed x section area
  ## .. not good in combination with manual xoff, since manual xpos stage had too much play
  domain = get_domain (pp);
##  lims_y = [min(lims_y_min) max(lims_y_max)];
##  lims_x_dom = [domain.xmin + 0*sf(1), domain.xmax - 1*sf(1)];
##  x_end_plus = 0;
##  for i_X = it_X
##    if (i_X == it_X(end))
##      x_end_plus = sf(1);
##    endif
##    lims_x{i_X} = lims_x_dom + ap.ids_X(i_X) + x_end_plus;
##  endfor
##  lims_x{1}(1) = lims_x{1}(1,1) - domain.border;
##  lims_x{end}(end) = lims_x{end}(end) + domain.border;

  ## .. rater use center of section overlap
  ## x
  lims_x{1} = [lims_x_min(1) mean([lims_x_max(1) lims_x_min(2)]) - sf(1)];
  for i_X = it_X(2:end-1)
    lims_x{i_X} = [mean([lims_x_max(i_X-1) lims_x_min(i_X)]) mean([lims_x_max(i_X) lims_x_min(i_X+1)]) - sf(1)];
  endfor
  lims_x{it_X(end)} = [mean([lims_x_max(it_X(end)-1) lims_x_min(it_X(end))]) lims_x_max(it_X(end))];
  ## ensure same extent of x for all ap.ids_M cases
  lims_x{1}(1) = domain.xmin + ap.ids_X(1) - domain.border;
  lims_x{it_X(end)}(2) = domain.xmax + ap.ids_X(it_X(end)) + domain.border;
  ## y
  lims_y = [min(lims_y_min) max(lims_y_max)];

  ## stitching mesh
  for i_X = it_X
    msh_gl_sec{i_X} = reg_mesh (lims_x{i_X}, lims_y, sf);
  endfor

  ## combined mesh
  msh_gl = stitch_x_cat (ncoords, it_X, msh_gl_sec);

endfunction
