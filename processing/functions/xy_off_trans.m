##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## translate mesh to minimize offset between Ic maps from wall contour
##
## Author: Sören J. Gerke
##

function [c_msh, xy_wall, pp] = xy_off_trans (c_msh, c_dat, pp, thrs)

  maxpass = 4;
  nswitch = 0;
  xoff = yoff = ones (1, numel (c_msh));
  xoff_c_fields = {"xoff_c", "xoff_c0", "xoff_c1"};
  yoff_c_fields = {"yoff_c", "yoff_c0", "yoff_c1"};

  ## init
  for i = 1 : numel (c_msh)

    if (isempty (pp.(xoff_c_fields{i}).data))
      pp.(xoff_c_fields{i}).data = 0.0;
    endif

    if isempty (pp.(yoff_c_fields{i}).data)
      pp.(yoff_c_fields{i}).data = 0.0;
    endif

  endfor

  ## minimize x and y offset
  k = 0;
  while ((max ([max(abs(xoff)), max(abs(yoff))]) > 1e-3) && (k <= maxpass))

    k++;
    if k > 1
       nswitch = 1;
    endif

    ## calculate y and x offset from detected wall contour and translate meshes
    xy_wall = update_wall_xy (c_msh, c_dat, pp, thrs);
    for i = 1 : numel (c_msh)
      [xoff(i), yoff(i)] = calc_offsets (xy_wall{i}, pp);

      ## update parameters
      pp.(xoff_c_fields{i}).data = pp.(xoff_c_fields{i}).data * nswitch + xoff(i);
      pp.(yoff_c_fields{i}).data = pp.(yoff_c_fields{i}).data * nswitch + yoff(i);

      ## translate meshes
      [c_msh{i}] = tform_mesh (c_msh{i}, pp, xoff_c_fields{i}, xoff(i));
      [c_msh{i}] = tform_mesh (c_msh{i}, pp, yoff_c_fields{i}, yoff(i));
    endfor
  endwhile

  ## final wall coordinate estimation
  xy_wall = update_wall_xy (c_msh, c_dat, pp, thrs);

endfunction
