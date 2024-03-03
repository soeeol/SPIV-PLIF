##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## draw the ideal wall geometry based on the cell identifier
##
## Author: Sören J. Gerke
##

function draw_cell (id_cell, xpos, is_rel_pos)

  opts = {"facecolor", 0.7 * ones(1,3), "facealpha", 0.75, "edgecolor", "none"};

  ww = 32; # in mm; full range of x scan with M13
  wms = 2; # in mm; structure base length

  if (isempty (is_rel_pos))
    is_rel_pos = true;
  endif

  if (isempty (xpos))
    if is_rel_pos
      xpos = 0;
    else
      xpos = 60;
    endif
  endif

  ## draw id_cell "flat" as base structure
  rectangle ("position", [-12+xpos -0.25 ww 0.25], opts{:})

  ms_type = [];
  xoff = 0;
  switch id_cell
    case {"2d-r10", "2d-r10-20", "2d-r10-40", "2d-r10-60"}
      ms_type = "r10";
    case {"2d-t10"}
      ms_type = "t10";
    case {"2d-c10"}
      ms_type = "c10";
  endswitch

  switch id_cell
    case {"2d-r10", "2d-r15", "2d-r20", "2d-c10", "2d-t10"}
      xoff = 0 + xpos;
    case {"2d-r10-20"}
      xoff = ((1:3) - 1) * (2 + wms) + xpos;
    case {"2d-r10-40"}
      xoff = ((1:3) - 1) * (4 + wms) + xpos;
    case {"2d-r10-60"}
      xoff = ((1:3) - 1) * (6 + wms) + xpos;
  endswitch

  ## draw structure
  draw_ms (ms_type, xoff, opts);

endfunction
