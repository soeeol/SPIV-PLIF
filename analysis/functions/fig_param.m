##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## figure parameter definition
##
## Author: Sören J. Gerke
##

function fp = fig_param (fig_def)
  fig_type = fig_def{1};
  fig_width = fig_def{2};
  publisher = fig_def{3};
  switch publisher
    case "elsevier"
      ## based on https://www.elsevier.com/authors/policies-and-guidelines/artwork-and-media-instructions/artwork-sizing
      ## simple adaption of above rules for rasterized plotting
      fp.fonts = {"Arial","Courier","Symbol","Times"};
      ## Line weights range from 0.10 pt to 1.5 pt
      fp.lw_min = 0.1;
      fp.lw_max = 1.5;
      ## font sizes sub 6 min 7 pt
      fp.fs_sub = 6;
      fp.fs_min = 7;
      ## printing resolution
      fp.dpi = [];
      switch fig_type
        case "combo"
          fp.dpi = 500;
        case "lineart"
          fp.dpi = 1000;
        case "image"
          fp.dpi = 300;
      endswitch
      ## standard sizes
      fp.xsize = [];
      switch fig_width
        case "min"
          fp.xsize = 3; # cm
        case "single"
          fp.xsize = 9;
        case "onehalf"
          fp.xsize = 14;
        case "double"
          fp.xsize = 19;
      endswitch
  endswitch
endfunction
