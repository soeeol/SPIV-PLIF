##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## draw structure geometry to currect axis
##
## Author: Sören J. Gerke
##

function draw_ms (ms_type, xoff, opts)
  switch ms_type
    case "r10"
      for i = 1:numel(xoff)
        rectangle ("position", [-1+xoff(i)-0.0 0.0 wms=2.0 hms=1.0], opts{:})
        ## 45-6 m higher structure from mould CNC machining offset?
        ## or on all sides of structure due to PDMS swelling?
      endfor
    case "r15"
      for i = 1:numel(xoff)
        rectangle ("position", [-1+xoff(i) 0 wms=2 hms=1.5], opts{:})
      endfor
    case "r20"
      for i = 1:numel(xoff)
        rectangle ("position", [-1+xoff(i) 0 wms=2 hms=2], opts{:})
      endfor
    case "t10"
      for i = 1:numel(xoff)
        hms = 1; wms = 2;
        patch ([-0.5*wms 0 .5*wms]' + xoff(i), [0 hms 0]', opts{:})
      endfor
    case "s10"
       hms = 2; wms = 2;
      for i = 1:numel(xoff)
        rectangle ("position", [-1+xoff(i) -1 wms hms], "curvature", [1 1], opts{:});
      endfor
  endswitch
endfunction
