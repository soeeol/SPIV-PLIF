##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## create "is liquid" mask covering field
##
## Author: Sören J. Gerke
##

function mask = masking (field, domain, size_map, y_min, y_if, sf, idx_off, out_val)

  mask = ones (size_map); # inside value

  switch domain

    case "gas"
      idx = 1 + int32 (ceil (abs (y_min - y_if) / sf(2)));
      for i = 1 : length (idx)
        mask(max([idx(i)+idx_off 1]):end,i) = out_val;
      endfor

    case "wall"
      idx = 1 + int32 (floor (abs (y_min - y_if) / sf(2)));
      for i = 1 : length (idx)
        mask(1:idx(i)-idx_off,i) = out_val;
      endfor

  endswitch

endfunction
