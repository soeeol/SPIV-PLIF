##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## create "is liquid" mask covering field
##
## Author: Sören J. Gerke
##

function mask = masking (field, domain, sm, ymin, hi, sf, offset, outside)

  mask = ones (sm); # inside

  switch domain
    case "gas"
      switch field
        case "u"
          idx = 1 + int32 (ceil (abs (ymin - hi) / sf(2)));
        case "c"
          idx = 1 + int32 (ceil (abs (ymin - hi) / sf(2)));
      endswitch
      for i = 1 : length (idx)
        mask(max([idx(i)+offset 1]):end,i) = outside;
      endfor

    case "wall"
      switch field
        case "u"
          idx = 1 + int32 (floor (abs (ymin - hi) / sf(2)));
        case "c"
          idx = 1 + int32 (floor (abs (ymin - hi) / sf(2)));
      endswitch
      for i = 1 : length (idx)
        mask(1:idx(i)-offset,i) = outside;
      endfor

    otherwise
      error (["masking: no masking for domain" domain]);
  endswitch

endfunction
