##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## read size of exported vector fields
##
## Author: Sören J. Gerke
##

function  [sx, sy] = get_vecmapsize (textdata)
  infostr = textdata{4};
  idx = strchr (infostr, "=,}");
  switch length (idx)
    case 4
      sx = str2num (infostr(idx(1)+1:idx(2)-1));
      sy = str2num (infostr(idx(3)+1:idx(4)-1));
    otherwise
      error ("get_vecmapsize: textdata header format changed?");
  endswitch
endfunction
