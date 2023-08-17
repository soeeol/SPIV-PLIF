##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## measurement domain definition
##
## Author: Sören J. Gerke
##

function domain = get_domain (pp)
  switch (pp.optset.data)
    case {"M13", "M13b", "M13c"}
      domain.border = 0.5;
      domain.xmax = 4;
      domain.xmin = -4;
    case "M26"
      domain.border = 0.1;
      domain.xmax = 2;
      domain.xmin = -2;
    otherwise
      error ("no matching optical setup identifier");
  endswitch
endfunction
