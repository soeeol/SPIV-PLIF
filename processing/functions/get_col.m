##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## get column index for matching header
##
## Author: Sören J. Gerke
##

function [idx] = get_col (header_filter, colid)
  idx = [];
  link_idx = find (ismember (header_filter{1}, colid) == 1);
  switch (numel (link_idx))
    case 0
      error (["colid " colid " was not found"]);
    case 1
      idx = header_filter{2}{link_idx};
    otherwise
      error (["colid " colid " was found " num2str(numel(link_idx)) " times"]);
  endswitch
endfunction
