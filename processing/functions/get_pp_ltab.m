##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## read processing parameters from linking table
##
## Author: Sören J. Gerke
##

function [pp, idx_measid, head] = get_pp_ltab (ltab, pp)
  ## pp pos in table header
  head{1} = head{2} = cell ();
  i = 0;
  for [val, key] = pp
    i++;
    head{1}(i) = val.id;
    head{2}(i) = {find(ismember(ltab(1,:),val.id(:))==1)};
  endfor
  ## find measurement
  idx_measid = get_measidx (ltab, head, pp);
  ## read parameters from table
  for [val, key] = pp
    col = get_col (head, val.id(:));
    if (! isempty (col))
      pp.(key).data = ltab(idx_measid,:){col};
    endif
  endfor
endfunction
