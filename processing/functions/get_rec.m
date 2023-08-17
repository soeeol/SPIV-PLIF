##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## get record by identifier string
##
## Author: Sören J. Gerke
##

function [rec] = get_rec (recs, recids, id)
  rec = recs{find(ismember(recids,id)==1)}{1};
endfunction
