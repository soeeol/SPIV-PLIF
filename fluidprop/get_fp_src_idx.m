##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## return index of datasets in db (fluidprop database) matching src
##
## Author: Sören J. Gerke
##

function idx = get_fp_src_idx (db, src)

  isin = false (1, numel (db));

  for i = 1:numel(db)
    if iscell (db(i).source)
      isin(i) = ismember (src, db(i).source{1});
    endif
  endfor

  if (sum (isin) > 0)
    idx = find (isin);
  else
    idx = 0;
  endif

endfunction
