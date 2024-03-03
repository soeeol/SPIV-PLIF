##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## return requested datasets from "fluidprop.v7" database. First lookup is for fluid name, then for property name, then for source name.
## property or source name can be empty.
##
## Author: Sören J. Gerke
##

function dataset = get_fp_dataset (fprop, fname, pname, sname)

  for i = 1 : numel(fprop)
    isin_f(i) = ismember (fname, fprop(i).fluid);
  endfor

  idx = find (isin_f);

  if (! isempty (pname))
    for i = 1 : numel (idx)
      for j = 1 : numel (fprop(idx(i)).prop)
        isin_p(i,j) = ismember (pname, fprop(idx(i)).prop{j});
      endfor
    endfor
    [i, j] = find (isin_p);
    idx = idx(i);
  endif

  if (! isempty (sname))
    idx_i = false (1, numel(idx));
    for i = 1 : numel (idx)
      idx_i(i) = get_fp_src_idx (fprop(idx(i)), sname{1});
    endfor
    idx = idx(idx_i);
  endif

  dataset = fprop(idx);

endfunction
