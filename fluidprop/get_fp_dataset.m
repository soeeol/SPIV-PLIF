##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## return requested data by fluid name, property name, source name from
## "fluidprop.v7" database
##
## Author: Sören J. Gerke
##

function dataset = get_fp_dataset (fprop, fname, pname, sname)
  ##i = 0
  ##for [key,val] = fprop
  ##  i++
  ##  key
  ##  val
  ##endfor
  for i = 1:numel(fprop)
    isin_f(i) = ismember (fname, fprop(i).fluid);
  endfor
  idx = find (isin_f);
  if !isempty(pname)
    for i = 1:numel(idx)
      for j = 1:numel(fprop(idx(i)).prop)
        isin_p(i,j) = ismember (pname, fprop(idx(i)).prop{j});
      endfor
    endfor
    [i, j] = find (isin_p);
    idx = idx(i);
  endif
  if !isempty(sname)
    for i = 1:numel(idx)
      for j = 1:numel(fprop(idx(i)).source)
        isin_s(i,j) = ismember (sname, fprop(idx(i)).source{j});
      endfor
    endfor
    [i, j] = find (isin_s);
    idx = idx(i);
  endif
  dataset = fprop(idx);

  function list = get_fp_avail (fprop)
    list = {fprop.fluid; fprop.prop};
  endfunction

endfunction

