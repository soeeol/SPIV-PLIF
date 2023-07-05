##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## method to read exported PLIF image
##
## Author: Sören J. Gerke
##

function [recc] = read_recc (filepath)
  ## tif images from recording are uint16 ranging from 0 to 2^16-1 = 65535
  ## further calculations in SINGLE precision showed to be sufficient
  ## accurate reprensentation of 8 significant digits
  files = glob (filepath);
  if ! isempty (files)
    for i = 1:numel (files)
      recc{i} = im2single (flip (imread (files{i}), 1));
    endfor
  else
    error (["no files found in " filepath]);
  endif
endfunction
