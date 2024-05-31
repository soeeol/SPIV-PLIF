##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## method to read exported PIV data
##
## Author: Sören J. Gerke
##

function  [recu] = read_recu (filepath)

  files = glob (filepath);
  if ! isempty (files)
    for i = 1 : numel (files)
      [dir, name, ~] = fileparts (files{i});
      tmp = importdata (files{i}, ",");
      recu{i}.textdata = tmp.textdata;
      recu{i}.data = double (tmp.data);
    endfor
  else
    error (["no files found in " filepath]);
  endif

endfunction
