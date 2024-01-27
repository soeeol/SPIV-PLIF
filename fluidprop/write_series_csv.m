##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## write header and data series seperately
##
## Author: Sören J. Gerke
##

function write_series_csv (fn, data, header, precision)
  head = cell(2);
  head(1,1) = mfilename;
  head(1,2) = date;
  if !isempty(header)
    for i = 1:numel(header)
      head(2,i) = header{i};
    endfor
  endif

  if isempty(precision)
    precision = "%.6e"
  endif
  cell2csv ([fn "_h.csv"], head);
  csvwrite ([fn "_d.csv"], data, "append", "off", "precision", precision)
endfunction
