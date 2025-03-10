##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## filter stream lines generated by "even_stream_data.m"
##
## Author: Sören J. Gerke
##

function xy_filtered = streamlines_filter (xy, skip, minLength)

  printf (["number of lines pre filter: " num2str(length(xy)) "\n"])

  ## remove single constant points and NaNs
  x = xy(:,1);
  y = xy(:,2);
  valid_lines = !isnan(x) & !isnan(y) & (diff([x; NaN]) != 0 | diff([y; NaN]) != 0);
  keep_indices = valid_lines | ([true; diff(valid_lines) != 0]);
  x_valid = x(keep_indices);
  y_valid = y(keep_indices);

  ## find the NaN separators
  nan_indices = find (isnan (x_valid));

  ## Initialize new vectors
  new_x = [];
  new_y = [];

  ## Process each line
  start_idx = 1;
  for i = 1:length (nan_indices)
      ## current line
      end_idx = nan_indices(i) - 1;
      line_x = x_valid(start_idx:end_idx);
      line_y = y_valid(start_idx:end_idx);
      [new_x new_y] = streamlines_reduce (line_x, line_y, new_x, new_y, skip, minLength);
      ## start index for next line
      start_idx = nan_indices(i) + 1;
  end

  ## last line if it doesn't end with NaN
  if start_idx <= length(x)
    ## last line
    line_x = x_valid(start_idx:end);
    line_y = y_valid(start_idx:end);
    [new_x new_y] = streamlines_reduce (line_x, line_y, new_x, new_y, skip, minLength);
  end

  ## only keep unique subsequent points
  unique_indices = [true; or(diff(new_x) != 0, diff(new_y) != 0)];
  new_x = new_x(unique_indices);
  new_y = new_y(unique_indices);

  xy_filtered = [vec(new_x) vec(new_y)];

  printf (["number of lines post filter: " num2str(length(xy_filtered)) "\n"])

endfunction

