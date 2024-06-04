##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## visual check of input recordings
##
## Author: Sören J. Gerke
##

function [is_valid, fh] = check_input_data (measid, cdata, udata)

  is_valid = [];
  cb_max = max (max (cdata{1}));

  titles = {"phi in a.u.", "phi des in a.u.", "phi sat in a.u.", "u mag in m/s"};

  n_c_maps = n_maps = numel (cdata);
  if (! isempty (udata))
    n_maps = n_c_maps + 1;
  endif

  fh = figure ();
  for i = 1:n_c_maps
    subplot (n_maps, 1, i);
    plot_map (cdata{i}, fh);
    axis tight;
    caxis ([0 max(max(cdata{1}))]);
    grid off;
    xlabel ("x in px");
    ylabel ("y in px");
    hax = colorbar ("location", "EastOutside");
    title (hax, titles{i})
  endfor
  if (! isempty (udata))
    subplot (n_maps, 1, n_maps);
    plot_map (udata{1}, fh);
    axis tight;
    colormap viridis;
    grid off;
    xlabel ("x in IA");
    ylabel ("y in IA");
    hax = colorbar ("location", "EastOutside");
    title (hax, titles{i})
  endif

  ##
  choice = questdlg ("does raw input fit measid?", "", "Yes", "No");
  if strcmp (choice, "Yes")
    is_valid = 1;
  else
    is_valid = 0;
    error ("check MeasID.csv to define correct combination of records");
  endif

endfunction
