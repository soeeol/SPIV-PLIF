##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## visual check of input recordings
##
## Author: Sören J. Gerke
##

function valid = check_input_data (measid, cdata, udata)
  valid = [];
  cb_max = max (max (cdata{1}));
  ## check raw data before processing
  fha = figure ();
  ## subimage to solve same colormap problem?
  subplot (2, 2, 1)
  surf (cdata{1})
  shading flat
  view ([0 0 1])
  axis image
  title ("Ic", "interpreter", "none");
  xlabel ("pixel"); ylabel ("pixel");
  caxis ([0 cb_max])
  ##
  subplot (2, 2, 3)
  surf (cdata{2})
  shading flat
  view ([0 0 1])
  axis image
  title ("Ic0", "interpreter", "none");
  xlabel ("pixel"); ylabel ("pixel");
  caxis ([0 cb_max])
  ##
  subplot (2, 2, 2)
  if (numel(cdata)>2)
    surf (cdata{3})
    shading flat
    view ([0 0 1])
    axis image
    caxis ([0 cb_max])
    xlabel ("pixel"); ylabel ("pixel");
  else
    axis off;
  endif
  title ("Ic1", "interpreter", "none");
  ##
  if !(isempty(udata))
    subplot (2, 2, 4)
    surf (udata{4})
    shading flat
    view ([0 0 1])
    axis image
    title ("u_mag","interpreter", "none");
    xlabel ("IA"); ylabel ("IA");
  endif
  axes ("visible", "off");
  title (measid, "interpreter", "none");
  ##
  set (fha, 'Position', get(0, 'Screensize'));
  hold off
  pause (2)
  choice = questdlg ("raw input fit measid?", "", "Yes", "No");
  if strcmp (choice, "Yes")
    valid = 1;
    close (fha)
  else
    valid = 0;
    error ("check MeadID.csv to define correct combination of records");
  endif
endfunction
