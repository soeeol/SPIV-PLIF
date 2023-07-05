##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## store per section processing result
##
## Author: Sören J. Gerke
##

function status = save_result (fullpath, workdir, c_msh, c_dat, c_masks, c_h, u_msh, u_dat, u_masks, u_h, pp)
  cd (fullpath);

  ##
  files = glob ([fullpath "/*.v7"]);
  if (isempty(files))
    try
      save -v7 "c.v7" c_msh c_dat c_masks c_h
      save -v7 "u.v7" u_msh u_dat u_masks u_h
      save -v7 "pp.v7" pp
      status = true;
    catch err
       status = false;
    end_try_catch
  else
    status = false;
    error (["found v7 files in " fullpath]);
  endif
  cd (workdir);
endfunction
