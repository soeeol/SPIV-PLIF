##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## assembly of per section data assuming variables defined in sd are temporal avg
##
## sd .. stitching descriptor
## ap .. analysis parameters
##
## Author: Sören J. Gerke
##

function [msh_gl dat_gl] = stitch_msh_dat (sd, ap)

  n_f = numel (sd.dat_fn);
  n_X = numel (ap.ids_X);
  it_X = 1:n_X;

  ## per section measurement ids
  for i_X = it_X
    ap.i_X = i_X;
    id_meas{i_X} = get_measid_ap (ap);
  endfor

  ##
  ## build combined mesh
  ##

  ## find mesh files
  msh_var = sd.msh_var;
  msh_file_x = {}
  for i_X = it_X
    msh_files = glob ([sd.a_dir id_meas{i_X} "/" "*_" sd.a_id_sec "/" sd.msh_fn]);
    msh_file_x{i_X} = msh_files{end}; # use latest result
  endfor

  ## load section meshes
  M = collect_sec_var (msh_file_x, {msh_var});

  ## correct inter section mesh offsets
  switch (msh_var)
    case {"c_msh"}
      if (isfield (ap, "yoff_X"))
        yoff = ap.yoff;
        for i_X = it_X
          M.(msh_var){i_X}{2} = M.(msh_var){i_X}{2} + ap.yoff_X(i_X);
        endfor
      endif
    otherwise
      ;
  endswitch

  ## build meshes
  [msh_xs msh_gl msh_gl_sec ~] = stitch_x_msh (M.(msh_var), ap, 3);

  ##
  ## combine data
  ##

  ## find relevant data files
  for i_f = 1:n_f
    for i_X = it_X
      dat_files = glob ([sd.a_dir id_meas{i_X} "/" "*_" sd.a_id_sec "/" sd.dat_fn{i_f}]);
      dat_file_X{i_f,i_X} = dat_files{end}; # use latest result
    endfor
  endfor

  ## prepare struct
  S = [];
  for i_f = 1:n_f
    vars{i_f} = sd.dat_var{i_f}
    n_v = numel (vars{i_f}) # number of variables
    for i_v = 1:n_v
      var = vars{i_f}{i_v}
      S.(var) = cell (1, n_X);
    endfor
  endfor

  ##  load all
  for i_f = 1:n_f
    vars = sd.dat_var{i_f}
    n_v = numel (vars) # number of variables
    for i_v = 1:n_v
      for i_X = 1:n_X
        AA = load (dat_file_X{i_f,i_X}, vars{i_v});
        S.(vars{i_v}){i_X} = AA.(vars{i_v});
      endfor
    endfor
  endfor

  ## init global struct for assembled data
  dat_gl = [];
  for i_f = 1 : numel (sd.dat_fn)
    vars = sd.dat_var{i_f}
    for i_v = 1 : numel (vars)
      for i_X = 1:n_X
        dat_gl.(vars{i_v}){i_X} = [];
      endfor
    endfor
  endfor

  ## ensure that yoffset is translated to y-dependent data variables
  if (isfield (ap, "yoff"))
    yoff = ap.yoff;
    for i_v = 1:n_v
      switch (data_vars{i_v})
        case {"delta_u_fit_avg"}
          for i_X = 1:n_X
            sec_avg_var = dat_avg.(data_vars{i_v}){i_X};
            dat_avg.(data_vars{i_v}){i_X} = sec_avg_var + yoff(i_X);
          endfor
      otherwise
      ;
      endswitch
    endfor
  endif

  ## x - sd each variable
  data_vars = fieldnames (S)
  for i_v = 1 : numel (data_vars)
    dat_gl.(data_vars{i_v}) = stitch_x_dim (msh_xs, S.(data_vars{i_v}), msh_gl_sec, it_X);
  endfor

endfunction
