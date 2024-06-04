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
  msh_var = sd.msh_var
  msh_file_x = {};
  for i_X = it_X
    msh_files = glob ([sd.a_dir id_meas{i_X} "/" "*" sd.a_id_sec "/" sd.msh_fn]);
    msh_file_x{i_X} = msh_files{end}; # use latest result
  endfor
  msh_file_x

  ## load section meshes
  M = collect_sec_var (msh_file_x, {msh_var});

  ## correct inter section mesh offsets
  ## xoff - applies to all types of section meshes so far
  if (isfield (ap.sd, "xoff_X"))
    xoff = ap.sd.xoff_X
    for i_X = 1 : numel (xoff)
      M.(msh_var){i_X}{1} = M.(msh_var){i_X}{1} + xoff(i_X);
    endfor
  endif
  ## yoff
  switch (msh_var)
    case {"c_msh", "u_msh"}
      if (isfield (ap.sd, "yoff_X"))
        yoff = ap.sd.yoff_X
        for i_X = 1 : numel (yoff)
          M.(msh_var){i_X}{2} = M.(msh_var){i_X}{2} + yoff(i_X);
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
    vars{i_f} = sd.dat_var{i_f};
    n_v = numel (vars{i_f}); # number of variables
    for i_v = 1:n_v
      var = vars{i_f}{i_v};
      S.(var) = cell (1, n_X);
    endfor
  endfor

  ##  load all
  for i_f = 1:n_f
    vars = sd.dat_var{i_f};
    n_v = numel (vars); # number of variables
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
    vars = sd.dat_var{i_f};
    for i_v = 1 : numel (vars)
      for i_X = 1:n_X
        dat_gl.(vars{i_v}){i_X} = [];
      endfor
    endfor
  endfor

  ## all variables
  data_vars = fieldnames (S)

  ## ensure that yoff is translated to y-dependent data variables
  if (isfield (ap.sd, "yoff_X"))
    yoff = ap.sd.yoff_X;
    for i_v = 1:n_v
      switch (data_vars{i_v})
        case {"delta_u_fit_avg"}
          for i_X = 1:n_X
            S.(data_vars{i_v}){i_X} = S.(data_vars{i_v}){i_X} + yoff(i_X);
          endfor
        case {"y_if_gas"}
          for i_X = 1:n_X
            for i_m = 1 : numel (S.(data_vars{i_v}){i_X})
              S.(data_vars{i_v}){i_X}{i_m} = S.(data_vars{i_v}){i_X}{i_m} + yoff(i_X);
            endfor
          endfor
      otherwise
      ;
      endswitch
    endfor
  endif

  ## x - stitch each variable
  for i_v = 1 : numel (data_vars)
    dat_gl.(data_vars{i_v}) = stitch_x_dim (msh_xs, S.(data_vars{i_v}), msh_gl_sec, it_X);
  endfor

endfunction
