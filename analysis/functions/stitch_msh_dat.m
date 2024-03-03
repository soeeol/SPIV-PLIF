##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## assembly
##
## var_list ..  variables to be loaded and assembled. first one has to be the name of mesh variable
##
##
## Author: Sören J. Gerke
##
##

function [msh_gl dat_gl] = stitch_msh_dat (fn_data, fn_sec, var_list, avg_method, aid)

    ## number of sections
    n_s = numel (fn_sec)

    ## load section meshes
    M = collect_sec_var (fn_sec, fn_data, var_list(1));

    ## build global mesh
    [msh_xs msh_gl msh_gl_sec ~] = stitch_x_msh (M.(var_list{1}), aid, 3);

    ## load section data
    dat = collect_sec_var (fn_sec, fn_data, var_list(2:end));
    var_names = fieldnames (dat);

    ## number of data variables
    n_v = numel (var_names)

    ## averaging method
    avgm = cell (n_v, 1);
    if (isempty (avg_method))
      avgm(:) = "median";
    elseif (abs (numel (avg_method) - n_v) > 0)
      avgm(:) = avg_method{1};
    else
      avgm = avg_method;
    endif

    ## stitching for sequential measurements only makes sense for avg values
    ## thus before stitching take average if time series is found in dat

    ## init avg data struct
    dat_avg = [];
    for i_v = 1:n_v
      dat_avg.(var_names{i_v}) = cell (n_s, 1);
    endfor

    ## avg mean or median per section if more than one time step is avail per variable and per section
    for i_v = 1:n_v
      for i_s = 1:n_s
        sec_var = dat.(var_names{i_v}){i_s};
        sec_avg_var = [];
        if (iscell (sec_var) && (numel (sec_var) > 1))
          if (isvector (sec_var{2}) || iscolumn (sec_var{2}))
            sec_avg_var = calc_vec_avg_cells (sec_var, avgm{i_v});
          else
            sec_avg_var = calc_im_avg_cells (sec_var, avgm{i_v});
          endif
          dat_avg.(var_names{i_v}){i_s} = sec_avg_var;
        else
          dat_avg.(var_names{i_v}){i_s} = sec_var;
        endif
      endfor
    endfor

    ## init global (assembled) struct
    dat_gl = [];
    for i_v = 1:n_v
      dat_gl.(var_names{i_v}) = [];
    endfor

    ## x - stitch each variable
    it_X = 1 : numel (aid.ids_X);
    for i_v = 1:n_v
      dat_gl.(var_names{i_v}) = stitch_x_dim (msh_xs, dat_avg.(var_names{i_v}), msh_gl_sec, it_X);
    endfor

endfunction
