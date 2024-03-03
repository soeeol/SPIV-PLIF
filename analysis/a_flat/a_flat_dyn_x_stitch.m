##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## template x-section assembly
##
## for a_flat_dyn (section wise analysis)
##
## Author: Sören J. Gerke
##

## (0) init
if 1
  aid.ids_L = {"WG141"};
  aid.ids_O = {"M13"};
  aid.ids_C = {"flat"};
  aid.ids_A = [60];
  aid.ids_M = [8 16 32 64];
  aid.ids_X = [-8 0 8 16];
  id_G = 2;
  id_T = 25
  id_Z = 0;

  a_type = "a_flat_dyn_stitch";
  c_method = "linear"; # "linear" "nonlin"
  c_if_method = "calib"; # "calib" "calib-if"

  ## iterators
  it_A = 1 : numel (aid.ids_A); # angles
  it_C = 1 : numel (aid.ids_C); # cells
  it_M = 1 : numel (aid.ids_M); # mass flow rates
  it_X = 1 : numel (aid.ids_X); # scanned x sections

  ## fixed
  i_L = i_O = 1; # liquid, optical setup
  i_A = 1;
  i_C = 1;
  ## overrides
##  i_M = it_M = 3; # single mass flow rate
##  i_X = it_X = 4; # single scanned x sections

  ## main result directories
  result_dir = [pdir.analyzed a_type "/"]

endif

## (_) x - stitch flat dyn avg result and save combined
if 1

  ## input data path of per section analysis
  sec_data_dir = [pdir.analyzed "a_flat_dyn/"]

  ## per section analysis method was
  a_id = ["cp-avg-" c_if_method "-" "cn-" c_method]

  ## file holding the data to be assembled
  fn_data{1} = "dyn_cn_cp.v7"
  fn_data{2} = "dyn_cn_cp_avg_fit.v7"

  ## variables to be loaded and assembled. first one has to be the name of mesh variable
  var_list{1} = {"p_msh" "cp_n" "cp_s" "cp_b"}
  var_list{2} = {"p_msh_o" "cp_nn" "delta_c" "h_g_mean"}
  n_f = numel (var_list);

  ## per data variable: method for time series averaging before stitching
  ## ("median" or "mean")
  ## or set one for all
  avg_method = {"median"}

  ## resulting stitching analysis identifier
  stitch_a_id = ["x-stitch_" a_id]

  for i_M = it_M
    for i_f = 1:n_f
      ## prepare result storage path
      id_meas = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X, id_Z)
      save_dir_id = [result_dir id_meas "/" date_str "_" stitch_a_id "/"];
      mkdir (save_dir_id);

      ## list analysis results to combine
      for i_X = it_X
        dyncase = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X(i_X), id_Z);
        fn_sec_data = glob ([sec_data_dir dyncase "/*_" a_id]);
        fn_sec_a{i_X} = fn_sec_data{end}; # use latest result
      endfor

      ## combine mesh and data of sections
      [msh_gl dat_gl] = stitch_msh_dat (fn_data{i_f}, fn_sec_a, var_list{i_f}, avg_method, aid);

      ## save data
      cd (save_dir_id)
      save ("-v7", ["x-stitch_msh_dat___" fn_data{i_f}], "msh_gl", "dat_gl")
    endfor
  endfor

endif

if 0

  figure ()
  plot (msh_gl{1}(1,:), dat_gl.delta_c, "k")

  figure ()
  plot (msh_gl{1}(1,:), dat_gl.h_g_mean, "k")

  plot_cp_norm (msh_gl, dat_gl.cp_nn)
  hold on
  plot (msh_gl{1}(1,:), dat_gl.delta_c, "r")

  figure ()
  hold on
  plot (dat_gl.cp_s, [";s;"])
  plot (dat_gl.cp_b, [";b;"])

endif
