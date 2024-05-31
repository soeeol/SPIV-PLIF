##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## single frame dynamic PLIF analysis
## interface statistics and boundary layer thickness
##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ## ap parameters defining the analysis
  ap = [];
  ap.a_type = "a_2DR10_dyn_cn"; # identifier of this analysis
  ap.p_type = "2d_dyn_Ic"; # for the analysis use the output of processing procedure "p_type"
  ap.c_method = "linear"; # method to transform fluorescence intensity to concentration ("linear" / "nonlinear" .. no impact on delta_c)
  ap.c_if_method = "calib"; # method to deal with fluorescence intensity decay at the interface ("calib" / "calib-if" .. high impact on delta_c)

  ## selection of experiments to be analyzed
  ap.ids_A = [60]; # [°] inlination IDs
  ap.ids_C = {"2d-r10"}; # cell IDs
  ap.ids_G = [2]; # [Nl/min] gas flow IDs
  ap.ids_O = {"M13"}; # optical setup IDs
  ap.ids_M = [8 16 32 64]; # [kg/h] mass flow IDs
  ap.ids_L = {"WG141"}; # liquid IDs
  ap.ids_T = [25]; # [°C] temperature IDs
  ap.ids_X = [-8 0 8 16]; # [mm] x section IDs
  ap.ids_Z = [0]; # [mm] z position of light sheet relative to center of cell

  ## iterators
  it_A = 1 : numel (ap.ids_A); # angles
  it_C = 1 : numel (ap.ids_C); # cells
  it_M = 1 : numel (ap.ids_M); # mass flow rates
  it_X = 1 : numel (ap.ids_X); # scanned x sections
  ## fixed
  i_A = 1; ap.i_A = i_A;
  i_C = 1; ap.i_C = i_C;
  i_G = 1; ap.i_G = i_G;
  i_L = 1; ap.i_L = i_L;
  i_O = 1; ap.i_O = i_O;
  i_T = 1; ap.i_T = i_T;
  i_Z = 1; ap.i_Z = i_Z;

  ## overrides
  i_M = it_M = 2
##  i_X = it_X = 2

  ap.dyn_cn_if_scmad_dev_max = 0.25; # used for threshold estimation for valid median interface deviation
  ap.dyn_cn_nt_max = 20; # limit number of valid single frame used for the analysis (valid: small deviation to median interface)
  ap.c_calib_sig_X = 0; # c calibration reference smoothing factor per ids_X (default)
  ap.c_isec_off_shift_lim = [0.25 0.1]; # [mm] inter section phi offset limits (default)
  ap.c_isec_rcurv_lim = 200; # [mm] threshold for "flat" interface curvature radius
  ## parameters for interface normal concentration profile estimation
  ap.cp_pd_M = [0.35 0.45 0.5 0.5]; # interface normal profile depth for each ids_M [mm]
  ## spline interface fit for stable interface normals (defaults)
  ap.cp_if_sfit_order = 3; # spline order
  ap.cp_if_sfit_sps = 9; # splines devisions per section (.. generally: increase for curved interface, decrease for flat)

  ## prepare directories
  ap.result_dir = [pdir.analyzed ap.a_type "/"];

endif

## [10] normalized concentration field per x-section
if 1

  ## analysis identifier
  sec_a_id = ["cn-" ap.c_method "_" ap.c_if_method]

  for i_X = it_X
    for i_M = it_M

      close all

      ap.i_M = i_M;
      ap.i_X = i_X;

      ## measurement identifier
      id_meas = get_measid_ap (ap)

      ## path to store per section analysis results and plots to
      ap.save_dir_id = [ap.result_dir id_meas "/" date_str "_" sec_a_id "/"];
      mkdir (ap.save_dir_id);

      c_calib_sig_X = []
      switch (i_M)
        case 1
          c_calib_sig_X = [0 0 0 0] # for i_M=1&i_X=1 saturation recorded film was slightly thinner and thus of lower fluorescence
        case 2
          c_calib_sig_X = [0 0 0 0] # TODO: first check for inter section offset
        case 3
          c_calib_sig_X = [0 0 0 0] # TODO: first check for inter section offset
        case 4
          c_calib_sig_X = [0 0 0 0] # TODO: first check for inter section offset
      endswitch
      if (! isempty (c_calib_sig_X))
        ap.c_calib_sig_X = c_calib_sig_X(i_M);
      endif

      [~] = a_dyn_cn___avg_calib (pdir, ap, id_meas);

    endfor
  endfor

endif


## [20] per x-section interface normal concentration profiles _avg
## shortcut to dyn statistics: interface normal profiles dyn analysis to
## analyze the time averaged peak cn shifted concentration profiles only
## ... using concentration fields only if gas-liquid interface is close to mean
if 1

  for i_M = it_M
    for i_X = it_X
      close all

      ## section analysis identifier
      sec_a_id = ["cp-avg-" c_if_method "-" "cn-" ap.c_method]

      id_meas = get_measid (ap.ids_C{i_C}, ap.ids_O{i_O}, ap.ids_A(i_A), ap.ids_T, ap.ids_L{i_L}, ap.ids_M(i_M), ap.ids_G, ap.ids_X(i_X), ap.ids_Z)

      ## path to store per section analysis results and plots to
      save_dir_id = [ap.result_dir id_meas "/" date_str "_" sec_a_id "/"];
      mkdir (save_dir_id);

      ## load normalized concentration field and interface data
      dirs = glob ([ap.result_dir id_meas "/" "*_cn-linear" "/"])
      files = glob ([dirs{end} "*.v7"]) # select latest result
      c_msh = cn_dyn = x = h_g_fit = []
      for i_f = 1 : numel (files)
        [~, fname, ~] = fileparts (files{i_f});
        switch (fname)
          case "cn_dyn"
            load ("-v7", files{i_f}, "c_msh", "cn_dyn", "cn_if_dyn");
          case "if_dyn_stat"
            load ("-v7", files{i_f}, "x", "h_g_fit");
          otherwise
            warning (["unexpected .v7 file: " fname]);
        endswitch
      endfor

      n_t = numel (cn_dyn) # number of time steps
      sf = get_sf (c_msh) # in mm
      n_p = numel (x) # number of concentration profiles
      profile_depth = ap.cp_pd_M(i_M) # in mm

      ## filter outlier interface height (film moving too much)
      [h_g_mean, ~] = calc_vec_avg_cells (h_g_fit, "mean");
      scmad = i_t_out = []
      for i_t = 1:n_t
        [~, ~, scmad(i_t)] = outlier_rm (h_g_fit{i_t}, h_g_mean);
      endfor

      ## gas-liquid interface outliers
      if (! (n_t == 10)) #
        lim_scmad = 1.25 * median (scmad);
      else
        lim_scmad = 10 * median (scmad);
      endif
      i_t_out = scmad > lim_scmad

      ##
      h_g_fit_valid = h_g_fit(! i_t_out);

      ## avg interface
      h_g_mean = calc_vec_avg_cells (h_g_fit_valid, "mean");

      fh = figure ();
      hold on
      plot (scmad, "kx")
      plot ([1 n_t], median (scmad) * [1 1], "--k;median;")
      plot ([1 n_t], lim_scmad * [1 1], "--r;limit;" )
      xlabel ("frame #")
      ylabel ("scmad")

      figure ()
      hold on;
      for i_t = 1:n_t
        plot (h_g_fit{i_t}, [";i_t = " num2str(i_t) ";"]);
      endfor
      plot (h_g_mean, "k", "linewidth", 2);

      ##
      ## per x-section interface normal concentration profile analysis
      ##

      ## concentration profile coordinates
      [snp, sf_p] = get_snp (profile_depth, sf, sn_idx_off=8);

      ## interpolate fields on interface normal line mesh
      switch (c_if_method)
        case {"calib"}
          cn_dyn_valid = cn_dyn (! i_t_out);
          [cp, p_msh, msh_n] = interp_if_norm (snp, h_g_fit_valid, c_msh, cn_dyn_valid);
        case {"calib-if"}
          cn_dyn_valid = cn_if_dyn (! i_t_out);
          [cp, p_msh, msh_n] = interp_if_norm (snp, h_g_fit_valid, c_msh, cn_dyn_valid);
      endswitch
      cn_dyn_valid_avg = calc_im_avg_cells (cn_dyn_valid, "median");

      ## test plot concentration field with interface normal lines
      fh = plot_cn_if_norm (x, h_g_fit{1}, msh_n, c_msh, cn_dyn{1}, ap.ids_C{i_C}, ap.ids_X(i_X));
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "if-normals"]);

      ## test plot interpolated concentration profiles
      fh = plot_cp_norm (p_msh, cp{1});
      title ("cp")
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "interp_if_norm"]);

      ## calculate bulk to interface concentration normalized profiles
      [cp_n, cp_b, cp_s] = cp_norm_field (snp, cp);

      ## test plot bulk and surface concentration
      fh = figure ();
      hold on;
      plot (x, cp_b{1}, "k;cn bulk;");
      plot (x, cp_s{1}, "b;cn surface;");
      xlabel ("x in mm");
      ylabel ("cn in -");
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_bulk_surface"]);

      ## test plot normalized concentration profiles
      fh = plot_cp_norm (p_msh, cp_n{1});
      title ("cp normalized")
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n"]);

      ## shift profile peak to snp = 0
      t1 = tic ();
      [snp_o, p_msh_o, cp_n_o] = cp_peak_shift (snp, sn_idx_off, p_msh, cp_n);
      toc (t1)

      ## test plot cp_n shifted
      fh = plot_cp_norm (p_msh_o, cp_n_o{1});
      title ("cp normalized, offset to max")
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n_o"]);

      ## compute temporal average concentration profile field
      cp_n_o_avg = calc_im_avg_cells (cp_n_o, "median");

      ## test plot averaged cp_n_o
      fh = plot_cp_norm (p_msh_o, cp_n_o_avg);
      title ("cp normalized, offset to max, median")
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n_o_avg"]);

      ##
      ## normalized concentration profile fit
      ##

      ## profile fit range parameters
      switch (i_M)
        case 1
          idx_r = [3 10];
          dcds_idx = 16;
          sig = 4;
        case 2
          idx_r = [2 8];
          dcds_idx = 16;
          sig = 4;
        case 3
          idx_r = [2 6];
          dcds_idx = 16;
          sig = 4;
        case 4
          idx_r = [2 4];
          dcds_idx = 16;
          sig = 4;
      endswitch

      ## erfc profile fit
      t_sgl = tic ();
      [a_fit, cp_nn] = erfc_profile_fit (snp_o, cp_n_o_avg, sig, idx_r, dcds_idx, [], false);
      toc (t_sgl)

      ## retrieve concentration boundary layer thickness
      delta_c = [];
      for i_p = 1:n_p
        delta_c(i_p) = cn_fit_delta_c (a_fit(i_p,:));
      endfor

      fh = figure ();
      hold on;
      plot (x, delta_c);
      plot (x, movmedian (delta_c, 41), "r");
      xlabel ("x in mm");
      ylabel ("delta_c in mm");
      ylim ([0 2.0*median(delta_c)])
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "delta_c"]);

      fh = plot_map_msh (p_msh_o, cp_nn);
      caxis ([0 1]);
      colorbar;
      hold on;
      plot3 (x, delta_c, ones (n_p, 1), "r");
      xlabel ("x in mm");
      ylabel ("s_n in mm");
      ylim ([0 profile_depth])
      set (gca (), "ydir", "reverse");
      legend ({"cp nn", "delta_c"});
      legend ("location", "southeast");
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "delta_c_cp_nn"]);

      ## store results
      cd (save_dir_id);
      save -v7 "dyn_cn_cp.v7" profile_depth sf_p snp p_msh cp_n cp_s cp_b i_t_out
      save -v7 "dyn_cn_cp_avg_fit.v7" snp_o p_msh_o cp_nn delta_c h_g_mean cp_n_o_avg a_fit # temporal avg result
      save -v7 "dyn_cn_cp_avg-cn.v7" c_msh cn_dyn_valid_avg
    endfor
  endfor

endif
## [30] x-stitch dyn avg data and store
it_X = 1 : numel (ap.ids_X) # all sections
if 1

  ## input data path of per section analysis
  sec_data_dir = [pdir.analyzed ap.a_type "/"]

  ## per section analysis method was
  a_id = ["cp-avg-" c_if_method "-" "cn-" ap.c_method]

  ## file holding the data to be assembled
  fn_data{1} = "dyn_cn_cp.v7"
  fn_data{2} = "dyn_cn_cp_avg_fit.v7"
  fn_data{3} = "dyn_cn_cp_avg-cn.v7"

  ## variables to be loaded and assembled. first one has to be the name of mesh variable
  var_list{1} = {"p_msh" "cp_n" "cp_s" "cp_b"}
  var_list{2} = {"p_msh_o" "cp_nn" "delta_c" "h_g_mean" "cp_n_o_avg" "a_fit"}
  var_list{3} = {"c_msh" "cn_dyn_valid_avg"}
  n_f = numel (var_list);

  ## per data variable: method for time series averaging before stitching
  ## ("median" or "mean")
  ## or set one for all
  avg_method = {"median"}

  ## resulting stitching analysis identifier
  stitch_a_id = ["x-stitch_" a_id]

  ##
  ## prepare inter section offset correction
  ##

  ## load needed section data
  for i_M = it_M
    for i_X = it_X
      dyncase = get_measid (ap.ids_C{i_C}, ap.ids_O{i_O}, ap.ids_A(i_A), ap.ids_T, ap.ids_L{i_L}, ap.ids_M(i_M), ap.ids_G, ap.ids_X(i_X), ap.ids_Z);
      ## per section analysis method was
      fn_M = glob ([pdir.analyzed ap.a_type "/" dyncase "/" "*_" ["cp-avg-" c_if_method "-" "cn-" ap.c_method] "/"]);
      fn_M = fn_M{end}; # use latest result

      ## {"p_msh_o" "cp_nn" "delta_c" "h_g_mean" "cp_n_o_avg" "a_fit"}
      load ([fn_M "/" fn_data{2}], "p_msh_o");
      load ([fn_M "/" fn_data{2}], "h_g_mean");
      x_MX{i_M,i_X} = p_msh_o{1}(1,:);
      h_MX{i_M,i_X} = h_g_mean;
    endfor
  endfor

  ## helper to find systematic offset
  xoff = []
  for i_M = it_M
    for i_X = it_X(1:end-1)
      x_l = x_MX{i_M,i_X};
      x_r = x_MX{i_M,i_X+1};
      xpos_l = ap.ids_X(i_X);
      xpos_r = ap.ids_X(i_X+1);
      delta_u_l = h_MX{i_M,i_X};
      delta_u_r = h_MX{i_M,i_X+1};
      xoff(i_M,i_X) = calc_xsec_if_offset_x (x_l, x_r, delta_u_l, delta_u_r, xpos_l, xpos_r);
    endfor
  endfor
  xoff

  ap.xoff = [+0.05 0 +0.2 +0.2] # mm
  yoff = zeros (numel(ap.ids_M), numel(ap.ids_X));
  yoff(1,4) = +0.005;
  yoff(2,1) = +0.02;
  yoff(2,4) = -0.005;
  yoff(3,3) = +0.015;
  yoff(3,4) = +0.015;
  yoff(4,1) = -0.015;
  yoff(4,3) = +0.01;
  yoff(4,4) = +0.005;

  for i_M = it_M
    ap.yoff = yoff(i_M,:);
    for i_f = 1:n_f
      ## prepare result storage path
      id_meas = get_measid (ap.ids_C{i_C}, ap.ids_O{i_O}, ap.ids_A(i_A), ap.ids_T, ap.ids_L{i_L}, ap.ids_M(i_M), ap.ids_G, ap.ids_X, ap.ids_Z)
      save_dir_id = [ap.result_dir id_meas "/" date_str "_" stitch_a_id "/"];
      mkdir (save_dir_id);

      ## list analysis results to combine
      for i_X = it_X
        dyncase = get_measid (ap.ids_C{i_C}, ap.ids_O{i_O}, ap.ids_A(i_A), ap.ids_T, ap.ids_L{i_L}, ap.ids_M(i_M), ap.ids_G, ap.ids_X(i_X), ap.ids_Z);
        fn_sec_data = glob ([sec_data_dir dyncase "/*_" a_id]);
        fn_sec_a{i_X} = fn_sec_data{end}; # use latest result
      endfor

      ## combine mesh and data of sections
      [msh_gl dat_gl] = stitch_msh_dat (fn_data{i_f}, fn_sec_a, var_list{i_f}, avg_method, ap);

      ## save data
      cd (save_dir_id)
      save ("-v7", ["x-stitch_msh_dat___" fn_data{i_f}], "msh_gl", "dat_gl")
    endfor
  endfor

endif
## [40] load dyn avg stitched data, filter, plot data overview and export
##       - export delta_c, snD and cp_nn
if 1

  ## per section analysis method was
  a_id = ["cp-avg-" c_if_method "-" "cn-" ap.c_method]

  ## resulting stitching analysis identifier
  stitch_a_id = ["x-stitch_" a_id]

  ##
  save_dir_id = [ap.result_dir stitch_a_id "/"];
  mkdir (save_dir_id);

  fn_data{1} = "dyn_cn_cp.v7"
  fn_data{2} = "dyn_cn_cp_avg_fit.v7"
  fn_data{3} = "dyn_cn_cp_avg-cn.v7"

  x_min = - 12.00 # mm
  x_max = + 20.00 # mm

  dat_gl = msh_gl = {};
  x_M = sf_M = {};
  cp_nn_M = delta_c_M = msh_p_M = {};
  cp_b_M = cp_s_M = {};
  cn_avg_M = msh_M = {};
  h_g_M = {};
  for i_M = it_M
    dyncase = get_measid (ap.ids_C{i_C}, ap.ids_O{i_O}, ap.ids_A(i_A), ap.ids_T, ap.ids_L{i_L}, ap.ids_M(i_M), ap.ids_G, ap.ids_X, ap.ids_Z);

    fn_M = glob ([ap.result_dir dyncase "/*_" stitch_a_id]);
    fn_M = fn_M{end}; # use latest result

    ## "p_msh" "cp_n" "cp_s" "cp_b"}
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{1}], "msh_gl");
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{1}], "dat_gl");

    x_M{i_M} = vec (msh_gl{1}(1,:));
    isin_x = (x_M{i_M} >= x_min) & (x_M{i_M} <= x_max);

    cp_b_M{i_M} = dat_gl.cp_b(isin_x);
    cp_s_M{i_M} = dat_gl.cp_s(isin_x);

    ## "p_msh_o" "cp_nn" "delta_c" "h_g_mean"
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{2}], "msh_gl");
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{2}], "dat_gl");

    msh_p_M{i_M} = msh_gl;
##    x_p_M{i_M} = x_M{i_M}; # x is the same
    sf_M{i_M} = get_sf (msh_p_M{i_M});

    cp_nn_M{i_M} = dat_gl.cp_nn;
    delta_c_M{i_M} = dat_gl.delta_c(isin_x);
    h_g_M{i_M} = dat_gl.h_g_mean;

    ## "c_msh" "cn_dyn_valid_avg"
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{3}], "msh_gl");
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{3}], "dat_gl");
    cn_avg_M{i_M} = dat_gl.cn_dyn_valid_avg;
    msh_M{i_M} = msh_gl;
    sf_M{i_M} = get_sf (msh_M{i_M});
  endfor

  x = x_M{it_M(1)}(isin_x);

##  i_M = 1;
  plot_map_msh (msh_M{i_M}, cn_avg_M{i_M}); hold on; plot (x_M{i_M}, h_g_M{i_M}, "r")

  ## filter for outliers
  delta_c_MM = delta_c_mov_M = {};
  for i_M = it_M
    delta_c_mov_M{i_M} = vec (movmedian (delta_c_M{i_M}, 81));
    delta_c_MM{i_M} = outlier_rm (delta_c_M{i_M}, delta_c_mov_M{i_M});
  endfor

  ## transform to diffusion front
  snD_M = snD_mov_M = {};
  for i_M = it_M
    snD_mov_M{i_M} = convert_deltac_snd (delta_c_mov_M{i_M});
    snD_M{i_M} = convert_deltac_snd (delta_c_MM{i_M});
  endfor

  ##
  ## print overview of data
  ##
  fh = figure ();
  hold on;
  for i_M = it_M
    plot (x, cp_b_M{i_M}, [";i_M = " num2str(i_M) ";"]);
    plot (x, cp_s_M{i_M}, [";i_M = " num2str(i_M) ";"]);
  endfor
  xlabel ("x in mm");
  ylabel ("cp in -");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "x_cp_b_s_M"]);

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (x, delta_c_MM{i_M}, [";i_M = " num2str(i_M) ";"]);
    plot (x, movmedian (delta_c_M{i_M}, 21), [";i_M = " num2str(i_M) ";"])
  endfor
  xlabel ("x in mm");
  ylabel ("delta_c in mm");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "x_delta_c_M"]);

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (x, 1e-6 * 0.5 * snD_M{i_M} .^ 2, [";i_M = " num2str(i_M) ";"]);
    plot (x, 1e-6 * 0.5 * snD_mov_M{i_M} .^ 2, [";i_M = " num2str(i_M) ";"]);
  endfor
  xlabel ("x in mm");
  ylabel ("snD^2 / 2 in m^2");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "x_snD_M"]);

  ##
  ## csv output delta_c vs. x, snD vs. x
  ##
  DC = cell2mat (delta_c_MM);
  DCM = cell2mat (delta_c_mov_M);
  fn_delta_c = [save_dir_id "x_deltac_M"]
  data = [x DC DCM];
  header = {"x in mm", "delta_c in mm for M1..M4", "delta_c smoothed in mm for M1..M4"};
  write_series_csv (fn_delta_c, data, header, []);
  ##
  SND = cell2mat (snD_M);
  SND_mov = cell2mat (snD_mov_M);
  fn_snd = [save_dir_id "x_snD_M"]
  data = [x SND SND_mov];
  header = {"x in mm", "sn_D in mm for M1..M4", "sn_D smoothed in mm for M1..M4"};
  write_series_csv (fn_snd, data, header, []);

  ## print cp_nn
  for i_M = it_M
    lim_x = [x_min x_max] # in mm
    lim_y = [0 0.1] # in mm
    lim_c = [0 1] # in -
    fn_cprint = [save_dir_id "cp_nn_M" num2str(i_M)]
    print_contour (fn_cprint, msh_p_M{i_M}{1}, max(lim_y) - msh_p_M{i_M}{2}, cp_nn_M{i_M}, lim_x, lim_y, sf_M{i_M}, lim_c, true);
  endfor

  ## print cn valid avg
  clims = [0.6 0.5 0.4 0.3]
  delta_u_avg = [0.42 0.58 0.74 0.93]
  yshift = -0.035
  for i_M = it_M
    mask = masking ("c", "gas", size(cn_avg_M{i_M}), min(min(msh_M{i_M}{2})), h_g_M{i_M}, get_sf(msh_M{i_M}), +10, 0.0);
    lim_x = [x_min x_max] # in mm
    lim_y = [0 2.5] # in mm
    lim_c = [0 clims(i_M)] # in -
    fn_cprint = [save_dir_id "cn_avg_M" num2str(i_M)]
    cprint = movmedian (cn_avg_M{i_M} .* mask, 3);
    print_contour (fn_cprint, msh_M{i_M}{1}, msh_M{i_M}{2}+yshift, cprint, lim_x, lim_y, sf_M{i_M}, lim_c, false);
  endfor

endif

## [50] comparison delta_c experimental vs. theoretical
if 1
  ## analysis per mass flow rate...##- delta_c (x)
  ## - sn_D^2 (t_c)
  ## - D_AB
  ## - INTEGRAL values ... what cn profile to use as input?!

  ## measured surface velocities from reference profile
  ref_prof = load ([pdir.analyzed "a_2DR10_reference_flow_profile/upstream/" "tab_meas_Re_deltau_us_deltac_mfr.txt"])
##  ref_prof = load ([pdir.analyzed "a_flat_reference_flow_profile/upstream/" "tab_meas_Re_deltau_us_deltac_mfr.txt"])

  u_s_meas = 1.02 * ref_prof.u_s

  T_K = ap.ids_T + 273.15;

  ## theoretical liquid mixture properties
  [w_lm, n_lm, rho_lm, eta_lm, ~, D_AB_lm] = get_fp_lm (pdir, ap.ids_L{1}, T_K);

  ## fluid properties experiment
  fp = get_fp_log (pdir, "2DR10_WG141");
  nu_exp = fp.eta / fp.rho

  nu_lm = eta_lm / rho_lm

  re_exp_M = nd_re_inlet (ref_prof.mfr_Nu, cell_width*1e-3, fp.eta)

  [delta_u, u_s, u_avg] = model_filmflow_laminar_u_profile_p (nu_exp, deg2rad(ap.ids_A), ref_prof.re_s);

  u_s = ref_prof.u_s

  ## contact time
  x_c = ((x + x_abs_meas)*1e-3 + x_off_inlet); # in m
  ## consider extra length of interface
  t_c = x_c ./ u_s;
  ## correct for curved interface: pseudo code: len_if (x) / mean (u_s(1:x))

  ## resulting contact time ranges
  x_sec(1,:)' ./ u_s

  ## diffusivity in m^2 / s
  D_AB = 0.807e-9
  ##D_AB = 0.486e-9

  x_eq = linspace (0, 0.12, 1000);

  figure ()
  hold on
  for i_M = it_M
    plot (x_c, 1e-3*delta_c_MM{i_M}, ["; delta_c meas i_M = " num2str(i_M) ";"]);
    plot (x_eq, model_filmflow_laminar_deltac (x_eq, D_AB_lm.PLIF2, u_s(i_M)), "--;delta_c eq;")
  endfor
  plot (x_sec + x_off_inlet, [0 0 0 0 0; 10e-5*[1 1 1 1 1]], "--k")
  xlabel ("x in m");
  ylabel ("delta_c in m");
  xlim ([0 0.1])
  ylim (1e-6 * [0 100])
  legend ("location", "northwest")

  figure ()
  hold on
  for i_M = it_M
    plot (t_c(:,i_M) , 1e-3*delta_c_MM{i_M}, ["k; delta_c meas i_M = " num2str(i_M) ";"]);
  endfor
  plot (x_eq ./ u_s(1), sqrt (pi * D_AB_lm.PLIF2*0.75 * (x_eq) ./ u_s(1)), "k--;eq 0.75*D;")
  plot (x_eq ./ u_s(1), - 5e-6 + sqrt (pi * D_AB_lm.PLIF2 * (x_eq) ./ u_s(1)), "-.b;eq;")
  plot (x_eq ./ u_s(1), sqrt (pi * D_AB_lm.PLIF2 * (x_eq) ./ u_s(1)), "k;eq;")
  plot (x_eq ./ u_s(1), + 5e-6 + sqrt (pi * D_AB_lm.PLIF2 * (x_eq) ./ u_s(1)), "-.b;eq;")
  plot (x_eq ./ u_s(1), sqrt (pi * D_AB_lm.PLIF2*1.25 * (x_eq) ./ u_s(1)), "k--;eq 1.25*D;")
  xlabel ("contact time s");
  ylabel ("delta_c in m");
  title ("D_AB = 0.807e-9");
  xlim ([0 1.25])
  ylim (1e-6 * [0 100])

##  ## compensate for delta_c overestimation
##  doff = 1e-3 * sf_M{i_M}(2) * [0 -1 -2 -3]
##
##  D_AB = D_AB_lm.PLIF2
##  figure ()
##  hold on
##  for i_M = it_M
##    plot (t_c(:,i_M) , 1e-3*delta_c_MM{i_M}, ["k; delta_c meas i_M = " num2str(i_M) ";"]);
##    plot (t_c(:,i_M) , doff(i_M) + 1e-3*delta_c_MM{i_M}, ["r; delta_c meas i_M = " num2str(i_M) ";"]);
##  endfor
##  plot (x_eq ./ u_s(1), sqrt (pi * D_AB*0.75 * (x_eq) ./ u_s(1)), "k--;eq 0.75*D;")
##  plot (x_eq ./ u_s(1), - 5e-6 + sqrt (pi * D_AB * (x_eq) ./ u_s(1)), "-.b;eq;")
##  plot (x_eq ./ u_s(1), sqrt (pi * D_AB * (x_eq) ./ u_s(1)), "k;eq;")
##  plot (x_eq ./ u_s(1), + 5e-6 + sqrt (pi * D_AB * (x_eq) ./ u_s(1)), "-.b;eq;")
##  plot (x_eq ./ u_s(1), sqrt (pi * D_AB*1.25 * (x_eq) ./ u_s(1)), "k--;eq 1.25*D;")
##  xlabel ("contact time s");
##  ylabel ("delta_c in m");
##  title ("D_AB = 0.807e-9");
##  xlim ([0 1.25])
##  ylim (1e-6 * [0 100])
##
##
##  figure ()
##  hold on
##  for i_M = it_M
##    plot (t_c(:,i_M), 1/2 * convert_deltac_snd (1e-3*delta_c_MM{i_M}) .^ 2, ["b; delta_c meas i_M = " num2str(i_M) ";"]);
##  endfor
##  D_AB = D_AB_lm.PLIF2
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, 0.75*D_AB, u_s(1))) .^ 2, "k--;eq 0.75*D;")
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, D_AB, u_s(1))).^ 2, "k;eq;")
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, 1.25*D_AB, u_s(1))).^ 2, "k--;eq 1.25*D;")
##  D_AB = D_AB_lm.PLIF1
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, 0.75*D_AB, u_s(1))) .^ 2, "r--;eq 0.75*D;")
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, D_AB, u_s(1))).^ 2, "r;eq;")
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, 1.25*D_AB, u_s(1))).^ 2, "r--;eq 1.25*D;")
##  xlabel ("contact time s");
##  ylabel ("snD^2 / 2 in m^2");
##  xlim ([0 1.25])
##  ylim (1e-9 * [0 1])
##


##
##  figure ()
##  hold on
##  for i_M = it_M
##    snd_sq = 1/2 * convert_deltac_snd (1e-3*delta_c_MM{i_M}) .^ 2;
##    snd_sq_mm = movmedian (snd_sq, 1111);
##    snd_sq_filter = outlier_rm (snd_sq, snd_sq_mm);
##    pfit{i_M} = polyfit (t_c(1:end,i_M), snd_sq_filter(1:end), 1)
##
##    plot (t_c(:,i_M), snd_sq, ["b; delta_c meas i_M = " num2str(i_M) ";"]);
##    plot (t_c(:,i_M), snd_sq_mm, ["g; delta_c meas i_M = " num2str(i_M) ";"]);
##
##    plot (t_c(:,i_M), polyval(pfit{i_M}, t_c(:,i_M)), ["-r; delta_c meas i_M = " num2str(i_M) ";"], "linewidth", 2);
##  endfor
##  D_AB = D_AB_lm.PLIF2
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, 0.75*D_AB, u_s(1))) .^ 2, "k--;eq 0.75*D;")
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, D_AB, u_s(1))).^ 2, "k;eq;")
##  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_deltac (x_eq, 1.25*D_AB, u_s(1))).^ 2, "k--;eq 1.25*D;")
##
##  D_fit = reshape (cell2mat (pfit), 2, 4)
##  D_fit = D_fit(1,:)
##  D_fit ./ D_AB_lm.PLIF2

  ## ! in contact time the measured surface velocity matters, probably lower than real
  ## ! in delta_c measurement probably overestimated for the very thin boundary layer (pixel blur i.e.)
  ## ! thus, rather compare in relative manner: inclination sn_D^2 / 2 vs. t_c

endif

