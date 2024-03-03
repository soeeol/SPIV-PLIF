##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## single frame dynamic PLIF analysis flat plate
## interface statistics and boundary layer thickness
##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ## select processed data to be analyzed
  aid.proc_type = "2d_dyn_Ic";
  aid.ids_L = {"WG141"};
  aid.ids_O = {"M13"};
  aid.ids_C = {"flat"};
  aid.ids_A = [60];
  aid.ids_M = [8 16 32 64];
  aid.ids_X = [-8 0 8 16];
  id_G = 2;
  id_T = 25
  id_Z = 0;

  a_type = "a_flat_dyn";
  c_method = "linear"; # "linear" "nonlinear" .. doesn't change delta_c
  c_if_method = "calib"; # "calib" "calib-if" .. high impact on delta_c

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
##  i_M = it_M = 2:4; # single mass flow rate
##  i_X = it_X = 2#[2 3 4]; # single scanned x sections

  ## set interface normal profile depth in mm
##  pd_M = [0.30 0.3 0.2 0.2]; # in mm, manually checked for bulk depth
##  pd_M = 0.1 + 3 * [0.051 0.037 0.029 0.022]; # in mm, setting three times the expected delta_c
  pd_M = [0.35 0.45 0.5 0.5];

  ## prepare directories
  result_dir = [pdir.analyzed a_type "/"]

endif

## [10] normalized concentration field per x-section
if 0

  for i_X = it_X
    for i_M = it_M

      close all

      ## section analysis identifier
      sec_a_id = ["cn-" c_method]

      ## load processed data (measured Ic(phi) images, PLIF data) per x section (aligned, interpolated on common grid)
      id_meas = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X(i_X), id_Z)
      filename = [pdir.processed id_meas "/*" aid.proc_type "*"];
      files = glob (filename);
      file = files{end} # select latest processing result
      [c_msh, c_dat, c_masks, c_h, ~, ~, ~, ~, pp] = load_uc (file, pdir.work);

      [c_h_g_mean, ~] = calc_vec_avg_cells (c_h.gas, "mean");

      x = c_msh{1}(1,:);
      x_abs = (x + aid.ids_X(i_X) + x_abs_meas) * 1e-3; # in m
      xmin = min (x);
      xmax = max (x);

      y = c_msh{2}(:,1);
      ymin = min (y);
      ymax = max (y);

      z = c_msh{3}(1,1)
      sf_c = get_sf (c_msh);
      n_t = numel (c_dat(1,:))

      dom = get_domain (pp);
      is_in_xdom = vec ((x >= dom.xmin) & (x <= dom.xmax));
      x_dom = x(is_in_xdom);

      ## path to store per section analysis results and plots to
      save_dir_id = [result_dir id_meas "/" date_str "_" sec_a_id "/"];
      mkdir (save_dir_id);

      ## per x-section calculate normalized concentration field
      run "a_flat_dyn_cn.m"

      ## store results
      cd (save_dir_id);
      save -v7 "cn_dyn.v7" c_msh cn_dyn cn_dyn_mean cn_if_dyn c_method
      save -v7 "if_dyn_stat.v7" x h_g_fit h_g_fit_mean x_dom if_vals max_dev mad_dev std_dev

    endfor
  endfor

endif

## [20] TODO: update script: per x-section dyn stat interface normal concentration profiles
if 0
##  run "a_flat_dyn_cprofile_fit.m"
endif

## [21] per x-section interface normal concentration profiles _avg
## shortcut to dyn statistics: interface normal profiles dyn analysis to
## analyze the time averaged peak cn shifted concentration profiles only
## ... using concentration fields only if gas-liquid interface is close to mean
if 0

  for i_M = it_M
    for i_X = it_X
      close all

      ## section analysis identifier
      sec_a_id = ["cp-avg-" c_if_method "-" "cn-" c_method]

      id_meas = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X(i_X), id_Z)

      ## path to store per section analysis results and plots to
      save_dir_id = [result_dir id_meas "/" date_str "_" sec_a_id "/"];
      mkdir (save_dir_id);

      ## load normalized concentration field and interface data
      dirs = glob ([result_dir id_meas "/" "*_cn-linear" "/"])
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
      sf_c = get_sf (c_msh) # in mm
      n_p = numel (x) # number of concentration profiles
      profile_depth = pd_M(i_M) # in mm

      ## filter outlier interface height (film moving too much)
      [h_g_mean, ~] = calc_vec_avg_cells (h_g_fit, "mean");

      scmad = i_t_out = []
      for i_t = 1:n_t
        [~, ~, scmad(i_t)] = outlier_rm (h_g_fit{i_t}, h_g_mean);
      endfor

      lim_scmad = 1.25 * median (scmad);
      i_t_out = scmad > lim_scmad;

      figure ();
      vnt = 1:n_t;
      hold on;
      for i_t = 1:n_t
        plot (h_g_fit{i_t},["; i_t = " num2str(i_t) " mad = " num2str(scmad(i_t)) ";"]);
      endfor
      plot (h_g_mean, "k;median;", "linewidth", 2)

      fh = figure ();
      hold on
      plot (scmad, "kx")
      plot ([1 n_t], median (scmad) * [1 1], "--k;median;")
      plot ([1 n_t], lim_scmad * [1 1], "--r;limit;" )
      xlabel ("frame #")
      ylabel ("scmad")

      ## gas-liquid interface outliers
      i_t_out = scmad > lim_scmad

      ## manual selection of outliers for measurements with slightly wobbly interface
      if ((i_M==1) && (i_X==1))
        i_t_out = ones (1, n_t);
##        i_t_out(1) = 0 # delta_c = 49 .. 43 µm
##        i_t_out(2) = 0 # delta_c = 53 .. 43 µm
##        i_t_out(3) = 0 # delta_c = 46 .. 53 µm
##        i_t_out(4) = 0 # delta_c = 49 .. 42 µm
##        i_t_out(5) = 0 # delta_c = 42 .. 52 µm
##        i_t_out(6) = 0 # delta_c = 45 .. 45 µm
##        i_t_out(7) = 0 # delta_c = 36 .. 48 µm
##        i_t_out(8) = 0 # delta_c = 36 .. 40 .. 55 µm
##        i_t_out(9) = 0 # delta_c = 43 .. 45 µm
##        i_t_out(10) = 0 # delta_c = 43 .. 37 .. 46 µm
##        i_t_out(11) = 0 # delta_c = 46 .. 40 .. 46 µm
##        i_t_out(12) = 0 # delta_c = 47 .. 38 .. 42 µm
##        i_t_out(13) = 0 # delta_c = 50 .. 45 .. 45 µm
##        i_t_out(14) = 0 # delta_c = 50 .. 40 µm
##        i_t_out(15) = 0 # delta_c = 50 .. 46 µm
##        i_t_out(16) = 0 # delta_c = 52 .. 42 µm
##        i_t_out(17) = 0 # delta_c = 40 .. 49 µm
##        i_t_out(18) = 0 # delta_c = 42 .. 49 µm
##        i_t_out(19) = 0 # delta_c = 38 .. 52 µm
##        i_t_out(20) = 0 # delta_c = 43 .. 50 µm
          ## combinations
##        i_t_out([7,8,9,10,12,13]) = 0 # same median deviation group, delta_c = 43 .. 38 .. 46 µm
        i_t_out([7,9,10,12,19]) = 0 # good, delta_c = 40 .. 39 .. 44 µm
      elseif ((i_M==2) && (i_X==4))
##        i_t_out = scmad > 0.9 * median (scmad) # 42 .. 49 µm
        i_t_out = ones (1, n_t);
##        i_t_out([9,10,14,15,17]) = 0 # delta_c = 42 .. 51 µm
##        i_t_out([1,5,16,20]) = 0 # delta_c = 40 .. 43 µm
##        i_t_out([2,4,8,11,12,13,18]) = 0 # delta_c = 38 .. 38 µm
        i_t_out([1,5,16,20,8,11,12,13,18]) = 0 # delta_c = 39 .. 42 µm
      endif

      plot (vnt(! i_t_out), scmad(! i_t_out), "go")
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "scmad"]);

      n_t = sum (! i_t_out)
      h_g_fit = h_g_fit (! i_t_out);
      ## avg interface
      h_g_mean = calc_vec_avg_cells (h_g_fit, "mean");

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
      [snp, sf_p] = get_snp (profile_depth, sf_c, sn_idx_off=4);

      ## interpolate fields on interface normal line mesh
      switch (c_if_method)
        case {"calib"}
          cn_dyn_valid = cn_dyn (! i_t_out);
          [cp, p_msh, msh_n] = interp_if_norm (snp, h_g_fit, c_msh, cn_dyn_valid);
        case {"calib-if"}
          cn_dyn_valid = cn_if_dyn (! i_t_out);
          [cp, p_msh, msh_n] = interp_if_norm (snp, h_g_fit, c_msh, cn_dyn_valid);
      endswitch
      cn_dyn_valid_avg = calc_im_avg_cells (cn_dyn_valid, "median");

      ## test plot concentration field with interface normal lines
      fh = plot_cn_if_norm (x, h_g_fit{1}, msh_n, c_msh, cn_dyn{1}, aid.ids_C{i_C});
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "if-normals"]);

      ## test plot interpolated concentration profiles
      fh = plot_cp_norm (p_msh, cp{i_t});
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
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n"]);

      ## shift profile peak to snp = 0
      [snp_o, p_msh_o, cp_n_o] = cp_peak_shift (snp, sn_idx_off, p_msh, cp_n);

      ## test plot cp_n shifted
      fh = plot_cp_norm (p_msh_o, cp_n_o{1});
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n_o"]);

      ## compute temporal average concentration profiles
      cp_n_o_avg = calc_im_avg_cells (cp_n_o, "median");
##      cp_n_o_avg = calc_im_avg_cells (cp_n_o, "mean");

      ## test plot averaged cp_n_o
      fh = plot_cp_norm (p_msh_o, cp_n_o_avg);
      print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n_o_avg"]);

      ##
      ## normalized concentration profile fit
      ##

      ## profile fit range parameters
      switch (i_M)
        case 1
          idx_r = [2 12];
          dcds_idx = 8;
          sig = 1;
        case 2
          idx_r = [2 8];
          dcds_idx = 8;
          sig = 1;
        case 3
          idx_r = [2 6];
          dcds_idx = 8;
          sig = 1;
        case 4
          idx_r = [1 4];
          dcds_idx = 6;
          sig = 1;
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
      plot (x, movmedian (delta_c, 21));
      xlabel ("x in mm");
      ylabel ("delta_c in µm");
      ylim ([0 1.5*median(delta_c)])
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
if 0

  ## input data path of per section analysis
  sec_data_dir = [pdir.analyzed "a_flat_dyn/"]

  ## per section analysis method was
  a_id = ["cp-avg-" c_if_method "-" "cn-" c_method]

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

  for i_M = it_M
    for i_f = 2#1:n_f
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
## [40] load dyn avg stitched data, filter, plot data overview and export
##       - export delta_c, snD and cp_nn
if 1

  ## per section analysis method was
  a_id = ["cp-avg-" c_if_method "-" "cn-" c_method]

  ## resulting stitching analysis identifier
  stitch_a_id = ["x-stitch_" a_id]

  ##
  save_dir_id = [result_dir stitch_a_id "/"];
  mkdir (save_dir_id);

  fn_data{1} = "dyn_cn_cp.v7"
  fn_data{2} = "dyn_cn_cp_avg_fit.v7"
  fn_data{3} = "dyn_cn_cp_avg-cn.v7"

  xmin = - 12.00 # mm
  xmax = + 20.00 # mm

  dat_gl = msh_gl = {};
  x_M = sf_M = {};
  cp_nn_M = delta_c_M = msh_p_M = {};
  cp_b_M = cp_s_M = {};
  cn_avg_M = msh_M = {};
  for i_M = it_M
    dyncase = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X, id_Z);

    fn_M = glob ([result_dir dyncase "/*_" stitch_a_id]);
    fn_M = fn_M{end}; # use latest result

    ## "p_msh" "cp_n" "cp_s" "cp_b"}
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{1}], "msh_gl");
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{1}], "dat_gl");

    x_M{i_M} = vec (msh_gl{1}(1,:));
    isin_x = (x_M{i_M} >= xmin) & (x_M{i_M} <= xmax);

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

    ## "c_msh" "cn_dyn_valid_avg"
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{3}], "msh_gl");
    load ([fn_M "/" "x-stitch_msh_dat___" fn_data{3}], "dat_gl");
    cn_avg_M{i_M} = dat_gl.cn_dyn_valid_avg;
    msh_M{i_M} = msh_gl;
    sf_M{i_M} = get_sf (msh_M{i_M});
  endfor

  x = x_M{1}(isin_x);

  ## filter for outliers
  delta_c_MM = delta_c_mov_M = {};
  for i_M = it_M
    delta_c_mov_M{i_M} = vec (movmedian (delta_c_M{i_M}, 41));
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
    lim_x = [xmin xmax] # in mm
    lim_y = [0 0.1] # in mm
    lim_c = [0 1] # in -
    fn_cprint = [save_dir_id "cp_nn_M" num2str(i_M)]
    print_contour (fn_cprint, msh_p_M{i_M}{1}, max(lim_y) - msh_p_M{i_M}{2}, cp_nn_M{i_M}, lim_x, lim_y, sf_M{i_M}, lim_c);
  endfor

  ## print cn valid avg
  clims = [0.6 0.5 0.4 0.3]
  delta_u_avg = [0.42 0.58 0.74 0.93]
  yshift = -0.035
  for i_M = it_M
    lim_x = [xmin xmax] # in mm
    lim_y = [0 2.5] # in mm
    lim_c = [0 clims(i_M)] # in -
    fn_cprint = [save_dir_id "cn_avg_M" num2str(i_M)]
    cprint = cn_avg_M{i_M};
    idx = msh_M{i_M}{2} > 1.25 * delta_u_avg(i_M);
    cprint(idx) = 0.0;
    print_contour (fn_cprint, msh_M{i_M}{1}, msh_M{i_M}{2}+yshift, cprint, lim_x, lim_y, sf_M{i_M}, lim_c);
  endfor

endif

## [50] comparison delta_c experimental vs. theoretical
if 0
  ## analysis per mass flow rate...##- delta_c (x)
  ## - sn_D^2 (t_c)
  ## - D_AB
  ## - INTEGRAL values ... what cn profile to use as input?!

  ## measured surface velocities from reference profile
  ref_prof = load ([pdir.analyzed "a_flat_reference_flow_profile/upstream/" "tab_meas_Re_deltau_us_deltac_mfr.txt"])

  u_s_meas = 1.02 * ref_prof.u_s

  T_K = id_T + 273.15;

  ## theoretical liquid mixture properties
  [w_lm, n_lm, rho_lm, eta_lm, ~, D_AB_lm] = get_fp_lm (pdir, aid.ids_L{1}, T_K);

  ## fluid properties experiment
  fp = get_fp_log (pdir, "flat_WG141");
  nu_exp = fp.eta / fp.rho

  nu_lm = eta_lm / rho_lm

  re_exp_M = nd_re_inlet (ref_prof.mfr_Nu, cell_width*1e-3, fp.eta)

  [delta_u, u_s, u_avg] = model_filmflow_laminar_u_profile_p (nu_exp, deg2rad(aid.ids_A), ref_prof.re_fp);

  u_s = ref_prof.u_s

  ## contact time
  x_off_inlet = 0.00 - 0.0050 + 0.0070 + 0.0075 # in m; longer contact than x_abs_meas with gas in inflow section
  x_c = ((x + x_abs_meas)*1e-3 + x_off_inlet); # in m
  t_c = x_c ./ u_s;

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
    plot (x_eq, model_filmflow_laminar_delta_x (x_eq, D_AB_lm.PLIF2, u_s(i_M)), ";delta_c eq;")
  endfor
  xlabel ("x in m");
  ylabel ("delta_c in m");
  xlim ([0 0.1])
  ylim (1e-6 * [0 100])

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

  ## compensate for delta_c overestimation
  doff = 1e-3 * sf_M{i_M}(2) * [0 -1 -2 -3]

  D_AB = D_AB_lm.PLIF2
  figure ()
  hold on
  for i_M = it_M
    plot (t_c(:,i_M) , 1e-3*delta_c_MM{i_M}, ["k; delta_c meas i_M = " num2str(i_M) ";"]);
    plot (t_c(:,i_M) , doff(i_M) + 1e-3*delta_c_MM{i_M}, ["r; delta_c meas i_M = " num2str(i_M) ";"]);
  endfor
  plot (x_eq ./ u_s(1), sqrt (pi * D_AB*0.75 * (x_eq) ./ u_s(1)), "k--;eq 0.75*D;")
  plot (x_eq ./ u_s(1), - 5e-6 + sqrt (pi * D_AB * (x_eq) ./ u_s(1)), "-.b;eq;")
  plot (x_eq ./ u_s(1), sqrt (pi * D_AB * (x_eq) ./ u_s(1)), "k;eq;")
  plot (x_eq ./ u_s(1), + 5e-6 + sqrt (pi * D_AB * (x_eq) ./ u_s(1)), "-.b;eq;")
  plot (x_eq ./ u_s(1), sqrt (pi * D_AB*1.25 * (x_eq) ./ u_s(1)), "k--;eq 1.25*D;")
  xlabel ("contact time s");
  ylabel ("delta_c in m");
  title ("D_AB = 0.807e-9");
  xlim ([0 1.25])
  ylim (1e-6 * [0 100])


  figure ()
  hold on
  for i_M = it_M
    plot (t_c(:,i_M), 1/2 * convert_deltac_snd (1e-3*delta_c_MM{i_M}) .^ 2, ["b; delta_c meas i_M = " num2str(i_M) ";"]);
  endfor
  D_AB = D_AB_lm.PLIF2
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, 0.75*D_AB, u_s(1))) .^ 2, "k--;eq 0.75*D;")
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, D_AB, u_s(1))).^ 2, "k;eq;")
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, 1.25*D_AB, u_s(1))).^ 2, "k--;eq 1.25*D;")
  D_AB = D_AB_lm.PLIF1
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, 0.75*D_AB, u_s(1))) .^ 2, "r--;eq 0.75*D;")
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, D_AB, u_s(1))).^ 2, "r;eq;")
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, 1.25*D_AB, u_s(1))).^ 2, "r--;eq 1.25*D;")
  xlabel ("contact time s");
  ylabel ("snD^2 / 2 in m^2");
  xlim ([0 1.25])
  ylim (1e-9 * [0 1])




  figure ()
  hold on
  for i_M = it_M
    snd_sq = 1/2 * convert_deltac_snd (1e-3*delta_c_MM{i_M}) .^ 2;
    snd_sq_mm = movmedian (snd_sq, 1111);
    snd_sq_filter = outlier_rm (snd_sq, snd_sq_mm);
    pfit{i_M} = polyfit (t_c(1:end,i_M), snd_sq_filter(1:end), 1)

    plot (t_c(:,i_M), snd_sq, ["b; delta_c meas i_M = " num2str(i_M) ";"]);
    plot (t_c(:,i_M), snd_sq_mm, ["g; delta_c meas i_M = " num2str(i_M) ";"]);

    plot (t_c(:,i_M), polyval(pfit{i_M}, t_c(:,i_M)), ["-r; delta_c meas i_M = " num2str(i_M) ";"], "linewidth", 2);
  endfor
  D_AB = D_AB_lm.PLIF2
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, 0.75*D_AB, u_s(1))) .^ 2, "k--;eq 0.75*D;")
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, D_AB, u_s(1))).^ 2, "k;eq;")
  plot (x_eq ./ u_s(1), 1/2 * convert_deltac_snd (model_filmflow_laminar_delta_x (x_eq, 1.25*D_AB, u_s(1))).^ 2, "k--;eq 1.25*D;")

  D_fit = reshape (cell2mat (pfit), 2, 4)
  D_fit = D_fit(1,:)
  D_fit ./ D_AB_lm.PLIF2

  ## ! in contact time the measured surface velocity matters, probably lower than real
  ## ! in delta_c measurement probably overestimated for the very thin boundary layer (pixel blur i.e.)
  ## ! thus, rather compare in relative manner: inclination sn_D^2 / 2 vs. t_c

  ## TODO: output for plotting

  ## contour plots of measured cn_dyn
##  lim_x_out = [-12 20];
##  imsc = 1; # print img px per meas px
##  xsc = 2; # aspect ratio
##  sf_p = 2.5e-3;
##  [XI, YI] = meshgrid ([lim_x_out(1):sf_p(1)*xsc/imsc:lim_x_out(2)], [0:sf_p(1)/imsc:0.1]);
##  for i_M = it_M
##    cprint = interp2 (msh_gl{i_M}{1}, msh_gl{i_M}{2}, cnn_gl{i_M}{1}, XI, YI, "pchip");
##    cprint(cprint<0) = 0.0;
##    cprint(cprint>1) = 1.0;
##    imwrite (ind2rgb(gray2ind(cprint), colormap("viridis")), [save_dir_p "/M" num2str(i_M) "_cn_norm.png"])
##    close all
##  endfor
##  ## measured dyn avg delta plus theoretical delta with D_AB.PLIF2 assumed
##  load ("-text", [pdir.analyzed "a_flat_x_stitch" "/inlet_local_Re.txt"]); # u_s_meas
##  for i_M = it_M
##    delta_dyn_mean_out = interp1 (msh_gl{i_M}{1}(1,:), delta_gl_mean{i_M}{1}, XI(1,:));
##    delta_dyn_std_out = interp1 (msh_gl{i_M}{1}(1,:), delta_gl_std{i_M}{1}, XI(1,:));
##    delta_eq = model_filmflow_laminar_delta_x ((XI(1,:)'+60)*1e-3, D_AB.PLIF2, u_s_meas(i_M));
##    plt_1_h = {mfilename, date, "", ""; "x in mm", "x_abs in m", "delta_c in m", "STD(delta_c) in m"};
##    cell2csv ([result_dir "/M" num2str(i_M) "_delta_dyn_mean.csv"], plt_1_h)
##    plt_1_d = [XI(1,:)' (XI(1,:)'+60)*1e-3 delta_dyn_mean_out'*1e3 delta_dyn_std_out'*1e3 delta_eq*1e3];
##    csvwrite ([result_dir "/M" num2str(i_M) "_delta_dyn_mean.csv"], plt_1_d, "append", "off", "precision", "%.4e")
##  endfor
endif

## fluid properties for liquid mixture
if 0

  ## mass fraction of actual liquid mixture in the tank during experiment
  ##
  ## logged fluid properties on 02.12.2021 (tracer + seeding inside)
  ## probe from 20 feeding liter tank
  T_test = 298.15;
  eta_test = 8.421e-3;
  nref_test = 1.41065;
  rho_test = 1147.7;
  ## mass fraction from test
  fname = "glycerol-water"; ext = [];
  n_PDMS = ri_PDMS_T (T_test, calib_w.fit_n_PDMS.c);
  w_exp = ri_matching_mf (nref_test, T_test, calib_w.fit_n_PT.c);
  ## compare to tabulated fp
  rho_exp = get_fp_tab (pdir, fname, pname="rho", T_test, w_exp, ext);
  eta_exp = get_fp_tab (pdir, fname, pname="eta", T_test, w_exp, ext);
  ##
  [~, ~, ~, ~, ~, D_AB] = get_fp_lm (pdir, aid.ids_L{i_L}, id_T+273.15);

endif

## (_) diffusivity estimate
if 0
  ## theoretical values
  run "a_flat_laminar_model.m"
  ##
  x_abs = (msh_gl{i_M}{1}(1,:) + 60) / 1e3;
  ##
  figure (); hold on
  ## all valid sections
  delta_c_nt = [delta_gl_mean{4}{1} delta_gl_mean{3}{1} delta_gl_mean{2}{1} delta_gl_mean{1}{1}(1:1700)];
  tc_all = [x_abs/u_s_meas(4) x_abs/u_s_meas(3) x_abs/u_s_meas(2) x_abs(1:1700)/u_s_meas(1)];
  ##delta_c_nt = [delta_gl_mean{4}{1} delta_gl_mean{3}{1} delta_gl_mean{1}{1}(1:1200)];
  ##tc_all = [x_abs/u_s_meas(4) x_abs/u_s_meas(3) x_abs(1:1200)/u_s_meas(1)];
  delta_c_nt = outlier_rm (delta_c_nt, movmedian(delta_c_nt,81));
  ##
  p_D_fit_all = polyfit (tc_all, (delta_c_nt*(sqrt (2 / pi))).^2/2, 1)
  plot ([0 1], polyval (p_D_fit_all, [0 1])-0*p_D_fit_all(2), "r;allall;")
  plot ([0 1], polyval (p_D_fit_all, [0 1])-1*p_D_fit_all(2), "r;allall;")
  ##
  delta_c_cmp = model_filmflow_laminar_delta_x (x_eq, D_eq, u_s_meas);
  for i_M = it_M
    plot (x_eq / u_s_meas(i_M), (delta_c_cmp(i_M,:)*(sqrt (2 / pi))).^2 / 2,  "m-;eq. abs.;")
  endfor
  plot ([0 1], polyval ([D_AB.PLIF1 0], [0 1]), "--;D_1 PLIF;")
  ##plot ([0 1], polyval ([D_AB.PLIF1 13e-11], [0 1]), "--;D_2 PLIF;")
  plot ([0 1], polyval ([D_AB.PLIF2 0], [0 1]), "--;D_2 PLIF;")
  ##plot ([0 1], polyval ([D_APLIF2F2 13e-11], [0 1]), "--;D_1 PLIF;")
  ##
  plot (tc_all, (delta_c_nt*sqrt(2/pi)).^2/2, "x")

  for i_M = it_M

    plot (x_eq / u_s_meas(i_M), (delta_c_cmp(i_M,:)*(sqrt (2 / pi))).^2 / 2,  "m-;eq. abs.;")
    tc = x_abs / u_s_meas(i_M);
  ##  switch i_M
  ##    case
  ##  endswitch
  ##  p_D_fit_1 = polyfit (tc, (delta_gl_mean{i_M}{1}*(sqrt (2 / pi))).^2/2, 1);
  ##  plot ([0 1], polyval (p_D_fit_1, [0 1]) - p_D_fit_all(2), "g;fit all;")
  ##  plot (tc, ((delta_gl_mean{i_M}{1}+1*-i_M*2.5e-6)*(sqrt (2 / pi))).^2/2, "x")
    plot (tc, ((delta_gl_mean{i_M}{1})*(sqrt (2 / pi))).^2/2-1*p_D_fit_all(2), "x")

    legend ("autoupdate", "off")
    plot (x_sec./ u_s_meas(i_M), [0 0 0 0 0; 1e-9*[1 1 1 1 1]], "k")
  endfor
  xlim ([0 0.6])
  ylim ([0 8e-10])

endif


