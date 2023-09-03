##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## single frame dynamic PLIF analysis flat plate
## interface statistics and boundary layer thickness
##
## per x-section alignment and analysis, then x-stitching
##
## Author: Sören J. Gerke
##

## (0) init
if 1
  proc_type = "2d_dyn_Ic";
  f_Hz = 10; # aquisition frequency
  testplots = testplots_fit = 0;
  ## select processed data to be part of this analysis
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
  ##
  a_type = "a_flat_dyn";
  c_method = "linear"; # "linear" "nonlin"
  c_if_method = "calib"; # "calib" "calib-if"
  ## iterators
  it_A = 1:numel(aid.ids_A); # angles
  it_C = 1:numel(aid.ids_C); # cells
##  it_M = 1:numel(aid.ids_M); # mass flow rate
  it_M = 1; # mass flow rate
  it_X = 1:numel(aid.ids_X); # scanned x sections
  it_X = 1; # scanned x sections
  ## fixed
  i_L = i_O = 1; # liquid, optical setup
  i_A = 1;
  i_C = 1;
  i_M = 1;
  i_X = 1;

  ## set interface normal profile depth in mm
  pd_M = [0.3 0.35 0.4 0.5];

  ## prepare directories
  save_dir = [pdir.analyzed a_type "/"]
  save_dir_p = [pdir.plot a_type "/"]
  mkdir (save_dir_p)
  ##
  for i_M = it_M
    measid_stitch{i_M} = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, [], id_Z);
    save_dirs{i_M} = [save_dir "c-" c_method "_" c_if_method "_" measid_stitch{i_M} "/"];
    mkdir (save_dirs{i_M});
  endfor
endif

## (1) load processed data (aligned, interpolated on common grid)
if 1
  id = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X(i_X), id_Z)
  filename = [pdir.processed id "/*" aid.proc_type "*"];
  files = glob (filename);
  file = files{end} # select newest processing result
  [c_msh, c_dat, c_masks, c_h, u_msh, u_dat, u_masks, u_h, pp] = load_uc (file, pdir.work);
  x = c_msh{1}(1,:);
  x_abs = ( x + aid.ids_X(i_X) + x_abs_meas ) * 1e-3;
  xmin = min (x);
  xmax = max (x);
  y = c_msh{2}(:,1);
  ymin = min (y);
  ymax = max (y);
  z = 0;
  sf_c = get_sf (c_msh);
  sf_u = get_sf (u_msh);
  n_t = numel (c_dat(1,:))
##  n_t = 1
endif

## (2) per section alignment of Ic recordings; concentration field, interface detection and statistics
if 1
  ## correct offsets betweens Ic recordings (per x-section)
  xoff = yoff = zeros (numel(it_X), 1);
  switch aid.ids_M(i_M)
    case {64,32,16,8} # same systematic x offset vs. Ic1 found. (spatial calibration coordinate system offset between Ic records?)
      xoff(1,1) = -50; # idx
  endswitch
  mxyoff = [xoff yoff];
  ##
  p_dat.c_msh{i_A, i_C, i_M, i_X} = c_msh;
  p_dat.c_dat{i_A, i_C, i_M, i_X} = {c_dat{1,1} c_dat{2,1} c_dat{3,1}};
  p_dat.c_h{i_A, i_C, i_M, i_X}.gas = c_h.gas{1};
  p_dat.c_h{i_A, i_C, i_M, i_X}.wall = c_h.wall;
  ## correct for systematic offset
  p_dat = corr_xy_min_offsets (p_dat, a_type, it_A, it_C, it_M, it_X, mxyoff);
  ## correct slight y offset between Ic recordings (per x-section)
  for i_t = 1:n_t
    pdat.c_dat{i_A, i_C, i_M, i_X}{1} = c_dat{1,i_t};
    pdat_n{i_t} = corr_y_offsets (p_dat, it_A, it_C, it_M, it_X);
  endfor

  ## TODO: adjust mesh

  ## concentration field
  ## linear 2 point ref calib, with intensity reduction towards interface in calib
  cn = cell (1, n_t);
  for i_t = 1:n_t
    cn{i_t} = calc_cn ({c_dat{1,i_t} c_dat{2,1} c_dat{3,1}}, [0 1], c_method, sig=2, testplots=false);
  endfor
  ## time average
  cn_mean = 0;
  for i_t = 1:n_t # t
    cn_mean = cn_mean + cn{i_t};
  endfor
  cn_mean = cn_mean / n_t;

  if testplots
    ## fluorescence maps
    plot_map (c_dat{2,1}); title ("Ic0");
    plot_map (c_dat{3,1}); title ("Ic1");
    fh1 = figure ()
    for i_t = 1:n_t
      clf (fh1);
      surf (c_dat{1,i});
      shading flat;
      view ([0 0 1]);
      colorbar;
      title (["Ic " num2str(i_t)]);
      pause (0.1);
    endfor
    ## concentration field
    fh2 = figure ()
    for i_t = 1:n_t
      clf (fh2);
      surf (c_msh{1}, c_msh{2}, c_msh{3}, cn{i_t});
      shading flat;
      view ([0 0 1]);
      colorbar;
      caxis ([0 1])
      title (["cn " num2str(i_t)]);
      pause (0.1);
    endfor
  endif

  ## interface detection
  ##
  ## uncomment next line to click select upstream inital interface start
  ##  pp_stitch.y0_if_c.data = [];
  [ifg_meas, ifg_meas_idx, ispeak] = interface_xy (c_msh, cn_mean, tol=10, "max", pp, []);
  ##
  fh3 = check_if_plot (ifg_meas, ispeak, c_msh, cn_mean);
  h_c_g_mean = ifg_meas(:,2);

  if testplots
    ## check cn y-profile
    x_idx = round (rand(1) * numel(cn_mean(1,:)))
    figure (); hold on;
    plot (cn{1}(:,x_idx), "b")
    plot (cn_mean(:,x_idx), "k")
    plot (ifg_meas_idx(x_idx,2)*[1 1], [0 1], "--")
  endif

  ## smooth spline fit representation of interface
  sps = 1;
  spf = splinefit (x(ispeak==1), h_c_g_mean(ispeak==1), round(numel(it_X)*sps*1), "order", 3, "beta", 0.75);
  h_c_g_fit_mean = ppval (spf, x);

  ## interface fit check
  fh = plot_map_msh (c_msh, cn_mean);
  hold on;
  draw_cell (aid.ids_C{i_C}, [], 1)
  plot (x, h_c_g_mean, "-k");
  plot (x, h_c_g_fit_mean, "b-", "linewidth", 2);
  plot (x, c_h.wall, "r-");
  xlabel ("x in mm")
  ylabel ("y in mm")
  xlim ([-4 4])
  ylim ([max(h_c_g_fit_mean)-0.5 1.1*max(h_c_g_fit_mean)])
##  axis image

  xy_if = xy_idx = ispeak = h_c_g_fit = cell (1, n_t);
  msg = {"interface detection ok?", "adjust tol, starting point or method"};
  for i_t = 1:n_t
    [xy_if{i_t}, xy_idx{i_t}, ispeak{i_t}] = interface_xy (c_msh, cn{i_t}, tol*2, "max", pp, ifg_meas_idx(:,2));
    printf ([">>> if detection in " num2str(i_t) " of " num2str(n_t) " c records \n"]);
    if (i==1)
      fh3 = check_if_plot (xy_if{i_t}, ispeak{i_t}, c_msh, cn{i_t});
      if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
        close (fh3)
      else
        error (msg{2});
      endif
    endif
  endfor

  ## smooth spline fit representation of interface
  for i_t = 1:n_t
    spf = splinefit (x(ispeak{i_t}==1), xy_if{i_t}(ispeak{i_t}==1,2), round(numel(it_X)*sps*1), "order", 3, "beta", 0.75);
    h_c_g_fit{i_t} = ppval (spf, x);
  endfor

  ##
  fh = figure ();
  for i_t = 1:n_t
    clf (fh);
    surf (c_msh{1}, c_msh{2}, cn{i_t});
    shading flat;
    view ([0 0 1]);
    colorbar;
    hold on;
    plot3 (x, h_c_g_fit{i_t}, ones(1,numel(x)), "-m");
    title (["cn t#" num2str(i_t)]);
    xlabel ("x in mm")
    ylabel ("y in mm")
    ylim ([max(h_c_g_fit_mean)-0.5 1.1*max(h_c_g_fit_mean)])
    zlim ([-1 1])
##    pause (0.1);
    print (fh, "-djpeg", "-color", ["-r" num2str(250)], [save_dir_p "/cn_t_" num2str(i_t) ".jpg"]);
  endfor

  x_idx = round (rand(1) * numel(cn_mean(1,:)))
  fh = figure ();
  for i_t = 1:n_t
    clf (fh);
    hold on;
    plot (y, cn{i_t}(:,x_idx), "b")
    plot (y,cn_mean(:,x_idx), "k")
    plot (h_c_g_fit{i_t}(x_idx)*[1 1], [0 1], "b-")
    plot (h_c_g_fit_mean(x_idx)*[1 1], [0 1], "k-")
    xlim ([max(h_c_g_fit_mean)-0.5 1.1*max(h_c_g_fit_mean)])
    xlabel ("y in mm")
    ylabel ("cn in -")
    title (["cn t#" num2str(i_t) " (for one random pixel row)"]);
##    pause (0.1);
    print (fh, "-djpeg", "-color", ["-r" num2str(250)], [save_dir_p "/cn_prof_t_" num2str(i_t) ".jpg"]);
  endfor

  ##
  ## interface statistics
  ##

  ## statistics
  vals = vals_fit = [];
  for i = 1:numel(h_c_g_fit{1})
    for j = 1:n_t
  ##    vals_fit(i,j) = h_c_g_fit{j}(i);
      vals(i,j) = xy_if{j}(i,2);
    endfor
  endfor

  max_dev = max (abs (vals - mean (vals,  2)), [], 2);
  mad_dev = mad (vals, 0, 2); # 0:mean 1:median
  std_dev = std (vals, [], 2);
  mmax_dev = median (max_dev)
  mmad_dev = median (mad_dev)
  mstd_dev = median (std_dev)

  ## output interface statistics
  fh4 = figure ();
  hold on;
  plot (x, max_dev, "k;max abs deviation;");
  plot (x, std_dev, "r;std deviation;");
  plot (x, mad_dev, "b;meam abs deviation;");
  lh = legend ("autoupdate", "off");
  plot ([-4 4], [1 1 ] .* mmax_dev, "k-.");
  plot ([-4 4], [1 1 ] .* mmad_dev, "b-.");
  plot ([-4 4], [1 1 ] .* mstd_dev, "r-.");
  hold off;
  title ("deviation from mean interface position")
  xlabel ("x in mm")
  ylabel ("dispersion measures in mm")
  print (fh4, "-djpeg", "-color", ["-r" num2str(250)], [save_dir_p "/if_stat.jpg"]);

  cd (save_dir);
  save -v7 "if_stat.v7" n_t x c_h mmax_dev mmad_dev mstd_dev max_dev mad_dev std_dev
endif

##
## interface normal concentration profile analysis
##
if 1
  sf = sf_c;
  sf_p = 0.5 * sf(1); # resolution along profile line
  px = x';
  sn_idx_off = 4;
  profile_depth = pd_M(i_M)
  ## concentration profile coordinates
  snp = sf_p * (-sn_idx_off:numel(linspace(0,profile_depth,round(profile_depth/sf_p+1)))-1) * 1e-3; # m
  ##
  p_msh = mesh_uc ([0 numel(snp)-1], [0 numel(px)-1], 0, xmin, [], "c", [sf_p sf(1)]);
  cp = cell (1, n_t);
  for i_t = 1:n_t
    printf ([">>> interpolation on interface normal lines: " num2str(i_t) " of " num2str(n_t) " records \n"]);
    ## interface normal line mesh
    h_g = h_c_g_fit{i_t};
    py = h_g';
    ap = angle2Points ([px(1:end-1) py(1:end-1)], [px(2:end) py(2:end)]);
    ## concentration profile lines normal to interface
    nlines = createLine ([px py], pi/2+[ap; ap(end)] .*ones (numel(px),1));
    ## create the line meshes
    msh_n_x = msh_n_y = zeros (length(nlines), numel(snp));
    for i_l = 1:length (nlines)
      edge = createEdge (nlines(i_l,:) .* ones(numel(snp), 1), -1e3*snp');
      msh_n_x(i_l,:) = edge(:,3);
      msh_n_y(i_l,:) = edge(:,4);
    endfor
    msh_n = {msh_n_x, msh_n_y, zeros(size(msh_n_x))};
    ## map with some profile lines - test line meshes
    if (i_t == 1)
      fh = plot_map_msh (c_msh, cn{1});
      hold on
      draw_cell (pp.cell.data, [], 1)
      for i_l = 1:50:size(msh_n_x,1)
        plot (msh_n_x(i_l,:), msh_n_y(i_l,:), "m")
      endfor
      axis image
      xlim ([-4 4])
      ylim ([0 1.1*max(h_c_g_fit_mean)])
      plot (x, h_g, "r")
      print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "profile-lines-test"]);
    endif
    ## interpolate fields on line meshes
    switch c_if_method
      case {"calib"}
        cp{i_t} = interp2 (c_msh{1}, c_msh{2}, cn{i_t}, msh_n_x, msh_n_y, "pchip", 0.0);
      case {"calib-if"}
        cp{i_t} = interp2 (c_msh{1}, c_msh{2}, cn_if{i_t}, msh_n_x, msh_n_y, "pchip", 0.0);
    endswitch
    if (i_t == 1)
      plot_map (cp{i_t});
    endif
  endfor

  ## per profile: measured bulk concentration, interface concentration, normalization
  cp_b = cp_s = cp_n = cell (1, n_t);
  for i_t = 1:n_t
##    cp_mm_b = movmean (cp{i_t}, 21, "Endpoints", 0.0);
    cp_b{i_t} = cp_s{i_t} = sn_max{i_t} = zeros (1, size(cp{i_t},1));
    for i_p = 1:size(cp{i_t},1)
      cp_b{i_t}(i_p) = mean (cp{i_t}(i_p,end-50:end));
##      cp_b{i_t}(i_p) = 0;
##      cp_s{i_t}(i_p) = cp(i_p,snp==0);
      [cp_s{i_t}(i_p), sn_max{i_t}(i_p)] = max (cp{i_t}(i_p,1:20));
    endfor
    cp_s_mm = movmedian (cp_s{i_t}, 81);
    [cp_s_r{i_t}, isout] = outlier_rm (cp_s{i_t}, cp_s_mm);
    if (i_t == 1)
      fh = figure (); hold on;
      plot (x, cp_s{i_t}, "k");
      plot (x, cp_s_r{i_t}, "g");
      plot (x, cp_s_mm, "b");
      plot (x, movmean (cp_s_r{i_t}, 21), "m");
      plot (x(isout==1), cp_s{i_t}(isout==1), "*r");
      xlabel ("x")
      ylabel ("surface concentration")
    endif
    cp_s{i_t} = cp_s_r{i_t};
    ##
    cp_b_mm = movmedian (cp_b{i_t}, 81);
    [cp_b_r{i_t}, isout] = outlier_rm (cp_b{i_t}, cp_b_mm);
    if (i_t == 1)
      fh = figure (); hold on;
      plot (x, cp_b{i_t}, "k");
      plot (x, cp_b_r{i_t}, "g");
      plot (x, cp_b_mm, "b");
      plot (x, movmean (cp_b_r{i_t}, 21), "m");
      plot (x(isout==1), cp_b{i_t}(isout==1), "*r");
      xlabel ("x")
      ylabel ("bulk concentration")
    endif

    ## normalize from measured cn_bulk to cn_interface
    cp_n{i_t} = [];
    for i_p = 1:size(cp{i_t},1)
      cp_n{i_t}(i_p,:) = (cp{i_t}(i_p,:) - cp_b{i_t}(i_p)) / (cp_s{i_t}(i_p) - cp_b{i_t}(i_p));
    endfor
    ## shift profiles to max
    for i_p = 1:size(cp_n{i_t},1)
      cp_n{i_t}(i_p,1:end-sn_max{i_t}(i_p)) = cp_n{i_t}(i_p,sn_max{i_t}(i_p):end-1);
    endfor
    cp_n{i_t} = movmedian (cp_n{i_t}, 21, 1, "Endpoints", 0.0);

    snp = snp + sn_idx_off*sf_p*1e-3;
    sn_idx_off = 0;

    if (testplots && (i_t ==1))
      plot_map (cp_n{i_t});
      xlabel ("sn"); ylabel ("st")
    endif
  endfor

  ## profile fit
  cp_nn = delta_fit = cn0 = p_c_fit = dcdsnp0 = {};
  dcds_idx = 16;
  switch i_M
    case 1
      idx_r = [2 3];
      dcds_idx = 10;
      sig = 0.5;
    case 2
      idx_r = [2 8];
    case 3
      idx_r = [2 6];
    case 4
      idx_r = [2 4];
      dcds_idx = 12;
  endswitch
  switch i_M
  case 1
    idx_r = [2 12];
    dcds_idx = 12;
    sig = 1;
  case 2
    idx_r = [2 8];
    dcds_idx = 16;
    sig = 1;
  case 3
    idx_r = [2 6];
    dcds_idx = 8;
    sig = 1;
  case 4
    idx_r = [1 5];
    dcds_idx = 8;
    sig = 1;
endswitch
  testplots_fit = 0
  ## profile fit singe thread
##  t1 = tic
##  for i_t = 1
##    printf ([">>> profile fit: " num2str(i_t) " of " num2str(n_t) " records \n"]);
##    [delta_fit{i_t} cp_nn{i_t} cn0{i_t} p_c_fit{i_t} p_sc ~] = erfc_profile_fit (snp, cp_n{i_t}, sf_p, 0, sig, idx_r, dcds_idx, testplots_fit);
##  endfor
##  toc (t1)
##  dt1 = toc (t1);
  nthreads = round (nproc/2); # no use of SMT for this
  printf ([">>> profile fit: " num2str(nthreads) " threads to process " num2str(n_t) " records \n"]);
  t2 = tic
  [delta_fit, cp_nn, cn0, p_c_fit, p_sc, ~] = parcellfun (nthreads, @(par_var) erfc_profile_fit(snp, par_var, sf_p, 0, sig, idx_r, dcds_idx, 0), cp_n, "UniformOutput", false);
  toc (t2)
##  dt2 = toc (t2);
##  dt2 / dt1
  ##
  delta_all = [];
  for i_t = 1:numel(delta_fit)
    for i = 1:numel(delta_fit{1})
      delta_all(i,i_t) = delta_fit{i_t}(i);
    endfor
  endfor
  delta_mean = median (delta_all, 2);
  delta_std = std (delta_all, [], 2);

  ##
  fh = figure (); hold on;
  for i_t = 1:numel(delta_fit)
    clf (fh)
     hold on;
    plot (x, delta_mean, ["k;median;"], "linewidth", 2)
##    plot (x, delta_fit{i_t}, [";i_M = " num2str(i_t) ";"])
    plot (x, outlier_rm(delta_fit{i_t}, movmedian(delta_fit{i_t},81)), [";i_M = " num2str(i_t) ";"])
    pause (0.5)
  endfor
  plot (x, delta_mean, ["k;median;"], "linewidth", 2)
  plot (x, delta_mean+2*delta_std, ["r--;+2 STD;"], "linewidth", 2)
  plot (x, delta_mean-2*delta_std, ["r--;-2 STD;"], "linewidth", 2)
  median (delta_mean)
  median (delta_std)
  100 * 2 * median (delta_std) ./ median (delta_mean)

  plot_map_msh (p_msh, cp_nn{1})
  hold on
  [mi, ixd_s] = min ( abs(cp_nn{1}-0.21),[],2);
  delta_i = ixd_s*sf_p;
  plot3 ((delta_i'), x, ones(1,numel(delta_mean)), "b")
  plot3 ((delta_fit{1}' * 1e3), x, ones(1,numel(delta_mean)), "g")
  plot3 ((delta_mean' * 1e3), x, ones(1,numel(delta_mean)), "r")

  plot (x, (delta_mean' * 1e3), "r")
  plot (x, (movmedian(delta_mean,81)' * 1e3), "k")
  plot (x, imsmooth(movmedian(delta_mean,81),8)' * 1e3, "g")

  valid = [];
  for i_t = 1:n_t
    p_test = polyfit (x(10:end-10), delta_fit{i_t}(10:end-10), 1);
    valid(i_t) = logical(p_test(1)>=0);
  endfor


  delta_all = [];
  k = 0
  for i_t = 1:numel(delta_fit)
    if (valid(i_t))
      k++;
      for i = 1:numel(delta_fit{1})
        delta_all(i,k) = delta_fit{i_t}(i);
      endfor
    endif
  endfor
  delta_mean = median (delta_all, 2);
  delta_std = std (delta_all, [], 2);

  delta_mean = outlier_rm (delta_mean, movmedian(delta_mean,81));
##  delta_mean = imsmooth (delta_mean,32);
  ## TODO: x_stitch result
  ## TODO: load single frame analysis in a_flat_main

  ## output result per section
  mkdir ([save_dir "/" id])
  cd ([save_dir "/" id])
  save -v7 "profiles_dyn_c.v7" cp cp_s cp_b cp_n
  save -v7 "profiles_dyn_msh.v7" profile_depth sf_p snp msh_n x_abs
  save -v7 "profiles_dyn_fit.v7" cn0 cp_nn delta_fit
  ##
  plt_1_h = {mfilename, date, "", ""; "x in mm", "x_abs in m", "delta_c in m", "STD(delta_c) in m"};
  cell2csv ([save_dir_p "delta_dyn_mean.csv"], plt_1_h)
  plt_1_d = [(x+aid.ids_X(i_X))' x_abs' delta_mean*1e3 delta_std*1e3];
  csvwrite ([save_dir_p "delta_dyn_mean.csv"], plt_1_d, "append", "off", "precision", "%.4e")

endif
##
