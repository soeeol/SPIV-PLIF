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
  it_X = 3; # scanned x sections
  ## fixed
  i_L = i_O = 1; # liquid, optical setup
  i_A = 1;
  i_C = 1;
  i_M = 1;
  i_X = 3;
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
  xmin = min (x);
  xmax = max (x);
  y = c_msh{2}(:,1);
  ymin = min (y);
  ymax = max (y);
  z = 0;
  sf_c = get_sf (c_msh);
  sf_u = get_sf (u_msh);
  n_t = numel (c_dat(1,:))
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
  for i = 1:n_t
    pdat.c_dat{i_A, i_C, i_M, i_X}{1} = c_dat{1,i};
    pdat_n{i} = corr_y_offsets (p_dat, it_A, it_C, it_M, it_X);
  endfor

  ## concentration field
  ## linear 2 point ref calib, with intensity reduction towards interface in calib
  cn = cell (1, n_t);
  for i = 1:n_t
    cn{i} = calc_cn ({c_dat{1,i} c_dat{2,1} c_dat{3,1}}, [0 1], c_method, sig=2, testplots=false);
  endfor
  ## time average
  cn_mean = 0;
  for i = 1:n_t # t
    cn_mean = cn_mean + cn{i};
  endfor
  cn_mean = cn_mean / n_t;

  if testplots
    ## fluorescence maps
    plot_map (c_dat{2,1}); title ("Ic0");
    plot_map (c_dat{3,1}); title ("Ic1");
    fh1 = figure ()
    for i = 1:n_t
      clf (fh1);
      surf (c_dat{1,i});
      shading flat;
      view ([0 0 1]);
      colorbar;
      title (["Ic " num2str(i)]);
      pause (0.1);
    endfor
    ## concentration field
    fh2 = figure ()
    for i = 1:n_t
      clf (fh2);
      surf (c_msh{1}, c_msh{2}, c_msh{3}, cn{i});
      shading flat;
      view ([0 0 1]);
      colorbar;
      caxis ([0 1])
      title (["cn " num2str(i)]);
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
  for i = 1:n_t
    [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (c_msh, cn{i}, tol*2, "max", pp, ifg_meas_idx(:,2));
    printf([">>> if detection in " num2str(i) " of " num2str(n_t) " c records \n"]);
    if (i==1)
      fh3 = check_if_plot (xy_if{i}, ispeak{i}, c_msh, cn{i});
      if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
        close (fh3)
      else
        error (msg{2});
      endif
    endif
  endfor

  ## smooth spline fit representation of interface
  for i = 1:n_t
    spf = splinefit (x(ispeak{i}==1), xy_if{i}(ispeak{i}==1,2), round(numel(it_X)*sps*1), "order", 3, "beta", 0.75);
    h_c_g_fit{i} = ppval (spf, x);
  endfor

  ##
  fh = figure ();
  for i = 1:n_t
    clf (fh);
    surf (c_msh{1}, c_msh{2}, cn{i});
    shading flat;
    view ([0 0 1]);
    colorbar;
    hold on;
    plot3 (x, h_c_g_fit{i}, ones(1,numel(x)), "-m");
    title (["cn t#" num2str(i)]);
    xlabel ("x in mm")
    ylabel ("y in mm")
    ylim ([max(h_c_g_fit_mean)-0.5 1.1*max(h_c_g_fit_mean)])
    zlim ([-1 1])
##    pause (0.1);
    print (fh, "-djpeg", "-color", ["-r" num2str(250)], [save_dir_p "/cn_t_" num2str(i) ".jpg"]);
  endfor

  x_idx = round (rand(1) * numel(cn_mean(1,:)))
  fh = figure ();
  for i = 1:n_t
    clf (fh);
    hold on;
    plot (y, cn{i}(:,x_idx), "b")
    plot (y,cn_mean(:,x_idx), "k")
    plot (h_c_g_fit{i}(x_idx)*[1 1], [0 1], "b-")
    plot (h_c_g_fit_mean(x_idx)*[1 1], [0 1], "k-")
    xlim ([max(h_c_g_fit_mean)-0.5 1.1*max(h_c_g_fit_mean)])
    xlabel ("y in mm")
    ylabel ("cn in -")
    title (["cn t#" num2str(i) " (for one random pixel row)"]);
##    pause (0.1);
    print (fh, "-djpeg", "-color", ["-r" num2str(250)], [save_dir_p "/cn_prof_t_" num2str(i) ".jpg"]);
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

  cd (save_dir_p);
  save -v7 "if_stat.v7" n_t x c_h mmax_dev mmad_dev mstd_dev max_dev mad_dev std_dev
endif

##
## TODO: interface normal concentration profile analysis
##
if 0
endif
