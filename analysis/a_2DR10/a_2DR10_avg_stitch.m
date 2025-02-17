##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## stitch section results
## - phi, cn, cp
## - velocity field
## - derived results: delta_u, delta_c_avg, u_s

##
## Author: Sören J. Gerke
##

##
## [00] init
if 1

  ## ap parameters defining the analysis
  ap = [];
  ap.a_type = "a_2DR10_avg_stitch"; # identifier of this analysis

  ## selection of experiments to be analyzed
  ap.ids_A = [60]; # [°] inlination IDs
  ap.ids_C = {"2d-r10"}; # cell IDs
  ap.ids_G = [2]; # [Nl/min] gas flow IDs
  ap.ids_L = {"WG141"}; # liquid IDs
  ap.ids_M = [8 16 32 64]; # [kg/h] mass flow IDs
  ap.ids_O = {"M13"}; # optical setup IDs
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
##  i_M = it_M = 4

  ## parameters concentration transformation
  ap.c_method = "linear"; # method to transform fluorescence intensity to concentration ("linear" / "nonlinear" .. no impact on delta_c_avg)
  ap.c_if_method = "calib"; # method to deal with fluorescence intensity decay at the interface ("calib" / "calib-if" .. high impact on delta_c_avg)

  ## per data variable: method for time series averaging before stitching
  ## ("median" or "mean")
  ## or set one for all
  ap.sd.avg_method = "median"

  ## spline interface fit for stable interface normals (defaults)
  ap.sd.if_sfit_order = 3; # spline order
  ap.sd.if_sfit_sps = 15; # splines devisions per section (.. generally: increase for curved interface, decrease for flat)

  ## prepare directories
  ap.date_str = datestr (now, "yyyymmdd");
  ap.result_dir = [pdir.analyzed ap.a_type "/"];

  ## prepare result storage path
  ap.i_M = it_M;
  ap.i_X = it_X;
  ap.id_meas = get_measid_ap (ap);
  ap.save_dir_id = [ap.result_dir ap.id_meas "/"];
  mkdir (ap.save_dir_id);

  ##
  ## stitching descriptors
  ##
  sd{1}.sd_name = "cn_dyn_avg";
  sd{1}.a_dir = [pdir.analyzed "a_2DR10_dyn_cn_cp/"];
  sd{1}.a_id_sec = ["cn-" ap.c_method "_" ap.c_if_method];
  sd{1}.msh_fn = "cn_dyn.v7";
  sd{1}.msh_var = "c_msh";
  sd{1}.dat_fn{1} = "delta_u_dyn.v7"; ##  "delta_u_dyn.v7" x y_wall delta_u delta_u_fit delta_u_fit_avg
  sd{1}.dat_var{1} = {"delta_u", "y_wall"};
  sd{1}.dat_fn{2} = "cn_dyn.v7"; ##  "cn_dyn.v7" c_msh cn_dyn cn_dyn_avg
  sd{1}.dat_var{2} = {"cn_dyn_avg"};
  sd{1}.dat_fn{3} = "phi_avg.v7"; ##  "phi_avg.v7" c_msh phi_avg phi_des phi_sat x phi_avg_s phi_des_s phi_sat_s
  sd{1}.dat_var{3} = {"phi_avg", "phi_des", "phi_sat"};
  ##
  sd{2}.sd_name = "cp_dyn_avg";
  sd{2}.a_dir = [pdir.analyzed "a_2DR10_dyn_cn_cp/"];
  sd{2}.a_id_sec = ["cp-avg-" ap.c_method "_" ap.c_if_method];
  sd{2}.msh_fn = "avg_cp.v7";
  sd{2}.msh_var = "p_msh_o";
  sd{2}.dat_fn{1} = "avg_cp.v7"; ##  "avg_cp.v7" snp_o p_msh_o cp_nn_avg cp_n_avg
  sd{2}.dat_var{1} = {"cp_n_avg", "cp_nn_avg"};
  sd{2}.dat_fn{2} = "avg_delta_c.v7"; ##  "avg_delta_c.v7" x delta_c a_fit
  sd{2}.dat_var{2} = {"delta_c", "a_fit"};
  sd{2}.dat_fn{3} = "dyn_cp.v7"; ##  "dyn_cp.v7" snp p_msh msh_n cp_n cp_s cp_b
  sd{2}.dat_var{3} = {"cp_s", "cp_b"};
  ##
  sd{3}.sd_name = "phi_avg";
  sd{3}.a_dir = [pdir.processed ""];
  sd{3}.a_id_sec = ["2d_avg_uIc1"];
  sd{3}.msh_fn = "c.v7";
  sd{3}.msh_var = "c_msh";
  sd{3}.dat_fn{1} = "c.v7"; ##  "c.v7" c_msh c_dat c_masks c_h delta_u y_wall
  sd{3}.dat_var{1} = {"c_dat", "y_wall", "delta_u"};
  ##
  sd{4}.sd_name = "u_avg";
  sd{4}.a_dir = [pdir.processed ""];
  sd{4}.a_id_sec = ["2d_avg_uIc1"];
  sd{4}.msh_fn = "u.v7";
  sd{4}.msh_var = "u_msh";
  sd{4}.dat_fn{1} = "u.v7"; ##  "u.v7" u_msh u_dat u_masks u_h
  sd{4}.dat_var{1} = {"u_dat"};
  ##  "y_if.v7" y_if_wall y_if_gas

endif

##
## [10] prepare inter section offset correction
if 1

  ## using gas-liquid interface as inter section offset indicator
  for i_M = it_M
    for i_X = it_X

      ap.i_M = i_M;
      ap.i_X = i_X;

      ## measurement identifier
      ap.id_meas = get_measid_ap (ap);

      ## dyn avg valid
      data_files = glob ([sd{1}.a_dir ap.id_meas "/" "*_" sd{1}.a_id_sec "/" sd{1}.dat_fn{1}]);
      data_file = data_files{end} # use latest result
      x = y_wall = delta_u = [];
      load (data_file, "x", "y_wall", "delta_u");
      x_MX{i_M,i_X} = x;
      y_wall_MX{i_M,i_X} = y_wall;
      delta_u_dyn_avg_MX{i_M,i_X} = calc_vec_avg_cells (delta_u, ap.sd.avg_method);

      ## avg all
      data_files = glob ([sd{3}.a_dir ap.id_meas "/" "*_" sd{3}.a_id_sec "/" sd{3}.dat_fn{1}]);
      data_file = data_files{end} # use latest result
      c_msh = y_wall = delta_u = [];
      load (data_file, "c_msh", "y_wall", "delta_u");
      x_avg_MX{i_M,i_X} = c_msh{1}(1,:);
      y_wall_avg_MX{i_M,i_X} = y_wall;
      delta_u_avg_MX{i_M,i_X} = delta_u;

      ## avg all
      data_files = glob ([sd{4}.a_dir ap.id_meas "/" "*_" sd{4}.a_id_sec "/" sd{4}.dat_fn{1}]);
      data_file = data_files{end} # use latest result
      u_msh = u_dat = [];
      load (data_file, "u_msh", "u_dat");
      u_msh_MX{i_M,i_X} = u_msh;
      ux_MX{i_M,i_X} = u_dat{1};

    endfor
  endfor

  if 0
    ## helper for systematic x offset
    xoff_test = []
    for i_M = it_M
      for i_X = it_X(1:end-1)
        x_l = x_MX{i_M,i_X};
        x_r = x_MX{i_M,i_X+1};
        xpos_l = ap.ids_X(i_X);
        xpos_r = ap.ids_X(i_X+1);
        delta_u_l = delta_u_avg_MX{i_M,i_X};
        delta_u_r = delta_u_avg_MX{i_M,i_X+1};
        xoff_test(i_M,i_X) = calc_xsec_if_offset_x (x_l, x_r, delta_u_l, delta_u_r, xpos_l, xpos_r);
      endfor
    endfor
    xoff_test
  endif

  xoff = zeros (numel(ap.ids_M), numel(ap.ids_X));
  xoff([1 2 4],1) = +0.05;
  xoff(1,3:4) = +0.1; # recorded on another day
##  xoff(2,1) = +0.025;
##  xoff(3,1) = +0.0;
  xoff(4,3:4) = +0.5; # recorded on another day

  yoff = zeros (numel(ap.ids_M), numel(ap.ids_X));
  yoff(1,3) = +0.005;
##  yoff(2,1) = -0.005;
##  yoff(3,3) = +0.005;
  yoff(4,1) = -0.01; yoff(4,3:4) = +0.005;

  if 0
    ## test plot offset correction
    figure ();
    hold on;
    for i_M = it_M
      for i_X = it_X
        plot (x_MX{i_M,i_X} + ap.ids_X(i_X) + xoff(i_M,i_X), delta_u_dyn_avg_MX{i_M,i_X} + yoff(i_M,i_X), "k");
        plot (x_avg_MX{i_M,i_X} + ap.ids_X(i_X) + xoff(i_M,i_X), delta_u_avg_MX{i_M,i_X} + yoff(i_M,i_X), "b");
        plot (x_MX{i_M,i_X} + ap.ids_X(i_X) + xoff(i_M,i_X), y_wall_MX{i_M,i_X} + yoff(i_M,i_X), "k");
        plot (x_avg_MX{i_M,i_X} + ap.ids_X(i_X) + xoff(i_M,i_X), y_wall_avg_MX{i_M,i_X} + yoff(i_M,i_X), "b");
      endfor
    endfor
    legend ({"filtered avg" "avg"})
    legend ("autoupdate", "off");
    plot ((x_sec-0.06)*1e3, [0 0 0 0 0; [1 1 1 1 1]], "--k");
    axis image;
    ylim ([0 max(delta_u_dyn_avg_MX{it_M(end),2})]);
  endif

  ## intra section u vs c offset - fine tuning
  x_u_off = zeros (numel(ap.ids_M), numel(ap.ids_X)); # + shifts u vs c in +x
  x_u_off(1,2) = + 0.02;
  x_u_off(2,2) = + 0.02;
  x_u_off(3,2) = + 0.02;
  x_u_off(4,2) = + 0.02;

  y_u_off = zeros (numel(ap.ids_M), numel(ap.ids_X)); # + shifts u vs c in +y
  y_u_off(1,1) = - 0.005;
  y_u_off(2,1) = - 0.01;
  y_u_off(3,1) = + 0.01;

endif


##
## [20] stich data and test continuity
if 1

  msh_cn = x_cn = y_cn = {};
  cn = delta_u_cn = {};
  phi_avg = phi_des = phi_sat = {};

  msh_cp = x_cp = snp = {};
  delta_c_avg = a_fit_avg = cp_n_avg = cp_nn_avg = {};

  msh_u = x_u = y_u = {};
  dat_u_avg = {};

  msh_c = x_c = y_c = {};
  delta_u_avg = y_wall_avg = {};

  for i_sd = 1 : numel (sd)
    for i_M = it_M

      ## inter section offset settings
      ap.sd.xoff_X = xoff(i_M,:);
      ap.sd.yoff_X = yoff(i_M,:);

      ap.sd.x_u_off_X = x_u_off(i_M,:);
      ap.sd.y_u_off_X = y_u_off(i_M,:);

      ap.i_M = i_M;
      ap.i_X = it_X;

      [msh_gl dat_gl] = stitch_msh_dat (sd{i_sd}, ap);

      switch (sd{i_sd}.sd_name)

        case {"cn_dyn_avg"}
          msh_cn{i_M} = msh_gl;
          x_cn{i_M} = msh_cn{i_M}{1}(1,:);
          y_cn{i_M} = msh_cn{i_M}{2}(:,1);
          ##
          cn_avg{i_M} = dat_gl.cn_dyn_avg;
          delta_u_avg{i_M} = calc_vec_avg_cells (dat_gl.delta_u, ap.sd.avg_method);
          phi_avg{i_M} = dat_gl.phi_avg;
          phi_des_avg{i_M} = dat_gl.phi_des;
          phi_sat_avg{i_M} = dat_gl.phi_sat;

        case {"cp_dyn_avg"}
          msh_cp{i_M} = msh_gl;
          x_cp{i_M} = msh_cp{i_M}{1}(1,:);
          snp{i_M} = msh_cp{i_M}{2}(:,1);
          ##
          cp_s_avg{i_M} = calc_vec_avg_cells (dat_gl.cp_s, ap.sd.avg_method);
          cp_b_avg{i_M} = calc_vec_avg_cells (dat_gl.cp_b, ap.sd.avg_method);
          ##
          delta_c_avg{i_M} = dat_gl.delta_c;
          a_fit_avg{i_M} = dat_gl.a_fit;
          cp_n_avg{i_M} = dat_gl.cp_n_avg;
          cp_nn_avg{i_M} = dat_gl.cp_nn_avg;

        case {"phi_avg"}
          msh_c{i_M} = msh_gl;
          x_c{i_M} = msh_c{i_M}{1}(1,:);
          y_c{i_M} = msh_c{i_M}{2}(:,1);
          ##
          y_wall_avg{i_M} = dat_gl.y_wall;
          delta_u_phi_avg{i_M} = dat_gl.delta_u;

        case {"u_avg"}
          msh_u{i_M} = msh_gl;
          x_u{i_M} = msh_u{i_M}{1}(1,:);
          y_u{i_M} = msh_u{i_M}{2}(:,1);
          ##
          dat_u_avg{i_M} = dat_gl.u_dat;

      endswitch

    endfor
  endfor

  ## check continuity of delta_u, y_wall and surface velocity estimate of stitching result
  fh = figure ();
  hold on;
  for i_M = it_M
    u_M_ip = interp2 (msh_u{i_M}{1}, msh_u{i_M}{2}, dat_u_avg{i_M}{4}, msh_cn{i_M}{1}, msh_cn{i_M}{2}, "pchip", 0.0); #
    u_s_ip = get_surface_val (msh_cn{i_M}, u_M_ip, delta_u_avg{i_M}, "min_y_dist");
    plot (x_cn{i_M}, u_s_ip/max(u_s_ip), "b;u_s norm;");
    plot (x_c{i_M}, delta_u_avg{i_M}/max(delta_u_avg{i_M}), "r;delta_u avg;");
    legend ("autoupdate", "off");
    plot ((x_sec-0.06)*1e3, [0 0 0 0 0; [1 1 1 1 1]], "--k");
  endfor
  xlabel ("x* in mm");
  ylabel ("");

  ## check cn field for continuity around section borders
  fh = figure ();
  plot_map_msh (msh_cn{i_M}, cn_avg{i_M}, fh)
  hold on;
  draw_cell (ap.ids_C{ap.i_C}, 0, 1);
  plot (x_c{i_M}, y_wall_avg{i_M}, "r");
  plot (x_c{i_M}, delta_u_avg{i_M}, "r");
  plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 2.5 * [1 1 1 1 1]], "--m");
  xlabel ("x* in mm");
  ylabel ("y in mm");
  ylim ([-0.1 2.5]);

  ## velocity profile per section - comparison with Nusselt profile
  ## check y u offset - valid for flat film
  ##

  yp_eq_nd = linspace (0, 1, 101);
  up_eq_nd = model_filmflow_laminar_u_profile (yp_eq_nd, 1, 1);
  for i_M = it_M

    fh = figure ();
    hold on;
    plot (up_eq_nd, yp_eq_nd, "-k;Nusselt normalized;", "linewidth", 2);

    yp = y_cn{i_M};

    for i_X = it_X

      ## select x section range for most flat film
      x_u_l = ap.ids_X(i_X) - 2;
      x_u_u = ap.ids_X(i_X) + 2;
      if ap.ids_X(i_X) < 0
        x_u_l = ap.ids_X(i_X) - 4;
        x_u_u = ap.ids_X(i_X) + 0;
      elseif ap.ids_X(i_X) > 0
        x_u_l = ap.ids_X(i_X) - 0;
        x_u_u = ap.ids_X(i_X) + 4;
      endif
      idx_sec = (x_u{i_M} >= x_u_l) & (x_u{i_M} <= x_u_u);

      up_mean = median (dat_u_avg{i_M}{1}(:,idx_sec), 2);
      up = interp1 (y_u{i_M}, up_mean, yp);
      ysp = median (delta_u_avg{i_M}(idx_sec));
      [~, idx_us] = min (abs (yp - ysp));
      usp = up (idx_us);

      plot (up / usp, yp / ysp, [".-;exp M" num2str(i_M) "X" num2str(i_X) " delta_u = " num2str(ysp) " mm;"]);

    endfor

    legend ("autoupdate", "off");
    plot ([0 1], 1 * [1 1], "--k");
    xlim ([0 1.1]);
    ylim ([0 1.1]);
    xlabel ("u / u_s");
    ylabel ("y / delta_u");
    legend ("location", "northwest");
    print (fh, "-djpeg", "-r500", [ap.save_dir_id "test_y_u_off___i_M=" num2str(i_M) ".jpg"]);
    close (fh);

  endfor


  ##
  ## check alignment of velocity vectors with interfaces and for continuity around section borders
  ##
  dx = 0.08; # mm double of measurement IA size
  dy = 0.08; # mm
  lim_um = 100e-6;
  for i_M = it_M
    mask_g{i_M} = masking ("gas", size (msh_u{i_M}{1}), min (y_u{i_M}), delta_u_avg{i_M}, get_sf (msh_u{i_M}), 0, nan);
    mask_w{i_M} = masking ("wall", size (msh_u{i_M}{1}), min (y_u{i_M}), y_wall_avg{i_M}, get_sf (msh_u{i_M}), 2, nan);
    mask = ones ( size (msh_u{i_M}{1}));

    lim_x = [min(x_u{i_M}) max(x_u{i_M})]; # in mm
    lim_y = [min(y_u{i_M}) max(delta_u_avg{i_M}) + 2*dy]; # in mm

    [x_v y_v ux_v uy_v um_v] = u_xy_vec (msh_u{i_M}, dat_u_avg{i_M}{1}, dat_u_avg{i_M}{2}, mask, dx, dy, lim_x, lim_y, lim_um);

    fh = figure ();
    plot_map_msh (msh_u{i_M}, dat_u_avg{i_M}{4}, fh)
    caxis ([ 0 median(max(dat_u_avg{i_M}{4}, [], 2), "omitnan") ]);
    caxis ([ 0 1.5*median(max(dat_u_avg{i_M}{4}, [], 2), "omitnan") ]);
    hold on;
    quiver (x_v, y_v, ux_v, uy_v, 1, "k");
    axis image;
    plot (x_cn{i_M}, delta_u_avg{i_M}, "r-", "linewidth", 1);
    plot (x_cn{i_M}, y_wall_avg{i_M}, "r-");
    xlabel ("x in mm");
    ylabel ("y in mm");
  endfor


  ## dirty fix for small section of velocity field of M 32: severe particle accumulation blocked PIV analysis
  ## replace area with moving median:
  ## -5.4 < x < -4.3
  ## y < 0.16
  for i_M = it_M
    if i_M == 3
      for i = 1 : numel (dat_u_avg{i_M})
        idx_x = (msh_u{i_M}{1} >= -5.3) & (msh_u{i_M}{1} <= -4.4);
        idx_y = (msh_u{i_M}{2} >= +0.02) & (msh_u{i_M}{2} <= +0.18);
        idx = idx_x & idx_y;
        switch (i)
          case {1, 4}
            median_wall_profile = repmat (max (dat_u_avg{i_M}{i} .* idx_y, [], 2), 1, size (idx_x, 2));
          otherwise
            median_wall_profile = repmat (median (dat_u_avg{i_M}{i} .* idx_y, 2, "omitnan"), 1, size (idx_x, 2));
        endswitch
        dat_u_avg{i_M}{i} = dat_u_avg{i_M}{i} .* (! idx) + median_wall_profile .* idx;
      endfor
    endif
  endfor

endif

##
## [30] interpolate on common mesh and store result
if 1

  msh = msh_cn;
  x = y = {};
  msh_p = {}; # to keep the finer profile resolution
  s_n = s_t = {};

  y_wall = delta_u = {};
  delta_u_fit = mask_g = mask_w = {};
  phi = phi_des = phi_sat = {};
  cn = {};
  cp_b = cp_s = cp_n = cp_nn = delta_c = a_fit_cp_scale = {};
  u_dat = u_x = u_y = u_z = u_m = {};
  u_s = incl_s = curve_s = l_s = y_h = vfr = {};

  ## stitched data
  for i_M = it_M

    x{i_M} = msh{i_M}{1}(1,:);
    y{i_M} = msh{i_M}{2}(:,1);
    z = msh{i_M}{3}(1,1);
    sf = get_sf (msh{i_M});
    sf_p = get_sf (msh_cp{i_M});
    s_n_min = min (snp{i_M})
    s_n_max = max (snp{i_M})
    s_n{i_M} = linspace (0, s_n_max, (s_n_max-s_n_min)/sf_p(2)+1);
    [XX, YY, ZZ] = meshgrid (x{i_M}, s_n{i_M}, z);
    msh_p{i_M} = {XX, YY, ZZ};

    ## wall
    y_wall{i_M} = interp1 (x_c{i_M}, y_wall_avg{i_M}, x{i_M}, "pchip", "extrap");
    ## gas-liquid interface
    delta_u{i_M} = interp1 (x_cn{i_M}, rm_ext (x_cn{i_M}, delta_u_avg{i_M}, 101), x{i_M}, "pchip", "extrap");

    ## using smooth gas-liquid interface to mask_g the reflective part but leaving some little extra tolerance
    spf{i_M} = splinefit (double (x{i_M}),  double (rm_ext (x{i_M}, delta_u{i_M}, 101)), ap.sd.if_sfit_sps * numel(it_X), "order", ap.sd.if_sfit_order, "beta", 0.75);
    delta_u_fit{i_M} = ppval (spf{i_M}, x{i_M});
    mask_g{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, get_sf(msh{i_M}), 0, nan);
    mask_w{i_M} = masking ("wall", size (msh{i_M}{1}), min (y{i_M}), y_wall{i_M}, get_sf(msh{i_M}), 0, nan);

    ## avg fluorescence recordings
    phi{i_M} = interp2 (msh_cn{i_M}{1}, msh_cn{i_M}{2}, phi_avg{i_M}, msh{i_M}{1}, msh{i_M}{2}, "pchip", 0.0);
    phi_des{i_M} = interp2 (msh_cn{i_M}{1}, msh_cn{i_M}{2}, phi_des_avg{i_M}, msh{i_M}{1}, msh{i_M}{2}, "pchip", 0.0);
    phi_sat{i_M} = interp2 (msh_cn{i_M}{1}, msh_cn{i_M}{2}, phi_sat_avg{i_M}, msh{i_M}{1}, msh{i_M}{2}, "pchip", 0.0);

    ## normalized concentration field
    cn{i_M} = interp2 (msh_cn{i_M}{1}, msh_cn{i_M}{2}, cn_avg{i_M}, msh{i_M}{1}, msh{i_M}{2}, "pchip", 0.0);

    ## bulk and surface normalized concentration
    cp_b{i_M} = interp1 (x_cp{i_M}, rm_ext (x{i_M}, cp_b_avg{i_M}, 101), x{i_M}, "pchip", "extrap");
    cp_s{i_M} = interp1 (x_cp{i_M}, rm_ext (x{i_M}, cp_s_avg{i_M}, 101), x{i_M}, "pchip", "extrap");
    ## concentration profiles
    cp_n{i_M} = interp2 (msh_cp{i_M}{1}, msh_cp{i_M}{2}, cp_n_avg{i_M}, msh_p{i_M}{1}, msh_p{i_M}{2}, "pchip", 0.0);
    cp_nn{i_M} = interp2 (msh_cp{i_M}{1}, msh_cp{i_M}{2}, cp_nn_avg{i_M}, msh_p{i_M}{1}, msh_p{i_M}{2}, "pchip", 0.0);

    ## concentration boundary layer thickness
    delta_c{i_M} = interp1 (x_cp{i_M}, rm_ext (x_cp{i_M}, delta_c_avg{i_M}, 101), x{i_M}, "pchip", "extrap");
    a_fit_cp_scale{i_M} = interp1 (x_cp{i_M}, a_fit_avg{i_M}(:,2), x{i_M}, "pchip", "extrap");

    ## velocity field
    for i_u = 1 : numel (dat_u_avg{i_M})
      u_dat{i_M}{i_u} = interp2 (msh_u{i_M}{1}, msh_u{i_M}{2}, dat_u_avg{i_M}{i_u}, msh{i_M}{1}, msh{i_M}{2}, "pchip", 0.0);
    endfor
    u_x{i_M} = u_dat{i_M}{1};
    u_y{i_M} = u_dat{i_M}{2};
    u_z{i_M} = u_dat{i_M}{3};
    u_m{i_M} = u_dat{i_M}{4};

    ## extract:
    ## - surface velocity
    u_s{i_M} = get_surface_val (msh{i_M}, u_dat{i_M}{4}, delta_u_fit{i_M}, "min_y_dist");
    u_s{i_M} = rm_ext (x{i_M}, u_s{i_M}, 401);
    ## - gas-liquid interface inclination vs. horizontal
    [~, ~, incl_s{i_M}] = calc_if_len (double(x{i_M}), delta_u_fit{i_M});
    ## - gas-liquid interface curvature estimate
    [curve_s{i_M}, ~] = calc_curvature_xy (double (x{i_M}), double (delta_u_fit{i_M}));
    ## - gas-liquid interface length and tangent coordinate
    [~, l_s{i_M}, ~, s_t{i_M}] = calc_if_len (double(x{i_M}), delta_u_fit{i_M});
    ## - local film thickness
    y_h{i_M} = delta_u_fit{i_M} - y_wall{i_M};

    ## flow profile specific flow rate along x over cell width
    vfr{i_M} = zeros (size (x{i_M}));
    u_x_masked = u_x{i_M} .* mask_w{i_M} .* mask_g{i_M};
    u_x_masked(isnan(u_x_masked)) = 0.0;
    for i_p = 1 : numel (x{i_M})
      u_x_yp = u_x_masked(:,i_p);
      vfr{i_M}(i_p) = sum (u_x_yp) * sf(2) / 1000 * cell_width / 1000 * 3600; # m^3 / s
    endfor

  endfor

  cd (ap.save_dir_id);
  save -v7 "sd_xy_map.v7" msh x y phi phi_des phi_sat cn u_x u_y u_z u_m
  save -v7 "sd_nt_map.v7" msh_p s_n s_t cp_n cp_nn
  save -v7 "sd_x_vec.v7" x cp_b cp_s l_s curve_s incl_s y_wall y_h delta_u delta_u_fit u_s delta_c a_fit_cp_scale
  save -text "sd_ap.txt" ap
  cd (pdir.work);

  ## stitching result overview plots
  if 1
    for i_M = it_M

      ## x vectors

      fh = figure ();
      hold on;
      plot (x{i_M}, cp_b{i_M}, ";cn bulk;");
      plot (x{i_M}, cp_s{i_M}, ";cn surface;");
      legend ("autoupdate", "off");
      plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 1.0 * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      ylabel ("cn in -");
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_cn_b___cn_s___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      fh = figure ();
      hold on;
      plot (x{i_M}, y_wall{i_M}, "k;wall;");
      plot (x{i_M}, delta_u{i_M}, "b;delta_u;");
      plot (x{i_M}, delta_u_fit{i_M}, "r;fit;");
      plot (x{i_M}, y_h{i_M}, "m;h;");
      legend ("autoupdate", "off");
      plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 2.0 * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      ylabel ("y in mm");
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_delta_u___y_h___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      fh = figure ();
      hold on;
      plot (x{i_M}, y_wall{i_M}, "k;y wall in mm;");
      plot (x{i_M}, delta_u_fit{i_M}, "k;delta_u fit in mm;");
      plot (x{i_M}, (incl_s{i_M}), ";inclination in rad;");
      plot (x{i_M}, curve_s{i_M}, ";curvature in 1/m;");
      legend ("autoupdate", "off");
      plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 2.0 * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_incl_s___curv_s___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      fh = figure ();
      hold on;
      plot (x{i_M}, y_wall{i_M}, "k;wall;");
      plot (x{i_M}, delta_u_fit{i_M}, "k;delta_u fit in mm;");
      plot (x{i_M}, u_s{i_M} / median (u_s{i_M}), ";u_s norm;");
      legend ("autoupdate", "off");
      plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 2.0 * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_u_s_norm___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      fh = figure ();
      hold on;
      plot (x{i_M}, vfr{i_M}, "b;VFR in m^3/s;");
      legend ("autoupdate", "off");
      legend ("location", "southeast");
      plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(vfr{i_M}) * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      ylabel ("VFR in m^3 / s");
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_vfr___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      fh = figure ();
      hold on;
      plot (x{i_M}, s_t{i_M} ./ (x{i_M}-min(x{i_M})), ";s_t;");
      legend ("autoupdate", "off");
      plot ((x_sec-0.06)*1e3, [min(s_t{i_M} ./ (x{i_M}-min(x{i_M}))) * [1 1 1 1 1] ; max(s_t{i_M} ./ (x{i_M}-min(x{i_M}))) * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_s_t_by_x___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      fh = figure ();
      hold on;
      plot (x{i_M}, delta_c{i_M}, ";delta_c;");
      legend ("autoupdate", "off");
      plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 0.1 * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      ylabel ("delta_c in mm");
      ylim ([0 0.1]);
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_delta_c___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      fh = figure ();
      hold on;
      plot (x{i_M}, a_fit_cp_scale{i_M}, ";a_fit_cp_scale;");
      legend ("autoupdate", "off");
      plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 2 * [1 1 1 1 1]], "--k");
      xlabel ("x* in mm");
      ylabel ("a_fit_cp_scale");
      print (fh, "-djpeg", "-r500", [ap.save_dir_id "vec_a_fit_cp_scale___i_M=" num2str(i_M) ".jpg"]);
      close (fh);

      ## maps
      if 1

        fh = plot_map_msh (msh{i_M}, phi{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_phi___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh{i_M}, phi_des{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_phi_des___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh{i_M}, phi_sat{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_phi_sat___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh{i_M}, cn{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_cn___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh_p{i_M}, cp_n{i_M}, []);
        caxis ([0 1]);
        set (gca (), "ydir", "reverse");
        ylim ([0 0.1]);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 0.1 * [1 1 1 1 1]], "--k");
        xlabel ("x* in mm");
        ylabel ("s_n in mm");
        print (fh, "-djpeg", "-r500", [ap.save_dir_id "map_cp_n___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh_p{i_M}, cp_nn{i_M}, []);
        caxis ([0 1]);
        set (gca (), "ydir", "reverse")
        ylim ([0 0.1]);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; 0.1 * [1 1 1 1 1]], "--k");
        xlabel ("x* in mm");
        ylabel ("s_n in mm");
        print (fh, "-djpeg", "-r500", [ap.save_dir_id "map_cp_nn___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh{i_M}, u_x{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        plot (x{i_M}, delta_u_fit{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_u_x___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh{i_M}, u_y{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        plot (x{i_M}, delta_u_fit{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_u_y___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh{i_M}, u_z{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        plot (x{i_M}, delta_u_fit{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_u_z___i_M=" num2str(i_M) ".jpg"]);
        close (fh);

        fh = plot_map_msh (msh{i_M}, u_m{i_M}, []);
        hold on;
        plot ((x_sec-0.06)*1e3, [0 0 0 0 0; max(y{i_M}) * [1 1 1 1 1]], "--k");
        plot (x{i_M}, y_wall{i_M}, "r");
        plot (x{i_M}, delta_u_fit{i_M}, "r");
        axis image;
        xlabel ("x* in mm");
        ylabel ("y in mm");
        print (fh, "-djpeg", "-r1000", [ap.save_dir_id "map_u_m___i_M=" num2str(i_M) ".jpg"]);
        close (fh);
      endif

    endfor
  endif

endif

