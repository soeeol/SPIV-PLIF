##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## sd section results
## - cn
## TODO: - velocity field

##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ## ap parameters defining the analysis
  ap = [];
  ap.a_type = "a_2DR10_avg_stitch"; # identifier of this analysis
  ap.c_method = "linear"; # method to transform fluorescence intensity to concentration ("linear" / "nonlinear" .. no impact on delta_c)
  ap.c_if_method = "calib"; # method to deal with fluorescence intensity decay at the interface ("calib" / "calib-if" .. high impact on delta_c)

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
  i_M = it_M = 3
##  i_M = it_M = 1:4

  ## prepare directories
  ap.date_str = datestr (now, "yyyymmdd");
  ap.result_dir = [pdir.analyzed ap.a_type "/"];


  ## per data variable: method for time series averaging before stitching
  ## ("median" or "mean")
  ## or set one for all
  ap.sd.avg_method = {"median"}

endif

## [10] x-sd dyn avg data and store
## - first sd cn map together with interface to estimate best intersection offset

if 1

  ## stitching descriptors

  sd{1}.sd_name = "cn_dyn_avg";
  sd{1}.a_dir = [pdir.analyzed "a_2DR10_dyn_cn_cp/"];
  sd{1}.a_id_sec = ["cn-" ap.c_method "_" ap.c_if_method];
  sd{1}.msh_fn = "cn_dyn.v7";
  sd{1}.msh_var = "c_msh";
  sd{1}.dat_fn{1} = "delta_u_dyn.v7";
  sd{1}.dat_var{1} = {"delta_u_fit_avg" "y_wall"};
  sd{1}.dat_fn{2} = "cn_dyn.v7";
  sd{1}.dat_var{2} = {"cn_dyn_avg"};
  sd{1}.dat_fn{3} = "phi_avg.v7";
  sd{1}.dat_var{3} = {"phi_avg" "phi_des" "phi_sat"};

  sd{2}.sd_name = "cp_dyn_avg";
  sd{2}.a_dir = [pdir.analyzed "a_2DR10_dyn_cn_cp/"];
  sd{2}.a_id_sec = ["cp-avg-" ap.c_method "_" ap.c_if_method];
  sd{2}.msh_fn = "avg_cp.v7";
  sd{2}.msh_var = "p_msh_o";
  sd{2}.dat_fn{1} = "avg_cp.v7";
  sd{2}.dat_var{1} = {"cp_n_avg"};
  sd{2}.dat_fn{2} = "avg_delta_c.v7";
  sd{2}.dat_var{2} = {"delta_c"};

  sd{3}.sd_name = "u_avg";
  sd{3}.a_dir = [pdir.processed ""];
  sd{3}.a_id_sec = ["2d_avg_uIc1"];
  sd{3}.msh_fn = "u.v7";
  sd{3}.msh_var = "u_msh";
  sd{3}.dat_fn{1} = "u.v7";
  sd{3}.dat_var{1} = {"u_dat"};

  sd{4}.sd_name = "c_avg";
  sd{4}.a_dir = [pdir.processed ""];
  sd{4}.a_id_sec = ["2d_avg_uIc1"];
  sd{4}.msh_fn = "c.v7";
  sd{4}.msh_var = "c_msh";
  sd{4}.dat_fn{1} = "c.v7";
  sd{4}.dat_var{1} = {"c_dat", "y_wall", "delta_u"};

  sd{5}.sd_name = "y_if";
  sd{5}.a_dir = [pdir.processed ""];
  sd{5}.a_id_sec = ["_2d_avg_uIc1"];
  sd{5}.msh_fn = "c.v7";
  sd{5}.msh_var = "c_msh";
  sd{5}.dat_fn{1} = "y_if.v7";
  sd{5}.dat_var{1} = {"y_if_gas" "y_if_wall"};

  ## TODO: per data variable: if is time series calc temporal average before stitching

  ##
  ## prepare inter section offset correction
  ##

  ## using gas-liquid interface for înter section offset indicator
  for i_M = it_M
    for i_X = it_X

      ap.i_M = i_M;
      ap.i_X = i_X;

      ## measurement identifier
      ap.id_meas = get_measid_ap (ap);

      ## per section
      data_files = glob ([sd{1}.a_dir ap.id_meas "/" "*_" sd{1}.a_id_sec "/" sd{1}.dat_fn{1}]);
      data_file = data_files{end} # use latest result

      load (data_file, "x", "y_wall", "delta_u_fit_avg");
      x_MX{i_M,i_X} = x;
      y_wall_MX{i_M,i_X} = y_wall;
      delta_u_avg_MX{i_M,i_X} = delta_u_fit_avg;

    endfor
  endfor


  xoff = zeros (numel(ap.ids_M), numel(ap.ids_X));
  xoff(:,1) = +0.005;
  xoff(:,3:4) = +0.1;
  xoff(4,3) = +0.3; xoff(4,4) = +0.3; # recorded on another day

  yoff = zeros (numel(ap.ids_M), numel(ap.ids_X));
  yoff(2,2) = +0.005;
  yoff(4,1) = -0.005;
  yoff(4,3:4) = +0.005;

  ## before inter section offset correction
  figure (); hold on;
  for i_M = it_M
    for i_X = it_X
      plot (x_MX{i_M,i_X} + ap.ids_X(i_X) + xoff(i_M,i_X), delta_u_avg_MX{i_M,i_X} + yoff(i_M,i_X), "k");
      plot (x_MX{i_M,i_X} + ap.ids_X(i_X) + xoff(i_M,i_X), y_wall_MX{i_M,i_X}, "k");
    endfor
  endfor
  axis image;
  ylim ([0 max(delta_u_avg_MX{it_M(end),2})]);

##  ## helper to find systematic offset
##  xoff_test = []
##  for i_M = it_M
##    for i_X = it_X(1:end-1)
##      x_l = x_MX{i_M,i_X};
##      x_r = x_MX{i_M,i_X+1};
##      xpos_l = ap.ids_X(i_X);
##      xpos_r = ap.ids_X(i_X+1);
##      delta_u_l = delta_u_avg_MX{i_M,i_X};
##      delta_u_r = delta_u_avg_MX{i_M,i_X+1};
##      xoff(i_M,i_X) = calc_xsec_if_offset_x (x_l, x_r, delta_u_l, delta_u_r, xpos_l, xpos_r);
##    endfor
##  endfor
##  xoff_test


###### c records to find connection of x positionings
#### i_M = 1
##20211125-M1-15 xoff25_1
##20211128-M2-28 xoff28_2 = 0
##20211128-M2-33 xoff28_3
##20211128-M2-34 xoff28_4

#### i_M = 2
##20211125-M1-14 xoff25_1 = 0.05
##20211128-M2-27 xoff28_2
##20211128-M2-32 xoff28_3 = 0.1;
##20211128-M2-35 xoff28_4 = 0.1;

#### i_M = 3
##20211125-M1-10 xoff25_1
##20211128-M2-25 xoff28_2
##20211129-M1-8
##20211129-M1-9

#### i_M = 4
##20211125-M1-12 xoff25_1
##20211128-M2-26 xoff28_2
##20211130-M1-6
##20211130-M1-4

## TODO: check that u records match ... no they wont exactly since they were usually recorded with phi_sat,
## if no huuge x offset its best not to offset x since resolution of PIV is much smaller than c


## x = 0 position is fixed by wall measured


##  yoff = zeros (numel(ap.ids_M), numel(ap.ids_X));
##  yoff(1,4) = +0.005;
##  yoff(2,1) = +0.02;
##  yoff(2,4) = -0.005;
##  yoff(3,3) = +0.015;
##  yoff(3,4) = +0.015;
##  yoff(4,1) = -0.015;
##  yoff(4,3) = +0.01;
##  yoff(4,4) = +0.005;

  msh_cn = x = y = {};
  cn = delta_u = y_wall = {};
  phi_avg = phi_des = phi_sat = {};

  msh_cp = snp = {};
  delta_c = cp_n = {};

  msh_u = x_u = y_u = {};
  dat_u = {};

  for i_sd = [1,2,3,4]#:2# : numel (sd)
    for i_M = 3#it_M

      ## inter section offset settings
      ap.sd.xoff_X = xoff(i_M,:);
      ap.sd.yoff_X = yoff(i_M,:);

      ap.i_M = i_M;
      ap.i_X = it_X;

      ## prepare result storage path
      ap.id_meas = get_measid_ap (ap);
      ap.save_dir_id = [ap.result_dir ap.id_meas "/" ap.date_str "_" "avg_stitch" "/"];
      mkdir (ap.save_dir_id);

      [msh_gl dat_gl] = stitch_msh_dat (sd{i_sd}, ap);


      switch (sd{i_sd}.sd_name)
        case {"cn_dyn_avg"}
          msh_cn{i_M} = msh_gl;
          x{i_M} = msh_gl{1}(1,:);
          y{i_M} = msh_gl{2}(:,1);
          cn{i_M} = dat_gl.cn_dyn_avg;
          delta_u{i_M} = dat_gl.delta_u_fit_avg;
          y_wall{i_M} = dat_gl.y_wall;
          phi_avg{i_M} = dat_gl.phi_avg;
          phi_des{i_M} = dat_gl.phi_des;
          phi_sat{i_M} = dat_gl.phi_sat;

        case {"cp_dyn_avg"}
          msh_cp{i_M} = msh_gl;
          snp{i_M} = msh_gl{2}(:,1);
          delta_c{i_M} = dat_gl.delta_c;
          cp_n{i_M} = dat_gl.cp_n_avg;
        case {"u_avg"}
          msh_u{i_M} = msh_gl;
          x_u{i_M} = msh_gl{1}(1,:);
          y_u{i_M} = msh_gl{2}(:,1);
          dat_u{i_M} = dat_gl.u_dat;
      endswitch

    endfor
  endfor

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (x{i_M}, delta_c{i_M}, [";i_M = " num2str(i_M) ";"]);
  endfor
  xlabel ("x* in mm");
  ylabel ("delta_c in mm");

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (x{i_M}, delta_u{i_M}, [";i_M = " num2str(i_M) ";"]);
    plot (x{i_M}, y_wall{i_M}, ["k;i_M = " num2str(i_M) ";"]);
  endfor
  axis image;
  xlabel ("x* in mm");
  ylabel ("delta_u in mm");


  fh = figure ();
  hold on;
  plot (x{i_M}, delta_c{i_M}/median(delta_c{i_M}), "r");
  plot (x{i_M}, delta_u{i_M}/median(delta_u{i_M}), "b");
  ## u_s
  xlabel ("x* in mm");
  ylabel ("y in mm");

  fh = figure ();
  plot_map_msh (msh_cn{i_M}, phi_sat{i_M}, fh)
  hold on;
  draw_cell (ap.ids_C{ap.i_C}, 0, 1);
  plot (x{i_M}, y_wall{i_M}, "r");
  plot (x{i_M}, delta_u{i_M}, "r");
  xlabel ("x* in mm");
  ylabel ("y in mm");
  axis image;
  ylim ([-0.1 2.5])


  fh = figure ();
  plot_map_msh (msh_gl, dat_gl.c_dat{3}, fh)
  hold on;
##  draw_cell (ap.ids_C{ap.i_C}, 0, 1);
  plot (msh_gl{1}(1,:), dat_gl.y_wall, "r");
  plot (msh_gl{1}(1,:), dat_gl.delta_u, "r");
  xlabel ("x* in mm");
  ylabel ("y in mm");
  axis image;
  ylim ([-0.1 2.5])


  fh = figure ();
  plot_map_msh (msh_cp{i_M}, cp_n{i_M}, fh)
  xlabel ("x* in mm");
  ylabel ("snp in mm");



      fh = figure ();
##      plot_map_msh (msh_cn{i_M}, phi_avg{i_M}, fh)
##      plot_map_msh (msh_cn{i_M}, phi_des{i_M}, fh)
      plot_map_msh (msh_cn{i_M}, phi_sat{i_M}, fh)
      hold on;
      draw_cell (ap.ids_C{ap.i_C}, 0, 1);
      plot (x{i_M}, y_wall{i_M}, "r");
      plot (x{i_M}, delta_u{i_M}, "r");
      xlabel ("x* in mm");
      ylabel ("y in mm");
      axis image;
      ylim ([-0.1 2.5]);


      fh = figure ();
      plot_map_msh (msh_cn{i_M}, cn{i_M}, fh)
      hold on;
      draw_cell (ap.ids_C{ap.i_C}, 0, 1);
      plot (x{i_M}, y_wall{i_M}, "r");
      plot (x{i_M}, delta_u{i_M}, "r");
      xlabel ("x* in mm");
      ylabel ("y in mm");
      axis image;
      ylim ([-0.1 2.5]);

      fh = figure ();
      plot_map_msh (msh_u{i_M}, dat_u{i_M}{4}, fh)
      caxis ([ 0 median(max(dat_u{i_M}{4}, [], 2), "omitnan") ])
      hold on;
      plot (msh_gl{1}(1,:), dat_gl.y_wall, "r");
      plot (msh_gl{1}(1,:), dat_gl.delta_u, "r");
      xlabel ("x* in mm");
      ylabel ("y in mm");
      axis image;
      ylim ([-0.1 2.5]);

      ##    if (i_sd==1)
      fh = figure ();
      hold on;
      for i = 1:3
        plot (msh_gl{1}(1,:), dat_gl.y_if_gas{i}, [";gas phi idx = " num2str(i) ";" ]);
        plot (msh_gl{1}(1,:), dat_gl.y_if_wall{i}, "k;wall;");
      endfor
      axis image;
      xlabel ("x* in mm");
      ylabel ("y in mm");


##      print (fh, "-dpng", "-color", ["-r" num2str(500)], [ap.save_dir_id sd{i_sd}.sd_name "_stitched.png"]);
##    endif

##    if (i_sd==2)
##
##      fh = plot_map_msh (msh_gl, dat_gl.cp_n_avg, [])
##      hold on;
##      plot (msh_gl{1}(1,:), dat_gl.delta_c, "r");
##      xlabel ("x* in mm");
##      ylabel ("s_n in mm");
##      print (fh, "-dpng", "-color", ["-r" num2str(500)], [ap.save_dir_id sd{i_sd}.sd_name "_stitched.png"]);
##
##    endif
##
##
##    if (i_sd==3)
##      yoff = 0.08 # offset u vs c msh
##
##      fh = plot_map_msh ({msh_gl{1} msh_gl{2}+yoff msh_gl{3}}, dat_gl.u_dat{4}, [])
##      hold on;
####      plot (msh_gl{1}(1,:), dat_gl.delta_c, "r");
##      axis image
####      caxis ([0 0.5])
##      xlabel ("x in mm");
##      ylabel ("y in mm");
##      for i_M = it_M
##        for i_X = it_X
##          plot (x_MX{i_M,i_X}+ap.ids_X(i_X), 0. + delta_u_avg_MX{i_M,i_X}, "r");
####          plot (x_MX{i_M,i_X}+ap.ids_X(i_X), 0.02 + delta_u_avg_MX{i_M,i_X}, "m");
####          plot (x_MX{i_M,i_X}+ap.ids_X(i_X), 0.04 + delta_u_avg_MX{i_M,i_X}, "c");
##          plot (x_MX{i_M,i_X}+ap.ids_X(i_X), y_wall_MX{i_M,i_X}, "r");
##        endfor
##      endfor
##
##      figure ();
##      hold on;
##      i_p = [20, 50, 600, 750, 800]
##      for k = 1:numel(i_p)
##
##        i = i_p(k)
##
##        u_p = dat_gl.u_dat{4}(:,i);
##
##        h_lim_l = 0.2;
##        h_lim_h = 1.04;
##
##        y_p = msh_gl{2}(:,1) + yoff;
##        idx_p = (y_p >= h_lim_l) & (y_p <= h_lim_h);
##
##        p_uy_fit = polyfit (y_p(idx_p), u_p(idx_p), 2);
##
##        uyoff = min (roots (p_uy_fit))
##
##        plot (msh_gl{2}(:,1) + yoff, u_p, [";" num2str(i) ";"])
##        plot (msh_gl{2}(:,1) + yoff, polyval (p_uy_fit, y_p), ["--;off " num2str(i) ";"])
##
##      endfor
##
##    endif
##  endfor

endif

## procedure:
## 1. stitch cn, y_wall and delta_u {best case: dyn result is available}
##    - defines global mesh with common resolution
## 2. relative to y_wall and delta_u check positioning of u_fields {avg result is good}
##    - adjust y offset if needed in processing step of that section
##    - stitch that at u resolution
##    - interpolate that on c resolution global mesh
##
## inter section offsets have to be correted before stitching! same x offsets for u and c mesh
##
##
## strategy for u field y wall offset
## find good y offset for flat region and adjust microstructure region to match nicely
##

## load avg velocity field per section (plus related mesh),
## find offset from velocity profile fit and cerrect offset before stitching.
## update avg processing and redo all the cases? or assume that x position of
## old and new processing is the same (which it probably is)
##

##

#### load processed data of all sections (from Ic and velocity records)
##pdat = load_all_2d (pdir, aid, it_A, it_C, it_M, it_X);
##
##
##  cd (fullpath);
##
##  c_tmp = load ("c.v7");
##  msh_c = c_tmp.c_msh;
##  cdata = c_tmp.c_dat;
##  masks_c = c_tmp.c_masks;
##  h_c = c_tmp.c_h;
##
##  u_tmp = load ("u.v7");
##  msh_u = u_tmp.u_msh;
##  udata = u_tmp.u_dat;
##  masks_u = u_tmp.u_masks;
##  h_u = u_tmp.u_h;
##
##  pp_tmp = load ("pp.v7");
##  pp = pp_tmp.pp;

#"2d_avg_uIc1"

##
## TODO: test for surface concentration influences
## vs_
## - phi
## - u_s
## - inclination of surface
## - film height
##

