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
  i_M = it_M = 2
##  i_M = it_M = 1:2

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

  ## stitch descriptor

  ## avilable cn data
##  save -text "if_stats.txt" x if_stats
##  save -v7 "phi_avg.v7" c_msh phi_avg phi_des phi_sat x phi_avg_s phi_des_s phi_sat_s
##  save -v7 "cn_dyn.v7" c_msh cn_dyn cn_dyn_avg
##  save -v7 "delta_u_dyn.v7" x y_wall delta_u delta_u_fit delta_u_fit_avg
  sd{1}.sd_name = "cn_dyn_avg"
  sd{1}.a_dir = [pdir.analyzed "a_2DR10_dyn_cn_cp/"]
  sd{1}.a_id_sec = ["cn-" ap.c_method "_" ap.c_if_method]
  sd{1}.msh_fn = "cn_dyn.v7"
  sd{1}.msh_var = "c_msh"
  sd{1}.dat_fn{1} = "delta_u_dyn.v7"
  sd{1}.dat_var{1} = {"delta_u_fit_avg" "y_wall"}
  sd{1}.dat_fn{2} = "cn_dyn.v7"
  sd{1}.dat_var{2} = {"cn_dyn_avg"}

    ## avilable cp data
##  save -v7 "dyn_cp.v7" snp p_msh msh_n cp_n cp_s cp_b
##  save -v7 "avg_cp.v7" snp_o p_msh_o cp_nn_avg cp_n_avg
##  save -v7 "avg_delta_c.v7" x delta_c a_fit
  sd{2}.sd_name = "cp_dyn_avg"
  sd{2}.a_dir = [pdir.analyzed "a_2DR10_dyn_cn_cp/"]
  sd{2}.a_id_sec = ["cp-avg-" ap.c_method "_" ap.c_if_method]
  sd{2}.msh_fn = "avg_cp.v7"
  sd{2}.msh_var = "p_msh_o"
  sd{2}.dat_fn{1} = "avg_cp.v7"
  sd{2}.dat_var{1} = {"cp_n_avg"}
  sd{2}.dat_fn{2} = "avg_delta_c.v7"
  sd{2}.dat_var{2} = {"delta_c"}

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

  ## before inter section offset correction
  figure (); hold on;
  for i_M = it_M
    for i_X = it_X
      plot (x_MX{i_M,i_X}+ap.ids_X(i_X), delta_u_avg_MX{i_M,i_X}, "k");
      plot (x_MX{i_M,i_X}+ap.ids_X(i_X), y_wall_MX{i_M,i_X}, "k");
    endfor
  endfor
  axis image
  ylim ([0 max(delta_u_avg_MX{it_M(end),2})])

  ## helper to find systematic offset
  xoff = []
  for i_M = it_M
    for i_X = it_X(1:end-1)
      x_l = x_MX{i_M,i_X};
      x_r = x_MX{i_M,i_X+1};
      xpos_l = ap.ids_X(i_X);
      xpos_r = ap.ids_X(i_X+1);
      delta_u_l = delta_u_avg_MX{i_M,i_X};
      delta_u_r = delta_u_avg_MX{i_M,i_X+1};
      xoff(i_M,i_X) = calc_xsec_if_offset_x (x_l, x_r, delta_u_l, delta_u_r, xpos_l, xpos_r);
    endfor
  endfor
  xoff

## TODO: inter section offset settings
##  ap.sd.xoff_X = [+0.05 0 +0.2 +0.2] # mm
##  yoff = zeros (numel(ap.ids_M), numel(ap.ids_X));
##  yoff(1,4) = +0.005;
##  yoff(2,1) = +0.02;
##  yoff(2,4) = -0.005;
##  yoff(3,3) = +0.015;
##  yoff(3,4) = +0.015;
##  yoff(4,1) = -0.015;
##  yoff(4,3) = +0.01;
##  yoff(4,4) = +0.005;

  ## x sd
  i_sd = 2
  for i_M = it_M

    ap.i_M = i_M;
    ap.i_X = it_X;

    ## prepare result storage path
    ap.id_meas = get_measid_ap (ap)
    ap.save_dir_id = [ap.result_dir ap.id_meas "/" ap.date_str "_" "avg_stitch" "/"];
    mkdir (ap.save_dir_id);

    [msh_gl dat_gl] = stitch_msh_dat (sd{i_sd}, ap);

    if (i_sd==1)
      fh = figure ("position", [0,0,1920,1080]);
      plot_map_msh (msh_gl, dat_gl.cn_dyn_avg, fh)
      hold on;
      draw_cell (ap.ids_C{ap.i_C}, 0, 1);
      plot (msh_gl{1}(1,:), dat_gl.y_wall, "r");
      plot (msh_gl{1}(1,:), dat_gl.delta_u_fit_avg, "r");
      xlabel ("x* in mm");
      ylabel ("y in mm");
      axis image
      ylim ([-0.1 2.5])
      print (fh, "-dpng", "-color", ["-r" num2str(500)], [ap.save_dir_id sd{i_sd}.sd_name "_stitched.png"]);
    endif

    if (i_sd==2)

      fh = plot_map_msh (msh_gl, dat_gl.cp_n_avg, [])
      hold on;
      plot (msh_gl{1}(1,:), dat_gl.delta_c, "r");
      xlabel ("x* in mm");
      ylabel ("s_n in mm");
      print (fh, "-dpng", "-color", ["-r" num2str(500)], [ap.save_dir_id sd{i_sd}.sd_name "_stitched.png"]);

    endif

  endfor

endif


##
## TODO: test for surface concentration influences
## vs_
## - phi
## - u_s
## - inclination of surface
## - film height
##

