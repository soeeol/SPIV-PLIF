##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## per section per frame PLIF analysis for
## - (dyn) normalized concentration field cn
## - (dyn) gas-liquid interface delta_u
## - (dyn) interface normal concentration profiles cp
## - (avg) concentration boundary layer thickness delta_c

##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ## ap parameters defining the analysis
  ap = [];
  ap.a_type = "a_2DR10_dyn_cn_cp"; # identifier of this analysis
  ap.p_type = "2d_dyn_Ic"; # for the analysis use the output of processing procedure "p_type"

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
##  it_M = 2:4
##  i_X = it_X = 1

  ap.dyn_cn_if_scmad_dev_max = 0.25; # used for threshold estimation for valid median interface deviation
  ap.dyn_cn_nt_max = 20; # limit number of valid single frame used for the analysis (valid: small deviation to median interface)
  ap.c_calib_sig_X = 0; # c calibration reference smoothing factor per ids_X (default)
  ap.c_isec_off_shift_lim = [0.25 0.1]; # [mm] intra section phi offset limits (default)
  ap.c_isec_rcurv_lim = 200; # [mm] threshold for "flat" interface curvature radius

  ## parameters for concentration transformation
  ap.c_method = "linear"; # method to transform fluorescence intensity to concentration ("linear" / "nonlinear" .. no impact on delta_c)
  ap.c_if_method = "calib"; # method to deal with fluorescence intensity decay at the interface ("calib" / "calib-if" .. high impact on delta_c)
  ap.c_calib_sig_X = 0; # c calibration reference smoothing factor per ids_X (default)
  ap.c_rm_bl = true; # remove cameras black level

  ## parameters for interface normal concentration profile estimation
  ap.cp_pd_M = [0.35 0.45 0.5 0.5]; # interface normal profile depth for each ids_M [mm]
  ## spline interface fit for stable interface normals (defaults)
  ap.cp_if_sfit_order = 3; # spline order
  ap.cp_if_sfit_sps = 9; # splines devisions per section (.. generally: increase for curved interface, decrease for flat)
  ap.cp_if_sn_idx_off = 20;

  ## prepare directories
  ap.date_str = datestr (now, "yyyymmdd");
  ap.result_dir = [pdir.analyzed ap.a_type "/"];

endif

## [10] normalized concentration field per x-section
if 0

  ## analysis identifier
  ap.sec_a_id = ["cn-" ap.c_method "_" ap.c_if_method]

  for i_X = it_X
    for i_M = it_M

      close all

      ap.i_M = i_M;
      ap.i_X = i_X;

      ## measurement identifier
      ap.id_meas = get_measid_ap (ap)

      ## path to store per section analysis results and plots to
      ap.save_dir_id = [ap.result_dir ap.id_meas "/" ap.date_str "_" ap.sec_a_id "/"];
      mkdir (ap.save_dir_id);

      c_calib_sig_X = [];
      switch (i_M)
        case 1
          c_calib_sig_X = [2 0 0 3] # for i_M=1&i_X=1 saturation recorded film was slightly thinner and thus of lower fluorescence
        case 2
          c_calib_sig_X = [0 0 0 0]
        case 3
          if (i_X==2)
            ap.cp_if_sfit_sps = 18;
          endif
          c_calib_sig_X = [0 0 0 0]
        case 4
          if (i_X==2)
            ap.cp_if_sfit_sps = 18;
          endif
          c_calib_sig_X = [0 0 0 0]
      endswitch
      if (! isempty (c_calib_sig_X))
        ap.c_calib_sig_X = c_calib_sig_X(i_X);
      endif

      [~] = a_dyn_cn_avg_calib (pdir, ap);

    endfor
  endfor

endif

## [20] per x-section interface normal concentration profiles
## shortcut to dyn statistics: interface normal profiles dyn analysis to
## analyze the time averaged peak cn shifted concentration profiles
if 1

  ## section analysis identifier
  ap.sec_a_id = ["cp-avg-" ap.c_method "_" ap.c_if_method]

  for i_M = it_M
    for i_X = it_X

      close all

      ap.i_M = i_M;
      ap.i_X = i_X;

      ## measurement identifier
      ap.id_meas = get_measid_ap (ap)

      ## path to store per section analysis results and plots to
      ap.save_dir_id = [ap.result_dir ap.id_meas "/" ap.date_str "_" ap.sec_a_id "/"];
      mkdir (ap.save_dir_id);

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
      ap.erfc_fit.dcds_idx = dcds_idx;
      ap.erfc_fit.idx_r = idx_r;
      ap.erfc_fit.sig = sig;

      [~] = a_dyn_cp_avg_delta_c (ap);

    endfor
  endfor

endif

