##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## return liquid mixuture properties
##
## input:
## pdir   .. projects base directory
## id_f   .. identifier of liquid mixture
## T_get  .. temperature in K
##
## output:
## w      .. mass fraction alcohol
## n      .. refractive index
## rho    .. density in kg / m^3
## eta    .. dynamic viscisity in Pa s
## c_sat  .. oxygen saturation concentration in mg / l
## D_AB   .. diffusivity in m^2 / s
##
## Author: Sören J. Gerke
##

function [w, n, rho, eta, c_sat, D_AB] = get_fp_lm (pdir, id_f, T_get)
  w = n = rho = eta = c_sat = D_AB = [];
  calib_w = load ([pdir.analyzed "a_ri-matching/ri_matching_calibration.txt"]);
  switch (id_f)
    case "WG141"
      fname = "glycerol-water"; ext = [];
      n = n_PDMS = ri_PDMS_T (T_get, calib_w.fit_n_PDMS.c); # see a_ri_matching.m
      w = ri_matching_mf (n_PDMS, T_get, calib_w.fit_n_PT.c);
      rho = get_fp_tab (pdir, fname, pname="rho", T_get, w, ext);
      eta = get_fp_tab (pdir, fname, pname="eta", T_get, w, ext);
      c_sat = 3.6399; # see "a_oxygen_diffusivity_solubility.m"
## TODO f(T, w)      D_AB = get_fp_tab (pdir, fname, pname="D", T_get, w, ext);
## TODO f(T, w)      c_sat = get_fp_tab (pdir, fname, pname="c_sat", T_get, w, ext);
##      D_AB.corr ... # see "a_oxygen_diffusivity_solubility.m"
      S = load ("-v7", [pdir.analyzed "a_2d_diff/DIFF_M13_A15_T25_WG141_M_G002_X_Z/D.v7"]);
      D_AB.PLIF1 = S.D_fit_1(1); # this studies result first lin. trend ("a_diff_main.m")
      D_AB.PLIF2 = S.D_fit_2(1); # this studies result 2nd lin. trend ("a_diff_main.m")
      ## diffusivity temperature correction
      eta_meas = get_fp_tab (pdir, fname="glycerol-water", pname="eta", T_meas=273.15+25.12, w, ext=[]);
      D_AB.PLIF1 = corr_D_T_eta (eta_meas, T_meas, D_AB.PLIF1, eta, T_get);
      D_AB.PLIF2 = corr_D_T_eta (eta_meas, T_meas, D_AB.PLIF2, eta, T_get);
    case "WG139"
      w = 0.44; # TODO
##      n = get_ri_matching_tab (pdir, fname="glycerol-water", pname="ri", T_get, w=w, ext=[]);
      n = ri_matching_ri (w, T_get, fit_n_PT.c);
      rho = get_fp_tab (pdir, fname, pname="rho", T_get, w, ext);
      eta = get_fp_tab (pdir, fname, pname="eta", T_get, w, ext);
      c_sat = 4.75; # see "a_oxygen_diffusivity_solubility.m"
    case "WP141"
      fname = "propylene glycol-water"; ext = [];
      n = n_PDMS = ri_PDMS_T (T_get, calib_w.fit_n_PDMS.c); # see a_ri_matching.m
      w = ri_matching_mf (n_PDMS, T_get, calib_w.fit_n_PD.c);
      rho = get_fp_tab (pdir, fname, pname="rho", T_get, w, ext);
      eta = get_fp_tab (pdir, fname, pname="eta", T_get, w, ext);
      c_sat = 7.1437; # see "a_oxygen_diffusivity_solubility.m"
      S = load ("-v7", [pdir.analyzed "a_2d_diff/DIFF_M26_A15_T25_WP141_M_G002_X_Z/D.v7"]);
      D_AB.PLIF1 = S.D_fit_1(1); # this studies result first lin. trend ("a_diff_main.m")
      D_AB.PLIF2 = S.D_fit_2(1); # this studies result 2nd lin. trend ("a_diff_main.m")
      ## diffusivity temperature correction
      eta_meas = get_fp_tab (pdir, fname, pname="eta", T_meas=273.15+24.98, w, ext=[]);
      D_AB.PLIF1 = corr_D_T_eta (eta_meas, T_meas, D_AB.PLIF1, eta, T_get);
      D_AB.PLIF2 = corr_D_T_eta (eta_meas, T_meas, D_AB.PLIF2, eta, T_get);
    otherwise
      error ("liquid parameters unknown");
  endswitch
endfunction
