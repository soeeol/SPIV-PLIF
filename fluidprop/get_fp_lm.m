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
## mf     .. mass fraction alcohol
## nref   .. refractive index
## rho    .. density in kg / m^3
## eta    .. dynamic viscisity in Pa s
## c_sat  .. oxygen saturation concentration in mg / l
## D_AB   .. diffusivity in m^2 / s
##
## Author: Sören J. Gerke
##

function [mf, nref, rho, eta, c_sat, D_AB] = get_fp_lm (pdir, id_f, T_get)
  mf = nref = rho = eta = csat = D_AB = [];
  switch (id_f)
    case "WG141"
      fname = "glycerol-water"; ext=[];
##      mf_PT_match = 5.8162e-01; # see ("ri_matching_regression.m")  mean of the 3 temperatures
##      mf_PT_match = 5.8131e-01; # see ("ri_matching_regression.m") 25°C
      ri_PDMS = ri_match_PDMS_ri_T (T_get);
      mf = ri_match_PT_mf_from_ri (ri_PDMS, T_get);
##      nref = get_ri_matching_tab (pdir, fname="glycerol-water", pname="ri", T_get, mf=mf_PT_match, ext=[]); # interpolated from experimental data
      nref = ri_match_PT_ri_from_mf (mf, T_get); # from exp. fit function
      rho = get_fp_tab (pdir, fname, pname="rho", T_get, mf, ext);
      eta = get_fp_tab (pdir, fname, pname="eta", T_get, mf, ext);
      c_sat = 3.6493e+00; # see "a_oxygen_diffusivity_solubility.m"
## TODO f(T, mf)      D_AB = get_fp_tab (pdir, fname, pname="D", T_get, mf, ext);
## TODO f(T, mf)      c_sat = get_fp_tab (pdir, fname, pname="c_sat", T_get, mf, ext);
##      D_AB.corr ... # see "a_oxygen_diffusivity_solubility.m"
      S = load ("-v7", [pdir.analyzed "a_2d_diff/DIFF_M13_A15_T25_WG141_M_G002_X_Z/D.v7"]);
      D_AB.PLIF1 = S.D_fit_1(1); # this studies result first lin. trend ("a_diff_main.m")
      D_AB.PLIF2 = S.D_fit_2(1); # this studies result 2nd lin. trend ("a_diff_main.m")
    case "WG139"
      w = 0.44; # TODO
      nref = get_ri_matching_tab (pdir, fname="glycerol-water", pname="ri", T_get, mf=w, ext=[]);
      rho = get_fp_tab (pdir, fname, pname="rho", T_get, mf, ext);
      eta = get_fp_tab (pdir, fname, pname="eta", T_get, mf, ext);
      c_sat = 4.75; # see "a_oxygen_diffusivity_solubility.m"
    case "WP141"
      mf_PD_match = 7.2519e-01; # see ("ri_matching_regression.m")
      ref = get_ri_matching_tab (pdir, fname="propylene glycol-water", pname="ri", T_get, mf=mf_PD_match, ext=[]);
      rho = get_fp_tab (pdir, fname, pname="rho", T_get, mf, ext);
      eta = get_fp_tab (pdir, fname, pname="eta", T_get, mf, ext);
      c_sat = 7.1437e+00; # see "a_oxygen_diffusivity_solubility.m"
      S = load ("-v7", [pdir.analyzed "a_2d_diff/DIFF_M26_A15_T25_WP141_M_G002_X_Z/D.v7"]);
      D_AB.PLIF1 = S.D_fit_1(1); # this studies result first lin. trend ("a_diff_main.m")
      D_AB.PLIF2 = S.D_fit_2(1); # this studies result 2nd lin. trend ("a_diff_main.m")
    otherwise
      error ("liquid parameters unknown");
  endswitch
endfunction
