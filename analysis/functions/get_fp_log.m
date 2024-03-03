##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## return fluid properties based on logged data
##
## output:
## w      .. mass fraction alcohol
## n      .. refractive index
## rho    .. density in kg / m^3
## eta    .. dynamic viscisity in Pa s
##
## Author: Sören J. Gerke
##

function fp = get_fp_log (pdir, id_exp)

  ## refractive index calibration
  calib_w = load ([pdir.analyzed "a_ri-matching/ri_matching_calibration.txt"]);

  switch (id_exp)
    case {"flat_WG141", "2DR10_WG141"}
      ## actual flat plate experimental values from logs
      ## logged fluid properties on 02.12.2021 (tracer + seeding inside)
      ## probe taken from 20 feeding liter tank
      T_test = 298.15;
      eta_test = 8.421e-3; # rotational viscometer
      nref_test = 1.41065; # Abbe refractometer, should be most reliable
      rho_test = 1147.7; # Coriolis MFM
      ## derive mass fraction from nref_test
      w_exp = ri_matching_mf (nref_test, T_test, calib_w.fit_n_PT.c);
      ##
      rho_tab = get_fp_tab (pdir, fname="glycerol-water", pname="rho", T_test, w_exp, ext=[]);
      eta_tab = get_fp_tab (pdir, fname, pname="eta", T_test, w_exp, ext);
      ##
      fp.n = nref_test;
      fp.w = w_exp;
      fp.rho = rho_tab;
      fp.eta = eta_tab;
      fp.eta = eta_test;
    otherwise
      error ("get_fp_exp: no data for id_exp");
  endswitch

endfunction
