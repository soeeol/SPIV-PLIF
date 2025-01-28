##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## laminar flat film theoretical solution


## [01] tabulate analytic solution for characetristic delta_u, u_s and delta_c versus inlet Reynolds number
ap.ids_A = [60]; # [°] inlination IDs
ap.a_type = "a_analytical_solutions";
ap.result_dir = [pdir.analyzed ap.a_type "/"];
ap.pos_ref_profile = "downstream";
ap.save_dir = [ap.result_dir "reference_flow_profile/" ap.pos_ref_profile "/"];
mkdir (ap.save_dir);

if 1

  switch (ap.pos_ref_profile)
    case {"microstructure"} # micro structure ref pos avg
      x_ref_prof = 0; # mm
    case {"full"} # whole film avg
      x_ref_prof = 8;
    case {"upstream"} # upstream profile
      x_ref_prof = -10;
    case {"downstream"} # downstream profile
      x_ref_prof = 10;
  endswitch

  ## fluid properties from exp log
  fp = get_fp_log (pdir, "flat_WG141");
  D_AB = [fp.D_AB_1 fp.D_AB_2]

  ## tabulate predicted delta_u_ref, u_s and delta_c depending on inlet Reynolds number
  ## Re | delta_u_ref in mm | u_s in m/s | delta_c in µm
  re = linspace (0, 50, 1001);
  delta_u = model_filmflow_laminar_deltau (fp.nu, deg2rad (ap.ids_A), re); # m
  u_s = model_filmflow_laminar_us (fp.nu, deg2rad (ap.ids_A), re); # m / s
  x_delta_c_tab = 1e-3 * (x_abs_meas + x_ref_prof) + x_off_inlet # in m; upstream profile position gas contact length
  delta_c = model_filmflow_laminar_deltac (x_delta_c_tab, D_AB(2), u_s); # m
  mfr = re * cell_width*1e-3 * fp.nu * fp.rho; # kg / s

  write_series_csv ([ap.save_dir "tab_eq_Re_deltau_us_deltac_mfr"], [re' 1e3*delta_u' u_s' 1e6*delta_c' 3600*mfr'], {"Re in -", "delta_u_ref in mm", "u_s in m/s", "delta_c in µm", "mfr in kg / h"}, []);

endif



if 1

  ## at the same Re number, influence of viscosity estimation
  visc_factor = 1.1

  u_s_2 = model_filmflow_laminar_us (fp.nu*visc_factor, deg2rad (ap.ids_A), re); # m / s
  delta_u_2 = model_filmflow_laminar_deltau (fp.nu*visc_factor, deg2rad (ap.ids_A), re); # m
  delta_c_2 = model_filmflow_laminar_deltac (x_delta_c_tab, D_AB(2), u_s_2); # m

  figure ()
  plot (re, 100 * (u_s_2-u_s) ./ u_s, "k")
  hold on
  plot (re, 100 * (delta_u_2-delta_u) ./ delta_u, "b")
  plot (re, 100 * (delta_c_2-delta_c) ./ delta_c, "r")

endif

