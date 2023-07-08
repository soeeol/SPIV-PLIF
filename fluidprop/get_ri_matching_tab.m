##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## load data table for fluid and property combination
## get data interpolated at T_get, mf_get from the refractive index calibration experiment
##
## Author: Sören J. Gerke
##

function out_ip = get_ri_matching_tab (pdir, fname, pname, T_get, mf_get, ext)
  if isempty(ext)
    switch ([fname "_" pname])
      case {"glycerol-water_eta"}
        tabdat = load ("-v7", [pdir.fptab "GerkeSJexp_glycerol-water_dynamic-viscosity_tab.v7"]);
      case {"glycerol-water_ri"}
        tabdat = load ("-v7", [pdir.fptab "GerkeSJexp_glycerol-water_refractive-index_tab.v7"]);
      case {"propylene glycol-water_eta"}
        tabdat = load ("-v7", [pdir.fptab "GerkeSJexp_propylene glycol-water_dynamic-viscosity_tab.v7"]);
      case {"propylene glycol-water_ri"}
        tabdat = load ("-v7", [pdir.fptab "GerkeSJexp_propylene glycol-water_refractive-index_tab.v7"]);
      otherwise
        error ("get_fp_tab: unknown")
    endswitch
    [XX, YY] = meshgrid (T_get, mf_get);
    out_ip = interp_dat_int (tabdat.T_ip, tabdat.mf_ip, tabdat.fp_ip, XX, YY);
  endif
endfunction
