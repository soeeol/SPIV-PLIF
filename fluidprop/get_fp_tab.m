##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## get fluid property of aqueous solution at (T, w)
##
## input:
## - fname .. fluid name {water, glycerol-water, propylene glycol-water}
## - pname .. property name {rho, eta}
##
## - T .. temperature in K
## - w .. mass fraction alcohol
##
## output:
## - out_ip .. density in kg / m^3 or dynamic viscosity in Pa s depending on pname
##
## Author: Sören J. Gerke
##

function out_ip = get_fp_tab (pdir, fname, pname, T, w, ext)

  if (isempty (ext))
    switch ([fname "_" pname])
      case {"water_rho"}
        out_ip = rho_PT_W_model (0, T);
      case {"water_eta"}
##        tabdat = load ("-v7", [pdir.fptab "IAPWS 2008_water_eta_tab.v7"]);
##        out_ip = interp1 (tabdat.T_ip, tabdat.eta_ip, T, "pchip", "extrap");
##        out_ip = eta_W_model (T);
        out_ip = eta_W_model_IAPWS (T);
      case {"glycerol-water_rho"}
        out_ip = rho_PT_W_model (w, T);
      case {"glycerol-water_eta"}
        tabdat = load ("-v7", [pdir.fptab "SegurJ1951_glycerol-water_eta_tab.v7"]);
        [XX, YY] = meshgrid (T, w);
        out_ip = interp_dat_ext (tabdat.T_ip, tabdat.mf_ip, tabdat.eta_ip, XX, YY);
      case "propylene glycol-water_rho"
        tabdat = load ("-v7", [pdir.fptab "GeorgeJ2003_propylene glycol-water_rho_tab.v7"]);
        [XX, YY] = meshgrid (T, w);
        out_ip = interp_dat_ext (tabdat.T_ip, tabdat.mf_ip, tabdat.rho_ip, XX, YY);
      case {"propylene glycol-water_eta"}
##        tabdat = load ("-v7", [pdir.fptab "KhattabIS2017_propylene glycol-water_eta_tab.v7"]); # values lower compared to literature, slightly different PD product?
##        tabdat = load ("-v7", [pdir.fptab "SunT2004_propylene glycol-water_eta_tab.v7"]);
        tabdat = load ("-v7", [pdir.fptab "GeorgeJ2003_propylene glycol-water_eta_tab.v7"]);
        [XX, YY] = meshgrid (T, w);
        out_ip = interp_dat_ext (tabdat.T_ip, tabdat.mf_ip, tabdat.eta_ip, XX, YY);
      otherwise
        error ("get_fp_tab: unknown / not implemented")
    endswitch
  endif

endfunction
