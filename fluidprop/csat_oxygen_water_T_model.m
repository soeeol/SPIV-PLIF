##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## two models describing the solubility of oxygen in water depending on the temperature
##
## Author: Sören J. Gerke
##

function c_sol = csat_oxygen_water_T_model (T, p_O2)
## BensonBB1979etal
##
## k_x in atm; with f = k_x * x, here assuming ideal gas: p = k_x * x;
##
##
## MiyamotoH2014etal fit to measured values of RettichTR2000etal
##
## k_x in Pa; with p = k * x_x;
##
## !T < 373 K!
##
## x mole fraction of dissolved oxygen
## p_O2 - partial pressure of the gas in the vapor phase in Pa
## c_sol - solubility of oxygen O2 in water in mg / l

  run ("fp_commons.m")

  rho_W = rho_PT_W_model (0, T);
  rho_O2 = 101325 ./ (const_ug ./ mm_O2 .* T);

  ##fit_fun_c_T = "BensonBB1979etal"
  fit_fun_c_T = "MiyamotoH2014etal"

  switch fit_fun_c_T
    case {"BensonBB1979etal"}
      f_k_x = @(T) (e .^ (3.71814 + 5596.17 ./ T - 1049668 ./ T .^ 2)); # atm # BensonBB1979etal
      mf_sol = fp_mf_mx ((p_O2/101325) ./ f_k_x(T), mm_O2, mm_W); # g/g
      c_sol = 1e3 * fp_mc_mf (rho_W, rho_O2, mf_sol); # mg/l
  case {"MiyamotoH2014etal"}
      f_k_x = @(T) e.^(14.989460 + 5.742622e3 ./ T - 1.070683e6 ./ (T.^2)); # Pa # MiyamotoH2014etal
      mf_sol = fp_mf_mx ((p_O2) ./ f_k_x(T), mm_O2, mm_W); # g/g
      c_sol = 1e3 * fp_mc_mf (rho_W, rho_O2, mf_sol); # mg/l
  endswitch

endfunction
