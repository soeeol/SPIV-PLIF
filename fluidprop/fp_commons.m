##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## some constants and some fluid properites of the pure components
##
## Author: Sören J. Gerke
##

## universal gas constant
const_ug = 8.31446261815324; # J / ( mol K )
## Avogadro constant
const_na = 6.02214076e23; # 1 / mol
## molar masses in kg/mol
mm_O2 = 31.999 * 1e-3;
mm_N2 = 28.013 * 1e-3;
mm_PT = 92.09 * 1e-3;
mm_PD = 76.095 * 1e-3;
mm_W = (15.9994 + 2*1.00794) * 1e-3;
mm_air = 28.9647 * 1e-3;
## total pressure
p_total = 101325; # Pa  = 1 atm
#### composition dry air
##n_N2 = 78.08e-2;
##n_O2 = 20.95e-2;
## composition synthetic air
n_N2 = 80e-2; # mol %
n_O2 = 20e-2; # mol %
mf_N2 = n_N2 * mm_N2; # wt-%
mf_O2 = n_O2 * mm_O2 / mm_air; # wt-%
## partial pressures
p_N2 = n_N2 * p_total;
p_O2 = n_O2 * p_total;

##
##
##
mv_O2 = 0.0256*1e-3; ## m^3 / mol from Le Bas method 1915 (see Reid et al. (1987)) see Gambil 1958 as well for values, group contribution
mv_O2 = 0.0279*1e-3; ## m^3 / mol exp.
##
## critical point data: CRC Handbook of Chemistry and Physics
##
## O2
TK_crit_O2 = 154.60; # K
p_crit_O2 = 50.46; # bar
rho_crit_O2 = 427; # kg/m^3
mv_crit_O2 = fp_mv_rho (mm_O2, rho_crit_O2); # m^3/mol
## C3H803
p_crit_PT = 66.9;
TK_crit_PT = 453+273.15;
rho_crit_PT = 361;
mv_crit_PT = 255e-6;
mv_crit_PT = fp_mv_rho (mm_PT, rho_crit_PT);
## H2O
p_crit_W = 220.5;
TK_crit_W = 374.1+273.15;
rho_crit_W = 322;
mv_crit_W = 56e-6;
mv_crit_W = fp_mv_rho (mm_W, rho_crit_W);
## C3H8O2
p_crit_PD = 60.8;
TK_crit_PD = 352+273.15;
rho_crit_PD = 321;
mv_crit_PD = 237e-6;
mv_crit_PD = fp_mv_rho (mm_PD, rho_crit_PD);
## latent heat of vaporization in kJ/kg
L_PD = 914;
L_PT = 974;
L_W = 2256;
L_O2 = 214;
