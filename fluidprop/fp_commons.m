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
w_N2 = n_N2 * mm_N2; # wt-%
w_O2 = n_O2 * mm_O2 / mm_air; # wt-%
## partial pressures
p_N2 = n_N2 * p_total;
p_O2 = n_O2 * p_total;

##
## critical point data
##
## CRC Handbook of Chemistry and Physics 95th Edition

## O2

## source?
##p_crit_O2 = 50.46; # bar
##TK_crit_O2 = 154.60; # K
##rho_crit_O2 = 427; # kg / m^3
##mv_crit_O2 = fp_mv_rho (mm_O2, rho_crit_O2); # m^3 / mol

## EOS, see CRC 6-22
p_crit_O2 = 5.043 * 1e6; # Pa
TK_crit_O2 = 154.58; # K
rho_crit_O2 = 436.1; # m^3 / mol
mv_crit_O2 = fp_mv_rho (mm_O2, rho_crit_O2); # 7.3375e-05 m^3 / mol

## experimental, see CRC 6-84, ref. 3
p_crit_O2 = 5.043 * 1e6; # Pa
TK_crit_O2 = 154.581; # K
mv_crit_O2 = 73 * 1e-6; # m^3 / mol
rho_crit_O2 = mm_O2 / mv_crit_O2; # 438.34 kg / m^3

## N2

## EOS, see CRC 6-22
p_crit_N2 = 3.396 * 1e6; # Pa
TK_crit_N2 = 126.19; # K
rho_crit_N2 = 313.3; # m^3 / mol
mv_crit_N2 = fp_mv_rho (mm_N2, rho_crit_N2); # 8.9413e-05 m^3 / mol

## experimental, see CRC 6-84, ref. 3
p_crit_N2 = 3.39 * 1e6; # Pa
TK_crit_N2 = 126.192; # K
mv_crit_N2 = 90 * 1e-6; # m^3 / mol
rho_crit_N2 = mm_N2 / mv_crit_N2; # 311.26 kg / m^3

## H2O - Water - W

## source?
##p_crit_W = 220.5; # bar
##TK_crit_W = 374.1; # K
##rho_crit_W = 322; # kg / m^3
##mv_crit_W = fp_mv_rho (mm_W, rho_crit_W); # 5.5948e-05 m^3 / mol

## CRC 6-9
p_crit_W = 22.064 * 1e6; # Pa
TK_crit_W = 373.946; # K
rho_crit_W = 322; # kg / m^3
mv_crit_W = fp_mv_rho (mm_W, rho_crit_W); # 5.5948e-05 m^3 / mol

## experimental, see CRC 6-85
p_crit_W = 22.06 * 1e6; # Pa
TK_crit_W = 373.12; # K
mv_crit_W = 56 * 1e-6; # m^3 / mol
rho_crit_W = mm_W / mv_crit_W; # 321.70 kg / m^3

## C3H803 - Glycerol - PT

## source?
##p_crit_PT = 66.9; # bar
##TK_crit_PT = 453; # °C
##rho_crit_PT = 361; # m^3 / mol
##mv_crit_PT = fp_mv_rho (mm_PT, rho_crit_PT); # 255e-6 m^3 / mol

## fit, CRC 6-160, ref. 6
TK_crit_PT = 800.00; # K
rho_crit_PT = 351.51; # m^3 / mol
mv_crit_PT = fp_mv_rho (mm_PT, rho_crit_PT); # 262e-6 m^3 / mol

## CRC 6-65, ref. 394
p_crit_PT = 7.6 * 1e6; # +/- 0.8 Pa
TK_crit_PT = 850; # +/- 9 K
mv_crit_PT = 251 * 1e-6; # +/- 15
rho_crit_PT = mm_PT / mv_crit_PT; # 366.89 kg / m^3

## C3H8O2 - 1,2-Propanediol - PD

## source?
##p_crit_PD = 60.8; # bar
##TK_crit_PD = 352; # °C
##rho_crit_PD = 321; # kg / m^3
##mv_crit_PD = fp_mv_rho (mm_PD, rho_crit_PD); # 237e-6; # m^3 / mol

## CRC 6-71, ref. 32
TK_crit_PD = 676; # +/- 1 K
p_crit_PD = 5.9 * 1e6; # +/- 0.8 Pa
mv_crit_PD = 237 * 1e-6; # +/- 15 m^3 / mol
rho_crit_PD = mm_PD / mv_crit_PD; # 303.17 kg / m^3

##
## molar volume at boiling point
##

##mv_O2 = 0.0256 * 1e-3; # m^3 / mol from Le Bas, see HimmelblauD1964
##mv_O2 = 0.0279 * 1e-3; # m^3 / mol experimental, see HimmelblauD1964
mv_O2 = fp_mv_bp_TC (mv_crit_O2);
mv_W = fp_mv_bp_TC (mv_crit_W);
mv_PD = fp_mv_bp_TC (mv_crit_PD);
mv_PT = fp_mv_bp_TC (mv_crit_PT);

##
## latent heat of vaporization in kJ / kg
##

L_PD = 914;
L_PT = 974;
L_W = 2256;
L_O2 = 214;
