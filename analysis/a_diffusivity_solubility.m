##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## oxygen diffusity and solubility estimates for the water-alcohol mixtures from
## correlations and experiments
##
## Author: Sören J. Gerke
##
load ("-v7", [pdir.fptab "fluidprop.v7"], "fp");
run "fp_commons.m"

## range of experimental temperatures in K
TK = [20, 22.5, 25, 28, 30] + 273.15;
T_ref = 25 + 273.15;
## range of mass fraction alcohol in g / g
mf = linspace (0, 1, 201);
## range of dynamic viscosity in Pa s
eta_W = get_fp_tab (pdir, fname="water", pname="eta", T_get=TK, mf_get=1, ext=[]);
eta_PT_W = get_fp_tab (pdir, fname="glycerol-water", pname="eta", T_get=T_ref, mf_get=mf, ext=[]);
eta_PD_W = get_fp_tab (pdir, fname="propylene glycol-water", pname="eta", T_get=T_ref, mf_get=mf, ext=[]);
## density in kg / m^3
rho_W = rho_PT_W_model (0, TK);
rho_PT = rho_PT_W_model (1, TK);
rho_PD = get_fp_tab (pdir, fname="propylene glycol-water", pname="rho", T_get=TK, mf_get=1, ext=[]);

## O2 ideal gas
rho_O2 = 101325 ./ (const_ug ./ mm_O2 .* T_ref);

##
## oxygen diffusivity
##

## Stokes-Einstein equation
f_D_SE      = @ (T, r_A, eta_B) const_ug .* T ./ const_na ./ (6 .* pi .* eta_B .* r_A);
f_D_SE_slip = @ (T, r_A, eta_B) const_ug .* T ./ const_na ./ (4 .* pi .* eta_B .* r_A);
r_O2 = 1/2 * (mv_O2 / const_na) ^ (1/3); # rough oxygen molecular radius estimate
D_PT_SE = f_D_SE_slip (T=T_ref, r_O2, [eta_W(TK==T_ref), eta_PT_W(2:end)]');
D_PD_SE = f_D_SE_slip (T=T_ref, r_O2, eta_PD_W');

## Wilke-Chan correlation
##
## Wilke and Chang (1955)
##
## original Wilke-Chang formulation for D in cm^2/s:
## D = 7.4e-8 * T * sqrt (phi_B * mm_B) / (eta_B * mv_A^(0.6))
## T:K, mm_B:g/mol, eta_B:centi poise, mv_A:cm^3/mol
##
## scaling factor to SI-units:
K_SI = (1e-4) * sqrt(1e3) / (1e3) / (1e6)^(0.6); # cm^2 to m^2 * sqrt(g/mol) to sqrt(kg/mol) / mPas to Pas / cm^3/mol^0.6 to m^3/mol^0.6
f_DAB_WC = @ (K, phi_B, mm_B, eta_B, mv_A, T) 7.4e-8 .* K .* T ./ (mv_A.^(0.6)) .* sqrt(phi_B .* mm_B) ./ (eta_B);
# _A solute
# _B solvent
# phi association factor (emperical)
phi_PD = 1.2;## glycol phi = 1.2  # NogamiH1962et
phi_PT = 1.4;## glycerol phi = 1.4
phi_W  = 2.6;## water phi = 2.6
## data for Wilke and Chang correlation:
mm_PTW_solvent = fp_mm_mf (mf, mm_PT, mm_W);
mm_PDW_solvent = fp_mm_mf (mf, mm_PD, mm_W);

## PerkinsLR1969et: apply Wilke and Chang for solvent mixture
## Perkins and Geankoplis (1969), doi: 10.1016/0009-2509(69)80075-8.
## PD:
x_PD = fp_mx_mf (mf, mm_PD, mm_PDW_solvent);
## with mixing rule:
phi_mm_PDW = x_PD*phi_PD*mm_PD + (1-x_PD)*phi_W*mm_W;
D_Am_PD_WC = f_DAB_WC (K_SI, phi_B=1, mm_B=phi_mm_PDW', eta_B=eta_PD_W', mv_A=mv_O2, T=T_ref);
## PT:
x_PT = fp_mx_mf (mf, mm_PT, mm_PTW_solvent);
## with mixing rule:
phi_mm_PTW = x_PT*phi_PT*mm_PT + (1-x_PT)*phi_W*mm_W;
D_Am_PT_WC = f_DAB_WC (K_SI, phi_B=1, mm_B=phi_mm_PTW', eta_B=eta_PT_W', mv_A=mv_O2, T=T_ref);

## effective diffusivity from binary diff coeff for gas dissolved in mix of two solvents
## HimmelblauDM1964, Tang, "crude estimate"
## Himmelblau (1964), doi: 10.1021/cr60231a002.
## Dim * sqrt(eta_m) = x_2 * D_W * sqrt (eta_2) + x_3 * D_PT * sqrt (eta_3)

## correlation of SitaramanR1963etal
## Sitaraman et al. (1963), doi: 10.1021/je60017a017.
f_DAB_S = @ (Ls, L, mm_B, eta_B, mv_A, T) 8.75e-17 * ((Ls^1/3 .* mm_B.^1/2 .* T) ./ (eta_B .* mv_A.^0.5 .* L^0.3)) .^ 0.93;
## Ls solvent , L solute latent heat of vaporization in kJ/kg
## O2 6.82 kJ/mol, glycerol 85.8 kJ/mol, water 44.2 kJ/mol
D_PD_S = f_DAB_S (Ls=L_PD, L=L_O2, mm_PD, eta_PD_W(end), mv_O2, T_ref)
D_PT_S = f_DAB_S (Ls=L_PT, L=L_O2, mm_PT, eta_3=eta_PT_W(end), mv_O2, T_ref)
D_W_S = f_DAB_S (Ls=L_W, L=L_O2, mm_W, eta_2=eta_W(TK==T_ref), mv_O2, T_ref)
## mixing rule of PerkinsLR1969et
mre_PT = 0.5;
x_2 = fp_mx_mf (1-mf, mm_W, mm_PTW_solvent);
x_3 = fp_mx_mf (mf, mm_PT, mm_PTW_solvent);
eta_2 = eta_W(TK==T_ref);
eta_3 = eta_PT_W(end);
eta_m = [eta_W(TK==T_ref) eta_PT_W(2:end)];
Dim_S_PT = (x_2*D_W_S*eta_2.^mre_PT + x_3*D_PT_S*eta_3.^mre_PT) ./ eta_m.^mre_PT;
mre_PD = 0.66;
x_2 = fp_mx_mf (1-mf, mm_W, mm_PDW_solvent);
x_3 = fp_mx_mf (mf, mm_PD, mm_PDW_solvent);
eta_2 = eta_W(TK==T_ref)
eta_3 = eta_PD_W(end)
eta_m = eta_PD_W;
Dim_S_PD = (x_2*D_W_S*eta_2.^mre_PD + x_3*D_PD_S*eta_3.^mre_PD) ./ eta_m.^mre_PD;

## fully emperical corrleation of DiazM1987etal @ 25 °C
## Díaz et al. (1987), doi: 10.1080/00986448708911872.
## viscosity proportionality is relaxed vs. Stokes-Einstein type correlations
##f_DAB_DM = @ (mv_A, mv_B, eta_B) 1e-4 * 6.02e-5 * ( (mv_B*1e6)^0.36 / ((eta_B*1e3)^0.61 * (mv_A*1e6)^0.64) ); # simplyfied option
f_DAB_DM = @ (mv_A, mv_B, eta_B) 1e-4 * 6.02e-5 * ( (mv_B*1e6)^0.36 / ((eta_B*1e3)^0.61 * (mv_A*1e6)^(0.43*(mv_B*1e6)^0.13)) );
D_W_DM  = f_DAB_DM (mv_A=mv_O2, mv_B=mv_W, eta_B=eta_W(TK==T_ref))
D_PT_DM = f_DAB_DM (mv_A=mv_O2, mv_B=mv_PT, eta_B=eta_PT_W(end))
D_PD_DM = f_DAB_DM (mv_A=mv_O2, mv_B=mv_PD, eta_B=eta_PD_W(end))
##D_PD_DM = 5.9e-10;
D_WPT_DM  = f_DAB_DM (mv_A=mv_W, mv_B=mv_PT, eta_B=eta_PT_W(end))
## mixing rule of PerkinsLR1969et
Dim_DM_PT = (x_2 * D_W_DM * eta_2.^mre_PT + x_3 * D_PT_DM * eta_3.^mre_PT) ./ eta_m.^mre_PT;
Dim_DM_PD = (x_2 * D_W_DM * eta_2.^mre_PD + x_3 * D_PD_DM * eta_3.^mre_PD) ./ eta_m.^mre_PD;

## correlation Tyn Calus 1975
## ... note that this should not be used for viscous solvents, they consider above 20-30 mPas as viscous
## should apply for other ~ T / eta correlations then too !

## proportionality viscosity: Hiss & Cussler 1973:  ~ eta^(-2/3)

## water in glycerol comparison
##
#### f_DAB_S: water in glycerol ## 3.87e-12 vs exp D ~ 1.4e-11
##D_WPT_S = f_DAB_S (Ls=L_PT, L=L_W, mm_PT, eta_PT_W(end), mv_W, T=T_ref)
#### f_DAB_WC: water in glycerol ## 4.1896e-12 vs exp D ~ 1.4e-11
##f_DAB_WC (K_SI, phi_B=phi_PT, mm_B=mm_PT, eta_B=eta_PT_W(end), mv_A=mv_W, T=T_ref)

##
## experimental studies
##

## difusivity for O2 in water is well established, still results ranging from 1.87e-9 to 2.6e-9 m^2/s @ 25 °C
##
## Diffusivity of oxygen in water by StDenisCE1971et
## St-Denis and Fell (1971), doi: 10.1002/cjce.5450490632.
##D_W_exp_mean = [2.31e-9]; # @ 25 °C
D_W_exp_mean = 6.92e-10 * T_ref / eta_W((TK==T_ref)) * 1e-5;

## NogamiH1962et, polarographic, T = 295.65 K
## Nogami and Kato (1962), doi: 10.1248/yakushi1947.82.1_120.
## PT
ds = get_fp_dataset (fp, {"PT_W"}, {"D_AB"}, {"NogamiH1962et"});
T_NK = ds.data{1};
vf_PT_NK = ds.data{2}; # vol fraction glycerol
D_PT_NK = ds.data{3}; # diffusivity m^2/s
ds = get_fp_dataset (fp, {"PT_W"}, {"c_sat"}, {"NogamiH1962et"});
csat_PT_NK = ds.data{3}; # mol/l
##phi_PT_NK = [2.6 2.57 2.45 2.3 2.06]; # association factor Wilke Chang
## conversion
mf_PT_NK = vf_PT_NK.*rho_PT(TK==T_NK) ./ (vf_PT_NK.*rho_PT(TK==T_NK) + (1-vf_PT_NK).*rho_W(TK==T_NK));
eta_PT_NK = get_fp_tab (pdir, fname="glycerol-water", pname="eta", T_get=T_NK, mf_get=mf_PT_NK, ext=[]);
eta_PT_NK_T25 = get_fp_tab (pdir, fname="glycerol-water", pname="eta", T_get=T_ref, mf_get=mf_PT_NK, ext=[]);
D_PT_NK_T25 = corr_D_T_eta (eta_PT_NK, T_NK, D_PT_NK, eta_PT_NK_T25, T_ref); # temperature correction 25°C ... 7 to 15 % increase

csat_PT_NK = csat_PT_NK * mm_O2 * 1e6; # mg/l
##
## PD
ds = get_fp_dataset (fp, {"PD_W"}, {"D_AB"}, {"NogamiH1962et"});
vf_PD_NK = ds.data{2}; # vol fraction glycerol
D_PD_NK = ds.data{3}; # diffusivity m^2/s
ds = get_fp_dataset (fp, {"PD_W"}, {"c_sat"}, {"NogamiH1962et"});
csat_PD_NK = ds.data{3}; # mol/l
## conversion
mf_PD_NK = vf_PD_NK.*rho_PD(TK==T_NK) ./ (vf_PD_NK.*rho_PD(TK==T_NK) + (1-vf_PD_NK).*rho_W(TK==T_NK));
eta_PD_NK = get_fp_tab (pdir, fname="propylene glycol-water", pname="eta", T_get=T_NK, mf_get=mf_PD_NK, ext=[]);
eta_PD_NK_T25 = get_fp_tab (pdir, fname="propylene glycol-water", pname="eta", T_get=T_ref, mf_get=mf_PD_NK, ext=[]);
D_PD_NK_T25 = corr_D_T_eta (eta_PD_NK, T_NK, D_PD_NK, eta_PD_NK_T25, T_ref); # temperature correction 25°C ... 6 to 12 % increase
csat_PD_NK = csat_PD_NK * mm_O2 * 1e6; # mg/l

## JordanJ1956etal, polarographic, T = 25°C
## Jordan et al. (1956), doi: 10.1021/ja01594a015.
## TODO: move to fluidprop collection
##
## PT
mf_org_Jetal = [0	0.124	0.2	0.349	0.526	0.68 0.8	0.843	0.875	0.925]; # mass fraction glycerol
mfa_Jetal = mf_org_Jetal(1:end-3); # mass fraction glycerol
rho_Jetal = [1006	1032	1058	1077	1126	1163	1197	1209	1215	1228]; # kg/m^3
eta_Jetal = [0.89	1.27	1.75	2.24	4.85	10.8	29	45.9	56.9	106]; # m Pa s  ## with electrolyte and maximum suppressor
D_org_Jetal = 1e-5*[2.12 3.31 2.99 1.21 1.15 0.87 0.51 0.48	0.43 0.24]; # Fick diffusivity in cm/s
Da_org_Jetal =1e-5*[1.80 1.60 1.00 0.30 0.2  0.1  0.05]; # "association" diffusivity in cm/s
c_org_Jetal = [1.18	0.889	0.735	0.649	0.478	0.378	0.302	0.28	0.264	0.238]; # mmol/l @ 1atm O2  ## by Winkler method +/- 1 %
#
mf_Jetal = mf_org_Jetal;
D_Jetal = D_org_Jetal * 1e-4; # m^2 / s
Da_Jetal = Da_org_Jetal * 1e-4; # m^2 / s
c_mg_Jetal = c_org_Jetal * mm_O2 * p_O2/p_total * 1e3; # mg/l

## checking the slope of D vs. viscosity
mre_deta_NK = 1./(eta_PD_NK_T25).^mre_PD;
mre_deta_NK_PT = 1./(eta_PT_NK_T25).^0.9;
mre_deta_Jetal = 1./(eta_Jetal*1e-3).^mre_PT;
fh = figure (); hold on;
loglog (eta_PD_NK_T25, D_PD_NK_T25, "d;PD exp. NogamiH1962et;")
loglog (eta_PT_NK_T25, D_PT_NK_T25, "x;PT exp. NogamiH1962et;")
loglog (eta_Jetal*1e-3, D_Jetal, "s;PT exp. JordanJ1956etal;")
loglog (eta_PD_NK_T25, D_PD_NK_T25(1).*mre_deta_NK./mre_deta_NK(1), ";slope mre=0.66;")
loglog (eta_PT_NK_T25, D_PT_NK_T25(1).*mre_deta_NK_PT./mre_deta_NK_PT(1), ";slope mre=0.8;")
loglog (eta_Jetal*1e-3, (D_Jetal(2)+D_Jetal(1))/2.*mre_deta_Jetal./mre_deta_Jetal(1), ";slope mre=0.5;")
xlabel ("dynamic viscosity in Pa s")
ylabel ("diffusivity in m^2 / s")
box on
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [pdir.plot "diffusivity/" "loglog_D_eta"]);

##
## PLIF diffusion front measurements
##

## master thesis; 25°C; considered as rough estimate, variance was extremely high with convective currents
mf_M_PT = 0.59;
D_M_PT = 0.344e-9;
mf_M_PD = 0.89;
D_M_PD = 0.232e-9;
## own SPIV-PLIF cell measurement; results of "a_diffusion.m"
[mf_SJG_PT, ~, rho_PT_ri_match, ~, ~, D_AB] = get_fp_lm (pdir, "WG141", T_ref); # from RI-matching
D_SJG_PT_1 = D_AB.PLIF1;
D_SJG_PT_2 = D_AB.PLIF2;
##
[mf_SJG_PD, ~, rho_PD_ri_match, ~, ~, D_AB] = get_fp_lm (pdir, "WP141", T_ref); # from RI-matching
D_SJG_PD_1 = D_AB.PLIF1;
D_SJG_PD_2 = D_AB.PLIF2;

## other PLIF experiments; solvent: water - ethanol - glycerol !
##
## JimenezM2012etal
## Jimenez et al. (2012), doi: 10.1002/aic.13805.
## 50 W 20 E 30 PT mass fraction
mf_PT_J = [0 0.2 0.3];
rho_J = [970 1029 1042];
eta_J = [1 1.7 2.4] * 1e-3;
T_J = 20 + 273.15;
D_J_3 = [19 12 9.47] * 1e-10; # case 3
## correction with assuming that viscosity decreases accordingly to that of water ~ 11 %
D_J_3_T25 = corr_D_T_eta (eta_J, T_J, D_J_3, eta_J*(eta_W(TK==T_ref)/eta_J(1)), T_ref); # estimate for 25°C ... ~ 14 % increase

## KapoustinaV2019etal
## Kapoustina et al. (2019), doi: 10.1016/j.cherd.2019.04.022.
## 20 W 30 E 50 PT mass fraction
mf_PT_K = 0.50;
eta_K_T25 = 9.4e-3; # measured at 25°C
rho_K = 1053; # at 25°C
D_K = 3.76e-10; # at 28°C
T_K = 28 + 273.15; # T = 28°C
## correction with assuming that viscosity decreases accordingly to that of water: ~ 7 %
eta_W_K = get_fp_tab (pdir, fname="water", pname="eta",T_get=T_K, mf_get=1, ext=[]);
eta_K = eta_K_T25/(eta_W(TK==T_ref)/eta_W_K);
D_K_T25 = corr_D_T_eta (eta_K, T_K, D_K, eta_K_T25, T_ref); # estimate for 25°C ... ~ 7 % decrease

## comparison solute water in solvent pure glycerol
##
## ErricoG2004etal
## D‘Errico et al. (2004), doi: 10.1021/je049917u.
## @25°C
D_PT_DErrico = D_water_glycerol = 1.4e-11;

##
## oxygen solubility
##

## Own csat measurements with UMS micro sensor in stirred beaker glas
## unknown how to measure reliably in viscous media with UMS micro sensor
## min flow for non decreasing reading needed to compensate oxygen consumption of probe
## possibly more flow needed for viscous liquids
mf_PT_SJG = [0 0.2915	0.4369 0.5109 0.5880 0.6994 0.7991 0.9097];
c_mg_PT_SJG_min = [8.025 5.9 4.85 4.15 3.42 2.93 2.5 2.12];
mf_PD_SJG = [0 0.2 0.4 0.725 0.85 0.975];
c_mg_PD_SJG_min = [8.2 7 6 6.5 7.15 8.5];

## KutscheI1984etal
## Kutsche et al. (1984), doi: 10.1021/je00037a018.
mf_Ketal = 1e-2 * [1.31 2.52 3.74 6.16 9.18 12.21 15.23 18.26]; # glycerol
aipera0 = [0.971 0.967 0.953 0.931 0.902 0.873 0.844 0.815]; # - @ 25°C
a0 = 0.0284; # lit.; Bunsen coefficient in water
a0 = 0.0279; # exp.; Bunsen coefficient in water
ai = aipera0 * a0;
Hscp_0 = a0 * 1 / (const_ug * 273.15);
##1e6 * mm_O2 .* p_O2 .* Hscp_0
Hscp_i = ai * 1 / (const_ug * 273.15);
c_mg_Ketal = 1e3 * mm_O2 .* p_O2 .* ([Hscp_0 Hscp_i]); # mg/l
## K_b = 4.74

## RischbieterE1996etal
## Rischbieter et al. (1996), doi: 10.1021/je960039c.
## experimental, T = 303.2 K
m_org_Retal = [0 23.4 46.8 92.7 138 209 492 691]; # kg/m^3, kg of glycerol per m^3
H_org = 1e3*[88.5 90 94.5 98.8 99.7 107 149.7 203.5]; # Henry's Const in Pa m^3 mol^-1
mf_Retal = fp_mf_mc (rho_W(4), rho_PT(4), m_org_Retal); # mass fraction in g/g
c_mg_Retal = mm_O2*1e3 * p_O2 ./ H_org; # x from p = H x in mg/l
## model
f_Kn = @(bn, bgH2O, bgT, T) bn + bgH2O + bgT .* (T - 298.15);
f_c_mg_Retal = @ (csat_W, Kn, cn) csat_W ./ (10.^(Kn .* (cn)));
#### Glycerol O2 original publication data:
bn_PD = 4.47 * 1e-4; # glycol m^3 kg^-1  #
bn_PT = 5.25 * 1e-4; # glycerol m^3 kg^-1; substance specific parameter
##
bgH2O = 0 * 1e-4; # m^3 kg^-1
bgT = - 0.044 * 1e-4; # m^3 kg^-1 K^-1; gas specific parameter
##
csat_W = csat_oxygen_water_T_model (T=TK, p_O2) # mg/l
Kn_PT = f_Kn (bn_PT, bgH2O, bgT, T=TK);
Kn_PD = f_Kn (bn_PD, bgH2O, bgT, T=TK);

## solubility comparison for JordanJ1956etal
##bn = 5.25 * 1e-4; # gylcerol m^3 kg^-1  # matches the slope for higer mass fraction

## YamamotoH1994etal
## Yamamoto et al. (1994), doi: 10.1252/jcej.27.455.
## @ 25°C and 1 atm
ds = get_fp_dataset (fp, {"PD_W"}, {"c_sat"}, {"YamamotoH1994etal"});
##mx_PD_Y = ds.data{2}; # mole fraction propylene glycol
##L = ds.data{3}; # Ostwald coeff
mf_PD_Y = fp_mf_mx (ds.data{2}, mm_PD, mm_W);
vf_PD_Y = fp_mf_mx (ds.data{2}, mm_PD/rho_PD(TK==T_ref), mm_W/rho_W(TK==T_ref));
##vsol = 1 .* ds.data{3} .* (p_O2 / p_total);
## conversion with ideal gas law
rho_ideal = p_total / (const_ug / mm_air * T_ref);
c_mg_Y = 1e3 * ds.data{3} * mf_O2 * rho_ideal;

## oxygen in water
T_fit = [270:1:330];
## BensonBB1979etal
## Benson et al. (1979), doi: 10.1007/bf01033696.
##
## f = k * x, assuming ideal gas for air here: p = k x
##
## T < 373 K
## pp_O2 - partial pressure of the gas in the vapor phase in Pa
## c_sol - solubility of oxygen O2 in water in mg / l
f_k_x_Benson = @(T) (e .^ (3.71814 + 5596.17 ./ T - 1049668 ./ T .^ 2)); # atm
mf_sat_Benson_fit = fp_mf_mx ((p_O2/p_total) ./ f_k_x_Benson(T_fit), mm_O2, mm_W); # g/g
c_sat_Benson_fit = 1e3 * fp_mc_mf (rho_PT_W_model(0,T_fit), rho_O2, mf_sat_Benson_fit); # mg/l

## MiyamotoH2014etal IUPAC solubility series:
## Miyamoto et al. (2014), doi: 10.1063/1.4883876.
## tabulated smoothed values based on RettichTR2000etal
## Rettich et al. (2000), doi: 10.1006/jcht.1999.0581.
## oxygen in water
T_iupac = [273.15 278.15 283.15 288.15 293.15 298.15 303.15 308.15 313.15 318.15 323.15 328.15];
H_iupac = 1e9 * [2.5591 2.9242 3.2966 3.6708 4.0415 4.4038 4.7533 5.0859 5.3985 5.6882 5.9531 6.1917];
##(1 .* p_total) ./ H_tab # mole fraction at 1 atm
mf_O2_iupac = fp_mf_mx (mx=(p_O2 ./ H_iupac), mm_A=mm_O2, mm_B=mm_W); # g/g
c_sat_iupac = 1e3 * fp_mc_mf (rho_PT_W_model(0,T_iupac), rho_O2, mf_O2_iupac); # mg/l
## reported fit function
##ln k_x = 14.989460 + 5.742622e3 ./ T - 1.070683e6 ./ (T.^2);
f_k_x = @(T) e.^(14.989460 + 5.742622e3 ./ T - 1.070683e6 ./ (T.^2)); # k_x in Pa
mf_O2_iupac_fit = fp_mf_mx (p_O2 ./ f_k_x(T_fit), mm_O2, mm_W); # g/g
c_sat_iupac_fit = 1e3 * fp_mc_mf (rho_PT_W_model(0,T_fit), rho_O2, mf_O2_iupac_fit); # mg/l

## O2 in 1,2-propanediol solubility @ 25°C experimental by OsborneAD1965et
## Osborne and Porter (1965), http://www.jstor.org/stable/2415099
##cx = @(p) 6/300 * p + 0; # 1e4 moles/litre, p in mmHg
cx = @(p) 2.25e-6 * p + 0; # moles/litre, p in mmHg
csat_PD_O = 1e6 * cx (160) * mm_O2  # mg/l

## correction of the polarographic diffusivity measurements by solubility based on Ilkovic equation
##
##csat_W = c_mg_Jetal(1); # mg/l @ 25°C
csat_Jetal_corr = f_c_mg_Retal (csat_W(TK==T_ref), Kn_PT(TK==T_ref), fp_mc_mf(rho_W(TK==T_ref), rho_PT(TK==T_ref), mf_Jetal));
Dcorr_Jetal = 1 .* (csat_Jetal_corr ./ c_mg_Jetal).^2;
##
csat_NK_PT_corr = f_c_mg_Retal (csat_W(TK==T_NK), Kn_PT(TK==T_NK), fp_mc_mf(rho_W(TK==T_NK), rho_PT(TK==T_NK), mf_PT_NK));
Dcorr_NK_PT = 1 .* (csat_NK_PT_corr ./ csat_PT_NK ).^2;
csat_NK_PD_corr = f_c_mg_Retal (csat_W(TK==T_NK), Kn_PD(TK==T_NK), fp_mc_mf(rho_W(TK==T_NK), rho_PD(TK==T_NK), mf_PD_NK));
Dcorr_NK_PD = 1 .* (csat_NK_PD_corr ./ csat_PD_NK ).^2;
Dcorr_NK_PD = [Dcorr_NK_PD(1:7) median(Dcorr_NK_PD(1:7))*[1 1 1 1]];
## ri matching values for this study
c_sat_PT_RI_match = f_c_mg_Retal(csat_W(TK==T_ref), Kn_PT(TK==T_ref), fp_mc_mf(rho_W(TK==T_ref), rho_PT(TK==T_ref), mf_SJG_PT))
c_sat_PD_RI_match = mean ([c_mg_Y(8) c_mg_PD_SJG_min(4)])
mf_PD_Y(8)
mf_PD_SJG(4)

##
## plots and output
##

##
save_dir = [pdir.plot "diffusivity/"]
cd (save_dir)
run "a_diffusivity_solubility_plots_diff.m"
close all
##
save_dir = [pdir.plot "solubility/"]
cd (save_dir)
run "a_diffusivity_solubility_plots_sol.m"
close all
