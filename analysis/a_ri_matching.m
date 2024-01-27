##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## refractive index matching analysis for liquid mixture mass fraction calibration
##
## Author: Sören J. Gerke
##

save_dir = [pdir.analyzed "a_ri-matching/"];
mkdir (save_dir)
cd (save_dir)

load ("-v7", [pdir.fptab "fluidprop.v7"], "fp");
sname = {"GerkeSJexp"};

## measurement temperatures
T = [20 25 30] + 273.15;

## H2O measurements
ds = get_fp_dataset (fp, fname={"water"}, pname={"refractive-index"}, sname);
n_W_T = ds.data{1,2};

##
## PDMS measurements
##
ds = get_fp_dataset (fp, fname={"PDMS"}, pname={"refractive-index"}, sname);
n_PDMS_T20 = ds.data{1,2}(1,:);
n_PDMS_T25 = ds.data{1,2}(2,:);
n_PDMS_T30 = ds.data{1,2}(3,:);
n_PDMS_T = median ([n_PDMS_T20', n_PDMS_T25', n_PDMS_T30'], 1);
n_PDMS_T = mean ([n_PDMS_T20', n_PDMS_T25', n_PDMS_T30'], 1);
n_PDMS_T_mad = mad ([n_PDMS_T20', n_PDMS_T25', n_PDMS_T30'], 1);
meas_PDMS_RI = [T', n_PDMS_T', n_PDMS_T_mad', [n_PDMS_T20;n_PDMS_T25;n_PDMS_T30]];
write_series_csv ([save_dir "meas_PDMS_RI"], meas_PDMS_RI, [], "%01.05f")
## linear fit for reference refractive index
T_vec = vec (T .* ones(3));
n_vec = [n_PDMS_T20 n_PDMS_T25 n_PDMS_T30];
c_n_PDMS = polyfit (T_vec, n_vec, 1);
n_PDMS_T_fit = ri_PDMS_T (T, c_n_PDMS)
## coefficient of determination linear regression
[Rxy, ~] = corrcoef (T_vec, n_vec);
RSQ = Rxy(1,2) * Rxy(1,2)
##
fit_n_PDMS = [];
fit_n_PDMS.descption = {"linear regression for n (T), fitting the measured refractive index n of PDMS"};
fit_n_PDMS.valid = {"273.15 <= T / K <= 293.15"};
fit_n_PDMS.material = {"PDMS, Momentive RTV615 1:10"};
fit_n_PDMS.fun = {"polyfit ()"};
fit_n_PDMS.c = c_n_PDMS;
fit_n_PDMS.RSQ = RSQ;

##
## 1,2,3-Propantriol calibration measurements "PT"
##
ds = get_fp_dataset (fp, fname={"glycerol"}, pname={"refractive-index"}, sname);
n_PT_T = ds.data{1,2};
ds = get_fp_dataset (fp, fname={"glycerol-water"}, pname={"refractive-index"}, sname);
w_PT = ds.data{1,2};
n_PT_T20 = ds.data{1,3}(1,:);
n_PT_T25 = ds.data{1,3}(2,:);
n_PT_T30 = ds.data{1,3}(3,:);
ds = get_fp_dataset (fp, fname={"glycerol-water"}, pname={"dynamic-viscosity"}, sname);
eta_PT_T20 = ds.data{1,3}(1,:);
eta_PT_T25 = ds.data{1,3}(2,:);
eta_PT_T30 = ds.data{1,3}(3,:);
##
meas_PT_RI = [w_PT', n_PT_T20', n_PT_T25', n_PT_T30'];
meas_PT_eta = [w_PT', eta_PT_T20', eta_PT_T25', eta_PT_T30'];
write_series_csv ([save_dir "meas_PT_RI"], meas_PT_RI, {"mass fraction","n @20°C","n @25°C","n @30°C"}, "%01.05f")
write_series_csv ([save_dir "meas_PT_eta"], meas_PT_eta, {"mass fraction","eta in Pa s @20°C","eta in Pa s @25°C","eta in Pa s @30°C"}, "%01.05f")

##
## 1,2-Propandiol measurements "PD"
##
ds = get_fp_dataset (fp, fname={"propylene glycol"}, pname={"refractive-index"}, sname);
n_PD_T = ds.data{1,2};
ds = get_fp_dataset (fp, fname={"propylene glycol-water"}, pname={"refractive-index"}, sname);
w_PD = ds.data{1,2};
n_PD_T20 = ds.data{1,3}(1,:);
n_PD_T25 = ds.data{1,3}(2,:);
n_PD_T30 = ds.data{1,3}(3,:);
ds = get_fp_dataset (fp, fname={"propylene glycol-water"}, pname={"dynamic-viscosity"}, sname);
eta_PD_T20 = ds.data{1,3}(1,:);
eta_PD_T25 = ds.data{1,3}(2,:);
eta_PD_T30 = ds.data{1,3}(3,:);
##
meas_PD_RI = [w_PD', n_PD_T20', n_PD_T25', n_PD_T30'];
meas_PD_eta = [w_PD', eta_PD_T20', eta_PD_T25', eta_PD_T30'];
write_series_csv ([save_dir "meas_PD_RI"], meas_PD_RI, {"mass fraction","n @20°C","n @25°C","n @30°C"}, "%01.05f")
write_series_csv ([save_dir "meas_PD_eta"], meas_PD_eta, {"mass fraction","eta in Pa s @20°C","eta in Pa s @25°C","eta in Pa s @30°C"}, "%01.05f")

##
## literature data refractive index
##
ds = get_fp_dataset (fp, fname={"glycerol-water"}, pname={"refractive-index"}, {"GP1963etal"});
w_PT_LIT = ds.data{1,2};
n_PT_LIT = ds.data{1,3};
ds = get_fp_dataset (fp, fname={"propylene glycol-water"}, pname={"refractive-index"}, {"MacBethG1951"});
w_PD_LIT = ds.data{1,2};
n_PD_LIT = ds.data{1,3};

##
## non-linear regression
##

## fit model in mf and T dimension for viscosity
## -> power law works well for glycerol - water and ok for propylene glycol - water:
##    eta(w, T) = P1(T) * w^P2(T,w) + P3
function eta = eta_model_fun_eta_mf (w, T, p)
  eta = (p(1).*T + p(2)) .* w.^(p(3).*T + p(4).*w) + p(5);
endfunction

## fit model in mf and T dimension for refractive index
## -> linear for temperature and quadratic for mass fraction works well to
##    describe the refractive index calibration measurements for both liquids in
##    the relevant range of mass fractions
##function ri = ri_model_fun_ri_mf (w, T, p)
##  ri = (p(1).*T + p(2)) .* w + (p(3).*T + p(4));
##endfunction
##function w = ri_model_fun_w_ri (ri, T, p)
##  w = ( ri - (p(3).*T + p(4)) ) / ( p(1).*T + p(2) );
##endfunction

## prepare regression input
init_eta_PT = [-1e3, 2e5, 1e-3, 3, 1]';
init_n_PT = ones (5, 1);
init_eta_PD = [-1e3, 2e3, 1e-3, 3, 1]';
init_n_PD = ones (5, 1);
## ... including pure components measurements:
##N_PT = ones (1, numel(w_PT) + 2);
##x_PT = [[1 w_PT 0], [1 w_PT 0], [1 w_PT 0]];
##y_PT = [T(1).*N_PT, T(2).*N_PT, T(3).*N_PT];
##obs_n_PT = [[n_PT_T(1) n_PT_T20 n_W_T(1)], [n_PT_T(2) n_PT_T25 n_W_T(2)], [n_PT_T(3) n_PT_T30 n_W_T(3)]];
##N_PD = ones (1, numel(w_PD) + 2);
##x_PD = [[1 w_PD 0], [1 w_PD 0], [1 w_PD 0]];
##y_PD = [T(1).*N_PD, T(2).*N_PD, T(3).*N_PD];
##obs_n_PD = [[n_PD_T(1) n_PD_T20 n_W_T(1)], [n_PD_T(2) n_PD_T25 n_W_T(2)], [n_PD_T(3) n_PD_T30 n_W_T(3)]];
## ... model works best for the range of interest:
N_PT = ones (1, numel(w_PT));
x_PT = [w_PT, w_PT, w_PT];
y_PT = [T(1).*N_PT, T(2).*N_PT, T(3).*N_PT];
obs_n_PT = [n_PT_T20, n_PT_T25, n_PT_T30];
obs_eta_PT = [eta_PT_T20, eta_PT_T25, eta_PT_T30];
N_PD = ones (1, numel(w_PD));
x_PD = [w_PD, w_PD, w_PD];
y_PD = [T(1).*N_PD, T(2).*N_PD, T(3).*N_PD];
obs_eta_PD = [eta_PD_T20, eta_PD_T25, eta_PD_T30];
obs_n_PD = [n_PD_T20, n_PD_T25, n_PD_T30];
##
eta_fit_fun = @ (p, x, y) eta_model_fun_eta_mf (x, y, p);
##n_fit_fun = @ (p, x, y) ri_model_fun_ri_mf (x, y, p);
n_fit_fun = @ (p, x, y) ri_matching_ri (x, y, p);
##
settings = optimset ("MaxIter", 100, "TolFun", 1e-9);

[p_eta_PT, residuals_eta_PT, cvg_eta_PT, outp_eta_PT] = nonlin_residmin (@(p_eta_PT) (eta_fit_fun(p_eta_PT, x_PT, y_PT) - obs_eta_PT), init_eta_PT, settings);
[c_n_PT, residuals_n_PT, cvg_n_PT, outp_n_PT] = nonlin_residmin (@(c_n_PT) (n_fit_fun(c_n_PT,x_PT,y_PT) - obs_n_PT), init_n_PT, settings);

[p_eta_PD, residuals_eta_PD, cvg_eta_PD, outp_eta_PD] = nonlin_residmin (@(p_eta_PD) (eta_fit_fun(p_eta_PD, x_PD, y_PD) - obs_eta_PD), init_eta_PD, settings);
[c_n_PD, residuals_n_PD, cvg_n_PD, outp_n_PD] = nonlin_residmin (@(c_n_PD) (n_fit_fun(c_n_PD,x_PD,y_PD) - obs_n_PD), init_n_PD, settings);

## standard error of regression
SER_eta_PT = sqrt (1 / (numel (residuals_eta_PT) - numel (p_eta_PT)) * sum (residuals_eta_PT.^2)) # Pas
SER_n_PT = sqrt (1 / (numel (residuals_n_PT) - numel (c_n_PT)) * sum (residuals_n_PT.^2)) # -
SER_eta_PD = sqrt (1 / (numel (residuals_eta_PD) - numel (p_eta_PD)) * sum (residuals_eta_PD.^2)) # Pas
SER_n_PD = sqrt (1 / (numel (residuals_n_PD) - numel (c_n_PD)) * sum (residuals_n_PD.^2)) # -

fit_n_PT = [];
fit_n_PT.descption = {"regression for n (T,w), fitting the measured refractive index depending on temperature T and mass fraction w"};
fit_n_PT.valid = {"273.15 <= T / K <= 293.15, 0.4 <= w <= 0.75"};
fit_n_PT.material = {"aqueous glycerol"};
fit_n_PT.fun = {"ri_matching_ri (), ri_matching_mf ()"};
fit_n_PT.c = c_n_PT
fit_n_PT.SER = SER_n_PT

fit_n_PD = [];
fit_n_PD.descption = {"regression for n (T,w), fitting the measured refractive index depending on temperature T and mass fraction w"};
fit_n_PD.valid = {"273.15 <= T / K <= 293.15, 0.5 <= w <= 0.95"};
fit_n_PD.material = {"aqueous propylene glycol"};
fit_n_PD.fun = {"ri_matching_ri (), ri_matching_mf ()"};
fit_n_PD.c = c_n_PD
fit_n_PD.SER = SER_n_PD

## save fit parameters
save -text "ri_matching_calibration.txt" fit_n_PDMS fit_n_PT fit_n_PD

##
## refractive index matching: mass fraction of liquid mixtures
##
w_match_T_PT = ri_matching_mf (n_PDMS_T_fit, T, c_n_PT)
w_match_PT = mean (w_match_T_PT) # 0.5833
mad (w_match_T_PT)
n_match_PT = ri_matching_ri (w_match_PT, T(2), c_n_PT)
##
w_match_T_PD = ri_matching_mf (n_PDMS_T_fit, T, c_n_PD)
w_match_PD = mean (w_match_T_PD) # 0.7241
mad (w_match_T_PD)
n_match_PD = ri_matching_ri (w_match_PD, T(2), c_n_PD)
## change of matching mass fraction for +/- 5 K deviation from 25 °C in %
w_match_T_PT / w_match_T_PT(2) * 100 - 100 # +/- 0.66 %
w_match_T_PD / w_match_T_PD(2) * 100 - 100 # +/- 0.1 %

##
## sensitivity of refractive index and of mass fraction estimate
##

function dndw = ri_mm_dw (w, T, p)
  dndw = (p(1) .*T + p(2)) .* (2 .* p(3) * w + p(4));
endfunction

function dndT = ri_mm_dT (w, p)
  dndT = p(1) .* (p(3).*w.^2 + p(4).*w + p(5));
endfunction

## sensitivity of liquid mixtures refractive index to temperature change
DndT_PT = ri_mm_dT (w_match_PT, c_n_PT)
DndT_PD = ri_mm_dT (w_match_PD, c_n_PD)
## sensitivity of PDMS refractive index to temperature change
DndT_PDMS = c_n_PDMS(1)

## sensitivtity of mass fraction estimate from calibration experiment at matching point
T_fit = [293.15:0.1:303.15];
n_fit = [-0.02:0.001:0.02] + n_PDMS_T_fit(2);
[DwDn DwDT] = gradient (ri_matching_mf (n_fit, T_fit', c_n_PT), n_fit, T_fit);
DwDn_PT = DwDn(T_fit==T(2),n_fit==n_PDMS_T_fit(2))
DwDT_PT = DwDT(T_fit==T(2),n_fit==n_PDMS_T_fit(2))
[DwDn DwDT] = gradient (ri_matching_mf (n_fit, T_fit', c_n_PD), n_fit, T_fit);
DwDn_PD = DwDn(T_fit==T(2),n_fit==n_PDMS_T_fit(2))
DwDT_PD = DwDT(T_fit==T(2),n_fit==n_PDMS_T_fit(2))

## temperature deviation
Delta_T = 0.1; # K

## precision of refractometer given is +/- 2e-5
Delta_n_meas = 1 * 2e-5;

100 * Delta_n_meas * DwDn_PT # +/- 0.014 %
100 * Delta_n_meas * DwDn_PD # +/- 0.023 %

## error of mass fraction estimate from liquid mixtures refractive index measurements
Delta_n_PT = mad (residuals_n_PT); #mean (mean (abs (meas_PT_RI(:,2:end) - ri_matching_ri (meas_PT_RI(:,1), T, c_n_PT))))
Delta_n_PD = mad (residuals_n_PD); #mean (mean (abs (meas_PD_RI(:,2:end) - ri_matching_ri (meas_PD_RI(:,1), T, c_n_PD))))
Delta_w_PT = 100 * DwDn_PT * Delta_n_PT + 100 * DwDT_PT * Delta_T # +/- 0.1 %
Delta_w_PD = 100 * DwDn_PD * Delta_n_PD + 100 * DwDT_PD * Delta_T # +/- 0.18 %

## error of mass fraction estimate from PDMS refractive index measurement
Delta_n_PDMS = mean (n_PDMS_T_mad)
Delta_w_PT_PDMS = 100 * DwDn_PT * Delta_n_PDMS + 100 * DwDT_PT * Delta_T # +/- 0.09 %
Delta_w_PD_PDMS = 100 * DwDn_PD * Delta_n_PDMS + 100 * DwDT_PD * Delta_T # +/- 0.16 %


## dyn. viscosity estimated from measurements
eta_PT_match_T = eta_model_fun_eta_mf (w_match_PT, T, p_eta_PT)
eta_PD_match_T = eta_model_fun_eta_mf (w_match_PD, T, p_eta_PD)
eta_PT_match_T = eta_model_fun_eta_mf (w_match_PT, T, p_eta_PT)
eta_PD_match_T = eta_model_fun_eta_mf (w_match_PD, T, p_eta_PD)
# test getting fluid properties for the ri-matching fluids
[w, n, rho, eta, c_sat, D_AB] = get_fp_lm (pdir, "WG141", 298.15)
[w, n, rho, eta, c_sat, D_AB] = get_fp_lm (pdir, "WP141", 298.15)


##
## plots
##

## PDMS n vs. T
fh = figure (); hold on;
plot (T_vec, n_vec, "k*;exp.;")
plot (T, n_PDMS_T, "rs-;mean;")
plot (T, n_PDMS_T_fit, "b-;fit;")
xlabel ("T in K")
ylabel ("refractive index PDMS in -")
print (fh, "-dpng", "-color", [save_dir "n_T_PDMS"]);
##
write_series_csv ([save_dir "fit_RI_PDMS_meas"], [vec([T(1) T(2) T(3)].*[1 1 1]') [n_PDMS_T20 n_PDMS_T25 n_PDMS_T30]'], {"T","n_PDMS_T20 n_PDMS_T25 n_PDMS_T30"}, "%01.05f")
write_series_csv ([save_dir "fit_RI_PDMS_fit"], [[T(1)-1 T(3)+1]' [ri_PDMS_T([T(1)-1 T(3)+1], c_n_PDMS)]'], {"T","n fit"}, "%01.05f")

## 2d visualization of ri-matching
w_fit = [0:0.001:1];
n_fit_PT = [n_fit_fun(c_n_PT, w_fit', T)];
n_fit_PD = [n_fit_fun(c_n_PD, w_fit', T)];
##
write_series_csv ([save_dir "RI_match"], [ w_match_PT 0.0 w_match_PD 0.0; w_match_PT n_match_PT  w_match_PD n_match_PD ], {"w_match_PT","n_match","w_match_PD","n_match"}, "%01.05f")
write_series_csv ([save_dir "RI_match_T"], [ w_match_T_PT' w_match_T_PD' n_PDMS_T_fit'], {"w_match_T_PT","w_match_T_PD","n_match_T"}, "%01.05f")
write_series_csv ([save_dir "meas_PDMS_RI_plot"], [[0 1]' n_PDMS_T_fit.*[1 1]'], {" ","n @20°C","n @25°C","n @30°C"}, "%01.05f")
write_series_csv ([save_dir "fit_PT_RI_plot"], [[w_fit' n_fit_PT]], {"mass fraction","n @20°C","n @25°C","n @30°C"}, "%01.05f")
write_series_csv ([save_dir "fit_PD_RI_plot"], [[w_fit' n_fit_PD]], {"mass fraction","n @20°C","n @25°C","n @30°C"}, "%01.05f")
##
fh = figure (); hold on;
plot ([0 w_PT 1], [n_W_T(1) n_PT_T20 n_PT_T(1)], "k*;exp. glycerol;")
plot ([0 w_PT 1], [n_W_T(2) n_PT_T25 n_PT_T(2)], "k*;exp. glycerol;")
plot ([0 w_PT 1], [n_W_T(3) n_PT_T30 n_PT_T(3)], "k*;exp. glycerol;")
plot (w_fit, n_fit_PT(:,1), "b-")
plot (w_fit, n_fit_PT(:,2), "b-")
plot (w_fit, n_fit_PT(:,3), "b-")
plot ([0 w_PD 1], [n_W_T(1) n_PD_T20 n_PD_T(1)], "kd;exp. propylene glycol;")
plot ([0 w_PD 1], [n_W_T(2) n_PD_T25 n_PD_T(2)], "kd;exp. propylene glycol;")
plot ([0 w_PD 1], [n_W_T(3) n_PD_T30 n_PD_T(3)], "kd;exp. propylene glycol;")
plot (w_fit, n_fit_PD(:,1), "b-")
plot (w_fit, n_fit_PD(:,2), "b-")
plot (w_fit, n_fit_PD(:,3), "b-")
plot ([0 1], n_PDMS_T_fit(1)*[1 1], "k-;PDMS @T = 20°C;")
plot ([0 1], n_PDMS_T_fit(2)*[1 1], "k-;PDMS @T = 25°C;")
plot ([0 1], n_PDMS_T_fit(3)*[1 1], "k-;PDMS @T = 30°C;")
plot (w_match_PT*[1 1], [1.33 n_PDMS_T_fit(2)], "b--")
plot (w_match_PD*[1 1], [1.33 n_PDMS_T_fit(2)], "b--")
plot (w_match_T_PT, n_PDMS_T_fit, "bx")
plot (w_match_T_PD, n_PDMS_T_fit, "bx")
plot (w_PD_LIT, n_PD_LIT, "rs;lit. propylene glycol 25 °C;")
plot (w_PT_LIT, n_PT_LIT, "rs;lit. glycerol 20 °C;")
xlabel("mass fraction alcohol")
ylabel("refractive index")
legend ("autoupdate", "off", "location", "northwest")
print (fh, "-dpng", "-color", [save_dir "ri_matching"]);

## 3d visualization of ri matching calibration
n_pts = 21;
T_eq = linspace (min(T)-1, max(T)+1, n_pts);
n_PDMS_T_eq = ri_PDMS_T (T_eq, c_n_PDMS);
w_PT_eq = ri_matching_mf (n_PDMS_T_eq, T_eq, c_n_PT);
w_PD_eq = ri_matching_mf (n_PDMS_T_eq, T_eq, c_n_PD);
##
[XX_PT, YY_PT] = meshgrid (linspace (min(w_PT)-0.025, max(w_PT)+0.01, n_pts), T_eq);
[XX_PT, YY_PT] = meshgrid (linspace (0.3, 1, n_pts), T_eq);
[XX_PD, YY_PD] = meshgrid (linspace (min(w_PD)-0.025, max(w_PD)+0.01, n_pts), T_eq);
[XX_PD, YY_PD] = meshgrid (linspace (0.4, 1, n_pts), T_eq);
##
write_series_csv ([save_dir "surf_meas_PT_n"], [x_PT' y_PT' obs_n_PT'], {"w","T","n"}, "%01.05f")
write_series_csv ([save_dir "surf_fit_PT_n"], [XX_PT(:) YY_PT(:) vec(n_fit_fun(c_n_PT,XX_PT,YY_PT))], {"w","T","n"}, "%01.05f")
write_series_csv ([save_dir "surf_meas_PD_n"], [x_PD' y_PD' obs_n_PD'], {"w","T","n"}, "%01.05f")
write_series_csv ([save_dir "surf_fit_PD_n"], [XX_PD(:) YY_PD(:) vec(n_fit_fun(c_n_PD,XX_PD,YY_PD))], {"w","T","n"}, "%01.05f")
write_series_csv ([save_dir "surf_match_PT"], [w_match_T_PT' T' n_PDMS_T_fit'], {"w_match_T_PT","T","n_match_T"}, "%01.05f")
write_series_csv ([save_dir "surf_match_PD"], [w_match_T_PD' T' n_PDMS_T_fit'], {"w_match_T_PD","T","n_match_T"}, "%01.05f")
write_series_csv ([save_dir "surf_eq_match_PT"], [w_PT_eq' T_eq' n_PDMS_T_eq'], {"w_PT_eq","T_eq","n_PDMS_T"}, "%01.05f")
write_series_csv ([save_dir "surf_eq_match_PD"], [w_PD_eq' T_eq' n_PDMS_T_eq'], {"w_PD_eq","T_eq","n_PDMS_T"}, "%01.05f")
##
fh = figure (); hold on; view ([-1.5 0.75 1])
stem3 (x_PT, y_PT, residuals_n_PT./obs_n_PT*100, "k;PT;");
stem3 (x_PD, y_PD, residuals_n_PD./obs_n_PD*100, "b;PD;"); # measurements
plot3 ([0.95 0.4], 293.15*[1 1], [0 0], "k-.")
plot3 ([0.95 0.4], 298.15*[1 1], [0 0], "k-.")
plot3 ([0.95 0.4], 303.15*[1 1], [0 0], "k-.")
legend ({"PT","PD"})
grid on
xlabel ("mass fraction alcohol")
ylabel ("temperature in K")
zlabel ("residuals n in %")
print (fh, "-dpng", "-color", [save_dir "ri_matching_residuals_n"]);
##
write_series_csv ([save_dir "fit_model_residulas_PT"], [[x_PT(1) y_PT(1) y_PT(9) y_PT(17) 0 0 0]; [[x_PT(1:8)]' reshape(y_PT,8,3) reshape([residuals_n_PT./obs_n_PT*100],8,3)]; [x_PT(end) y_PT(1) y_PT(9) y_PT(17) 0 0 0];], {"w PT","T","residuals n in %"}, "%01.06f")
write_series_csv ([save_dir "fit_model_residulas_PD"], [[x_PD(1) y_PD(1) y_PD(8) y_PD(15) 0 0 0]; [[x_PD(1:7)]' reshape(y_PD,7,3) reshape([residuals_n_PD./obs_n_PD*100],7,3)]; [x_PD(end) y_PD(1) y_PD(8) y_PD(15) 0 0 0];], {"w PD","T","residuals n in %"}, "%01.06f")


fh = figure (); hold on; view ([-1 1 1]);
plot3 (x_PT, y_PT, residuals_eta_PT./obs_eta_PT*100, "kx", "MarkerSize", 5, "LineWidth", 2);
plot3 (x_PD, y_PD, residuals_eta_PD./obs_eta_PD*100, "bx", "MarkerSize", 5, "LineWidth", 2);
grid on;
xlabel ("mass fraction alcohol")
ylabel ("temperature in K")
zlabel ("residuals eta in %")
print (fh, "-dpng", "-color", [save_dir "ri_matching_residuals_eta"]);


n_z0 = 1.37;
fh = figure (); hold on;
plot3 (x_PT, y_PT, obs_n_PT, "k*", "MarkerSize", 5, "LineWidth", 1.5); # measurements
view ([-1 1 1])
mesh (XX_PT, YY_PT, n_fit_fun(c_n_PT,XX_PT,YY_PT), "facecolor", "none", "edgecolor", [0 0 0], "LineWidth", 1);
plot3 (w_match_PT.*[1 1 1], T, n_PDMS_T, "b*", "MarkerSize", 5, "LineWidth", 1.5); # fit
T_line = linspace (min(T)-1, max(T)+1, 30);
w_line = w_match_PT .* ones (1, numel(T_line));
plot3 (w_line, T_line, n_fit_fun(c_n_PT,w_line,T_line), "g-", "LineWidth", 1.5); # fit
legend ("exp.", "eq. RI (w, T)", "matched w", "eq. RI (T) matched");#
legend ("autoupdate", "off");
grid on;
xlim ([0.35 0.8])
ylim ([292 305])
zlim ([n_z0 1.44])
plot3 (w_match_PT*[1 1], [306 T(1)], n_z0*[1 1], "b--")
plot3 ([0.35 w_match_T_PT(1)], T(1)*[1 1], n_z0*[1 1], "b--")
plot3 ([0.35 w_match_T_PT(2)], T(2)*[1 1], n_z0*[1 1], "b--")
plot3 ([0.35 w_match_T_PT(3)], T(3)*[1 1], n_z0*[1 1], "b--")
plot3 (w_match_PT*[1 1], T(1)*[1 1], [n_z0 n_PDMS_T(1)], "b--")
plot3 (w_match_PT*[1 1], T(2)*[1 1], [n_z0 n_PDMS_T(2)], "b--")
plot3 (w_match_PT*[1 1], T(3)*[1 1], [n_z0 n_PDMS_T(3)], "b--")
xlabel ("mass fraction glycerol")
ylabel ("temperature in K")
zlabel ("refractive index")
print (fh, "-dpng", "-color", [save_dir "n_PT_matching"]);


n_z0 = 1.38;
fh = figure (); hold on;
plot3 (x_PD, y_PD, obs_n_PD, "k*", "MarkerSize", 5, "LineWidth", 1.5); # measurements
view ([-1 1 1])
mesh (XX_PD, YY_PD, n_fit_fun(c_n_PD,XX_PD,YY_PD), "facecolor", "none", "edgecolor", [0 0 0], "LineWidth", 1);
plot3 (w_match_PD.*[1 1 1], T, n_PDMS_T, "b*", "MarkerSize", 5, "LineWidth", 1.5); # fit
T_line = linspace (min(T)-1, max(T)+1, 30);
w_line = w_match_PD .* ones (1, numel(T_line));
plot3 (w_line, T_line, n_fit_fun(c_n_PD,w_line,T_line),"g-","LineWidth",1.5); # fit
legend ("exp.", "eq. RI (w, T)", "matched w", "eq. RI (T) matched");
legend ("autoupdate", "off");
grid on;
xlim ([0.45 1])
ylim ([292 305])
zlim ([n_z0 1.44])
plot3 (w_match_PD*[1 1], [306 T(1)], n_z0*[1 1], "b--")
plot3 ([0.35 w_match_T_PD(1)], T(1)*[1 1], n_z0*[1 1], "b--")
plot3 ([0.35 w_match_T_PD(2)], T(2)*[1 1], n_z0*[1 1], "b--")
plot3 ([0.35 w_match_T_PD(3)], T(3)*[1 1], n_z0*[1 1], "b--")
plot3 (w_match_PD*[1 1], T(1)*[1 1], [n_z0 n_PDMS_T(1)], "b--")
plot3 (w_match_PD*[1 1], T(2)*[1 1], [n_z0 n_PDMS_T(2)], "b--")
plot3 (w_match_PD*[1 1], T(3)*[1 1], [n_z0 n_PDMS_T(3)], "b--")
xlabel ("mass fraction propylene glycol")
ylabel ("temperature in K")
zlabel ("refractive index")
print (fh, "-dpng", "-color", [save_dir "n_PD_matching"]);


fh = figure (); hold on;
plot3 (x_PT, y_PT, obs_eta_PT, "k*", "MarkerSize", 5, "LineWidth", 1.5); # measurements
view ([-1 1 1])
mesh (XX_PT, YY_PT, eta_fit_fun(p_eta_PT,XX_PT,YY_PT),"facecolor", "none", "edgecolor", [0 0 0], "LineWidth", 1);
plot3 (w_match_PT.*[1 1 1], T, eta_PT_match_T, "b*", "MarkerSize", 5, "LineWidth", 1.5); # fit
T_line = linspace (min(T)-1, max(T)+1, 30);
w_line = w_match_PT .* ones (1, numel(T_line));
plot3 (w_line, T_line, eta_fit_fun(p_eta_PT,w_line,T_line),"g-","LineWidth",1.5); # fit
legend ("meas", "dynvisc(w,T) - fit", "matched w", "eta match (T)");legend ("autoupdate", "off");
grid on;
xlim ([0.35 0.8])
ylim ([292 305])
zlim ([0 40e-3])
xlabel ("mass fraction glycerol")
ylabel ("temperature in K")
zlabel ("dyn. viscosity in Pa s")
plot3 (w_match_PT*[1 1], [306 T(1)], [0 0], "b--")
plot3 ([0.35 w_match_PT(1)], T(1)* [1 1], [0 0], "b--")
plot3 ([0.35 w_match_T_PT(2)], T(2)*[1 1], [0 0], "b--")
plot3 ([0.35 w_match_T_PT(3)], T(3)*[1 1], [0 0], "b--")
plot3 (w_match_PT*[1 1], T(1)*[1 1], [0 eta_PT_match_T(1)], "b--")
plot3 (w_match_PT*[1 1], T(2)*[1 1], [0 eta_PT_match_T(2)], "b--")
plot3 (w_match_PT*[1 1], T(3)*[1 1], [0 eta_PT_match_T(3)], "b--")
print (fh, "-dpng", "-color", [save_dir "eta_PT_matching"]);


fh = figure (); hold on;
plot3 (x_PD, y_PD, obs_eta_PD, "k*", "MarkerSize", 5, "LineWidth", 1.5); # measurements
view ([-1 1 1])
mesh (XX_PD, YY_PD, eta_fit_fun(p_eta_PD,XX_PD,YY_PD), "facecolor", "none", "edgecolor", [0 0 0], "LineWidth", 1);
plot3 (w_match_PD.*[1 1 1], T, eta_PD_match_T, "b*", "MarkerSize", 5, "LineWidth", 1.5); # fit
T_line = linspace (min(T)-1, max(T)+1, 30);
w_line = w_match_PD .* ones (1, numel(T_line));
plot3 (w_line, T_line, eta_fit_fun(p_eta_PD,w_line,T_line), "g-", "LineWidth", 1.5); # fit
legend ("meas", "dynvisc(w,T) - fit", "matched w", "eta match (T)");legend ("autoupdate", "off");
grid on;
caxis ([0 35e-3])
xlim ([0.45 1])
ylim ([292 305])
zlim ([0 50e-3])
xlabel ("mass fraction propylene glycol")
ylabel ("temperature in K")
zlabel ("dyn. viscosity in Pa s")
plot3 (w_match_PD*[1 1], [306 T(1)], [0 0], "b--")
plot3 ([0.35 w_match_T_PD(1)], T(1)*[1 1], [0 0], "b--")
plot3 ([0.35 w_match_T_PD(2)], T(2)*[1 1], [0 0], "b--")
plot3 ([0.35 w_match_T_PD(3)], T(3)*[1 1], [0 0], "b--")
plot3 (w_match_PD*[1 1], T(1)*[1 1], [0 eta_PD_match_T(1)], "b--")
plot3 (w_match_PD*[1 1], T(2)*[1 1], [0 eta_PD_match_T(2)], "b--")
plot3 (w_match_PD*[1 1], T(3)*[1 1], [0 eta_PD_match_T(3)], "b--")
print (fh, "-dpng", "-color", [save_dir "eta_PD_matching"]);

## result table
write_series_csv ([save_dir "results_ri-matching"], [[T' n_PDMS_T_fit' w_match_T_PT' w_match_T_PD'];], {"T","n PDMS","w PT","w PD"}, "%01.06f")

