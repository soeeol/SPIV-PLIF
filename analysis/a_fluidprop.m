## fluid properties
## measured and literature

save_dir = [pdir.analyzed "a_fluidprop/"];
mkdir (save_dir)
cd (save_dir)

load ("-v7", [pdir.fptab "fluidprop.v7"], "fp");
run ("fp_commons.m")

w_plot = 0:0.0125:1;
[w_PT_match, n_match, rho_PT_match, eta_PT_match, ~, ~] = get_fp_lm (pdir, "WG141", 298.15)
[w_PD_match, ~      , rho_PD_match, eta_PD_match, ~, ~] = get_fp_lm (pdir, "WP141", 298.15)
calib_w = load ([pdir.analyzed "a_ri-matching/ri_matching_calibration.txt"]);



##
## density aq. propylene glycol
##
fname = {"propylene glycol-water"}
pname = {"density"}
ds_rho_PD = get_fp_dataset (fp, fname, pname, [])
## lit. exp.
i = id_R1_PD = get_fp_src_idx (ds_rho_PD, "KhattabIS2017");
T_PD_R1 = ds_rho_PD(i).data{1};
w_PD_R1 = fp_mf_mx (ds_rho_PD(i).data{2}, mm_PD, mm_W);
rho_PD_R1_T20 = ds_rho_PD(i).data{3}(1,:);
rho_PD_R1_T25 = ds_rho_PD(i).data{3}(2,:);
rho_PD_R1_T30 = ds_rho_PD(i).data{3}(3,:);
##
i = id_R2_PD = get_fp_src_idx (ds_rho_PD, "NakanishiK1967");
T_PD_R2 = ds_rho_PD(i).data{1};
w_PD_R2 = ds_rho_PD(i).data{2};
rho_PD_R2_T25 = ds_rho_PD(i).data{3}(1,:);
##
i = id_R3_PD = get_fp_src_idx (ds_rho_PD, "GeorgeJ2003");
T_PD_R3 = ds_rho_PD(i).data{1};
w_PD_R3 = fp_mf_mx (ds_rho_PD(i).data{2}, mm_PD, mm_W);
rho_PD_R3_T25 = ds_rho_PD(i).data{3}(1,:);
rho_PD_R3_T35 = ds_rho_PD(i).data{3}(2,:);
##
i = id_R4_PD = get_fp_src_idx (ds_rho_PD, "SunT2004");
ds_rho_PD_ext = get_fp_dataset (fp, {"propylene glycol"}, pname, {"SunT2004"}); # pure
T_PD_R4 = mean (ds_rho_PD(i).data{1}');
w_PD_R4 = fp_mf_mx (ds_rho_PD(i).data{2}, mm_PD, mm_W);
rho_PD_R4_T24 = ds_rho_PD(i).data{3}(1,:);
T_PD_R4_ext = ds_rho_PD_ext.data{1}(1);
w_PD_R4_ext = 1;
rho_PD_R4_ext = ds_rho_PD_ext.data{2}(1);
##
i = id_R5_PD = get_fp_src_idx (ds_rho_PD, "MacBethG1951");
T_PD_R5 = ds_rho_PD(i).data{1};
w_PD_R5 = ds_rho_PD(i).data{2};
rho_PD_R5_T35 = ds_rho_PD(i).data{3};
##
fh = figure (); hold on
plot (w_PD_R1, rho_PD_R1_T20, ["x" ";" ds_rho_PD(id_R1_PD).source{1}{1} " T = " num2str(T_PD_R1(1)) " K" ";"])
plot (w_PD_R1, rho_PD_R1_T25, ["kx" ";" ds_rho_PD(id_R1_PD).source{1}{1} " T = " num2str(T_PD_R1(2)) " K" ";"])
plot (w_PD_R1, rho_PD_R1_T30, ["x" ";" ds_rho_PD(id_R1_PD).source{1}{1} " T = " num2str(T_PD_R1(3)) " K" ";"])
plot (w_PD_R2, rho_PD_R2_T25, ["ko" ";" ds_rho_PD(id_R2_PD).source{1}{1} " T = " num2str(T_PD_R2(1)) " K" ";"])
plot (w_PD_R3, rho_PD_R3_T25, ["k*" ";" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(T_PD_R3(1)) " K" ";"])
plot (w_PD_R3, rho_PD_R3_T35, ["*" ";" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(T_PD_R3(2)) " K" ";"])
plot (1, rho_PD_R4_ext, ["k>" ";" ds_rho_PD(id_R4_PD).source{1}{1} " T = " num2str(T_PD_R4_ext) " K" ";"])
plot (w_PD_R4, rho_PD_R4_T24, ["^" ";" ds_rho_PD(id_R4_PD).source{1}{1} " T = " num2str(T_PD_R4(1)) " K" ";"])
plot (w_PD_R5, rho_PD_R5_T35, ["+" ";" ds_rho_PD(id_R5_PD).source{1}{1} " T = " num2str(T_PD_R5) " K" ";"])
plot (w_plot, get_fp_tab (pdir, fname{1}, "rho", 273.15+20, w_plot, []), ["k-.;" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(273.15+20) " K" " extrapolated;"])
plot (w_plot, get_fp_tab (pdir, fname{1}, "rho", 273.15+25, w_plot, []), ["k-;" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(273.15+25) " K" " interpolated;"])
plot (w_plot, get_fp_tab (pdir, fname{1}, "rho", 273.15+30, w_plot, []), ["k--;" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(273.15+30) " K" " interpolated;"])
plot (w_PD_match * [1 1 1 0], rho_PD_match * [0 1 1 1], ["r-.;ri-matching PDMS @ T = " num2str(298.15) " K" ";"])
legend ("location", "southeast")
xlabel ("mass fraction propylene glycol in g / g")
ylabel ("density in kg / m^3")
ylim ([990 1050])
print (fh, "-dpng", "-color", [save_dir "rho_PD"]);

write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_rho_PD(id_R1_PD).source{1}{1}], [w_PD_R1' rho_PD_R1_T20' rho_PD_R1_T25' rho_PD_R1_T30'], {"w PD", "rho in kg / m^3 @ 20 °C", "rho in kg / m^3 @ 25 °C", "rho in kg / m^3 @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_rho_PD(id_R2_PD).source{1}{1}], [w_PD_R2' rho_PD_R2_T25'], {"w PD", "rho in kg / m^3 @ 25 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_rho_PD(id_R3_PD).source{1}{1}], [w_PD_R3' rho_PD_R3_T25' rho_PD_R3_T35'], {"w PD", "rho in kg / m^3 @ 25 °C", "rho in kg / m^3 @ 35 °C"}, [])
##write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_rho_PD(id_R4_PD).source{1}{1}], [w_PD_R4' rho_PD_R4_T24'], {"w PD", "rho in kg / m^3 @ 24 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_rho_PD(id_R5_PD).source{1}{1}], [w_PD_R5' rho_PD_R5_T35'], {"w PD", "rho in kg / m^3 @ 35 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_rho_PD(id_R3_PD).source{1}{1} "_interpolated"], [w_plot' [get_fp_tab(pdir, fname{1}, "rho", 273.15+[20 25 30], w_plot, [])]], {"w PD", "rho in kg / m^3 @ 20 °C", "rho in kg / m^3 @ 25 °C", "rho in kg / m^3 @ 30 °C"}, [])



##
## density aq. glycerol
##
fname = {"glycerol-water"}
pname = {"density"}
ds_rho_PT = get_fp_dataset (fp, fname, pname, [])
## lit. exp.
i = id_R1_PT = get_fp_src_idx (ds_rho_PT, "BosartL1927");
T_PT_R1 = ds_rho_PT(i).data{1};
w_PT_R1 = ds_rho_PT(i).data{2};
rho_PT_R1_T20 = ds_rho_PT(i).data{3}(1,:);
rho_PT_R1_T25 = ds_rho_PT(i).data{3}(2,:);
rho_PT_R1_T30 = ds_rho_PT(i).data{3}(3,:);
## lit. mod.
rho_PT_R2_T20 = get_fp_tab (pdir, fname{1}, "rho", 273.15+20, w_plot, []);
rho_PT_R2_T25 = get_fp_tab (pdir, fname{1}, "rho", 273.15+25, w_plot, []);
rho_PT_R2_T30 = get_fp_tab (pdir, fname{1}, "rho", 273.15+30, w_plot, []);
##
fh = figure (); hold on
plot (w_PT_R1, rho_PT_R1_T20, ["x" ";" ds_rho_PT(id_R1_PT).source{1}{1} " T = " num2str(T_PT_R1(1)) " K" ";"])
plot (w_PT_R1, rho_PT_R1_T25, ["kx" ";" ds_rho_PT(id_R1_PT).source{1}{1} " T = " num2str(T_PT_R1(2)) " K" ";"])
plot (w_PT_R1, rho_PT_R1_T30, ["x" ";" ds_rho_PT(id_R1_PT).source{1}{1} " T = " num2str(T_PT_R1(3)) " K" ";"])
plot (w_plot, rho_PT_R2_T20, ["k-" ";" "VolkA2018" " T = " num2str(273.15+20) " K" ";"])
plot (w_plot, rho_PT_R2_T25, ["k-" ";" "VolkA2018" " T = " num2str(273.15+25) " K" ";"])
plot (w_plot, rho_PT_R2_T30, ["k-" ";" "VolkA2018" " T = " num2str(273.15+30) " K" ";"])
plot (w_PT_match * [1 1 1 0], rho_PT_match * [0 1 1 1], ["r-.;ri-matching PDMS @ T = " num2str(298.15) " K" ";"])
legend ("location", "northwest")
xlabel ("mass fraction glycerol in g / g")
ylabel ("density in kg / m^3")
ylim ([990 1275])
print (fh, "-dpng", "-color", [save_dir "rho_PT"]);

write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_rho_PT(id_R1_PT).source{1}{1}], [w_PT_R1' rho_PT_R1_T20' rho_PT_R1_T25' rho_PT_R1_T30'], {"w PT", "rho in kg / m^3 @ 20 °C", "rho in kg / m^3 @ 25 °C", "rho in kg / m^3 @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" "VolkA2018"], [w_plot' rho_PT_R2_T20' rho_PT_R2_T25' rho_PT_R2_T30'], {"w PT", "rho in kg / m^3 @ 20 °C", "rho in kg / m^3 @ 25 °C", "rho in kg / m^3 @ 30 °C"}, [])



##
## refractive index aq. propylene glycol
##
fname = {"propylene glycol-water"}
pname = {"refractive-index"}
ds_n_PD = get_fp_dataset (fp, fname, pname, [])
## lit
i = id_N1_PD = get_fp_src_idx (ds_n_PD, "MacBethG1951");
w_PD_N1 = ds_n_PD(i).data{1,2};
n_PD_N1 = ds_n_PD(i).data{1,3};
## exp water
ds_n_W = get_fp_dataset (fp, {"water"}, pname, []);
i = get_fp_src_idx (ds_n_PD, "GerkeSJexp");
n_W_T = ds_n_W(i).data{1,2};
## exp propylene glycol
ds_n_PD_pure = get_fp_dataset (fp, {"propylene glycol"}, pname, []);
i = get_fp_src_idx (ds_n_PD_pure, "GerkeSJexp");
n_PD_T = ds_n_PD_pure(i).data{1,2};
## exp aqueous propylene glycol
i = id_N2_PD = get_fp_src_idx (ds_n_PD, "GerkeSJexp");
w_PD_N2 = [0 ds_n_PD(i).data{1,2} 1];
n_PD_N2_T20 = [n_W_T(1) ds_n_PD(i).data{1,3}(1,:) n_PD_T(1)];
n_PD_N2_T25 = [n_W_T(2) ds_n_PD(i).data{1,3}(2,:) n_PD_T(2)];
n_PD_N2_T30 = [n_W_T(3) ds_n_PD(i).data{1,3}(3,:) n_PD_T(3)];
## model fitted to exp
w_plot_PD = [0.5:0.01:0.95];
n_PD_N3_T20 = ri_matching_ri (w_plot_PD, 273.15+20, calib_w.fit_n_PD.c);
n_PD_N3_T25 = ri_matching_ri (w_plot_PD, 273.15+25, calib_w.fit_n_PD.c);
n_PD_N3_T30 = ri_matching_ri (w_plot_PD, 273.15+30, calib_w.fit_n_PD.c);
##
fh = figure (); hold on
plot (w_PD_N1, n_PD_N1, ["kx" ";" ds_n_PD(id_N1_PD).source{1}{1} " T = " num2str(273.15+25) " K" ";"])
plot (w_PD_N2, n_PD_N2_T20, ["^" ";" ds_n_PD(id_N2_PD).source{1}{1} " T = " num2str(273.15+20) " K" ";"])
plot (w_PD_N2, n_PD_N2_T25, ["k^" ";" ds_n_PD(id_N2_PD).source{1}{1} " T = " num2str(273.15+25) " K" ";"])
plot (w_PD_N2, n_PD_N2_T30, ["^" ";" ds_n_PD(id_N2_PD).source{1}{1} " T = " num2str(273.15+30) " K" ";"])
plot (w_plot_PD, n_PD_N3_T20, ["k-" ";" " matching model = " num2str(273.15+20) " K" ";"])
plot (w_plot_PD, n_PD_N3_T25, ["k-" ";" " matching model = " num2str(273.15+25) " K" ";"])
plot (w_plot_PD, n_PD_N3_T30, ["k-" ";" " matching model = " num2str(273.15+25) " K" ";"])
plot (w_PD_match * [1 1 1 0], n_match * [0 1 1 1], ["r-.;ri-matching PDMS @ T = " num2str(298.15) " K" ";"])
legend ("location", "northwest")
xlabel ("mass fraction glycerol in g / g")
ylabel ("refractive index in -")
ylim ([1.32 1.44])
print (fh, "-dpng", "-color", [save_dir "n_PD"]);

write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_n_PD(id_N1_PD).source{1}{1}], [w_PD_N1' n_PD_N1'], {"w PD", "n in - @ 25 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_n_PD(id_N2_PD).source{1}{1}], [w_PD_N2' n_PD_N2_T20' n_PD_N2_T25' n_PD_N2_T30'], {"w PD", "n in - @ 20 °C", "n in - @ 25 °C", "n in - @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_n_PD(id_N2_PD).source{1}{1} "matching model"], [w_plot_PD' n_PD_N3_T20' n_PD_N3_T25' n_PD_N3_T30'], {"w PD", "n in - @ 20 °C", "n in - @ 25 °C", "n in - @ 30 °C"}, [])



##
## refractive index aq. glycerol
##
fname = {"glycerol-water"}
pname = {"refractive-index"}
ds_n_PT = get_fp_dataset (fp, fname, pname, [])
## lit
i = id_N1_PT = get_fp_src_idx (ds_n_PT, "GP1963etal");
w_PT_N1 = ds_n_PT(i).data{1,2};
n_PT_N1 = ds_n_PT(i).data{1,3};
## exp propylene glycol
ds_n_PT_pure = get_fp_dataset (fp, {"glycerol"}, pname, []);
i = get_fp_src_idx (ds_n_PT_pure, "GerkeSJexp");
n_PT_T = ds_n_PT_pure(i).data{1,2};
## exp
i = id_N2_PT = get_fp_src_idx (ds_n_PT, "GerkeSJexp");
w_PT_N2 = [0 ds_n_PT(i).data{1,2} 1];
n_PT_N2_T20 = [n_W_T(1) ds_n_PT(i).data{1,3}(1,:) n_PT_T(1)];
n_PT_N2_T25 = [n_W_T(2) ds_n_PT(i).data{1,3}(2,:) n_PT_T(2)];
n_PT_N2_T30 = [n_W_T(3) ds_n_PT(i).data{1,3}(3,:) n_PT_T(3)];
## model fitted to exp
w_plot_PT = [0.4:0.01:0.75];
n_PT_N3_T20 = ri_matching_ri (w_plot_PT, 273.15+20, calib_w.fit_n_PT.c);
n_PT_N3_T25 = ri_matching_ri (w_plot_PT, 273.15+25, calib_w.fit_n_PT.c);
n_PT_N3_T30 = ri_matching_ri (w_plot_PT, 273.15+30, calib_w.fit_n_PT.c);
##
fh = figure (); hold on
plot (w_PT_N1, n_PT_N1, ["kx" ";" ds_n_PT(id_N1_PT).source{1}{1} " T = " num2str(273.15+20) " K" ";"])
plot (w_PT_N2, n_PT_N2_T20, ["^" ";" ds_n_PT(id_N2_PT).source{1}{1} " T = " num2str(273.15+20) " K" ";"])
plot (w_PT_N2, n_PT_N2_T25, ["k^" ";" ds_n_PT(id_N2_PT).source{1}{1} " T = " num2str(273.15+25) " K" ";"])
plot (w_PT_N2, n_PT_N2_T30, ["^" ";" ds_n_PT(id_N2_PT).source{1}{1} " T = " num2str(273.15+30) " K" ";"])
plot (w_plot_PT, n_PT_N3_T20, ["k-" ";" " matching model = " num2str(273.15+20) " K" ";"])
plot (w_plot_PT, n_PT_N3_T25, ["k-" ";" " matching model = " num2str(273.15+25) " K" ";"])
plot (w_plot_PT, n_PT_N3_T30, ["k-" ";" " matching model = " num2str(273.15+30) " K" ";"])
plot (w_PT_match * [1 1 1 0], n_match * [0 1 1 1], ["r-.;ri-matching PDMS @ T = " num2str(298.15) " K" ";"])
legend ("location", "northwest")
xlabel ("mass fraction glycerol in g / g")
ylabel ("refractive index in -")
ylim ([1.32 1.48])
print (fh, "-dpng", "-color", [save_dir "n_PT"]);

write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_n_PT(id_N1_PT).source{1}{1}], [w_PT_N1' n_PT_N1'], {"w PT", "n in - @ 20 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_n_PT(id_N2_PT).source{1}{1}], [w_PT_N2' n_PT_N2_T20' n_PT_N2_T25' n_PT_N2_T30'], {"w PT", "n in - @ 20 °C", "n in - @ 25 °C", "n in - @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_n_PT(id_N2_PT).source{1}{1} "matching model"], [w_plot_PT' n_PT_N3_T20' n_PT_N3_T25' n_PT_N3_T30'], {"w PT", "n in - @ 20 °C", "n in - @ 25 °C", "n in - @ 30 °C"}, [])



##
## dynamic viscosity aq. glycerol
##
fname = {"glycerol-water"}
pname = {"dynamic-viscosity"}
ds_eta_PT = get_fp_dataset (fp, fname, pname, [])
## lit
i = id_E1_PT = get_fp_src_idx (ds_eta_PT, "SegurJ1951");
T_PT_E1 = ds_eta_PT(i).data{1,1};
w_PT_E1 = ds_eta_PT(i).data{1,2};
eta_PT_E1_T20 = ds_eta_PT(i).data{1,3}(2,:);
eta_PT_E1_T30 = ds_eta_PT(i).data{1,3}(3,:);
## exp
i = id_E2_PT = get_fp_src_idx (ds_eta_PT, "GerkeSJexp");
T_PT_E2 = ds_eta_PT(i).data{1,1};
w_PT_E2 = ds_eta_PT(i).data{1,2};
eta_PT_E2_T20 = ds_eta_PT(i).data{1,3}(1,:);
eta_PT_E2_T25 = ds_eta_PT(i).data{1,3}(2,:);
eta_PT_E2_T30 = ds_eta_PT(i).data{1,3}(3,:);
## lit model of VolkA2018
w_PT_E3_T20 = eta_PT_W_model (w_plot, 273.15+20)';
w_PT_E3_T25 = eta_PT_W_model (w_plot, 273.15+25)';
w_PT_E3_T30 = eta_PT_W_model (w_plot, 273.15+30)';
##
fh = figure (); hold on
plot (w_PT_E1, eta_PT_E1_T20, ["x" ";" ds_eta_PT(id_E1_PT).source{1}{1} " T = " num2str(T_PT_E1(2)) " K" ";"])
plot (w_PT_E1, eta_PT_E1_T30, ["x" ";" ds_eta_PT(id_E1_PT).source{1}{1} " T = " num2str(T_PT_E1(3)) " K" ";"])
plot (w_PT_E2, eta_PT_E2_T20, ["^" ";" ds_eta_PT(id_E2_PT).source{1}{1} " T = " num2str(T_PT_E2(1)) " K" ";"])
plot (w_PT_E2, eta_PT_E2_T25, ["k^" ";" ds_eta_PT(id_E2_PT).source{1}{1} " T = " num2str(T_PT_E2(2)) " K" ";"])
plot (w_PT_E2, eta_PT_E2_T30, ["^" ";" ds_eta_PT(id_E2_PT).source{1}{1} " T = " num2str(T_PT_E2(3)) " K" ";"])
##plot (w_plot, get_fp_tab (pdir, fname{1}, "eta", 273.15+20, w_plot, []), ["b--" ";" ds_eta_PT(id_E1_PT).source{1}{1} " T = " num2str(273.15+20) " K interpolated" ";"])
plot (w_plot, get_fp_tab (pdir, fname{1}, "eta", 273.15+25, w_plot, []), ["b--" ";" ds_eta_PT(id_E1_PT).source{1}{1} " T = " num2str(273.15+25) " K interpolated" ";"])
##plot (w_plot, get_fp_tab (pdir, fname{1}, "eta", 273.15+30, w_plot, []), ["b--" ";" ds_eta_PT(id_E1_PT).source{1}{1} " T = " num2str(273.15+30) " K interpolated" ";"])
plot (w_plot, w_PT_E3_T20, ["k-." ";" "model VolkA2018 T = " num2str(273.15+20) " K" ";"])
plot (w_plot, w_PT_E3_T25, ["k-" ";" "model VolkA2018 T = " num2str(273.15+25) " K" ";"])
plot (w_plot, w_PT_E3_T30, ["k--" ";" "model VolkA2018 T = " num2str(273.15+30) " K" ";"])
plot (w_PT_match * [1 1 1 0], eta_PT_match * [0 1 1 1], ["r-.;ri-matching PDMS @ T = " num2str(298.15) " K" ";"])
legend ("location", "northwest")
xlabel ("mass fraction glycerol in g / g")
ylabel ("dynamic viscosity in Pa s")
print (fh, "-dpng", "-color", [save_dir "eta_PT"]);

write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PT(id_E1_PT).source{1}{1}], [w_PT_E1' eta_PT_E1_T20' eta_PT_E1_T30'], {"w PT", "eta in Pa s @ 20 °C", "eta in Pa s @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PT(id_E2_PT).source{1}{1}], [w_PT_E2' eta_PT_E2_T20' eta_PT_E2_T25' eta_PT_E2_T30'], {"w PT", "eta in Pa s @ 20 °C", "eta in Pa s @ 25 °C", "eta in Pa s @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" "VolkA2018"], [w_plot' w_PT_E3_T20' w_PT_E3_T25' w_PT_E3_T30'], {"w PT", "eta in Pa s @ 20 °C", "eta in Pa s @ 25 °C", "eta in Pa s @ 30 °C"}, [])



##
## dynamic viscosity aq. propylene glycol
##
fname = {"propylene glycol-water"}
pname = {"dynamic-viscosity"}
ds_eta_PD = get_fp_dataset (fp, fname, pname, [])
## lit
i = id_E1_PD = get_fp_src_idx (ds_eta_PD, "GeorgeJ2003");
T_PD_E1 = ds_eta_PD(i).data{1,1};
w_PD_E1 = fp_mf_mx (ds_eta_PD(i).data{1,2}, mm_PD, mm_W);
eta_PD_E1_T25 = ds_eta_PD(i).data{1,3}(1,:);
eta_PD_E1_T35 = ds_eta_PD(i).data{1,3}(2,:);
##
i = id_E2_PD = get_fp_src_idx (ds_eta_PD, "TanakaY1988");
T_PD_E2 = ds_eta_PD(i).data{1,1};
w_PD_E2 = ds_eta_PD(i).data{1,2};
eta_PD_E2_T25 = ds_eta_PD(i).data{1,3}(1,:);
##
i = id_E3_PD = get_fp_src_idx (ds_eta_PD, "SunT2004");
T_PD_E3 =  mean (ds_eta_PD(i).data{1}');
w_PD_E3 = fp_mf_mx (ds_eta_PD(i).data{1,2}, mm_PD, mm_W);
eta_PD_E3_T24 = ds_eta_PD(i).data{1,3}(1,:);
ds_eta_PD_ext = get_fp_dataset (fp, {"propylene glycol"}, pname, {"SunT2004"}); # pure
T_PD_E3_ext = ds_eta_PD_ext.data{1}(1);
w_PD_E3_ext = 1;
eta_PD_E3_ext = ds_eta_PD_ext.data{2}(1);
##
i = id_E4_PD = get_fp_src_idx (ds_eta_PD, "KhattabIS2017");
T_PD_E4 = ds_eta_PD(i).data{1,1};
w_PD_E4 = fp_mf_mx (ds_eta_PD(i).data{1,2}, mm_PD, mm_W);
eta_PD_E4_T20 = ds_eta_PD(i).data{1,3}(1,:);
eta_PD_E4_T25 = ds_eta_PD(i).data{1,3}(2,:);
eta_PD_E4_T30 = ds_eta_PD(i).data{1,3}(3,:);
## exp
i = id_E5_PD = get_fp_src_idx (ds_eta_PD, "GerkeSJexp");
T_PD_E5 = ds_eta_PD(i).data{1,1};
w_PD_E5 = ds_eta_PD(i).data{1,2};
eta_PD_E5_T20 = ds_eta_PD(i).data{1,3}(1,:);
eta_PD_E5_T25 = ds_eta_PD(i).data{1,3}(2,:);
eta_PD_E5_T30 = ds_eta_PD(i).data{1,3}(3,:);
##
fh = figure (); hold on
plot (w_PD_E1, eta_PD_E1_T25, ["kx" ";" ds_eta_PD(id_E1_PD).source{1}{1} " T = " num2str(T_PD_E1(1)) " K" ";"])
plot (w_PD_E1, eta_PD_E1_T35, ["x" ";" ds_eta_PD(id_E1_PD).source{1}{1} " T = " num2str(T_PD_E1(2)) " K" ";"])
plot (w_PD_E2, eta_PD_E2_T25, ["kd" ";" ds_eta_PD(id_E2_PD).source{1}{1} " T = " num2str(T_PD_E2(1)) " K" ";"])
plot (w_PD_E3, eta_PD_E3_T24, ["^" ";" ds_eta_PD(id_E3_PD).source{1}{1} " T = " num2str(T_PD_E3(1)) " K" ";"])
plot (w_PD_E3_ext, eta_PD_E3_ext, ["^" ";" ds_eta_PD(id_E3_PD).source{1}{1} " T = " num2str(T_PD_E3_ext(1)) " K" ";"])
plot (w_PD_E4, eta_PD_E4_T20, ["o" ";" ds_eta_PD(id_E4_PD).source{1}{1} " T = " num2str(T_PD_E4(1)) " K" ";"])
plot (w_PD_E4, eta_PD_E4_T25, ["ko" ";" ds_eta_PD(id_E4_PD).source{1}{1} " T = " num2str(T_PD_E4(2)) " K" ";"])
plot (w_PD_E4, eta_PD_E4_T30, ["o" ";" ds_eta_PD(id_E4_PD).source{1}{1} " T = " num2str(T_PD_E4(3)) " K" ";"])
plot (w_PD_E5, eta_PD_E5_T20, ["*" ";" ds_eta_PD(id_E5_PD).source{1}{1} " T = " num2str(T_PD_E5(1)) " K" ";"])
plot (w_PD_E5, eta_PD_E5_T25, ["k*" ";" ds_eta_PD(id_E5_PD).source{1}{1} " T = " num2str(T_PD_E5(2)) " K" ";"])
plot (w_PD_E5, eta_PD_E5_T30, ["*" ";" ds_eta_PD(id_E5_PD).source{1}{1} " T = " num2str(T_PD_E5(3)) " K" ";"])
plot (w_plot, get_fp_tab (pdir, fname{1}, "eta", 273.15+20, w_plot, []), ["k-." ";" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(273.15+20) " K interpolated" ";"])
plot (w_plot, get_fp_tab (pdir, fname{1}, "eta", 273.15+25, w_plot, []), ["k-" ";" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(273.15+25) " K interpolated" ";"])
plot (w_plot, get_fp_tab (pdir, fname{1}, "eta", 273.15+30, w_plot, []), ["k--" ";" ds_rho_PD(id_R3_PD).source{1}{1} " T = " num2str(273.15+30) " K interpolated" ";"])
plot (w_PD_match * [1 1 1 0], eta_PD_match * [0 1 1 1], ["r-.;ri-matching PDMS @ T = " num2str(298.15) " K" ";"])
legend ("location", "northwest")
xlabel ("mass fraction propylene glycol in g / g")
ylabel ("dynamic viscosity in Pa s")
print (fh, "-dpng", "-color", [save_dir "eta_PD"]);

write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PD(id_E1_PD).source{1}{1}], [w_PD_E1' eta_PD_E1_T25' eta_PD_E1_T35'], {"w PD", "eta in Pa s @ 25 °C", "eta in Pa s @ 35 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PD(id_E2_PD).source{1}{1}], [w_PD_E2' eta_PD_E2_T25'], {"w PD", "eta in Pa s @ 25 °C"}, [])
##write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PD(id_E3_PD).source{1}{1}], [w_PD_E3' eta_PD_E3_T24'], {"w PD", "eta in Pa s @ 24 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PD(id_E4_PD).source{1}{1}], [w_PD_E4' eta_PD_E4_T20' eta_PD_E4_T25' eta_PD_E4_T30'], {"w PD", "eta in Pa s @ 20 °C", "eta in Pa s @ 25 °C", "eta in Pa s @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PD(id_E5_PD).source{1}{1}], [w_PD_E5' eta_PD_E5_T20' eta_PD_E5_T25' eta_PD_E5_T30'], {"w PD", "eta in Pa s @ 20 °C", "eta in Pa s @ 25 °C", "eta in Pa s @ 30 °C"}, [])
write_series_csv ([save_dir pname{1} "_" fname{1} "_" ds_eta_PD(id_E1_PD).source{1}{1} "_interpolated"], [w_plot' [get_fp_tab(pdir, fname{1}, "eta", 273.15+[20 25 30], w_plot, [])]], {"w PD", "eta in Pa s @ 20 °C", "eta in Pa s @ 25 °C", "eta in Pa s @ 30 °C"}, [])




## table comparing fluid properties for 20 25 and 30 °C for both refractive index matching liquid mixtures
T = 273.15 + [20 25 30]'
[w_L1_T, n_L1_T, rho_L1_T, eta_L1_T, c_sat_L1_T, D_AB_L1_T] = get_fp_lm (pdir, "WG141", T)
[w_L2_T, n_L2_T, rho_L2_T, eta_L2_T, c_sat_L2_T, D_AB_L2_T] = get_fp_lm (pdir, "WP141", T)
write_series_csv ([save_dir "WG141_fluidprop"], [T n_L1_T w_L1_T rho_L1_T eta_L1_T], {"T in K", "n in -", "w in g/g", "rho in kg/m^3", "eta in Pa s"}, [])
write_series_csv ([save_dir "WP141_fluidprop"], [T n_L2_T w_L2_T rho_L2_T eta_L2_T], {"T in K", "n in -", "w in g/g", "rho in kg/m^3", "eta in Pa s"}, [])
write_series_csv ([save_dir "WG141_WP141_fluidprop"], [[T n_L1_T w_L1_T rho_L1_T eta_L1_T];[T n_L2_T w_L2_T rho_L2_T eta_L2_T]], {"T in K", "n in -", "w in g/g", "rho in kg/m^3", "eta in Pa s"}, [])


## TODO: sensitivity of fluid properties to mass fraction estimation by refractive index measurement



