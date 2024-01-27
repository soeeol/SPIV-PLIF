##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## trying out several plots and output data to csv
## to be run from a_diffusivity_solubility.m
##
## Author: Sören J. Gerke
##

figp = fig_param ({"lineart", "double", "elsevier"});
fig_size_x = 20;
fig_size_y = 15;

##
## Water - Glycerol
##
fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
title("Water-Glycerol O_2 solubility from air @ 1 atm")
hold on;
grid on;
xlabel("mass fraction glycerol")
ylabel("saturation concentration O_2 in mg/l")
plot(w_PT_SJG, c_mg_PT_SJG_min, "k-.d;exp. @25°C UMS min. stirred GerkeSJxexp;");
plot(w_PT_NK, csat_PT_NK,"-.x;exp. @22.5°C NogamiH1962et;");
plot(w_Jetal, c_mg_Jetal,"-.o;exp. @25°C JordanJ1956etal;");
##plot(w_Jetal, 2.1/1.6 *c_mg_Jetal,"-.o;JordanJ1956etal;");
##plot([0], csat_oxygen_water_T_model (T=TK(3), p_O2), "x; 25°C BensonBB1979etal;")
plot(w_Retal, c_mg_Retal,"-.s;exp. @30°C RischbieterE1996etal;");
plot([0 w_Ketal], c_mg_Ketal,"-.*;exp. @25°C KutscheI1984etal;");
##
plot(fp_mf_mc (rho_W(1), rho_PT(1), [0:10:rho_PT(1)]), f_c_mg_Retal(csat_W(1), Kn_PT(1), [0:10:rho_PT(1)]),"k-;corr. @20°C RischbieterE1996etal;", "linewidth", 1);
##plot(fp_mf_mc (rho_W(2),  rho_PT(2), [0:10: rho_PT(2)]), f_c_mg_Retal(csat_W, Kn_PT(2), [0:10: rho_PT(2)]),"k-;corr. @22.5°C;", "linewidth", 1);
plot(fp_mf_mc (rho_W(3), rho_PT(3), [0:10:rho_PT(3)]), f_c_mg_Retal(csat_W(3), Kn_PT(3), [0:10:rho_PT(3)]),"k-;corr. @25°C RischbieterE1996etal;", "linewidth", 1);
plot(fp_mf_mc (rho_W(4), rho_PT(4), [0:10:rho_PT(4)]), f_c_mg_Retal(csat_W(4), Kn_PT(4), [0:10:rho_PT(4)]),"k-;corr. @30°C RischbieterE1996etal;", "linewidth", 1);
legend ("location", "eastoutside")
##
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir "fig_csat_water-glycerol_mf"]);
## output
plt_1_h = {mfilename, date; "w_PT", "c_sat in mg/l; exp. @25°C UMS min. stirred GerkeSJxexp"};
plt_1_d = [w_PT_SJG' c_mg_PT_SJG_min'];
cell2csv (["csat_water-glycerol_SJG_h.csv"], plt_1_h)
csvwrite (["csat_water-glycerol_SJG_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
plt_2_h = {mfilename, date; "w_PT", "c_sat in mg/l; exp. @22.5°C NogamiH1962et"};
plt_2_d = [w_PT_NK' csat_PT_NK'];
cell2csv (["csat_water-glycerol_NogamiH1962et_h.csv"], plt_2_h)
csvwrite (["csat_water-glycerol_NogamiH1962et_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
plt_3_h = {mfilename, date; "w_PT", "c_sat in mg/l; exp. @25°C JordanJ1956etal"};
plt_3_d = [w_Jetal' c_mg_Jetal'];
cell2csv (["csat_water-glycerol_JordanJ1956etal_h.csv"], plt_3_h)
csvwrite (["csat_water-glycerol_JordanJ1956etal_d.csv"], plt_3_d, "append", "off", "precision","%.4e")
plt_4_h = {mfilename, date; "w_PT", "c_sat in mg/l; exp. @30°C RischbieterE1996etal"};
plt_4_d = [w_Retal' c_mg_Retal'];
cell2csv (["csat_water-glycerol_RischbieterE1996etal_h.csv"], plt_4_h)
csvwrite (["csat_water-glycerol_RischbieterE1996etal_d.csv"], plt_4_d, "append", "off", "precision","%.4e")
plt_5_h = {mfilename, date; "w_PT", "c_sat in mg/l; exp. @25°C KutscheI1984etal"};
plt_5_d = [[0 w_Ketal]' c_mg_Ketal'];
cell2csv (["csat_water-glycerol_KutscheI1984etal_h.csv"], plt_5_h)
csvwrite (["csat_water-glycerol_KutscheI1984etal_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_6_h = {mfilename, date; "w_PT", "c_sat in mg/l;corr. @20°C RischbieterE1996etal"};
plt_6_d = [[fp_mf_mc(rho_W(1),rho_PT(1),[0:10:rho_PT(1)])]' f_c_mg_Retal(csat_W(1),Kn_PT(1),[0:10:rho_PT(1)])'];
cell2csv (["csat_water-glycerol_RischbieterE1996etal_T20_h.csv"], plt_6_h)
csvwrite (["csat_water-glycerol_RischbieterE1996etal_T20_d.csv"], plt_6_d, "append", "off", "precision","%.4e")
plt_7_h = {mfilename, date; "w_PT", "c_sat in mg/l;corr. @25°C RischbieterE1996etal"};
plt_7_d = [[fp_mf_mc(rho_W(3),rho_PT(3),[0:10:rho_PT(3)])]' f_c_mg_Retal(csat_W(3),Kn_PT(3),[0:10:rho_PT(3)])'];
cell2csv (["csat_water-glycerol_RischbieterE1996etal_T25_h.csv"], plt_7_h)
csvwrite (["csat_water-glycerol_RischbieterE1996etal_T25_d.csv"], plt_7_d, "append", "off", "precision","%.4e")
plt_8_h = {mfilename, date; "w_PT", "c_sat in mg/l;corr. @30°C RischbieterE1996etal"};
plt_8_d = [[fp_mf_mc(rho_W(4),rho_PT(4),[0:10:rho_PT(4)])]' f_c_mg_Retal(csat_W(4),Kn_PT(4),[0:10:rho_PT(4)])'];
cell2csv (["csat_water-glycerol_RischbieterE1996etal_T30_h.csv"], plt_8_h)
csvwrite (["csat_water-glycerol_RischbieterE1996etal_T30_d.csv"], plt_8_d, "append", "off", "precision","%.4e")

##
## Water - Propylene Glycol
##
fh = figure(); hold on
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
title("Water - Propylene Glycol \n O_2 solubility from air @ 1 atm")
hold on;
grid on;
xlabel("mass fraction Propylene Glycol")
ylabel("saturation concentration O_2 in mg/l")
plot (w_PD_NK, csat_PD_NK,"-.x;exp. @22.5°C NogamiH1962et;");
plot (w_PD_Y, c_mg_Y,"-.x;exp. @25°C YamamotoH1994etal;");
plot([0], csat_oxygen_water_T_model (T=TK(1), p_O2), "s;Water @20°C BensonBB1979etal;")
plot([0], csat_oxygen_water_T_model (T=TK(3), p_O2), "s;Water @25°C BensonBB1979etal;")
plot ([1], csat_PD_O, "k*;exp. @25°C OsborneAD1965et;")
plot (w_PD_SJG, c_mg_PD_SJG_min, "-.ko;exp. @25°C UMS min. stirred GerkeSJxexp;")
plot(fp_mf_mc (rho_W(3), rho_PD(3), [0:10:rho_PD(3)]), f_c_mg_Retal(csat_W(3), Kn_PD(3), [0:10:rho_PD(3)]),"k-;corr. @25°C RischbieterE1996etal;", "linewidth", 1);
legend ("location", "eastoutside")
##
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir "fig_csat_water-propylene glycol_mf"]);
## output
plt_1_h = {mfilename, date; "w_PD", "c_sat in mg/l; exp. @25°C UMS min. stirred GerkeSJxexp"};
plt_1_d = [w_PD_SJG' c_mg_PD_SJG_min'];
cell2csv (["csat_water-propylene glycol_SJG_h.csv"], plt_1_h)
csvwrite (["csat_water-propylene glycol_SJG_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
plt_2_h = {mfilename, date; "w_PD", "c_sat in mg/l; exp. @22.5°C NogamiH1962et"};
plt_2_d = [w_PD_NK' csat_PD_NK'];
cell2csv (["csat_water-propylene glycol_NogamiH1962et_h.csv"], plt_2_h)
csvwrite (["csat_water-propylene glycol_NogamiH1962et_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
plt_3_h = {mfilename, date; "w_PD", "c_sat in mg/l; exp. @25°C YamamotoH1994etal"};
plt_3_d = [w_PD_Y' c_mg_Y'];
cell2csv (["csat_water-propylene glycol_YamamotoH1994etal_h.csv"], plt_3_h)
csvwrite (["csat_water-propylene glycol_YamamotoH1994etal_d.csv"], plt_3_d, "append", "off", "precision","%.4e")
plt_4_h = {mfilename, date; "w_PD", "c_sat in mg/l; exp. @25°C OsborneAD1965et"};
plt_4_d = [1' csat_PD_O'];
cell2csv (["csat_water-propylene glycol_OsborneAD1965et_h.csv"], plt_4_h)
csvwrite (["csat_water-propylene glycol_OsborneAD1965et_d.csv"], plt_4_d, "append", "off", "precision","%.4e")
plt_5_h = {mfilename, date; "w_PD", "c_sat in mg/l; water @20°C BensonBB1979etal"};
plt_5_d = [0' [csat_oxygen_water_T_model(T=TK(1),p_O2)]'];
cell2csv (["csat_water_BensonBB1979etal_T20_h.csv"], plt_5_h)
csvwrite (["csat_water_BensonBB1979etal_T20_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_6_h = {mfilename, date; "w_PD", "c_sat in mg/l; water @25°C BensonBB1979etal"};
plt_6_d = [0' [csat_oxygen_water_T_model(T=TK(3),p_O2)]'];
cell2csv (["csat_water_BensonBB1979etal_T25_h.csv"], plt_6_h)
csvwrite (["csat_water_BensonBB1979etal_T25_d.csv"], plt_6_d, "append", "off", "precision","%.4e")
plt_7_h = {mfilename, date; "w_PD", "c_sat in mg/l;Glycol corr. @25°C RischbieterE1996etal"};
plt_7_d = [[fp_mf_mc(rho_W(3),rho_PD(3),[0:10:rho_PD(3)])]' f_c_mg_Retal(csat_W(3),Kn_PD(3),[0:10:rho_PD(3)])'];
cell2csv (["csat_water-propylene glycol_RischbieterE1996etal_T25_h.csv"], plt_7_h)
csvwrite (["csat_water-propylene glycol_RischbieterE1996etal_T25_d.csv"], plt_7_d, "append", "off", "precision","%.4e")


##
## Water
##
fh = figure(); hold on
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
##plot (Ttab, cOOtab, "o;tab;")
##plot ([min(T_Tan):1:max(T_Tan)], p_O2 * 1e5 * mm_O2 * 1e3 * f_Hoxy_T (Href, Coxy, [min(T_Tan):1:max(T_Tan)]),"-;Sander 2015;")
##plot(T_Tan, c_sat_Tan, "s;Tan;")
plot(T_iupac, c_sat_iupac, "d;MiyamotoH2014etal;")
plot(T_fit, c_sat_iupac_fit, "-;fit MiyamotoH2014etal;")
##plot([min(T_Tan):1:max(T_Tan)], csat_oxygen_water_T_model([min(T_Tan):1:max(T_Tan)], p_O2), "k-;Benson et al. 1979;")
plot(T_fit, c_sat_Benson_fit, "k-;fit BensonBB1979etal;")
xlabel ("T in K")
ylabel ("c_s_a_t in mg/l")
title ("oxygen solubility in water for air (20% O_2, 80% N_2) @ 25°C")
##
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir "fig_csat_water_T"]);
## output
plt_1_h = {mfilename, date; "T_iupac", "(MiyamotoH2014etal, smoothed exp. of RettichTR2000etal) oxygen solubility in water for air (20 % O2) @ 25°C"};
plt_1_d = [T_iupac' c_sat_iupac'];
cell2csv (["csat_water_IUPAC_h.csv"], plt_1_h)
csvwrite (["csat_water_IUPAC_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
plt_2_h = {mfilename, date; "T_iupac", "(BensonBB1979etal, fit model"};
plt_2_d = [T_fit' c_sat_Benson_fit'];
cell2csv (["csat_water_BensonBB1979etal_h.csv"], plt_2_h)
csvwrite (["csat_water_BensonBB1979etal_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
plt_3_h = {mfilename, date; "T_iupac", "MiyamotoH2014etal: fit of oxygen solubility in water for air (20 % O2) @ 25°C"};
plt_3_d = [T_fit' c_sat_iupac_fit'];
cell2csv (["csat_water_IUPAC_fit_h.csv"], plt_3_h)
csvwrite (["csat_water_IUPAC_fit_d.csv"], plt_3_d, "append", "off", "precision","%.4e")


