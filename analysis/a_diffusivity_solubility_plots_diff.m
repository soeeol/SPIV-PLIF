##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## trying out several plots and output data to csv
## to be run from a_diffusivity_solubility.m
##
## Author: Sören J. Gerke
##

figp = fig_param ({"lineart", "double", "elsevier"});
fig_size_x = 20;
fig_size_y = 10;

##
## Water - Glycerol
##

## output
plt_0_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. water @25°C StDenisCE1971et"};
plt_0_d = [0 D_W_exp_mean];
cell2csv (["D_water-StDenisCE1971et_h.csv"], plt_0_h)
csvwrite (["D_water-StDenisCE1971et_d.csv"], plt_0_d, "append", "off", "precision","%.4e")
plt_8_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. water diffusivity in glycerol @25°C ErricoG2004etal"};
plt_8_d = [1 D_PT_DErrico];
cell2csv (["D_water-in-glycerol_ErricoG2004etal_h.csv"], plt_8_h)
csvwrite (["D_water-in-glycerol_ErricoG2004etal_d.csv"], plt_8_d, "append", "off", "precision","%.4e")
##
plt_1_h = {mfilename, date; "w_PT", "D_AB in m^2/s; model Stokes-Einstein eq. (4*Pi)"};
plt_1_d = [w' D_PT_SE];
cell2csv (["D_water-PT_Stokes-Einstein_h.csv"], plt_1_h)
csvwrite (["D_water-PT_Stokes-Einstein_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
plt_2_h = {mfilename, date; "w_PT", "D_AB in m^2/s; corr. Wilke-Chang-mix"};
plt_2_d = [w' D_Am_PT_WC];
cell2csv (["D_water-PT_PerkinsLR1969et_h.csv"], plt_2_h)
csvwrite (["D_water-PT_PerkinsLR1969et_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
plt_3_h = {mfilename, date; "w_PT", "D_AB in m^2/s; corr. SitaramanR1963etal with mixing rule eta^0.8"};
plt_3_d = [w' Dim_S_PT'];
cell2csv (["D_water-PT_SitaramanR1963etal_h.csv"], plt_3_h)
csvwrite (["D_water-PT_SitaramanR1963etal_d.csv"], plt_3_d, "append", "off", "precision","%.4e")
plt_4_h = {mfilename, date; "w_PT", "D_AB in m^2/s; corr. DiazM1987etal with mixing rule eta^0.8"};
plt_4_d = [w' Dim_DM_PT'];
cell2csv (["D_water-PT_DiazM1987etal_h.csv"], plt_4_h)
csvwrite (["D_water-PT_DiazM1987etal_d.csv"], plt_4_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. NogamiH1962et"};
plt_5_d = [w_PT_NK' D_PT_NK'];
cell2csv (["D_water-PT_NogamiH1962et_h.csv"], plt_5_h)
csvwrite (["D_water-PT_NogamiH1962et_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_5_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. NogamiH1962et, corr to 25°C"};
plt_5_d = [w_PT_NK' D_PT_NK_T25' D_PT_NK_T25'./Dcorr_NK_PT'];
cell2csv (["D_water-PT_NogamiH1962et_T25_h.csv"], plt_5_h)
csvwrite (["D_water-PT_NogamiH1962et_T25_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. this study, first trend"};
plt_5_d = [w_SJG_PT' D_SJG_PT_1'];
cell2csv (["D_water-PT_SJG_1_h.csv"], plt_5_h)
csvwrite (["D_water-PT_SJG_1_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_6_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. this study, second trend"};
plt_6_d = [w_SJG_PT' D_SJG_PT_2'];
cell2csv (["D_water-PT_SJG_2_h.csv"], plt_6_h)
csvwrite (["D_water-PT_SJG_2_d.csv"], plt_6_d, "append", "off", "precision","%.4e")
plt_7_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. JordanJ1956etal"};
plt_7_d = [w_Jetal' D_Jetal' D_Jetal'./Dcorr_Jetal'];
cell2csv (["D_water-PT_JordanJ1956etal_h.csv"], plt_7_h)
csvwrite (["D_water-PT_JordanJ1956etal_d.csv"], plt_7_d, "append", "off", "precision","%.4e")
plt_9_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. PLIF JimenezM2012etal-1; corr to 25°C"};
plt_9_d = [w_PT_J' D_J_3_T25'];
cell2csv (["D_water-PT-E_JimenezM2012etal-1_h.csv"], plt_9_h)
csvwrite (["D_water-PT-E_JimenezM2012etal-1_d.csv"], plt_9_d, "append", "off", "precision","%.4e")
plt_10_h = {mfilename, date; "w_PT", "D_AB in m^2/s; exp. PLIF KapoustinaV2019etal; corr to 25°C"};
plt_10_d = [w_PT_K' D_K_T25'];
cell2csv (["D_water-PT-E_KapoustinaV2019etal_h.csv"], plt_10_h)
csvwrite (["D_water-PT-E_KapoustinaV2019etal_d.csv"], plt_10_d, "append", "off", "precision","%.4e")

##
fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (w, D_PT_SE, "k.-;Stokes-Einstein eq. (4 Pi);");
semilogy (w, D_Am_PT_WC, "k-;Wike-Chang eq. mod. of Perkins and Geankoplis (1969);");
semilogy (w, Dim_DM_PT, "k--;corr. w. MM (M. Diaz 1987);");
semilogy (w, Dim_S_PT, "k-.;corr. w. MM (Sitaraman et al. 1963);");
semilogy (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
semilogy (w_PT_NK, D_PT_NK_T25,"ks;exp. (Nogami and Kato 1961);");
semilogy (w_Jetal, D_Jetal, "ko;exp. D_F_i_c_k (Jordan et al. 1956);");
semilogy (w_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
semilogy (w_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction Glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T corrected")
grid on
box on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [pdir.plot "diffusivity/" "overview_PT_T25"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (w, D_PT_SE, "k.-;Stokes-Einstein eq. (4 Pi);");
semilogy (w, D_Am_PT_WC, "k-;Wike-Chang eq. mod. of Perkins and Geankoplis (1969);");
semilogy (w, Dim_DM_PT, "k--;corr. w. MM (M. Diaz 1987);");
semilogy (w, Dim_S_PT, "k-.;corr. w. MM (Sitaraman et al. 1963);");
semilogy (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
semilogy (w_PT_NK, D_PT_NK_T25./Dcorr_NK_PT,"ks;exp. (Nogami and Kato 1961);");
semilogy (w_Jetal, D_Jetal./Dcorr_Jetal, "ko;exp. D_F_i_c_k (Jordan et al. 1956);");
semilogy (w_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
semilogy (w_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction Glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T and C corrected")
grid on
box on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [pdir.plot "diffusivity/" "overview_PT_T25_Ccorr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (w, D_PT_SE*D_W_exp_mean/D_PT_SE(1), "k.-;Stokes-Einstein eq. (4 Pi);");
semilogy (w, D_Am_PT_WC*D_W_exp_mean/D_Am_PT_WC(1), "k-;Wike-Chang eq. mod. of Perkins and Geankoplis (1969);");
semilogy (w, Dim_DM_PT*D_W_exp_mean/Dim_DM_PT(1), "k--;corr. w. MM (M. Diaz 1987);");
semilogy (w, Dim_S_PT*D_W_exp_mean/Dim_S_PT(1), "k-.;corr. w. MM (Sitaraman et al. 1963);");
semilogy (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
semilogy (w_PT_NK, D_PT_NK_T25./Dcorr_NK_PT*D_W_exp_mean/(D_PT_NK_T25(1)./Dcorr_NK_PT(1)),"ks;exp. (Nogami and Kato 1961);");
semilogy (w_Jetal, D_Jetal./Dcorr_Jetal * D_W_exp_mean /(D_Jetal(1) ./ Dcorr_Jetal(1)), "ko;exp. D_F_i_c_k (Jordan et al. 1956);");
semilogy (w_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
semilogy (w_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction Glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T, C and DW corrected")
grid on
box on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [pdir.plot "diffusivity/" "overview_PT_T25_Ccorr_DWcorr"]);


##
## Propylene Glycol
##

## output
plt_1_h = {mfilename, date; "w_PD", "D_AB in m^2/s; model Stokes-Einstein eq. (4*Pi)"};
plt_1_d = [w' D_PD_SE];
cell2csv (["D_water-PD_Stokes-Einstein_h.csv"], plt_1_h)
csvwrite (["D_water-PD_Stokes-Einstein_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
plt_2_h = {mfilename, date; "w_PD", "D_AB in m^2/s; corr. Wilke-Chang-mix"};
plt_2_d = [w' D_Am_PD_WC];
cell2csv (["D_water-PD_PerkinsLR1969et_h.csv"], plt_2_h)
csvwrite (["D_water-PD_PerkinsLR1969et_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
plt_3_h = {mfilename, date; "w_PD", "D_AB in m^2/s; corr. SitaramanR1963etal with mixing rule eta^0.8"};
plt_3_d = [w' Dim_S_PD'];
cell2csv (["D_water-PD_SitaramanR1963etal_h.csv"], plt_3_h)
csvwrite (["D_water-PD_SitaramanR1963etal_d.csv"], plt_3_d, "append", "off", "precision","%.4e")
plt_4_h = {mfilename, date; "w_PD", "D_AB in m^2/s; corr. DiazM1987etal with mixing rule eta^0.8"};
plt_4_d = [w' Dim_DM_PD'];
cell2csv (["D_water-PD_DiazM1987etal_h.csv"], plt_4_h)
csvwrite (["D_water-PD_DiazM1987etal_d.csv"], plt_4_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "w_PD", "D_AB in m^2/s; exp. NogamiH1962et"};
plt_5_d = [w_PD_NK' D_PD_NK'];
cell2csv (["D_water-PD_NogamiH1962et_h.csv"], plt_5_h)
csvwrite (["D_water-PD_NogamiH1962et_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_5_h = {mfilename, date; "w_PD", "D_AB in m^2/s; exp. NogamiH1962et, corr to 25°C"};
plt_5_d = [w_PD_NK' D_PD_NK_T25' D_PD_NK_T25'./Dcorr_NK_PD'];
cell2csv (["D_water-PD_NogamiH1962et_T25_h.csv"], plt_5_h)
csvwrite (["D_water-PD_NogamiH1962et_T25_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "w_PD", "D_AB in m^2/s; exp. this study, first trend"};
plt_5_d = [w_SJG_PD' D_SJG_PD_1'];
cell2csv (["D_water-PD_SJG_1_h.csv"], plt_5_h)
csvwrite (["D_water-PD_SJG_1_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_6_h = {mfilename, date; "w_PD", "D_AB in m^2/s; exp. this study, second trend"};
plt_6_d = [w_SJG_PD' D_SJG_PD_2'];
cell2csv (["D_water-PD_SJG_2_h.csv"], plt_6_h)
csvwrite (["D_water-PD_SJG_2_d.csv"], plt_6_d, "append", "off", "precision","%.4e")

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (w, D_PD_SE, "k.-;Stokes-Einstein eq. (4 Pi);");
semilogy (w, D_Am_PD_WC, "k-;Wike-Chang eq. mod. of Perkins and Geankoplis (1969);");
semilogy (w, Dim_DM_PD, "k--;corr. w. MM (M. Diaz 1987);");
semilogy (w, Dim_S_PD, "k-.;corr. w. MM (Sitaraman et al. 1963);");
semilogy (w_PD_NK, D_PD_NK_T25,"ks;exp. (Nogami and Kato 1961);");
semilogy (w_SJG_PD, D_SJG_PD_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_SJG_PD, D_SJG_PD_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
xlabel ("mass fraction Propylene Glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n T corrected")
grid on
box on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [pdir.plot "diffusivity/" "overview_PD_T25"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (w, D_PD_SE, "k.-;Stokes-Einstein eq. (4 Pi);");
semilogy (w, D_Am_PD_WC, "k-;Wike-Chang eq. mod. of Perkins and Geankoplis (1969);");
semilogy (w, Dim_DM_PD, "k--;corr. w. MM (M. Diaz 1987);");
semilogy (w, Dim_S_PD, "k-.;corr. w. MM (Sitaraman et al. 1963);");
semilogy (w_PD_NK, D_PD_NK_T25./Dcorr_NK_PD,"ks;exp. (Nogami and Kato 1961);");
semilogy (w_SJG_PD, D_SJG_PD_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_SJG_PD, D_SJG_PD_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
xlabel ("mass fraction Propylene Glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n T and C corrected")
grid on
box on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [pdir.plot "diffusivity/" "overview_PD_T25_Ccorr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (w, D_PD_SE*D_W_exp_mean/D_PD_SE(1), "k.-;Stokes-Einstein eq. (4 Pi);");
semilogy (w, D_Am_PD_WC*D_W_exp_mean/D_Am_PD_WC(1), "k-;Wike-Chang eq. mod. of Perkins and Geankoplis (1969);");
semilogy (w, Dim_DM_PD*D_W_exp_mean/Dim_DM_PD(1), "k--;corr. w. MM (M. Diaz 1987);");
semilogy (w, Dim_S_PD*D_W_exp_mean/Dim_S_PD(1), "k-.;corr. w. MM (Sitaraman et al. 1963);");
semilogy (w_PD_NK, D_PD_NK_T25./Dcorr_NK_PD*D_W_exp_mean/(D_PD_NK_T25(1)./Dcorr_NK_PD(1)),"ks;exp. (Nogami and Kato 1961);");
semilogy (w_SJG_PD, D_SJG_PD_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (w_SJG_PD, D_SJG_PD_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
xlabel ("mass fraction Propylene Glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n T, C and DW corrected")
grid on
box on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [pdir.plot "diffusivity/" "overview_PD_T25_Ccorr_DWcorr"]);
