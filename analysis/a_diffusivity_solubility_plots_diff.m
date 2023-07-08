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
fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
##plot (0, D_W_exp_mean, "kd;exp. mean (St. Denis Fell 1971);");
plot (mf, D_PT_SE*D_W_exp_mean/D_PT_SE(1), "k.-;Stokes-Einstein eq. (D_W corrected) (4*Pi);");
##plot (mf, D_PT_SE, "b-.;Stokes-Einstein eq. (4*Pi);");
plot (mf, D_Am_PT_WC*D_W_exp_mean/D_Am_PT_WC(1), "k-;Wike-Chang eq. mixture (Perkins and Geankoplis (1969);");
##plot (mf, D_Am_PT_WC, "k-.;Wike-Chang corr. for mixed solvent (Perkins and Geankoplis (1969);");
plot (mf, Dim_DM_PT*D_W_exp_mean/Dim_DM_PT(1), "k--;corr. with mixture model (D_W corrected) (M Diaz 1987);");
##plot (mf, Dim_DM_PT, "r-.;emperical corr. with mixture model (D_W corrected) (M Diaz 1987);");
##plot (mf, Dim_PT, "k-;S mixture model;");
##plot (1, D_PT_DM, "x-;DM;");
plot (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
##plot (mf_PT_NK, D_PT_NK_T25,"-.x;PT exp NK;");
plot (mf_PT_NK, D_PT_NK_T25*D_W_exp_mean/D_PT_NK_T25(1),"k-.s;exp. (Nogami and Kato 1961);");
##plot (mf_PT_NK, D_PT_NK,"k-.;exp. NK;");
plot (mf_Jetal, D_Jetal ./ Dcorr_Jetal *D_W_exp_mean /(D_Jetal(1) ./ Dcorr_Jetal(1)) ,"k-.o;exp. D_F_i_c_k (Jordan et al. 1956);");
##plot (mf_Jetal, D_Jetal,"k-.o;exp. Jetal;");
##plot (mfa_Jetal, Da_Jetal*D_W_exp_mean/Da_Jetal(1), "k-.x;exp. D_a (Jordan et al. 1956);");
##plot (mfa_Jetal, Da_Jetal, "k-.x;exp. Da Jetal;");
##plot (mf_M_PT, D_M_PT,"-.*;PT exp master thesis;");
plot (mf_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
plot (mf_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol")
grid on
box on
legend ("location", "southoutside")
legend ("location", "eastoutside")
##
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [save_dir "overview_PT"]);
## output
plt_0_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. water @25°C StDenisCE1971et"};
plt_0_d = [0 D_W_exp_mean];
cell2csv (["D_water-StDenisCE1971et_h.csv"], plt_0_h)
csvwrite (["D_water-StDenisCE1971et_d.csv"], plt_0_d, "append", "off", "precision","%.4e")
plt_8_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. water diffusivity in glycerol @25°C ErricoG2004etal"};
plt_8_d = [1 D_PT_DErrico];
cell2csv (["D_water-in-glycerol_ErricoG2004etal_h.csv"], plt_8_h)
csvwrite (["D_water-in-glycerol_ErricoG2004etal_d.csv"], plt_8_d, "append", "off", "precision","%.4e")
##
plt_1_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; model Stokes-Einstein eq. (4*Pi)"};
plt_1_d = [mf' D_PT_SE];
cell2csv (["D_water-PT_Stokes-Einstein_h.csv"], plt_1_h)
csvwrite (["D_water-PT_Stokes-Einstein_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
plt_2_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; corr. Wilke-Chang-mix"};
plt_2_d = [mf' D_Am_PT_WC];
cell2csv (["D_water-PT_PerkinsLR1969et_h.csv"], plt_2_h)
csvwrite (["D_water-PT_PerkinsLR1969et_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
plt_3_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; corr. SitaramanR1963etal with mixing rule eta^0.8"};
plt_3_d = [mf' Dim_S_PT'];
cell2csv (["D_water-PT_SitaramanR1963etal_h.csv"], plt_3_h)
csvwrite (["D_water-PT_SitaramanR1963etal_d.csv"], plt_3_d, "append", "off", "precision","%.4e")
plt_4_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; corr. DiazM1987etal with mixing rule eta^0.8"};
plt_4_d = [mf' Dim_DM_PT'];
cell2csv (["D_water-PT_DiazM1987etal_h.csv"], plt_4_h)
csvwrite (["D_water-PT_DiazM1987etal_d.csv"], plt_4_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. NogamiH1962et"};
plt_5_d = [mf_PT_NK' D_PT_NK'];
cell2csv (["D_water-PT_NogamiH1962et_h.csv"], plt_5_h)
csvwrite (["D_water-PT_NogamiH1962et_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_5_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. NogamiH1962et, corr to 25°C"};
plt_5_d = [mf_PT_NK' D_PT_NK_T25'];
cell2csv (["D_water-PT_NogamiH1962et_T25_h.csv"], plt_5_h)
csvwrite (["D_water-PT_NogamiH1962et_T25_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. this study, first trend"};
plt_5_d = [mf_SJG_PT' D_SJG_PT_1'];
cell2csv (["D_water-PT_SJG_1_h.csv"], plt_5_h)
csvwrite (["D_water-PT_SJG_1_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_6_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. this study, second trend"};
plt_6_d = [mf_SJG_PT' D_SJG_PT_2'];
cell2csv (["D_water-PT_SJG_2_h.csv"], plt_6_h)
csvwrite (["D_water-PT_SJG_2_d.csv"], plt_6_d, "append", "off", "precision","%.4e")
plt_7_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. JordanJ1956etal"};
plt_7_d = [mf_Jetal' D_Jetal'];
cell2csv (["D_water-PT_JordanJ1956etal_h.csv"], plt_7_h)
csvwrite (["D_water-PT_JordanJ1956etal_d.csv"], plt_7_d, "append", "off", "precision","%.4e")
plt_9_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. PLIF JimenezM2012etal-1; corr to 25°C"};
plt_9_d = [mf_PT_J' D_J_3_T25'];
cell2csv (["D_water-PT-E_JimenezM2012etal-1_h.csv"], plt_9_h)
csvwrite (["D_water-PT-E_JimenezM2012etal-1_d.csv"], plt_9_d, "append", "off", "precision","%.4e")
plt_10_h = {mfilename, date; "mf_PT", "D_AB in m^2/s; exp. PLIF KapoustinaV2019etal; corr to 25°C"};
plt_10_d = [mf_PT_K' D_K_T25'];
cell2csv (["D_water-PT-E_KapoustinaV2019etal_h.csv"], plt_10_h)
csvwrite (["D_water-PT-E_KapoustinaV2019etal_d.csv"], plt_10_d, "append", "off", "precision","%.4e")


fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (0, D_W_exp_mean, "kd;exp. mean (St. Denis Fell 1971);");
semilogy (mf, D_PT_SE, "b-.;Stokes-Einstein eq. (4*Pi);");
semilogy (mf, D_Am_PT_WC, "k-.;Wike-Chang corr. for mixed solvent (Perkins and Geankoplis (1969);");
semilogy (mf, Dim_DM_PT, "r-.;emperical corr. with mixture model (D_W corrected) (M Diaz 1987);");
semilogy (1, D_PT_DM, "x-;M Diaz;");
semilogy (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
semilogy (mf_PT_NK, D_PT_NK_T25,"-.x;PT exp NK;");
semilogy (mf_Jetal, D_Jetal,"k-.o;exp. Jetal;");
semilogy (mf_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
semilogy (mf_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T corrected")
grid on
box on
legend ("location", "southoutside")
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PT_log_T-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
plot (0, D_W_exp_mean, "kd;exp. mean (St. Denis Fell 1971);");
plot (mf, D_PT_SE, "b-.;Stokes-Einstein eq. (4*Pi);");
plot (mf, D_Am_PT_WC, "k-.;Wike-Chang corr. for mixed solvent (Perkins and Geankoplis (1969);");
plot (mf, Dim_DM_PT, "r-.;emperical corr. with mixture model (D_W corrected) (M Diaz 1987);");
plot (1, D_PT_DM, "x-;M Diaz;");
plot (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
plot (mf_PT_NK, D_PT_NK_T25,"-.x;PT exp NK;");
plot (mf_Jetal, D_Jetal,"k-.o;exp. Jetal;");
plot (mf_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
plot (mf_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T corrected")
grid on
box on
legend ("location", "southoutside")
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PT_T-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (0, D_W_exp_mean, "kd;exp. mean (St. Denis Fell 1971);");
semilogy (mf, D_PT_SE, "b-.;Stokes-Einstein eq. (4*Pi);");
semilogy (mf, D_Am_PT_WC, "k-.;Wike-Chang corr. for mixed solvent (Perkins and Geankoplis (1969);");
semilogy (mf, Dim_DM_PT, "r-.;emperical corr. with mixture model (D_W corrected) (M Diaz 1987);");
semilogy (1, D_PT_DM, "x-;M Diaz;");
semilogy (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
semilogy (mf_PT_NK, D_PT_NK_T25./Dcorr_NK,"-.x;PT exp NK;");
semilogy (mf_Jetal, D_Jetal./Dcorr_Jetal,"k-.o;exp. Jetal;");
semilogy (mf_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
semilogy (mf_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T and C corrected")
grid on
box on
legend ("location", "southoutside")
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PT_log_T-C-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
plot (0, D_W_exp_mean, "kd;exp. mean (St. Denis Fell 1971);");
plot (mf, D_PT_SE, "b-.;Stokes-Einstein eq. (4*Pi);");
plot (mf, D_Am_PT_WC, "k-.;Wike-Chang corr. for mixed solvent (Perkins and Geankoplis (1969);");
plot (mf, Dim_DM_PT, "r-.;emperical corr. with mixture model (D_W corrected) (M Diaz 1987);");
plot (1, D_PT_DM, "x-;M Diaz;");
plot (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
plot (mf_PT_NK, D_PT_NK_T25./Dcorr_NK,"-.x;PT exp NK;");
plot (mf_Jetal, D_Jetal./Dcorr_Jetal,"k-.o;exp. Jetal;");
plot (mf_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
plot (mf_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T and C corrected")
grid on
box on
legend ("location", "southoutside")
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PT_T-C-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
plot (0, D_W_exp_mean, "kd;exp. mean (St. Denis Fell 1971);");
plot (mf, D_PT_SE, "b-.;Stokes-Einstein eq. (4*Pi);");
plot (mf, D_Am_PT_WC, "k-.;Wike-Chang corr. for mixed solvent (Perkins and Geankoplis (1969);");
plot (mf, Dim_DM_PT, "r-.;emperical corr. with mixture model (D_W corrected) (M Diaz 1987);");
plot (1, D_PT_DM, "x-;DM;");
plot (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
plot (mf_PT_NK, D_PT_NK_T25,"-.x;PT exp NK;");
plot (mf_Jetal, D_Jetal,"k-.o;exp. Jetal;");
##plot (mfa_Jetal, Da_Jetal, "k-.x;exp. Da Jetal;");
plot (mf_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
plot (mf_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
plot (mf_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n temperature corrected")
grid on
box on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PT_T-C-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (mf, D_PT_SE*D_W_exp_mean/D_PT_SE(1), "k.-;Stokes-Einstein eq. (D_W corrected) (4*Pi);");
semilogy (mf, D_Am_PT_WC*D_W_exp_mean/D_Am_PT_WC(1), "k-;Wike-Chang eq. mixture (Perkins and Geankoplis (1969);");
semilogy (mf, Dim_DM_PT*D_W_exp_mean/Dim_DM_PT(1), "k--;corr. with mixture model (D_W corrected) (M Diaz 1987);");
semilogy (1, D_water_glycerol, "kd;exp. water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
semilogy (mf_PT_NK, D_PT_NK_T25./Dcorr_NK*D_W_exp_mean/(D_PT_NK_T25(1)./Dcorr_NK(1)),"k-.s;exp. (Nogami and Kato 1961);");
semilogy (mf_Jetal, D_Jetal./Dcorr_Jetal * D_W_exp_mean /(D_Jetal(1) ./ Dcorr_Jetal(1)), "k-.o;exp. D_F_i_c_k (Jordan et al. 1956);");
semilogy (mf_SJG_PT, D_SJG_PT_1,"bs*;exp. this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_SJG_PT, D_SJG_PT_2,"b+;exp. this study, second trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_PT_J, D_J_3_T25,"m*;exp. with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
semilogy (mf_PT_K, D_K_T25,"mo;exp. with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n T, D_W and C corrected")
grid on
box on
legend ("location", "southoutside")
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PT_corr_log"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (mf, D_PT_SE, "k.-;model @25°C, Stokes-Einstein eq. (4*Pi);");
semilogy (mf, D_Am_PT_WC, "k-;corr. @25°C, Wike-Chang mix (Perkins and Geankoplis (1969);");
semilogy (mf, Dim_DM_PT, "k--;corr. @25°C, with mixture model (M Diaz 1987);");
semilogy (1, D_water_glycerol, "kd;exp. @25°C, water (as solute) in glycerol (D'Errico et al. 2004);", "markersize", 5,"linewidth", 2)
semilogy (mf_PT_NK, D_PT_NK, "k-.s;exp. @22.5°C, (Nogami and Kato 1961);");
semilogy (mf_Jetal, D_Jetal, "ko-.;exp. @25°C, D_F_i_c_k (Jordan et al. 1956);");
semilogy (mf_SJG_PT, D_SJG_PT_1, "bs*;exp. @25°C, this study, fist trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_SJG_PT, D_SJG_PT_2, "b+;exp. @25°C, this study, second trend;", "markersize", 8, "linewidth", 1.5);
semilogy (mf_PT_J, D_J_3, "m*;exp. @20°C, with 20% Ethanol (Jimenez et al. 2012) ;", "markersize", 5,"linewidth", 1.5);
semilogy (mf_PT_K, D_K, "mo;exp. @28°C, with 30% Ethanol (Kapoustina et al. 2019);", "markersize", 5,"linewidth", 2);
xlabel ("mass fraction glycerol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Glycerol \n original data")
grid on
box on
legend ("location", "southoutside")
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PT_log_org"]);



##
## Propylene Glycol
##

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
plot (0, D_W_exp_mean, "kd;exp. mean St. Denis Fell 1971;");
plot (mf, D_PD_SE, "k-..;Stokes-Einstein eq. (4*Pi);");
plot (mf, D_Am_PD_WC, "m-;corr. @25°C, Wike-Chang mix (Perkins and Geankoplis (1969);");
plot (mf, Dim_S_PD, "g-;corr. @25°C, SitaramanR1963etal;");
plot (mf, Dim_DM_PD, "k-;corr. @25°C, DiazM1987etal;");
##plot (1, D_PD_DM, "x-;@25°C, DM;");
plot (mf_PD_NK, D_PD_NK,"-.*;exp. @22.5°C NogamiH1962et;");
plot (mf_SJG_PD, D_SJG_PD_1,"ms;exp. this study @25°C,;");
plot (mf_SJG_PD, D_SJG_PD_2,"ms;exp. this study @25°C,;");
xlabel ("mass fraction propylene glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n original data")
grid on
legend ("location", "eastoutside")
##
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [save_dir "overview_PD_org"]);
## output
plt_1_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; model Stokes-Einstein eq. (4*Pi)"};
plt_1_d = [mf' D_PD_SE];
cell2csv (["D_water-PD_Stokes-Einstein_h.csv"], plt_1_h)
csvwrite (["D_water-PD_Stokes-Einstein_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
plt_2_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; corr. Wilke-Chang-mix"};
plt_2_d = [mf' D_Am_PD_WC];
cell2csv (["D_water-PD_PerkinsLR1969et_h.csv"], plt_2_h)
csvwrite (["D_water-PD_PerkinsLR1969et_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
plt_3_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; corr. SitaramanR1963etal with mixing rule eta^0.8"};
plt_3_d = [mf' Dim_S_PD'];
cell2csv (["D_water-PD_SitaramanR1963etal_h.csv"], plt_3_h)
csvwrite (["D_water-PD_SitaramanR1963etal_d.csv"], plt_3_d, "append", "off", "precision","%.4e")
plt_4_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; corr. DiazM1987etal with mixing rule eta^0.8"};
plt_4_d = [mf' Dim_DM_PD'];
cell2csv (["D_water-PD_DiazM1987etal_h.csv"], plt_4_h)
csvwrite (["D_water-PD_DiazM1987etal_d.csv"], plt_4_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; exp. NogamiH1962et"};
plt_5_d = [mf_PD_NK' D_PD_NK'];
cell2csv (["D_water-PD_NogamiH1962et_h.csv"], plt_5_h)
csvwrite (["D_water-PD_NogamiH1962et_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_5_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; exp. NogamiH1962et, corr to 25°C"};
plt_5_d = [mf_PD_NK' D_PD_NK_T25'];
cell2csv (["D_water-PD_NogamiH1962et_T25_h.csv"], plt_5_h)
csvwrite (["D_water-PD_NogamiH1962et_T25_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
##
plt_5_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; exp. this study, first trend"};
plt_5_d = [mf_SJG_PD' D_SJG_PD_1'];
cell2csv (["D_water-PD_SJG_1_h.csv"], plt_5_h)
csvwrite (["D_water-PD_SJG_1_d.csv"], plt_5_d, "append", "off", "precision","%.4e")
plt_6_h = {mfilename, date; "mf_PD", "D_AB in m^2/s; exp. this study, second trend"};
plt_6_d = [mf_SJG_PD' D_SJG_PD_2'];
cell2csv (["D_water-PD_SJG_2_h.csv"], plt_6_h)
csvwrite (["D_water-PD_SJG_2_d.csv"], plt_6_d, "append", "off", "precision","%.4e")



fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (0, D_W_exp_mean, "kd;exp. mean St. Denis Fell 1971;");
semilogy (mf, D_PD_SE, "k-..;Stokes-Einstein eq. (4*Pi);");
semilogy (mf, D_Am_PD_WC, "m-;corr. @25°C, Wike-Chang mix (Perkins and Geankoplis (1969);");
semilogy (mf, Dim_S_PD, "g-;corr. @25°C, SitaramanR1963etal;");
semilogy (mf, Dim_DM_PD, "k-;corr. @25°C, DiazM1987etal;");
##semilogy (1, D_PD_DM, "x-;@25°C, DM;");
semilogy (mf_PD_NK, D_PD_NK,"-.*;exp. @22.5°C NK;");
semilogy (mf_SJG_PD, D_SJG_PD_1,"ms;exp. this study @25°C,;");
semilogy (mf_SJG_PD, D_SJG_PD_2,"ms;exp. this study @25°C,;");
xlabel ("mass fraction propylene glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n original data")
grid on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PD_log_org"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (0, D_W_exp_mean, "kd;exp. mean St. Denis Fell 1971;");
semilogy (mf, D_PD_SE, "k-..;Stokes-Einstein eq. (4*Pi);");
semilogy (mf, D_Am_PD_WC, "m-;corr. @25°C, Wike-Chang mix (Perkins and Geankoplis (1969);");
semilogy (mf, Dim_S_PD, "g-;corr. @25°C, SitaramanR1963etal;");
semilogy (mf, Dim_DM_PD, "k-;corr. @25°C, DiazM1987etal;");
##semilogy (1, D_PD_DM, "x-;@25°C, DM;");
semilogy (mf_PD_NK, D_PD_NK_T25,"-.*;exp. @22.5°C NK;");
semilogy (mf_SJG_PD, D_SJG_PD_1,"ms;exp. this study @25°C,;");
semilogy (mf_SJG_PD, D_SJG_PD_2,"ms;exp. this study @25°C,;");
xlabel ("mass fraction propylene glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n T corrected")
grid on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PD_log_T-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
plot (0, D_W_exp_mean, "kd;exp. mean St. Denis Fell 1971;");
plot (mf, D_PD_SE, "k-..;Stokes-Einstein eq. (4*Pi);");
plot (mf, D_Am_PD_WC, "m-;corr. @25°C, Wike-Chang mix (Perkins and Geankoplis (1969);");
plot (mf, Dim_S_PD, "g-;corr. @25°C, SitaramanR1963etal;");
plot (mf, Dim_DM_PD, "k-;corr. @25°C, DiazM1987etal;");
##plot (1, D_PD_DM, "x-;@25°C, DM;");
plot (mf_PD_NK, D_PD_NK_T25,"-.*;exp. @22.5°C NK;");
plot (mf_SJG_PD, D_SJG_PD_1,"ms;exp. this study @25°C,;");
plot (mf_SJG_PD, D_SJG_PD_2,"ms;exp. this study @25°C,;");
xlabel ("mass fraction propylene glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n T corrected")
grid on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PD_T-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
semilogy (0, D_W_exp_mean, "kd;exp. mean St. Denis Fell 1971;");
semilogy (mf, D_PD_SE, "k-..;Stokes-Einstein eq. (4*Pi);");
semilogy (mf, D_Am_PD_WC, "m-;corr. @25°C, Wike-Chang mix (Perkins and Geankoplis (1969);");
semilogy (mf, Dim_S_PD, "g-;corr. @25°C, SitaramanR1963etal;");
semilogy (mf, Dim_DM_PD, "k-;corr. @25°C, DiazM1987etal;");
##semilogy (1, D_PD_DM, "x-;@25°C, DM;");
semilogy (mf_PD_NK, D_PD_NK_T25 * D_W_exp_mean / D_PD_NK_T25(1),"-.*;exp. @22.5°C NK;");
semilogy (mf_SJG_PD, D_SJG_PD_1,"ms;exp. this study @25°C,;");
semilogy (mf_SJG_PD, D_SJG_PD_2,"ms;exp. this study @25°C,;");
xlabel ("mass fraction propylene glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n T corrected")
grid on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PD_log_T_D_W-corr"]);

fh = figure (); hold on;
set (fh, "DefaultAxesFontSize", figp.fs_min-1, "DefaultAxesFontName", figp.fonts{1});
set (fh, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
set (fh, "PaperSize", [fig_size_x fig_size_y]);
set (fh, "renderer", "painters");
plot (0, D_W_exp_mean, "kd;exp. mean St. Denis Fell 1971;");
plot (mf, D_PD_SE, "k-..;Stokes-Einstein eq. (4*Pi);");
plot (mf, D_Am_PD_WC, "m-;corr. @25°C, Wike-Chang mix (Perkins and Geankoplis (1969);");
plot (mf, Dim_S_PD, "g-;corr. @25°C, SitaramanR1963etal;");
plot (mf, Dim_DM_PD, "k-;corr. @25°C, DiazM1987etal;");
##plot (1, D_PD_DM, "x-;@25°C, DM;");
plot (mf_PD_NK, D_PD_NK_T25 * D_W_exp_mean / D_PD_NK_T25(1),"-.*;exp. @22.5°C NK;");
plot (mf_SJG_PD, D_SJG_PD_1,"ms;exp. this study @25°C,;");
plot (mf_SJG_PD, D_SJG_PD_2,"ms;exp. this study @25°C,;");
xlabel ("mass fraction propylene glycol in g/g")
ylabel ("D in m^2/s")
title ("Diffusity of Oxygen in Aqueous Solutions of Propylene Glycol \n T corrected")
grid on
legend ("location", "eastoutside")
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [pdir.plot "diffusivity/" "overview_PD_T-D_W-corr"]);

