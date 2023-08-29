##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Collection of fluid properties relevant to this study. Script to generate
## "fluidprop.v7" database.
## Original data points of own measurements and of available literature.
## Credit to the original authors is given in the field "source" of each
## dataset "ds"
##
## Author: Sören J. Gerke
##

##
## datasets
##

##
tables_dir = pdir.fptab;

##
fp = init_fp ();

##
## Water
##
ds = init_fp ();
ds.fluid = {"water", "W", "H2O"};
ds.prop{1} = {"molar mass", "mm"};
ds.unit(1) = {"kg / mol"};
ds.data{1} = 1e-3 * 18.0153;
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"water", "W", "H2O"};
ds.prop{1} = {"refractive-index", "RI"};
ds.source{1} = {"GerkeSJexp"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30] + 273.15};
ds.unit(2) = {"RI in -"};
ds.data(2) = {[1.33299 1.33265 1.33217]};
ds.desc = {"3 temperatures for destilled water, measured with Abbe refractometer"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"water", "W", "H2O"};
ds.prop{1} = {"density", "rho"};
ds.source{1} = {"IAPWS 2008 (https://wiki.anton-paar.com/us-en/water/)"};
ds.unit(1) = {"T in K"};
ds.data(1) = {273.15 + [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 45 50 55 60 65 70 75 80]};
ds.unit(2) = {"density in kg/m^3"};
ds.data(2) = {1e3 * [0.9999 1.0 1.0 1.0 0.9999 0.9999 0.9999 0.9998 0.9997 0.9996 0.9995 0.9994 0.9992 0.9991 0.9989 0.9988 0.9986 0.9984 0.9982 998 0.9978 0.9975 0.9973 997 0.9968 0.9965 0.9962 0.9959 0.9956 0.9953 995 0.9947 0.9944 994 0.9937 0.9933 993 0.9926 0.9922 0.9902 988 0.9857 0.9832 0.9806 0.9778 0.9748 0.9718]};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"water", "W", "H2O"};
ds.prop{1} = {"dynamic viscosity", "eta"};
ds.source{1} = {"IAPWS 2008 (https://wiki.anton-paar.com/us-en/water/)"};
ds.unit(1) = {"T in K"};
ds.data(1) = {273.15 + [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 45 50 55 60 65 70 75 80]};
ds.unit(2) = {"dynamic viscosity in Pa s"};
ds.data(2) = {1e-3 * [1.6735 1.6190 1.5673 1.5182 1.4715 1.4271 1.3847 1.3444 1.3059 1.2692 1.2340 1.2005 1.1683 1.1375 1.1081 1.0798 1.0526 1.0266 1.0016 0.9775 0.9544 0.9321 0.9107 0.8900 0.8701 0.8509 0.8324 0.8145 0.7972 0.7805 0.7644 0.7488 0.7337 0.7191 0.7050 0.6913 0.6780 0.6652 0.6527 0.5958 0.5465 0.5036 0.4660 0.4329 0.4035 0.3774 0.3540]};
[fp, idx] = add_fp_dataset (fp, ds);

##
## Glycerol
##
ds = init_fp ();
ds.fluid =  {"glycerol", "PT", "1,2,3-Propantriol"};
ds.prop{1} = {"molar mass", "mm"};
ds.unit(1) = {"kg / mol"};
ds.data{1} = 1e-3 * 92.094;
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"glycerol", "PT", "1,2,3-Propantriol"};
ds.prop{1} = {"refractive-index", "RI"};
ds.source{1} = {"GerkeSJexp"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30] + 273.15};
ds.unit(2) = {"RI in -"};
ds.data(2) = {[1.47428 1.47313 1.47193]};
ds.desc = {"3 temperatures for glycerol (99.5%), measured with Abbe refractometer, mean of 3 repetitions"};
[fp, idx] = add_fp_dataset (fp, ds);

##
## Propylene Glycol
##
ds = init_fp ();
ds.fluid = {"propylene glycol", "PD", "1,2-Propandiol"};
ds.prop{1} = {"molar mass", "mm"};
ds.unit(1) = {"kg / mol"};
ds.data{1} = 1e-3 * 76.095;
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol", "PD", "1,2-Propandiol"};
ds.prop{1} = {"refractive-index", "RI"};
ds.source{1} = {"GerkeSJexp"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30] + 273.15};
ds.unit(2) = {"RI in -"};
ds.data(2) = {[1.43298	1.43134	1.42963]};
ds.desc = {"3 temperatures for glycerol (99.5%), measured with Abbe refractometer, mean of 3 repetitions"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol", "PD"};
ds.prop{1} = {"dynamic-viscosity", "eta"};
ds.source{1} = {"SunT2004et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[297.25 302.55 317.10 332.30]};
ds.unit(2) = {"dynamic viscosity in Pa s"};
ds.data(2) = {1e-3 *[46.4 34.0 16.4 8.75]};
[fp, idx] = add_fp_dataset (fp, ds);

##
## Glycerol - Water
##
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"refractive-index", "RI"};
ds.source{1} = {"GerkeSJexp"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30]' + 273.15};
ds.unit(2) = {"mass fraction glycerol in g / g"};
ds.data(2) = {1e-2 * [75 70 65 60 55 50 45 40]};
ds.unit(3) = {"RI in -"};
ds.data(3) = {[[1.43536 1.42792 1.42054 1.41315 1.40565 1.39850 1.39174 1.38447];
               [1.43426 1.42682 1.41953 1.41221 1.40459 1.39758 1.39086 1.38369];
               [1.43314 1.42571 1.41843 1.41100 1.40346 1.39669 1.38973 1.38300]]};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"dynamic-viscosity", "eta"};
ds.source{1} = {"GerkeSJexp"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30]' + 273.15};
ds.unit(2) = {"mass fraction glycerol in g/g"};
ds.data(2) = {1e-2 * [75 70 65 60 55 50 45 40]};
ds.unit(3) = {"dynamic viscosity in Pa s"};
ds.data(3) = {1e-3 * [[33.7 21.5 15.3 10.9 7.68 6.22 4.95 4.01];
                      [26.3 17.5 12.4 9 6.45 5.3 4.22 3.48];
                      [21.1 14.2 10.3 7.58 5.47 4.55 3.68 3.04];]};
ds.desc = {"glycerol (99.5%) - water(dest.), measured with rotational cylinder viscometer"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"density", "rho"};
ds.source{1} = {"BosartLW1927et", "BOSART,L. W., and SNODDY,A. O., Ind. Eng. Chem., 19, 506-510 (1927)"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30]' + 273.15};
ds.unit(2) = {"mass fraction glycerol in g / g"};
ds.data(2) = {1:-0.1:0};
ds.unit(3) = {"density in kg / m^3"};
ds.data(3) = {1e-3 * [[1.26108 1.23510 1.20850 1.18125 1.15380 1.12630 1.09930 1.07270 1.04690 1.02210 0.99823];
                      [1.25802 1.23200 1.20545 1.17840 1.15105 1.12375 1.09710 1.07070 1.04525 1.02070 0.99708];
                      [1.25495 1.22890 1.20240 1.17565 1.14830 1.12110 1.09475 1.06855 1.04350 1.01905 0.99568];]};
ds.desc = {"literature"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"viscosity", "eta"};
ds.source{1} = {"Segur1951et", "SEGUR.J. B.. and O ~ ~ R S T AHELEN. R. Ind. & Enn. Chem.,~.43. 2117-2120 (1951)"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 30]' + 273.15};
ds.unit(2) = {"mass fraction glycerol in g / g"};
ds.data(2) = {[0:0.1:0.6 0.65:0.05:1]};
ds.unit(3) = {"dynamic viscosity in Pa s"};
ds.data(3) = {1e-3 * [[1.005 1.31 1.76 2.50 3.72 6.00 10.8 15.2 22.5 35.5 60.1 109 219 523 1410];
                      [.8007 1.03 1.35 1.87 2.72 4.21 7.19 9.85 14.1 21.2 33.9 58 109 237 612];]};
ds.desc = {"literature"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"refractive-index", "RI"};
ds.source{1} = {"HoytLF1933", "HOYT,L. F., Oil & Soap, 10, 43-47 (1933); Ind. Eng. Chem., 26.329-332 (1934)"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20]' + 273.15};
ds.unit(2) = {"mass fraction glycerol in g / g"};
ds.data(2) = {1:-0.05:0};
ds.unit(3) = {"RI in -"};
ds.data(3) = {[1.47399 1.46597 1.45839 1.45085 1.44290 1.43534 1.42789 1.42044 1.41299 1.40554 1.39809 1.39089 1.38413 1.37740 1.37070 1.36404 1.35749 1.35106 1.34481 1.33880 1.33303];};
ds.desc = {"literature, polyfit(mf,RI, 3) works well"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"density", "rho"};
ds.source{1} = {"VolkA2018et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[0 100]' + 273.15};
ds.unit(2) = {"mass fraction glycerol in g / g"};
ds.data(2) = {[0 1]};
ds.unit(3) = {"density in kg / m^3"};
##ds.data(3) = {"rho_PT_W_model.m"};
ds.desc = {"model, open literature"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"oxygen diffusivity", "D_AB"};
ds.source{1} = {"NogamiH1962et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[22.5]' + 273.15};
ds.unit(2) = {"volume fraction glycerol in l / l"};
ds.data(2) = {1e-4 * [0 980 3840 5755 7750]};
ds.unit(3) = {"diffusivity in m^2 / s"};
ds.data(3) = {1e-4 * 1e-5 * 1e-2 * [165 116 42 14 5]};
ds.desc = {"literature; polarographic measurements"};
[fp, idx] = add_fp_dataset (fp, ds);
ds = init_fp ();
ds.fluid = {"glycerol-water", "PT_W"};
ds.prop{1} = {"oxygen solubility", "c_sat"};
ds.source{1} = {"NogamiH1962et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[22.5]' + 273.15};
ds.unit(2) = {"volume fraction glycerol in l / l"};
ds.data(2) = {1e-4 * [0 980 3840 5755 7750]};
ds.unit(3) = {"dissolved oxygen in mmol / l"};
ds.data(3) = {1e-3 * 1e-2 * [28 26 16 12 7]};
ds.desc = {"literature"};
[fp, idx] = add_fp_dataset (fp, ds);

##
## Propylene Glycol - Water
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"refractive-index", "RI"};
ds.source{1} = {"GerkeSJexp"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30]' + 273.15};
ds.unit(2) = {"mass fraction propylene glycol in g / g"};
ds.data(2) = {1e-2 * [95 90 85 80 70 60 50]};
ds.unit(3) = {"RI in -"};
ds.data(3) = {[[1.42918 1.42573 1.4217 1.41774 1.40905 1.39915 1.38933];
               [1.42763 1.42394 1.42004 1.4162 1.40751 1.39784 1.38809];
               [1.42608 1.42232 1.4185 1.41468 1.40604 1.39645 1.38671]]};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"dynamic-viscosity", "eta"};
ds.source{1} = {"GerkeSJexp"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30]' + 273.15};
ds.unit(2) = {"mass fraction propylene glycol in g / g"};
ds.data(2) = {1e-2 * [95 90 85 80 70 60 50]};
ds.unit(3) = {"dynamic viscosity in Pa s"};
ds.data(3) = {1e-3 *[[41 32.9 25.5 20.8 13.8 9.7 6.65];
                     [32.6 24.5 20.1 16.4 11.2 7.92 5.52];
                     [25.4 20.3 16 13.3 9.14 6.6 4.62]]};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"dynamic-viscosity", "eta"};
ds.source{1} = {"SunT2004et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[[297];
               [313];]};
ds.unit(2) = {"mole fraction propylene glycol in mol / mol"};
ds.data(2) = [0.25 0.5 0.75];
ds.unit(3) = {"dynamic viscosity in Pa s"};
ds.data(3) = {1e-3 *[[7.39 17.5 32.4];
                     [3.90 8.47 14.4];]};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"dynamic-viscosity", "eta"};
ds.source{1} = {"KhattabIS2017etal"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30]' + 273};
ds.unit(2) = {"mole fraction propylene glycol in mol / mol"};
ds.data(2) = [0.000 0.027 0.058 0.095 0.141 0.197 0.269 0.364 0.495 0.688 1.000];
ds.unit(3) = {"dynamic viscosity in Pa s"};
ds.data(3) = {1e-3 *[[1.003	1.434	1.975	2.78 3.826	5.325	9.235	12.021 18.332 29.494 57.571];
                     [0.976	1.17 1.569	2.346	3.331	4.529	6.757	9.864	14.508 22.908 39.436];
                     [0.937	1.026	1.396	1.976	2.727	3.513	5.086	7.424	10.524 16.043 26.852]]};
ds.desc = {"literature"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"density", "rho"};
ds.source{1} = {"KhattabIS2017etal"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[20 25 30]' + 273};
ds.unit(2) = {"mole fraction propylene glycol in mol / mol"};
ds.data(2) = [0.000 0.027 0.058 0.095 0.141 0.197 0.269 0.364 0.495 0.688 1.000];
ds.unit(3) = {"density in kg / m^3"};
ds.data(3) = {[[997.8	1005.1	1013.6	1022.7	1031.2	1036.5	1040.8	1042.7	1042.5	1040.1	1035.3];
               [995.8	1003.1	1011.3	1020.2	1028.2	1033.3	1037.6	1039.3	1038.9	1036.7	1032.3];
               [993.8	1000.7	1008.6	1016.6	1024.2	1029.4	1032.9	1034.7	1034.2	1033.5	1027.6];]};
ds.desc = {"literature"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"density", "rho"};
ds.source{1} = {"NakanishiK1967etal", "Koichiro Nakanishi, Nobuyuki Kato, Masato Maruyama: J. physic. Chem. 71(1967) 814/B"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[25]' + 273.15};
ds.unit(2) = {"mass fraction propylene glycol in g / g"};
ds.data(2) = {1e-2 * [0	3.617	8.03 14.39 15.949 18.081 19.846	21.782 28.554	35.798 40.26 44.649	47.054 54.121 56.879 60.453 64.032 69.897	80.688 82.536	90.89	100]};
ds.unit(3) = {"density in kg / m^3"};
ds.data(3) = {[997 999.39 1002.65 1007.8 1009.06 1010.98 1012.49 1014.16 1019.95 1025.58 1028.95 1031.56 1033.24 1036.61 1037.55 1038.81 1039.54 1040.18 1039.83 1039.46 1037.06 1032.41]};
ds.desc = {"literature; pycnometer measurements"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"dynamic-viscosity", "eta"};
ds.source{1} = {"TanakaY1988etal"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[25]' + 273.15};
ds.unit(2) = {"mass fraction propylene glycol in g / g"};
ds.data(2) = {[0 0.5136 0.7379 0.8637 0.9441 1.0000]};
ds.unit(3) = {"dynamic viscosity in Pa s"};
ds.data(3) = {1e-3 *[0.891 5.452 12.77 21.63 32.85 44.39]};
ds.desc = {"literature; falling cylinder measurements"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"dynamic-viscosity", "eta"};
ds.source{1} = {"GeorgeJ2003et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[25 35]' + 273.15};
ds.unit(2) = {"mole fraction propylene glycol in mol / mol"};
ds.data(2) = {1 - [0 0.0440 0.0987 0.1484 0.2498 0.3494 0.4494 0.5029 0.5499  0.5597 0.6501 0.7499 0.8500 0.9003  0.9115 0.9154 0.9270  0.9367 0.9504 0.9551 0.9600  0.9696 0.9797 0.9902 0.9923 0.9948 1]};
ds.unit(3) = {"dynamic viscosity in Pa s"};
ds.data(3) = {1e-3 *[[43.428 41.249 38.410 35.744 30.202 24.829 19.708 17.139 14.997 14.565 10.839 7.287 4.325 3.049 2.782 2.692 2.425 2.208 1.907 1.806 1.701 1.501 1.294 1.083 1.041 0.992 0.890];
                     [24.247 23.027 21.431 19.930 16.803 13.774 10.900 9.465 8.274 8.034 5.979 4.049 2.473 1.810 1.673 1.626 1.491 1.379 1.227 1.177 1.124 1.023 0.919 0.815 0.794 0.769 0.719];]};
ds.desc = {"literature; Ubbelohde viscometer measurements"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"density", "rho"};
ds.source{1} = {"GeorgeJ2003et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[25 35]' + 273.15};
ds.unit(2) = {"mole fraction propylene glycol in mol / mol"};
ds.data(2) = {1 - [0.0440 0.0987 0.1484 0.2498 0.3494 0.4494 0.5029 0.5499  0.5597 0.6501 0.7499 0.8500 0.9003  0.9115 0.9154 0.9270  0.9367 0.9504 0.9551 0.9600  0.9696 0.9797 0.9902 0.9923 0.9948]};
ds.unit(3) = {"density in kg / m^3"};
ds.data(3) = {1e3 * [[1.03339 1.03415 1.03484 1.03630 1.03780 1.03924 1.03987 1.04027 1.04033 1.04020  1.03773 1.03031 1.02326 1.02126 1.02051 1.01816 1.01602 1.01269 1.01146  1.01013 1.00736 1.00421 1.00065 0.99990 0.99899];
                     [1.02601 1.02676 1.02745 1.02892 1.03043 1.03191 1.03259 1.03305 1.03312 1.03319 1.03114 1.02453 1.01814 1.01630 1.01563 1.01347 1.01151 1.00845 1.00732 1.00610 1.00355 1.00065 0.99736 0.99667 0.99583];]};
ds.desc = {"literature; vibrating tube digital densimeter measurements"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"density", "rho"};
ds.source{1} = {"MacBethG1951et", "MacBeth, G., Thompson, A.R., 1951. Densities and refractive indexes for propylene glycol–water solutions. Anal. Chem. 23, 618–619"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[35]' + 273.15};
ds.unit(2) = {"mass fraction propylene glycol in g / g"};
ds.data(2) = {0:.1:1};
ds.unit(3) = {"density in kg / m^3"};
ds.data(3) = {1e3 * [0.99406 1.00089 1.00852 1.01615 1.02295 1.02818 1.03154 1.03299 1.03240 1.02982 1.02510]};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"refractive-index", "RI"};
ds.source{1} = {"MacBethG1951et", "MacBeth, G., Thompson, A.R., 1951. Densities and refractive indexes for propylene glycol–water solutions. Anal. Chem. 23, 618–619"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[25]' + 273.15};
ds.unit(2) = {"mass fraction propylene glycol in g / g"};
ds.data(2) = {0:.1:1};
ds.unit(3) = {"RI in -"};
ds.data(3) = {[1.3325 1.3432 1.3544 1.3658 1.3770 1.3878 1.3980 1.4075 1.4162 1.4241 1.4316]};
ds.desc = {"literature, polyfit(mf,RI, 3) works well"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"oxygen diffusivity", "D_AB"};
ds.source{1} = {"NogamiH1962et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[22.5]' + 273.15};
ds.unit(2) = {"volume fraction propylene glycol in l / l"};
ds.data(2) = {1e-4 * [0 975 1980 2960 3900 4940 5960 6890 7885 8905 10000]};
ds.unit(3) = {"diffusivity in m^2 / s"};
ds.data(3) = {1e-4 * 1e-5 * 1e-2 * [183 136 101 82 68 56 43 37 33 34 40]};
ds.desc = {"literature; polarographic measurements"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"oxygen solubility", "c_sat"};
ds.source{1} = {"NogamiH1962et"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[22.5]' + 273.15};
ds.unit(2) = {"volume fraction propylene glycol in l / l"};
ds.data(2) = {1e-4 * [0 975 1980 2960 3900 4940 5960 6890 7885 8905 10000]};
ds.unit(3) = {"dissolved oxygen in mol / l"};
ds.data(3) = {1e-3 * 1e-3 * [295 272 248 219 201 176 159 211 251 302 366]};
ds.desc = {"literature"};
[fp, idx] = add_fp_dataset (fp, ds);
##
ds = init_fp ();
ds.fluid = {"propylene glycol-water", "PD_W"};
ds.prop{1} = {"oxygen solubility", "c_sat"};
ds.source{1} = {"YamamotoH1994etal"};
ds.unit(1) = {"T in K"};
ds.data(1) = {[25]' + 273.15};
ds.unit(2) = {"mole fraction propylene glycol in mol / mol"};
ds.data(2) = {[0 0.029 0.0606 0.0961 0.1398 0.1996 0.2702 0.3784 0.4991 0.6854 1]};
ds.unit(3) = {"Ostwald coefficient in m^3/m^3"};
ds.data(3) = {[0.03161 0.02911 0.027 0.02443 0.02341 0.02364 0.02575 0.02977 0.03585 0.04694 0.05816]};
ds.desc = {"literature"};
[fp, idx] = add_fp_dataset (fp, ds);

##
##
##
save ("-v7", [tables_dir "fluidprop.v7"], "fp");

