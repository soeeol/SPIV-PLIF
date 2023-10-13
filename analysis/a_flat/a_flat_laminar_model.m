##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## laminar flat film theoretical solution
##
##
## Author: Sören J. Gerke
##

## input values:
##
## theoretical concentration field for measurement domain
x_eq = 1e-3 * (x_abs_meas + [-12+5e-3-20:5e-3:20+20]);
y_eq = 1e-3 * [0:1e-3:max(pd_M)];
[XX_eq, YY_eq] = meshgrid (x_eq, y_eq);
## fluid properties:
##  D_eq = D_AB.PLIF1;
D_eq = D_AB.PLIF2;
eta_eq = eta_exp;
nu_eq = eta_eq / rho_exp;
Re_eq = re_l_exp;
Sc_eq = nd_sc (nu_eq, D_eq);

## u and c distribution calculation:
y_u_eq = u_py_eq = cell ();
[h_eq, u_s_eq, u_avg_eq] = model_filmflow_laminar_u_profile_p (nu_eq, deg2rad(aid.ids_A), Re_eq)
for i_M = it_M
  y_u_eq{i_M} = linspace (0, h_eq(i_M), 100);
  u_py_eq{i_M} = model_filmflow_laminar_u_profile (y_u_eq{i_M}, u_s_eq(i_M), h_eq(i_M));
endfor
c_eq = dcdy_eq = [];
for i_M = it_M
  [c_eq{i_M} dcdy_eq{i_M}] = model_filmflow_laminar_c_field (x_eq, y_eq, c_i=1, c_b=0, u_s_eq(i_M), D_eq);
endfor

## integral values for the section x_abs = 0 .. L_x_eq
L_x_eq  = 1e-3*[48:0.1:80]; # m
L_x_nd_eq = nd_x_plate (Re_eq, Sc_eq, L_x_eq', h_eq);
Sh_h_x_nd_eq = model_filmflow_laminar_nd_sh (L_x_nd_eq);
beta_eq = model_filmflow_laminar_beta (u_s_eq, D_eq, L_x_eq');
Sh_h_eq = model_filmflow_laminar_sh (u_s_eq, D_eq, L_x_eq', h_eq);
##  Sh_h_eq = nd_sh (beta_eq, h_eq, D_eq);
Sh_L_eq = nd_sh (beta_eq, L_x_eq', D_eq);
Pe_h_eq = nd_pe (h_eq, u_avg_eq, D_eq);
Pe_L_eq = nd_pe (L_x_eq', u_avg_eq, D_eq);
delta_c_x_eq = model_filmflow_laminar_delta_x (x_eq, D_eq, u_s_eq);
beta_x_eq = model_filmflow_laminar_beta_x (x_eq, D_eq, u_s_eq);

## integral mass transfer nondimensional equation
x_nd_out = 0.00001:0.0001:1;
Sh_nd_out = model_filmflow_laminar_nd_sh (x_nd_out);

