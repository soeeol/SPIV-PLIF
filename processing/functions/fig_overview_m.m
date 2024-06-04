##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## section processing result overview: vel. field, film height and max. velocity
##
## Author: Sören J. Gerke
##

function fh = fig_overview_m (pp, pdir, c_msh, c_dat, c_h, c_masks, u_msh, u_dat, u_h, u_masks, massflow, HU)

  [Lmax, ~] = max (u_dat{4} .* u_masks.gas .* u_masks.wall, [], 1); # TODO: this is generally not surface velocity
  mf_avg = median (massflow(!(massflow==0)));

  ## nusselt_film
  [~, ~, rho, eta, ~, ~] = get_fp_lm (pdir, pp.liquid.data, pp.T.data+273.15);
  mfr = 1.0 / 3600 * mf_avg; # kg / s
  nu = eta / rho;
  alpha = deg2rad (pp.alpha.data);
  re_l = nd_re_inlet (mfr, 50*1e-3, eta);

  [~, umax_N] = model_filmflow_laminar_u_profile_p (nu, alpha, re_l);

  skip = 2*6; # vectors

  fh = figure ();

  subplot (2, 1, 1);
  hold on;
  title ([pp.measid.data " \n" "Mavg = " num2str(mf_avg) "kg/h --- setpoint: M = " num2str(pp.M.data) " kg/h \n liquid in domain: " num2str(HU) " mm^2"]);

  plot (c_msh{1}(1,:), c_h.gas, "r-");
  plot (c_msh{1}(1,:), c_h.wall, "k-");

  axis image;

  quiver3 (u_msh{1}(1:skip:end), u_msh{2}(1:skip:end), 1+u_msh{3}(1:skip:end),
            u_dat{1}(1:skip:end) .* u_masks.gas(1:skip:end) .* u_masks.wall(1:skip:end),
            u_dat{2}(1:skip:end) .* u_masks.gas(1:skip:end) .* u_masks.wall(1:skip:end),
            u_dat{3}(1:skip:end) .* u_masks.gas(1:skip:end) .* u_masks.wall(1:skip:end),
            scale=1.0, color="k");

  ah = annotation ("textbox", [0.75 0.85 2 1]);
  ylabel ("y in mm");
  grid on;
  hold off;

  subplot (2, 1, 2);
  hold on;
  [ax, h1, h2] = plotyy (u_msh{1}(1,:), massflow / pp.M.data, u_msh{1}(1,:), Lmax/umax_N);
  grid on;
  xlabel ("x in mm");
  ylabel (ax(1), "norm. mass flow rate");
  ylabel (ax(2), "norm. surface velocity \n u_surf / u_Nusselt_surf", "interpreter", "none");
  hold off;

endfunction
