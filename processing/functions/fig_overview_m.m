##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## section processing result overview: vel. field, film height and max. velocity
##
## Author: Sören J. Gerke
##

function fh = fig_overview_m (pp, pdir, c_msh, c_dat, c_h, c_masks, u_msh, u_dat, u_h, u_masks, massflow, HU)
  fp = fig_param ({"combo", "double", "elsevier"});
  domain = get_domain (pp);
  [Lmax, ~] = max (u_dat{4}.*u_masks.gas.*u_masks.wall, [], 1);
  mf_avg = median (massflow(!(massflow==0)));
  ## nusselt_film
  [~, ~, rho, eta, ~, ~] = get_fp_lm (pdir, pp.liquid.data, pp.T.data+273.15);
  mfr = 1.0 / 3600 * mf_avg; # kg / s
  nu = eta / rho;
  alpha = deg2rad (pp.alpha.data);
  re_l = nd_re_inlet (mfr, 50*1e-3, eta)
  [~, umax_N] = model_filmflow_laminar_u_profile_p (nu, alpha, re_l);
  skip = 1*6; # vectors
  fig_size_x = 14; # cm
  fig_size_y = 9;
  fh = figure ("DefaultAxesFontSize", fp.fs_min, "DefaultAxesFontName", fp.fonts{1});
  set (gcf, "PaperPositionMode", "manual");
  set (gcf, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
  set (gcf, "PaperSize", [fig_size_x fig_size_y]);
  set (gcf, "renderer", "opengl");
  subplot (2, 1, 1); hold on;
  title ([pp.measid.data " \n" "Mavg = " num2str(mf_avg) "kg/h --- setpoint: M = " num2str(pp.M.data) " kg/h \n liquid in domain: " num2str(HU) " mm^2"], "interpreter", "none");
  plot (c_msh{1}(1,:), c_h.gas,"r-","LineWidth",2);
  plot (c_msh{1}(1,:), c_h.wall,"k-","LineWidth",2)
  ##plot (c_msh{1}(1,:), c_h.gas - c_h.wall,"r-","LineWidth",2)
  axis image
  ##surf(u_msh{1},u_msh{2},u_msh{3},u_dat{4}); shading flat;
  quiver3 (u_msh{1}(1:skip:end),u_msh{2}(1:skip:end),1+u_msh{3}(1:skip:end),
          u_dat{1}(1:skip:end).*u_masks.gas(1:skip:end).*u_masks.wall(1:skip:end),
          u_dat{2}(1:skip:end).*u_masks.gas(1:skip:end).*u_masks.wall(1:skip:end),
          u_dat{3}(1:skip:end).*u_masks.gas(1:skip:end).*u_masks.wall(1:skip:end),
          scale=1.1, color='k');
  ah = annotation ("textbox", [0.75 0.85 2 1]);
  set (ah, "string", ["showing every\n" num2str(skip) "th vector"]);
  set (ah, "fontsize", fp.fs_min);
  set (ah, "linestyle", "none");
  xlim ([domain.xmin domain.xmax]);
  ylabel ("y in mm");
  grid on
  hold off
  subplot (2, 1, 2); hold on;
  [ax, h1, h2] = plotyy (u_msh{1}(1,:), massflow / pp.M.data, u_msh{1}(1,:), Lmax/umax_N);
  set ([h1 h2],"linewidth",1);
  grid on
  xlim ([domain.xmin domain.xmax]);
  ylim (ax(2),[min(Lmax/umax_N) max(Lmax/umax_N)]);
  xlabel ("x in mm");
  ylabel (ax(1),"norm. mass flow rate");
  ylabel (ax(2),"norm. surface velocity \n u_surf / u_Nusselt_surf", "interpreter", "none");
  hold off
endfunction
