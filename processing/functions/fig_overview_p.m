##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## section processing result overview: some fluorescence and velocity profiles
##
## Author: Sören J. Gerke
##

function fh = fig_overview_p (pp, pdir, c_msh, c_dat, c_masks, u_msh, u_dat, u_masks, massflow)
  fp = fig_param ({"combo", "double", "elsevier"});
  domain = get_domain(pp);
  ##xpos = [domain.xmin, domain.xmin+1, 0, domain.xmax-1, domain.xmax];
  xpos = [domain.xmin, -0.5, 0, 0.5, domain.xmax];
  ##xpos = [0];
  mf_avg = median (massflow(!(massflow==0)));
  fig_size_x = fp.xsize; # cm
  fig_size_y = 8;
  fh = figure ("DefaultAxesFontSize", fp.fs_min, "DefaultAxesFontName", fp.fonts{1});
  set (gcf, "PaperPositionMode", "manual");
  set (gcf, "PaperUnits", "centimeters", "PaperPosition", [0 0 fig_size_x fig_size_y]);
  set (gcf, "PaperSize", [fig_size_x fig_size_y]);
  set (gcf, "renderer", "opengl");
  ## Ic profiles
  subplot (1, 2, 1)
  title ("Ic profiles")
  hold on
  xlabel ("y in mm"); ylabel ("Ic in a.u.");
  styles = {"r-", "m-", "c-", "g-", "b-"};
  for i = 1:numel(xpos)
    [ ~, idx] = min (abs ( - xpos(i) + c_msh{1}(1,:)));
    plot (c_msh{2}(:,idx),
            c_dat{1}(:,idx).*c_masks.wall(:,idx).*c_masks.gas(:,idx),
            [styles{i} "x"],"LineWidth",1)
    l_str{i} = ["x = " num2str(xpos(i)) " mm"];
    hold on
  endfor
  for i = 1:numel(xpos)
    [ ~, idx] = min (abs ( - xpos(i) + c_msh{1}(1,:)));
    for j = 1:numel(c_dat)
      plot (c_msh{2}(:,idx), c_dat{1}(:,idx), "k-", "LineWidth", 0.25)
    endfor
  endfor
  legend (l_str, "location", "southeast");
  grid on;
  hold off
  ## u_mag profiles
  if (!isempty(u_msh))
    subplot (1, 2, 2)
    title ("velocity profiles")
    hold on
    xlabel ("y in mm"); ylabel ("um / u_Nusselt_surf", "interpreter", "none");
    ## nusselt_film
    [~, ~, rho, eta, ~, ~] = get_fp_lm (pdir, pp.liquid.data, pp.T.data+273.15);
    mfr = 1.0 / 3600 * mf_avg; # kg / s
    nu = eta / rho;
    alpha = deg2rad (pp.alpha.data);
    re_l = nd_re_inlet (mfr, 50*1e-3, eta)
    [~, umax_N] = model_filmflow_laminar_u_profile_p (nu, alpha, re_l);
    for i = 1:numel(xpos)
      [ ~, idx_u] = min (abs ( - xpos(i) + u_msh{1}(1,:)));
      plot (u_msh{2}(:,idx_u), u_dat{4}(:,idx_u)/umax_N, "k-.", "LineWidth", 0.25)
      plot (u_msh{2}(:,idx_u),
            u_dat{4}(:,idx_u)/umax_N.*u_masks.wall(:,idx_u).*u_masks.gas(:,idx_u),
            [styles{i} "x"], "LineWidth", 1)
      idxyyy = (u_masks.wall(:,idx_u).*u_masks.gas(:,idx_u))==1;
      idxyyy(1:int32(sum(idxyyy))/2)=false;
      p = polyfit (u_msh{2}(idxyyy,idx_u), u_dat{4}(idxyyy,idx_u)/umax_N, 2);
      x = linspace (-0.1, max (u_msh{2}(:,idx_u)), 50);
      y = polyval (p, x);
      plot (x, y, styles{i}, "LineWidth", 0.25)
    endfor
    ylim ([0 1.5])
    grid on;
    ah = annotation ("textbox", [0.75 0.75 2 1]);
    set (ah, "string", ["colored lines:\n parabola fit"]);
    set (ah, "fontsize", fp.fs_min);
  endif
  hold off
endfunction
