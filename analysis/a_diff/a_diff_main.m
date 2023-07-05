##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## script to analyze the results of the PLIF diffusion measurements
##
## "Interactive" script with following steps:
## (0) init: parameters defining the analysis
## (1) concentration field, profiles normal to gas-liquid interface, diffusion front estimate
## (2) diffusivity estimate and outputs
##
## background: erfc fitted to measured normalized concentration profile
## normal to gas-liquid interface to estimate diffusivity
## c* (sn, a) = erfc (n / a)
## where a = 2 * sqrt (D * t); sn in m, t in s, D in m^2/s
## second derivative is 4 * sn * exp (-sn^2 / a^2) / (sqrt (pi) * a^3)
## max. of second derivative is located at sn_D = a / sqrt(2)
## Fick's diffusivity then is D = sn_D^2 / (2 * t)
##
## The fit function is utilized for a stable position estimate of maximal second
## derivative, the "diffusion front" sn_D.
##
## Author: Sören J. Gerke
##

## (0) init
if 1
  proc_type = "a_2d_diff";
  measid = "DIFF_M13_A15_T25_WG141_M_G002_X_Z"
##  measid = "DIFF_M26_A15_T25_WP141_M_G002_X_Z"

  save_dir_p = [pdir.plot proc_type "/" measid "/"];
  save_dir = [pdir.analyzed proc_type "/" measid "/"];
  mkdir (save_dir_p)
  mkdir (save_dir)
  cd (save_dir)

  testplots = false

  ## c* at max. 2nd. derivative: erfc (n_front / a) =
  cddsnm = erfc (1 / sqrt(2))

  ## available processed measurements
  switch measid
    case {"DIFF_M13_A15_T25_WG141_M_G002_X_Z"}
      ## absorbtion of O2 from air to desorbed liquid WG141
      profile_depth = 1400e-3; # mm
      lims_x = [-1, 5];
      lims_y = [ 0, 3];
      xrange = [1, 2.5];
      ## logged times @ frames recorded
      t_s{1} = [0]; # in s
      t_s{2} = [1:159];
      t_s{3} = t_s{2}(end) + 1*60;
      t_s{4} = t_s{3}(end) + 10*60;
      t_s{5} = t_s{4}(end) + 13*60;
      t_s{6} = 9999;
    case {"DIFF_M26_A15_T25_WP141_M_G002_X_Z"}
      ## absorbtion of O2 from air to desorbed liquid WP141
      profile_depth = 1000e-3; # mm
      lims_x = [-2,+2];
      lims_y = [ 0, 3];
      xrange = [-1.7, 1.7];
      ## logged times @ frames recorded
      t_s{1} = [0]; # s ... desorbed
      t_s{2} = [1:79];
      t_s{3} = 260;
      t_s{4} = 545;
      t_s{5} = 955;
      t_s{6} = 1425;
      t_s{7} = 9999;
  endswitch
endif

## (1) cn, cp, sn_D
if 0
  ## load processed data (aligned, interpolated on common grid)
  ldirs = glob ([pdir.processed measid "/" "*___" proc_type "/"]);
  ldir = ldirs{end} # newest processing result
  [msh, c_dat, ~, h_g, ~, ~, ~, ~, pp] = load_uc (ldir, pdir.work);
  x = msh{1}(1,:);
  y = msh{2}(:,1);

  ## transform quenching signal to normalized concentration
  Ic0 = c_dat{1}{1};
  Ic1 = c_dat{end}{end}; ## last image: of the oxygen saturated
  cn = cell ();
  k = 0;
  for i = 1:numel (t_s)
    for j = 1:numel (c_dat{i})
      k++;
      [cn{k}] = calc_cn ({c_dat{i}{j}, Ic0, Ic1}, [0, 1], "linear", 3, false);
    endfor
  endfor

  ## time vector
  tn = [];
  k = 0;
  for i = 1:numel (t_s)
    for j = 1:1:numel (c_dat{i})
      k++;
      tn(k) = t_s{i}(j);
    endfor
  endfor

  if testplots
    fh = figure ();
    for i = 1:numel (cn)
      clf
      surf (msh{1}, msh{2}, msh{3}, cn{i})
      hold on
      grid off
      shading flat
      caxis ([0 1])
      colormap viridis
      axis image
      xlabel ("x in mm")
      ylabel ("y in mm")
      title ([ "cn, t = " num2str(tn(i)) " s"])
      plot3 (x, h_g, ones (1, numel (h_g)), "m", "linewidth", 2)
      view ([0 0 1])
      pause (0.5)
    endfor
  endif

  if testplots
    idx = 800
    figure ()
    hold on
    ylim ([-.1 1])
    xlim ([0.5 3])
    for i = 1:numel (cn)
      plot (y, mean (cn{i}(:,idx-20:idx+20), 2))
      plot (mean(h_g(idx-20:idx+20))*[1 1], [0 1], "r")
      title ([ "idx = " num2str((i)) "; " "nt = " num2str(tn(i)) " s"])
      pause (0.5)
    endfor
  endif

  ## create mesh of lines normal to gas-liquid interface
  sf_c = get_sf (msh);
  skip = 1;
  sf_p = 1 * sf_c(1);
  px = x';
  py = h_g';
  ap = angle2Points ([px(1:skip:end-1) py(1:skip:end-1)], [px(2:skip:end) py(2:skip:end)]);
  clines = createLine ([px(1:skip:end) py(1:skip:end)], pi/2 + [ap(1:skip:end); ap(end)] .* ones (numel (px(1:skip:end)),1)); # line from angle and position
  nt = - linspace (0, profile_depth, round(profile_depth/sf_p+1));
  ## concentration profile interface coordinate
  snp = sf_p * (0:numel (nt) - 1) * 1e-3; # m
  for i = 1:length (clines)
    clines_edges{i} = createEdge (clines(i,:) .* ones (numel (nt),1), nt');
  endfor
  msh_n{1} = msh_n{2} = zeros (length(clines), numel (nt));
  for i = 1:length (clines)
    edge = createEdge (clines(i,:) .* ones (numel (nt),1), nt');
    msh_n{1}(i,:) = edge(:,3);
    msh_n{2}(i,:) = edge(:,4);
  endfor
  msh_n = {msh_n{1}, msh_n{2}, zeros(size(msh_n{1}))};

  ##
  run ("a_diff_fig_chamber.m")

  ## local interace curvature
  [curvC, curvR] = calc_curvature_xy (x, h_g);

  ## interpolate cn values on line meshes
  cn_s = [];
  for k = 1:numel (cn)
    cn_s{k} = interp2 (msh{1}, msh{2}, cn{k}, msh_n{1}, msh_n{2}, "pchip", 0.0);
    cn_s{k}(isnan(cn_s{k})) = 0;
  endfor

  ## c field with some profile lines
  plot_map_msh (msh, imsmooth(cn{12}))
  hold on
  draw_cell (pp.cell.data, 1)
  for i = 1:100:size(cn_s{1},1)
    plot (msh_n{1}(i,:), msh_n{2}(i,:), "m")
  endfor
  axis equal
  plot3 (x, py, ones(1,numel (h_g)), "r", "LineWidth", 1.5)
  title ("normalized concentration from c0 to csat and excemplaric interface normal profiles")
  xlabel ("x in mm")
  ylabel ("y in mm")
  xlim ([min(x) max(x)])

  [~, x_idx_min] = min (abs ( x - xrange(1)));
  [~, x_idx_max] = min (abs ( x - xrange(2)));
  x_idx_r = x_idx_min:x_idx_max;

  cp = [];
  for i = 1:1:numel (cn_s)
    cp(i,:) = median (cn_s{i}(x_idx_r,:), 1);
    cp(i,(cp(i,:)>1)) = 1;
  endfor

  figure ()
  surf (snp*1e3, tn, cp)
  view ([0 0 1])
  shading flat
  colormap viridis
  colorbar
  xlabel ("s in mm")
  ylabel ("t in s")
  xlim ([0 profile_depth]);
  ylim ([0 90]);
  caxis ([0 1])

  if testplots
    figure (); hold on;
    ylim ([0 1]);
    xlim ([0 profile_depth])
    xlabel ("s in mm")
    ylabel ("cn");
    for i = 2:numel (cn_s)
      plot (snp*1e3, cp(i,:))
      title (["idx = " num2str(i) " , nt = " num2str(tn(i)) ]);
      pause (0.5)
    endfor
  endif

  ## bulk and interface estimate for c profile normalization
  cp_b = cp_s = [];
  for i = 1:numel (cn_s)
    cp_b(i) = min (movmean (cp(i,:), 21));
    cp_s(i) = cp(i,1);
  endfor
  ##
  figure ()
  hold on
  plot (tn, cp_b)
  plot (tn, cp_s)
  xlim ([0 160]);
  ylim ([0 1]);
  ##
  cp_n = [];
  for i = 1:numel (tn)
      cp_n(i,:) = (cp(i,:) - cp_b(i)) / (cp_s(i) - cp_b(i));
  endfor

  figure ()
  surf (snp*1e3, tn, cp_n)
  view ([0 0 1])
  shading flat
  colormap viridis
  xlabel ("s in mm")
  ylabel ("t in s")
  xlim ([0 profile_depth]);
  ylim ([2 90]);
  caxis ([0 1])
  title ("normalized c profile vs. time")

  if testplots
    figure (); hold on;
    ylim ([0 1]); xlim ([0 profile_depth]); xlabel ("s in mm"); ylabel ("cn");
    for i = 2:numel (tn)
      plot (snp*1e3, cp_n(i,:))
      title (["idx = " num2str(i) " , t = " num2str(tn(i)) ]);
      pause (0.5)
    endfor
  endif

  cp_n(isnan(cp_n)) = 0.0;

  ##
  ## estimate D from position of maximum of second derivative of cp
  ##
  c_fitf = @ (p, s) fitfn_cn_diff (p, s);
  settings = optimset ("lbound", [0.01; 1], "ubound", [100; 10], "TolFun", 1e-16);
  init = [2.0 1.0]'
  ##
  p_fit = fitrange = cell (1, numel (tn));
  p_fit(:) = init;
  a_fit = sn_D = zeros (1, numel (tn));
  ##
  cp_nn = cp_n;
  cn0 = ones (1, numel (tn));
  for i = 1:numel (tn)
    p_fit{i} = init;
    spf = splinefit (snp(1:end), cp_n(i,1:end), 20, "order", 3);
    d1s = ppval (ppder(spf), snp(1:end));
    [PKS, LOC, ~] = findpeaks (abs (d1s));
    [~, LOCm] = max (PKS);
    idx_maxd1s = LOC(LOCm);
  ##  fitrange{i} = [idx_maxd1s-1:min(idx_maxd1s+1+tn(i), numel (snp))];
    fitrange{i} = [idx_maxd1s-1:numel(snp)];
    if i == 1
      fitrange{i} = [1:10];
    endif
  ##  try
      [p_fit{i}, ~] = nonlin_curvefit (c_fitf, init, snp(fitrange{i})', (cp_n(i,fitrange{i}))', settings);
      cn0(i) = c_fitf (p_fit{i}, 0);
      cp_nn(i,1:idx_maxd1s+1) = c_fitf (p_fit{i}, snp(1:idx_maxd1s+1));
      cp_nn(i,:) = cp_nn(i,:) / cn0(i);
      a_fit(i) = p_fit{i}(1) * 1e-4;
      sn_D(i) = a_fit(i) / sqrt(2);
  ##  catch
  ##    warning (["fit issue at i = " num2str(i)]);
  ##  end
  endfor

  figure ()
  surf (snp*1e3, tn, cp_nn)
  view ([0 0 1])
  shading flat
  colormap viridis
  xlabel ("s in mm")
  ylabel ("t in s")
  xlim ([0 profile_depth]);
  ylim ([2 90]);
  ##zlim ([0 1])
  caxis ([0 1])
  hold on
  plot3 (sn_D*1e3, tn, 1.1*ones(1,length(tn))*cddsnm, "r", "linewidth", 1)

  ## sn_D from intersection with cddsnm
  [~, idx_D] = min ( abs (cp_nn - cddsnm), [], 2)
  sn_Di = snp(idx_D);

  if testplots
    fh = figure ();
    for i = 4:numel (tn)-5
      subplot (1, 2, 1)
      plot (snp*1e3, cp_nn(i,:), "-", "linewidth", 1)
      hold on
      plot ([1]*sn_D(i)*1e3, [cddsnm], "b*")
      xlim ([0 0.4]);
      ylim ([0 1]);
      xlabel ("sn in mm")
      ylabel ("c*")
      box on
      title ("pos. of max. ( d c*(sn) / d^2 sn )", "fontsize", 11);
      ##
      subplot (1, 2, 2)
      title ("= diffusion front vs. time", "fontsize", 11);
      plot (tn(i), sn_D(i)*1e3, "b*")
    ##  plot (tn(i), sn_Di(i)*1e3, "r*")
      hold on
      xlim ([0 90]);
      ylim ([0 0.4]);
      xlabel ("t in s")
      ylabel ("sn_D in mm")
      ax_main = axes ("title", ["t = " num2str(tn(i)) " s"], "fontsize", 11, "visible", "off");
      pause(0.5)
    endfor
  endif

  ##
  cd (save_dir)
  save -v7 "pp.v7" pp
  save -v7 "msh.v7" tn msh x y h_g curvC curvR x_idx_r
  save -v7 "c.v7" cn
  save -v7 "profiles_msh.v7" profile_depth sf_p snp msh_n ap
  save -v7 "profiles_c.v7" cn_s cp cp_s cp_b cp_n
  save -v7 "profiles_fit.v7" p_fit fitrange cn0 cp_nn sn_D sn_Di
endif


## (2) D, out
if 1
  ##
  cd (save_dir)
  load -v7 "pp.v7"
  load -v7 "msh.v7"
  load -v7 "c.v7"
  load -v7 "profiles_msh.v7"
  load -v7 "profiles_c.v7"
  load -v7 "profiles_fit.v7"
  c_fitf = @ (p, s) fitfn_cn_diff (p, s);

  ## diffusion front vs time
  switch (pp.liquid.data)
    case {"WG141"}
      idx_ref = 2;
      t_range_1 = 3:14;
      t_range_2 = 40:160;
    case {"WP141"}
      idx_ref = 2;
      t_range_1 = [3:14];
      t_range_2 = 15:70;
  endswitch
  ## slope of  sn_D(t)^2 / 2 is proportional to diffusivity
  sfsq = (sn_D.^2) / 2; # D = sn_D^2 / (2 * t)
  [D_fit_0, S0] = polyfit (tn(t_range_1), (sn_D(t_range_1).^2)/2, 1);
  tnull = - D_fit_0(2) / D_fit_0(1)
  tvec = tn - tnull;
  ##
  [D_fit_1, S1] = polyfit (tvec(t_range_1), sfsq(t_range_1), 1);
  [D_fit_2, S2] = polyfit (tvec(t_range_2), sfsq(t_range_2), 1);
  D_fit_1
  D_fit_2
  ## standard deviation of polynomial fit coefficients
  std_D1 = abs (sqrt (diag (S1.C)/S1.df) * S1.normr);
  std_D2 = abs (sqrt (diag (S2.C)/S2.df) * S2.normr);
  ## +/- 2 sigma for 95 %
  ## standard deviation of polynomial coefficients in percent
  std_D1_rel = abs (std_D1 ./ D_fit_1') * 100
  std_D2_rel = abs (std_D2 ./ D_fit_2') * 100

  cd (save_dir)
  save -v7 "D.v7" D_fit_1 D_fit_2 std_D1 std_D2 std_D1_rel std_D2_rel

  ##
  fh = figure (); hold on;
  plot (tvec(2:end), sfsq(2:end), "k*")
  plot (tvec, polyval (D_fit_1, tvec), "k-")
  plot ([0,tvec(t_range_1)], polyval (D_fit_1, [0,tvec(t_range_1)]), "-b")
  plot (tvec(t_range_1), polyval (D_fit_1+2*std_D1', tvec(t_range_1)), "-.b")
  plot (tvec(t_range_1), polyval (D_fit_1-2*std_D1', tvec(t_range_1)), "-.b")
  plot (tvec, polyval (D_fit_2, tvec), "k-")
  plot (tvec(t_range_2), polyval (D_fit_2, tvec(t_range_2)), "-r")
  plot (tvec(t_range_2), polyval (D_fit_2+2*std_D2', tvec(t_range_2)), "-.r")
  plot (tvec(t_range_2), polyval (D_fit_2-2*std_D2', tvec(t_range_2)), "-.r")
  ##plot (tvec, polyval ([D_fit_2(1)-2*std_D2(1) D_fit_2(2)], tvec), "-.r")
  switch (pp.liquid.data)
    case {"WG141"}
      xlim([0 tn(t_range_2(end))-tnull]);
      ylim([0 sfsq(t_range_2(end))]);
    case {"WP141"}
      xlim([0 tn(t_range_2(end))-tnull+10]);
      ylim([0 1.25*sfsq(t_range_2(end))]);
  endswitch
  xlabel("t in s")
  ylabel("sn_D^2 / 2 in m^2")

  ##
  ## write data for tikz plots
  ##
  print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [save_dir_p "diffusion-front"]);
  plt_1_h = {mfilename, date; "t in s (offset with tnull)", "S^2 in m^2"};
  plt_1_d = [tvec' sfsq'];
  plt_2_h = {mfilename, date, "", "", "", ""; "Dfit1", "Dfit2", "sigma Dfit1", "sigma Dfit2", "sigma Dfit1 in %", "sigma Dfit2 in %"};
  plt_2_d = [D_fit_1' D_fit_2' std_D1 std_D2 std_D1_rel std_D2_rel];
  plt_3_h = {mfilename, date; "extrapolated time of the start of the diffustion in s", ""};
  plt_3_d = [tnull];
  cd (save_dir_p)
  cell2csv (["tikz_diff-front_h.csv"], plt_1_h)
  csvwrite (["tikz_diff-front_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
  cell2csv (["tikz_D_fit_h.csv"], plt_2_h)
  csvwrite (["tikz_D_fit_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
  cell2csv (["tikz_tnull_h.csv"], plt_3_h)
  csvwrite (["tikz_tnull_d.csv"], plt_3_d, "append", "off", "precision","%.4e")

  ##
  ## cn field for printing
  ##
  S = snp * 1e3; # depth
  X = msh{1}(1,x_idx_r);
  [SS, XX] = meshgrid (S, X);
  switch (pp.liquid.data)
    case {"WG141"}
      [SI, XI] = meshgrid ([0:2e-3:1], [1.9:2e-3:2.1]);
    case {"WP141"}
      [SI, XI] = meshgrid ([0:2e-3:1], [0.4:2e-3:0.6]); # exemplaric section
  endswitch
  ##
  fh = figure ();
  for i = 1:numel (tn)
    cprint = cn_s{i}(x_idx_r,:);
    cprint = (cprint - cp_b(i)) / (cp_s(i) - cp_b(i));
    cprint = interp2(SS,XX,cprint,SI,XI);
    cprint(cprint<0) = 0;
    cprint(cprint>1) = 1;
  ##
    clf
    sph = surf (SI, XI, imsmooth(cprint,1))
    hold on
    view([0 0 1])
    shading flat
    grid off
    colormap viridis
    plot3([1 1] * (sn_D(i)*1e3), [min(XI(:,1)) max(XI(:,1))], [1 1], "r", "linewidth", 1.5)
    plot3([1 1] * (sn_D(1)*1e3), [min(XI(:,1)) max(XI(:,1))], [1 1], "m", "linewidth", 1.5)
    axis image
    title (["t = " num2str(tn(i)) " s; sn_D = " num2str(sn_D(i)*1e3) " mm; sn_D ref = " num2str(sn_D(1)*1e3) "mm"])
    xlabel ("sn in mm"); ylabel ("st in mm")
    print (fh, "-djpeg", "-color", ["-r" num2str(250)], [save_dir_p "fig_" num2str(i,"%03i")]);
    imwrite ((ind2rgb (gray2ind(cprint), colormap("viridis"))), [save_dir_p "tikz_" num2str(tn(i)) ".png"])
  endfor

  cd (save_dir_p)
  system ('magick mogrify -trim "./fig_*.jpg"');
  system ('magick mogrify -resize 1280x408! "./fig_*.jpg"');
  system ('rm out.mp4')
  system ('ffmpeg -framerate 10 -pattern_type glob -i "./fig_*.jpg" -c:v libx264 out.mp4');
  system ('rm out.gif')
  system ('ffmpeg -framerate 10 -pattern_type glob -i "./fig_*.jpg" out.gif');

  fh = figure ();
  subplot (3, 1, 1)
    i = 10
    sph = surf (rgb2gray(imread([save_dir_p "tikz_" num2str(tn(i)) ".png"])))
    shading flat
    view ([0 0 1])
    axis image
    axis off
    hold on
    plot3 ([1 1] * (sn_D(i)*1e3)/2e-3, [0 0.2/2e-3], 2^16*[1 1], "r", "linewidth", 1.5)
    plot3 ([1 1] * (sn_D(2)*1e3)/2e-3, [0 0.2/2e-3], 2^16*[1 1], "m", "linewidth", 1.5)
    title (["t = " num2str(tn(i)) " s; sn_D = " num2str(sn_D(i)*1e3) " mm; sn_D ref = " num2str(sn_D(1)*1e3) "mm"])
  subplot(3, 1, 2)
    i = 37
    sph = surf (rgb2gray(imread([save_dir_p "tikz_" num2str(tn(i)) ".png"])))
    shading flat
    view([0 0 1])
    axis image
    axis off
    hold on
    plot3 ([1 1] * (sn_D(i)*1e3)/2e-3, [0 0.2/2e-3], 2^16*[1 1], "r", "linewidth", 1.5)
    plot3 ([1 1] * (sn_D(2)*1e3)/2e-3, [0 0.2/2e-3], 2^16*[1 1], "m", "linewidth", 1.5)
    title (["t = " num2str(tn(i)) " s; sn_D = " num2str(sn_D(i)*1e3) " mm; sn_D ref = " num2str(sn_D(1)*1e3) "mm"])
  subplot(3, 1, 3)
    i = 80
    sph = surf (rgb2gray(imread([save_dir_p "tikz_" num2str(tn(i)) ".png"])))
    shading flat
    view([0 0 1])
    axis image
    axis off
    hold on
    plot3 ([1 1] * (sn_D(i)*1e3)/2e-3, [0 0.2/2e-3], 2^16*[1 1], "r", "linewidth", 1.5)
    plot3 ([1 1] * (sn_D(2)*1e3)/2e-3, [0 0.2/2e-3], 2^16*[1 1], "m", "linewidth", 1.5)
    title (["t = " num2str(tn(i)) " s; sn_D = " num2str(sn_D(i)*1e3) " mm; sn_D ref = " num2str(sn_D(1)*1e3) "mm"])
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_p "3diff_maps_2tikz"]);

  ## 3 profiles
  switch (pp.liquid.data)
    case {"WG141"}
      i_print = [10,37,79] - round(4.8)
    case {"WP141"}
      i_print = [10,37,79]
  endswitch
  ##
  fh = figure ();
  hold on
  i = i_print(1);
  k = 0;
  styles = {"r", "g", "b"};
  for i = i_print
    k++;
    plot (snp(fitrange{i})*1e3, cp_nn(i,fitrange{i}), [styles{k} "x;t = " num2str(tn(i)) " s;"])
##      plot(snp(:)*1e3, cp_nn(i,:), [styles{k} "-;1;"])
  endfor
  legend ("autoupdate", "off")
  for i = i_print
    plot([1 1] * (sn_D(i)*1e3), [0 c_fitf(p_fit{i}, sn_D(i))/c_fitf(p_fit{i}, snp(1))], "k-.")
    plot([1] * (sn_D(i)*1e3), [cddsnm], "k*")
    plot(snp*1e3, c_fitf(p_fit{i}, snp)/c_fitf(p_fit{i}, snp(1)),"k-", "linewidth", 1)
  endfor
  xlabel("sn in mm")
  ylabel("c*")
  ylim ([0 1])
  print (fh, "-djpeg", "-color", ["-r" num2str(1000)], [save_dir_p "diffusion-front-profiles"]);

  ## for tikz
  i = i_print(1);
  plt_1_h = {mfilename, date, ""; "s in mm", "c* exp. t=9", "c* fit"};
  plt_1_d = [snp'*1e3 cp_nn(i,:)' (c_fitf(p_fit{i},snp)/c_fitf(p_fit{i},snp(1)))'];
  i = i_print(2);
  plt_2_h = {mfilename, date, ""; "s in mm", "c* exp. t=36", "c* fit"};
  plt_2_d = [snp'*1e3 cp_nn(i,:)' (c_fitf(p_fit{i},snp)/c_fitf(p_fit{i},snp(1)))'];
  i = i_print(3);
  plt_3_h = {mfilename, date, ""; "s in mm", "c* exp. t=79", "c* fit"};
  plt_3_d = [snp'*1e3 cp_nn(i,:)' (c_fitf(p_fit{i},snp)/c_fitf(p_fit{i},snp(1)))'];
  cell2csv (["tikz_diff-profile-1_h.csv"], plt_1_h)
  csvwrite (["tikz_diff-profile-1_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
  cell2csv (["tikz_diff-profile-2_h.csv"], plt_2_h)
  csvwrite (["tikz_diff-profile-2_d.csv"], plt_2_d, "append", "off", "precision","%.4e")
  cell2csv (["tikz_diff-profile-3_h.csv"], plt_3_h)
  csvwrite (["tikz_diff-profile-3_d.csv"], plt_3_d, "append", "off", "precision","%.4e")
  plt_4_h = {mfilename, date; "t in s", "sn_D in mm"};
  plt_4_d = [tvec(i_print)' sn_D(i_print)'*1e3];
  cell2csv (["tikz_diff-s_fronts_h.csv"], plt_4_h)
  csvwrite (["tikz_diff-s_fronts_d.csv"], plt_4_d, "append", "off", "precision","%.4e")

  ##
  ## estimate D from direct fit, assuming knowledge of starting time
  ##
  ##  init = [1.0 1.0]'
  ##  c_fit_D_f = @ (p, s, t) p(2) * erfc (s / (2 * sqrt (p(1)*1e-10 * t)));
  ##  for i = 1:numel (tn)
  ##    ...
  ##  endfor
  ##  figure ()
  ##  plot (tn, D_fit, "*")
  ##  hold on
  ##  plot ([0 160], [1 1] .* D_fit_2(1),"r-")
  ##  plot ([0 160], [1 1] .* D_fit_1(1),"r-")
  ##  ylim ([0 1e-9])
  ##  xlim ([0 160])

  ## exposure time
  ## exp time 1: 500 ms  ..Glycerol-W 100 µm / 25 s ... 2 µm per exposure... 1px = 5 µm
  ## exp time 2: 500 ms  ..Propdiol-W 80 µm / 25 s ... 1.6 µm per exposure... 1px = 2 µm
  ## ... good

  RR = outlier_rm (curvR, movmedian(curvR,21));
  figure (); hold on;
  plot (x', RR, "k*")
  title ("interface curvature in mm")
endif
