##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fit data tables and plot data overview for the refractive index matching
## experiment the two binary mixtures
##
## Author: Sören J. Gerke
##

tables_dir = pdir.fptab;
load ("-v7", [tables_dir "fluidprop.v7"], "fp");
fnames = {"PT_W", "PD_W"};
pnames = {"eta", "n"};
for j = 1:numel(fnames)
  for k = 1:numel(pnames)
    clear rho_ip eta_fit eta_ip T_dat mf_dat fp_fit fp_ip p_fit
    ds = get_fp_dataset (fp, fnames(j), pnames(k), {"GerkeSJexp"});
    ds_str = [ds.source{1}{1} "_" ds.fluid{1} "_" ds.prop{1}{1}];
    fh1 = figure ();
##    pause (0.2); set (fh1, "Position", [1 1 1000 600])
    for i = 1:numel(ds.data)
      hax(i) = subplot (numel(ds.data), 1, i);
      plot (ds.data{i},"-x")
      ylabel (ds.unit{i})
      xlabel ("idx")
    endfor
    title (hax(1), [ds_str "_dataset"], "interpreter", "none")
    save_fn = [tables_dir ds_str "_dataset"];
    print (fh1, "-dpng", save_fn)
    close (fh1)
    mf_ip = min(ds.data{2}):0.0125:max(ds.data{2});
    T_ip = min(ds.data{1}):0.5:max(ds.data{1});
    for i = 1:numel(ds.data{1})
      p_fit{i} = polyfit (ds.data{2}, ds.data{3}(i,:), 4);
      fp_fit(i,:) = polyval (p_fit{i}, mf_ip);
    endfor
    [T_dat, mf_dat] = meshgrid (ds.data{1}, ds.data{2});
    [T_fit, mf_fit] = meshgrid (ds.data{1}, mf_ip);
    [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
    fp_ip = interp_dat_int (T_fit, mf_fit, fp_fit', T_ip, mf_ip);
    save_fn = [tables_dir ds_str "_tab.v7"];
    save ("-v7", save_fn, "fp_ip", "T_ip", "mf_ip")
##
    fh2 = figure ();
##    pause (0.2); set (fh2, "Position", [1 1 1000 600])
    surf (T_ip, mf_ip, fp_ip, "facecolor", "none");
    xlabel (ds.unit{1})
    ylabel (ds.unit{2})
    zlabel (ds.unit{3})
    hold on
    plot3 (T_dat, mf_dat, ds.data{3}', "b*")
    title ([ds_str "_fit"], "interpreter", "none")
    legend ("fit", "data")
    save_fn = [tables_dir ds_str "_fit"];
    print (fh2, "-dpng", save_fn)
    close (fh2)
  endfor
endfor

