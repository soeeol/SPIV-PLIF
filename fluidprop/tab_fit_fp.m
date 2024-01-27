##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## preparation of data tables for get_tp_tab.m function
##
## Author: Sören J. Gerke
##

load ("-v7", [pdir.fptab "fluidprop.v7"], "fp");
mmds = get_fp_dataset (fp, {"water"}, {"molar mass"}, []);
mm_W = mmds.data{1};
mmds = get_fp_dataset (fp, {"propylene glycol"}, {"molar mass"}, []);
mm_PD = mmds.data{1};

##
function plot_tab (pdir, ds, ds_str, T_ds, mf_ds, fp_ds, T_ip, mf_ip, fp_ip)

  if !(isempty(ds_str))
    fh1 = figure ();
##    pause (0.2); set (fh1, "Position", [1 1 1000 600])
    for i = 1:numel(ds.data)
      hax(i) = subplot (numel (ds.data), 1, i);
      plot (ds.data{i}, "-x")
      ylabel (ds.unit{i})
      xlabel ("idx")
    endfor
    title (hax(1), [ds_str "_dataset"], "interpreter", "none")
    save_fn = [pdir.fptab ds_str "_dataset"];
    print (fh1, "-dpng", save_fn)
    close (fh1)

    if numel(ds.unit)>2
      fh2 = figure ();
##      pause (0.2); set (fh2, "Position", [1 1 1000 600])
      surf (T_ip, mf_ip, fp_ip, "facecolor", "none");
      xlabel ("T in K")
      ylabel ("mass fraction in g / g")
      zlabel (ds.unit{3})
      if !isempty(fp_ds) & !isempty(T_ds)
        hold on
        [T_dat, mf_dat] = meshgrid (T_ds, mf_ds);
        plot3 (T_dat, mf_dat, fp_ds, "b*")
      endif
      title ([ds_str "_fit"],"interpreter", "none")
      legend ("fit", "data")
      save_fn = [pdir.fptab ds_str "_fit"];
      print (fh2, "-dpng", save_fn)
      close (fh2)
    endif

  endif
endfunction

##
fnames = {"glycerol-water", "propylene glycol-water", "water"};
pnames = {"rho", "eta", "n"};
for j = 1:numel(fnames)
  for k = 1:numel(pnames)
    ds_str = [];
    T_dat = mf_dat = T_ip = mf_ip = fp_ip = T_fit = mf_fit = fp_ds = T_ds = mf_ds = [];
    ## select datasets for fit range
    switch ([fnames{j} "_" pnames{k}])
      case "glycerol-water_rho"
        ds = get_fp_dataset (fp, fnames(j), pnames(k), {"VolkA2018"});
        ds_str = [ds.source{1}{1} "_" fnames{j} "_" pnames{k}];
        T_ds = ds.data{1};
        mf_ds = ds.data{2};
        fp_ds = []; # model
        T_ip = 293.15:0.5:303.15;
        mf_ip = 0:0.0125:1;
        [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
        fp_ip = rho_ip = rho_PT_W_model (mf_ip, T_ip);
        save_fn = [pdir.fptab ds_str "_tab.v7"];
        save ("-v7", save_fn, "rho_ip", "T_ip", "mf_ip")
      case "propylene glycol-water_rho"
        ds = get_fp_dataset (fp, fnames(j), pnames(k), {"GeorgeJ2003"});
        ds_str = [ds.source{1}{1} "_" fnames{j} "_" pnames{k}];
        T_ds = ds.data{1};
        mf_ds = fp_mf_mx (ds.data{2}, mm_PD, mm_W);
        fp_ds = ds.data{3}';
        T_ip = 298.15:0.5:308.15;
        mf_ip = 0:0.0125:1;
        [T_fit, mf_fit] = meshgrid (T_ds, mf_ds);
        [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
        fp_ip = rho_ip = interp_dat_ext (T_fit, mf_fit, fp_ds, T_ip, mf_ip);
        save_fn = [pdir.fptab ds_str "_tab.v7"];
        save ("-v7", save_fn, "rho_ip", "T_ip", "mf_ip")
      case "water_eta"
        ds = get_fp_dataset (fp, fnames(j), pnames(k), {"IAPWS 2008"});
        ds_str = [ds.source{1}{1} "_" fnames{j} "_" pnames{k}];
        T_ds = ds.data{1};
        mf_ds = 1;
        fp_ds = ds.data{2}';
        T_ip = 275.15:0.5:353.15;
        mf_ip = 1;
        [T_fit, mf_fit] = meshgrid (T_ds, mf_ds);
        [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
        fp_ip = eta_ip = interp1 (T_fit', fp_ds, T_ip, "pchip");
        save_fn = [pdir.fptab ds_str "_tab.v7"];
        save ("-v7", save_fn, "eta_ip", "T_ip", "mf_ip")
      case "glycerol-water_eta"
        ds = get_fp_dataset (fp, fnames(j), pnames(k), {"SegurJ1951"});
        ds_str = [ds.source{1}{1} "_" fnames{j} "_" pnames{k}];
        T_ds = ds.data{1};
        mf_ds = ds.data{2};
        fp_ds = ds.data{3}';
        T_ip = 293.15:0.5:303.15;
        mf_ip = 0:0.0125:1;
        [T_fit, mf_fit] = meshgrid (T_ds, mf_ds);
        [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
        fp_ip = eta_ip = interp_dat_int (T_fit, mf_fit, fp_ds, T_ip, mf_ip);
        save_fn = [pdir.fptab ds_str "_tab.v7"];
        save ("-v7", save_fn, "eta_ip", "T_ip", "mf_ip")
      case "propylene glycol-water_eta"
        ## KhattabIS2017
        ds = get_fp_dataset (fp, fnames(j), pnames(k), {"KhattabIS2017"});
        ds_str = [ds.source{1}{1} "_" fnames{j} "_" pnames{k}];
        T_ds = ds.data{1};
        mf_ds = fp_mf_mx (ds.data{2}, mm_PD, mm_W);
        fp_ds = ds.data{3}';
        T_ip = 293:0.5:303;
        mf_ip = 0:0.0125:1;
        [T_fit, mf_fit] = meshgrid (T_ds, mf_ds);
        [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
        fp_ip = eta_ip = interp_dat_int (T_fit, mf_fit, fp_ds, T_ip, mf_ip);
        save_fn = [pdir.fptab ds_str "_tab.v7"];
        save ("-v7", save_fn, "eta_ip", "T_ip", "mf_ip")
        plot_tab (pdir, ds, ds_str, T_ds, mf_ds, fp_ds, T_ip, mf_ip, fp_ip);
        ## GeorgeJ2003
        ds = get_fp_dataset (fp, fnames(j), pnames(k), {"GeorgeJ2003"});
        ds_str = [ds.source{1}{1} "_" fnames{j} "_" pnames{k}];
        T_ds = ds.data{1};
        mf_ds = fp_mf_mx (ds.data{2}, mm_PD, mm_W);
        fp_ds = ds.data{3}';
        T_ip = 298.15:0.5:308.15;
        mf_ip = 0:0.0125:1;
        [T_fit, mf_fit] = meshgrid (T_ds, mf_ds);
        [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
        fp_ip = eta_ip = interp_dat_int (T_fit, mf_fit, fp_ds, T_ip, mf_ip);
        save_fn = [pdir.fptab ds_str "_tab.v7"];
        save ("-v7", save_fn, "eta_ip", "T_ip", "mf_ip")
        plot_tab (pdir, ds, ds_str, T_ds, mf_ds, fp_ds, T_ip, mf_ip, fp_ip)
        ## SunT2004 - combine data
        ## pure water
        ds = get_fp_dataset (fp, {"water"}, {"eta"}, {});
        T_ip = [297 313];
        eta_mf0 = interp1 (ds.data{1}, ds.data{2}, T_ip, "spline");
        ## pure glycol
        ds = get_fp_dataset (fp, {"propylene glycol"}, {"eta"}, {"SunT2004"});
        eta_mf1 =  interp1 (ds.data{1}, ds.data{2}, T_ip, "spline", "extrap");
        ##
        ds = get_fp_dataset (fp, fnames(j), pnames(k), {"SunT2004"});
        ds_str = [ds.source{1}{1} "_" fnames{j} "_" pnames{k}];
        T_ds = mean (ds.data{1}');
        mf_ds = fp_mf_mx (ds.data{2}, mm_PD, mm_W);
        mf_ds = [0 mf_ds 1]';
        fp_ds =  [eta_mf0; ds.data{3}'; eta_mf1];
        mf_ip = 0:0.0125:1;
        T_ip = 297:0.5:313;
        [T_fit, mf_fit] = meshgrid (T_ds, mf_ds);
        [T_ip, mf_ip] = meshgrid (T_ip, mf_ip);
        fp_ip = eta_ip = interp_dat_int (T_fit, mf_fit, fp_ds, T_ip, mf_ip);
        save_fn = [pdir.fptab ds_str "_tab.v7"];
        save ("-v7", save_fn, "eta_ip", "T_ip", "mf_ip")
      otherwise
    endswitch
    plot_tab (pdir, ds, ds_str, T_ds, mf_ds, fp_ds, T_ip, mf_ip, fp_ip);
  endfor
endfor

