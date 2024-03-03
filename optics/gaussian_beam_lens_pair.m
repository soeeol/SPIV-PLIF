##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## Gaussian beam calculation for the basic design of a beam expander and a
## light sheet optic
##
## Author: Sören J. Gerke 421016+git@mailbox.org
##

## bp  .. beam parameter product
## div .. beam divergence full angle in rad
## f   .. focal length in m
## msq .. "M squared" beam quality factor in -
## n   .. refractive index of free space, air in -
## r   .. beam radius in m
## rw  .. beam waist radius in m
## wl  .. wave length in m
## zr  .. Rayleigh range in m
## z   .. coordinate beam axis in m
## zo  .. z position of optical element in m

##
## functions
##

## beam parameter
func_bp = @ (z, zr) complex (z, zr);

## divergence in rad
function div = func_div (msq, wl, rw, n)
   div = 2 * msq * wl / (pi * n * rw);
endfunction
func_gb_div = @ (gb, n) func_div (gb.msq, gb.wl, gb.rw, n);

## beam radius (1/e^2) along z
function r = func_r (z, rw, zr)
   r = rw .* sqrt (1 + (z ./ zr).^2);
endfunction

## beam waist radius
function rw = func_rw (msq, wl, zr, n)
  rw = sqrt (zr * msq * wl / ( n * pi));
endfunction
func_gb_rw = @ (gb, n) func_rw (gb.msq, gb.wl, gb.zr, n);

## Rayleigh range
func_zr = @ (msq, wl, rw, n) (pi * rw^2 * n) / (msq * wl);
func_gb_zr = @ (gb, n) func_zr (gb.msq, gb.wl, gb.rw, n);

function rtm = get_rtm (element, p)
## ray transfer matrix definitions
  switch element
    case {"thin lens"}
      rtm = [[1 0];
      			 [-1/p(1) 1]];
    case {"free space"}
      rtm = [[1 p(1)];
             [0 1]];
    otherwise
      error ([element " ... unknown"]);
    endswitch
endfunction

function gbo = beam_propagate (optics, gbi, n)
  rtm_g = eye (2);
  for i = 1 : numel (optics)
    rtm_g = optics{i} * rtm_g;
    rtm{i} = rtm_g;
    ## beam parameter product propagation through i RTMs
    gbo.bp{i} = (rtm{i}(1,1)*gbi.bp + rtm{i}(1,2)) / (rtm{i}(2,1)*gbi.bp + rtm{i}(2,2));
    gbo.zr{i} = imag (gbo.bp{i});
    gbo.zp{i} = real (gbo.bp{i});
    gbo.msq{i} = gbi.msq;
    gbo.wl{i} = gbi.wl;
    gbo.rw{i} = func_rw (gbo.msq{i}, gbo.wl{i}, gbo.zr{i}, n);
    gbo.div{i} = func_div (gbo.msq{i}, gbo.wl{i}, gbo.rw{i}, n);
  endfor
endfunction

function [optics, zo] = optics_thin_lens_pair (f_1, f_2, dz_in, dz_lp, dz_out, z_0)
  dz(1) = dz_in;
  zo(1) = z_0 + dz(1);
  optics{1} = get_rtm ("free space", dz(1));
  dz(2) = 0.0;
  zo(2) = zo(1) + dz(2);
  optics{2} = get_rtm ("thin lens", f_1);
  dz(3) = dz_lp;
  zo(3) = zo(2) + dz(3);
  optics{3} = get_rtm ("free space", dz(3));
  dz(4) = 0.0;
  zo(4) = zo(3) + dz(4);
  optics{4} = get_rtm ("thin lens", f_2);
  dz(5) = dz_out;
  zo(5) = zo(4) + dz(5);
  optics{5} = get_rtm ("free space", dz(5));
endfunction

function [r, z] = beam_radius_propagation (gb, zo, z_0, npoints)
  zz = linspace (z_0, zo(end), npoints);
  z_sec{1} = zz((zz>=z_0)&(zz<=zo(1)));
  r_sec{1} = func_r (z_sec{1}-zo(1)+gb.zp{1}, gb.rw{1}, gb.zr{1});
  for i = 2 : numel (zo)
    z_sec{i} = zz((zz>=zo(i-1))&(zz<=zo(i)));
    r_sec{i} = func_r (z_sec{i}-zo(i)+gb.zp{i}, gb.rw{i}, gb.zr{i});
  endfor
  z = cell2mat (z_sec);
  r = cell2mat (r_sec);
endfunction

function fh = plot_beam (r, z, gb)
  fh = figure (); hold on; grid on;
  plot (z, +1e3 * r, "-b");
  ylabel ("beam radius in mm");
  xlabel ("z in m");
  title (["M^2 = " num2str(gb.msq{1}) ", wl = " num2str(1e9*gb.wl{1}) " nm"]);
endfunction

##
## calculation
##
close all

## input parameters
n = 1; # air
z_0 = 0; # z start position
npoints = 1000; # z resolution in points

## PLIF laser parameters
## model CNI, MBL-III-473-50-3-LED, parameters from protocol:
## major axis
gb_1.wl = 473e-9;
gb_1.msq = 1.07;
gb_1.rw = 0.5 * 955e-6;
gb_1.zr = func_gb_zr (gb_1, n);
gb_1.bp = func_bp (z_0, gb_1.zr);
gb_1.div = func_gb_div (gb_1, n);
gb_1
## minor axis
gb_2.wl = 473e-9;
gb_2.msq = 1.11;
gb_2.rw = 0.5 * 818e-6;
gb_2.zr = func_gb_zr (gb_2, n);
gb_2.bp = func_bp (z_0, gb_2.zr);
gb_2.div = func_gb_div (gb_2, n);
gb_2

## basic Galilean beam expander
dz_in = 0.5;
f_1 = -0.030; # Thorlabs LC1060-A
f_2 = +0.125; # Thorlabs LA1986-A
dz_lp = f_1 + f_2;
dz_out = 1.0;
[o_be, zo] = optics_thin_lens_pair (f_1, f_2, dz_in, dz_lp, dz_out, z_0);
gbo = beam_propagate (o_be, gb_2, n)
[r, z] = beam_radius_propagation (gbo, zo, z_0, npoints);
fh_1 = plot_beam (r, z, gbo);
sin (cell2mat(gbo.div)) ./ cell2mat(gbo.div)
write_series_csv ([pdir.result "r_z_beam-expander"], [z' r'], {"z in m", "r in m"}, [])

## basic light sheet optics
dz_in = 0.5;
f_1 = -0.1; # Thorlabs LC1120-A, 25.4 mm
f_2_cyl_1 = -0.0254; # Thorlabs LK1900L1-A, 16 x 18 mm
f_2_cyl_2 = +0.150; # Thorlabs LJ1629L1-A, 30 x 32 mm
dz_out = 1.0;

## y
dz_lp = 0.135;
[o_be, zo] = optics_thin_lens_pair (f_1, f_2_cyl_2, dz_in, dz_lp, dz_out, z_0);
gbo = beam_propagate (o_be, gb_2, n)
[r, z] = beam_radius_propagation (gbo, zo, z_0, npoints);
fh_22 = plot_beam (r, z, gbo);
sin (cell2mat(gbo.div)) ./ cell2mat(gbo.div)
write_series_csv ([pdir.result "r_z_light-sheet_y"], [z' r'], {"z in m", "r in m"}, [])

## x
dz_lp = 0.125;
[o_be, zo] = optics_thin_lens_pair (f_1, f_2_cyl_1, dz_in, dz_lp, dz_out, z_0);
gbo = beam_propagate (o_be, gb_2, n)
[r, z] = beam_radius_propagation (gbo, zo, z_0, npoints);
fh_21 = plot_beam (r, z, gbo);
sin (cell2mat(gbo.div)) ./ cell2mat(gbo.div)
write_series_csv ([pdir.result "r_z_light-sheet_x"], [z' r'], {"z in m", "r in m"}, [])

