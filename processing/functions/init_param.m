##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## initialize processing parameter structure
##
## Author: Sören J. Gerke
##

function [param] = init_param ()
  fields = {
    "rot_c"; "rot_u";
    "xoff_c"; "xoff_c0"; "xoff_c1"; "xoff_u";
    "yoff_c_ini"; "yoff_c"; "yoff_c0"; "yoff_c1";
    "yoff_u"; "y0_if_c";
    "tol_wall"; "tol_if";
    "measid"; "recid_Ic"; "recid_Ic0"; "recid_Ic1"; "recid_u";
    "valid";
    "optset"; "cell"; "liquid";
    "M"; "G"; "T";
    "X"; "Z"; "alpha";
    "type"; "savedate"
  };
  ids = {
    "dtheta_c in rad"; "dtheta_u in rad";
    "xoff_c in mm"; "xoff_c0 in mm"; "xoff_c1 in mm"; "xoff_u in mm";
    "yoff_c_ini in mm"; "yoff_c in mm"; "yoff_c0 in mm"; "yoff_c1 in mm";
    "yoff_u in mm"; "y0_if_c in mm";
    "tol_wall_detect"; "tol_if_detect";
    "Meas ID"; "RecordID Ic"; "RecordID Ic0"; "RecordID Ic1"; "RecordID u";
    "checked records?";
    "Optical Setup"; "Cell"; "Liquid Mixture";
    "M in kg/h"; "G in Nl/min"; "T in °C";
    "X* in mm"; "Z* in mm"; "Inclination in °";
    "type of processing"; "last save"
  };
  for i = 1:numel (fields)
    param.(fields{i}).id = ids(i,:);
    param.(fields{i}).data = [];
  endfor
endfunction
