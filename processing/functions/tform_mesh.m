##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## translate mesh
##
## Author: Sören J. Gerke
##

function [mesh_out] = tform_mesh (mesh_in, pp, id, param)
  if isempty (param)
    param = getfield (pp, id).data;
  endif
  mesh_out = mesh_in;
  if abs (param) > 0
    switch (id)
      case {"rot_c", "rot_u"}
        [mesh_out{1}, mesh_out{2}]  = mesh_rotate (mesh_in{1}, mesh_in{2}, param);
      case {"xoff_c", "xoff_c0", "xoff_c1", "xoff_u"}
        mesh_out{1} = mesh_in{1} + param;
      case {"yoff_c_ini", "yoff_c", "yoff_c0", "yoff_c1", "yoff_u"}
        mesh_out{2} = mesh_in{2} + param;
      otherwise
        error (["no transform for" id]);
    endswitch
  elseif param == 0
    ;
  else
    error (["false param: " num2str(param)]);
  endif
endfunction
