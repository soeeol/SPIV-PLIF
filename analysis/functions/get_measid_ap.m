##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## build measurement identifier
##
## Author: Sören J. Gerke
##

function measid = get_measid_ap (ap)

  id_A = ap.ids_A(ap.i_A);
  id_C = ap.ids_C{ap.i_C};
  id_G = ap.ids_G(ap.i_G);
  id_L = ap.ids_L{ap.i_L};
  id_M = ap.ids_M(ap.i_M);
  id_O = ap.ids_O{ap.i_O};
  id_T = ap.ids_T(ap.i_T);
  id_X = ap.ids_X(ap.i_X);
  id_Z = ap.ids_Z(ap.i_Z);

  measid = [toupper(id_C) "_" (id_O) "_A" num2str(id_A,"%02d") "_T" ...
            num2str(id_T,"%02d") "_" id_L "_M" num2str(id_M,"%02d") "_G" ...
            num2str(id_G,"%03d") "_X" num2str(id_X,"%+03d") "_Z" num2str(id_Z,"%+03d")];

endfunction
