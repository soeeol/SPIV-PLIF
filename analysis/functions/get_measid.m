##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## build measurement identifier
##
## Author: Sören J. Gerke
##

function measid = get_measid (id_C, id_O, id_A, id_T, id_L, id_M, id_G, id_X, id_Z)

  measid = [toupper(id_C) "_" (id_O) "_A" num2str(id_A,"%02d") "_T" ...
            num2str(id_T,"%02d") "_" id_L "_M" num2str(id_M,"%02d") "_G" ...
            num2str(id_G,"%03d") "_X" num2str(id_X,"%+03d") "_Z" num2str(id_Z,"%+03d")];

endfunction
