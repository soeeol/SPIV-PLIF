##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## wall indicator map for u  field
##
## Author: Sören J. Gerke
##

function map_out = ind_wall_u (map_in, cut_off)
	map_out = abs (map_in);
	maxxy = max (max (map_out));
	for i = 1:numel(map_out(1,:))
		maxx = max (map_out(:,i));
		map_out(map_out(:,i)>cut_off*maxx,i) = cut_off * maxxy;
	endfor
endfunction
