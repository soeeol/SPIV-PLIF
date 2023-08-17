##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## return path for record identificator
##
## recordid: <date>-<project>-<run>
## sample recordid: 20220211-M1-1
##
## Author: Sören J. Gerke
##

function idpath = id_path (recordid)
	idx = strchr (recordid, "-");
	switch (length (idx))
		case 2
		  idpath = [recordid(1:idx(1)-1) "/" recordid(idx(1)+1:idx(2)-1) "/" ...
	    						recordid(idx(2)+1:end) "/"];
		otherwise
		  error ("id2path: wrong recordid structure?");
	endswitch
endfunction
