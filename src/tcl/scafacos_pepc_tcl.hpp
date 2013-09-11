/*
  Copyright (C) 2012,2013 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef _TCL_SCAFACOS_PEPC_TCL_H
#define _TCL_SCAFACOS_PEPC_TCL_H

#include "config.hpp"
#include <tcl.h>

#ifdef SCAFACOS_PEPC
int tclprint_to_result_scafacos_pepc(Tcl_Interp *interp);

/// parse the basic direct parameters
int tclcommand_inter_coulomb_parse_scafacos_pepc(Tcl_Interp * interp, int argc, char ** argv);

#endif /* of ifdef SCAFACOS_PEPC */
#endif
