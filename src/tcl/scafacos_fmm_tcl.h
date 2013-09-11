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


#include "config.h"
#include <tcl.h>

#include "interaction_data_tcl.h"

#ifdef SCAFACOS_FMM

/// parse the basic fmm parameters
int tclcommand_inter_coulomb_parse_scafacos_fmm(Tcl_Interp * interp, int argc, char ** argv);

int tclprint_to_result_scafacos_fmm(Tcl_Interp *interp);
#endif /* of ifdef SCAFACOS_FMM */

