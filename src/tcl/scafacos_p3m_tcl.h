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

#ifndef _TCL_SCAFACOS_P3M_TCL_H
#define _TCL_SCAFACOS_P3M_TCL_H

#include "config.h"
#include <tcl.h>

#ifdef SCAFACOS_P3M

/// parse the basic scafacos P3M parameters
int tclcommand_inter_coulomb_parse_scafacos_p3m(Tcl_Interp * interp, int argc, char ** argv);
int tclcommand_inter_coulomb_parse_scafacos_p3m_conventional_tune(Tcl_Interp * interp, int argc, char ** argv);
int tclcommand_inter_coulomb_parse_scafacos_p3m_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive);
int tclprint_to_result_scafacos_p3m(Tcl_Interp *interp);

#endif
#endif
