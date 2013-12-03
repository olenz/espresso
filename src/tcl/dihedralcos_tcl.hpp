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
#ifndef _DIHEDRALCOS_TCL_HPP
#define _DIHEDRALCOS_TCL_HPP
/** \file dihedralcos_tcl.hpp
 * Tcl interface for \ref dihedralcos.hpp
 */

#include "parser.hpp"
#include "dihedralcos.hpp"
#include "interaction_data.hpp"

///
int tclprint_to_result_dihedralcosIA(Tcl_Interp *interp,
                                     Bonded_ia_parameters *params);

/// parse parameters for the dihedral potential
inline int tclcommand_inter_parse_dihedralcos(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, bend1, bend2;

  if (argc < 4 ) {
    Tcl_AppendResult(interp, "dihedralcos needs 3 parameters: "
		     "<bend> <bend1> <bend2>", (char *) NULL);
    return (ES_ERROR);
  }
  if ( !ARG_IS_D(1, bend) || !ARG_IS_D(2, bend1) || !ARG_IS_D(3, bend2) ) {
    Tcl_AppendResult(interp, "dihedralcos needs 3 parameters of types DOUBLE DOUBLE DOUBLE: "
		     "<bend> <bend1> <bend2> ", (char *) NULL);
    return ES_ERROR;
  }
  
  CHECK_VALUE(dihedralcos_set_params(bond_type, bend, bend1, bend2), "bond type must be nonnegative");
}

#endif

