/*
  Copyright (C) 2012 The ESPResSo project
  
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

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "scafacos_pepc_tcl.hpp"
#include "utils.hpp"
#include "integrate.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "fft.hpp"
#include "thermostat.hpp"
#include "cells.hpp"
#include "tuning.hpp"
#include "elc.hpp"
#include "global_tcl.hpp"
#include "scafacos.hpp"

#ifdef SCAFACOS_PEPC


int tclcommand_inter_coulomb_parse_scafacos_pepc(Tcl_Interp * interp, int argc, char ** argv)
{
  double epsilon = -1, theta = -1;
  int dipole_correction;
  scafacos.virial = 0;
  scafacos.short_range_flag = 1 ;
  
  switch(coulomb.method){
    case COULOMB_SCAFACOS_DIRECT:
    case COULOMB_SCAFACOS_EWALD:
    case COULOMB_SCAFACOS_FMM:
    case COULOMB_SCAFACOS_MEMD:
    case COULOMB_SCAFACOS_MMM1D:
    case COULOMB_SCAFACOS_MMM2D:
    case COULOMB_SCAFACOS_P2NFFT:
    case COULOMB_SCAFACOS_P3M:
    case COULOMB_SCAFACOS_PEPC:
    case COULOMB_SCAFACOS_PP3MG:
    case COULOMB_SCAFACOS_VMG:
      fcs_destroy(fcs_handle);
      break;
    default:
      break;
  }
  coulomb.method = COULOMB_SCAFACOS_PEPC;

  if(PERIODIC(0) == 0 || PERIODIC(1) == 0 || PERIODIC(2) == 0)
    {
      Tcl_AppendResult(interp, "Coulomb scafacos_pepc needs periodic boundaries", (char *) NULL);
      return TCL_ERROR;  
    }

  while(argc > 0) {
    if(ARG0_IS_S("epsilon")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_pepc.epsilon) && scafacos_pepc.epsilon > -1)) {
	Tcl_AppendResult(interp, "epsilon expects a positive double", (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("theta")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_pepc.theta) && scafacos_pepc.theta  > 0 )) {
	Tcl_AppendResult(interp, "theta a positive double", (char *) NULL);
	return TCL_ERROR;
      }
    }else if(ARG0_IS_S("dipole_correction")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_pepc.dipole_correction) && (scafacos_pepc.dipole_correction == 0 || scafacos_pepc.dipole_correction == 1))) {
	Tcl_AppendResult(interp, "dipole_correction expects 0 or 1", (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("virial")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial == 0 || scafacos.virial == 1))) {
	Tcl_AppendResult(interp, "virial expects 0 or 1", (char *) NULL);
	return TCL_ERROR;
      }
    } 
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }

  Tcl_AppendResult(interp, "scafacos_pepc is set up \n ", (char *) NULL);
  
  return TCL_OK;
}


int tclprint_to_result_scafacos_pepc(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, scafacos_pepc.epsilon, buffer);
  Tcl_AppendResult(interp, "scafacos_pepc epsilon ", buffer, " ", (char *) NULL);

  Tcl_PrintDouble(interp, scafacos_pepc.theta, buffer);
  Tcl_AppendResult(interp, "theta ", buffer, " ", (char *) NULL);
 
  sprintf(buffer, "%d", scafacos_pepc.dipole_correction);
  Tcl_AppendResult(interp, "dipole_correction ", buffer, " ", (char *) NULL);
  
  sprintf(buffer, "%d", scafacos.virial);
  Tcl_AppendResult(interp, "virial ", buffer, " ", (char *) NULL);
  return TCL_OK;
} 
#endif
