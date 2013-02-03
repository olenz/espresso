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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "integrate.h"
#include "global.h"
#include "grid.h"
#include "domain_decomposition.h"
#include "particle_data.h"
#include "communication.h"
#include "fft.h"
#include "thermostat.h"
#include "scafacos_mmm1d_tcl.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "scafacos.h"

#ifdef SCAFACOS
int tclcommand_inter_coulomb_parse_scafacos_mmm1d(Tcl_Interp * interp, int argc, char ** argv)
{
  double far_switch_radius = 0, maxPWerror = 0;
  int besselcutoff = 0;
  scafacos.virial = 0;
  scafacos.short_range_flag = 1;
  
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

  coulomb.method = COULOMB_SCAFACOS_MMM1D;
    
  
/*periodicity is set hard by scafacos_common_set for scafacos_mmm1d and mmm2d*/
    

  while(argc > 0) {
    if(ARG0_IS_S("far_switch_radius")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_mmm1d.far_switch_radius) && scafacos_mmm1d.far_switch_radius >= 0)) {
	Tcl_AppendResult(interp, "far_switch_radius expects a positive double", (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("besselcutoff")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_mmm1d.besselcutoff) && scafacos_mmm1d.besselcutoff > 0 )) {
	Tcl_AppendResult(interp, "besselcutoff expects positive integer", (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("maxPWerror")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_mmm1d.maxPWerror) && scafacos_mmm1d.maxPWerror > 0)) {
	Tcl_AppendResult(interp, "maxPWerror expects a positive double", (char *) NULL);
	return TCL_ERROR;
      }

    }else if(ARG0_IS_S("virial")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial == 0 || scafacos.virial == 1))) {
	Tcl_AppendResult(interp, "virial expects 0 or 1", (char *) NULL);
	return TCL_ERROR;
      }

    } 
    /* unknown parameter. */
    else break;
    
    argc -= 2;
    argv += 2;
  }
  
  Tcl_AppendResult(interp, "Scafacos mmm1d Solver is set up. \n ", (char *) NULL);
  
  return TCL_OK;
}


int tclprint_to_result_scafacos_mmm1d(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];
  

  Tcl_PrintDouble(interp, scafacos_mmm1d.far_switch_radius, buffer);
  Tcl_AppendResult(interp, "scafacos_mmm1d farswitchrad ", buffer, " ", (char *) NULL);
 
  Tcl_PrintDouble(interp, scafacos_mmm1d.besselcutoff, buffer);  
  Tcl_AppendResult(interp,"besselcutoff ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_mmm1d.maxPWerror, buffer);
  Tcl_AppendResult(interp, "maxPWerror ", buffer, " ", (char *) NULL);
  
  sprintf(buffer, "%d", scafacos.virial);  
  Tcl_AppendResult(interp,"virial " ,buffer, " ", (char *) NULL);

  return TCL_OK;
} 

#endif