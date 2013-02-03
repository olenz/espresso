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
#include "scafacos_mmm2d_tcl.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "scafacos.h"

#ifdef SCAFACOS
int tclcommand_inter_coulomb_parse_scafacos_mmm2d(Tcl_Interp * interp, int argc, char ** argv)
{
  double far_cutoff = -1, maxPWerror = -1, delta_bot = -1, delta_top = -1, skin = -1;
  int layers_per_node = -1, require_total_energy;
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

 coulomb.method = COULOMB_SCAFACOS_MMM2D;

 
 /* periodicity is set hard by scafacos_common_set for scafacos_mmm1d and mmm2d */

  while(argc > 0) {
    if(ARG0_IS_S("far_cutoff")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_mmm2d.far_cutoff) && scafacos_mmm2d.far_cutoff >= 0)) {
	Tcl_AppendResult(interp, "far_cutoff expects a positive double ", (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("layers_per_node")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_mmm2d.layers_per_node) && scafacos_mmm2d.layers_per_node > -1)) {
	Tcl_AppendResult(interp, "layers_per_node expects an nonnegative integer ", (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("maxPWerror")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_mmm2d.maxPWerror) && scafacos_mmm2d.maxPWerror > 0 )) {
	Tcl_AppendResult(interp, "maxPWerror expects apositive double", (char *) NULL);
	return TCL_ERROR;
      } 
    } else if(ARG0_IS_S("delta_bot")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_mmm2d.delta_bot) && scafacos_mmm2d.delta_bot > 0)) {
	Tcl_AppendResult(interp, "delta_bot expects a positive double",	 (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("delta_top")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_mmm2d.delta_top) && scafacos_mmm2d.delta_top > 0)) {
	Tcl_AppendResult(interp, "delta_top expects a positive double",	 (char *) NULL);
	return TCL_ERROR;
      }
    } else if (ARG0_IS_S("skin")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_mmm2d.skin) && scafacos_mmm2d.skin >= 0)) {
	Tcl_AppendResult(interp, "skin expects an nonnegative double", (char *) NULL);
	return TCL_ERROR;
      }
    }else if(ARG0_IS_S("virial")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial == 0 || scafacos.virial == 1))) {
	Tcl_AppendResult(interp, "virial expects 0 or 1", (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("total_energy")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_mmm2d.require_total_energy) && (scafacos_mmm2d.require_total_energy == 0 || scafacos_mmm2d.require_total_energy == 1))) {
	Tcl_AppendResult(interp, "total_energy expects 0 or 1", (char *) NULL);
	return TCL_ERROR;
      }
    } 
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }

  Tcl_AppendResult(interp, "Scafacos mmm2d Solver is set up. \n ", (char *) NULL);
  
  return TCL_OK;
}




int tclprint_to_result_scafacos_mmm2d(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];

  sprintf(buffer,"%d",scafacos_mmm2d.layers_per_node);  
  Tcl_AppendResult(interp, "scafacos_mmm2d layers_per_node ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_mmm2d.maxPWerror, buffer);
  Tcl_AppendResult(interp, "maxPWerror ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_mmm2d.far_cutoff, buffer);
  Tcl_AppendResult(interp, "far_cutoff ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_mmm2d.delta_top, buffer);
  Tcl_AppendResult(interp, "delta_top ",buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, scafacos_mmm2d.delta_bot, buffer);
  Tcl_AppendResult(interp, "delta_bot ",buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_mmm2d.skin, buffer);  
  Tcl_AppendResult(interp,"skin ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_mmm2d.total_energy, buffer);   
  Tcl_AppendResult(interp, "total_energy ", buffer, " ", (char *) NULL);

  sprintf(buffer,"%d",virial);
  Tcl_AppendResult(interp, "virial ", buffer, " ", (char *) NULL);
 
  return TCL_OK;
} 

#endif