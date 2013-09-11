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
#include "scafacos_memd_tcl.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "interaction_data.h"
#include "scafacos.h"

#ifdef SCAFACOS_MEMD
int tclcommand_inter_coulomb_parse_scafacos_memd(Tcl_Interp * interp, int argc, char ** argv){
  double lightspeed = 0,permittivity = 0;
  int mesh_size =0, init_flag =0;
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
  coulomb.method = COULOMB_SCAFACOS_MEMD;

  while(argc > 0) {
    if(ARG0_IS_S("lightspeed")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_memd.lightspeed) && scafacos_memd.lightspeed > 0)) {
	Tcl_AppendResult(interp, "lightspeed expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("mesh_size")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_memd.mesh_size) && scafacos_memd.mesh_size > 0)) {
	Tcl_AppendResult(interp, "mesh_size expects an integer > 0",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("permittivity")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_memd.permittivity) && scafacos_memd.permittivity > 0)) {
	Tcl_AppendResult(interp, "permittivity expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      } 

    } else if(ARG0_IS_S("init_flag")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_memd.init_flag) && (scafacos_memd.init_flag == 0 || scafacos_memd.init_flag == 1))) {
	Tcl_AppendResult(interp, "init_flag expects 0 or 1",
			 (char *) NULL);
	return TCL_ERROR;
      }

    } else if(ARG0_IS_S("virial")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial == 0 || scafacos.virial == 1))) {
	Tcl_AppendResult(interp, "init_flag expects 0 or 1",
			 (char *) NULL);
	return TCL_ERROR;
      }

    } 
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }

  Tcl_AppendResult(interp, "scafacos_memd is set up  ", (char *) NULL);

  return TCL_OK;
}

int tclprint_to_result_scafacos_memd(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];
  /*
  fcs_float lightspeed;
  fcs_memd_get_speed_of_light(fcs_handle, &lightspeed);
  Tcl_PrintDouble(interp, (double)lightspeed, buffer);
  Tcl_AppendResult(interp, "scafacos_p2nfft", buffer, " ", (char *) NULL);
  
  fcs_float permittivity;
  fcs_memd_get_permittivity(fcs_handle, &permittivity);
  Tcl_PrintDouble(interp, permittivity, buffer);
  Tcl_AppendResult(interp,"permittivity ", buffer, " ", (char *) NULL);
  
  fcs_int mesh_size;
  fcs_memd_get_mesh_size_1D(fcs_handle, &mesh_size);
  sprintf(buffer, "%d", (int)mesh_size);  
  Tcl_AppendResult(interp, "mesh_size ",buffer, " ", (char *) NULL);
  
  fcs_int init_flag;
  fcs_memd_get_init_flag(fcs_handle, &init_flag);
  sprintf(buffer,"%d", (int)init_flag);  
  Tcl_AppendResult(interp, "init_flag ", buffer, " ", (char *) NULL);
  
  fcs_int virial;
  fcs_get_virial(fcs_handle, &virial);
  sprintf(buffer, "%d", (int)virial);  
  Tcl_AppendResult(interp, "virial ", buffer, " ", (char *) NULL);
  */
  return TCL_OK;
} 
#endif
