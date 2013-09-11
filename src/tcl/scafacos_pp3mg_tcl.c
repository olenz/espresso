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
#include "scafacos_pp3mg_tcl.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "scafacos.h"

#ifdef SCAFACOS_PP3MG

int tclcommand_inter_coulomb_parse_scafacos_pp3mg(Tcl_Interp * interp, int argc, char ** argv){
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
  coulomb.method = COULOMB_SCAFACOS_PP3MG;

  if(PERIODIC(0) == 0 || PERIODIC(1) == 0 || PERIODIC(2) == 0){
    Tcl_AppendResult(interp, "Coulomb scafacos_pp3mg needs periodic boundaries", (char *) NULL);
    return TCL_ERROR;  
  }

  while(argc > 0) {
    if(ARG0_IS_S("ghosts")) {
      if (! (argc > 1 && ARG1_IS_I(scafacos_pp3mg.ghosts) && scafacos_pp3mg.ghosts > -1)) {
	Tcl_AppendResult(interp, "ghosts expects a positive double", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } else if(ARG0_IS_S("tolerance")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_pp3mg.tolerance) && scafacos_pp3mg.tolerance > -1 )) {
	Tcl_AppendResult(interp, "tolerance expects nonnegative double", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } else if(ARG0_IS_S("degree")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_pp3mg.degree) && scafacos_pp3mg.degree > 0)) {
	Tcl_AppendResult(interp, "degree expects a positive integer", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } else if(ARG0_IS_S("cells")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_pp3mg.cells_x) && ARG_IS_I(2,scafacos_pp3mg.cells_y) && ARG_IS_I(3, scafacos_pp3mg.cells_z) && 
	scafacos_pp3mg.cells_x >0 && scafacos_pp3mg.cells_y >0 && scafacos_pp3mg.cells_z >0)) {
	Tcl_AppendResult(interp, "cells expects 3 positive integers", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 4;
      argv += 4;
    } else if(ARG0_IS_S("max_iterations")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_pp3mg.max_iterations) && scafacos_pp3mg.max_iterations > 0)) {
	Tcl_AppendResult(interp, "max_iterations expects a positive integer", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } else if(ARG0_IS_S("max_particles")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_pp3mg.max_particles) && scafacos_pp3mg.max_particles > 0)) {
	Tcl_AppendResult(interp, "max_particles expects a positive integer", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } else if(ARG0_IS_S("virial")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial == 0 || scafacos.virial == 1))) {
	Tcl_AppendResult(interp, "virial expects 0 or 1", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } 
    /* unknown parameter. Probably one of the optionals */
    else break;
       
  }

  Tcl_AppendResult(interp, "scafacos_pp3mg is set up \n ", (char *) NULL);
  
  return TCL_OK;
}



int tclprint_to_result_scafacos_pp3mg(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  
  sprintf(buffer,"%d",scafacos_pp3mg.cells_x);
  Tcl_AppendResult(interp, "scafacos_pp3mg cells ", buffer, " ", (char *) NULL);
 
  sprintf(buffer,"%d",scafacos_pp3mg.cells_y);
  Tcl_AppendResult(interp,  buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_pp3mg.cells_z);
  Tcl_AppendResult(interp,  buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_pp3mg.tolerance, buffer);
  Tcl_AppendResult(interp, "tolerance ",buffer, " ", (char *) NULL);
 
  sprintf(buffer,"%d",scafacos_pp3mg.ghosts);
  Tcl_AppendResult(interp, "ghosts ",buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_pp3mg.degree);
  Tcl_AppendResult(interp, "degree ",buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_pp3mg.max_iterations);
  Tcl_AppendResult(interp, "max_iterations ",buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_pp3mg.max_particles);
  Tcl_AppendResult(interp, "max_particles ",buffer, " ", (char *) NULL);
 
  return TCL_OK;
} 
#endif
