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
#include "scafacos_ewald_tcl.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "scafacos.h"

#ifdef SCAFACOS_EWALD

int tclcommand_inter_coulomb_parse_scafacos_ewald(Tcl_Interp * interp, int argc, char ** argv)
{
  double cutoff = -1, alpha = -1, tolerance_field = -1;
  int kmax = -1, maxkmax = -1 ;
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
      
   coulomb.method = COULOMB_SCAFACOS_EWALD;

  if(PERIODIC(0) == 0 || PERIODIC(1) == 0 || PERIODIC(2) == 0) {
    Tcl_AppendResult(interp, "Coulomb scafacos_ewald needs periodic boundaries", (char *) NULL);
    return TCL_ERROR;  
  } 

  while(argc > 0) {
    if(ARG0_IS_S("cutoff")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_ewald.cutoff) && scafacos_ewald.cutoff >= -1)) {
	Tcl_AppendResult(interp, "cutoff expects a positive double", (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("alpha")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_ewald.alpha) && scafacos_ewald.alpha >= -1)) {
	Tcl_AppendResult(interp, "alpha expects a positive double >= -1", (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("kmax")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_ewald.kmax) && scafacos_ewald.kmax > 0)) {
	Tcl_AppendResult(interp, "kmax expects a positive integer", (char *) NULL);
	return TCL_ERROR;
      } 

    } else if(ARG0_IS_S("tolerance_field")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_ewald.tolerance_field) && scafacos_ewald.tolerance_field > 0)) {
	Tcl_AppendResult(interp, "tolerance_field expects a positive double",	 (char *) NULL);
	return TCL_ERROR;
      }

    } else if (ARG0_IS_S("maxkmax")) {
      if (! (argc > 1 && ARG1_IS_I(scafacos_ewald.maxkmax) && scafacos_ewald.maxkmax > 0)) {
	Tcl_AppendResult(interp, "maxkmax expects a positive integer", (char *) NULL);
	return TCL_ERROR;
      }
    }else if(ARG0_IS_S("virial")) {
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

  scafacos.short_range_flag = 1;
    
  Tcl_AppendResult(interp, "Scafacos ewald Solver is set up. \n", (char *) NULL);
  
  return TCL_OK;
}




int tclprint_to_result_scafacos_ewald(Tcl_Interp *interp)
{

  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, scafacos_ewald.cutoff, buffer);
  Tcl_AppendResult(interp, "scafacos_ewald cutoff ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_ewald.alpha, buffer);
  Tcl_AppendResult(interp, "alpha ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_ewald.tolerance_field, buffer);  
  Tcl_AppendResult(interp, "tolerance_field " ,buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d", scafacos_ewald.kmax);
  Tcl_AppendResult(interp,"kmax ", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d", scafacos_ewald.maxkmax); 
  Tcl_AppendResult(interp,"maxkmax ", buffer, " ", (char *) NULL);

  sprintf(buffer,"%d",scafacos.virial);
  Tcl_AppendResult(interp, "virial ", buffer, " ", (char *) NULL);
  
  return TCL_OK;
} 
#endif
