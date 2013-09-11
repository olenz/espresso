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
#include "scafacos_p2nfft_tcl.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "scafacos.h"
#include "interaction_data.h"

#ifdef SCAFACOS_P2NFFT

int tclcommand_inter_coulomb_parse_scafacos_p2nfft(Tcl_Interp * interp, int argc, char ** argv){ 
  double tolerance_value = -1, alpha = -1, r_cut = -1;;
  int  m = -1, n0 = -1, n1 = -1, n2 = -1, N0 = -1, N1 = -1, N2 = -1, tolerance_type = -1;
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

  coulomb.method = COULOMB_SCAFACOS_P2NFFT;


#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 || PERIODIC(1) == 0 || PERIODIC(2) == 0) {
      Tcl_AppendResult(interp, "Coulomb scafacos_p2nfft needs periodic boundaries", (char *) NULL);
      return TCL_ERROR;  
  }
#endif
    while(argc > 0) {
    if(ARG0_IS_S("r_cut")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_p2nfft.cutoff) && scafacos_p2nfft.cutoff >= -1)) {
	Tcl_AppendResult(interp, "r_cut expects a positive double ", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } 
    else if(ARG0_IS_S("tolerance_type")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_p2nfft.tolerance_type) && scafacos_p2nfft.tolerance_type > -1)) {
	Tcl_AppendResult(interp, "type expects an nonnegative integer ", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } 
    else if(ARG0_IS_S("tolerance_value")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_p2nfft.tolerance_value) && scafacos_p2nfft.tolerance_value > -1 )) {
	Tcl_AppendResult(interp, "tolerance_value expects apositive double", (char *) NULL);
	return TCL_ERROR;
      } 
      argc -= 2;
      argv += 2;
    } 
    else if(ARG0_IS_S("alpha")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_p2nfft.alpha) && scafacos_p2nfft.alpha > 0)) {
	Tcl_AppendResult(interp, "alpha expects a positive double",	 (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } 
    else if(ARG0_IS_S("m")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_p2nfft.m) && scafacos_p2nfft.m > 0)) {
	Tcl_AppendResult(interp, "m expects a positive integer",	 (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } 
    else if (ARG0_IS_S("grid")) {
      if (! (argc > 1 && ARG1_IS_I(scafacos_p2nfft.N0) && ARG_IS_I(2,scafacos_p2nfft.N1) && ARG_IS_I(3,scafacos_p2nfft.N2) && 
	scafacos_p2nfft.N0 > 0 && scafacos_p2nfft.N1 > 0 && scafacos_p2nfft.N2 > 0)) {
	Tcl_AppendResult(interp, "grid expects positive integers", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 4;
      argv += 4;
    }
    else if (ARG0_IS_S("oversampled_grid")) {
      if (! (argc > 1 && ARG1_IS_I(scafacos_p2nfft.n0) && ARG_IS_I(2,scafacos_p2nfft.n1) && ARG_IS_I(3,scafacos_p2nfft.n2) && 
	scafacos_p2nfft.n0 > 0 && scafacos_p2nfft.n1 > 0 && scafacos_p2nfft.n2 > 0)) {
	Tcl_AppendResult(interp, "oversampled_grid expects positive integers", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 4;
      argv += 4;
    }
    else if(ARG0_IS_S("srf")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.short_range_flag) && (scafacos.short_range_flag ==1 || scafacos.short_range_flag == 0))) {
	Tcl_AppendResult(interp, "srf expects 0 or 1",	 (char *) NULL);
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
 
  Tcl_AppendResult(interp ,"scafacos p2nfft is set up \n  ",(char*) NULL);
  return TCL_OK;
}


int tclprint_to_result_scafacos_p2nfft(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];

  sprintf(buffer, "%d", scafacos_p2nfft.tolerance_type);
  Tcl_AppendResult(interp, "scafacos_p2nfft tolerance_type ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, scafacos_p2nfft.tolerance_value, buffer);
  Tcl_AppendResult(interp, "tolerance_value ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_p2nfft.cutoff, buffer);
  Tcl_AppendResult(interp, "r_cut " ,buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_p2nfft.N0);
  Tcl_AppendResult(interp, "grid ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",scafacos_p2nfft.N1);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",scafacos_p2nfft.N2);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_p2nfft.alpha, buffer);
  Tcl_AppendResult(interp, "alpha ", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_p2nfft.m);
  Tcl_AppendResult(interp,"m ", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos.short_range_flag);
  Tcl_AppendResult(interp, "srf ", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos.virial);
  Tcl_AppendResult(interp, "virial ", buffer, " ", (char *) NULL);
  
  return TCL_OK;
} 
#endif
