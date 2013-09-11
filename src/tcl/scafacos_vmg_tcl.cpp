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
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "utils.hpp"
#include "integrate.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "fft.hpp"
#include "thermostat.hpp"
#include "scafacos_vmg_tcl.hpp"
#include "cells.hpp"
#include "tuning.hpp"
#include "elc.hpp"
#include "global_tcl.hpp"
#include "scafacos.hpp"
#include "interaction_data.hpp"

#ifdef SCAFACOS_VMG


int tclcommand_inter_coulomb_parse_scafacos_vmg(Tcl_Interp * interp, int argc, char ** argv){
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
  coulomb.method = COULOMB_SCAFACOS_VMG;

  if(PERIODIC(0) == 0 || PERIODIC(1) == 0 || PERIODIC(2) == 0) {
      Tcl_AppendResult(interp, "Coulomb scafacos_vmg needs periodic boundaries", (char *) NULL);
      return TCL_ERROR;  
    }

  while(argc > 0) {
    if(ARG0_IS_S("cycle_type")) {
      if (! (argc > 1 && ARG1_IS_I(scafacos_vmg.cycle_type) && (scafacos_vmg.cycle_type == 1 || scafacos_vmg.cycle_type == 2))) {
	Tcl_AppendResult(interp, "cycle_type expects 1 or 2", (char *) NULL);
	return TCL_ERROR;
      }

    } else if(ARG0_IS_S("precision")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_vmg.precision) && scafacos_vmg.precision > -1 )) {
	Tcl_AppendResult(interp, "precision expects nonnegative double", (char *) NULL);
	return TCL_ERROR;
      }

    } else if(ARG0_IS_S("max_level")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_vmg.max_level) && scafacos_vmg.max_level > 0)) {
	Tcl_AppendResult(interp, "max_level expects a positive integer", (char *) NULL);
	return TCL_ERROR;
      }

    } else if(ARG0_IS_S("smoothing_steps")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_vmg.smoothing_steps) && scafacos_vmg.smoothing_steps >0 )) {
	Tcl_AppendResult(interp, "smoothing_steps expects positive integer", (char *) NULL);
	return TCL_ERROR;
      }

    } else if(ARG0_IS_S("max_iterations")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_vmg.max_iterations) && scafacos_vmg.max_iterations > 0)) {
	Tcl_AppendResult(interp, "max_iterations expects a positive integer", (char *) NULL);
	return TCL_ERROR;
      }

    } else if(ARG0_IS_S("near_field_cells")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_vmg.near_field_cells) && scafacos_vmg.near_field_cells > 0)) {
	Tcl_AppendResult(interp, "near_field_cells expects a positive integer", (char *) NULL);
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

  Tcl_AppendResult(interp, "scafacos_vmg is set up \n ", (char *) NULL);

  return TCL_OK;
}



int tclprint_to_result_scafacos_vmg(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];
  
  Tcl_PrintDouble(interp, scafacos_vmg.precision, buffer);
  Tcl_AppendResult(interp,"scafacos_vmg precision ", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_vmg.max_level);
  Tcl_AppendResult(interp, "max_level ", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_vmg.max_iterations);
  Tcl_AppendResult(interp, "max_iterations", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_vmg.smoothing_steps);
  Tcl_AppendResult(interp, "smoothing_steps ",buffer, " ", (char *) NULL);

  sprintf(buffer,"%d",scafacos_vmg.cycle_type);
  Tcl_AppendResult(interp,"cycle_type ", buffer, " ", (char *) NULL);

  sprintf(buffer,"%d",scafacos_vmg.near_field_cells);
  Tcl_AppendResult(interp, "near_field_cells ",buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos.virial);
  Tcl_AppendResult(interp, "virial ",buffer, " ", (char *) NULL);
  return TCL_OK;
} 
#endif
