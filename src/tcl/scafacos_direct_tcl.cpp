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

#include "scafacos_direct_tcl.hpp"

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
#include "interaction_data.hpp"
#include "scafacos.hpp"


#ifdef SCAFACOS_DIRECT

int tclcommand_inter_coulomb_parse_scafacos_direct_images(Tcl_Interp * interp, int argc, char ** argv);

int tclcommand_inter_coulomb_parse_scafacos_direct(Tcl_Interp * interp, int argc, char ** argv){
  
  double cutoff = 0;
  int images[3]={0,0,0};
  int i;
  scafacos.virial =0;
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
      
  coulomb.method = COULOMB_SCAFACOS_DIRECT;

 
  while(argc > 0) {
    if(ARG0_IS_S("images")) {
      if (! (argc > 1 && ARG1_IS_I(scafacos_direct.periodic_images[0]) && ARG_IS_I(2,scafacos_direct.periodic_images[1])&& ARG_IS_I(3,scafacos_direct.periodic_images[2]) && 
	scafacos_direct.periodic_images[0] >= 0 && scafacos_direct.periodic_images[1] >= 0 && scafacos_direct.periodic_images[2] >= 0)){
	Tcl_AppendResult(interp, "images expects a nonnegative integer", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 4;
      argv += 4;
    }
    else if (ARG0_IS_S("cutoff")){
      if(! (argc > 1 && ARG1_IS_D(scafacos_direct.cutoff) && scafacos_direct.cutoff >= 0)){
	Tcl_AppendResult(interp, "cutoff must be double", (char *) NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    }
    else if(ARG0_IS_S("virial")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial == 0 || scafacos.virial == 1))) {
	 Tcl_AppendResult(interp, "virial expects 0 or 1", (char *) NULL);
	 return TCL_ERROR;
       }

    } 

    
    else break;
  }
      
  Tcl_AppendResult(interp, "Scafacos Direct Solver is set up. \n", (char *) NULL);
  
  return TCL_OK;
}



int tclprint_to_result_scafacos_direct(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  

  Tcl_PrintDouble(interp, scafacos_direct.cutoff, buffer);
  Tcl_AppendResult(interp, "scafacos_direct cutoff ", buffer, " ", (char *) NULL);

  
  Tcl_AppendResult(interp, "images ", (char *) NULL);
  sprintf(buffer,"%d",scafacos_direct.periodic_images[0]);
  Tcl_AppendResult(interp,"", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",scafacos_direct.periodic_images[1]);
  Tcl_AppendResult(interp,"", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",scafacos_direct.periodic_images[2]);
  Tcl_AppendResult(interp,"", buffer, " ", (char *) NULL);
  return TCL_OK;
}

#endif
