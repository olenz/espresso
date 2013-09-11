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
#include "scafacos_fmm_tcl.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "interaction_data.h"
#include "scafacos.h"

#ifdef SCAFACOS_FMM
int tclcommand_inter_coulomb_parse_scafacos_fmm(Tcl_Interp * interp, int argc, char ** argv){
  double tolerance_energy = 0.001, cuspradius = 0;
  int absrel = 0, dipole_correction = 0, potential = 0, internal_tuning = 0;
  scafacos.short_range_flag = 1;
  scafacos.virial = 0;
  
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
  coulomb.method = COULOMB_SCAFACOS_FMM; 

#ifdef PARTIAL_PERIODIC
    if(PERIODIC(0) == 0 || PERIODIC(1) == 0 || PERIODIC(2) == 0)    {
      Tcl_AppendResult(interp, "Coulomb scafacos_fmm needs periodic boundaries!\n", (char *) NULL);
      return TCL_ERROR;  
    } 
#endif 
  
  while(argc > 0) {
    if(ARG0_IS_S("tolerance_energy")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_fmm.tolerance_energy) && scafacos_fmm.tolerance_energy >= 0)) {
	Tcl_AppendResult(interp, "tolerance_energy expects a nonegative double", (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("absrel")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_fmm.absrel) && absrel >= 0 && scafacos_fmm.absrel <= 2)) {
	Tcl_AppendResult(interp, "absrel expects an integer 0,1,2", (char *) NULL);
	return TCL_ERROR;
      }     
    } else if(ARG0_IS_S("dipole_correction")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_fmm.dipole_correction) && scafacos_fmm.dipole_correction >= -1 && scafacos_fmm.dipole_correction <= 1)) {
	Tcl_AppendResult(interp, "dipole_correction expects -1, 0, 1", (char *) NULL);
	return TCL_ERROR;
      } 
    }else if(ARG0_IS_S("virial")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial == 0 || scafacos.virial == 1))) {
	Tcl_AppendResult(interp, "virial expects 0 or 1", (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("internal_tuning")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_fmm.internal_tuning) && (scafacos_fmm.internal_tuning == 0 || scafacos_fmm.internal_tuning == 1))) {
	Tcl_AppendResult(interp, "dipole_correction expects 0 (homogenous) or 1 (inhomogenous)", (char *) NULL);
	return TCL_ERROR;
      } 
    }else if(ARG0_IS_S("potential")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_fmm.potential) && (scafacos_fmm.potential == 64 || scafacos_fmm.potential == 65))) {
	Tcl_AppendResult(interp, "potential expects 64 (Coulomb) or 65 (Cusp)", (char *) NULL);
	return TCL_ERROR;
      } 
    }else if(ARG0_IS_S("cuspradius")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_fmm.cuspradius) && scafacos_fmm.cuspradius >= 0)) {
	Tcl_AppendResult(interp, "cuspradius expects a nonnegative double", (char *) NULL);
	return TCL_ERROR;
      } 
    }
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }
    

/*
#define FCS_FMM_COULOMB 64
#define FCS_FMM_CUSP 65

#define FCS_FMM_NO_DIPOLE_CORRECTION -1
#define FCS_FMM_STANDARD_DIPOLE_CORRECTION 0
#define FCS_FMM_ACTIVE_DIPOLE_CORRECTION 1

#define FCS_FMM_STANDARD_ERROR 0
#define FCS_FMM_CUSTOM_ABSOLUTE 1
#define FCS_FMM_CUSTOM_RELATIVE 2

#define FCS_FMM_INHOMOGENOUS_SYSTEM 1LL
#define FCS_FMM_HOMOGENOUS_SYSTEM 0LL
*/

  Tcl_AppendResult(interp, "Scafacos fmm Solver is set up. \n", (char *) NULL);
  
  return TCL_OK;
}




int tclprint_to_result_scafacos_fmm(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];
  
  sprintf(buffer,"%d",scafacos_fmm.absrel);
  Tcl_AppendResult(interp, "scafacos_fmm absrel ", buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_fmm.tolerance_energy, buffer);
  Tcl_AppendResult(interp, "tolerance_energy ", buffer, " ", (char *) NULL);  
  
  sprintf(buffer,"%d",scafacos_fmm.dipole_correction);
  Tcl_AppendResult(interp, "dipole_correction " ,buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d",scafacos_fmm.potential);
  Tcl_AppendResult(interp, "potential " ,buffer," ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_fmm.cuspradius, buffer);
  Tcl_AppendResult(interp, "cuspradius", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d", scafacos_fmm.internal_tuning);
  Tcl_AppendResult(interp, "internal_tuning " ,buffer, " ",(char *) NULL);
  
  sprintf(buffer,"%d",(int)scafacos.virial);
  Tcl_AppendResult(interp,"virial ", buffer , " ", (char *) NULL);

  return TCL_OK;
} 
#endif
