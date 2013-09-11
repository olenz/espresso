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

#include "scafacos_p3m_tcl.hpp"
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
#include "initialize.hpp"

#include "scafacos.hpp"


/************************************************************/
#ifdef SCAFACOS
//int tclcommand_inter_coulomb_parse_scafacos_p3m_tune(Tcl_Interp * interp, int argc, char **argv);
int tclcommand_inter_coulomb_parse_scafacos_p3m(Tcl_Interp * interp, int argc, char ** argv);
int tclprint_to_result_scafacos_p3m(Tcl_Interp *interp);

// int tclcommand_inter_coulomb_parse_scafacos_p3m_tune(Tcl_Interp * interp, int argc, char **argv){
  
//   double tolerance_field = 0;
//   scafacos.short_range_flag = 1;
 
//   while(argc > 0) {
//   if(ARG0_IS_S("tolerance_field")) {
//       if(! (argc > 1 && ARG1_IS_D(scafacos_p3m.params.tolerance_field) && scafacos_p3m.params.tolerance_field > 0)) {
// 	Tcl_AppendResult(interp, "tolerance_field expects a positive double",
// 			 (char *) NULL);
// 	return TCL_ERROR;
//       }
//   } else if(ARG0_IS_S("srf")) {
//     if(! (argc > 1 && ARG1_IS_I(scafacos.short_range_flag) && (scafacos.short_range_flag ==1 || scafacos.short_range_flag == 0))) {
//       Tcl_AppendResult(interp, "srf expects 0 or 1",	 (char *) NULL);
//       return TCL_ERROR;
//     }
//   }  else if(ARG0_IS_S("virial")) {
//       if(! (argc > 1 && ARG1_IS_I(scafacos.virial) && (scafacos.virial ==1 || scafacos.virial == 0))) {
//         Tcl_AppendResult(interp, "virial expects 0 or 1",	 (char *) NULL);
//         return TCL_ERROR;
//       }
//     } 
//   /* unknown parameter. Probably one of the optionals */
//     else break;
    
//     argc -= 2;
//     argv += 2;
//   }
  
//   /* check for optional parameters */
//   if (argc > 0) {
//       return TCL_ERROR;
//   }
  
//     /* do the tuning */  


//   if (scafacos_p3m_tuning() == ES_ERROR) {  
//     Tcl_AppendResult(interp, "\nfailed to tune scafacos p3m parameters to required accuracy \n", (char *) NULL);
//     return TCL_ERROR;
//   }
//   Tcl_AppendResult(interp,"tuning finished \n", (char *) NULL);
//   return TCL_OK;
  
// }


// int tclcommand_inter_coulomb_parse_scafacos_p3m_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive)
// { /*
//   int grid = 0, cao = 0, n_interpol = -1,  virial = 0;
//   double r_cut = 0, alpha = -1, tolerance = 0;
//   coulomb.method = COULOMB_SCAFACOS_P3M;
//   short_range_flag = 1;
//   tune_scafacos= 1;
  
//   if(fcs_handle){
//     fcs_destroy(fcs_handle);
//   }

//   while(argc > 0) {
//     if(ARG0_IS_S("r_cut")) {
//       if (! (argc > 1 && ARG1_IS_D(r_cut) && r_cut >= -1)) {
// 	Tcl_AppendResult(interp, "r_cut expects a positive double",
// 			 (char *) NULL);
// 	return TCL_ERROR;
//       }
      
//     } else if(ARG0_IS_S("grid")) {
//       if(! (argc > 1 && ARG1_IS_I(grid) && grid >= -1)) {
// 	Tcl_AppendResult(interp, "grid expects an integer >= -1",
// 			 (char *) NULL);
// 	return TCL_ERROR;
//       }
      
//     } else if(ARG0_IS_S("cao")) {
//       if(! (argc > 1 && ARG1_IS_I(cao) && cao >= -1 && cao <= 7)) {
// 	Tcl_AppendResult(interp, "cao expects an integer between -1 and 7",
// 			 (char *) NULL);
// 	return TCL_ERROR;
//       } 

//     } else if(ARG0_IS_S("tolerance")) {
//       if(! (argc > 1 && ARG1_IS_D(tolerance) && tolerance > 0)) {
// 	Tcl_AppendResult(interp, "tolerance expects a positive double",
// 			 (char *) NULL);
// 	return TCL_ERROR;
//       }

//     } else if (ARG0_IS_S("alpha")) {
//       if (! (argc > 1 && ARG1_IS_D(alpha) && alpha >= 0)) {
// 	Tcl_AppendResult(interp, "alpha expects a nonnegative double", (char *) NULL);
// 	return TCL_ERROR;
//       }
//     }else if(ARG0_IS_S("virial")) {
//       if(! (argc > 1 && ARG1_IS_I(virial) && (virial == 0 || virial == 1))) {
// 	Tcl_AppendResult(interp, "virial expects 0 or 1", (char *) NULL);
// 	return TCL_ERROR;
//       }

//     }
//     else if(ARG0_IS_S("srf")) {
//       if(! (argc > 1 && ARG1_IS_I(short_range_flag) && (short_range_flag ==1 || short_range_flag == 0))) {
// 	Tcl_AppendResult(interp, "srf expects 0 or 1",	 (char *) NULL);
// 	return TCL_ERROR;
//       }
//     } 
//    //  unknown parameter. Probably one of the optionals 
//     else break;
    
//     argc -= 2;
//     argv += 2;
//   }

//   scafacos_p3m_pre_init();
//   scafacos_p3m_set_tune_params( r_cut, grid, cao, -1.0, tolerance, n_interpol);


//  //  check for optional parameters 
//   if (argc > 0) {
//       return TCL_ERROR;
//   }

//   // do the tuning 
//   char *log = NULL;
//   if (scafacos_p3m_adaptive_tune(&log) == ES_ERROR) {  
//     Tcl_AppendResult(interp, log, "\nfailed to tune SCAFACOS P3M parameters to required accuracy", (char *) NULL);
//     if (log)
//       free(log);
//     return TCL_ERROR;
//   }
//   // Tell the user about the tuning outcome 
//   Tcl_AppendResult(interp, log, (char *) NULL);
//   if (log)
//     free(log);

//   Tcl_AppendResult(interp, "\nscafacos_p3m is set up \n", (char *) NULL);
  
//   return TCL_OK;
//   */
// }


int tclcommand_inter_coulomb_parse_scafacos_p3m(Tcl_Interp * interp, int argc, char ** argv)
{

  int grid = -1, cao = -1;
  double tolerance_field = -1, alpha = -1, cutoff = -1;
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
  coulomb.method = COULOMB_SCAFACOS_P3M;
    
  
  
#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 || PERIODIC(1) == 0 || PERIODIC(2) == 0) {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with Coulomb scafacos_P3M", (char *) NULL);
      return TCL_ERROR;  
   }
#endif

  while(argc > 0) {
    // if (ARG0_IS_S("tune"))
    //   return tclcommand_inter_coulomb_parse_scafacos_p3m_tune(interp, argc-1, argv+1);

    if(ARG0_IS_S("cutoff")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_p3m.params.cutoff) && scafacos_p3m.params.cutoff > 0)) {
	Tcl_AppendResult(interp, "cutoff expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("grid")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_p3m.params.grid) && scafacos_p3m.params.grid >= -1)) {
	Tcl_AppendResult(interp, "grid expects an integer >= -1",
			 (char *) NULL);
	return TCL_ERROR;
      }
    } else if(ARG0_IS_S("cao")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos_p3m.params.cao) && scafacos_p3m.params.cao >= -1 && scafacos_p3m.params.cao <= 7)) {
	Tcl_AppendResult(interp, "cao expects an integer between -1 and 7",
			 (char *) NULL);
	return TCL_ERROR;
      } 
    } else if(ARG0_IS_S("tolerance_field")) {
      if(! (argc > 1 && ARG1_IS_D(scafacos_p3m.params.tolerance_field) && scafacos_p3m.params.tolerance_field > 0)) {
	Tcl_AppendResult(interp, "tolerance_field expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }
    } else if (ARG0_IS_S("alpha")) {
      if (! (argc > 1 && ARG1_IS_D(scafacos_p3m.params.alpha) && scafacos_p3m.params.alpha_L >= 0 && scafacos_p3m.params.alpha_L <= 1)) {
	Tcl_AppendResult(interp, "alpha expects a double between 0 and 1", (char *) NULL);
	return TCL_ERROR;
      }
    }
    else if(ARG0_IS_S("srf")) {
      if(! (argc > 1 && ARG1_IS_I(scafacos.short_range_flag) && (scafacos.short_range_flag ==1 || scafacos.short_range_flag == 0))) {
	Tcl_AppendResult(interp, "srf expects 0 or 1",	 (char *) NULL);
	return TCL_ERROR;
      }
    } 
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }

  Tcl_AppendResult(interp, "scafacos_p3m is set up \n", (char *) NULL);

  return TCL_OK;
}


/*********************** miscelanea of functions *************************************/

 int tclprint_to_result_scafacos_p3m(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, scafacos_p3m.params.cutoff, buffer);
  Tcl_AppendResult(interp, "scafacos_p3m cutoff ", buffer, "  ", (char *) NULL);
  

  sprintf(buffer,"%d",scafacos_p3m.params.grid);
  Tcl_AppendResult(interp, "grid ", buffer, " ", (char *) NULL);
  Tcl_AppendResult(interp,  buffer, " ", (char *) NULL);
  Tcl_AppendResult(interp,  buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d", scafacos_p3m.params.cao);
  Tcl_AppendResult(interp, "cao ",buffer, " ", (char *) NULL);
  
  Tcl_PrintDouble(interp, scafacos_p3m.params.alpha, buffer);
  Tcl_AppendResult(interp, "alpha ",buffer, " ", (char *) NULL);

  Tcl_PrintDouble(interp, scafacos_p3m.params.tolerance_field, buffer);
  Tcl_AppendResult(interp, "tolerance_field ", buffer, " ", (char *) NULL);


  sprintf(buffer,"%d",scafacos.short_range_flag);
  Tcl_AppendResult(interp, "srf ", buffer, " ", (char *) NULL);
  
  sprintf(buffer,"%d", scafacos.virial);
  Tcl_AppendResult(interp, "virial ", buffer, " ", (char *) NULL);


  return TCL_OK;
} 

#endif /* SCAFACOS */
