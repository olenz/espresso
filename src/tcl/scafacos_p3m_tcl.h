/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

#include "config.h"
#include <tcl.h>

#include "interaction_data_tcl.h"
#include "scafacos.h"
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
#include "cells.h"
#include "tuning.h"
#include "elc.h"
#include "global_tcl.h"
#include "interaction_data.h"
#include "initialize.h"

#include "scafacos.h"
#include "scafacos_p3m_tuning.h"

#ifdef SCAFACOS

/// parse the basic scafacos P3M parameters
int tclcommand_inter_coulomb_parse_scafacos_p3m(Tcl_Interp * interp, int argc, char ** argv);
int tclcommand_inter_coulomb_parse_scafacos_p3m_conventional_tune(Tcl_Interp * interp, int argc, char ** argv);
int tclcommand_inter_coulomb_parse_scafacos_p3m_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive);
int tclprint_to_result_scafacos_p3m(Tcl_Interp *interp);

#endif