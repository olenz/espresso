/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef LJGEN_H
#define LJGEN_H

/** \file ljgen.h
 *  Routines to calculate the generalized lennard jones energy and/or  force 
 *  for a particle pair. "Generalized" here means that the LJ energy is of the
 *  form
 *
 *  eps * [ b1 * (sigma/(r-r_offset))^a1 - b2 * (sigma/(r-r_offset))^a2 + shift]
 *
 *  \ref forces.c
*/

/* These headers are needed to define types used in this header, hence they
 * are included here.  */
#include "particle_data.h"
#include "interaction_data.h"

int tclprint_to_result_ljgenIA(Tcl_Interp *interp, int i, int j);

int tclcommand_inter_parse_ljgen(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv);


/** Calculate lennard Jones force between particle p1 and p2 */
void add_ljgen_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3]);

/** calculate Lennard jones energy between particle p1 and p2. */
double ljgen_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist);

/** calculate lj_capradius from lj_force_cap */
void calc_ljgen_cap_radii(double force_cap);

/* LJGEN_H */
#endif 
