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
#ifndef _SCAFACOS_H
#define _SCAFACOS_H

#include "utils.hpp"

#ifdef SCAFACOS
#include <fcs.h>
#include <mpi.h>

extern FCS fcs_handle;

typedef struct{
  //char  method[];
  int needs_initialize;
  int needs_solver_specific_set;
  int needs_tuning;
  int short_range_flag ;
  int virial;
  double global_r_cut_fcs;
   
  double box_a[3]; 
  double box_c[3]; 
  double box_b[3]; 
  double offset[3];
  int periodicity[3];   

  int n_total_particles;
} scafacos_data_struct;
extern scafacos_data_struct scafacos;

#ifdef SCAFACOS_DIRECT
typedef struct{
  double cutoff;
  int periodic_images[3];
} scafacos_direct_parameter_struct;

extern scafacos_direct_parameter_struct scafacos_direct;
#endif

#ifdef SCAFACOS_EWALD
typedef struct{
  double cutoff;
  int kmax;
  int maxkmax;
  double alpha;
  double tolerance_field;
} scafacos_ewald_parameter_struct;

extern scafacos_ewald_parameter_struct scafacos_ewald;
#endif

#ifdef SCAFACOS_FMM
typedef struct{
  int absrel;
  double tolerance_energy;
  int dipole_correction;
  int internal_tuning;
  double cuspradius;
  int potential; 
} scafacos_fmm_parameter_struct;

extern scafacos_fmm_parameter_struct scafacos_fmm;
#endif

#ifdef SCAFACOS_MEMD
typedef struct{
  double length_x;
  double length_y;
  double length_z;
  double timestep;
  int local_number_of_particles;
  int mesh_size;
  double lightspeed;
  double temperature;
  double permittivity;
  int init_flag;
  
} scafacos_memd_parameter_struct;

extern scafacos_memd_parameter_struct scafacos_memd;
#endif

#ifdef SCAFACOS_MMM1D
typedef struct{
  double far_switch_radius; 
  int besselcutoff;
  double maxPWerror;
} scafacos_mmm1d_parameter_struct;

extern scafacos_mmm1d_parameter_struct scafacos_mmm1d;
#endif

#ifdef SCAFACOS_MMM2D
typedef struct{
  double maxPWerror;
  double far_cutoff;
  double delta_top;
  double delta_bot;
  int layers_per_node;
  double skin;
  int require_total_energy;
} scafacos_mmm2d_parameter_struct;

extern scafacos_mmm2d_parameter_struct scafacos_mmm2d;
#endif

#ifdef SCAFACOS_P2NFFT
typedef struct{
  double cutoff;
  int tolerance_type;
  double tolerance_value;
  double alpha;
  int m;
  int N0;
  int N1;
  int N2;
  int n0;
  int n1;
  int n2;
} scafacos_p2nfft_parameter_struct;

extern scafacos_p2nfft_parameter_struct scafacos_p2nfft;
#endif

#ifdef SCAFACOS_PEPC
typedef struct{
  double epsilon; 
  double theta;  
  int dipole_correction;
} scafacos_pepc_parameter_struct;

extern scafacos_pepc_parameter_struct scafacos_pepc;
#endif

#ifdef SCAFACOS_PP3MG
typedef struct{
  int cells_x; 
  int cells_y; 
  int cells_z; 
  int ghosts; 
  int degree; 
  int max_particles; 
  int max_iterations; 
  double tolerance;
} scafacos_pp3mg_parameter_struct;

extern scafacos_pp3mg_parameter_struct scafacos_pp3mg;
#endif

#ifdef SCAFACOS_VMG
typedef struct{
  int max_level;
  int max_iterations;
  int smoothing_steps; 
  double precision; 
  int near_field_cells;
  int interpolation_order;
  int discretization_order;
  int cycle_type;
} scafacos_vmg_parameter_struct;

extern scafacos_vmg_parameter_struct scafacos_vmg;
#endif

#ifdef SCAFACOS_P3M
typedef struct {
/** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double alpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double r_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  int    grid;
  /** charge assignment order ([0,7]). */
  int    cao;
  /** Accuracy of the actual parameter set. */
  double tolerance_field;
  /** unscaled \ref alpha_L for use with fast inline functions only */
  double alpha;
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  double cutoff;
} scafacos_p3m_parameter_struct;

extern scafacos_p3m_parameter_struct scafacos_p3m;
#endif

int scafacos_run();
void scafacos_tune();
void mpi_bcast_coulomb_method();
void mpi_scafacos_bcast_common_params();
void mpi_scafacos_bcast_solver_specific();
void mpi_scafacos_init();
void mpi_scafacos_set_common();
void mpi_scafacos_solver_specific_set();

void fcs_p3m_count_charged_particles();

#endif

#endif
