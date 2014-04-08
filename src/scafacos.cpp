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
#include "scafacos.hpp"

#ifdef SCAFACOS

#include <mpi.h>
#include "grid.hpp"
#include "cells.hpp"
#include "initialize.hpp"
#include "interaction_data.hpp"
#include "tcl/interaction_data_tcl.hpp"
#include "domain_decomposition.hpp"
#include "verlet.hpp"
#include "tcl/global_tcl.hpp"

FCS fcs_handle;
scafacos_data_struct scafacos;

#ifdef SCAFACOS_DIRECT
scafacos_direct_parameter_struct scafacos_direct;
#endif
#ifdef SCAFACOS_EWALD
scafacos_ewald_parameter_struct scafacos_ewald;
#endif
#ifdef SCAFACOS_FMM
scafacos_fmm_parameter_struct scafacos_fmm;
#endif
#ifdef SCAFACOS_MEMD
scafacos_memd_parameter_struct scafacos_memd;
#endif
#ifdef SCAFACOS_MMM1D
scafacos_mmm1d_parameter_struct scafacos_mmm1d;
#endif
#ifdef SCAFACOS_MMM2D
scafacos_mmm2d_parameter_struct scafacos_mmm2d;
#endif
#ifdef SCAFACOS_P2NFFT
scafacos_p2nfft_parameter_struct scafacos_p2nfft;
#endif
#ifdef SCAFACOS_PP3MG
scafacos_pp3mg_parameter_struct scafacos_pp3mg;
#endif
#ifdef SCAFACOS_PEPC
scafacos_pepc_parameter_struct scafacos_pepc;
#endif
#ifdef SCAFACOS_VMG
scafacos_vmg_parameter_struct scafacos_vmg;
#endif
#ifdef SCAFACOS_P3M
scafacos_p3m_parameter_struct scafacos_p3m;
#endif

int scafacos_run (){
  Cell *cell;
  Particle *p;
  int np,c,i,j, dir;
  int count = 0, n_local_particles =0;
  
  // count local_particles
  for (c = 0; c < local_cells.n; c++) {   
    cell = local_cells.cell[c];
    np = cell->n;
    
    n_local_particles += np;
  }
  fcs_float positions[3*n_local_particles];
  fcs_float field[3*n_local_particles];
  fcs_float charges[n_local_particles];
  fcs_float potentials[n_local_particles];
  
  int image_box[3];
  //transform for ScaFaCoS's interface
  //loop over all cells
  for (c = 0; c < local_cells.n; c++) {   
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
       
    //loop over all particles in a cell
    for(i=0;i<np;i++) {  
      //transform into arrays needed by ScaFaCoS
      for (dir=0; dir<3; dir++) {
        //fold each coordinate of each particle back into primary box
        fold_coordinate(p[i].r.p , image_box, dir);
        positions [3*(i+count) + dir] = p[i].r.p[dir];
        //            field [3*(i+count) + j] = p[i].f.f[j];
      }
      charges[i+count] = p[i].p.q;
      //          potentials[i+count] = 0;
    }
    count += np;
  }
  // end transform

 // fprintf(stderr, "scafacos.c: run_scafacos(): fcs_run is being called \n");
  
  if (fcs_run(fcs_handle, (fcs_int)n_local_particles, positions, charges, field, potentials ) != NULL  )
    fprintf(stderr,"Warning: fcs_run exited with an error \n");

  count = 0;
  
   //retransform for Espresso's interface
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
       
        for(i=0;i<np;i++) {   
          for(j=0;j<3;j++){
            p[i].f.f[j] +=  field [3*(i+count) + j] * p[i].p.q * coulomb.bjerrum;
          }
        }
	count += np;
      }
      //printf("node: %d, count: %d \n", this_node, count);

   // end retransform 
  return ES_OK;
}


void scafacos_tune(){
    
  Cell *cell;
  Particle *p;
  int np,c,i,j;
  int count = 0, n_local_particles = 0;
  
  for (c = 0; c < local_cells.n; c++) {   
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    
    n_local_particles += np;
  }
  fcs_float positions[3*n_local_particles];
  fcs_float charges[n_local_particles];
  //printf("node: %d, n_local_particles: %d \n", this_node, n_local_particles);
  

    //loop over all cells
      for (c = 0; c < local_cells.n; c++) {   
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
       
        //loop over all particles in a cell
        for(i=0;i<np;i++) {   
          for (j=0; j<3; j++) {
            positions [3*(i+count) + j] = p[i].r.p[j];
          }
          charges[i+count] = p[i].p.q;
        }
        
        count += np;
      }
     // printf("node: %d, count: %d \n", this_node, count);

  fcs_tune(fcs_handle, (fcs_int)n_local_particles, positions , charges);
  return;
}


void mpi_scafacos_set_common(){
  if (fcs_set_common(fcs_handle, scafacos.short_range_flag, 
                     scafacos.box_a, scafacos.box_b, scafacos.box_c, 
                     scafacos.offset, scafacos.periodicity, 
                     scafacos.n_total_particles) != 0){
    char *errtext = runtime_error(128);
    ERROR_SPRINTF(errtext,"{scafacos.c : scafacos_set_common : fcs_set_common has failed} \n");
  }

  return;
}


void mpi_bcast_coulomb_method(){
  MPI_Bcast(&coulomb.method, 1, MPI_INT, 0, comm_cart);
  return; 
}


void  mpi_scafacos_update_params()
{
    int i;
   // fprintf(stderr, "mpi_scafacos_bcast_common_params is running \n");
    scafacos.box_a[0] = box_l[0];
    scafacos.box_a[1] = 0;
    scafacos.box_a[2] = 0;
    scafacos.box_b[0] = 0;
    scafacos.box_b[1] = box_l[1];
    scafacos.box_b[2] = 0;
    scafacos.box_c[0] = 0;
    scafacos.box_c[1] = 0;
    scafacos.box_c[2] = box_l[2];
    for(i=0;i<3;i++){
      scafacos.offset[i] = 0;
      scafacos.periodicity[i] = 1;
    }

    switch(coulomb.method){
      case COULOMB_SCAFACOS_MMM1D:
        scafacos.periodicity[0] = 0;
        scafacos.periodicity[1] = 0;
        scafacos.periodicity[2] = 1;
        printf( " MMM1D periodicity \n");
        break;
      case COULOMB_SCAFACOS_MMM2D:
        scafacos.periodicity[0] = 1;
        scafacos.periodicity[1] = 1;
        scafacos.periodicity[2] = 0;
        printf( " MMM2D periodicity \n");
        break;
      default:
        break;
    }
    scafacos.n_total_particles = n_part;
}

void mpi_scafacos_bcast_common_params(){

  /*
  MPI_Bcast( &scafacos.short_range_flag, 1, MPI_INT, 0, comm_cart); 
  MPI_Bcast( scafacos.box_a, 3, MPI_DOUBLE, 0, comm_cart); 
  MPI_Bcast( scafacos.box_b, 3, MPI_DOUBLE, 0, comm_cart); 
  MPI_Bcast( scafacos.box_c, 3, MPI_DOUBLE, 0, comm_cart); 
  MPI_Bcast( scafacos.offset, 3, MPI_DOUBLE, 0, comm_cart); 
  MPI_Bcast( scafacos.periodicity, 3, MPI_INT, 0, comm_cart); 
  MPI_Bcast( &scafacos.n_total_particles, 1, MPI_INT, 0, comm_cart); 
  MPI_Bcast(&scafacos.needs_initialize, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&scafacos.needs_solver_specific_set, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&scafacos.needs_tuning, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&scafacos.virial, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&scafacos.global_r_cut_fcs, 1, MPI_DOUBLE, 0 ,comm_cart);

  MPI_Bcast(&coulomb.bjerrum, 1 ,MPI_DOUBLE, 0, comm_cart);*/
  MPI_Bcast(&scafacos, sizeof(scafacos_data_struct), MPI_BYTE, 0, comm_cart);
  return;
}


void mpi_scafacos_bcast_solver_specific(){
 // fprintf(stderr, "mpi_scafacos_bcast_solver_specific is running %d\n", this_node);
  switch(coulomb.method){
#ifdef SCAFACOS_DIRECT
    case COULOMB_SCAFACOS_DIRECT:
      MPI_Bcast(&scafacos_direct, sizeof(scafacos_direct_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_EWALD
    case COULOMB_SCAFACOS_EWALD:
      MPI_Bcast(&scafacos_ewald, sizeof(scafacos_ewald_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_FMM
    case COULOMB_SCAFACOS_FMM:
      MPI_Bcast(&scafacos_fmm, sizeof(scafacos_fmm_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_MEMD
    case COULOMB_SCAFACOS_MEMD:
      MPI_Bcast(&scafacos_memd, sizeof(scafacos_memd_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_MMM1D
    case COULOMB_SCAFACOS_MMM1D:
      MPI_Bcast(&scafacos_mmm1d, sizeof(scafacos_mmm1d_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_MMM2D
    case COULOMB_SCAFACOS_MMM2D:
      MPI_Bcast(&scafacos_mmm2d, sizeof(scafacos_mmm2d_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_P2NFFT
    case COULOMB_SCAFACOS_P2NFFT:
      MPI_Bcast(&scafacos_p2nfft, sizeof(scafacos_p2nfft_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_P3M
    case COULOMB_SCAFACOS_P3M:
      MPI_Bcast(&scafacos_p3m, sizeof(scafacos_p3m_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_PP3MG
    case COULOMB_SCAFACOS_PP3MG:
      MPI_Bcast(&scafacos_pp3mg, sizeof(scafacos_pp3mg_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_VMG
    case COULOMB_SCAFACOS_VMG:
      MPI_Bcast(&scafacos_vmg, sizeof(scafacos_vmg_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
#ifdef SCAFACOS_PEPC
    case COULOMB_SCAFACOS_PEPC:
      MPI_Bcast(&scafacos_pepc, sizeof(scafacos_pepc_parameter_struct), MPI_BYTE, 0, comm_cart);
      break;
#endif
/*    case COULOMB_SCAFACOS_DIRECT:
      MPI_Bcast(&scafacos_direct.cutoff, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(scafacos_direct.periodic_images, 3, MPI_INT, 0, comm_cart);
      break;
      
    case COULOMB_SCAFACOS_EWALD:
      MPI_Bcast(&scafacos_ewald.cutoff, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_ewald.kmax, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_ewald.maxkmax, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_ewald.alpha, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_ewald.tolerance_field, 1, MPI_DOUBLE, 0, comm_cart);
      break;
      
    case COULOMB_SCAFACOS_FMM:
      MPI_Bcast(&scafacos_fmm.absrel, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_fmm.tolerance_energy, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_fmm.dipole_correction, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_fmm.internal_tuning, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_fmm.cuspradius, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_fmm.potential, 1, MPI_INT, 0, comm_cart);
      break;
      
    case COULOMB_SCAFACOS_MEMD:
      MPI_Bcast(&temperature, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&coulomb.bjerrum, 1 , MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_memd.length_x, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_memd.length_y, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_memd.length_z, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_memd.mesh_size, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_memd.lightspeed, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_memd.permittivity, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_memd.init_flag, 1, MPI_INT, 0, comm_cart);
      break;
    
    case COULOMB_SCAFACOS_MMM1D:
      MPI_Bcast(&scafacos_mmm1d.far_switch_radius, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm1d.besselcutoff, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm1d.maxPWerror, 1, MPI_DOUBLE, 0, comm_cart);
      break;
    case COULOMB_SCAFACOS_MMM2D:
      MPI_Bcast(&scafacos_mmm2d.far_cutoff, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm2d.maxPWerror, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm2d.layers_per_node, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm2d.delta_bot, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm2d.delta_top, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm2d.skin, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_mmm2d.require_total_energy, 1, MPI_INT, 0, comm_cart);
      break;
    case COULOMB_SCAFACOS_P2NFFT:
      MPI_Bcast(&scafacos_p2nfft.tolerance_type, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.tolerance_value, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.cutoff, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.alpha, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.m, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.N0, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.N1, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.N2, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.n0, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.n1, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p2nfft.n2, 1, MPI_INT, 0, comm_cart);
      break;
    case COULOMB_SCAFACOS_P3M:
      MPI_Bcast(&scafacos_p3m.params.cutoff, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_p3m.params.alpha, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_p3m.params.grid, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p3m.params.cao, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_p3m.params.tolerance_field, 1, MPI_DOUBLE, 0, comm_cart);
      break;
    case COULOMB_SCAFACOS_PEPC:
      MPI_Bcast(&scafacos_pepc.epsilon, 1 , MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_pepc.theta, 1 , MPI_DOUBLE, 0, comm_cart);
      break;
    case COULOMB_SCAFACOS_PP3MG:
      MPI_Bcast(&scafacos_pp3mg.cells_x, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_pp3mg.cells_y, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_pp3mg.cells_z, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_pp3mg.ghosts, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_pp3mg.max_iterations, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_pp3mg.degree, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_pp3mg.tolerance, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_pp3mg.max_particles, 1, MPI_INT, 0, comm_cart);
      break;
    case COULOMB_SCAFACOS_VMG:
      MPI_Bcast(&scafacos_vmg.max_level, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_vmg.max_iterations, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_vmg.precision, 1, MPI_DOUBLE, 0, comm_cart);
      MPI_Bcast(&scafacos_vmg.near_field_cells, 1, MPI_INT, 0, comm_cart);
      MPI_Bcast(&scafacos_vmg.smoothing_steps, 1, MPI_INT, 0, comm_cart);
      break;*/
    default:
      break;
  }
  return;
}


void mpi_scafacos_solver_specific_set(){

  
 // printf("mpi_scafacos_solver_specific_set is running \n");
  switch(coulomb.method){
#ifdef SCAFACOS_DIRECT
    case COULOMB_SCAFACOS_DIRECT:
      fcs_direct_set_cutoff(fcs_handle, (fcs_float) scafacos_direct.cutoff);
      fcs_direct_set_periodic_images(fcs_handle, scafacos_direct.periodic_images);
      break;
#endif
#ifdef SCAFACOS_EWALD
    case COULOMB_SCAFACOS_EWALD:
      if(scafacos_ewald.cutoff == 0){
        fcs_ewald_set_r_cut_tune(fcs_handle);
      }
      if(scafacos_ewald.cutoff != 0){
        fcs_ewald_set_r_cut(fcs_handle, (fcs_float)scafacos_ewald.cutoff);
      }
      if(scafacos_ewald.alpha == 0){
        fcs_ewald_set_alpha_tune(fcs_handle);
      }
      if(scafacos_ewald.alpha != 0){
        fcs_ewald_set_alpha(fcs_handle, (fcs_float) scafacos_ewald.alpha);
      }
      if(scafacos_ewald.tolerance_field == 0){
        fcs_ewald_set_tolerance_field_tune(fcs_handle);
      }
      if(scafacos_ewald.tolerance_field != 0){
        fcs_ewald_set_tolerance_field(fcs_handle, (fcs_float)scafacos_ewald.tolerance_field);
      }
      if(scafacos_ewald.kmax == 0){
        fcs_ewald_set_kmax_tune(fcs_handle);
      }
      if(scafacos_ewald.kmax != 0){
        fcs_ewald_set_kmax(fcs_handle, (fcs_int) scafacos_ewald.kmax);
      }
      if(scafacos_ewald.maxkmax == 0){
        fcs_ewald_set_maxkmax_tune(fcs_handle);
      }
      if(scafacos_ewald.maxkmax != 0){
        fcs_ewald_set_maxkmax(fcs_handle, (fcs_int) scafacos_ewald.maxkmax);
      }
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1); 
      break;
#endif
#ifdef SCAFACOS_FMM
    case COULOMB_SCAFACOS_FMM:
      fcs_fmm_set_absrel(fcs_handle, (fcs_int)scafacos_fmm.absrel);
      fcs_fmm_set_tolerance_energy(fcs_handle, (fcs_float)scafacos_fmm.tolerance_energy);
      fcs_fmm_set_dipole_correction(fcs_handle, (fcs_int)scafacos_fmm.dipole_correction);
      if(scafacos_fmm.internal_tuning == 0)
        fcs_fmm_set_internal_tuning(fcs_handle, 0LL);
      if(scafacos_fmm.internal_tuning == 1)
        fcs_fmm_set_internal_tuning(fcs_handle, 1LL);
      if(scafacos_fmm.cuspradius != 0)
        fcs_fmm_set_cusp_radius(fcs_handle, (fcs_float)scafacos_fmm.cuspradius);
      if(scafacos_fmm.potential != 0)
        fcs_fmm_set_potential(fcs_handle, (fcs_int)scafacos_fmm.potential);
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1);
      break;
#endif
#ifdef SCAFACOS_MEMD
    case COULOMB_SCAFACOS_MEMD:
/*
      fcs_memd_set_box_size(fcs_handle, box_l[0],  box_l[1], box_l[2]);
      fcs_memd_set_total_number_of_particles(fcs_handle, (fcs_int) n_total_particles);
      fcs_memd_set_bjerrum_length(fcs_handle, coulomb.bjerrum);
      fcs_memd_set_time_step(fcs_handle, (fcs_float) time_step);
      fcs_memd_set_temperature(fcs_handle, (fcs_float) temperature);
      Cell *cell;
      int np,c;
      int count = 0;
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        np = cell->n;
  
	count += np;
      }
      fcs_memd_set_local_number_of_particles(fcs_handle, (fcs_int) count);
      
      if(scafacos_memd.init_flag != 0)
        fcs_memd_set_init_flag(fcs_handle, (fcs_int) scafacos_memd.init_flag);
      if(scafacos_memd.mesh_size != 0)
        fcs_memd_set_mesh_size_1D(fcs_handle, (fcs_int) scafacos_memd.mesh_size);
      if(scafacos_memd.lightspeed !=0)
        fcs_memd_set_speed_of_light(fcs_handle, (fcs_float) scafacos_memd.lightspeed);
      if (scafacos_memd.permittivity != 0)
        fcs_memd_set_permittivity(fcs_handle, (fcs_float) scafacos_memd.permittivity);
        
*/
      break; 
#endif
#ifdef SCAFACOS_MMM1D
    case COULOMB_SCAFACOS_MMM1D:

      if(scafacos_mmm1d.far_switch_radius != 0){
        fcs_mmm1d_set_far_switch_radius(fcs_handle, (fcs_float)scafacos_mmm1d.far_switch_radius);
      }
      if(scafacos_mmm1d.besselcutoff != 0){
        fcs_mmm1d_set_bessel_cutoff(fcs_handle, (fcs_int)scafacos_mmm1d.besselcutoff);
      }
      if(scafacos_mmm1d.maxPWerror != 0){
        fcs_mmm1d_set_maxPWerror(fcs_handle, (fcs_float)scafacos_mmm1d.maxPWerror);
      }
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1);
      break;
#endif
#ifdef SCAFACOS_MMM2D
    case COULOMB_SCAFACOS_MMM2D:

      if(scafacos_mmm2d.far_cutoff != 0){
        fcs_mmm2d_set_far_cutoff(fcs_handle, (fcs_float)scafacos_mmm2d.far_cutoff);
        }
      if(scafacos_mmm2d.maxPWerror != 0){
       fcs_mmm2d_set_maxPWerror(fcs_handle, (fcs_float)scafacos_mmm2d.maxPWerror);
      }
      if(scafacos_mmm2d.layers_per_node != 0){
        fcs_mmm2d_set_layers_per_node(fcs_handle, (fcs_int)scafacos_mmm2d.layers_per_node);
      }
      if(scafacos_mmm2d.delta_top != 0 && scafacos_mmm2d.delta_bot != 0){
        fcs_mmm2d_set_dielectric_contrasts(fcs_handle, (fcs_float)scafacos_mmm2d.delta_top, (fcs_float)scafacos_mmm2d.delta_bot);
      }
      if(scafacos_mmm2d.skin != 0){
        fcs_mmm2d_set_skin(fcs_handle, (fcs_float)scafacos_mmm2d.skin);
      }
      if(scafacos_mmm2d.require_total_energy != 0){
        fcs_mmm2d_require_total_energy(fcs_handle, (fcs_int) scafacos_mmm2d.require_total_energy);
      }
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1);

      break;
#endif
#ifdef SCAFACOS_P2NFFT
    case COULOMB_SCAFACOS_P2NFFT:

      if(scafacos_p2nfft.tolerance_type != 0){
        fcs_p2nfft_set_tolerance(fcs_handle, (fcs_int)scafacos_p2nfft.tolerance_type, (fcs_float)scafacos_p2nfft.tolerance_value);
	int buffer1;
	double buffer2;
	fcs_p2nfft_get_tolerance(fcs_handle, &buffer1, &buffer2);
	printf("tolerance_type is: %d, tolerance_value is: %f \n", buffer1, buffer2);
      }
      if(scafacos_p2nfft.tolerance_type == 0){
        fcs_p2nfft_set_ignore_tolerance(fcs_handle, 1);
		printf("tolerance_type tune \n");
      }
      if(scafacos_p2nfft.cutoff != 0){
        fcs_p2nfft_set_r_cut(fcs_handle, (fcs_float)scafacos_p2nfft.cutoff);
      }
      if(scafacos_p2nfft.cutoff == 0){
        fcs_p2nfft_set_r_cut_tune(fcs_handle);
      }
      if(scafacos_p2nfft.alpha != 0){
        fcs_p2nfft_set_alpha(fcs_handle, (fcs_float)scafacos_p2nfft.alpha);
      }
      if(scafacos_p2nfft.alpha == 0){
        fcs_p2nfft_set_alpha_tune(fcs_handle);
      }
      if(scafacos_p2nfft.m != 0){
        fcs_p2nfft_set_m(fcs_handle, (fcs_int)scafacos_p2nfft.m);
      }
      if(scafacos_p2nfft.m == 0){
        fcs_p2nfft_set_m_tune(fcs_handle);
      }
      if(scafacos_p2nfft.N0 != 0 && scafacos_p2nfft.N1 != 0 && scafacos_p2nfft.N2 != 0){
        fcs_p2nfft_set_gridsize(fcs_handle, (fcs_int)scafacos_p2nfft.N0, (fcs_int)scafacos_p2nfft.N1, (fcs_int)scafacos_p2nfft.N2);
      }
      if (scafacos_p2nfft.N0 == 0 || scafacos_p2nfft.N1 == 0 || scafacos_p2nfft.N2 == 0){
        fcs_p2nfft_set_gridsize_tune(fcs_handle);
      }
      if(scafacos_p2nfft.n0 != 0 && scafacos_p2nfft.n1 != 0 && scafacos_p2nfft.n2 != 0){
        fcs_p2nfft_set_oversampled_gridsize(fcs_handle, (fcs_int) scafacos_p2nfft.n0, (fcs_int) scafacos_p2nfft.n1, (fcs_int) scafacos_p2nfft.n2);
      }
      if(scafacos_p2nfft.n0 == 0 && scafacos_p2nfft.n1 == 0 && scafacos_p2nfft.n2 == 0){
        fcs_p2nfft_set_oversampled_gridsize_tune(fcs_handle);
      }
      break;
#endif
#ifdef SCAFACOS_P3M
    case COULOMB_SCAFACOS_P3M:

      if(scafacos_p3m.cutoff == 0){
        fcs_p3m_set_r_cut_tune(fcs_handle);
      }
      if(scafacos_p3m.cutoff!= 0){
       fcs_p3m_set_r_cut(fcs_handle, (fcs_float)scafacos_p3m.cutoff);
      }
      if(scafacos_p3m.alpha == 0){
        fcs_p3m_set_alpha_tune(fcs_handle);
      }
      if(scafacos_p3m.alpha != 0){
        fcs_p3m_set_alpha(fcs_handle, (fcs_float)scafacos_p3m.alpha);
      }
      if(scafacos_p3m.grid == 0){
        fcs_p3m_set_grid_tune(fcs_handle);
      }
      if(scafacos_p3m.grid != 0){
        fcs_p3m_set_grid(fcs_handle, (fcs_int)scafacos_p3m.grid);
      }
      if(scafacos_p3m.cao == 0){
        fcs_p3m_set_cao_tune(fcs_handle);
      }
      if(scafacos_p3m.cao != 0){
        fcs_p3m_set_cao(fcs_handle, (fcs_int)scafacos_p3m.cao);
      }
      if(scafacos_p3m.tolerance_field == 0){
        fcs_p3m_set_tolerance_field_tune(fcs_handle);
      }
      if(scafacos_p3m.tolerance_field != 0){
        fcs_p3m_set_tolerance_field(fcs_handle, (fcs_float)scafacos_p3m.tolerance_field);
      }
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1);
      break;
#endif
#ifdef SCAFACOS_PEPC
    case COULOMB_SCAFACOS_PEPC:
      if(scafacos_pepc.epsilon != 0){
       fcs_pepc_set_epsilon(fcs_handle, (fcs_float)scafacos_pepc.epsilon);
      }
      if(scafacos_pepc.theta != 0){
        fcs_pepc_set_theta(fcs_handle, (fcs_float)scafacos_pepc.theta);
      }
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1);

      break;
#endif
#ifdef SCAFACOS_PP3MG
    case COULOMB_SCAFACOS_PP3MG:
      if(scafacos_pp3mg.cells_x != 0){
        fcs_pp3mg_set_cells_x(fcs_handle, (fcs_int) scafacos_pp3mg.cells_x);
      }
      if(scafacos_pp3mg.cells_y != 0){
        fcs_pp3mg_set_cells_y(fcs_handle, (fcs_int) scafacos_pp3mg.cells_y);
      }
      if(scafacos_pp3mg.cells_z != 0){
        fcs_pp3mg_set_cells_z(fcs_handle, (fcs_int) scafacos_pp3mg.cells_z);
      }  
      if(scafacos_pp3mg.ghosts != 0){
        fcs_pp3mg_set_ghosts(fcs_handle, (fcs_int) scafacos_pp3mg.ghosts);
      }
      if(scafacos_pp3mg.max_iterations != 0){
        fcs_pp3mg_set_max_iterations(fcs_handle, (fcs_int) scafacos_pp3mg.max_iterations);
      }
      if(scafacos_pp3mg.degree != 0){
        fcs_pp3mg_set_degree(fcs_handle, (fcs_int) scafacos_pp3mg.degree);
      }
      if(scafacos_pp3mg.tolerance != 0){
        fcs_pp3mg_set_tol(fcs_handle, (fcs_float) scafacos_pp3mg.tolerance);
      }
      if(scafacos_pp3mg.max_particles != 0){
        fcs_pp3mg_set_max_particles(fcs_handle, (fcs_int) scafacos_pp3mg.max_particles);
      }
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1);
      break;
#endif
#ifdef SCAFACOS_VMG
    case COULOMB_SCAFACOS_VMG:
      if(scafacos_vmg.max_level != 0) {
        fcs_vmg_set_max_level( fcs_handle, (fcs_int)scafacos_vmg.max_level );
      }
      if(scafacos_vmg.max_iterations != 0) {
        fcs_vmg_set_max_iterations( fcs_handle, (fcs_int)scafacos_vmg.max_iterations );
      }
      if(scafacos_vmg.cycle_type != 0) {
        fcs_vmg_set_cycle_type( fcs_handle, (fcs_int)scafacos_vmg.cycle_type );
      }
      printf("vmg.precision: %f \n", scafacos_vmg.precision);
      if(scafacos_vmg.precision != 0) {
        fcs_vmg_set_precision( fcs_handle, (fcs_float)scafacos_vmg.precision );
      }
      if(scafacos_vmg.near_field_cells != 0) {
        fcs_vmg_set_near_field_cells( fcs_handle, (fcs_int)scafacos_vmg.near_field_cells );
      }
      if(scafacos_vmg.smoothing_steps != 0) {
        fcs_vmg_set_smoothing_steps( fcs_handle, (fcs_int)scafacos_vmg.smoothing_steps );
      }
      if(scafacos.virial != 0)
        fcs_require_virial(fcs_handle, 1);
      break;
#endif
    default:
      break;
  }
  return;
}

void mpi_scafacos_init(){
  MPI_Bcast(&coulomb.method, 1, MPI_INT, 0, comm_cart);
  switch(coulomb.method){
    case COULOMB_SCAFACOS_DIRECT:
      //scafacos.method = "direct";
      //MPI_Bcast(&scafacos.method, 6, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "direct", comm_cart);
      break;
    case COULOMB_SCAFACOS_EWALD:
      //scafacos.method = "ewald";
      //MPI_Bcast(&scafacos.method, 5, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "ewald", comm_cart);
      break;
    case COULOMB_SCAFACOS_FMM:
      //scafacos.method = "fmm";
      //MPI_Bcast(&scafacos.method, 3, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "fmm", comm_cart);
      break;
    case COULOMB_SCAFACOS_MEMD:
      //scafacos.method = "memd";
      //MPI_Bcast(&scafacos.method, 4, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "memd", comm_cart);
      break;
    case COULOMB_SCAFACOS_MMM1D:
      //scafacos.method = "mmm1d";
      //MPI_Bcast(&scafacos.method, 5, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "mmm1d", comm_cart);
      break;
    case COULOMB_SCAFACOS_MMM2D:
      //scafacos.method = "mmm2d";
      //MPI_Bcast(&scafacos.method, 5, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "mmm2d", comm_cart);
      break;
    case COULOMB_SCAFACOS_P2NFFT:
      //scafacos.method = "p2nfft";
      //MPI_Bcast(&scafacos.method, 6, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "p2nfft", comm_cart);
      break;
    case COULOMB_SCAFACOS_P3M:
      //scafacos.method = "p3m";
      //MPI_Bcast(&scafacos.method, 3, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "p3m", comm_cart);
      break;
    case COULOMB_SCAFACOS_PEPC:
      //scafacos.method = "pepc";
      //MPI_Bcast(&scafacos.method, 4, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "pepc", comm_cart);
      break;
    case COULOMB_SCAFACOS_PP3MG:
      //scafacos.method = "pp3mg";
      //MPI_Bcast(&scafacos.method, 5, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "pp3mg", comm_cart);
      break;
    case COULOMB_SCAFACOS_VMG:
      //scafacos.method = "vmg";
      //MPI_Bcast(&scafacos.method, 3, MPI_CHAR, 0, comm_cart);
      fcs_init(&fcs_handle, "vmg", comm_cart);
      break;
    default:
      break;
  }
  return;
}
#endif /*ifdef SCAFACOS */
