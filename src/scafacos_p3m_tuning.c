 

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>



#include "utils.h"
#include "integrate.h"
#include "initialize.h"
#include "global.h"
#include "grid.h"
#include "domain_decomposition.h"
#include "particle_data.h"
#include "communication.h"
#include "fft.h"
#include "thermostat.h"
#include "cells.h"
#include "elc.h"
#include "cells.h"
#include "initialize.h"
#include "interaction_data.h"
#include "errorhandling.h"

#include "tuning.h"
#include "p3m-common.h"
#include "scafacos.h"
#include "scafacos_p3m_tuning.h"
#include "communication.h"


#ifdef SCAFACOS
int scafacos_p3m_tuning();


int scafacos_p3m_tuning(){
  
  int  mesh[3] = {0, 0, 0}, tmp_mesh_points; 
  int tmp_mesh[3];
  double r_cut_iL_min, r_cut_iL_max, r_cut_iL = -1, tmp_r_cut_iL=0.0;
  int cao;
  double                             alpha_L  = -1, tmp_alpha_L=0.0;
  double                             accuracy = -1, tmp_accuracy=0.0;
  double                            time_best=1e20, tmp_time;
  char b[3*ES_INTEGER_SPACE + 3*ES_DOUBLE_SPACE + 128];
  int buffer_cao, buffer_grid;
  double buffer_alpha_L, buffer_r_cut_iL, buffer_tolerance_field;
  
  if (skin == -1) {
//    *log = strcat_alloc(*log, "p3m cannot be tuned, since the skin is not yet set");
    return ES_ERROR;
  }

  
  /* preparation */
  mpi_bcast_event(FCS_P3M_COUNT_CHARGES);
  /* Print Status */
  sprintf(b, "P3M tune parameters: Accuracy goal = %.5e\n", scafacos_p3m.params.tolerance_field);
 // *log = strcat_alloc(*log, b);
  sprintf(b, "System: box_l = %.5e # charged part = %d Sum[q_i^2] = %.5e\n", box_l[0], scafacos_p3m.sum_qpart, scafacos_p3m.sum_q2);
//  *log = strcat_alloc(*log, b);

  if (scafacos_p3m.sum_qpart == 0) {
  //  *log = strcat_alloc(*log, "no charged particles in the system, cannot tune P3M");
    return ES_ERROR;
  }
  double density, r_cut;
  density = n_total_particles /(box_l[0]* box_l[1]* box_l[2]);
  int number = 100;
  r_cut= pow(number/ density/3.14 *3 /4 , 1/3.0);
  if(scafacos_p3m.params.r_cut_iL == 0.0) {
    // WARNING: r_cut_iL_min should not be small !

    r_cut_iL_min = dmin(min_local_box_l, min_box_l/2);
    r_cut_iL_min *= 0.18;
    r_cut_iL_min = r_cut *0.8;
    r_cut_iL_max = dmin(min_local_box_l, min_box_l/2);
    r_cut_iL_max *= 0.28;
    r_cut_iL_max = r_cut *1.5;
    if(r_cut_iL_max >= 4.5)
      r_cut_iL_max = 4.5;
//    r_cut_iL_max = dmin(min_local_box_l, min_box_l/2);// - skin;

//    r_cut_iL_min *= box_l_i[0];
//    r_cut_iL_max *= box_l_i[0];
    fprintf(stderr, "r_cut_iL_max: %f r_cut_iL_min: %f \n", r_cut_iL_max, r_cut_iL_min);
    fprintf(stderr, "box_l is: %f \n", box_l[0]);
  }
  else{
    fprintf(stderr, "Scafacos_p3m tuning varies the cutoff radius; you have to leave cutoff unspecified");
    return ES_ERROR;
  }
  
 // *log = strcat_alloc(*log, "mesh cao r_cut_iL     alpha_L      err          rs_err     ks_err     time [ms]\n");



  for(tmp_r_cut_iL = r_cut_iL_min; tmp_r_cut_iL <= r_cut_iL_max; tmp_r_cut_iL += 0.1){
    scafacos_p3m.params.cutoff = tmp_r_cut_iL;
    
    tmp_time = time_force_calc(1);
//    *log = strcat_alloc(*log, " %f ms ");
    fprintf(stderr, "time: %f \n", tmp_time);
//    *log = strcat_alloc(*log, b);
    
    if ((tmp_time < time_best) && (tmp_time > 0)) {
      time_best = tmp_time;
      fcs_p3m_get_grid(fcs_handle, &buffer_grid);
      fcs_p3m_get_cao(fcs_handle, &buffer_cao);
      fcs_p3m_get_r_cut(fcs_handle, &buffer_r_cut_iL);
      fcs_p3m_get_alpha(fcs_handle, &buffer_alpha_L);
      fcs_p3m_get_tolerance_field(fcs_handle, &buffer_tolerance_field);
      // target accuracy remains the same
    }
    /* no hope of further optimisation */
    //else if (tmp_time > time_best + P3M_TIME_GRAN) {
    //  break;
    //}
  }
  
  
  if(time_best == 1e20) {
    fprintf(stderr, "failed to tune P3M parameters to required accuracy\n");
    return ES_ERROR;
  }

  /* prepare tuned p3m parameters for broadcast*/
  scafacos_p3m.params.cao = buffer_cao;
  scafacos_p3m.params.alpha_L = buffer_alpha_L;
  scafacos_p3m.params.grid = buffer_grid;
  scafacos_p3m.params.cutoff = buffer_r_cut_iL;
  scafacos_p3m.params.alpha_L = buffer_alpha_L;
  /* Tell the user about the outcome */
  fprintf(stderr, "\nresulting parameters:\n%-4d %-3d %.5e %.5e %.5e %-8d\n", buffer_grid, buffer_cao, buffer_r_cut_iL, buffer_alpha_L, buffer_tolerance_field, (int)time_best);
  return ES_OK;
}
#endif /* SCAFACOS */