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
#ifndef CONSTRAINT_H
#define CONSTRAINT_H
#include "statistics.h"
#include "energy.h"
#include "forces.h"
#include "grid.h"
#include "errorhandling.h"
#include "tunable_slip.h"

/** \file constraint.h
 *  Routines for handling of constraints.
 *  Only active if the feature CONSTRAINTS is activated.
 *  see also \ref interaction_data.h
 */

#ifdef CONSTRAINTS

// for the charged rod "constraint"
#define C_GAMMA   0.57721566490153286060651209008

extern int reflection_happened ;

MDINLINE double sign(double x) {
  if (x > 0)
    return 1.;
  else
    return -1;
}
MDINLINE void calculate_wall_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_wall *c, double *dist, double *vec)
{
  int i;

  *dist = -c->d;
  for(i=0;i<3;i++) *dist += ppos[i]*c->n[i];
  
  for(i=0;i<3;i++) vec[i] = c->n[i] * *dist;
  
}


MDINLINE void calculate_sphere_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_sphere *c, double *dist, double *vec)
{
  int i;
  double fac,  c_dist;

  c_dist=0.0;
  for(i=0;i<3;i++) {
    vec[i] = c->pos[i] - ppos[i];
    c_dist += SQR(vec[i]);
  }
  c_dist = sqrt(c_dist);
  
  if ( c->direction == -1 ) {
  /* apply force towards inside the sphere */
    *dist = c->rad - c_dist;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= fac;
  } else {
   /* apply force towards outside the sphere */
    *dist = c_dist - c->rad;
    fac = *dist / c_dist;
    for(i=0;i<3;i++) vec[i] *= -fac;
  }
}


MDINLINE void calculate_maze_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_maze *c, double *dist, double *vec)
{
  int i,min_axis,cursph[3],dim;
  double diasph,fac,c_dist,sph_dist,cyl_dist,temp_dis;
  double sph_vec[3],cyl_vec[3];

  dim=(int) c->dim;
  diasph = box_l[0]/c->nsphere;

  /* First determine the distance to the sphere */
  c_dist=0.0;
  for(i=0;i<3;i++) {
    cursph[i] = (int) (ppos[i]/diasph);
    sph_vec[i] = (cursph[i]+0.5) * diasph  - (ppos[i]);
    c_dist += SQR(sph_vec[i]);
  }
  c_dist = sqrt(c_dist);
  sph_dist = c->sphrad - c_dist;
  fac = sph_dist / c_dist;
  for(i=0;i<3;i++) cyl_vec[i] = sph_vec[i];
  for(i=0;i<3;i++) sph_vec[i] *= fac;
  
  /* Now calculate the cylinder stuff */
  /* have to check all 3 cylinders */
  min_axis=2;
  cyl_dist=sqrt(cyl_vec[0]*cyl_vec[0]+cyl_vec[1]*cyl_vec[1]);
  
  if(dim > 0 ){
    temp_dis=sqrt(cyl_vec[0]*cyl_vec[0]+cyl_vec[2]*cyl_vec[2]);
    if ( temp_dis < cyl_dist) {
        min_axis=1;
        cyl_dist=temp_dis;
    }

    if(dim > 1 ){
        temp_dis=sqrt(cyl_vec[1]*cyl_vec[1]+cyl_vec[2]*cyl_vec[2]);
        if ( temp_dis < cyl_dist) {
            min_axis=0;
            cyl_dist=temp_dis;
        }
    }
  }
  cyl_vec[min_axis]=0.;
  
  c_dist=cyl_dist;
  cyl_dist = c->cylrad - c_dist;
  fac = cyl_dist / c_dist;
  for(i=0;i<3;i++) cyl_vec[i] *= fac;
  
  /* Now decide between cylinder and sphere */
  if ( sph_dist > 0 ) {
    if ( sph_dist>cyl_dist ) {
        *dist = sph_dist;
        for(i=0;i<3;i++) vec[i] = sph_vec[i];
    } else {
        *dist = cyl_dist;
        for(i=0;i<3;i++) vec[i] = cyl_vec[i];  
    }
  } else {
    *dist = cyl_dist;
    for(i=0;i<3;i++) vec[i] = cyl_vec[i];  
  }
}

MDINLINE void calculate_cylinder_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_cylinder *c, double *dist, double *vec)
{
  int i;
  double d_per,d_par,d_real,d_per_vec[3],d_par_vec[3],d_real_vec[3];

  d_real = 0.0;
  for(i=0;i<3;i++) {
    d_real_vec[i] = ppos[i] - c->pos[i];
    d_real += SQR(d_real_vec[i]);
  }
  d_real = sqrt(d_real);
    
  d_par=0.;
  for(i=0;i<3;i++) {
    d_par += (d_real_vec[i] * c->axis[i]);
  }
    
  for(i=0;i<3;i++) {
    d_par_vec[i] = d_par * c->axis[i] ;
    d_per_vec[i] = ppos[i] - (c->pos[i] + d_par_vec[i]) ;
  }
		
  d_per=sqrt(SQR(d_real)-SQR(d_par));
  d_par = fabs(d_par) ;

  if ( c->direction == -1 ) {
    /*apply force towards inside cylinder */
    d_per = c->rad - d_per ;
    d_par = c->length - d_par;
    if (d_per < d_par )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= -d_per_vec[i] * d_per /  (c->rad - d_per) ;
      }
    } else {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= -d_par_vec[i] * d_par /  (c->length - d_par) ;
      }
    }
  } else {
    /*apply force towards outside cylinder */
    d_per = d_per - c->rad ;
    d_par = d_par - c->length ;
    if (d_par < 0 )  {
      *dist = d_per ;   
      for (i=0; i<3;i++) {
	vec[i]= d_per_vec[i] * d_per /  (d_per + c->rad) ;
      }
    } else if ( d_per < 0) {
      *dist = d_par ;
      for (i=0; i<3;i++) {
	vec[i]= d_par_vec[i] * d_par /  (d_par + c->length) ;
      }
    } else {
      *dist = sqrt( SQR(d_par) + SQR(d_per)) ;
      for (i=0; i<3;i++) {
	vec[i]=
	  d_per_vec[i] * d_per /  (d_per + c->rad) +
	  d_par_vec[i] * d_par /  (d_par + c->length) ;
      }	
    }
  }
}

MDINLINE void calculate_pore_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_pore *c, double *dist, double *vec)
{
  int i;
  double c_dist[3];           /* cartesian distance from pore center */
  double z , r;             /* cylindrical coordinates, coordinate system parallel to
                               pore, origin at pore centera */
  double z_vec[3], r_vec[3]; /* cartesian vectors that correspond to these coordinates */
  double e_z[3], e_r[3];    /* unit vectors in the cylindrical coordinate system */
//  double d_per, d_per_sh, d_par, d_par_sh, d_real, d_real_2,
//    d_per_vec[3], d_par_vec[3], d_real_vec[3];
  /* helper variables, for performance reasons should the be move the the
   * constraint struct*/
     double slope, average_radius, z_left, z_right;
  /* and another helper that is hopefully optmized out */
     double norm, help_dist, z_0;
     double c1_r, c1_z, c2_r, c2_z;

     
     c->smoothing_radius = 1.;

      
     slope = (c->rad_right - c->rad_left)/2./c->length;
     average_radius = 0.5*(c->rad_left + c->rad_right);

  /* compute the position relative to the center of the pore */
  for(i=0;i<3;i++) {
    c_dist[i] = ppos[i] - c->pos[i]; 
  } 
  
  /* compute the component parallel to the pore axis */
  z =0.; 
  for(i=0;i<3;i++) {
    z += (c_dist[i] * c->axis[i]);
  }
  
  /* decompose the position into parallel and perpendicular to the axis */
  r = 0.;
  for(i=0;i<3;i++) {
    z_vec[i] = z * c->axis[i]; 
    r_vec[i] = c_dist[i] - z_vec[i];
    r += r_vec[i]*r_vec[i];
  }
  r = sqrt(r);


  /* calculate norm and unit vectors for both */
  norm = 0;
  for(i=0;i<3;i++) 
    norm += z_vec[i]*z_vec[i]; 
  norm = sqrt(norm);
  for(i=0;i<3;i++) 
    e_z[i] = z_vec[i]/norm;
  norm = 0;
  for(i=0;i<3;i++) 
    norm += r_vec[i]*r_vec[i]; 
  norm = sqrt(norm);
  for(i=0;i<3;i++) 
    e_r[i] = r_vec[i]/norm;
  
  /* c?_r/z and are the centers of the circles that are used to smooth 
   * the entrance of the pore in cylindrical coordinates*/
  c1_z = - (c->length - c->smoothing_radius);
  c2_z = + (c->length - c->smoothing_radius);
  z_left = c1_z - sign(slope) * sqrt(slope*slope/(1+slope*slope))*c->smoothing_radius;
  z_right = c2_z + sign(slope) * sqrt(slope*slope/(1+slope*slope))*c->smoothing_radius;

  c1_r = c->rad_left + slope * ( z_left + c->length ) +
      sqrt( c->smoothing_radius * c->smoothing_radius  - SQR( z_left - c1_z ) );
  c2_r = c->rad_left + slope * ( z_right + c->length ) +
      sqrt( c->smoothing_radius * c->smoothing_radius  - SQR( z_right - c2_z ) );
//  printf("r %f z %f c1_r %f c1_z %f c2_r %f c2_z%f\n", r, z, c1_r, c1_z, c2_r, c2_z);
 
  /* Check if we are in the region of the left wall */
  if ( r > c1_r && z < c1_z ) {
    *dist = -z - c->length;
    for (i=0; i<3;i++) 
      vec[i] = -e_z[i]*(*dist) ;
    return;
  }
  /* Check if we are in the region of the right wall */
  if ( r > c2_r && z > c2_z ) {
    *dist = +z - c->length;
    for (i=0; i<3;i++) 
      vec[i] = +e_z[i]*(*dist) ;
    return;
  }

  /* check if the particle sould feel the smoothed ends or the middle of the pore */
  /* calculate aufpunkt in z direction first.   */

  /* the distance of the particle from the pore cylinder/cone calculated by projection on the
   * cone normal. Should be > 0 if particle is inside the pore */
  help_dist = (- ( r - average_radius ) + (z * slope)); ///sqrt(1+slope*slope);
  z_0 = z + help_dist * slope / sqrt(1+slope*slope);

  /* Check if we are in the range of the cone in the center */
  if (z_0 > z_left && z_0 < z_right) { 
    *dist = help_dist;
    for (i=0; i<3;i++) 
      vec[i]= - e_r[i] * help_dist; 
    return;
  }
  /* Check if we are in the range of the left smoothing circle */
  if ( z_0 < z_left ) {
    /* distance from the smoothing center */
    norm = sqrt( (z - c1_z)*(z - c1_z) + (r - c1_r)*(r - c1_r) );
    *dist = norm - c->smoothing_radius;
    for (i=0; i<3;i++) 
      vec[i] = (c->smoothing_radius/norm - 1)*(z - c1_z)*e_z[i] - (c->smoothing_radius/norm -1)*(r - c1_r)*e_r[i];
    return;
  }
  /* Check if we are in the range of the right smoothing circle */
  if ( z_0 > z_right ) {
    norm = sqrt( (z - c2_z)*(z - c2_z) + (r - c2_r)*(r - c2_r) );
    *dist = norm - c->smoothing_radius;
    for (i=0; i<3;i++) 
      vec[i] = -(c->smoothing_radius/norm - 1)*(z - c2_z)*e_z[i] - (c->smoothing_radius/norm -1)*(r - c2_r)*e_r[i];
//    if (*dist < -0.3)
//      printf("warning - reflection w length %f in right smoo area\n", *dist);
    return;
  }
  exit(printf("should never be reached, z %f, r%f\n",z, r));

}

MDINLINE void calculate_plane_dist(Particle *p1, double ppos[3], Particle *c_p, Constraint_plane *c, double *dist, double *vec)
{
  int i;
  double c_dist_sqr,c_dist;
  
  c_dist_sqr=0.0;
  for(i=0;i<3;i++) {
    if(c->pos[i] >= 0) {
      vec[i] = c->pos[i] - ppos[i];
      c_dist_sqr += SQR(vec[i]);
    }else{
      vec[i] = 0.0;
      c_dist += SQR(vec[i]);
    }
  }
  c_dist = sqrt(c_dist_sqr);
  *dist = c_dist;

  
  for(i=0;i<3;i++) {
    vec[i] *= -1;
  }
}

MDINLINE void add_rod_force(Particle *p1, double ppos[3], Particle *c_p, Constraint_rod *c)
{
#ifdef ELECTROSTATICS
  int i;
  double fac, vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = ppos[i] - c->pos[i];
    c_dist_2 += SQR(vec[i]);
  }

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->lambda != 0.0) {
    fac = 2*coulomb.prefactor*c->lambda*p1->p.q/c_dist_2;
    p1->f.f[0]  += fac*vec[0];
    p1->f.f[1]  += fac*vec[1];
    c_p->f.f[0] -= fac*vec[0];
    c_p->f.f[1] -= fac*vec[1];
  }
#endif
}

MDINLINE double rod_energy(Particle *p1, double ppos[3], Particle *c_p, Constraint_rod *c)
{
#ifdef ELECTROSTATICS
  int i;
  double vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = ppos[i] - c->pos[i];
    c_dist_2 += SQR(vec[i]);
  }

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->lambda != 0.0)
    return coulomb.prefactor*p1->p.q*c->lambda*(-log(c_dist_2*SQR(box_l_i[2])) + 2*(M_LN2 - C_GAMMA));
#endif
  return 0;
}

MDINLINE void add_plate_force(Particle *p1, double ppos[3], Particle *c_p, Constraint_plate *c)
{
#ifdef ELECTROSTATICS
  double f;

  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->sigma != 0.0) {
    f = 2*M_PI*coulomb.prefactor*c->sigma*p1->p.q;
    if (ppos[2] < c->pos)
      f = -f;
    p1->f.f[2]  += f;
    c_p->f.f[2] -= f;
  }
#endif
}

MDINLINE double plate_energy(Particle *p1, double ppos[3], Particle *c_p, Constraint_plate *c)
{
#ifdef ELECTROSTATICS
  if (coulomb.prefactor != 0.0 && p1->p.q != 0.0 && c->sigma != 0.0)
    return -2*M_PI*coulomb.prefactor*c->sigma*p1->p.q*fabs(ppos[2] - c->pos);
#endif
  return 0;
}

//ER
MDINLINE void add_ext_magn_field_force(Particle *p1, Constraint_ext_magn_field *c)
{
#ifdef ROTATION
#ifdef DIPOLES
  p1->f.torque[0] += p1->r.dip[1]*c->ext_magn_field[2]-p1->r.dip[2]*c->ext_magn_field[1];
  p1->f.torque[1] += p1->r.dip[2]*c->ext_magn_field[0]-p1->r.dip[0]*c->ext_magn_field[2];
  p1->f.torque[2] += p1->r.dip[0]*c->ext_magn_field[1]-p1->r.dip[1]*c->ext_magn_field[0];
#endif
#endif
}

MDINLINE double ext_magn_field_energy(Particle *p1, Constraint_ext_magn_field *c)
{
#ifdef DIPOLES
     return -1.0 * scalar(c->ext_magn_field,p1->r.dip);
#endif
  return 0;
}
//end ER

MDINLINE void reflect_particle(Particle *p1, double *distance_vec, int reflecting) {
  double vec[3];
  double norm; 
  memcpy(vec, distance_vec, 3*sizeof(double));
//      printf("distance vector %f %f %f\n", vec[0], vec[1], vec[2]);
//      printf("reflecting particle %d at position %f %f %f new position %f %f %f\n", p1->p.identity, p1->r.p[0], p1->r.p[1], p1->r.p[2    ], p1->r.p[0]-2*vec[0], p1->r.p[1]-2*vec[1], p1->r.p[2]-2*vec[2]);
      reflection_happened = 1;
       p1->r.p[0] = p1->r.p[0]-2*vec[0];
       p1->r.p[1] = p1->r.p[1]-2*vec[1];
       p1->r.p[2] = p1->r.p[2]-2*vec[2];
       /* vec seams to be the vector that points from the wall to the particle*/
       /* now normalize it */ 
       norm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
       if ( reflecting==1 ) {
         vec[0] /= norm;
         vec[1] /= norm;
         vec[2] /= norm;
         /* calculating scalar product - reusing var norm */
         norm = vec[0] *  p1->m.v[0] + vec[1] * p1->m.v[1] + vec[2] * p1->m.v[2];
         /* now add twice the normal component to the velcity */
   //      printf("normal vector %f %f %f\n", vec[0], vec[1], vec[2]);
   //      printf("velocity before %f %f %f ", p1->m.v[0], p1->m.v[1], p1->m.v[2]);
          p1->m.v[0] = p1->m.v[0]-2*vec[0]*norm; /* norm is still the scalar product! */
          p1->m.v[1] = p1->m.v[1]-2*vec[1]*norm;
          p1->m.v[2] = p1->m.v[2]-2*vec[2]*norm;
  //       printf("velocity after %f %f %f\n", p1->m.v[0], p1->m.v[1], p1->m.v[2]);
       } else if (reflecting == 2) {
         /* if bounce back, invert velocity */
          p1->m.v[0] =-p1->m.v[0]; /* norm is still the scalar product! */
          p1->m.v[1] =-p1->m.v[1];
          p1->m.v[2] =-p1->m.v[2];

       }

}

MDINLINE void add_constraints_forces(Particle *p1)
{
  int n, j;
  double dist, vec[3], force[3], torque1[3], torque2[3];

  IA_parameters *ia_params;
  char *errtxt;
  double folded_pos[3];
  int img[3];

  /* fold the coordinate[2] of the particle */
  memcpy(folded_pos, p1->r.p, 3*sizeof(double));
  memcpy(img, p1->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);

  for(n=0;n<n_constraints;n++) {
    ia_params=get_ia_param(p1->p.type, (&constraints[n].part_rep)->p.type);
    dist=0.;
    for (j = 0; j < 3; j++) {
      force[j] = 0;
#ifdef ROTATION
      torque1[j] = torque2[j] = 0;
#endif
    }

    switch(constraints[n].type) {
    case CONSTRAINT_WAL: 
      if(checkIfInteraction(ia_params)) {
	calculate_wall_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.wal.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
    if(constraints[n].c.wal.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.wal.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{061 wall constraint %d violated by particle %d} ", n, p1->p.identity);
    }
	}
      }
      break;

    case CONSTRAINT_SPH:
      if(checkIfInteraction(ia_params)) {
	calculate_sphere_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.sph.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
    if(constraints[n].c.sph.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.sph.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{062 sphere constraint %d violated by particle %d} ", n, p1->p.identity);
    }
	}
      }
      break;
    
    case CONSTRAINT_CYL: 
      if(checkIfInteraction(ia_params)) {
	calculate_cylinder_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.cyl, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.cyl.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
    if(constraints[n].c.cyl.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.cyl.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
    }
        }
      }
      break;
	
    case CONSTRAINT_MAZE: 
      if(checkIfInteraction(ia_params)) {
	calculate_maze_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else if ( dist <= 0 && constraints[n].c.maze.penetrable == 1 ) {
	  if ( dist < 0 ) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,-1.0*dist,dist*dist, force,
				     torque1, torque2);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{064 maze constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_PORE: 
      if(checkIfInteraction(ia_params)) {
	calculate_pore_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.pore, &dist, vec); 
	if ( dist > 0 ) {
	  calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				     ia_params,vec,dist,dist*dist, force,
				     torque1, torque2);
	}
	else {
    if(constraints[n].c.pore.reflecting){
      reflect_particle(p1, &(vec[0]), constraints[n].c.pore.reflecting);
    } else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 pore constraint %d violated by particle %d} ", n, p1->p.identity);
        }
      }
      }
      break;
	
      /* electrostatic "constraints" */
    case CONSTRAINT_ROD:
      add_rod_force(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.rod);
      break;

    case CONSTRAINT_PLATE:
      add_plate_force(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.plate);
      break;
      
#ifdef MAGNETOSTATICS 
    case CONSTRAINT_EXT_MAGN_FIELD:
      add_ext_magn_field_force(p1, &constraints[n].c.emfield);
      break;
#endif
    
    case CONSTRAINT_PLANE:
     if(checkIfInteraction(ia_params)) {
	calculate_plane_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.plane, &dist, vec); 
	if (dist > 0) {
	    calc_non_bonded_pair_force(p1, &constraints[n].part_rep,
				       ia_params,vec,dist,dist*dist, force,
				       torque1, torque2);
#ifdef TUNABLE_SLIP
	    add_tunable_slip_pair_force(p1, &constraints[n].part_rep,ia_params,vec,dist,force);
#endif
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{063 plane constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
    }
    for (j = 0; j < 3; j++) {
      p1->f.f[j] += force[j];
      constraints[n].part_rep.f.f[j] -= force[j];
#ifdef ROTATION
      p1->f.torque[j] += torque1[j];
      constraints[n].part_rep.f.torque[j] += torque2[j];
#endif
    }
  }
}

MDINLINE double add_constraints_energy(Particle *p1)
{
  int n, type;
  double dist, vec[3];
  double nonbonded_en, coulomb_en, magnetic_en;
  IA_parameters *ia_params;
  char *errtxt;
  double folded_pos[3];
  int img[3];

  /* fold the coordinate[2] of the particle */
  memcpy(folded_pos, p1->r.p, 3*sizeof(double));
  memcpy(img, p1->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
  for(n=0;n<n_constraints;n++) { 
    ia_params = get_ia_param(p1->p.type, (&constraints[n].part_rep)->p.type);
    nonbonded_en = 0.;
    coulomb_en   = 0.;
    magnetic_en = 0.;

    dist=0.;
    switch(constraints[n].type) {
    case CONSTRAINT_WAL: 
      if(checkIfInteraction(ia_params)) {
	calculate_wall_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else if ( dist <= 0 && constraints[n].c.wal.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{065 wall constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_SPH: 
      if(checkIfInteraction(ia_params)) {
	calculate_sphere_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else if ( dist <= 0 && constraints[n].c.sph.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{066 sphere constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;
	
    case CONSTRAINT_CYL: 
      if(checkIfInteraction(ia_params)) {
	calculate_cylinder_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.cyl, &dist , vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);

	}
	else if ( dist <= 0 && constraints[n].c.cyl.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{067 cylinder constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_MAZE: 
      if(checkIfInteraction(ia_params)) {
	calculate_maze_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);
	}
	else if ( dist <= 0 && constraints[n].c.maze.penetrable == 1 ) {
	  if ( dist < 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, -1.0*dist, dist*dist);
	  }
	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{068 maze constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_PORE: 
      if(checkIfInteraction(ia_params)) {
	calculate_pore_dist(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.pore, &dist , vec); 
	if ( dist > 0 ) {
	  nonbonded_en = calc_non_bonded_pair_energy(p1, &constraints[n].part_rep,
						     ia_params, vec, dist, dist*dist);

	}
	else {
	  errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt, "{067 pore constraint %d violated by particle %d} ", n, p1->p.identity);
	}
      }
      break;

    case CONSTRAINT_ROD:
      coulomb_en = rod_energy(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.rod);
      break;

    case CONSTRAINT_PLATE:
      coulomb_en = plate_energy(p1, folded_pos, &constraints[n].part_rep, &constraints[n].c.plate);
      break;
    case CONSTRAINT_EXT_MAGN_FIELD:
      magnetic_en = ext_magn_field_energy(p1, &constraints[n].c.emfield);
      break;
    }

    if (energy.n_coulomb > 0)
      energy.coulomb[0] += coulomb_en;
    
    if (energy.n_dipolar > 0)
      energy.dipolar[0] += magnetic_en;

    type = (&constraints[n].part_rep)->p.type;
    if (type >= 0)
      *obsstat_nonbonded(&energy, p1->p.type, type) += nonbonded_en;
  }
  return 0.;
}

MDINLINE void init_constraint_forces()
{
  int n, i;
  
  for (n = 0; n < n_constraints; n++)
    for (i = 0; i < 3; i++)
      constraints[n].part_rep.f.f[i] = 0;
}
#endif

#endif
