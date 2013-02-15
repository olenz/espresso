/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
#ifndef _GRID_H
#define _GRID_H
/** \file grid.h   Domain decomposition for parallel computing.
 *
 *  The primary simulation box is divided into orthogonal rectangular
 *  subboxes which are assigned to the different nodes (or processes
 *  or threads if you want). This grid is described in \ref
 *  node_grid. Each node has a number \ref this_node and a position
 *  \ref node_pos in that grid. Each node has also 6 nearest neighbors
 *  \ref node_neighbors which are necessary for the communication
 *  between the nodes (see also \ref ghosts.c and \ref p3m.c for more
 *  details about the communication.
 *
 *  For the 6 directions \anchor directions we have the following convention:
 *
 *  \image html directions.gif "Convention for the order of the directions"
 *
 *  The Figure illustrates the direction convetion used for arrays
 *  with 6 (e.g. \ref node_neighbors, \ref #boundary) and 3 entries
 *  (e.g \ref node_grid, \ref box_l , \ref my_left,...).
 *  
 *
 *  For more information on the domain decomposition, see \ref grid.c "grid.c". 
*/
#include "lees_edwards.h"
#include "utils.h"
#include <limits.h>
#include "communication.h"
#include "errorhandling.h"

/** Macro that tests for a coordinate being periodic or not. */
#ifdef PARTIAL_PERIODIC
#define PERIODIC(coord) (periodic & (1L << coord))
#else
#define PERIODIC(coord) 1
#endif

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** The number of nodes in each spatial dimension. */
extern int node_grid[3];
/** position of node in node grid */
extern int node_pos[3];
/** the six nearest neighbors of a node in the node grid. */
extern int node_neighbors[6];
/** where to fold particles that leave local box in direction i. */
extern int boundary[6];
/** Flags for all three dimensions wether pbc are applied (default).
    The first three bits give the periodicity */
extern int periodic;

/** Simulation box dimensions. */ 
extern double box_l[3];
/** 1 / box dimensions. */ 
extern double box_l_i[3];
/** Smallest simulation box dimension (\ref box_l). 
    Remark: with PARTIAL_PERIODIC, only the periodic directions 
    are taken into account! */
extern double min_box_l;
/** Dimensions of the box a single node is responsible for. */ 
extern double local_box_l[3];
/** Smallest local simulation box dimension (\ref local_box_l).
    Remark: with PARTIAL_PERIODIC, only the periodic directions 
    are taken into account! */
extern double min_local_box_l;
/** Left (bottom, front) corner of this nodes local box. */ 
extern double my_left[3];
/** Right (top, back) corner of this nodes local box. */ 
extern double my_right[3];

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Make sure that the node grid is set, eventually
    determine one automatically. */
void setup_node_grid();

/** return wether node grid was set. */
int node_grid_is_set();

/** node mapping: array -> node. 
 *
 * \param node   rank of the node you want to know the position for.
 * \param pos    position of the node in node grid.        
*/
MDINLINE void map_node_array(int node, int pos[3])
{
  MPI_Cart_coords(comm_cart, node, 3, pos);
}

/** node mapping: node -> array. 
 *
 * \return      rank of the node at position pos.
 * \param pos   position of the node in node grid.        
*/
MDINLINE int map_array_node(int pos[3])
{
  int rank;
  MPI_Cart_rank(comm_cart, pos, &rank);
  return rank;
}

/** map a spatial position to the node grid */
int map_position_node_array(double pos[3]);

/** fill neighbor lists of node. 
 *
 * Calculates the numbers of the 6 nearest neighbors for a node and
 * stors them in \ref node_neighbors.
 *
 * \param node number of the node.  */
void calc_node_neighbors(int node);

/** called from \ref mpi_bcast_parameter . */
void grid_changed_n_nodes();

/** called from \ref mpi_bcast_parameter . */
void grid_changed_box_l();

/** Calculates the smallest box and local box dimensions for periodic
 * directions.  This is needed to check if the interaction ranges are
 * compatible with the box dimensions and the node grid.  
 * see also \ref box_l, \ref local_box_l, \ref min_box_l 
 * and \ref min_local_box_l.
 * Remark: In the apreiodic case min_box_l is set to 
 * 2 * \ref MAX_INTERACTION_RANGE . */
void calc_minimal_box_dimensions();

/** calculate most square 2d grid. */
void calc_2d_grid(int n, int grid[3]);

/** calculate 'best' mapping between a 2d and 3d grid.
 *  This we need for the communication from 3d domain decomposition 
 *  to 2d row decomposition. 
 *  The dimensions of the 2d grid are resorted, if necessary, in a way
 *  that they are multiples of the 3d grid dimensions.
 *  \param g3d      3d grid.
 *  \param g2d      2d grid.
 *  \param mult     factors between 3d and 2d grid dimensions
 *  \return         index of the row direction [0,1,2].
*/ 
int map_3don2d_grid(int g3d[3],int g2d[3], int mult[3]);

/** rescales the box in dimension 'dir' to the new value 'd_new', and rescales the particles accordingly */
void rescale_boxl(int dir, double d_new);

/** get the minimal distance vector of two vectors in the current bc.
    @param a the vector to subtract from
    @param b the vector to subtract
    @param res where to store the result
*/
MDINLINE void get_mi_vector(double res[3], double a[3], double b[3])
{
#ifdef LEES_EDWARDS
  int y_img_count;
  double delta_x;

  y_img_count = (int)floor(b[1]*box_l_i[1]) - (int)floor(a[1]*box_l_i[1]);
  delta_x     = y_img_count * lees_edwards_offset;   
  
  res[0]  = (a[0] - ( b[0] + delta_x) );
  res[0] -= dround(res[0]*box_l_i[0])*box_l[0];
  res[1]  = (a[1] - b[1]); 
  res[1] -= dround(res[1]*box_l_i[1])*box_l[1];
  res[2]  = (a[2] - b[2]); 
  res[2] -= dround(res[2]*box_l_i[2])*box_l[2];
#else
  int i;

  for(i=0;i<3;i++) {
    res[i] = a[i] - b[i];
#ifdef PARTIAL_PERIODIC
    if (PERIODIC(i))
#endif
      res[i] -= dround(res[i]*box_l_i[i])*box_l[i];
  }
#endif

}

/** fold a coordinate to primary simulation box.
    \param pos         the position...
    \param image_box   and the box
    \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
MDINLINE void fold_coordinate(double pos[3], int image_box[3], int dir)
{
  int tmp;
#ifdef PARTIAL_PERIODIC
  if (PERIODIC(dir))
#endif
    {
#ifdef LEES_EDWARDS
      if( dir == 0 ){
          int y_img_count;
          y_img_count   = (int)floor(pos[1]*box_l_i[1]);
          pos[0]       -= (lees_edwards_offset * y_img_count); 
      }
#endif

      image_box[dir] += (tmp = (int)floor(pos[dir]*box_l_i[dir]));
      pos[dir]        = pos[dir] - tmp*box_l[dir];    
      if(pos[dir] < 0 || pos[dir] >= box_l[dir]) {
	/* slow but safe */
	if (fabs(pos[dir]*box_l_i[dir]) >= INT_MAX/2) {
	  char *errtext = runtime_error(128 + ES_INTEGER_SPACE + ES_DOUBLE_SPACE);
	  ERROR_SPRINTF(errtext,"{086 particle coordinate out of range, pos = %g, image box = %d} ", pos[dir], image_box[dir]);
	  image_box[dir] = 0;
	  pos[dir] = 0;
	  return;
	}
      }
    }
}

/** fold particle coordinates to primary simulation box.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
MDINLINE void fold_position(double pos[3],int image_box[3])
{
  int i;
  for(i=0;i<3;i++)
    fold_coordinate(pos, image_box, i);
}

/** unfold coordinates to physical position.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
    afterwards.
*/
MDINLINE void unfold_position(double pos[3],int image_box[3])
{
#ifdef LEES_EDWARDS

  int y_img_count;
  y_img_count   = (int)floor( pos[1]*box_l_i[1] + image_box[1] );

  pos[0] = pos[0] + image_box[0]*box_l[0] + y_img_count*lees_edwards_offset;
  pos[1] = pos[1] + image_box[1]*box_l[1];
  pos[2] = pos[2] + image_box[2]*box_l[2];
  image_box[0] = image_box[1] = image_box[2] = 0;
#else
  int i;
  for(i=0;i<3;i++) {
    pos[i]       = pos[i] + image_box[i]*box_l[i];    
    image_box[i] = 0;
  }
#endif
}



/*************************************************************/
/** \name Distance calculations.  */
/*************************************************************/
/*@{*/

/** returns the distance between two position. 
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
MDINLINE double distance(double pos1[3], double pos2[3])
{
#ifndef LEES_EDWARDS
  return sqrt( SQR(pos1[0]-pos2[0]) + SQR(pos1[1]-pos2[1]) + SQR(pos1[2]-pos2[2]) );
#else
  int y_img_count;
  y_img_count = (int)floor(pos1[1]*box_l_i[1]) - (int)floor(pos2[1]*box_l_i[1]);
 
  return( sqrt(SQR(pos1[0]+y_img_count*lees_edwards_offset-pos2[0]) 
             + SQR(pos1[1]-pos2[1]) 
             + SQR(pos1[2]-pos2[2])) );
  
#endif
}

/** returns the distance between two positions squared.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
MDINLINE double distance2(double pos1[3], double pos2[3])
{
#ifndef LEES_EDWARDS
  return SQR(pos1[0]-pos2[0]) + SQR(pos1[1]-pos2[1]) + SQR(pos1[2]-pos2[2]);
#else
  int y_img_count;
  y_img_count = (int)floor(pos1[1]*box_l_i[1]) - (int)floor(pos2[1]*box_l_i[1]);
 
  return( SQR(pos1[0]+y_img_count*lees_edwards_offset-pos2[0]) 
             + SQR(pos1[1]-pos2[1]) 
             + SQR(pos1[2]-pos2[2]) );
  
#endif
}

/** Returns the distance between two positions squared and stores the
    distance vector pos1-pos2 in vec.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 *  \param vec  vecotr pos1-pos2.
 *  \return distance squared
*/
MDINLINE double distance2vec(double pos1[3], double pos2[3], double vec[3])
{
#ifdef LEES_EDWARDS
  int y_img_count;
  y_img_count = (int)floor(pos1[1]*box_l_i[1]) - (int)floor(pos2[1]*box_l_i[1]);
  vec[0]      = pos1[0]+y_img_count*lees_edwards_offset - pos2[0];
  while(vec[0] >  0.5*box_l[0]){vec[0]-=box_l[0];}
  while(vec[0] < -0.5*box_l[0]){vec[0]+=box_l[0];}
#else
  vec[0] = pos1[0]-pos2[0];
#endif
  vec[1] = pos1[1]-pos2[1];
  vec[2] = pos1[2]-pos2[2];
  return SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]);
}

/** returns the distance between the unfolded coordintes of two particles. 
 *  \param pos1       Position of particle one.
 *  \param image_box1 simulation box index of particle one .
 *  \param pos2       Position of particle two.
 *  \param image_box2 simulation box index of particle two .
 *  \param box_l      size of simulation box.
*/
MDINLINE double unfolded_distance(double pos1[3], int image_box1[3], 
                  double pos2[3], int image_box2[3], double box_l[3])
{
  double dist;
  double lpos1[3],lpos2[3];


  /*unrolling the loop so can neatly add Lees-Edwards: 
   *compiler probably unrolls anyway*/
  lpos1[0]  = pos1[0] + image_box1[0]*box_l[0];
  lpos2[0]  = pos2[0] + image_box2[0]*box_l[0];
#ifdef LEES_EDWARDS
  lpos1[0] += image_box1[1] * lees_edwards_offset;
  lpos2[0] += image_box2[1] * lees_edwards_offset;
#endif
  dist      = SQR(lpos1[0]-lpos2[0]);

  lpos1[1]  = pos1[1] + image_box1[1]*box_l[1];
  lpos2[1]  = pos2[1] + image_box2[1]*box_l[1];
  dist     += SQR(lpos1[1]-lpos2[1]);

  lpos1[2]  = pos1[2] + image_box1[2]*box_l[2];
  lpos2[2]  = pos2[2] + image_box2[2]*box_l[2];
  dist     += SQR(lpos1[2]-lpos2[2]);

  return sqrt(dist);
}
/*@}*/


/*@}*/
#endif
