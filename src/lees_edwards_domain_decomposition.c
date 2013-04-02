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
/** \file lees_edwards_domain_decomposition.c
 *
 *  This file contains modifications to the domain decomposition functions 
 *  which are specific to Lees-Edwards boundary conditions.
 *  See also \ref lees_edwards_domain_decomposition.h and  \ref domain_decomposition.h
 */


#include "domain_decomposition.h"
#include "lees_edwards_domain_decomposition.h"
#include "errorhandling.h"
#include "forces.h"
#include "pressure.h"
#include "energy.h"
#include "constraint.h"
#include "initialize.h"
#include "lees_edwards.h"
#include "grid.h"

#ifdef LEES_EDWARDS

/** Init cell interactions for the Lees-Edwards cell system.
 * initializes the interacting neighbor cell list of a cell. The
 * created list of interacting neighbor cells is used by the verlet
 * algorithm (see verlet.c) to build the verlet lists.
 */
void le_dd_init_cell_interactions()
{
  int m,n,o,p,q,r,ind1,ind2,c_cnt=0,n_cnt;
  int my_ncount, extra_cells = 0;
 

  /* initialize cell neighbor structures */
  dd.cell_inter = (IA_Neighbor_List *) realloc(dd.cell_inter,local_cells.n*sizeof(IA_Neighbor_List));
  for(m=0; m<local_cells.n; m++) { 
    dd.cell_inter[m].nList = NULL; 
    dd.cell_inter[m].n_neighbors=0; 
  }

  /* loop over non-ghost cells */
  for(o=1; o<=dd.cell_grid[2]; o++) \
    for(n=1; n<=dd.cell_grid[1]; n++) \
      for(m=1; m<=dd.cell_grid[0]; m++) {

    /* plenty for most cases */
    dd.cell_inter[c_cnt].nList = (IA_Neighbor *) realloc(dd.cell_inter[c_cnt].nList, 14*sizeof(IA_Neighbor));
    
    n_cnt=0;
    ind1 = get_linear_index(m,n,o,dd.ghost_cell_grid);

    /* loop all 'conventional' neighbor cells */
    for(p=o-1; p<=o+1; p++)	{      /*z-loop*/
      for(q=n-1; q<=n+1; q++) {    /*y-loop*/
        for(r=m-1; r<=m+1; r++) {  /*x-loop*/

            /* special conditions at the top and bottom y-surfaces are dealt with later */
            if( q == dd.cell_grid[1] + 1 && boundary[3] == -1)
                continue;
            if( q == 0 && boundary[2] == 1)
                continue;

            /* skip the un-ghosted corners: imaging in x for these happens at interaction time. */
            /* also skip ghost interactions which are duped in the LE-extra ghosts */
            if( q == 0 && boundary[2] == 1 )           
                if( r == 0 || r == dd.cell_grid[0] + 1 )
                    continue;
            if( q == dd.cell_grid[1] + 1 && boundary[3] == -1 )           
                if( r == 0 || r == dd.cell_grid[0] + 1 )
                    continue;
            
	        ind2 = get_linear_index(r,q,p,dd.ghost_cell_grid);
	        if(ind2 >= ind1) {
	            dd.cell_inter[c_cnt].nList[n_cnt].cell_ind = ind2;
	            dd.cell_inter[c_cnt].nList[n_cnt].pList    = &cells[ind2];
	            init_pairList(&dd.cell_inter[c_cnt].nList[n_cnt].vList);
	
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[0] = my_left[0] + r * dd.cell_size[0];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[1] = my_left[1] + q * dd.cell_size[1];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[2] = my_left[2] + p * dd.cell_size[2];

	            n_cnt++;
   

           }
        }
      }
	}

    /* special conditions at the top y-surface */
    if( n == dd.cell_grid[1] && boundary[3] == -1){ 
        q = n + 1;

        /* ask for more neighbor cells */
        extra_cells = 3 * dd.cell_grid[0]; /*strip at fixed y, 3 wide in z. This is a slight over-alloc. */
        dd.cell_inter[c_cnt].nList = 
           (IA_Neighbor *)realloc(dd.cell_inter[c_cnt].nList,(n_cnt+extra_cells)*sizeof(IA_Neighbor));

        for(p=o-1; p<=o+1; p++)	{ /* usual loop over z */
            for( r=1; r <= dd.cell_grid[0]; r++){ /* ANY x could be a neighbor */

//                if( abs(r-m) <= 1) continue; /* skip dupes */
//                if( r == 1 && m == dd.cell_grid[0] ) continue;
//                if( r == dd.cell_grid[0] && m == 1 ) continue;

	            ind2 = get_linear_index(r,q,p,dd.ghost_cell_grid);
	            if(ind2 >= ind1) {
	                dd.cell_inter[c_cnt].nList[n_cnt].cell_ind = ind2;
	                dd.cell_inter[c_cnt].nList[n_cnt].pList    = &cells[ind2];
	                init_pairList(&dd.cell_inter[c_cnt].nList[n_cnt].vList);
	
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[0] = my_left[0] + r * dd.cell_size[0];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[1] = my_left[1] + q * dd.cell_size[1];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[2] = my_left[2] + p * dd.cell_size[2];

	                n_cnt++;
	            }
            }
        }
    }
    /* at the bottom y-surface */
    if( n == 1 && boundary[2] == 1){ 
        q = 0;
        /* ask for more neighbor cells */
        extra_cells = 3 * dd.cell_grid[0];
        dd.cell_inter[c_cnt].nList = 
           (IA_Neighbor *)realloc(dd.cell_inter[c_cnt].nList,(n_cnt+extra_cells)*sizeof(IA_Neighbor));
        for(p=o-1; p<=o+1; p++)	{ /* usual loop over z */
            for( r=1; r <= dd.cell_grid[0]; r++){ /* ANY x could be a neighbor */

//                if( abs(r-m) <= 1) continue; /* skip dupes */
//                if( r == 1 && m == dd.cell_grid[0] ) continue;
//                if( r == dd.cell_grid[0] && m == 1 ) continue;

	            ind2 = get_linear_index(r,q,p,dd.ghost_cell_grid);
	            if(ind2 >= ind1) {
	                dd.cell_inter[c_cnt].nList[n_cnt].cell_ind = ind2;
	                dd.cell_inter[c_cnt].nList[n_cnt].pList    = &cells[ind2];
	                init_pairList(&dd.cell_inter[c_cnt].nList[n_cnt].vList);

    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[0] = my_left[0] + r * dd.cell_size[0];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[1] = my_left[1] + q * dd.cell_size[1];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[2] = my_left[2] + p * dd.cell_size[2];


	                n_cnt++;
	            }
            }
        }
    }

#if 1
    /* Also need ghost cells from other nodes, if running in parallel.*/
    if( n == dd.cell_grid[1] && boundary[3] == -1 ){

        extra_cells  = dd.ghost_cell_grid_le_extra[0] * dd.ghost_cell_grid_le_extra[2];

        dd.cell_inter[c_cnt].nList = 
           (IA_Neighbor *)realloc(dd.cell_inter[c_cnt].nList, 
                                  (extra_cells+n_cnt)*sizeof(IA_Neighbor));

        /* ghosted lower surface */
        q = dd.ghost_cell_grid_le_extra[1] - 1; /* this q (y-index to the extra grid) is 0 or 1 */
        for( p = o-1; p <= o+1; p++ ){          /* same as the usual loop over z */
            for( r = 0; r < dd.ghost_cell_grid_le_extra[0]; r++ ){ /* same as the usual loop over x */

                /* these extra cells correspond to ghost cells which already exist on other nodes */
                int pp, qq, rr, neighbor_pos_x;
        
                neighbor_pos_x = r / dd.cell_grid[0];
                if( neighbor_pos_x >= node_pos[0] ){neighbor_pos_x++;}

                /* enforce no two-way links simply, based only on x-index,
                   because x-index is never the same for an le_extra and its parent */
                if( neighbor_pos_x >= node_pos[0] ) continue;


                rr = neighbor_pos_x * dd.cell_grid[0] + (r % dd.cell_grid[0]) + 1; 
                qq = (1-q) *  dd.ghost_cell_grid[1];
                pp = p;

                ind2  = get_linear_index(r,q,p,dd.ghost_cell_grid_le_extra);
                ind2 += (n_cells - n_le_extra_cells);

	            dd.cell_inter[c_cnt].nList[n_cnt].cell_ind = ind2;
	            dd.cell_inter[c_cnt].nList[n_cnt].pList    = &cells[ind2];
	            init_pairList(&dd.cell_inter[c_cnt].nList[n_cnt].vList);

    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[0] = rr * dd.cell_size[0];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[1] = my_left[1] + (qq-1) * dd.cell_size[1];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[2] = pp * dd.cell_size[2];

	            n_cnt++;
            }
        }
    }
    if( n == 1 && boundary[2] == 1 ){

        extra_cells  = dd.ghost_cell_grid_le_extra[0] * dd.ghost_cell_grid_le_extra[2];

        dd.cell_inter[c_cnt].nList = 
           (IA_Neighbor *)realloc(dd.cell_inter[c_cnt].nList, 
                                  (extra_cells+n_cnt)*sizeof(IA_Neighbor));

        /* ghosted upper surface */
        q = 0;
        for( p = o-1; p <= o+1; p++ ){
            for( r = 0; r < dd.ghost_cell_grid_le_extra[0]; r++ ){

                /* these extra cells correspond to ghost cells which already exist on other nodes */
                int pp, qq, rr, neighbor_pos_x;
        
                neighbor_pos_x = r / dd.cell_grid[0];
                if( neighbor_pos_x >= node_pos[0] ){neighbor_pos_x++;}

                /* enforce no two-way links simply, based only on x-index,
                   because x-index is never the same for le_extras */
                if( neighbor_pos_x >= node_pos[0] ) continue;

                rr = neighbor_pos_x * dd.cell_grid[0] + (r % dd.cell_grid[0]) + 1; 
                qq = (1-q) *  dd.ghost_cell_grid[1];
                pp = p;

	            ind2  = get_linear_index(r,q,p,dd.ghost_cell_grid_le_extra);
                ind2 += (n_cells - n_le_extra_cells);
	            dd.cell_inter[c_cnt].nList[n_cnt].cell_ind = ind2;
	            dd.cell_inter[c_cnt].nList[n_cnt].pList    = &cells[ind2];
	            init_pairList(&dd.cell_inter[c_cnt].nList[n_cnt].vList);

    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[0] = rr * dd.cell_size[0];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[1] = 0.0;
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[2] = pp * dd.cell_size[2];

	            n_cnt++;
            }
        }
    }
#endif

    dd.cell_inter[c_cnt].n_neighbors = n_cnt; 
    c_cnt++;
  }


  FILE *cells_fp;
  char cLogName[64];
  int  c,nn,this_n;
  double myPos[3], nPos[3];
  sprintf(cLogName, "cells_map%i.dat", this_node);
  cells_fp = fopen(cLogName,"w");

  n_le_extra_cells =  ( dd.ghost_cell_grid_le_extra[0] *
                        dd.ghost_cell_grid_le_extra[1] * 
                        dd.ghost_cell_grid_le_extra[2] );

  for(c=0;c<c_cnt;c++){
     myPos[0] = my_left[0] + dd.cell_size[0] * ( 1 + c % dd.cell_grid[0] );  
     myPos[1] = my_left[1] + dd.cell_size[1] * ( 1 + (c / dd.cell_grid[0]) % dd.cell_grid[1]);  
     myPos[2] = my_left[2] + dd.cell_size[2] * ( 1 + (c / (dd.cell_grid[0] * dd.cell_grid[1])));  

     for(nn=0;nn<dd.cell_inter[c].n_neighbors;nn++){
        
        this_n = dd.cell_inter[c].nList[nn].cell_ind;


        fprintf(cells_fp,"%i %i %i %f %f %f %f %f %f\n",c,nn,this_n,
            myPos[0], myPos[1], myPos[2], 
            dd.cell_inter[c].nList[nn].my_pos[0], 
            dd.cell_inter[c].nList[nn].my_pos[1], 
            dd.cell_inter[c].nList[nn].my_pos[2]);
          
     }
  }  
  fclose(cells_fp);

}


/** Fill a communication cell pointer list. Fill the cell pointers of
    all cells which are inside the grid of additional ghost cells which need to
    be tracked for Lees-Edwards.
    \param part_lists      List of cell pointers to store the result.
    \param cell_array_base position of the first 'extra' cells in the cell array
    \param lc              lower left corner of the subgrid.
    \param hc              high up corner of the subgrid.
 */
int le_dd_fill_comm_cell_lists(Cell **part_lists, Cell *cell_array_base, int lc[3], int hc[3])
{
  int i,m,n,o,c=0;

  /* sanity check */
  for(i=0; i<3; i++) {
    if(lc[i]<0     || lc[i] >= dd.ghost_cell_grid_le_extra[i]) return 0;
    if(hc[i]<lc[i] || hc[i] >= dd.ghost_cell_grid_le_extra[i]) return 0;
  }

  for(o=lc[0]; o<=hc[0]; o++) 
    for(n=lc[1]; n<=hc[1]; n++) 
      for(m=lc[2]; m<=hc[2]; m++) {
	i = get_linear_index(o,n,m,dd.ghost_cell_grid_le_extra);
	//CELL_TRACE(fprintf(stderr,"%d: dd_fill_comm_cell_list: added ghost cell %d based on crds %i %i %i\n",this_node,i,o,n,m));
	part_lists[c] = &cell_array_base[i];
	c++;
      }

  return c;
}

/** update the 'shift' member of those GhostCommunicators, which use
    that value to speed up the folding process of its ghost members
    (see \ref dd_prepare_comm for the original), i.e. all which have
    GHOSTTRANS_POSSHFTD or'd into 'data_parts' upon execution of \ref
    dd_prepare_comm. */
void le_dd_update_communicators_w_boxl()
{
  int cnt=0;
  int neighbor_index, neighbor_rank, dir, lr, i;


  for( i = 0; i < my_neighbor_count; i++ ){
     neighbor_index = i;

     /* find some geometry about this comm */
     dir = 1;
     if( i == 0 || i == 1 )      { dir = 0; }     
     else if ( i == 2 || i == 3 ){ dir = 2; neighbor_index += 2;}
     else if ( i == 4 || i == 5 ){ neighbor_index -= 2;}

     lr = node_neighbor_lr[neighbor_index];
     neighbor_rank = node_neighbors[neighbor_index];

     if( neighbor_rank == this_node ) { 
	    /* prepare folding of ghost positions */
	    if(node_neighbor_wrap[neighbor_index] != 0) {
	      cell_structure.exchange_ghosts_comm.comm[cnt].shift[dir]  = node_neighbor_wrap[neighbor_index]*box_l[dir]; 
	      cell_structure.update_ghost_pos_comm.comm[cnt].shift[dir] = node_neighbor_wrap[neighbor_index]*box_l[dir];
          if( dir == 1 ){
           cell_structure.exchange_ghosts_comm.comm[cnt].shift[0]   = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
           cell_structure.update_ghost_pos_comm.comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
          }
	    }
        GHOST_TRACE(fprintf(stderr, "%d: Ghost_shift z %i, dir %i, wrap:%i: %f,%f,%f  %f,%f,%f \n",this_node,cnt,dir,
                node_neighbor_wrap[neighbor_index],
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[0], 
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[1],
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[2],
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[0], 
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[1],
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[2]);)
	    cnt++;	  
     }
     else {
	 /* send/recv loop: easiest to just set the shift for both of them, only sends will actually be shifted. */
	for(int send_rec=0; send_rec<2; send_rec++) {  

	      /* prepare folding of ghost positions */
	      if(node_neighbor_wrap[neighbor_index] != 0) {
		cell_structure.exchange_ghosts_comm.comm[cnt].shift[dir]  = node_neighbor_wrap[neighbor_index]*box_l[dir];
		cell_structure.update_ghost_pos_comm.comm[cnt].shift[dir] = node_neighbor_wrap[neighbor_index]*box_l[dir];
            if( dir == 1 ){
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[0]   = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
            }
	      }

        GHOST_TRACE(fprintf(stderr, "%d: Ghost_shift %i, dir %i, wrap:%i: %f,%f,%f  %f,%f,%f \n",this_node,cnt,dir,
                node_neighbor_wrap[neighbor_index],
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[0], 
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[1],
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[2],
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[0], 
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[1],
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[2]);)
	      cnt++;
	    }

    }
  }
}

/** Create communicators for cell structure domain decomposition. (see \ref GhostCommunicator) */

void  le_dd_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  int dir,lr,i,cnt, num, n_comm_cells[3], send_rec;
  int lc[3],hc[3], neighbor_index;
  int have_le[2];

  have_le[0] = 0;
  have_le[1] = 0;
   
  /* prepare communicator... need 2 comms for non-trivial exchange, 1 comm for local exchange */
  num = my_neighbor_count * 2;
  for(i = 0; i < my_neighbor_count; i++) {
        if( node_neighbors[i] == this_node ){
            num--;
        }
  }

  CELL_TRACE(fprintf(stderr,"%d Create Communicator: prep_comm data_parts %d num_neighbor_comms: %d\n",this_node,data_parts,num));
  prepare_comm(comm, data_parts, num);

  /* number of cells owned by this node to communicate in each direction:
     note that  dd.ghost_cell_grid[1] doesn't get a mention. */
  n_comm_cells[0] = dd.cell_grid[1]       * dd.cell_grid[2];
  n_comm_cells[1] = dd.ghost_cell_grid[2] * dd.ghost_cell_grid[0]; /* order of comms is x-z-y */
  n_comm_cells[2] = dd.ghost_cell_grid[0] * dd.cell_grid[1];

  cnt=0;

  CELL_TRACE(fprintf(stderr,"%d: neighbours:", this_node));  
  CELL_TRACE(for(i = 0; i < my_neighbor_count; i++)fprintf(stderr," %d",node_neighbors[i]);) 
  CELL_TRACE(fprintf(stderr,"\n"));  

  /* loop over neighbours */
  for(i = 0; i < my_neighbor_count; i++) {
  
    int buf, neighbor_coords[3], neighbor_rank;
    
    /* must loop on order x-z-y or (z-x-y): y imaging must carry ghost cells already genned in the 
       plane perpendicular to y. */
    neighbor_index = i;
    if( i == 0 || i == 1 ){
        dir = 0; /* sharing an x-plane of the cube */
        lc[1] = 1;
        lc[2] = 1;
        hc[1] = dd.cell_grid[1];
        hc[2] = dd.cell_grid[2];
    }
    else if( i == 2 || i == 3 ){ 
        dir   = 2; /* sharing a z-plane of the cube */
        neighbor_index += 2;
        lc[0] = 0;
        lc[1] = 1;
        hc[0] = dd.cell_grid[0] + 1;
        hc[1] = dd.cell_grid[1];
    } 
    else  { /*( i >= 4 ) */ 
        dir   = 1; /* sharing a y-plane of the cube */
        if( i == 4 || i == 5 ) {
            neighbor_index -= 2;
        }
        lc[0] = 0;
        lc[2] = 0;
        hc[0] = dd.cell_grid[0] + 1;
        hc[2] = dd.cell_grid[2] + 1;
    }
    lr            = node_neighbor_lr[neighbor_index]; 


    /* find out where this neighbor is */
    neighbor_rank = node_neighbors[neighbor_index];
    map_node_array(neighbor_rank, neighbor_coords);
 
    /* transfers along a given axis have lc[axis] == hc[axis] == fixed: you are sending/receiving
     * a layer of cells defined by that value of x,y or z. 
     *
     * For example, a receive from the cell above you in y would the range:
     * (??, dd.cell_grid[1] + 1, ??) for both lc and hc
     * while a send to the cell above you would be the range:
     * (??, dd.cell_grid[1], ??)
     * because the arriving cells go into the "ghost" plane, while the sent ones are 'real'.
     *
     * For example, a (send to /receive from) the cell below in y would be:
     * send: (??, 1, ??)
     * recv: (??, 0, ??)
     * because the arriving cells go into the "ghost" plane, while the sent ones are 'real'.
     */
     
    
    if( neighbor_rank == this_node ) { /* if copying cells on a single node, then comms are a bit simpler... */
        
        comm->comm[cnt].type          = GHOST_LOCL;
	    comm->comm[cnt].node          = this_node;


	    /* prepare folding of ghost positions */
        data_parts |= GHOSTTRANS_POSSHFTD;
        comm->comm[cnt].shift[dir] = node_neighbor_wrap[neighbor_index]*box_l[dir];  

        /* special case: sending through a y-wrap */
        if( dir == 1 ){
            comm->comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;  

            /* we don't send corner-ghosts through a LE y-wrap */
            lc[0]                     = 1;
            hc[0]                     = dd.cell_grid[0];
             
            n_comm_cells[1] = dd.cell_grid[0] * dd.ghost_cell_grid[2];

        }

        /*if top plane */
        if( lr == 1 ) {
                    lc[dir] = dd.cell_grid[dir];
                    hc[dir] = dd.cell_grid[dir];
        }
        /*bottom plane*/
        else {
                    lc[dir] = 1;
                    hc[dir] = 1;
        }

        
	    /* Buffer has to contain both Send and Recv cells -> factor 2 */
	    comm->comm[cnt].part_lists    = malloc(2*n_comm_cells[dir]*sizeof(ParticleList *));
	    comm->comm[cnt].n_part_lists  = 2*n_comm_cells[dir];


        if( dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc) != comm->comm[cnt].n_part_lists / 2){
            fprintf(stderr, "ERROR failed loading comm cell lists for node-local send.\n");
            exit( 1 );
        }
 

	    CELL_TRACE(fprintf(stderr,"%d: prep_comm %d copy %d from grid (%d,%d,%d)-(%d,%d,%d) shift: %f-%f-%f\n",this_node,cnt, 
                   comm->comm[cnt].n_part_lists/2,
			       lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]));
	    
        /* fill recv comm cells */
        if( lr == 1 ) {
                    lc[dir] = 0;
                    hc[dir] = 0;
        }
        else {
                    lc[dir] = dd.cell_grid[dir] + 1;
                    hc[dir] = dd.cell_grid[dir] + 1;
        }
        
	    /* place receive cells after send cells */
        if( dd_fill_comm_cell_lists(&comm->comm[cnt].part_lists[n_comm_cells[dir]],lc,hc) != comm->comm[cnt].n_part_lists / 2){
            fprintf(stderr, "ERROR failed loading comm cell lists for node-local receive.\n");
            exit( 1 );
        }
	    CELL_TRACE(fprintf(stderr,"%d: prep_comm %d copy %d into grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,comm->comm[cnt].n_part_lists/2,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	    cnt++;
         

    }else{
        /* send/recv loop */
        for(send_rec=0; send_rec<2; send_rec++) {  
        
            comm->comm[cnt].node          = neighbor_rank;
 
            if( (send_rec == 0  && neighbor_rank > this_node) 
              ||(send_rec == 1  && neighbor_rank < this_node)){
                comm->comm[cnt].type          = GHOST_SEND;
          
                /* prepare fold-on-send of ghost positions */
                if(neighbor_index < 6 ){
                    if(node_neighbor_wrap[neighbor_index] != 0 ){                    
                        data_parts |= GHOSTTRANS_POSSHFTD;
                        comm->comm[cnt].shift[dir] = node_neighbor_wrap[neighbor_index]*box_l[dir];
                        if( dir == 1 ){
                            comm->comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
                        }
                    }
                }

                /*choose the plane to send*/
                if( lr == 1 ) {
                        lc[dir] = dd.cell_grid[dir];
                        hc[dir] = dd.cell_grid[dir];
                }else {
                        lc[dir] = 1;
                        hc[dir] = 1;
                }
 

                if( neighbor_index >= 6){ /* send to extra Lees-Edwards ghost cells */
                    int node_off = neighbor_coords[0] - node_pos[0];
                    int n_lists;

                    data_parts |= GHOSTTRANS_POSSHFTD;
                    comm->comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index] * lees_edwards_offset;
                    comm->comm[cnt].shift[1]  = node_neighbor_wrap[neighbor_index] * box_l[1]; 

                    /* need a different indexing for send to LE ghosts */              
                    lc[0] = 1;
                    hc[0] = dd.cell_grid[0];

                    if( hc[0] >= 0 && lc[0] >= 0 ){
                        /* number of lists is (x-range * z-range * 1) */
                        n_lists = (hc[0]-lc[0]+1) * (dd.ghost_cell_grid[2]);
                    }else{
                        n_lists = 0;
                    }
                    if( n_lists > 0 ){                    

                        comm->comm[cnt].part_lists    = malloc(n_lists*sizeof(ParticleList *));
                        comm->comm[cnt].n_part_lists  = n_lists;

                        dd_fill_comm_cell_lists(comm->comm[cnt].part_lists, lc, hc);

                    }else{            
                        comm->comm[cnt].part_lists    = NULL;
                        comm->comm[cnt].n_part_lists  = 0;
                    }
                }else{


                    /* we don't send corner-ghosts through a LE y-wrap */
                    if( dir == 1 ){
                        if(  node_neighbor_wrap[neighbor_index] != 0 ){
                            lc[0]     = 1;
                            hc[0]     = dd.cell_grid[0];
                            n_comm_cells[1] = dd.cell_grid[0] * dd.ghost_cell_grid[2];
                        }else{
                            n_comm_cells[1] = (1+hc[0]-lc[0]) * dd.ghost_cell_grid[2];
                        }
                    }



                    comm->comm[cnt].part_lists    = malloc(n_comm_cells[dir]*sizeof(ParticleList *));
                    comm->comm[cnt].n_part_lists  = n_comm_cells[dir];


                    if( dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc) != comm->comm[cnt].n_part_lists ){
                            fprintf(stderr, "%i: ERROR failed loading comm cell lists for send lr=%i in dir %i wanted %i and got %i.\n",
                               this_node, lr, dir, comm->comm[cnt].n_part_lists, dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc));

                            fprintf(stderr,"%d: grid was (%d,%d,%d)-(%d,%d,%d)\n",this_node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]);
                            fprintf(stderr,"%d: full grid was (%d,%d,%d)\n",this_node,dd.ghost_cell_grid[0],dd.ghost_cell_grid[1],dd.ghost_cell_grid[2]);
                            exit( 1 );
                    }
                }



                int pc_index, part_count = 0;
                CELL_TRACE(for( pc_index = 0; pc_index < comm->comm[cnt].n_part_lists; pc_index++) part_count+=comm->comm[cnt].part_lists[pc_index]->n;);
                CELL_TRACE(fprintf(stderr,"%d: prep_comm %d send %d lists %d parts to   node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d shift (%f,%f,%f)\n",this_node,cnt,comm->comm[cnt].n_part_lists,part_count,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],dir,comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]));

            }else {
                comm->comm[cnt].type          = GHOST_RECV;
          
                if( lr == 1 ) {
                    /* if top plane */
                    lc[dir] = dd.cell_grid[dir] + 1;
                    hc[dir] = dd.cell_grid[dir] + 1;
                }
                else {
                /* bottom plane */
                    lc[dir] = 0;
                    hc[dir] = 0;
                }
          
                if( neighbor_index < 6 ){

                    /* we don't send corner-ghosts through a LE y-wrap */
                    if( dir == 1 ){
                        if( node_neighbor_wrap[neighbor_index] != 0 ){               
                            lc[0]     = 1;
                            hc[0]     = dd.cell_grid[0];
                        }
                        n_comm_cells[1] = (1+hc[0]-lc[0]) * dd.ghost_cell_grid[2];
                    }
                    comm->comm[cnt].part_lists    = malloc(n_comm_cells[dir]*sizeof(ParticleList *));
                    comm->comm[cnt].n_part_lists  = n_comm_cells[dir];
                    dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);

                }else{ /* receive into Lees-Edwards extra ghost cells */
                    int node_off = neighbor_coords[0] - node_pos[0];
                    int n_lists, receive_from_this;

                    /* the neighbor is to my right, then already have some*/
                    if( node_off > 0){
                        lc[0] = (neighbor_coords[0] - 1) * dd.cell_grid[0];
                    }else{
                        lc[0] =  neighbor_coords[0] * dd.cell_grid[0];
                    }
                    hc[0] = lc[0] + dd.cell_grid[0] - 1;


                    /* y and z grid sections fixed more as normal */
                    /* there are max of two layers in y */
                    if( lr == 0 ){
                        lc[1] = 0;
                        hc[1] = 0;
                    }else{
                        lc[1] = dd.ghost_cell_grid_le_extra[1] - 1;
                        hc[1] = dd.ghost_cell_grid_le_extra[1] - 1;
                    }                    
                      lc[2] = 0;
                    hc[2] = dd.ghost_cell_grid_le_extra[2] - 1;


                    /* number of lists is (x-range * z-range * 1) */
                    if( hc[0] >= lc[0] && lc[0] >= 0 ){ 
                        n_lists = (hc[0]-lc[0]+1) * (hc[2]-lc[2]+1);
                    }else{
                        n_lists = 0;
                    }

                    if( n_lists > 0 ){                    
                        comm->comm[cnt].part_lists    = malloc(n_lists*sizeof(ParticleList *));
                        comm->comm[cnt].n_part_lists  = n_lists;
                        
                        /* need to address into the Lees-Edwards 'extra' ghost array */
                        le_dd_fill_comm_cell_lists(comm->comm[cnt].part_lists, 
                                               &cells[n_cells - n_le_extra_cells], 
                                               lc, hc);
                      
                    }else{        
                        comm->comm[cnt].part_lists    = NULL;
                        comm->comm[cnt].n_part_lists  = 0;
                    }
                 }
          
                CELL_TRACE(fprintf(stderr,"%d: prep_comm %d recv %d lists from node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d\n",this_node,cnt,comm->comm[cnt].n_part_lists,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2], dir));
            }
            cnt++;
        }
      }
    }

}




/************************************************************/
void  le_dd_exchange_and_sort_particles(int global_flag)
{
  int           neighbor_index, dir, c, p, i, finished=0, count = 0;
  ParticleList *cell,*sort_cell, send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;
  Particle     *part;

  CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles(%d):\n",this_node,global_flag));

  init_particlelist(&send_buf_l);
  init_particlelist(&send_buf_r);
  init_particlelist(&recv_buf_l);
  init_particlelist(&recv_buf_r);

  while(finished == 0 ) {
    finished=1; count++;

    /* loop over x,y,z */  
    for( dir = 0; dir < 2; dir++ ){    
      if(node_grid[dir] > 1) {

	/* Communicate particles that have left the node domain */
	/* particle loop */
	for(c=0; c<local_cells.n; c++) {
	  cell = local_cells.cell[c];
	  for (p = 0; p < cell->n; p++) {
	    part = &cell->part[p];
       
	    /* Move particles to the left side */
	    if(part->r.p[dir] - my_left[dir] < -ROUND_ERROR_PREC*box_l[dir]) {
		  CELL_TRACE(fprintf(stderr,"%d: dd_ex_and_sort_p: send part left %d\n",this_node,part->p.identity));
		  local_particles[part->p.identity] = NULL;
		  move_indexed_particle(&send_buf_l, cell, p);
		  if(p < cell->n) p--;
		}
	    /* Move particles to the right side */
	    else if(part->r.p[dir] - my_right[dir] >= ROUND_ERROR_PREC*box_l[dir]) {
		  CELL_TRACE(fprintf(stderr,"%d: dd_ex_and_sort_p: send part right %d\n",this_node,part->p.identity));
		  local_particles[part->p.identity] = NULL;
		  move_indexed_particle(&send_buf_r, cell, p);
		  if(p < cell->n) p--;
	    }
	    /* Sort particles in cells of this node during last direction */
	    else if(dir==2) {
	      sort_cell = dd_save_position_to_cell(part->r.p);
	      if(sort_cell != cell) {
		if(sort_cell==NULL) {
		  CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles: Take another loop",this_node));
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: CP1 Particle %d (%f,%f,%f) not inside node domain.\n",
				     this_node,part->p.identity,part->r.p[0],part->r.p[1],part->r.p[2]));		 
		  finished=0;
		  sort_cell = local_cells.cell[0];
		  if(sort_cell != cell) {
		    move_indexed_particle(sort_cell, cell, p);
		    if(p < cell->n) p--;
		  }
		}
		else {
		  move_indexed_particle(sort_cell, cell, p);
		  if(p < cell->n) p--;
		}
	      }
	    }
	  }
    }

	/* Exchange particles */
	if(node_pos[dir]%2==0) {
	  send_particles(&send_buf_l, node_neighbors[2*dir]);
	  recv_particles(&recv_buf_r, node_neighbors[2*dir+1]);
	  send_particles(&send_buf_r, node_neighbors[2*dir+1]);
	  recv_particles(&recv_buf_l, node_neighbors[2*dir]);
	}
	else {
	  recv_particles(&recv_buf_r, node_neighbors[2*dir+1]);
	  send_particles(&send_buf_l, node_neighbors[2*dir]);
	  recv_particles(&recv_buf_l, node_neighbors[2*dir]);
	  send_particles(&send_buf_r, node_neighbors[2*dir+1]);
	}
	/* sort received particles to cells, folding of coordinates also happens in here. */
	if(dd_append_particles(&recv_buf_l, 2*dir  ) && dir == 2) finished = 0;
	if(dd_append_particles(&recv_buf_r, 2*dir+1) && dir == 2) finished = 0; 
	/* reset send/recv buffers */
	send_buf_l.n = 0;
	send_buf_r.n = 0;
	recv_buf_l.n = 0;
	recv_buf_r.n = 0;
      }
      else {
	/* Single node direction case (no communication) */
	/* Fold particles that have left the box */
	/* particle loop */
	for(c=0; c<local_cells.n; c++) {
	  cell = local_cells.cell[c];
	  for (p = 0; p < cell->n; p++) {
	    part = &cell->part[p];
		fold_coordinate(part->r.p, part->l.i, dir);
	    if (dir==2) {
	      sort_cell = dd_save_position_to_cell(part->r.p);
	      if(sort_cell != cell) {
		if(sort_cell==NULL) {
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: CP2 Particle %d (%f,%f,%f) not inside node domain.\n",
				     this_node,part->p.identity,part->r.p[0],part->r.p[1],part->r.p[2]));
		  finished=0;
		  sort_cell = local_cells.cell[0];
		  if(sort_cell != cell) {
		    move_indexed_particle(sort_cell, cell, p);
		    if(p < cell->n) p--;
		  }      
		}
		else {
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: move particle id %d\n", this_node,part->p.identity));
		  move_indexed_particle(sort_cell, cell, p);
		  if(p < cell->n) p--;
		}
	      }
	    }
	  }
	}
      }
    }

    /* Communicate whether particle exchange is finished */
    if(global_flag == CELL_GLOBAL_EXCHANGE) {
      if(this_node==0) {
        int sum;
   
        MPI_Reduce(&finished, &sum, 1, MPI_INT, MPI_SUM, 0, comm_cart);
        if( sum < n_nodes ) finished=0; else finished=sum;
        
    fprintf(stderr, "domdec: %f , try %i , %i nodes have all particles in correct cells\n", lees_edwards_offset, count, sum);
    
    
    
      } else {
	MPI_Reduce(&finished, NULL, 1, MPI_INT, MPI_SUM, 0, comm_cart);
      }
      MPI_Bcast(&finished, 1, MPI_INT, 0, comm_cart);
      
    
      
    } else {
      if(finished == 0) {
	char *errtext = runtime_error(128);
	ERROR_SPRINTF(errtext,"{004 some particles moved more than min_local_box_l, reduce the time step} ");
	/* the bad guys are all in cell 0, but probably their interactions are of no importance anyways.
	   However, their positions have to be made valid again. */
	finished = 1;
	/* all out of range coordinates in the left overs cell are moved to (0,0,0) */
	cell = local_cells.cell[0];
	for (p = 0; p < cell->n; p++) {
	  part = &cell->part[p];
	  if(dir < 3 && dir >= 0 && (part->r.p[dir] < my_left[dir] || part->r.p[dir] > my_right[dir]))
	    for (i = 0; i < 3; i++)
	      part->r.p[i] = 0;
	}
      }
    }
    CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles: finished value: %d\n",this_node,finished));
  }

  realloc_particlelist(&send_buf_l, 0);
  realloc_particlelist(&send_buf_r, 0);
  realloc_particlelist(&recv_buf_l, 0);
  realloc_particlelist(&recv_buf_r, 0);

#ifdef ADDITIONAL_CHECKS
  check_particle_consistency();
#endif

  CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles finished\n",this_node));
}

#endif //LEES_EDWARDS
