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
/** \file lees_edwards_domain_decomposition.h
    Data and methods for Lees-Edwards periodic boundary conditions. Most of the domain decomposition is handled
    as normal, (see \ref domain_decomposition.h ) however imaging in the Y-axis is neccessarily a little difference.
*/
#ifndef _LEES_EDWARDS_DD_HPP
#define _LEES_EDWARDS_DD_HPP

#include "ghosts.hpp"

#ifdef LEES_EDWARDS
/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Create communicators for cell structure domain decomposition. (see \ref GhostCommunicator and also \ref dd_prepare_comm) */
void  le_dd_prepare_comm(GhostCommunicator *comm, int data_parts);

/** update the 'shift' member of those GhostCommunicators, which use
    that value to speed up the folding process of its ghost members
    (see \ref dd_prepare_comm for the original), i.e. all which have
    GHOSTTRANS_POSSHFTD or'd into 'data_parts' upon execution of \ref
    dd_prepare_comm. */
void le_dd_update_communicators_w_boxl();
/** Init cell interactions for the Lees-Edwards cell system.
 * initializes the interacting neighbor cell list of a cell. The
 * created list of interacting neighbor cells is used by the verlet
 * algorithm (see verlet.c) to build the verlet lists.
 */
void le_dd_init_cell_interactions();
#endif //LEES_EDWARDS
#endif //LEES_EDWARDS_DD_H
