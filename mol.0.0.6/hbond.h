/*
Copyright (c) 2009-2012, Structural Bioinformatics Laboratory, Boston University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- Neither the name of the author nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific
  prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef _HBOND_H_
#define _HBOND_H_

#include <math.h>
#include <limits.h>

//#define USE_LONG_DOUBLE
#define NON_POSITIVE_POLY

#ifdef USE_LONG_DOUBLE
   #define FLOAT long double
#else
   #define FLOAT double
#endif   

#define LINEAR_FADE  1
#define LOG_FADE     2
#define BSPLINE_FADE 3

#ifndef FADING_FUNCTION 
  #define FADING_FUNCTION BSPLINE_FADE
#endif  

#define HB_ENG_MAX    0.0     // Maximum energy value to consider as hbond

#define MIN_AH  1.4	// Acceptor-Hydrogen distance range 
#define MAX_AH	3.0

#define MIN_xH  0.0	// Psi cutoff range 
#define MAX_xH  1.0

#define MIN_xD  0.0 	// Theta cutoff range
#define MAX_xD  1.0

#define MIN_AH2 ( MIN_AH * MIN_AH )	// Not used
#define MAX_AH2 ( MAX_AH * MAX_AH )	// Not used

#define INTPOL_EDGE_DIST	2.1
#define INTPOL_EDGE_ANGL	0.05
#define INTPOL_RANGE   		0.2	// Not Used

#define DIST_SWITCH             2.1	// Hbond Distance cutoff for sidechains
#define INTPOL_MIN ( DIST_SWITCH - INTPOL_RANGE )	// minimum AH to run interpolation
#define INTPOL_MAX ( DIST_SWITCH + INTPOL_RANGE )	// maximum AH to run interpolation

#define D_WP_DIFFERENCE_TOLERANCE 0.3  /* for water-mediated hydrogen bonding the maximum allowed difference (in angstroms) between dist( water-oxygen, donor ) and dist( water-oxygen, acceptor ) */


struct dvector
{
   FLOAT X; /**< translation in the X dimension */
   FLOAT Y; /**< translation in the Y dimension */
   FLOAT Z; /**< translation in the Z dimension */
};


#define ARY( a, n, i, j ) a[ ( i ) * ( n ) + ( j ) ]

typedef struct
{
    int hydro_id, acc_id;
    double en;
} FLOW_EDGE;


typedef struct
{
    int id, cap;
} FLOW_MAP;


typedef struct
{
    int max_n_fedge, max_n;
    FLOW_EDGE *flow_edge;
    FLOW_MAP *flow_map;
    int *cap;
    FLOAT *cost;
    int *deg;
    int *par;
    int *q;
    int *inq;
    FLOAT *pi;
    FLOAT *d;
    int *adj;
    int *fnet;
} FLOW_STRUCT;


int init_flow_struct( FLOW_STRUCT *flow_struct );
void free_flow_struct( FLOW_STRUCT *flow_struct );
int reinit_flow_struct( FLOW_STRUCT *flow_struct, int max_n_fedge, int max_n );

void mark_hbond_donors (
		 struct atomgrp* ag,
		 struct prm* prm
		);

void mark_hbond_acceptors (
		 struct atomgrp* ag,
		 struct prm* prm
		);

void fix_acceptor_bases( struct atomgrp *ag, struct prm *prm );

void hbondeng( struct atomgrp *ag, double *energy, struct nblist *nblst );
void hbondeng_weighted(struct atomgrp *ag, double *energy, struct nblist *nblst, double weight);

void hbondeng_split(struct atomgrp *ag, double *energy, struct nblist *nblst, int atom_split, double weight);

void hbondeng_split_atom_weighted(struct atomgrp *ag, double *energy, struct nblist *nblst, int atom_split, double weight);

void hbondengcat( struct atomgrp *ag, double *energy, struct nblist *nblst );

void hbondeng_all( struct atomgrp *ag, double *energy, struct nblist *nblst );

void hbondeng_bbexc( struct atomgrp *ag, double *energy, struct nblist *nblst );

void water_mediated_hbondeng( struct atomgrp *ag, double *energy );

void flow_hbondeng( struct atomgrp *ag, double *energy, struct nblist *nblst );

void hbondeng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy );

double *alloc_categorized_hbondeng( void );

void init_categorized_hbondeng( double *engcat );

void print_categorized_hbondeng( double *engcat );

void get_categorized_hbondeng( double *engcat, double *bb_bb_sr, double *bb_bb_lr, double *bb_sc, double *sc_sc );

void set_categorized_hbondeng( double *engcat, double bb_bb_sr, double bb_bb_lr, double bb_sc, double sc_sc );

void check_fading_funcs( void );

double get_pairwise_hbondeng_nblist( mol_atom *atoms_hydro, int hydro_id, mol_atom *atoms_acc, int acc_id, double *engcat, double rc2, int comp_grad, double weight );

void mol_set_hbond_bases(struct atomgrp *ag);

int get_hbe_type(mol_atom * atoms, mol_atom * hydro, mol_atom * acc);

#endif
