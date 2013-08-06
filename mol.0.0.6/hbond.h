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

enum HBondProp { UNKNOWN_HPROP = 0, HBOND_DONOR = 1, HBOND_ACCEPTOR = 2, DONATABLE_HYDROGEN = 4, ROTATABLE_HYDROGEN = 8, FIXED_HYDROGEN = 16 };

enum Hybridization_State { UNKNOWN_HYBRID = 0, SP1_HYBRID, SP2_HYBRID, SP3_HYBRID, RING_HYBRID };
enum HB_Donor_Type { hbdon_NO = 0, hbdon_BB, hbdon_SC, hbdon_SM, hbdon_MAX }; // Small-mol added		// Donor Classes
enum HB_Acceptor_Type { hbacc_NO = 0, hbacc_BB, hbacc_SP2,  hbacc_SP3, hbacc_RING, hbacc_SP2_SM,  hbacc_SP3_SM, hbacc_RING_SM ,hbacc_MAX }; // Acceptor Classes
enum HB_Type_Evaluation { hbe_NONE = 0, hbe_BB, hbe_BBTURN, hbe_BBHELIX, hbe_BBOTHER, hbe_SP2B, hbe_SP3B, hbe_RINGB, hbe_BSC, hbe_SP2SC, hbe_SP3SC, hbe_RINGSC, /* Small Molecule Section -> */ hbe_BSM, hbe_SP2SCSM, hbe_SP3SCSM, hbe_RINGSCSM, hbe_SP2SMSC, hbe_SP3SMSC, hbe_RINGSMSC, hbe_SP2SMB, hbe_SP3SMB, hbe_RINGSMB,hbe_SP2SM, hbe_SP3SM, hbe_RINGSM, hbe_MAX };
enum HB_Weight_Type { hbw_NONE = -1,  hbw_TOTAL, hbw_SR_BB, hbw_LR_BB, hbw_BB_SC, hbw_SC, hbw_BB_SM, hbw_SC_SM, hbw_SM, hbw_H2O};

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

double get_pairwise_hbondeng_nblist( mol_atom *atoms_hydro, int hydro_id, mol_atom *atoms_acc, int acc_id, double *engcat, double rc2, int comp_grad );


/* --------------------------------------------------------------------------

Hydrogen Bonding Energy Functions in presence of Small Molecules

Receptor - Small Molecule Complex

By: Mohammad Moghadasi          mohamad@bu.edu

March 2013

-------------------------------------------------------------------------- */

struct tnode{
        int rig_labl;       // the label for this rigid cluster
        int atm_size;      //number of atoms in this rigid cluster
        int *atm_num;      // atoms which are in this rigid cluster
};

struct rig_tree {
        int v;          // number of nodes in this tree, as it is a tree we know that it has v-1 edges
        struct tnode* vertex;     // pointers to the nodes of the tree
        int *deg;       // The degree for each node
        int **Adj;      // The Adjecent list for the tree
        int *par;       // for each node keep its parent
        int *e;         // list of edges in the tree
        int *atom_e;    // for each edge, index of the atom for the rotatable bond
        int total_atm_size; // total number of atoms in this tree
        int *atm_rig_label; // for each atom we keep it belongs to which rigid cluster
        int *bfs_qu;       // It is the queue of BFS
        double *theta;     // It keeps for each rotatable bond, what is the corresponding rotation
        int *e_index;      // it keeps for each edge according to the bfs queue which edge is responsible from the list of edges and save the index of that edge
        //double *orig;      // it keeps the original coordinate of the atoms
        double *rot;      // 3 parameters for the exponential coordinate of the whole body
        double *trans;    // 3 parameters for the translation of the whole body
        double *cent_t;    // the center of rotation for whole body
        double *der;
        int *tree_atoms;   // keep the list of the atoms which are in this tree and the size of it is total_atm_size
};

struct rigid_body {
        int rig_size;           // number of atoms in this rigid block
        int *rigid_atoms;       // list of the atom numbers which are in this block
        //double *rigid_orig;
        // original coordinate of the atoms which are in this rigid block

double *rigid_der;      // derivative of the energy function ( 6 parameters ) with respect to translation and rotation
        double *trans;          // 3 parameters for the translation of this rigid body
        double *rot;            // 3 parameters for the rotation of this rigid body
        double *cent_t;         // the center of rotation for whole body, by default it is the center of mass for rigid bodies
};

struct rig_forest {
        int numf;
        int numt;
        int numr;
        int numtor;

        int *free_atoms;
        double *free_cor;
        double *free_der;

        double *orig;         // the original coordinate of the atoms

        struct rig_tree** tr;
        struct rigid_body** r_b;
        struct rig_tree** tor;
};

// Functions

void hbondengcat_smallmol( struct atomgrp *ag, double *energy, struct nblist *nblst , struct rig_forest* prot );

void get_categorized_hbondeng_smallmol( double *engcat, double *bb_bb_sr, double *bb_bb_lr, double *bb_sc, double *sc_sc, double *bb_sm, double *sc_sm, double *sm_sm );

void set_categorized_hbondeng_smallmol( double *engcat, double bb_bb_sr, double bb_bb_lr, double bb_sc, double sc_sc, double bb_sm, double sc_sm, double sm_sm );

void init_categorized_hbondeng_smallmol( double *engcat );

double *alloc_categorized_hbondeng_smallmol( void );


#endif
