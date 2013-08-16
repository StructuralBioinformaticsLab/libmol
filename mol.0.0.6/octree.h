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
#ifndef _MOL_DYN_OCTREE_H_
#define _MOL_DYN_OCTREE_H_

/** sqrt( 3 ) / 2 */
#define HALF_SQRT_THREE 0.86602540378443864676372317075294

/** initial number of nodes in the free nodes pool */
#define INIT_NUM_OCTREE_NODES 100

/** number of bits to shift right to extract node id from octree pointer */
#define LOW_BITS  14
/** bit mask to extract index in node from octree pointer */
#define LOW_MASK  0x3FFF

/** initial size of the migration array per node required for batch update */
#define INIT_MIGRATION_ARRAY_SIZE 8

/** combine node id ('a') and index in node ('b') to store them as a single integer */
#define create_octree_ptr( a, b ) ( ( ( a ) << LOW_BITS ) + b )
/** extract node id from octree pointer 'c' */
#define get_node_id( c ) ( ( c ) >> LOW_BITS )
/** extract index in node from octree pointer 'c' */
#define get_index_in_node( c ) ( ( c ) & LOW_MASK )

#ifdef zeroIfLess
   #undef zeroIfLess
#endif
/** return 0 if a < b, 1 otherwise */
#define zeroIfLess( a, b ) ( ( ( a ) < ( b ) ) ? 0 : 1 )

/*
#ifndef RECURSION_DEPTH
   #define RECURSION_DEPTH 50
#endif
*/

//#define ADD_ATTR

/** stores the properties of an octree node */
typedef struct
{
        /** minimum x co-ordinate of the cube */
        double lx;     
        /** minimum y co-ordinate of the cube */        
        double ly;
        /** minimum z co-ordinate of the cube */        
        double lz;
        /** dimension of the cube */
        double dim;            

#ifdef ADD_ATTR
        /** sums of the x, y and z co-ordinates of the atoms under this node */
	double sx, sy, sz;     
	/** sum of charges of atoms under this node */
	double sq;             
#endif

        /** number of items/atoms under this node */
	int n;                 
	/** number of fixed items/atomd under this node */
	int nfixed;            

        /** index of (ptr to) the parent node in the octree (-1 for the parent of the root) */
	int p_ptr;             
	/** indices of (ptrs to) child nodes */
	int c_ptr[ 8 ];        

        /** 1 if this node is a leaf (i.e., no further subdivision), 0 otherwise */
	int leaf;              
	/** capacity of the indices array */
	int id_cap;            
	/** number of items in the indices array (used for temporary storage in non-leafs) */
	int id_num;            
	/** indices to the items/atoms in this node (fixed items are listed first) */
	int *indices;          

} OCTREE_NODE;


/** stores the properties of an octree */
typedef struct
{
        /** nodes of the octree with nodes[ 0 ] being the root */
	OCTREE_NODE *nodes;    

        /** maximum number of nodes in a leaf */
	int max_leaf_size;  
	/** maximum dimension of a leaf node */
 	double max_leaf_dim;   

        /** pointer to the list of atoms */
 	mol_atom *atoms;       
 	/** number of atoms in the list */
 	int natoms;            

        /** number of nodes in the nodes array */
	int num_nodes;         
	/** index of the first free node in the nodes array, -1 if array if full */
	int free_node_ptr;     
} OCTREE;


/** stores parameters for evaluating various energy functions using octrees */
typedef struct OPAR
{
        /** pointer to structure containing properties and setup parameters for the given group of atoms */
        struct agsetup *ags;
        /** parameter for scaling charges during electrostatics computation */
        double eps;
        
        /** pointer to the octree containing the static molecule */
        OCTREE *octree_static;
        /** index of the current node in the static octree */
        int node_static;
        /** index of the current node in the moving octree */        
        int node_moving;

        /** pointer to the octree containing the moving molecule */        
        OCTREE *octree_moving;
        
        /** distance cutoff for pairwise interaction evaluation */        
        double dist_cutoff;
        /** distance cutoff for hydrogen bonding pairwise interaction evaluation */        
        double hdist_cutoff;
        /** distance cutoff for approximate pairwise interaction evaluation */                
        double approx_cutoff;
        
        /** 3 x 4 transformation matrix */
        double *trans;
        
        /** array for storing categorized hydrogen bonding energy */
        double *engcat;
        
        /** pointer to data to be passed as parameters to the interaction evaluation function */
        void *proc_func_params;
        /** user-defined pairwise interaction evaluation function */
        void ( * processing_function )( struct OPAR *, double * );
        
        /** if set to 1, ignore computations involving node pairs containing only fixed atoms */        
        int fixed_cull;
} OCTREE_PARAMS;


/*
typedef struct
{
   int node_static, static_cid;
   int node_moving, moving_cid;
} REC_PARAMS;
*/

#ifdef _WIN32
#define _mol_sinline static
#else
#define _mol_sinline static inline
#endif

/**
   Given a 3 x 4 transformation matrix pointed to by trans_mat, this function transforms
   point < x, y, z > to point < *nx, *ny, *nz >.
*/
_mol_sinline void transform_point( double x, double y, double z, double *trans_mat, double *nx, double *ny, double *nz )
{
   *nx = trans_mat[  0 ] * x + trans_mat[  1 ] * y + trans_mat[  2 ] * z + trans_mat[  3 ];
   *ny = trans_mat[  4 ] * x + trans_mat[  5 ] * y + trans_mat[  6 ] * z + trans_mat[  7 ];
   *nz = trans_mat[  8 ] * x + trans_mat[  9 ] * y + trans_mat[ 10 ] * z + trans_mat[ 11 ];
}

/**
   Given a point < x, y, z > and a box of dimension dim x dim x dim with its minimum X, Y and Z 
   coordinates given by bx, by and bz, respectively, this function returns the squared distance
   from < x, y, z > to the closest point on the box.
*/
_mol_sinline double min_pt2bx_dist2( double bx, double by, double bz, double dim, double x, double y, double z )
{
   double dx, dy, dz;
   double dl, dr;
   
   dl = bx - x; 
   dr = - dim - dl;
   dx = ( dl > 0 ) ? dl : ( ( dr > 0 ) ? dr : 0 );   

   dl = by - y; 
   dr = - dim - dl;
   dy = ( dl > 0 ) ? dl : ( ( dr > 0 ) ? dr : 0 );      
   
   dl = bz - z; 
   dr = - dim - dl;
   dz = ( dl > 0 ) ? dl : ( ( dr > 0 ) ? dr : 0 );   
      
   return ( dx * dx + dy * dy + dz * dz );
}

/**
   Returns 0 if the two boxes given by < x_1, y_1, z_1, dim_1 > (where x_1, y_1 and z_1 are the smallest
   coordinates of the box, and dim1 is its dimension) and < x_2, y_2, z_2, dim_2 > apart by a distance
   of at least 'ext'. If not, returns 1 provided the separation between the smallest spheres bounding 
   the two boxes is less than 'ext'.    
*/
_mol_sinline int within_distance_cutoff( double x_1, double y_1, double z_1, double dim_1, double x_2, double y_2, double z_2, double dim_2, double ext )
{
   double d_1 = dim_1 + ext, d_2 = - ( dim_2 + ext );

   if ( ( x_2 - x_1 >= d_1 ) || ( x_2 - x_1 <= d_2 )
     || ( y_2 - y_1 >= d_1 ) || ( y_2 - y_1 <= d_2 )
     || ( z_2 - z_1 >= d_1 ) || ( z_2 - z_1 <= d_2 ) ) return 0;   
   
   d_1 = 0.5 * ( dim_2 - dim_1 );
  
   d_2 = ( x_2 - x_1 + d_1 ) * ( x_2 - x_1 + d_1 ) 
      + ( y_2 - y_1 + d_1 ) * ( y_2 - y_1 + d_1 )
      + ( z_2 - z_1 + d_1 ) * ( z_2 - z_1 + d_1 );            

   d_1 = ( HALF_SQRT_THREE * ( dim_1 + dim_2 ) + ext );
   d_1 *= d_1;
   
   return ( d_2 < d_1 );   
}

/**
   Returns nonzero if the center of *atom is inside the cubic volume spanned
   by *node (an octree node), otherwise return zero.  
*/
_mol_sinline int inside_node( OCTREE_NODE *node, mol_atom *atom )
{
   return ( ( atom->X - node->lx >= 0 ) && ( atom->X - node->lx < node->dim )
         && ( atom->Y - node->ly >= 0 ) && ( atom->Y - node->ly < node->dim ) 
         && ( atom->Z - node->lz >= 0 ) && ( atom->Z - node->lz < node->dim ) );
}


/**
   Given a mol_atom_group pointed to by 'ag', build an octree with each leaf containing at most 'max_leaf_size' atoms
   and spanning a cubic volume of at most 'max_leaf_dim' x 'max_leaf_dim' x 'max_leaf_dim' (cubic angstroms). 
   The pointer to the octree is assigned to the parameter 'octree'.
*/
int build_octree( OCTREE *octree, int max_leaf_size, double max_leaf_dim, double slack_factor, mol_atom_group *ag );

/**
   Given a mol_atom_group pointed to by 'ag', build an octree from the fixed atoms of *ag, with each leaf containing 
   at 'most max_leaf_size' atoms and spanning a cubic volume of at most 'max_leaf_dim' x 'max_leaf_dim' x 'max_leaf_dim' 
   (cubic angstroms). The pointer to the octree is assigned to the parameter 'octree'.
*/
int build_octree_excluding_fixed_atoms( OCTREE *octree, int max_leaf_size, double max_leaf_dim, double slack_factor, mol_atom_group *ag );

/**
   Print the tree pointed to by 'octree'. 
*/
void print_octree( OCTREE *octree );

/**
   Get the space used (in bytes) by the OCTREE data structure pointed to by 'octree'. 
*/
int get_octree_size( OCTREE *octree );

/**
   Update the OCTREE data structure pointed to by 'octree' with the (possibly) new position
   of the atom pointed to by 'atom'.
*/
void update_octree( OCTREE *octree, mol_atom* atom );

/**
   Batch update the OCTREE data structure pointed to by 'octree' with the (possibly) new positions
   of the atoms it is constructed from. 
*/
int reorganize_octree( OCTREE *octree, int batch_update );

/**
   This function returns the sum of the pairwise interaction values between the atoms stored in 'octree_static' 
   and those in 'octree_moving'. All atom pairs with pairwise distance larger than 'dist_cutoff'
   are ignored. If set to a value smaller than 'dist_cutoff', the parameter 'approx_cutoff' 
   speeds up computation by allowing some false positives during the determination of whether two octree 
   nodes are within 'dist_cutoff'. If 'fixed_cull' is set to nonzero, ignores computations involving node 
   pairs containing only fixed atoms. The parameter 'trans' points to a 3 x 4 transformation matrix
   for transforming the molecule stored in 'octree_moving' (no transformation if 'trans' is
   set to NULL). The parameter 'proc_func_params' points to data to be passed to the user-defined
   pairwise interaction computation function pointed to by 'processing_function'.
*/
double octree_accumulation_excluding_far( OCTREE *octree_static, OCTREE *octree_moving,
                                          double dist_cutoff, double approx_cutoff, int fixed_cull, double *trans,
                                          void *proc_func_params,
                                          void ( * processing_function )( OCTREE_PARAMS *, double * ) );
                                          
/**
   Free the memory allocated to the OCTREE data structure pointed to by 'octree'.
*/
void destroy_octree( OCTREE *octree );

#endif
