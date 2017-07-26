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
#ifndef _MOL_ATOM_H_
#define _MOL_ATOM_H_
#include <stdbool.h>

#include "enums.h"
#include "prms.h"

typedef struct atom mol_atom;


// deprecate: rename atom to mol_atom
struct atom
{
	int ingrp ; /**< atom index in the atomgroup */
	int atom_typen; /**< atom type number */
	int acp_type; /**< atom type number */
	int atom_ftypen;/**< atom type number in the forcefield */
        int octree_ptr; /**< index (ptr) to octree leaf node to which this atoms belongs */           
	char *name;
	char *residue_name;
	char *ftype_name;
	char *element;

	int sa; /**< solvent accessible: 1 => solvent accessible, 0 => !1, -1 => undefined */
	int fixed; /**< =1 if atom is immovable , 0 otherwise*/
	int mask; /**< mask this atom: 1=>mask */
	float attl; /**< attraction level: 0 => repel,
				  1 => standard sa atom,
				  >1 => attract at this level (up to 9)
				  -1 => undefined
				 */

	double X,Y,Z; /**< X,Y,Z coordinates of atom */

	double GX,GY,GZ; /**< forces acting on an atom */

	double eps, rminh, eps03, rminh03; /**< vdw parameters of atom */
    
        double acevolume;

	double chrg; /**< electrostatic charge of atom */

        double B; /**< b-factor of an atom */

        double fft_repvdw_radius; /**< b-factor of an atom */
        
        int res_num;/** global residue number based on libmol numbering */
        int res_seq; /**< residue sequence number */
        int comb_res_seq; /**< single sequence of residue numbers combining all chains */
        
        bool backbone; /**< 1 if this atom is part of the backbone, 0 otherwise */
        bool fft_repvdw_filter; /**< 1 if this atom is part of the backbone, 0 otherwise */
        bool fft_no_attvdw; /**< 1 if this atom is part of the backbone, 0 otherwise */

	// deprecated
	struct atombond** bonds; /**< first level bonds of this atom */
	// end deprecated

	// bond indices of this atom's bonds
	int* bondis;

	struct atomangle** angs; /**< angles this atom is involved in */

	struct atomtorsion** tors; /**< torsions this atom is involved in */

	struct atomimproper** imps; /**< impropers this atom is involved in */
	
	int nbonds;
	int nbondis;
	int nangs;
	int ntors;
	int nimps;
	enum HBondProp hprop; /**< properties related to hydrogen bonding */
	enum Hybridization_State hybridization; /**< hybridization state of the hydrogen bond acceptor atom */
        int base, base2; /**< indices to the base atoms of the hbond acceptor and the donor atom (in base) for hydrogens */                   
        char icode;  /**< insertion code */
	enum mol_yeti yeti_type; /* enum mol_yeti */
	double hbond_weight;
};

// warning: this function does not handle all of mol_atom
// - this function mallocs
void
mol_atom_create (mol_atom* a, int nbondis);

// warning: this function does not handle all of mol_atom
// - this function frees
void
mol_atom_destroy (mol_atom* a);

// warning: this function does not handle all of mol_atom
// - this function neither mallocs nor frees
void
mol_atom_copy (mol_atom* asrc, mol_atom* adst);


/**********************************************************
the functions and types below should be deprecated or
moved to separate files
**********************************************************/



/**
  	angle struct
	a0
	 \
	  \
	   a1---a2

	e = k * (th - th0)^2
*/
struct atomangle
{
	struct atom *a0, *a1, *a2; /**< atoms involved in the angle */
	float th; /**< angle theta */

	float th0; /**< equilibrium angle theta */
	float k; /**< spring constant */
};

/**
  	torsion struct
	a0         a3
	 \         /
	  \       /
	   a1---a2

	e = k * (1 + cos (n*chi - d))
*/
struct atomtorsion
{
	struct atom *a0, *a1, *a2, *a3; /**< atoms involved in the torsion */
	float chi; /**< angle chi: angle between the planes a0,a1,a2 and a1,a2,a3 */

	float k; /**< k constant */
	float d; /**< delta constant */
	int n; /**< number of bonds made by atoms a1 and a2 */
};

/**
  	improper struct
	 a1---a0---a3
	      |
	      |
	      a2

	e = k * (psi - psi0)^2
*/
struct atomimproper
{
	struct atom *a0, *a1, *a2, *a3; /**< atoms involved in the improper */
	float psi; /**< angle psi: angle between the planes a0,a1,a2 and a1,a2,a3 */

	float k; /**< spring constant */
	float psi0; /**< equilibrium angle psi */
};

/**
        spring
*/
struct spring
{
        struct atomgrp *agp;  /**< affected atomgroup */
        int *laspr; /**< list of atoms */
        double fkspr; /**< force constant */
        double X0,Y0,Z0; /**< anchor point */
        int naspr ; /**< number of affected atoms */
};

struct springset
{
        struct spring *springs;  /**< array of springs */
        int nsprings; /**< number of springs */
};

/**
        frees a spring set
*/
void free_springset(struct springset *sprst);

/**
  Copies contents of atom struct src to dest.
  The pointers must be to previously allocated memory.
*/
void copy_atom (struct atom* src, struct atom* dest);

/**
  Copies contents of atom struct src to dest.
  The pointers must be to previously allocated memory.
*/
void fullcopy_atom (struct atom* src, struct atom* dest);

/**
	Prints to stderr the contents of the atom struct
*/
void fprintf_stderr_atom (struct atom* atom, struct prm* prm);

/**
 Free bonds, angles, tors, impropers pointers
*/
void free_atom(struct atom* atm);

#endif
