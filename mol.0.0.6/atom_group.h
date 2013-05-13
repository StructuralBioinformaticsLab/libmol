/*
Copyright (c) 2009-2012, Structural Bioinformatics Laboratory, Boston University
Copyright (c) 2013, Acpharis Inc
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
#ifndef _MOL_ATOM_GROUP_H_
#define _MOL_ATOM_GROUP_H_
#include <stdbool.h>

#include "atom.h"
#include "bond.h"
#include "matrix.h"
#include "vector.h"

typedef struct atomgrp mol_atom_group;

enum mol_res_type {
	UNK,
	ALA,
	ARG,
	ARGN,
	ASN,
	ASP,
	ASPH,
	CYS,
	GLN,
	GLU,
	GLUH,
	GLY,
	HIS,
	ILE,
	LEU,
	LYS,
	LYSN,
	MET,
	PHE,
	PRO,
	SER,
	SEP,
	THR,
	TRP,
	TYR,
	VAL,
	HSC,
	HSD,
	ACE,
	PCA,
	HYL,
	HYP,
	HSE,
	ORN,
	PEN,
	ALB,
	ABU,
	ST2,
	TIP3,
	OH2,
	HOH,
	WAT,
};

struct atom_type_name_num {

  int num;
  char name[5];

};

/**
	Holds list of atoms and
	the number of atoms in
	the list.
*/
struct atomgrp {

	int natoms; /**< number of atoms in the group */
        int nactives;/**< number of active atoms in the group */
	mol_atom* atoms; /**< pointer to the array of atoms in the group */

        struct atom_type_name_num* atypenn;

        int *activelist;/**< list of active atoms in the group */

	int nbonds;
        int nbact;
	mol_bond* bonds; /**< all first level bonds of this atomgrp */
        
        struct atombond** bact; /**< array of pointers to the bonds of active atom/s */

	int nangs;
        int nangact;
	struct atomangle* angs; /**< all angs of this atomgrp */

        struct atomangle** angact; /**< all active angs of this atomgrp */

	int ntors;
        int ntoract;
	struct atomtorsion* tors; /**< all tors of this atomgrp */

        struct atomtorsion** toract; /**< all active tors of this atomgrp */

	int nimps;
        int nimpact;
	struct atomimproper* imps; /**< all imps of this atomgrp */
        
        struct atomimproper** impact; /**< all active imps of this atomgrp */

        int num_atom_types;
        int nres;
        int *iares;
        int *rot;
        char **idres;
        enum mol_res_type *res_type;

        struct prm *prm;        
        
        void *flow_struct; // for netfork-flow based hydrogen bonding
	char *atom_group_name;
        bool is_psf_read; //psf has been read in
};


// warning: this function does not handle all of mol_atom_group
// - this function mallocs
void
mol_atom_group_create_from_template (
		mol_atom_group* ag, mol_atom_group* agtmplt);


// warning: this function does not handle all of mol_atom_group
// - this function frees
void
mol_atom_group_destroy (mol_atom_group* ag);


// warning: this function does not handle all of mol_atom_group
// - this function neither mallocs nor frees
void
mol_atom_group_copy (mol_atom_group* agsrc, mol_atom_group* agdst);


// input req: minmaxs should be sizeof 6 floats
// output: minmaxs filled with min and max coords
//         of ag for every dimension
void
mol_atom_group_minmaxs (mol_atom_group* ag, float *minmaxs);


/**
	Frees all atoms in the atomgrp and the atomgrp itself.
*/
void full_free_atomgrp (struct atomgrp* ag);

/**
        Frees atomgroup light way 
*/
void free_atomgrp (struct atomgrp* ag);

/**
	Creates a copy of srcag and returns it.
*/
struct atomgrp* copy_atomgrp (struct atomgrp* srcag);

/**
        Creates a copy of srcag and returns it.
*/
struct atomgrp* fullcopy_atomgrp (struct atomgrp* srcag);

/**
  extract all atoms of type and return them in atomgrp
*/
struct atomgrp* extract_type (struct atomgrp* ag, const char* type, struct prm* prm);

/**
  remove all atoms of type and return the remaining atoms in atomgrp
*/
struct atomgrp* rm_type (struct atomgrp* ag, const char* type, struct prm* prm);

struct atomgrp* exrm_type (struct atomgrp* ag, const char* type, struct prm* prm, int direction);

/**
  typemaj extraction and removal
*/
struct atomgrp* extract_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm);
struct atomgrp* rm_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm);
struct atomgrp* exrm_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm, int direction);

/**
  make all ag atoms sa
*/
void full_sa (struct atomgrp* ag);

/**
  join ags into a single atom group ag and return ag
*/
struct atomgrp* join_atomgrps (struct atomgrp** ags);
/**
  join ag1 & ag2 into a single atom group ag and return ag
*/
struct atomgrp* join_2atomgrps (struct atomgrp* ag1, struct atomgrp* ag2);

/**
	Prints the contents of the atomgrp struct to stdout
*/
void print_atomgrp (struct atomgrp* ag, struct prm* prm);

/**
  check the properties of atomgrp for consistency with prm
*/
void check_atomgrp (struct atomgrp* ag, struct prm* prm);

/**
  fill the atomgroup index in atoms
*/
void fill_ingrp (struct atomgrp* ag);

/**
  Copy coordinates from atom group to an array
*/

void ag2array( double* array, struct atomgrp* ag);

/**
  Copy coordinates from an array to the atom group
*/

void array2ag ( double* array, struct atomgrp* ag);

/**
  Replace coordinates in an atom group with those from a pdb file
*/

void replace_coordinates(struct atomgrp* ag,const char* pdb_path);


void transform_atomgrpf(struct atomgrp* ag, struct mol_matrix3f rotation, struct mol_vector3f translation);
struct atomgrp* join_rec_lig_ff(struct atomgrp* rec, struct atomgrp* lig);

#endif
