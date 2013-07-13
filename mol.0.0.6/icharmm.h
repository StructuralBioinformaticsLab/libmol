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
#ifndef _MOL_ICHARMM_H_
#define _MOL_ICHARMM_H_

/**
	routines for reading charmm forcefield parameters for the
	particular atom structure 
*/

/**
	return an array of charmm atom type names for the
	particular atom structure
*/
char* get_atnam_rtf(char* rtffile, int* atoms, int nat);
/**
        read charmm psf file: 
	o_nat-number of atoms, 		o_atoms-atom type numbers,
	o_nb-number of bonds, 		o_bonds-2atom numbers for each bond,
	o_nang-number of angles, 	o_ang-3atom numbers for each angle,
	o_ndih-number of dihedrals, 	o_dihs-4atom numbers for each dihedral,
	o_nimp-number of impropers, 	o_imps-4atom numbers for each improper,
	o_cha-charges for atoms
*/
void read_psf(const char* psffile, int* o_nat,  int** o_atoms, char*** o_atom_names, int* o_nb,   int** o_bonds,
                                   int* o_nang, int** o_angs,  int* o_ndih, int** o_dihs,
                                   int* o_nimp, int** o_imps, double** o_cha, 
                                   int* o_nres, int** o_ires, char** o_dres,
                                   int* o_ndon, int** o_dons, int* o_nacc, int** o_accs);
/**
	read forcefield constants from charmm papameter file
	**kbond, **lbond:-bonds
	**o_kangle, **o_fangle:-angles
	**o_kdih, **o_fdih, **o_pdih:-dihedrals (pdih-periodicity)
	**o_kimp,**o_fimp:-impropers
	**o_epa, **o_siga:-vdw parameters
*/
void read_par(char* prmfile, int nat,  char* atnam,
                int nb,   int *bonds, double **kbond, double **lbond,
                int nang, int *angles, double **o_kangle, double **o_fangle,
                int *ndih, int *dihs, double **o_kdih, double **o_fdih, int **o_pdih, int **o_tdih,
                int nimp, int *imps, double **o_kimp, double **o_fimp,
                double **o_epa, double **o_siga, double **o_acevolumes);
/**
        return *pptext, a pointer to the n-th word (starting from 0) in a string ptext
	with a delimiter ch
*/
int nthv(char **pptext, char *ptext, char ch, int n);
/**
	get bond parameters from a string buffer
*/
void chkbond(char* buffer, int nb, int* bonds,
        char* atnam, double* kbond, double* lbond);
/**
	based on existing bond i generate a string tag
	to be matched by a parameter record
*/
void bondtag(int i, int* bonds, char* atnam, char* tag);
/**
        get angle parameters from a string buffer
*/
void chkangle(char* buffer, int nang, int* angles,
        char* atnam, double* kangle, double* fangle);
/**
        based on existing angle i generate a string tag
        to be matched by a parameter record
*/
void angletag(int i, int* angles, char* atnam, char* tag);
/**
        get torsion parameters from a string buffer
*/
void chkdih(char* buffer, int *ndih, int ndih0, int* dihs,
                char* atnam, double* kdih, double* fdih, int* pdih, int *tdih, int *wdih);
/**
        based on existing torsion or improper i generate a string tag
        to be matched by a parameter record (mode is a way of dealing
	with wild cards)
*/
void dihtag(int i, int* dihs, char* atnam, int mode, char* tag);
/**
        get improper parameters from a string buffer
*/
void chkimp(char* buffer, int nimp, int* imps,
                char* atnam, double* kimp, double* fimp);
/** 
	put "X   " replacement into string p starting with character n
*/
void putXtoP(char* p, int n);
/**
        get vdw parameters from a string buffer
*/
void chkvdw(char* buffer, int nat, char* atnam, double* epa, double* siga);
/**
        get ACE volumes parmaeter from string
*/
void chkacevolumes(char* buffer, int nat, char* atnam, double* acevolumes);
/**
        main routine to read all charmm parameters for an atom group ag
*/

void read_ff_charmm(const char* psffile, char* prmfile, char* rtffile, struct atomgrp* ag);
/**
       update fixed atom index and bonded active structures. list contains indices of fixed atoms
*/
void fixed_update(struct atomgrp* ag, int nlist, int* list);
void fixed_update_nolist(struct atomgrp *ag);
/**
	initialize active structures
*/
void fixed_init(struct atomgrp* ag);
/**
         routine to read bond parameters for an atom group ag from psf
*/
void read_ff_charmm_light(const char* psffile, struct atomgrp* ag);

void read_atom_type_name_num(char *rtffile, int *num_atom_types, struct atom_type_name_num* atypenn);

#endif
