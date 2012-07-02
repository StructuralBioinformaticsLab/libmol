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
#ifndef _MOL_PRMS_H_
#define _MOL_PRMS_H_

/**
	\file prm.h
	This file contains structures and functions
	that are used to initialize prm from
	the atom parameter file.
*/

struct prm
{
	// atom
	int natoms; /**< number of atoms */
	struct prmatom* atoms;
	int nsubatoms; /**< number of subatom types */

	//float* rs;
	//float* chrgs;

	// pairwise potential
	struct prmpwpot* pwpot;

	// bond
	int nbonds;
	struct prmbond* bonds;

};

struct prmatom
{
	int id; // atom id
	char* typemaj; // e.g. GLY
	char* typemin; // e.g. CA
	int subid; // subatom id

	float r; // vdw radius
	float q; // partial charge
};

/**
  pairwise potential
*/
struct prmpwpot
{
	float r1,r2; /**< r1,r2 for potential. */
	int k; /**< number of eigenvalues */
	float* lambdas; /**< pointer to array of eigenvalues */
	float* Xs; /**< pointer to matrix of eigenvectors */
        float **eng;/** matrix of eps_ij */
};

struct prmbond
{
	int i,j; // subatom i,j
	float k; // spring constant
	float l0; // equilibrium length
};

/**
	Current line type being read.
*/
enum ereadstate
{
	VERSION,
	ATOM,
	HYDROGEN,
	RADIUS,
	POTENTIAL,
	BOND,
};

enum ereaderr
{
	ERR_VERSION,
	ERR_ATOM,
	ERR_NAMELEN,
	ERR_HYDROGEN,
	ERR_RADIUS,
	ERR_POTENTIAL,
};

/**
	Defines atoms types, pairwise energies,
	and charges.
*/
struct prm* read_prm (const char* path, const char* bin_version);
void read_prmversion (const char* path, const char* bin_version);
void read_prmatom (struct prm* prm, const char* path);
void read_prmpwpot (struct prm* prm, const char* path);
void read_prmbond (struct prm* prm, const char* path);

/**
	Read the PDB file, and populate only the typemin and typemaj fields 
	of prm->atoms with only a single entry for each unique <typemaj, typemin> 
	tuple in the PDB. This is mainly for compatibility with other parts
	of the code already written. However, for large PDB's this will also 
	save some space.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]
*/
void read_typeinfo_from_pdb (struct prm* prm, const char* path);

/**
  Return the atom id of the atypemaj, atypemin key
*/
int atomid (struct prm* prm, const char* atypemaj, const char* atypemin);

/**
	Returns a malloced copy of iprm,
	currently only copies prmatom

	\param iprms parameter struct to be copied
	\return pointer to copy of iprms
*/
struct prm* copy_prm (struct prm* iprm);
void copy_prmatom (struct prmatom* iatom, struct prmatom* oatom);

/**
	Multiply radii in prm by k.
	\param prm parameter struct to be modified
	\param k value with which to multiply radii
	\return void
*/
void modify_prms_radii (struct prm* prm, float k);

/**
	Frees the contents of prm and prm itself.
*/
/*
void free_prms (struct prm* prm);
*/

/**
	Prints an error message based on readerr.
*/
void print_readerr (enum ereaderr readerr, const char* path, char* line);

/**
	Prints the prm struct to stdout
*/
void print_prm (struct prm* prm);
void print_prmatom (struct prm* prm);
void print_prmpwpot (struct prm* prm);
void print_prmbond (struct prm* prm);

void free_prm (struct prm* prm);

#endif
