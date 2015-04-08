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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>

#include _MOL_INCLUDE_

// warning: this function does not handle all of mol_atom_group
void
mol_atom_group_create_from_template (
		mol_atom_group* ag, mol_atom_group* agtmplt)
{
	int i;

	//nulls out any unset pointers, zeroes out counts
	memset(ag, 0, sizeof(mol_atom_group));

	assert (ag != NULL);
	assert (agtmplt != NULL);

	// create atoms
	ag->natoms = agtmplt->natoms;
	if (ag->natoms > 0)
	{
		ag->atoms = _mol_malloc (ag->natoms * sizeof (mol_atom));
		for (i=0; i < ag->natoms; i++)
			mol_atom_create (&ag->atoms[i], agtmplt->atoms[i].nbondis);
	}
	else
		ag->atoms = NULL;

	// create bonds
	ag->nbonds = agtmplt->nbonds;
	if (ag->nbonds > 0)
	{
		ag->bonds = _mol_malloc (ag->nbonds * sizeof (mol_bond));
	}
	else
		ag->bonds = NULL;

    ag->nres = 0;
    ag->nangs = 0;
    ag->ntors = 0;
    ag->nimps = 0;

	return;
}

// warning: this function does not handle all of mol_atom_group
void
mol_atom_group_destroy (mol_atom_group* ag)
{
	int i;

	assert (ag != NULL);

	// destroy atoms
	if (ag->natoms > 0)
	{
		for (i=0; i < ag->natoms; i++)
			mol_atom_destroy (&ag->atoms[i]);
		free (ag->atoms);
	}

	// destroy bonds
	if (ag->nbonds > 0)
	{
		free (ag->bonds);
	}

    if (ag->nres > 0)
	{
		for (i=0; i < ag->nres; i++)
			free (ag->idres[i]);
		free (ag->idres);
		free (ag->iares);
	}

    if (ag->nangs > 0)
	{
		free (ag->angs);
    }

    if (ag->ntors > 0)
	{
		free (ag->tors);
    }

    if (ag->nimps > 0)
	{
		free (ag->imps);
    }

    if (ag->nactives > 0)
	{
		free (ag->activelist);
    }

    if (ag->nbact > 0)
	{
		free (ag->bact);
    }

    if (ag->nangact > 0)
	{
		free (ag->angact);
    }

    if (ag->ntoract > 0)
	{
		free (ag->toract);
    }

    if (ag->nimpact > 0)
	{
		free (ag->impact);
    }

	return;
}

// warning: this function does not handle all of mol_atom_group
void
mol_atom_group_copy (mol_atom_group* agsrc, mol_atom_group* agdst)
{
	int i;

	// check integrity
	assert (agsrc != NULL);
	assert (agdst != NULL);

	assert (agsrc->natoms == agdst->natoms);
	if (agdst->natoms > 0)
		assert (agdst->atoms != NULL);

	assert (agsrc->nbonds == agdst->nbonds);
	if (agdst->nbonds > 0)
		assert (agdst->bonds != NULL);

	// copy atoms
	for (i = 0; i < agdst->natoms; i++)
		mol_atom_copy (&agsrc->atoms[i], &agdst->atoms[i]);

	// copy bonds
	for (i = 0; i < agdst->nbonds; i++)
		mol_bond_copy (&agsrc->bonds[i], &agdst->bonds[i]);

	return;
}

void
mol_atom_group_minmaxs (mol_atom_group* ag, float *minmaxs)
{
	float minx, maxx, miny, maxy, minz, maxz;
	int i;

	assert (minmaxs != NULL);
	assert (ag->natoms > 0);

	minx = ag->atoms[0].X;
	maxx = ag->atoms[0].X;
	miny = ag->atoms[0].Y;
	maxy = ag->atoms[0].Y;
	minz = ag->atoms[0].Z;
	maxz = ag->atoms[0].Z;

	for (i = 0; i < ag->natoms; i++)
	{
		if (ag->atoms[i].X < minx)
			minx = ag->atoms[i].X;
		if (ag->atoms[i].X > maxx)
			maxx = ag->atoms[i].X;

		if (ag->atoms[i].Y < miny)
			miny = ag->atoms[i].Y;
		if (ag->atoms[i].Y > maxy)
			maxy = ag->atoms[i].Y;

		if (ag->atoms[i].Z < minz)
			minz = ag->atoms[i].Z;
		if (ag->atoms[i].Z > maxz)
			maxz = ag->atoms[i].Z;
	}

	// minmaxs packed like this:
	minmaxs[0] = minx;
	minmaxs[1] = maxx;
	minmaxs[2] = miny;
	minmaxs[3] = maxy;
	minmaxs[4] = minz;
	minmaxs[5] = maxz;
}


void free_atomgrp (struct atomgrp* ag)
{
        free (ag->atoms); // free all the atoms in the ag
        free (ag); // free the ag itself
}

void full_free_atomgrp (struct atomgrp* ag)
{
        int i;
        for(i=0; i<ag->natoms; i++)
              free_atom(ag->atoms+i); // free bonded pointers of each atom
	free (ag->atoms); // free all the atoms in the ag
        free (ag->activelist);// free list of active atoms
	if (ag->nbact > 0)
	        free (ag->bact); // free active bond pointers
        free (ag->bonds); // free bonds
	if (ag->nangact > 0)
	        free (ag->angact); // free active angle pointers
        free (ag->angs); // free angles
	if (ag->ntoract > 0)
	        free (ag->toract); // free active torsion pointers
        free (ag->tors); // free torsions
	if (ag->nimpact > 0)
	        free (ag->impact); // free active improper pointers
        free (ag->imps); // free impropers
        for(i=0; i<ag->nres; i++)
              free (ag->idres[i]);    //free residue names
        free (ag->idres);  //free pointers to residue names
        free (ag->iares);  //free array of first atom in residues
        free (ag->rot);
        free (ag->res_type);
        free (ag->atypenn);
        free (ag->atom_group_name);

	free (ag); // free the ag itself
}

struct atomgrp* fullcopy_atomgrp (struct atomgrp* srcag)
{
	int i;
	int j;
	struct atomgrp* destag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));

	destag->natoms = srcag->natoms;
	destag->num_atom_types = srcag->num_atom_types;
	destag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * destag->natoms);
	memcpy(destag->atoms, srcag->atoms, sizeof(struct atom) * destag->natoms);


        //destag->nactives = srcag->nactives ;
        //destag->activelist = _mol_malloc(destag->nactives*sizeof(int));
        //for (i = 0; i < destag->nactives; i++)
        //{
        //     destag->activelist[i]=srcag->activelist[i];
        //}


        destag->nbonds=srcag->nbonds;
	if (srcag->nbonds > 0) {
		destag->bonds= _mol_malloc (destag->nbonds*sizeof(struct atombond));
		memcpy(destag->bonds, srcag->bonds, sizeof(struct atombond) * destag->nbonds);
	}


 	//destag->nbact=srcag->nbact;
        //destag->bact= _mol_malloc (destag->nbact*sizeof(struct atombond*));
        //for (i = 0; i < destag->nbact; i++)
        //{
        //     destag->bact[i]=srcag->bact[i];
        //}

	destag->nangs=srcag->nangs;
	if (srcag->nangs > 0) {
		destag->angs= _mol_malloc (destag->nangs*sizeof(struct atomangle));
		memcpy(destag->angs, srcag->angs, sizeof(struct atomangle) * destag->nangs);
	}

	//destag->nangact=srcag->nangact;
        //destag->angact= _mol_malloc (destag->nangact*sizeof(struct atomangle*));
        //for (i = 0; i < destag->nangact; i++)
        //{
        //     destag->angact[i]=srcag->angact[i];
        //}

	destag->ntors=srcag->ntors;
	if (srcag->ntors > 0) {
		destag->tors= _mol_malloc (destag->ntors*sizeof(struct atomtorsion));
		memcpy(destag->tors, srcag->tors, sizeof(struct atomtorsion) * destag->ntors);
	}

	//destag->ntoract=srcag->ntoract;
        //destag->toract= _mol_malloc (destag->ntoract*sizeof(struct atomtorsion*));
        //for (i = 0; i < destag->ntoract; i++)
        //{
        //     destag->toract[i]=srcag->toract[i];
        //}

        destag->nimps=srcag->nimps;
	if (srcag->nimps > 0) {
		destag->imps= _mol_malloc (destag->nimps*sizeof(struct atomimproper));
		memcpy(destag->imps, srcag->imps, sizeof(struct atomimproper) * destag->nimps);
	}

	//destag->nimpact=srcag->nimpact;
        //destag->impact= _mol_malloc (destag->nimpact*sizeof(struct atomimproper*));
        //for (i = 0; i < destag->nimpact; i++)
        //{
        //     destag->impact[i]=srcag->impact[i];
        //}

	//point to correct atoms
	for (i = 0; i < destag->nbonds; i++) {
		ptrdiff_t a0_offset = destag->bonds[i].a0 - srcag->atoms;
		ptrdiff_t a1_offset = destag->bonds[i].a1 - srcag->atoms;
		destag->bonds[i].a0 = destag->atoms + a0_offset;
		destag->bonds[i].a1 = destag->atoms + a1_offset;
	}
	for (i = 0; i < destag->nangs; i++) {
		ptrdiff_t a0_offset = destag->angs[i].a0 - srcag->atoms;
		ptrdiff_t a1_offset = destag->angs[i].a1 - srcag->atoms;
		ptrdiff_t a2_offset = destag->angs[i].a2 - srcag->atoms;
		destag->angs[i].a0 = destag->atoms + a0_offset;
		destag->angs[i].a1 = destag->atoms + a1_offset;
		destag->angs[i].a2 = destag->atoms + a2_offset;
	}
	for (i = 0; i < destag->ntors; i++) {
		ptrdiff_t a0_offset = destag->tors[i].a0 - srcag->atoms;
		ptrdiff_t a1_offset = destag->tors[i].a1 - srcag->atoms;
		ptrdiff_t a2_offset = destag->tors[i].a2 - srcag->atoms;
		ptrdiff_t a3_offset = destag->tors[i].a3 - srcag->atoms;
		destag->tors[i].a0 = destag->atoms + a0_offset;
		destag->tors[i].a1 = destag->atoms + a1_offset;
		destag->tors[i].a2 = destag->atoms + a2_offset;
		destag->tors[i].a3 = destag->atoms + a3_offset;
	}
	for (i = 0; i < destag->nimps; i++) {
		ptrdiff_t a0_offset = destag->imps[i].a0 - srcag->atoms;
		ptrdiff_t a1_offset = destag->imps[i].a1 - srcag->atoms;
		ptrdiff_t a2_offset = destag->imps[i].a2 - srcag->atoms;
		ptrdiff_t a3_offset = destag->imps[i].a3 - srcag->atoms;
		destag->imps[i].a0 = destag->atoms + a0_offset;
		destag->imps[i].a1 = destag->atoms + a1_offset;
		destag->imps[i].a2 = destag->atoms + a2_offset;
		destag->imps[i].a3 = destag->atoms + a3_offset;
	}

	//point to atoms to correct bonds, angs, tors, imps
	for (i = 0; i < destag->natoms; i++) {
//		struct atom local_atom = destag->atoms[i];
		destag->atoms[i].name = strdup(srcag->atoms[i].name);
		destag->atoms[i].residue_name = strdup(srcag->atoms[i].residue_name);
		destag->atoms[i].ftype_name = strdup(srcag->atoms[i].ftype_name);

		if (srcag->nbonds > 0) {
			destag->atoms[i].bonds = _mol_malloc( destag->atoms[i].nbonds*sizeof(struct atombond*));
			memcpy(destag->atoms[i].bonds, srcag->atoms[i].bonds, sizeof(struct atombond*) * destag->atoms[i].nbonds);
			for (j = 0; j < destag->atoms[i].nbonds; j++) {
				ptrdiff_t bond_offset = destag->atoms[i].bonds[j] - srcag->bonds;
				destag->atoms[i].bonds[j] = destag->bonds + bond_offset;
			}
		}
		if (srcag->nangs > 0) {
			destag->atoms[i].angs = _mol_malloc( destag->atoms[i].nangs*sizeof(struct atomangle*));
			memcpy(destag->atoms[i].angs, srcag->atoms[i].angs, sizeof(struct atomangle*) * destag->atoms[i].nangs);
			for (j = 0; j < destag->atoms[i].nangs; j++) {
				ptrdiff_t ang_offset = destag->atoms[i].angs[j] - srcag->angs;
				destag->atoms[i].angs[j] = destag->angs + ang_offset;
			}
		}
		if (srcag->ntors > 0) {
			destag->atoms[i].tors = _mol_malloc( destag->atoms[i].ntors*sizeof(struct atomtorsion*));
			memcpy(destag->atoms[i].tors, srcag->atoms[i].tors, sizeof(struct atomtorsion*) * destag->atoms[i].ntors);
			for (j = 0; j < destag->atoms[i].ntors; j++) {
				ptrdiff_t tor_offset = destag->atoms[i].tors[j] - srcag->tors;
				destag->atoms[i].tors[j] = destag->tors + tor_offset;
			}
		}
		if (srcag->nimps > 0) {
			destag->atoms[i].imps = _mol_malloc( destag->atoms[i].nimps*sizeof(struct atomimproper*));
			memcpy(destag->atoms[i].imps, srcag->atoms[i].imps, sizeof(struct atomimproper*) * destag->atoms[i].nimps);
			for (j = 0; j < destag->atoms[i].nimps; j++) {
				ptrdiff_t imp_offset = destag->atoms[i].imps[j] - srcag->imps;
				destag->atoms[i].imps[j] = destag->imps + imp_offset;
			}
		}
	}

	//destag->nres=srcag->nres;
	//destag->iares=srcag->iares;
	//destag->idres=srcag->idres;
	//destag->is_psf_read=srcag->is_psf_read;


	return destag;
}

struct atomgrp* copy_atomgrp (const struct atomgrp* srcag)
{
        int atomn;
        struct atomgrp* destag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
        destag->natoms = srcag->natoms;
        destag->atoms = (struct atom*) _mol_calloc (destag->natoms, sizeof (struct atom));
        for (atomn = 0; atomn < destag->natoms; atomn++)
        {
		printf("atom %d\n", atomn);
		printf("residue %s\n", srcag->atoms[atomn].residue_name);
                copy_atom (&srcag->atoms[atomn], &destag->atoms[atomn]);
        }
        return destag;
}

void ag2array( double* array, struct atomgrp* ag)
{

	int i,j;

	for(j=0; j<ag->nactives; j++){
		i=ag->activelist[j];
		array[3*j]=ag->atoms[i].X;
		array[3*j+1]=ag->atoms[i].Y;
		array[3*j+2]=ag->atoms[i].Z;
	}

}

void array2ag ( double* array, struct atomgrp* ag)
{

        int i,j;
        for(j=0; j<ag->nactives; j++){
		i= ag->activelist[j];
                ag->atoms[i].X=array[3*j];
                ag->atoms[i].Y=array[3*j+1];
                ag->atoms[i].Z=array[3*j+2];
        }

}

struct atomgrp* extract_type (struct atomgrp* ag, const char* type, struct prm* prm)
{
	return exrm_type (ag, type, prm, 0);
}

struct atomgrp* rm_type (struct atomgrp* ag, const char* type, struct prm* prm)
{
	return exrm_type (ag, type, prm, 1);
}

struct atomgrp* exrm_type (struct atomgrp* ag, const char* type, struct prm* prm, int direction)
{
	int atomn; // loop var
	int natoms; // type atom number

	struct atomgrp* exag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp)); // resultant extracted atomgrp
	exag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms); // allocate memory for the array of atoms, realloc later to make it tight

	natoms = 0;
	for (atomn = 0; atomn < ag->natoms; atomn++)
	{
		if (
				(direction == 0 && strcmp (type, prm->atoms[ag->atoms[atomn].atom_typen].typemin) == 0) // extract type
				||
				(direction == 1 && ! (strcmp (type, prm->atoms[ag->atoms[atomn].atom_typen].typemin) == 0)) // rm type
		   )
		{
			exag->atoms[natoms].atom_typen = ag->atoms[atomn].atom_typen;
			exag->atoms[natoms].sa = ag->atoms[atomn].sa;
			exag->atoms[natoms].X = ag->atoms[atomn].X;
			exag->atoms[natoms].Y = ag->atoms[atomn].Y;
			exag->atoms[natoms].Z = ag->atoms[atomn].Z;
			/*
			exag->atoms[natoms].bonds[0] = ag->atoms[atomn].bonds[0];
			exag->atoms[natoms].bonds[1] = ag->atoms[atomn].bonds[1];
			exag->atoms[natoms].bonds[2] = ag->atoms[atomn].bonds[2];
			exag->atoms[natoms].bonds[3] = ag->atoms[atomn].bonds[3];
			*/
			natoms++;
		}
	}

	if (natoms == 0)
	{
		if (direction == 0)
		{
			fprintf (stderr, "there are no atoms of type %s to extract\n", type);
		}
		if (direction == 1)
		{
			fprintf (stderr, "removing atoms of type %s leaves no remaining atoms\n", type);
		}
		exit (EXIT_FAILURE);
	}

	exag->atoms = (struct atom*) _mol_realloc (exag->atoms, sizeof (struct atom) * natoms); // realloc exag->atoms to make it tight
	exag->natoms = natoms; // attach number of atoms

	return exag;
}

struct atomgrp* extract_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm)
{
	return exrm_typemaj (ag, typemaj, prm, 0);
}

struct atomgrp* rm_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm)
{
	return exrm_typemaj (ag, typemaj, prm, 1);
}

struct atomgrp* exrm_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm, int direction)
{
	int atomn; // loop var
	int natoms; // type atom number

	struct atomgrp* exag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp)); // resultant extracted atomgrp
	exag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms); // allocate memory for the array of atoms, realloc later to make it tight

	natoms = 0;
	for (atomn = 0; atomn < ag->natoms; atomn++)
	{
		if (
				(direction == 0 && strcmp (typemaj, prm->atoms[ag->atoms[atomn].atom_typen].typemaj) == 0) // extract typemaj
				||
				(direction == 1 && ! (strcmp (typemaj, prm->atoms[ag->atoms[atomn].atom_typen].typemaj) == 0)) // rm typemaj
		   )
		{
			exag->atoms[natoms].atom_typen = ag->atoms[atomn].atom_typen;
			exag->atoms[natoms].sa = ag->atoms[atomn].sa;
			exag->atoms[natoms].X = ag->atoms[atomn].X;
			exag->atoms[natoms].Y = ag->atoms[atomn].Y;
			exag->atoms[natoms].Z = ag->atoms[atomn].Z;
			natoms++;
		}
	}

	if (natoms == 0)
	{
		if (direction == 0)
		{
			fprintf (stderr, "there are no atoms of typemaj %s to extract\n", typemaj);
		}
		if (direction == 1)
		{
			fprintf (stderr, "removing atoms of typemaj %s leaves no remaining atoms\n", typemaj);
		}
		exit (EXIT_FAILURE);
	}

	exag->atoms = (struct atom*) _mol_realloc (exag->atoms, sizeof (struct atom) * natoms); // realloc exag->atoms to make it tight
	exag->natoms = natoms; // attach number of atoms

	return exag;
}

void full_sa (struct atomgrp* ag)
{
	int atomi;
	for (atomi = 0; atomi < ag->natoms; atomi++)
	{
		ag->atoms[atomi].sa = 1;
	}
}

struct atomgrp* join_atomgrps (struct atomgrp** ags)
{
	int agai = 0; // ag atom index
	int agi = 0; // index into ags
	struct atomgrp* ag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
	ag->natoms = 100; // just a guess, realloc as necessary
	ag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms);

	while (ags[agi] != NULL)
	{
		int agsai; // ags atom index
		for (agsai = 0; agsai < ags[agi]->natoms; agsai++, agai++)
		{
			if (agai+1 > ag->natoms)
			{
				ag->natoms *= 2;
				ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);
			}

			copy_atom (&ags[agi]->atoms[agsai], &ag->atoms[agai]);
		}

		agi++;
	}


	// final realloc of the arrays to make them tight
	ag->natoms = agai;
	ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);

	return ag;
}

struct atomgrp* join_2atomgrps (struct atomgrp* ag1, struct atomgrp* ag2)
{
	struct atomgrp* ag;
	struct atomgrp** ags = (struct atomgrp**) _mol_malloc (3 * sizeof (struct atomgrp*));
	ags[0] = ag1;
	ags[1] = ag2;
	ags[2] = NULL;
	ag = join_atomgrps  (ags);
	free (ags);
	return ag;
}

void print_atomgrp (struct atomgrp* ag, struct prm* prm)
{
	int i;
	printf ("number of atoms: %d\n\n", ag->natoms);
	for (i = 0; i < ag->natoms; i++)
	{
		printf ("atom index: %d\n", i);
		printf ("\tatom type number: %d\n", ag->atoms[i].atom_typen);
		printf ("\tatom type name prefix: %s\n", prm->atoms[ag->atoms[i].atom_typen].typemaj);
		printf ("\tatom type name suffix: %s\n", prm->atoms[ag->atoms[i].atom_typen].typemin);
		printf ("\tradius: %.3f\n", prm->atoms[ag->atoms[i].atom_typen].r);
		printf ("\tcharge: %.3f\n", prm->atoms[ag->atoms[i].atom_typen].q);
		printf ("\tsa: %d\n", ag->atoms[i].sa);
		printf ("\tmask: %d\n", ag->atoms[i].mask);
		printf ("\tattl: %8.3f\n", ag->atoms[i].attl);
		printf ("\tX: %8.3f\n", ag->atoms[i].X);
		printf ("\tY: %8.3f\n", ag->atoms[i].Y);
		printf ("\tZ: %8.3f\n", ag->atoms[i].Z);
		printf ("\n");
	}
}

void check_atomgrp (struct atomgrp* ag, struct prm* prm)
{
	int atomi;
	for (atomi = 0; atomi < ag->natoms; atomi++)
	{
		if (ag->atoms[atomi].atom_typen < 0 || ag->atoms[atomi].atom_typen >= prm->natoms)
		{
			fprintf (stderr, "error: atom type number %d not defined in prm\n", ag->atoms[atomi].atom_typen);
			exit (EXIT_FAILURE);
		}
	}
}

void fill_ingrp (struct atomgrp* ag)
{
	int i;
	for (i=0; i<ag->natoms; i++)
	{
		ag->atoms[i].ingrp=i;
	}
}

void
replace_coordinates(struct atomgrp* ag, const char* pdb_path)
{
    int i;
    int *list = _mol_malloc(sizeof(int) * ag->natoms);
    struct atomgrp* ag_new = read_pdb_nopar(pdb_path);

    if (ag_new->natoms != ag->natoms) {
        _mol_error("new atomgroup has %d atoms, old has %d atoms", ag_new->natoms, ag->natoms);
        exit(EXIT_FAILURE);
    }

    for (i=0; i<ag->natoms; i++) {
        list[i] = i;
    }

    setup_subag(ag, ag_new, ag->natoms, list);

    free(list);
    free_atomgrp(ag_new);
}

void transform_atomgrpf(struct atomgrp* ag, struct mol_matrix3f rotation, struct mol_vector3f translation) {
	int i;
	for (i=0; i < ag->natoms; i++) {
		double X = ag->atoms[i].X;
		double Y = ag->atoms[i].Y;
		double Z = ag->atoms[i].Z;
		ag->atoms[i].X = X*rotation.m11 + Y*rotation.m12 + Z*rotation.m13;
		ag->atoms[i].Y = X*rotation.m21 + Y*rotation.m22 + Z*rotation.m23;
		ag->atoms[i].Z = X*rotation.m31 + Y*rotation.m32 + Z*rotation.m33;
		ag->atoms[i].X += translation.X;
		ag->atoms[i].Y += translation.Y;
		ag->atoms[i].Z += translation.Z;
	}
}

struct atomgrp* join_rec_lig_ff(struct atomgrp* rec, struct atomgrp* lig)
{
	int i;
	int j;
	struct atomgrp* ag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));

	//sum counts
	ag->natoms   = rec->natoms  + lig->natoms;
	ag->nbonds   = rec->nbonds  + lig->nbonds;
	ag->nangs    = rec->nangs   + lig->nangs;
	ag->ntors    = rec->ntors   + lig->ntors;
	ag->nimps    = rec->nimps   + lig->nimps;
	ag->num_atom_types = rec->num_atom_types + lig->num_atom_types;
	//don't think we need actives yet
	//ag->nactives = rec->nactives + lig->nactives;
	//ag->nbact    = rec->nbact   + lig->nbact;
	//ag->nangact  = rec->nangact + lig->nangact;
	//ag->ntoract  = rec->ntoract + lig->ntoract;
	//ag->nimpact  = rec->nimpact + lig->nimpact;

	//allocate first level
	ag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms);
	ag->bonds = (struct atombond*) _mol_malloc (sizeof (struct atombond) * ag->nbonds);
	ag->angs = (struct atomangle*) _mol_malloc (sizeof (struct atomangle) * ag->nangs);
	ag->tors = (struct atomtorsion*) _mol_malloc (sizeof (struct atomtorsion) * ag->ntors);
	ag->imps = (struct atomimproper*) _mol_malloc (sizeof (struct atomimproper) * ag->nimps);

	//copy rec items
	memcpy(ag->atoms, rec->atoms, sizeof(struct atom) * rec->natoms);
	memcpy(ag->bonds, rec->bonds, sizeof(struct atombond) * rec->nbonds);
	memcpy(ag->angs, rec->angs, sizeof(struct atomangle) * rec->nangs);
	memcpy(ag->tors, rec->tors, sizeof(struct atomtorsion) * rec->ntors);
	memcpy(ag->imps, rec->imps, sizeof(struct atomimproper) * rec->nimps);

	//point to correct atoms
	for (i = 0; i < rec->nbonds; i++) {
		ptrdiff_t a0_offset = ag->bonds[i].a0 - rec->atoms;
		ptrdiff_t a1_offset = ag->bonds[i].a1 - rec->atoms;
		ag->bonds[i].a0 = ag->atoms + a0_offset;
		ag->bonds[i].a1 = ag->atoms + a1_offset;
	}
	for (i = 0; i < rec->nangs; i++) {
		ptrdiff_t a0_offset = ag->angs[i].a0 - rec->atoms;
		ptrdiff_t a1_offset = ag->angs[i].a1 - rec->atoms;
		ptrdiff_t a2_offset = ag->angs[i].a2 - rec->atoms;
		ag->angs[i].a0 = ag->atoms + a0_offset;
		ag->angs[i].a1 = ag->atoms + a1_offset;
		ag->angs[i].a2 = ag->atoms + a2_offset;
	}
	for (i = 0; i < rec->ntors; i++) {
		ptrdiff_t a0_offset = ag->tors[i].a0 - rec->atoms;
		ptrdiff_t a1_offset = ag->tors[i].a1 - rec->atoms;
		ptrdiff_t a2_offset = ag->tors[i].a2 - rec->atoms;
		ptrdiff_t a3_offset = ag->tors[i].a3 - rec->atoms;
		ag->tors[i].a0 = ag->atoms + a0_offset;
		ag->tors[i].a1 = ag->atoms + a1_offset;
		ag->tors[i].a2 = ag->atoms + a2_offset;
		ag->tors[i].a3 = ag->atoms + a3_offset;
	}
	for (i = 0; i < rec->nimps; i++) {
		ptrdiff_t a0_offset = ag->imps[i].a0 - rec->atoms;
		ptrdiff_t a1_offset = ag->imps[i].a1 - rec->atoms;
		ptrdiff_t a2_offset = ag->imps[i].a2 - rec->atoms;
		ptrdiff_t a3_offset = ag->imps[i].a3 - rec->atoms;
		ag->imps[i].a0 = ag->atoms + a0_offset;
		ag->imps[i].a1 = ag->atoms + a1_offset;
		ag->imps[i].a2 = ag->atoms + a2_offset;
		ag->imps[i].a3 = ag->atoms + a3_offset;
	}

	//point to atoms to correct bonds, angs, tors, imps
	for (i = 0; i < rec->natoms; i++) {
		ag->atoms[i].name = strdup(rec->atoms[i].name);
		ag->atoms[i].residue_name = strdup(rec->atoms[i].residue_name);
		ag->atoms[i].ftype_name = strdup(rec->atoms[i].ftype_name);

		ag->atoms[i].bonds = _mol_malloc( rec->atoms[i].nbonds*sizeof(struct atombond*));
		memcpy(ag->atoms[i].bonds, rec->atoms[i].bonds, sizeof(struct atombond*) * rec->atoms[i].nbonds);
		for (j = 0; j < ag->atoms[i].nbonds; j++) {
			ptrdiff_t bond_offset = ag->atoms[i].bonds[j] - rec->bonds;
			ag->atoms[i].bonds[j] = ag->bonds + bond_offset;
		}
		ag->atoms[i].angs = _mol_malloc( rec->atoms[i].nangs*sizeof(struct atomangle*));
		memcpy(ag->atoms[i].angs, rec->atoms[i].angs, sizeof(struct atomangle*) * rec->atoms[i].nangs);
		for (j = 0; j < ag->atoms[i].nangs; j++) {
			ptrdiff_t ang_offset = ag->atoms[i].angs[j] - rec->angs;
			ag->atoms[i].angs[j] = ag->angs + ang_offset;
		}
		ag->atoms[i].tors = _mol_malloc( rec->atoms[i].ntors*sizeof(struct atomtorsion*));
		memcpy(ag->atoms[i].tors, rec->atoms[i].tors, sizeof(struct atomtorsion*) * rec->atoms[i].ntors);
		for (j = 0; j < ag->atoms[i].ntors; j++) {
			ptrdiff_t tor_offset = ag->atoms[i].tors[j] - rec->tors;
			ag->atoms[i].tors[j] = ag->tors + tor_offset;
		}
		ag->atoms[i].imps = _mol_malloc( rec->atoms[i].nimps*sizeof(struct atomimproper*));
		memcpy(ag->atoms[i].imps, rec->atoms[i].imps, sizeof(struct atomimproper*) * rec->atoms[i].nimps);
		for (j = 0; j < ag->atoms[i].nimps; j++) {
			ptrdiff_t imp_offset = ag->atoms[i].imps[j] - rec->imps;
			ag->atoms[i].imps[j] = ag->imps + imp_offset;
		}
	}


	//copy lig items
	memcpy(ag->atoms + rec->natoms, lig->atoms, sizeof(struct atom) * lig->natoms);
	memcpy(ag->bonds + rec->nbonds, lig->bonds, sizeof(struct atombond) * lig->nbonds);
	memcpy(ag->angs + rec->nangs, lig->angs, sizeof(struct atomangle) * lig->nangs);
	memcpy(ag->tors + rec->ntors, lig->tors, sizeof(struct atomtorsion) * lig->ntors);
	memcpy(ag->imps + rec->nimps, lig->imps, sizeof(struct atomimproper) * lig->nimps);

	//point to correct atoms
	struct atom *newlig_atoms = ag->atoms + rec->natoms;
	for (i = rec->nbonds; i < ag->nbonds; i++) {
		ptrdiff_t a0_offset = ag->bonds[i].a0 - lig->atoms;
		ptrdiff_t a1_offset = ag->bonds[i].a1 - lig->atoms;
		ag->bonds[i].a0 = newlig_atoms + a0_offset;
		ag->bonds[i].a1 = newlig_atoms + a1_offset;
		ag->bonds[i].ai += rec->natoms;
		ag->bonds[i].aj += rec->natoms;
	}
	for (i = rec->nangs; i < ag->nangs; i++) {
		ptrdiff_t a0_offset = ag->angs[i].a0 - lig->atoms;
		ptrdiff_t a1_offset = ag->angs[i].a1 - lig->atoms;
		ptrdiff_t a2_offset = ag->angs[i].a2 - lig->atoms;
		ag->angs[i].a0 = newlig_atoms + a0_offset;
		ag->angs[i].a1 = newlig_atoms + a1_offset;
		ag->angs[i].a2 = newlig_atoms + a2_offset;
	}
	for (i = rec->ntors; i < ag->ntors; i++) {
		ptrdiff_t a0_offset = ag->tors[i].a0 - lig->atoms;
		ptrdiff_t a1_offset = ag->tors[i].a1 - lig->atoms;
		ptrdiff_t a2_offset = ag->tors[i].a2 - lig->atoms;
		ptrdiff_t a3_offset = ag->tors[i].a3 - lig->atoms;
		ag->tors[i].a0 = newlig_atoms + a0_offset;
		ag->tors[i].a1 = newlig_atoms + a1_offset;
		ag->tors[i].a2 = newlig_atoms + a2_offset;
		ag->tors[i].a3 = newlig_atoms + a3_offset;
	}
	for (i = rec->nimps; i < ag->nimps; i++) {
		ptrdiff_t a0_offset = ag->imps[i].a0 - lig->atoms;
		ptrdiff_t a1_offset = ag->imps[i].a1 - lig->atoms;
		ptrdiff_t a2_offset = ag->imps[i].a2 - lig->atoms;
		ptrdiff_t a3_offset = ag->imps[i].a3 - lig->atoms;
		ag->imps[i].a0 = newlig_atoms + a0_offset;
		ag->imps[i].a1 = newlig_atoms + a1_offset;
		ag->imps[i].a2 = newlig_atoms + a2_offset;
		ag->imps[i].a3 = newlig_atoms + a3_offset;
	}

	//point to atoms to correct bonds, angs, tors, imps
	struct atombond *newlig_bonds = ag->bonds + rec->nbonds;
	struct atomangle *newlig_angs = ag->angs + rec->nangs;
	struct atomtorsion *newlig_tors = ag->tors + rec->ntors;
	struct atomimproper *newlig_imps = ag->imps + rec->nimps;
	for (i = rec->natoms; i < ag->natoms; i++) {
//		struct atom local_atom = ag->atoms[i];
		ag->atoms[i].atom_ftypen += rec->num_atom_types;
		ag->atoms[i].name = strdup(lig->atoms[i-(rec->natoms)].name);
		ag->atoms[i].residue_name = strdup(lig->atoms[i-(rec->natoms)].residue_name);
		ag->atoms[i].ftype_name = strdup(lig->atoms[i-(rec->natoms)].ftype_name);

		ag->atoms[i].bonds = _mol_malloc( lig->atoms[i-(rec->natoms)].nbonds*sizeof(struct atombond*));
		memcpy(ag->atoms[i].bonds, lig->atoms[i-(rec->natoms)].bonds, sizeof(struct atombond*) * lig->atoms[i-(rec->natoms)].nbonds);
		for (j = 0; j < ag->atoms[i].nbonds; j++) {
			ptrdiff_t bond_offset = ag->atoms[i].bonds[j] - lig->bonds;
			ag->atoms[i].bonds[j] = newlig_bonds + bond_offset;
		}
		ag->atoms[i].angs = _mol_malloc( lig->atoms[i-(rec->natoms)].nangs*sizeof(struct atomangle*));
		memcpy(ag->atoms[i].angs, lig->atoms[i-(rec->natoms)].angs, sizeof(struct atomangle*) * lig->atoms[i-(rec->natoms)].nangs);
		for (j = 0; j < ag->atoms[i].nangs; j++) {
			ptrdiff_t ang_offset = ag->atoms[i].angs[j] - lig->angs;
			ag->atoms[i].angs[j] = newlig_angs + ang_offset;
		}
		ag->atoms[i].tors = _mol_malloc( lig->atoms[i-(rec->natoms)].ntors*sizeof(struct atomtorsion*));
		memcpy(ag->atoms[i].tors, lig->atoms[i-(rec->natoms)].tors, sizeof(struct atomtorsion*) * lig->atoms[i-(rec->natoms)].ntors);
		for (j = 0; j < ag->atoms[i].ntors; j++) {
			ptrdiff_t tor_offset = ag->atoms[i].tors[j] - lig->tors;
			ag->atoms[i].tors[j] = newlig_tors + tor_offset;
		}
		ag->atoms[i].imps = _mol_malloc( lig->atoms[i-(rec->natoms)].nimps*sizeof(struct atomimproper*));
		memcpy(ag->atoms[i].imps, lig->atoms[i-(rec->natoms)].imps, sizeof(struct atomimproper*) * lig->atoms[i-(rec->natoms)].nimps);
		for (j = 0; j < ag->atoms[i].nimps; j++) {
			ptrdiff_t imp_offset = ag->atoms[i].imps[j] - lig->imps;
			ag->atoms[i].imps[j] = newlig_imps + imp_offset;
		}
	}

	//copy into ryan's structure
//	_mol_atom_group_copy_from_deprecated(ag);

	fill_ingrp(ag);
	
	return ag;
}
