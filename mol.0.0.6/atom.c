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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "atom.h"
#include "mem.h"

// warning: this function does not handle all of mol_atom
void mol_atom_create(mol_atom * a, int nbondis)
{
	assert(a != NULL);

	a->nbondis = nbondis;
	if (a->nbondis > 0)
		a->bondis = _mol_malloc(a->nbondis * sizeof(int));
	else
		a->bondis = NULL;

	a->nbonds = 0;
	a->nangs = 0;
	a->ntors = 0;
	a->nimps = 0;

	return;
}

// warning: this function does not handle all of mol_atom
void mol_atom_destroy(mol_atom * a)
{
	assert(a != NULL);

	if (a->nbondis > 0)
		free(a->bondis);
	if (a->nbonds > 0)
		free(a->bonds);
	if (a->nangs > 0)
		free(a->angs);
	if (a->ntors > 0)
		free(a->tors);
	if (a->nimps > 0)
		free(a->imps);

	return;
}

// warning: this function does not handle all of mol_atom
void mol_atom_copy(mol_atom * asrc, mol_atom * adst)
{
	int i;

	assert(asrc != NULL);
	assert(adst != NULL);
	assert(asrc->nbondis == adst->nbondis);
	if (adst->nbondis > 0)
		assert(adst->bondis != NULL);

	adst->ingrp = asrc->ingrp;
	adst->atom_typen = asrc->atom_typen;
	adst->atom_ftypen = asrc->atom_ftypen;

	adst->sa = asrc->sa;
	adst->fixed = asrc->fixed;
	adst->mask = asrc->mask;
	adst->attl = asrc->attl;

	adst->X = asrc->X;
	adst->Y = asrc->Y;
	adst->Z = asrc->Z;

	adst->acp_type = asrc->acp_type;
	adst->rminh = asrc->rminh;
	adst->chrg = asrc->chrg;

	for (i = 0; i < adst->nbondis; i++)
		adst->bondis[i] = asrc->bondis[i];

	return;
}

void free_atom(struct atom *atm)
{
	free(atm->name);
	free(atm->residue_name);
	free(atm->ftype_name);
	//free(atm->element);
	free(atm->bonds);
	free(atm->angs);
	free(atm->tors);
	free(atm->imps);
}

void fullcopy_atom(struct atom *src, struct atom *dest)
{
	int i;
	if (src == NULL || dest == NULL) {
		fprintf(stderr,
			"err in function copy_atom: src or dest arg is NULL\n");
		exit(EXIT_FAILURE);
	}

	dest->ingrp = src->ingrp;
	dest->atom_typen = src->atom_typen;
	dest->atom_ftypen = src->atom_ftypen;
	dest->sa = src->sa;
	dest->fixed = src->fixed;
	dest->attl = src->attl;
	dest->mask = src->mask;
	dest->name = strdup(src->name);
	dest->residue_name = strdup(src->residue_name);

	dest->X = src->X;
	dest->Y = src->Y;
	dest->Z = src->Z;

	dest->GX = src->GX;
	dest->GY = src->GY;
	dest->GZ = src->GZ;

	dest->eps = src->eps;
	dest->rminh = src->rminh;
	dest->eps03 = src->eps03;
	dest->rminh03 = src->rminh03;
	dest->acevolume = src->acevolume;
	dest->chrg = src->chrg;
	dest->B = src->B;
	dest->nbonds = src->nbonds;
	dest->bonds = _mol_malloc(dest->nbonds * sizeof(struct atombond *));
	for (i = 0; i < dest->nbonds; i++)
		dest->bonds[i] = src->bonds[i];
	dest->nangs = src->nangs;
	dest->angs = _mol_malloc(dest->nangs * sizeof(struct atomangle *));
	for (i = 0; i < dest->nangs; i++)
		dest->angs[i] = src->angs[i];
	dest->ntors = src->ntors;
	dest->tors = _mol_malloc(dest->ntors * sizeof(struct atomtorsion *));
	for (i = 0; i < dest->ntors; i++)
		dest->tors[i] = src->tors[i];
	dest->nimps = src->nimps;
	dest->imps = _mol_malloc(dest->nimps * sizeof(struct atomimproper *));
	for (i = 0; i < dest->nimps; i++)
		dest->imps[i] = src->imps[i];
}

void copy_atom(struct atom *src, struct atom *dest)
{
	if (src == NULL || dest == NULL) {
		fprintf(stderr,
			"err in function copy_atom: src or dest arg is NULL\n");
		exit(EXIT_FAILURE);
	}

	dest->atom_typen = src->atom_typen;
	dest->sa = src->sa;
	dest->attl = src->attl;
	dest->mask = src->mask;
	dest->fixed = src->fixed;

	dest->X = src->X;
	dest->Y = src->Y;
	dest->Z = src->Z;
	if (src->name != NULL) {
		dest->name = strdup(src->name);
	}
	if (src->residue_name != NULL) {
		dest->residue_name = strdup(src->residue_name);
	}

	/*
	   dest->bonds[0] = src->bonds[0];
	   dest->bonds[1] = src->bonds[1];
	   dest->bonds[2] = src->bonds[2];
	   dest->bonds[3] = src->bonds[3];
	 */
}

void fprintf_stderr_atom(struct atom *atom, struct prm *prm)
{
	fprintf(stderr, "\tatom type number: %d\n", atom->atom_typen);
	fprintf(stderr, "\tatom type name prefix: %s\n",
		prm->atoms[atom->atom_typen].typemaj);
	fprintf(stderr, "\tatom type name suffix: %s\n",
		prm->atoms[atom->atom_typen].typemin);
	fprintf(stderr, "\tcharge: %.2f\n", prm->atoms[atom->atom_typen].q);
	fprintf(stderr, "\tsa: %d\n", atom->sa);
	fprintf(stderr, "\tX: %8.3f\n", atom->X);
	fprintf(stderr, "\tY: %8.3f\n", atom->Y);
	fprintf(stderr, "\tZ: %8.3f\n", atom->Z);
	fprintf(stderr, "\n");
}

void free_springset(struct springset *sprst)
{
	int i;
	for (i = 0; i < sprst->nsprings; i++)
		free(sprst->springs[i].laspr);
	if (sprst->springs != NULL)
		free(sprst->springs);
	if (sprst != NULL)
		free(sprst);
}
