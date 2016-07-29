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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include _MOL_INCLUDE_

struct atomgrp **extract_nitrogen_residues(struct atomgrp *ag, struct prm *prm)
{
	char nitrogen[2] = "N";
	struct atomgrp **ress;	// array of residues result
	int nresidues = 0;
	int atomi;		// atom index
	int natoms = 0;		// natoms in the residue
	int maxnatoms = 0;	// maximum number of atoms in a residue

	// if no atoms or first atom is not N, return ag as index 0 and NULL as index 1 of ress
	if (ag->natoms < 1
	    || strcmp(nitrogen,
		      prm->atoms[ag->atoms[0].atom_typen].typemin) != 0) {
		ress =
		    (struct atomgrp **)_mol_malloc(2 *
						   sizeof(struct atomgrp *));
		ress[0] = copy_atomgrp(ag);
		ress[1] = NULL;	// flag the end with NULL
		return ress;
	}

	for (atomi = 0; atomi < ag->natoms; atomi++)	// first count the number of residues
	{
		natoms++;
		if (strcmp
		    (nitrogen,
		     prm->atoms[ag->atoms[atomi].atom_typen].typemin) == 0) {
			nresidues++;
			if (natoms > maxnatoms) {
				maxnatoms = natoms - 1;
			}
			natoms = 1;
		}
	}
	if (natoms > maxnatoms)	// check natoms in last residue
		maxnatoms = natoms;

	if (nresidues == 0)	// no res, return ag as index 0 and NULL as index 1 of ress
	{
		ress =
		    (struct atomgrp **)_mol_malloc(2 *
						   sizeof(struct atomgrp *));
		ress[0] = copy_atomgrp(ag);
		ress[1] = NULL;	// flag the end with NULL
	} else {
		int resn = -1;	// residue number
		int resatomn = 0;	// residue atom number
		ress =
		    (struct atomgrp **)_mol_malloc((nresidues + 1) *
						   sizeof(struct atomgrp *));
		ress[nresidues] = NULL;	// flag the end with NULL

		for (atomi = 0; atomi < ag->natoms; atomi++)	// separate and store the ress
		{
			if (strcmp
			    (nitrogen,
			     prm->atoms[ag->atoms[atomi].atom_typen].typemin) ==
			    0) {
				resn++;
				if (resn != 0)
					ress[resn - 1]->natoms = resatomn;	// attach the natoms for the prev res
				resatomn = 0;	// reset

				ress[resn] =
				    (struct atomgrp *)_mol_calloc(1,
								  sizeof(struct
									 atomgrp));
				ress[resn]->atoms =
				    (struct atom *)_mol_malloc(maxnatoms *
							       sizeof(struct
								      atom));
			}

			copy_atom(&ag->atoms[atomi],
				  &ress[resn]->atoms[resatomn]);

			resatomn++;
		}
		ress[resn]->natoms = resatomn;	// attach the natoms for the last res
	}

	return ress;
}

void free_atomgrps(struct atomgrp **ress)
{
	int resn = 0;
	while (ress[resn] != NULL) {
		free_atomgrp(ress[resn]);
		resn++;
	}
	free(ress);
}

void print_residues(struct atomgrp **ress, struct prm *prm)
{
	int resi = 0;
	while (ress[resi] != NULL) {
		printf("residue index: %d\n", resi);
		print_atomgrp(ress[resi], prm);
		resi++;
	}
}

struct atomgrp *around(int nlist, int *list, struct atomgrp *ag, double distup)
{
	int i, j;
	int natoms = ag->natoms;
	double x1, y1, z1, x2, y2, z2, d2;
	double du2 = distup * distup;
	int *tmplist;
	int nout;
	struct atomgrp *destag;
	if (natoms == 0) {
		fprintf(stderr, "around> ERROR: no atoms initially");
		exit(EXIT_FAILURE);
	}
	tmplist = _mol_malloc(natoms * sizeof(int));

	for (i = 0; i < natoms; i++)
		tmplist[i] = 1;
	for (i = 0; i < nlist; i++)
		tmplist[list[i]] = 0;

	nout = 0;
	for (i = 0; i < nlist; i++) {
		x1 = ag->atoms[list[i]].X;
		y1 = ag->atoms[list[i]].Y;
		z1 = ag->atoms[list[i]].Z;
		for (j = 0; j < natoms; j++) {
			if (tmplist[j] != 1)
				continue;
			x2 = ag->atoms[j].X - x1;
			y2 = ag->atoms[j].Y - y1;
			z2 = ag->atoms[j].Z - z1;
			d2 = x2 * x2 + y2 * y2 + z2 * z2;
			if (d2 <= du2) {
				tmplist[j] = 2;
				nout++;
			}
		}
	}
	destag = (struct atomgrp *)_mol_calloc(1, sizeof(struct atomgrp));
	destag->natoms = nout;
	destag->atoms =
	    (struct atom *)_mol_malloc(sizeof(struct atom) * destag->natoms);
	j = 0;
	for (i = 0; i < natoms; i++) {
		if (tmplist[i] == 2)
			copy_atom(&ag->atoms[i], &destag->atoms[j++]);
	}
	free(tmplist);
	return destag;
}
