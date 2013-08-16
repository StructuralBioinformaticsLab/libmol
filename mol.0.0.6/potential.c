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
#include <math.h>
#include <errno.h>

#ifdef _WIN32
#include "../mol.0.0.6.h"
#else
#include _MOL_INCLUDE_
#endif

struct matrix2df* potential_matrix2df_ncontacts_bin (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm, float r1, float r2, int only_sab)
{
	int Aatomi, Batomi; // loop iters
	float r1sq;
	float r2sq;
	struct matrix2df *M;

	if (r1 < 0.0 || r2 < 0.0)
	{
		fprintf (stderr, "begin error\n");
		fprintf (stderr, "in function potential_matrix2df_ncontacts_bin\n");
		fprintf (stderr, "at least one of the bin limits is less than 0\n");
		fprintf (stderr, "end error\n");
		exit (EXIT_FAILURE);
	}


	// squared vals for euclidean dist
	r1sq = _mol_sq(r1);
	r2sq = _mol_sq(r2);

	M = matrix2df_create (prm->natoms, prm->natoms);
	matrix2df_init (M, 0); // init all matrix vals to 0

	// loop through every atom in agA
	for (Aatomi = 0; Aatomi < agA->natoms; Aatomi++)
	{
		int Atypen;
		float AX;
		float AY;
		float AZ;

		if (only_sab && ! agA->atoms[Aatomi].sa)
			continue;

		Atypen = agA->atoms[Aatomi].atom_typen;
		if (Atypen < 0 || Atypen > prm->natoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "in function potential_matrix2df_ncontacts_bin\n");
			fprintf (stderr, "in the first atom group\n");
			fprintf (stderr, "atom type number of atom index %d is not defined in the argument atom prm\n", Aatomi);
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		AX = agA->atoms[Aatomi].X;
		AY = agA->atoms[Aatomi].Y;
		AZ = agA->atoms[Aatomi].Z;

		// loop through every atom in agB
		for (Batomi = 0; Batomi < agB->natoms; Batomi++)
		{
			int Btypen;
			float BX;
			float BY;
			float BZ;
			float rsq;
			if (only_sab && ! agB->atoms[Batomi].sa)
				continue;

			Btypen = agB->atoms[Batomi].atom_typen;

			if (Atypen < 0 || Atypen > prm->natoms-1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function potential_matrix2df_ncontacts_bin\n");
				fprintf (stderr, "in the second atom group\n");
				fprintf (stderr, "atom type number of atom index %d is not defined in the argument atom prm\n", Batomi);
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}

			BX = agB->atoms[Batomi].X;
			BY = agB->atoms[Batomi].Y;
			BZ = agB->atoms[Batomi].Z;

			// calculate euclidean distance
			rsq = (_mol_sq(AX - BX) +
			       _mol_sq(AY - BY) +
			       _mol_sq(AZ - BZ));

			if (rsq >= r1sq && rsq < r2sq) // atom distance is within the bin
			{
				// (symmetric matrix)
				M->vals[Atypen][Btypen]++;
				M->vals[Btypen][Atypen]++;
			}
		}
	}

	return M;
}

struct matrix2df* potential_matrix2df_rowcol_subatom_join (struct matrix2df* A, struct prm* prm)
{
	int i, j;
	// create subatom matrix
	struct matrix2df* B = matrix2df_create (prm->nsubatoms, prm->nsubatoms);
	matrix2df_init (B, 0); // init all matrix vals to 0

	for (i = 0; i < A->ni; i++)
	{
		int subi;
		if (i < 0 || i > prm->natoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
			fprintf (stderr, "the row %d does not have a corresponding atom type", i);
			fprintf (stderr, "in the argument atom prm\n");
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		subi = prm->atoms[i].subid; // get subatom mapping

		if (subi < 0)
			continue; // ignore this subatom type

		if (subi > prm->nsubatoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
			fprintf (stderr, "in the first atom group\n");
			fprintf (stderr, "the argument atom prm subatom mapping of atom (row) %d", i);
			fprintf (stderr, "is greater than the maximum subatom type index\n");
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		for (j = 0; j < A->nj; j++)
		{
			int subj;
			if (j < 0 || j > prm->natoms-1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
				fprintf (stderr, "the column %d does not have a corresponding atom type", j);
				fprintf (stderr, "in the argument atom prm\n");
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}

			subj = prm->atoms[j].subid; // get subatom mapping

			if (subj < 0)
				continue; // ignore this subatom type

			if (subj > prm->nsubatoms-1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
				fprintf (stderr, "in the first atom group\n");
				fprintf (stderr, "the argument atom prm subatom mapping of atom (column) %d", j);
				fprintf (stderr, "is greater than the maximum subatom type index\n");
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}

			B->vals[subi][subj] += A->vals[i][j];
		}
	}

	return B;
}
