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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include _MOL_INCLUDE_

struct atomgrp* read_ms (const char* path, struct prm* prm)
{
	FILE* fp = myfopen (path, "r");

	struct atomgrp* ag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));

	char* line = NULL;
	size_t len = 0;

	// tmp scanf vals
	char atypemaj[5];
	char atypemin[5];

	// read every line of the pdb file
	int atomi = 0;
	ag->natoms = 100; // just a guess, realloc as necessary
	ag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms);
	while (getline (&line, &len, fp) != -1)
	{
		if (strncmp (line, "ATOM  ", 6) != 0) // check for ATOM line
			continue;

		if (atomi+1 > ag->natoms)
		{
			ag->natoms *= 2;
			ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);
		}

		/*
		// init bonds
		ag->atoms[atomi].bonds[0] = -1;
		ag->atoms[atomi].bonds[1] = -1;
		ag->atoms[atomi].bonds[2] = -1;
		ag->atoms[atomi].bonds[3] = -1;
		*/

		// init sa
		ag->atoms[atomi].sa = -1;
		// init mask
		ag->atoms[atomi].mask = 0;
		ag->atoms[atomi].nbondis = 0;
		ag->atoms[atomi].nbonds = 0;
		ag->atoms[atomi].nangs = 0;
		ag->atoms[atomi].ntors = 0;
		ag->atoms[atomi].nimps = 0;

		if (sscanf (line, "%*s %*d %4s %4s", atypemin, atypemaj) < 2)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "in function read_ms\n");
			fprintf (stderr, "incorrect format for ATOM line\n");
			fprintf (stderr, "at file:\n");
			fprintf (stderr, "%s\n", path);
			fprintf (stderr, "at line:\n");
			fprintf (stderr, "%s\n", line);
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		if (strncmp (atypemin, "HT", 2) == 0)
			continue; // ignore terminal hydrogens
		if (strncmp (atypemin, "OCT1", 4) == 0)
		{
			atypemin[1] = '\0';
		}
		if (strncmp (atypemin, "OCT2", 4) == 0)
			continue;


		ag->atoms[atomi].atom_typen = atomid (prm, atypemaj, atypemin);
		if (ag->atoms[atomi].atom_typen == -1) // val not found
		{
			if (atypemin[0] == 'H') // try one more time for hydrogen
			{
				atypemin[1] = '\0';
				ag->atoms[atomi].atom_typen = atomid (prm, atypemaj, atypemin);
			}
			if (ag->atoms[atomi].atom_typen == -1) // val still not found
			{
				fprintf (stderr, "error: in file %s line %s: atom type number of %s %s not defined in prm\n", path, line, atypemaj, atypemin);
				exit (EXIT_FAILURE);
			}
		}

		ag->atoms[atomi].X = atof (&line[30]);
		ag->atoms[atomi].Y = atof (&line[38]);
		ag->atoms[atomi].Z = atof (&line[46]);
		//if (line[57] != '0' && line[57] != '1') // verify that sa is 1 or 0
		if (! (line[57] >= '0' && line[57] <= '9')) // verify that sa is 1 or 0
		{
			fprintf (stderr, "error: file %s does not appear to be an ms file\n", path);
			exit (EXIT_FAILURE);
		}
		//ag->atoms[atomi].sa = atoi (&line[57]);
		ag->atoms[atomi].attl = atof (&line[57]);
		if (ag->atoms[atomi].attl > 0.0)
		{
			ag->atoms[atomi].sa = 1;
		}
		else
		{
			ag->atoms[atomi].sa = 0;
		}

		atomi++;
	}
	if (line)
		free (line);
	myfclose (fp);

	// final realloc of the arrays to make them tight
	ag->natoms = atomi;
	assert (ag->natoms > 0);
	ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);

	check_atomgrp (ag, prm);

	return ag;
}

struct atomgrp* read_ms_nopar (const char* path)
{
	FILE* fp = myfopen (path, "r");

	struct atomgrp* ag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));

	char* line = NULL;
	size_t len = 0;

	// read every line of the pdb file
	int atomi = 0;
	ag->natoms = 100; // just a guess, realloc as necessary
	ag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms);
	while (getline (&line, &len, fp) != -1)
	{
		if (strncmp (line, "ATOM  ", 6) != 0) // check for ATOM line
			continue;

		if (atomi+1 > ag->natoms)
		{
			ag->natoms *= 2;
			ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);
		}

		/*
		// init bonds
		ag->atoms[atomi].bonds[0] = -1;
		ag->atoms[atomi].bonds[1] = -1;
		ag->atoms[atomi].bonds[2] = -1;
		ag->atoms[atomi].bonds[3] = -1;
		*/

		// init sa
		ag->atoms[atomi].sa = -1;
		// init mask
		ag->atoms[atomi].mask = 0;

        ag->atoms[atomi].atom_typen = 1;
		ag->atoms[atomi].X = atof (&line[30]);
		ag->atoms[atomi].Y = atof (&line[38]);
		ag->atoms[atomi].Z = atof (&line[46]);
		//if (line[57] != '0' && line[57] != '1') // verify that sa is 1 or 0
		if (! (line[57] >= '0' && line[57] <= '9')) // verify that sa is 1 or 0
		{
			fprintf (stderr, "error: file %s does not appear to be an ms file\n", path);
			exit (EXIT_FAILURE);
		}
		//ag->atoms[atomi].sa = atoi (&line[57]);
		ag->atoms[atomi].attl = atof (&line[57]);
		if (ag->atoms[atomi].attl > 0.0)
		{
			ag->atoms[atomi].sa = 1;
		}
		else
		{
			ag->atoms[atomi].sa = 0;
		}

		atomi++;
	}
	if (line)
		free (line);
	myfclose (fp);

	// final realloc of the arrays to make them tight
	ag->natoms = atomi;
	ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);

	return ag;
}

void fprint_ms (struct atomgrp* ag, struct prm* prm, const char* path)
{
	FILE* fp = myfopen (path, "w");

	//struct atomgrp** ags = (struct atomgrp**) _mol_malloc (2 * sizeof (struct atomgrp*));
	struct atomgrp** ags = extract_nitrogen_residues (ag, prm);
	//ags[0] = copy_atomgrp (ag);
	//ags[1] = NULL; // flag the end with NULL

	int sum_atomi = 0;
	int agi = 0;
	while (ags[agi] != NULL)
	{
		int atomi;
		for (atomi = 0; atomi < ags[agi]->natoms; atomi++, sum_atomi++)
		{
			float attl;
            char atomname[5];
			fprintf (fp, "%-6s", "ATOM"); // atom number
			fprintf (fp, "%5d", sum_atomi+1); // atom number
			fprintf (fp, " "); // 1 space

            if (strlen(prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemin) == 4) {
                strcpy(atomname, prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemin);
            } else {
                atomname[0] = ' ';
                strcpy(atomname+1, prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemin);
            }
			fprintf (fp, "%-4.4s", atomname); // atom typemin
			fprintf (fp, " "); // Alternate location indicator

			if (ags[agi]->atoms[atomi].mask)
			{
				fprintf (fp, "%-3.3s", "DED"); // atom typemaj
			}
			else
			{
				fprintf (fp, "%-3.3s", prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemaj); // atom typemaj
			}

			fprintf (fp, " "); // 1 space
			fprintf (fp, "%1s", "A"); // chain id
			fprintf (fp, "%4d", agi+1); // residue number

			fprintf (fp, " "); // code for insertion of residues
			fprintf (fp, "   "); // 3 spaces

			fprintf (fp, "%8.3f", ags[agi]->atoms[atomi].X); // X coordinate
			fprintf (fp, "%8.3f", ags[agi]->atoms[atomi].Y); // Y coordinate
			fprintf (fp, "%8.3f", ags[agi]->atoms[atomi].Z); // Z coordinate

			fprintf (fp, "   "); // 3 spaces
			//fprintf (fp, "%1d", ags[agi]->atoms[atomi].sa); // sa boolean
			attl=0.0;
			if (ags[agi]->atoms[atomi].attl > 0.0)
			{
				attl += ags[agi]->atoms[atomi].attl;
			}
			if (ags[agi]->atoms[atomi].sa > 0)
			{
				attl += (float) ags[agi]->atoms[atomi].sa;
			}
			if (attl > 9.0)
			{
				fprintf (stderr, "warning: attraction level is > 9.0 for an atom\n");
			}
			fprintf (fp, "%.1f", attl); // attraction level

			//fprintf (fp, "%8d", ags[agi]->atoms[atomi].atom_typen); // atom type number

			fprintf (fp, "   "); // 3 spaces
			/*
			// bonding
			int bi = 0;
			while (bi < 4 && ags[agi]->atoms[atomi].bonds[bi] != -1)
			{
				fprintf (fp, "%6d", ags[agi]->atoms[atomi].bonds[bi]+1); // atom type number
				bi++;
			}
			*/

			fprintf (fp, "\n"); // new line
		}
		agi++;
	}

	//free (ags);
	free_atomgrps (ags);

	myfclose (fp);
}

void fprint_ms_nopar (struct atomgrp* ag, const char* inf, const char* ouf)
{
	FILE* fp = myfopen (inf, "r");
	FILE* fop = myfopen (ouf, "w");

	char* line = NULL;
	size_t len = 0;
	
	
	// read every line of the pdb file
	int atomi = 0;
	while (getline (&line, &len, fp) != -1)
	{
		float attl;
		if (strncmp (line, "ATOM  ", 6) != 0 && strncmp (line, "HETATM", 6) != 0)
		{
			fprintf (fop,"%s",line);
			continue;
		}
		sprintf(line+30,"%8.3f",ag->atoms[atomi].X );
		sprintf(line+38,"%8.3f",ag->atoms[atomi].Y );
		sprintf(line+46,"%8.3f",ag->atoms[atomi].Z );
		attl=0.0;
		if (ag->atoms[atomi].attl > 0.0)
		{
			attl += ag->atoms[atomi].attl;
		}
		if (ag->atoms[atomi].sa > 0)
		{
			attl += (float) ag->atoms[atomi].sa;
		}
		if (attl > 9.0)
		{
			fprintf (stderr, "warning: attraction level is > 9.0 for an atom\n");
		}
		sprintf(line+54,"%6.1f   \n",attl );
		atomi++;
		fprintf (fop,"%s",line);
	}
	if (line)
	        free (line);
	myfclose (fp);
	myfclose (fop);
}
