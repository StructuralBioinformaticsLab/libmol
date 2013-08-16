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

#ifdef _WIN32
#include "../mol.0.0.6.h"
#else
#include _MOL_INCLUDE_
#endif

#define MAXSLEN 200

char* main_dir ()
{
	_PRINT_DEPRECATED_

	size_t slen; // string length

	char here[] = ".";

	char* mdir = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (here);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", here);
		exit (EXIT_FAILURE);
	}

	mdir = strcpy (mdir, here);
	//mdir = strcat (mdir, "/.mol");
	mdir = strcat (mdir, "");

	return mdir;
}

char* prms_dir ()
{
	_PRINT_DEPRECATED_

	size_t slen; // string length

	char* mdir = main_dir ();
	char* pdir = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (mdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", mdir);
		exit (EXIT_FAILURE);
	}

	pdir = strcpy (pdir, mdir);
	free (mdir);
	pdir = strcat (pdir, "/prm");

	return pdir;
}

char* atom_prm_file (const char* atom_prm)
{
	_PRINT_DEPRECATED_

	size_t slen; // string length

	char* pdir = prms_dir ();
	char* atom_pfile = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (pdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", pdir);
		exit (EXIT_FAILURE);
	}
	slen = strlen (pdir) + strlen (atom_prm);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", atom_prm);
		exit (EXIT_FAILURE);
	}

	atom_pfile = strcpy (atom_pfile, pdir);
	free (pdir);
	atom_pfile = strcat (atom_pfile, "/");
	atom_pfile = strcat (atom_pfile, atom_prm);

	return atom_pfile;
}

char* coeffs_prm_file (const char* coeffs_prm)
{
	_PRINT_DEPRECATED_

	size_t slen; // string length

	char* pdir = prms_dir ();
	char* coeffs_pfile = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (pdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", pdir);
		exit (EXIT_FAILURE);
	}
	slen = strlen (pdir) + strlen (coeffs_prm);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", coeffs_prm);
		exit (EXIT_FAILURE);
	}

	coeffs_pfile = strcpy (coeffs_pfile, pdir);
	free (pdir);
	coeffs_pfile = strcat (coeffs_pfile, "/");
	coeffs_pfile = strcat (coeffs_pfile, coeffs_prm);

	return coeffs_pfile;
}

char* rots_prm_file (const char* rots_prm)
{
	_PRINT_DEPRECATED_

	size_t slen; // string length

	char* pdir = prms_dir ();
	char* rots_pfile = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (pdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", pdir);
		exit (EXIT_FAILURE);
	}
	slen = strlen (pdir) + strlen (rots_prm);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", rots_prm);
		exit (EXIT_FAILURE);
	}

	rots_pfile = strcpy (rots_pfile, pdir);
	free (pdir);
	rots_pfile = strcat (rots_pfile, "/");
	rots_pfile = strcat (rots_pfile, rots_prm);

	return rots_pfile;
}
