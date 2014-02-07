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
#include <errno.h>

#include _MOL_INCLUDE_

enum File_Type
{
	FILE_PDB,
	FILE_MS, // marksur pdb
	FILE_JSON,
	FILE_UNKNOWN
};

static enum File_Type file_ext (const char* path)
{
	char* ri = strrchr(path, '.');

	if (ri == NULL)
		return FILE_UNKNOWN;
	if (strncmp (ri, ".pdb", 4) == 0)
		return FILE_PDB;
	if (strncmp (ri, ".ms", 3) == 0)
		return FILE_MS;
	if (strncmp (ri, ".json", 5) == 0)
		return FILE_JSON;

    return FILE_UNKNOWN;
}

struct atomgrp* read_file_atomgrp (const char* path, struct prm* prm, float msur_k)
{
	struct atomgrp* ag = NULL;
	if (file_ext (path) == FILE_PDB)
	{
		ag=read_pdb (path, prm);
		msur (ag, prm, msur_k);
		return ag;
	}

	if (file_ext (path) == FILE_MS)
		return read_ms (path, prm);

#ifndef _NO_JANSSON_
	if (file_ext (path) == FILE_JSON)
	{
		ag=read_json_ag(path);
		msur2(ag, msur_k);
		return ag;
	}
#endif /* _NO_JANSSON_ */

	// file type unknown
	fprintf (stderr, "file ext of %s is not a recognized file ext\n", path);
	exit (EXIT_FAILURE);

	return NULL;
}

void fprint_file_atomgrp (const char* path, struct atomgrp* ag, struct prm* prm)
{
	if (file_ext (path) == FILE_PDB)
		fprint_pdb (ag, prm, path);
		return;

	if (file_ext (path) == FILE_MS)
		fprint_ms (ag, prm, path);
		return;

	// file type unknown
	fprintf (stderr, "file ext of %s is not a recognized file ext\n", path);
	exit (EXIT_FAILURE);
}

void read_mod_vdw(char *mfile, int *nmod, int **mod, double **modeps, double **modrminh)
{
   int linesz=91;
   int na;
   char *buffer=mymalloc(sizeof(char)*linesz);
   FILE* fp = myfopen (mfile, "r");
   *nmod=0;
   while (fgets(buffer, linesz-1, fp)!=NULL)
   {
      if(!strncmp(buffer,"ATOM",4))(*nmod)++;
   }
   fclose(fp);
   *mod=mymalloc(*nmod*sizeof(int));
   *modeps=mymalloc(*nmod*sizeof(double));
   *modrminh=mymalloc(*nmod*sizeof(double));
   fp = myfopen (mfile, "r");
   na=0;
   while(fgets(buffer, linesz-1, fp)!=NULL)
   {
      if(!strncmp(buffer,"ATOM",4))
      {
         (*mod)[na]=atoi(buffer+4)-1;
	 (*modeps)[na]=atof(buffer+54);
	 (*modrminh)[na]=atof(buffer+60);
         na++;
      }
   }
   free(buffer);
   fclose(fp);
}
