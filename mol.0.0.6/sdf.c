/*
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
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include _MOL_INCLUDE_

/* Documentation for the sdf file format can be found at
 * http://accelrys.com/products/informatics/cheminformatics/ctfile-formats/no-fee.php
 */

struct atomgrp **read_sdf(const char *path, int *rmodels)
{
	FILE *fp = myfopen(path, "r");
	char *line = NULL;
	size_t len = 0;
	getline(&line, &len, fp);
	getline(&line, &len, fp);
	getline(&line, &len, fp);
	int line_length = getline(&line, &len, fp);
	if (line_length < 39) {
		fprintf(stderr, "No SDF version found\n");
		return NULL;
	}
	char version[5];
	strncpy(version, line+34, 5);
	fclose(fp);
	if (line)
		free(line);

	if ( strncmp(version, "V3000", 5) == 0) {
		return read_sdf_v3000(path, rmodels);
	} else { //assume v2000
		return read_sdf_v2000(path, rmodels);
	}
	return NULL;
}
		
struct atomgrp **read_sdf_v2000(const char *path, int *rmodels)
{
	FILE *fp = myfopen(path, "r");
	int nmodels = 100;	//Initial guess
	struct atomgrp **ag_models =
	    _mol_malloc(nmodels * sizeof(struct atomgrp *));
	char *line = NULL;
	size_t len = 0;
	int modeli = -1;
	while (getline(&line, &len, fp) != -1) {
		if ((modeli + 1 > (*rmodels - 1)) && ((*rmodels) != -1))	//If were over requested models and we don't read everyhting
			break;
		if (modeli + 1 > nmodels) {
			nmodels *= 2;
			ag_models =
			    _mol_realloc(ag_models,
					 nmodels * sizeof(struct atomgrp *));
		}
		modeli++;
		ag_models[modeli] = _mol_calloc(1, sizeof(struct atomgrp));

		getline(&line, &len, fp);
		getline(&line, &len, fp);
		getline(&line, &len, fp);

		

		char tmp[4]; //temporary holding for variables
		tmp[3] = '\0';
		strncpy(tmp, line, 3); //copy the number of atoms into tmp

		errno = 0;
		ag_models[modeli]->natoms = atoi(tmp);
		if (errno) {
			perror("atoi");
			exit(EXIT_FAILURE);
		}
		ag_models[modeli]->atoms =
		    _mol_malloc(sizeof(struct atom) *
				ag_models[modeli]->natoms);
                ag_models[modeli]->prm = NULL;

		for (int atomi = 0; atomi < ag_models[modeli]->natoms; atomi ++) {
			if (getline(&line, &len, fp) == -1) {
				fprintf(stderr, "Not enough atom lines in sdf file\n");
			}
			ag_models[modeli]->atoms[atomi].name = _mol_calloc(4, sizeof(char));
			sscanf(line, "%10lf%10lf%10lf %3s", 
			       &(ag_models[modeli]->atoms[atomi].X),
			       &(ag_models[modeli]->atoms[atomi].Y),
			       &(ag_models[modeli]->atoms[atomi].Z),
			       ag_models[modeli]->atoms[atomi].name );

		}
		int status;
		do {
			status=getline(&line, &len, fp);
		} while ( (status != -1) && (strncmp(line, "$$$$", 4) != 0) );
	}
	if (line)
		free(line);
	myfclose(fp);
	nmodels = modeli + 1;
	ag_models = _mol_realloc(ag_models, nmodels * sizeof(struct atomgrp));
	// final realloc of the arrays to make them tight

	*rmodels = modeli + 1;
	return ag_models;
}

struct atomgrp **read_sdf_v3000(const char *path, int *rmodels)
{
	FILE *fp = myfopen(path, "r");
	int nmodels = 100;	//Initial guess
	struct atomgrp **ag_models =
	    _mol_malloc(nmodels * sizeof(struct atomgrp *));
	char *line = NULL;
	size_t len = 0;
	int modeli = -1;
	while (getline(&line, &len, fp) != -1) {
		if ((modeli + 1 > (*rmodels - 1)) && ((*rmodels) != -1))	//If were over requested models and we don't read everyhting
			break;
		if (modeli + 1 > nmodels) {
			nmodels *= 2;
			ag_models =
			    _mol_realloc(ag_models,
					 nmodels * sizeof(struct atomgrp *));
		}
		modeli++;
		ag_models[modeli] = _mol_calloc(1, sizeof(struct atomgrp));

		getline(&line, &len, fp);
		getline(&line, &len, fp);
		getline(&line, &len, fp);
		int status;
		do {
			status=getline(&line, &len, fp);
		} while ( (status != -1) && (strncmp(line, "M  V30 BEGIN CTAB", 17) != 0) );

		//next should be counts line, eg M  V30 COUNTS 50 50 0 0 1
		getline(&line, &len, fp);
		sscanf(line, "M  V30 COUNTS %d", &(ag_models[modeli]->natoms));

		ag_models[modeli]->atoms =
		    _mol_malloc(sizeof(struct atom) * ag_models[modeli]->natoms);
                ag_models[modeli]->prm = NULL;
		
		do {
			status=getline(&line, &len, fp);
		} while ( (status != -1) && (strncmp(line, "M  V30 BEGIN ATOM", 17) != 0) );

		for (int atomi = 0; atomi < ag_models[modeli]->natoms; atomi ++) {
			if (getline(&line, &len, fp) == -1) {
				fprintf(stderr, "Not enough atom lines in sdf file\n");
			}
			sscanf(line, "M  V30 %*d %ms %lf %lf %lf", 
			       &(ag_models[modeli]->atoms[atomi].name),
			       &(ag_models[modeli]->atoms[atomi].X),
			       &(ag_models[modeli]->atoms[atomi].Y),
			       &(ag_models[modeli]->atoms[atomi].Z) );

		}
		do {
			status=getline(&line, &len, fp);
		} while ( (status != -1) && (strncmp(line, "$$$$", 4) != 0) );
	}
	if (line)
		free(line);
	myfclose(fp);
	nmodels = modeli + 1;
	ag_models = _mol_realloc(ag_models, nmodels * sizeof(struct atomgrp));
	// final realloc of the arrays to make them tight

	*rmodels = modeli + 1;
	return ag_models;
}
