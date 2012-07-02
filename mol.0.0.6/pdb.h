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
#ifndef _MOL_PDB_H_
#define _MOL_PDB_H_

/**
  generate a single residue sequence combining all chains
*/
void assign_combined_residue_sequence( struct atomgrp *ag );

/**
  generate a single residue number corresponding to libmol numbering [0;nres)
*/
void assign_global_residue_number( struct atomgrp *ag );

/**
  read a pdb file
*/
struct atomgrp* read_pdb (const char* path, struct prm* prm);

/**
  read a pdb file
*/
struct atomgrp* read_pdb_nopar (const char* path);

/**
	Read the PDB file, but assign an integer id to each atom
	based on <residue-name, atom-name> instead of explicitly 
	storing the name of the atom and the residue it is part of.
	This is mainly for compatibility with other parts of the 
	code already written. However, for large PDB's this will also 
	save some space.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]
*/
struct atomgrp* read_pdb_with_compressed_typeinfo (const char* path, struct prm* prm);

/**
  read pdb model files
  rmodels specifies how much models is required, if set to -1    
*/
struct atomgrp** read_pdb_models (const char* path, struct prm* prm, int* rmodels);

/**
 * read pdb model files without requiring the atoms prms file
 * \return an array of pointers to atomgrp with models
 * @param[in] path location of the models file
 * @param[in,out] rmodels number of models required from file, if -1, read all, 
 * updated to number of models read in
 */
struct atomgrp **read_pdb_modelsnopar(const char *path, int *rmodels);

/**
  print a pdb file
*/
void fprint_pdb (struct atomgrp* ag, struct prm* prm, const char* path);

/**
  print a pdb file
*/
void write_pdb_nopar (struct atomgrp* ag, const char* inf, const char* ouf);

#endif
