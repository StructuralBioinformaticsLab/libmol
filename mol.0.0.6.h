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
#ifndef _MOL_H_
#define _MOL_H_

#ifndef _PRINT_DEPRECATED_
#define _PRINT_DEPRECATED_ fprintf (stderr, "%s: Deprecated function\n", __func__);
#endif
#ifndef _mol_error
#define _mol_error(format,...) fprintf (stderr, "%s in %s@%d: " format "\n", __func__, __FILE__, __LINE__, __VA_ARGS__)
#endif
#ifndef strequal
#define strequal(s1,s2) (!strcmp(s1,s2))
#endif
#ifndef strnequal
#define strnequal(s1,s2,n) (!strncmp(s1,s2,n))
#endif

#include "mol.0.0.6/mem.h"
#include "mol.0.0.6/myhelpers.h"
#include "mol.0.0.6/prms.h"
#include "mol.0.0.6/io.h"
#include "mol.0.0.6/bond.h"
#include "mol.0.0.6/atom.h"
#include "mol.0.0.6/atom_group.h"
#include "mol.0.0.6/_atom_group_copy_from_deprecated.h"
#include "mol.0.0.6/icharmm.h"
#include "mol.0.0.6/init.h"
#include "mol.0.0.6/protein.h"
#include "mol.0.0.6/pdb.h"
#include "mol.0.0.6/ms.h"
#include "mol.0.0.6/octree.h"
#include "mol.0.0.6/matrix.h"
#include "mol.0.0.6/sasa.h"
#include "mol.0.0.6/potential.h"
#include "mol.0.0.6/energy.h"
#include "mol.0.0.6/benergy.h"
#include "mol.0.0.6/nbenergy.h"
#include "mol.0.0.6/minimize.h"
#include "mol.0.0.6/compare.h"
#include "mol.0.0.6/subag.h"
#include "mol.0.0.6/hbond.h"
#include "mol.0.0.6/version.h"
#include "mol.0.0.6/mol2.h"
#endif
