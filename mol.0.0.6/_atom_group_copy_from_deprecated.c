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
#include <assert.h>
#include <stdlib.h>

#include _MOL_INCLUDE_

void
_mol_atom_create_bond_indices (mol_atom* a, int nbondis)
{
	assert (a != NULL);

	a->nbondis = nbondis;
	if (a->nbondis > 0)
		a->bondis = _mol_malloc (a->nbondis * sizeof (int));
	else
		a->bondis = NULL;

	return;
}

void
_mol_atom_copy_bond_ptrs_to_bond_indices (mol_atom* atom, int nbonds, mol_bond* bonds)
{
	int bondssi, bondsi;
	mol_bond** bondss;

	atom->nbondis = atom->nbonds;
	bondss = atom->bonds;

	for (bondssi = 0; bondssi < atom->nbonds; bondssi++)
	{
		int target_count = 0;
		mol_bond* bond_target = bondss[bondssi];

		for (bondsi = 0; bondsi < nbonds; bondsi++)
		{
			if (&bonds[bondsi] == bond_target)
			{
				atom->bondis[bondssi] = bondsi;
				target_count++;
			}
		}

		assert (target_count == 1);
	}
}

void
_mol_bond_copy_atom_ptrs_to_atom_indices (mol_bond* b)
{
	b->ai = b->a0->ingrp;
	b->aj = b->a1->ingrp;
}

void
_mol_atom_group_copy_from_deprecated (mol_atom_group* ag)
{
	int atomsi, bondsi;

	for (atomsi = 0; atomsi < ag->natoms; atomsi++)
	{
		_mol_atom_copy_bond_ptrs_to_bond_indices (&ag->atoms[atomsi],
				ag->nbonds, ag->bonds);
	}

	for (bondsi = 0; bondsi < ag->nbonds; bondsi++)
	{
		_mol_bond_copy_atom_ptrs_to_atom_indices (&ag->bonds[bondsi]);
	}
}
