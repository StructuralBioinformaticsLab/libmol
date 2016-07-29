/*
Copyright (c) 2013, Acpharis
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

#include <math.h>

/* This is based on Lone-Pair Directionality in Hydrogen Bond Potential Functions
   by Vedani and Dunitz
   J. Am. Chem. Soc. 1985, 107, 7653-7658
   See Tables I, II, and footnote 14
*/

static double mol_yeti_pairwise_hbondeng(struct atomgrp *ag, int donor_h,
					 int acceptor, double rc2)
{
	mol_atom *hydro = &(ag->atoms[donor_h]);
	mol_atom *acc = &(ag->atoms[acceptor]);

	double dx = hydro->X - acc->X;
	double dy = hydro->Y - acc->Y;
	double dz = hydro->Z - acc->Z;

	double d2 = dx * dx + dy * dy + dz * dz;

	if (d2 > rc2)
		return 0;

	double d_h_a_component;
	double h_a_aa_component;
	double displacement_component = 1.0;

	/* calculate j-i coefficients */
	double A;
	double C;

	/*
	   for this case, we do i=12, j= 10, based on yeti paper
	   (A/r^i - C/r^j)
	   general form:
	   A = (j/(j-i)) * E_0 * (r_0^i)
	   C = (i/(j-i)) * E_0 * (r_0^j)

	   In this case:
	   A = -5 * E_0 * (r_0^i)
	   C = -6 * E_0 * (r_0^j)
	 */

}

void mol_yeti_hbondeng(struct atomgrp *ag, double *energy, struct nblist *nblst)
{
	for (i = 0; i < nblst->nfat; i++) {
		int ai = nblst->ifat[i];
		mol_atom *atom_i = &(ag->atoms[ai]);
		int n2;
		int *p;
		int j;

		if (!(atom_i->hprop & DONATABLE_HYDROGEN)
		    && (atom_i->yeti_type == MOL_YETI_NONE))
			continue;

		n2 = nblst->nsat[i];
		p = nblst->isat[i];

		for (j = 0; j < n2; j++) {
			int aj = p[j];
			mol_atom *atom_j = &(ag->atoms[aj]);

			if (atom_i->hprop & DONATABLE_HYDROGEN) {
				if (atom_j->yeti_type == MOL_YETI_NONE)
					continue;

				(*energy) +=
				    mol_yeti_pairwise_hbondeng(ag, ai, aj, rc2);
			} else {
				if (!(atom_j->hprop & DONATABLE_HYDROGEN))
					continue;

				(*energy) +=
				    mol_yeti_pairwise_hbondeng(ag, aj, ai, rc2);
			}
		}
	}
}
