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
#include <math.h>
#include <float.h>

#ifndef _WIN32
#include <unistd.h>
#endif
#include _MOL_INCLUDE_

/*Recursive algorithm for detecting symmetry based on bonds and atom types of the molecule.
Basic idea is in the fact that algorithm on depth k try to find atom in the potential permutation which behaves exactly like
atom k in original numbering this is achieved by comparing bonds of potential candidate with all atoms to the bonds of atom k
. When all atoms are arranged reuired permutation is obtained, however extra check is need to make sure that atoms are connected the same way as in original arrnagement
 */
void detectsymmetry_submit2(struct list *molecule,
							int nftypes,
							struct atomgrp *ag,
							struct pointlist *symmetry) {
	int i, j;
	struct list templist;
	struct list bondstemp;
	struct list chklist;
	templist.list = (int *)_mol_malloc(sizeof(int) * ag->natoms);
	bondstemp.list = (int *)_mol_malloc(sizeof(int) * nftypes);
	for (i = 0; i < ag->natoms; i++) {
		templist.list[i] = 0;
	}
	for (i = 0; i < nftypes; i++) {
		bondstemp.list[i] = 0;
	}
	for (i = 0; i < molecule->n; i++) {
		templist.list[molecule->list[i]] = 2;
	}
	//Describe all bonds at level n (molecule->n) for original order molecule in bondstem
	//Check for all atoms not used yet, as feasible for level n
	for (i = 0; i < ag->natoms; i++){
		struct atombond* bond;
		struct atom*  a1;
		struct atom*  a0;
		//We check wether atoms is used, wether it has the same type as we need, and the same number of bonds
		if ((templist.list[ag->atoms[i].ingrp]==0)&&(ag->atoms[i].atom_ftypen==ag->atoms[molecule->n].atom_ftypen) && (ag->atoms[i].nbonds==ag->atoms[molecule->n].nbonds))
		{
			//printf("At level %d ; Atom %d looks promising\n",molecule->n,ag->atoms[i].ingrp);
			struct atom* ac=&(ag->atoms[i]);
			int k;
			int flag=0;
			for (j=0;j<nftypes;j++){
				bondstemp.list[j]=0;
			}
			//Here we check that it is bonded to the same atom type as the nth level atom
			for (j=0;j<ag->atoms[molecule->n].nbonds;j++){
				bond=  ag->atoms[molecule->n].bonds[j];
				a1 =bond->a1;
				a0=bond->a0;
				if (a1->ingrp==molecule->n) {a1=a0;}
				bondstemp.list[a1->atom_ftypen]++;
			}
			for (j=0;j<ac->nbonds;j++){
				bond= ac->bonds[j];
				a1 =bond->a1;
				a0=bond->a0;
				if (a1->ingrp==ac->ingrp) {a1=a0;}
				bondstemp.list[a1->atom_ftypen]--;
			}
			for (k=0;k<nftypes;k++){
				if  (bondstemp.list[k]!=0) flag=1;
			}
			// printf("At level %d; Flag is %d\n",molecule->n,flag); 
			if (flag==0){
				int flag2;
				molecule->list[molecule->n]=ac->ingrp;
				if (molecule->n==(ag->natoms-1))
				{
					//We have good permutation where all atoms has bonding pattern as original, however it could happen that they are not connected the same way, since atom order is arbitrary
					//we can check only now once full solution is available. We check that atoms our atom at position n is bonded, are the same atoms in the mapping.
					//printf("Checking Solution\n");
					chklist.list=(int*)_mol_malloc(sizeof(int)*ag->natoms);
					flag2=0;
					for (j=0;j<ag->natoms;j++){
						chklist.list[j]=0;
					}
					for (j=0;j<ag->natoms;j++){
						for (k=0;k<ag->atoms[j].nbonds;k++){
							bond=  ag->atoms[j].bonds[k];
							a1 =bond->a1;
							a0=bond->a0;
							if (a1->ingrp==j) {a1=a0;}
							chklist.list[molecule->list[a1->ingrp]]++;
						}
						for (k=0;k<ag->atoms[molecule->list[j]].nbonds;k++){
							bond=  ag->atoms[molecule->list[j]].bonds[k];
							a1 =bond->a1;
							a0=bond->a0;
							if (a1->ingrp==molecule->list[j]) {a1=a0;}
							chklist.list[a1->ingrp]--;
						} 
						for (k=0;k<ag->natoms;k++){
							if  (chklist.list[k]!=0) flag2=1;
						}
					}
					free(chklist.list);
					if (flag2==0){
						//printf("Solution found %d \n",flag2);
						symmetry->list[symmetry->n]=(int*)_mol_malloc(sizeof(int)*ag->natoms);
						int* solution=(int*)symmetry->list[symmetry->n];
						for (j=0;j<ag->natoms;j++) solution[j]=molecule->list[j];
						symmetry->n++;
					}
				} else {
					molecule->n++;
					// printf("Entering level %d\n",molecule->n);
					detectsymmetry_submit2(molecule, nftypes, ag, symmetry);
					molecule->n--;
					//printf("Exit level %d\n",molecule->n);
				}
			}
		}
	}
	//	printf("Exit level 2 %d\n",molecule->n);
	free(templist.list);
	free(bondstemp.list);
}
struct pointlist* detectsymmetry(struct atomgrp *ag)
{
	struct pointlist *symmetry;
	struct list molecule;
	int nftypes = -1;
	int i;
	symmetry = _mol_malloc(sizeof(struct pointlist));
	symmetry->list =
	    _mol_malloc(sizeof(int **) * ag->natoms * ag->natoms * ag->natoms);
	symmetry->n = 0;
	molecule.n = 0;
	molecule.list = (int *)_mol_malloc(sizeof(int *) * ag->natoms);
	for (i = 0; i < ag->natoms; i++) {
		if (ag->atoms[i].atom_ftypen>nftypes) {
			nftypes = ag->atoms[i].atom_ftypen;
		}
	}
	nftypes++;
	//printf ("Maxftype %d\n",nftypes);
	if (ag->natoms > 1) {
		for (i = 0; i < ag->natoms; i++) {
			if   (ag->atoms[i].atom_ftypen ==
			      ag->atoms[0].atom_ftypen) {
				//printf("Types %d %d %d\n",ag->atoms[i].atom_ftypen,i,ag->atoms[i].nbonds);
				molecule.list[molecule.n++] = i;
				detectsymmetry_submit2(&molecule, nftypes, ag, symmetry);
				molecule.n--;
			}
		}
	} else {
//Default solution
		int *solution;
		symmetry->list[symmetry->n] =
		    (int *)_mol_malloc(sizeof(int) * ag->natoms);
		solution = (int *)symmetry->list[symmetry->n];
		solution[0] = 0;
		symmetry->n++;
	}
	symmetry->list =
	    (int **)_mol_realloc(symmetry->list, symmetry->n * sizeof(int **));
	return symmetry;
}

float rmsd_sym(struct atomgrp *pA, struct atomgrp *pB, struct pointlist *sym)
{
	float sum;
	int nis = pA->natoms;
	int i, j;
	int *curlist;
	float rmsd_min = 1E+6;
	float rmsd_val;
	if (nis == 0) {
		fprintf(stderr, "error: no indices for rmsd calculation\n");
		exit(EXIT_FAILURE);
	}
	for (j = 0; j < sym->n; j++) {
		curlist = (int *)sym->list[j];
		sum = 0.0;
		for (i = 0; i < pA->natoms; i++) {
			int ii = curlist[i];	// indices index
			float dev_squared =
			    ((pA->atoms[i].X -
			      pB->atoms[ii].X) * (pA->atoms[i].X -
						  pB->atoms[ii].X)) +
			    ((pA->atoms[i].Y -
			      pB->atoms[ii].Y) * (pA->atoms[i].Y -
						  pB->atoms[ii].Y)) +
			    ((pA->atoms[i].Z -
			      pB->atoms[ii].Z) * (pA->atoms[i].Z -
						  pB->atoms[ii].Z));

			sum += dev_squared;
		}
		rmsd_val = sqrt(sum / (float)nis);
		if (rmsd_val < rmsd_min)
			rmsd_min = rmsd_val;
	}
	return rmsd_min;
}

float rmsd_sym_no_bb(struct atomgrp *pA, struct atomgrp *pB,
		     struct pointlist *sym)
{
	int nis = pA->natoms - 5;
	float rmsd_min = FLT_MAX;
	int j;
	if (nis == 0) {
		fprintf(stderr, "error: no indices for rmsd calculation\n");
		exit(EXIT_FAILURE);
	}
	for (j = 0; j < sym->n; j++) {
		int *curlist = (int *)sym->list[j];
		float sum = 0.0;
		int i;
		float rmsd_val;
		for (i = 3; i < pA->natoms - 2; i++) {
			int ii = curlist[i];	// indices index
			float dev_squared =
			    ((pA->atoms[i].X -
			      pB->atoms[ii].X) * (pA->atoms[i].X -
						  pB->atoms[ii].X)) +
			    ((pA->atoms[i].Y -
			      pB->atoms[ii].Y) * (pA->atoms[i].Y -
						  pB->atoms[ii].Y)) +
			    ((pA->atoms[i].Z -
			      pB->atoms[ii].Z) * (pA->atoms[i].Z -
						  pB->atoms[ii].Z));

			sum += dev_squared;
		}
		rmsd_val = sqrtf(sum / (float)nis);
		if (rmsd_val < rmsd_min)
			rmsd_min = rmsd_val;
	}
	return rmsd_min;
}
