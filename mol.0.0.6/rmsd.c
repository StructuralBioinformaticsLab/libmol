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
#include <unistd.h>

#include _MOL_INCLUDE_

/*Recursive algorithm for detecting symmetry based on bonds and atom types of the molecule.
 Input molecule should have atom numbers arranged the way that atom number k has at least one connection with previous atoms.
Basic idea is in the fact that algorithm on depth k try to find atom in the potential permutation which behaves exactly like
atom k in original numbering this is achieved by comparing bonds of potential candidate with previous atoms to the bonds of atom k
with its previous numbers. When all atoms are arranged reuired permutation is obtained
 */
struct pointlist* detectsymmetry(struct atomgrp* ag){
    struct pointlist* symmetry;
    void submit2(struct list* molecule){
	int i,j,count;
	struct list templist;
	templist.list=(int*)_mol_malloc(sizeof(int)*ag->natoms);
	for (i=0;i<ag->natoms;i++){
	    templist.list[i]=0;
	}
	for (i=0;i<molecule->n;i++){
	    templist.list[molecule->list[i]]=2;
	}
	int value=-1;
	count=0;
	for (i=0;i<ag->atoms[molecule->n].nbonds;i++){
	struct atombond* bond=  ag->atoms[molecule->n].bonds[i];
	 struct atom*  a1 =bond->a1;
	 struct atom*  a0=bond->a0;
         if (a1->ingrp==molecule->n) {a1=a0;}
	 if (a1->ingrp<molecule->n){
	 if ((a1->ingrp>value)) value=a1->ingrp;
	 templist.list[molecule->list[(a1->ingrp)]]=3;
	 count++;
         }
	}
	for (i=0;i<ag->atoms[molecule->list[value]].nbonds;i++)
	{
        struct atombond* bond=ag->atoms[molecule->list[value]].bonds[i];
	struct atom*  a1=bond->a1;
	struct atom*  a0=bond->a0;
	if (a1->ingrp==molecule->list[value]) {a1=a0;}
	if ((templist.list[a1->ingrp]==0)&&(a1->atom_ftypen==ag->atoms[molecule->n].atom_ftypen)){
	    struct atom* ac=a1;
	    int c=0;
	    int flag=0;
	    for (j=0;j<ac->nbonds;j++){
		bond= ac->bonds[j];
		a1 =bond->a1;
		a0=bond->a0;
		if (a1->ingrp==ac->ingrp) {a1=a0;}
		int k;
		for (k=0;k<molecule->n;k++) if (molecule->list[k]==a1->ingrp) {
			if (templist.list[molecule->list[k]]==3) {c++;}
			else {flag=1;}
		    }
	    }
	    if ((!flag)&& (c==count)){
		molecule->list[molecule->n]=ac->ingrp;
		if (molecule->n==(ag->natoms-1)){
		    symmetry->list[symmetry->n]=(int*)_mol_malloc(sizeof(int)*ag->natoms);
		    int* solution=(int*)symmetry->list[symmetry->n];
		    for (j=0;j<ag->natoms;j++) solution[j]=molecule->list[j];
		symmetry->n++;
	    } else {
		    molecule->n++;
		    submit2(molecule);
	     	    molecule->n--;
		}
	    }
	}
	}
	free(templist.list);
}
int i;
symmetry=_mol_malloc(sizeof(struct pointlist));
symmetry->list=_mol_malloc(sizeof(int**)*ag->natoms*ag->natoms*ag->natoms);
symmetry->n=0;
struct list molecule;
molecule.n=0;
molecule.list=(int*)_mol_malloc(sizeof(int*)*ag->natoms);
if (ag->natoms>1){
for (i=0;i<ag->natoms;i++){
    if   (ag->atoms[i].atom_ftypen==ag->atoms[0].atom_ftypen){
	//printf("Types %d %d %d\n",ag->atoms[i].atom_ftypen,i,ag->atoms[i].nbonds);
	molecule.list[molecule.n++]=i;
	submit2(&molecule);
	molecule.n--;
    }
}
} else
{
//Default solution
symmetry->list[symmetry->n]=(int*)_mol_malloc(sizeof(int)*ag->natoms);
int* solution=(int*)symmetry->list[symmetry->n];
solution[0]=0;
symmetry->n++;
}
symmetry->list=(int**)_mol_realloc(symmetry->list,symmetry->n*sizeof(int**));
return symmetry;
}

float rmsd_sym (struct atomgrp* pA, struct atomgrp* pB, struct pointlist* sym)
{
	float sum;
	int nis = pA->natoms;
	if (nis == 0)
	{
		fprintf (stderr, "error: no indices for rmsd calculation\n");
		exit (EXIT_FAILURE);
	}
	int i,j;
	int *curlist;
	float rmsd_min=1E+6;
	float rmsd_val;
	for (j=0;j<sym->n;j++){
	    curlist=(int*)sym->list[j];
	    sum=0.0;
	    for (i = 0; i < pA->natoms; i++)
	    {
		int ii = curlist[i]; // indices index
		float dev_squared =
		    ((pA->atoms[i].X - pB->atoms[ii].X)*(pA->atoms[i].X - pB->atoms[ii].X)) +
		    ((pA->atoms[i].Y - pB->atoms[ii].Y)*(pA->atoms[i].Y - pB->atoms[ii].Y)) +
		    ((pA->atoms[i].Z - pB->atoms[ii].Z)*(pA->atoms[i].Z - pB->atoms[ii].Z));

		sum += dev_squared;
	}
	    rmsd_val = sqrt (sum / (float) nis);
	    if (rmsd_val<rmsd_min) rmsd_min=rmsd_val;
	}
	return rmsd_min;
}

float rmsd_sym_no_bb(struct atomgrp *pA, struct atomgrp *pB,
		     struct pointlist *sym)
{
	int nis = pA->natoms - 5;
	if (nis == 0) {
		fprintf(stderr, "error: no indices for rmsd calculation\n");
		exit(EXIT_FAILURE);
	}
	float rmsd_min = INFINITY;
	for (int j = 0; j < sym->n; j++) {
		int *curlist = (int *)sym->list[j];
		float sum = 0.0;
		for (int i = 3; i < pA->natoms - 2; i++) {
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
		float rmsd_val = sqrtf(sum / (float)nis);
		if (rmsd_val < rmsd_min)
			rmsd_min = rmsd_val;
	}
	return rmsd_min;
}
