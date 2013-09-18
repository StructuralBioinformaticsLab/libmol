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
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
#include "../mol.0.0.6.h"
#else
#include _MOL_INCLUDE_
#endif

#define ACE_NONPOLAR 1 // 01
#define ACE_POLAR    2 // 10
#define ACE_ALL   (ACE_NONPOLAR | ACE_POLAR) // 11

typedef unsigned char ACE_ENERGY_TYPE;

//Init ACE
void ace_ini(struct atomgrp* ag,struct acesetup* ac_s)
{
	double minr=0.6;
	int i,j,k;
	double  sqrtpi=sqrt(M_PI);
	double  sqrttwo=sqrt(2);
	double   c1 = 4.0 / (3.0*M_PI);
	double   c2 = 77.0 * M_PI *sqrttwo/ 512.0;
	double   c3 = 2.0 * M_PI * sqrtpi;
	double   pi2 = 1.0 / (M_PI*M_PI);
	double width=1.2;
	double ri,ri2,rk,vk,rk2,alpha,alpha2,alpha4,prod2,prod4,ratio,tik2,qik,qterm,temp,fik,omik,s2ik,s3ik,uik,omgik;
	int nbsize;
	ac_s->ntypes=ag->num_atom_types+1;
	ac_s->rsolv=(double*)_mol_malloc(sizeof(double)*ac_s->ntypes);
	ac_s->vsolv=(double*)_mol_malloc(sizeof(double)*ac_s->ntypes);
	ac_s->lwace=(double*)_mol_malloc(sizeof(double)*ac_s->ntypes);
	ac_s->s2ace=(double*)_mol_malloc(sizeof(double)*ac_s->ntypes*ac_s->ntypes);
	ac_s->uace=(double*)_mol_malloc(sizeof(double)*ac_s->ntypes*ac_s->ntypes);
	ac_s->wace=(double*)_mol_malloc(sizeof(double)*ac_s->ntypes*ac_s->ntypes);
	ac_s->hydr=(double*)_mol_malloc(sizeof(double)*ac_s->ntypes);
	for (i=0;i<ac_s->ntypes;i++){
		ac_s->vsolv[i]=-1.0;
		ac_s->rsolv[i]=0.0;
	}
	for (i=0;i<ag->natoms;i++){
		if (ac_s->vsolv[ag->atoms[i].atom_ftypen]<=0.0){
			struct atombond bp;
			ac_s->vsolv[ag->atoms[i].atom_ftypen]=ag->atoms[i].acevolume;
			ac_s->rsolv[ag->atoms[i].atom_ftypen]=_mol_max(ag->atoms[i].rminh,minr);
			//Hydrogen correction
			if (strcmp("H", ag->atoms[i].ftype_name) == 0) {
				ac_s->rsolv[ag->atoms[i].atom_ftypen]=_mol_max(0.85, ac_s->rsolv[ag->atoms[i].atom_ftypen]);
			} else if (strcmp("HC", ag->atoms[i].ftype_name) == 0) {
				ac_s->rsolv[ag->atoms[i].atom_ftypen]=_mol_max(0.83, ac_s->rsolv[ag->atoms[i].atom_ftypen]);
			}
			//	printf(" %d %.3f \n",ag->atoms[i].atom_ftypen,ag->atoms[i].rminh);
			for (j=0;j<ag->nbonds;j++){
				bp=ag->bonds[j];
				if (bp.a0->atom_ftypen == ag->atoms[i].atom_ftypen) {
					ac_s->rsolv[ag->atoms[i].atom_ftypen]=_mol_max(ac_s->rsolv[ag->atoms[i].atom_ftypen],bp.a1->rminh - bp.l0);
				}
				if (ag->bonds[j].a1->atom_ftypen==ag->atoms[i].atom_ftypen){
					ac_s->rsolv[ag->atoms[i].atom_ftypen]=_mol_max(ac_s->rsolv[ag->atoms[i].atom_ftypen],bp.a0->rminh - bp.l0);
				}
			}
		}
	}
	//Hydrogen correction
//	if (ac_s->rsolv[1]<0.85) ac_s->rsolv[1]=0.85; //H   ?
//	if (ac_s->rsolv[2]<0.83) ac_s->rsolv[2]=0.83; //HC  ?
	for (i=0;i<ac_s->ntypes;i++){
		if (ac_s->vsolv[i]>=0.0){
			//	printf ("%d %.3f %.3f\n",i,ac_s->vsolv[i],ac_s->rsolv[i]);
		}
	}
	for (i=0;i<ac_s->ntypes;i++){
		ac_s->hydr[i]=-4*M_PI*0.003*ac_s->rsolv[i]*(ac_s->rsolv[i]+1.4)*(ac_s->rsolv[i]+1.4);
		ri = ac_s->rsolv[i];
		ri2 = ri * ri;
		if (ri>0){
			for (k=0;k<ac_s->ntypes;k++){
				rk = ac_s->rsolv[k];
				vk = ac_s->vsolv[k];
				if (rk>0){
					rk2 = rk * rk;
					alpha = _mol_max(width,ri/rk);
					alpha2 = alpha * alpha;
					alpha4 = alpha2 * alpha2;
					prod2 = alpha2 * rk2;
					prod4 = prod2 * prod2;
					ratio = alpha2 * rk2 / ri2;
					tik2 = 0.5 * M_PI * ratio;
					temp = 1.0 / (1.0+2.0*tik2);
					fik = 2.0/(1.0+tik2) - temp;
					qik = tik2 * sqrt(temp);
					qterm = qik - atan(qik);
					if (k == i) {
						omik = -((vk*pi2)-c1 *(ri2*ri))* qterm / (alpha4 * ri2*ri2);
						ac_s->lwace[i]   = omik;
					}
					if (vk > 1.0e-8){
						omgik = vk * qterm * pi2 / prod4;
						s2ik = 3.0 * qterm * prod2/((3.0+fik)*qik-4.0*atan(qik));
						s3ik = s2ik * sqrt(s2ik);
						uik = c2 * ri / (1.0-(c3*s3ik*ri*omgik/vk));
						ac_s->wace[i*ac_s->ntypes+k] = omgik;
						ac_s->s2ace[i*ac_s->ntypes+k] = s2ik;
						ac_s->uace[i*ac_s->ntypes+k] = uik;
					}
					else
					{
						ac_s->wace[i*ac_s->ntypes+k] = 0.0;
						ac_s->s2ace[i*ac_s->ntypes+k]= 1.0e-6;
						ac_s->uace[i*ac_s->ntypes+k] = 1.0e+6;
					}
				}
			}
		}
	}
	for (i=0;i<ac_s->ntypes;i++){
		if (ac_s->vsolv[i]>=0.0){
			for (k=0;k<ac_s->ntypes;k++){
				if (ac_s->vsolv[k]>=0.0){
					//	printf ("%d %d %.3f %.3f\n",i,k,ac_s->wace[i*ac_s->ntypes+k],ac_s->lwace[i]);
				}
			}
		}
	}
	ac_s->nbsize=1;
	nbsize=1;
	ac_s->eself=_mol_calloc(sizeof(double), ag->natoms);
	ac_s->rborn=_mol_calloc(sizeof(double), ag->natoms);
	ac_s->swarr=_mol_calloc(sizeof(double), nbsize);
	ac_s->dswarr=_mol_calloc(sizeof(double), nbsize);
	ac_s->darr=_mol_calloc(sizeof(double), nbsize);
	//d(rb)/d(eself)
	ac_s->dbrdes=_mol_calloc(sizeof(double), ag->natoms);
	ac_s->xsf=_mol_calloc(sizeof(double), 2*nbsize);
	ac_s->ysf=_mol_calloc(sizeof(double), 2*nbsize);
	ac_s->zsf=_mol_calloc(sizeof(double), 2*nbsize);
	ac_s->xf=_mol_calloc(sizeof(double), ag->natoms);
	ac_s->yf=_mol_calloc(sizeof(double), ag->natoms);
	ac_s->zf=_mol_calloc(sizeof(double), ag->natoms);
	ac_s->diarr=_mol_calloc(sizeof(double), ag->natoms);
	ac_s->list0123=_mol_calloc(sizeof(int), 1);
}

static int*  compute_0123_list(struct atomgrp* ag, int* n0123, int *list03,int n03,int *list02,int n02,int *na01,int **pna01)
{
	int i,j,n2,ind;
	int* p;
	int *list0123;
	*n0123=(ag->nbonds+n03+n02);
	list0123=(int*)_mol_malloc(*n0123*2*sizeof(int));
	ind=0;
	//list03 goes first
	for(i=0; i<n03; i++)
	{
		list0123[ind++]=list03[2*i];
		list0123[ind++]=list03[2*i+1];
	}
	for(i=0; i<ag->natoms; i++)
	{
		n2=na01[i];
		p=pna01[i];
		for (j=0;j<n2;j++){
			if (p[j]>i) {
				list0123[ind++]=i;
				list0123[ind++]=p[j];
			}
		}
	}
	for(i=0; i<n02; i++)
	{
		list0123[ind++]=list02[2*i];
		list0123[ind++]=list02[2*i+1];
	}


	return list0123;
}

static void ace_eselfupdate(int i1,int i2, int it, int kt,int ij,double dx, double dy, double dz, struct acesetup* ac_s, double *eself,double* swarr, double* dswarr,double* darr, double* xf, double* yf, double* zf, double* xsf, double* ysf, double* zsf, double rul3, double rul12, double nb2cot, double nb2cof,int nbsize)
{
	double sw,dsw, r,r2,r3,r4,expterm,temp,rmu,term,rl,ru,u4ace,ffk;
	r2 = dx*dx + dy*dy + dz*dz;
	darr[ij]=-1;

	sw=1.0;
	dsw=0;
	if (r2>nb2cot) {
		rl  = (nb2cot) - r2;
		ru  = (nb2cof) - r2;
		sw  = ru*ru*(ru-3*rl)*rul3;
		dsw = rl*ru*rul12;
	}

	swarr[ij]=sw;
	dswarr[ij]=dsw;
	r=sqrt(r2);
	darr[ij]=r;
	r3 = r2 * r;
	r4 = r2 * r2;
	if (ac_s->vsolv[kt]>0){
		expterm = ac_s->wace[it*ac_s->ntypes+kt] * exp(-r2/ac_s->s2ace[it*ac_s->ntypes+kt]);
		u4ace=pow(ac_s->uace[it*ac_s->ntypes+kt],4);
		rmu = r4 + u4ace;
		term = (ac_s->vsolv[kt]/(8.0*M_PI)) * pow((r3/rmu),4);

		temp = 2.0*(expterm+term);
		eself[i1] -= temp*sw;
		ffk = ((-8*term*(3*u4ace-r4)/(r2*rmu)) + (4*expterm/ac_s->s2ace[it*ac_s->ntypes+kt]))*sw - temp*dsw;

		xsf[ij]=ffk*dx;
		ysf[ij]=ffk*dy;
		zsf[ij]=ffk*dz;

		xf[i1]-=xsf[ij];
		yf[i1]-=ysf[ij];
		zf[i1]-=zsf[ij];
	} else {
		xsf[ij]=0;
		ysf[ij]=0;
		zsf[ij]=0;
	}
	if (ac_s->vsolv[it]>0){
		int ind;
		expterm = ac_s->wace[kt*ac_s->ntypes+it] * exp(-r2/ac_s->s2ace[kt*ac_s->ntypes+it]);
		u4ace=pow(ac_s->uace[kt*ac_s->ntypes+it],4);
		rmu = r4 + u4ace;
		term = (ac_s->vsolv[it]/(8.0*M_PI)) * pow((r3/rmu),4);
		temp = 2.0*(expterm+term);
		eself[i2] -= temp*sw;
		ffk = ((-8*term*(3*u4ace-r4)/(r2*rmu)) + (4*expterm/ac_s->s2ace[kt*ac_s->ntypes+it]))*sw - temp*dsw;

		ind=ij+nbsize;
		xsf[ind]=-ffk*dx;
		ysf[ind]=-ffk*dy;
		zsf[ind]=-ffk*dz;

		xf[i2]-=xsf[ind];
		yf[i2]-=ysf[ind];
		zf[i2]-=zsf[ind];
	} else {
		xsf[ij+nbsize]=0;
		ysf[ij+nbsize]=0;
		zsf[ij+nbsize]=0;
	}
}


//Fixed update
void ace_fixedupdate(struct atomgrp* ag, struct agsetup* ags, struct acesetup* ac_s)
{
	int *list01=_mol_malloc(2*(ag->nbonds)*sizeof(int));
	int *na01=_mol_malloc((ag->natoms)*sizeof(int));
	int **pna01=_mol_malloc((ag->natoms)*sizeof(int*));
	int* list0123;
	int n0123;
	comp_list01(ag, list01, na01, pna01);
	if (ac_s->list0123!=NULL) free(ac_s->list0123);
	list0123=compute_0123_list(ag, &n0123, ags->listf03,ags->nf03,ags->list02,ags->n02,na01,pna01);
	ac_s->list0123=list0123;
	ac_s->n0123=n0123;
	free(pna01);
	free(na01);
	free(list01);
}
void ace_updatenblst(struct agsetup* ags, struct acesetup* ac_s)
{
	int nbsize=0;
	int i;
	for(i=0; i<ags->nblst->nfat; i++)
	{
		nbsize+=ags->nblst->nsat[i];
	}
	nbsize+=ac_s->n0123;
	ac_s->nbsize=nbsize;
	ac_s->swarr=_mol_realloc(ac_s->swarr,nbsize*sizeof(double));
	ac_s->dswarr=_mol_realloc(ac_s->dswarr,nbsize*sizeof(double));
	ac_s->darr=_mol_realloc(ac_s->darr,nbsize*sizeof(double));
	ac_s->xsf=_mol_realloc(ac_s->xsf,2*nbsize*sizeof(double));
	ac_s->ysf=_mol_realloc(ac_s->ysf,2*nbsize*sizeof(double));
	ac_s->zsf=_mol_realloc(ac_s->zsf,2*nbsize*sizeof(double));

}

//Free ace
void destroy_acesetup(struct acesetup* ac_s)
{
	free(ac_s->list0123);
	free(ac_s->eself);
	free(ac_s->rborn);
	free(ac_s->swarr);
	free(ac_s->dswarr);
	free(ac_s->darr);
	free(ac_s->dbrdes);
	free(ac_s->xsf);
	free(ac_s->ysf);
	free(ac_s->zsf);
	free(ac_s->xf);
	free(ac_s->yf);
	free(ac_s->zf);
	free(ac_s->diarr);
	free(ac_s->lwace);
	free(ac_s->rsolv);
	free(ac_s->vsolv);
	free(ac_s->s2ace);
	free(ac_s->uace);
	free(ac_s->wace);
	free(ac_s->hydr);
}
void free_acesetup(struct acesetup* ac_s)
{
	destroy_acesetup(ac_s);
	free(ac_s);
}

void test_acegrads(struct atomgrp *ag,struct agsetup* ags, struct acesetup* acs,double d)
{
	int n=ag->natoms, i;
	double en, en1, t;
	double *fs=_mol_malloc(3*n*sizeof(double));
	//en0
	en=0;
	aceeng(ag, &en,acs,ags);
	//        vdweng03(f, ag, &en, n03, list03);
	//        eleng03(f, ag, 1.0, &en, n03, list03);
	//        elengs03(f, nblst->nbcof, ag, 1.0, &en, n03, list03);
	//        vdweng(ag, &en, nblst);
	//        eleng(ag, 1.0, &en, nblst);

	for(i=0; i<n; i++)
	{
		//x
		en1=0;
		t=ag->atoms[i].X;
		ag->atoms[i].X=d+t;
		aceeng(ag, &en1,acs,ags);
		//vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
		//                vdweng03(f, ag, &en1, n03, list03);
		//                eleng03(f, ag, 1.0, &en1, n03, list03);
		//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
		//                vdweng(ag, &en1, nblst);
		//                eleng(ag, 1.0, &en1, nblst);
		ag->atoms[i].X=t;
		fs[3*i]=(en-en1)/d;
		//y
		en1=0;
		t=ag->atoms[i].Y;
		ag->atoms[i].Y=d+t;
		aceeng(ag, &en1,acs,ags);
		// vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
		//                vdweng03(f, ag, &en1, n03, list03);
		//                eleng03(f, ag, 1.0, &en1, n03, list03);
		//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
		//                vdweng(ag, &en1, nblst);
		//                eleng(ag, 1.0, &en1, nblst);
		ag->atoms[i].Y=t;
		fs[3*i+1]=(en-en1)/d;
		//z
		en1=0;
		t=ag->atoms[i].Z;
		ag->atoms[i].Z=d+t;
		aceeng(ag, &en1,acs,ags);
		//vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
		//                vdweng03(f, ag, &en1, n03, list03);
		//                eleng(f, ag, 1.0, &en1, n03, list03);
		//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
		//                vdweng(ag, &en1, nblst);
		//                eleng(ag, 1.0, &en1, nblst);
		ag->atoms[i].Z=t;
		fs[3*i+2]=(en-en1)/d;
	}
	en=0;
	zero_grads(ag);
	aceeng(ag, &en,acs,ags);
	// vdwengs03(f, nblst->nbcof, ag, &en, n03, list03);
	//        vdweng03(f, ag, &en, n03, list03);
	//        eleng(f, ag, 1.0, &en, n03, list03);
	//        elengs03(f, nblst->nbcof, ag, 1.0, &en, n03, list03);
	//        vdweng(ag, &en, nblst);
	//        eleng(ag, 1.0, &en, nblst);
	for(i=0; i<n; i++)
	{
		/*	    ag->atoms[i].GX=fs[3*i];
			    ag->atoms[i].GY=fs[3*i+1];
			    ag->atoms[i].GZ=fs[3*i+2];
		 */
		printf("PREF %d calculated: %.5lf %.5lf %.5lf\n",
				i,ag->atoms[i].GX,ag->atoms[i].GY,ag->atoms[i].GZ);
		printf("PREF %d numerical : %.5lf %.5lf %.5lf\n",
				i,fs[3*i],fs[3*i+1],fs[3*i+2]);
		//ag->atoms[i].GX=fs[3*i];
		//  ag->atoms[i].GY=fs[3*i+1];
		// ag->atoms[i].GZ=fs[3*i+2];

	}
	printf("Etotal %.9lf\n",en);
	free(fs);
}

void test_acegrads_nonpolar(struct atomgrp *ag,struct agsetup* ags, struct acesetup* acs,double d)
{
	int n=ag->natoms, i;
	double en, en1, t;
	double *fs=_mol_malloc(3*n*sizeof(double));
	//en0
	en=0;
	aceeng_nonpolar(ag, &en,acs,ags);
	//        vdweng03(f, ag, &en, n03, list03);
	//        eleng03(f, ag, 1.0, &en, n03, list03);
	//        elengs03(f, nblst->nbcof, ag, 1.0, &en, n03, list03);
	//        vdweng(ag, &en, nblst);
	//        eleng(ag, 1.0, &en, nblst);

	for(i=0; i<n; i++)
	{
		//x
		en1=0;
		t=ag->atoms[i].X;
		ag->atoms[i].X=d+t;
		aceeng_nonpolar(ag, &en1,acs,ags);
		//vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
		//                vdweng03(f, ag, &en1, n03, list03);
		//                eleng03(f, ag, 1.0, &en1, n03, list03);
		//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
		//                vdweng(ag, &en1, nblst);
		//                eleng(ag, 1.0, &en1, nblst);
		ag->atoms[i].X=t;
		fs[3*i]=(en-en1)/d;
		//y
		en1=0;
		t=ag->atoms[i].Y;
		ag->atoms[i].Y=d+t;
		aceeng_nonpolar(ag, &en1,acs,ags);
		// vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
		//                vdweng03(f, ag, &en1, n03, list03);
		//                eleng03(f, ag, 1.0, &en1, n03, list03);
		//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
		//                vdweng(ag, &en1, nblst);
		//                eleng(ag, 1.0, &en1, nblst);
		ag->atoms[i].Y=t;
		fs[3*i+1]=(en-en1)/d;
		//z
		en1=0;
		t=ag->atoms[i].Z;
		ag->atoms[i].Z=d+t;
		aceeng_nonpolar(ag, &en1,acs,ags);
		//vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
		//                vdweng03(f, ag, &en1, n03, list03);
		//                eleng(f, ag, 1.0, &en1, n03, list03);
		//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
		//                vdweng(ag, &en1, nblst);
		//                eleng(ag, 1.0, &en1, nblst);
		ag->atoms[i].Z=t;
		fs[3*i+2]=(en-en1)/d;
	}
	en=0;
	zero_grads(ag);
	aceeng_nonpolar(ag, &en,acs,ags);
	// vdwengs03(f, nblst->nbcof, ag, &en, n03, list03);
	//        vdweng03(f, ag, &en, n03, list03);
	//        eleng(f, ag, 1.0, &en, n03, list03);
	//        elengs03(f, nblst->nbcof, ag, 1.0, &en, n03, list03);
	//        vdweng(ag, &en, nblst);
	//        eleng(ag, 1.0, &en, nblst);
	for(i=0; i<n; i++)
	{
		/*	    ag->atoms[i].GX=fs[3*i];
			    ag->atoms[i].GY=fs[3*i+1];
			    ag->atoms[i].GZ=fs[3*i+2];
		 */
		printf("PREF %d calculated: %.5lf %.5lf %.5lf\n",
				i,ag->atoms[i].GX,ag->atoms[i].GY,ag->atoms[i].GZ);
		printf("PREF %d numerical : %.5lf %.5lf %.5lf\n",
				i,fs[3*i],fs[3*i+1],fs[3*i+2]);
		//ag->atoms[i].GX=fs[3*i];
		//  ag->atoms[i].GY=fs[3*i+1];
		// ag->atoms[i].GZ=fs[3*i+2];

	}
	printf("Etotal %.9lf\n",en);
	free(fs);
}

static void aceeng_internal(struct atomgrp* ag,double* en,struct acesetup* ac_s,
                            struct agsetup* ags, ACE_ENERGY_TYPE ace_energy_type)
{
	double x1,y1,z1,dx,dy,dz;
	int i1,i2,it,j,i,n2,kt,ij=0;
	int *p;
	double b0=0;
	double nb2cot=8.0*8.0;//Switching start
	double nb2cof=ags->nblst->nbcof*ags->nblst->nbcof;
	int nbsize=ac_s->nbsize;
	double etotal=0;
	double ecoul=0;
	int* list0123=ac_s->list0123;
	int n0123=ac_s->n0123;
	double* eself=ac_s->eself;
	double* rborn=ac_s->rborn;
	double* swarr=ac_s->swarr;
	double* dswarr=ac_s->dswarr;
	double* darr=ac_s->darr;
	//d(rb)/d(eself)
	double* dbrdes=ac_s->dbrdes;
	double* xsf=ac_s->xsf;
	double* ysf=ac_s->ysf;
	double* zsf=ac_s->zsf;
	double* xf=ac_s->xf;
	double* yf=ac_s->yf;
	double* zf=ac_s->zf;
	double* diarr=ac_s->diarr;
	double rul3=1.0/pow((nb2cof-nb2cot),3.0);
	double rul12=12.0*rul3;
	double kelec;
	double factor_E;
	double tau;
	double ehydr = 0;
	double fac1;
	double facc1;
	//NBLST RBORN
	for (i=0;i<ag->natoms;i++){
		double ri;
		it =ag->atoms[i].atom_ftypen;
		ri=ac_s->rsolv[it];
		eself[i] = 1.0/ri + 2.0*ac_s->lwace[it];
		b0 = b0 +ag->atoms[i].acevolume;
		xf[i]=0;
		yf[i]=0;
		zf[i]=0;
	}
	b0 = pow((0.75*b0/M_PI),1.0/3.0);
	ij=0;
	//Loop through non bonded atoms
	for(i=0; i<ags->nblst->nfat; i++) {
		i1=ags->nblst->ifat[i];
		x1=ag->atoms[i1].X;
		y1=ag->atoms[i1].Y;
		z1=ag->atoms[i1].Z;
		n2=ags->nblst->nsat[i];

		p=ags->nblst->isat[i];
		it=ag->atoms[i1].atom_ftypen;
		for(j=0; j<n2; j++)
		{
			i2=p[j];
			dx=x1-ag->atoms[i2].X;
			dy=y1-ag->atoms[i2].Y;
			dz=z1-ag->atoms[i2].Z;
			kt = ag->atoms[i2].atom_ftypen;
			ace_eselfupdate(i1,i2, it, kt,ij++,dx, dy, dz, ac_s, eself,swarr, dswarr,darr, xf, yf, zf,xsf,ysf,zsf,rul3,rul12,nb2cot,nb2cof,nbsize);
		}
	}
	//Loop through 1-2-3-4 list
	for(i=0; i<n0123; i++) {
		i1=list0123[2*i];
		i2=list0123[2*i+1];
		x1=ag->atoms[i1].X;
		y1=ag->atoms[i1].Y;
		z1=ag->atoms[i1].Z;
		it=ag->atoms[i1].atom_ftypen;
		dx=x1-ag->atoms[i2].X;
		dy=y1-ag->atoms[i2].Y;
		dz=z1-ag->atoms[i2].Z;
		kt = ag->atoms[i2].atom_ftypen;
		ace_eselfupdate(i1,i2, it, kt,ij++,dx, dy, dz, ac_s, eself,swarr, dswarr,darr, xf, yf, zf,xsf,ysf,zsf,rul3,rul12,nb2cot,nb2cof,nbsize);
	}
	//Electrostatic constant need to carry over to constants
	kelec=332.0716;
	factor_E=-kelec/2.0;
	//change to epsilons;
	tau=((1/4.0)-(1/78.0));

	for (i=0;i<ag->natoms;i++) {
		double c2;
		diarr[i]=0;
		if (ace_energy_type & ACE_NONPOLAR) {
			it=ag->atoms[i].atom_ftypen;
			ehydr+=-ac_s->hydr[it]*eself[i];
			diarr[i]+=-(ac_s->hydr[it]/factor_E);
		}
		if (ace_energy_type & ACE_POLAR) {
			if (eself[i] >= 1.0/b0) {
				rborn[i] = 1.0 / eself[i];
				dbrdes[i]=rborn[i]/(factor_E*eself[i]);
			}
			else {
				rborn[i] = 1.0 * b0 * (2.0-b0*eself[i]);
				dbrdes[i]=b0*b0/factor_E;
			}
			c2=ag->atoms[i].chrg*ag->atoms[i].chrg*tau;
			eself[i]=factor_E*eself[i]*c2;
			etotal+=eself[i];
			diarr[i]+=c2;
		}
	}
	if (ace_energy_type & ACE_POLAR) {
		//Second loop calculating energy nb
		ij=0;
		fac1=-kelec*tau;
		facc1=kelec/(4.0);
		for(i=0; i<ags->nblst->nfat; i++) {
			i1=ags->nblst->ifat[i];
			n2=ags->nblst->nsat[i];
			p=ags->nblst->isat[i];
			x1=ag->atoms[i1].X;
			y1=ag->atoms[i1].Y;
			z1=ag->atoms[i1].Z;
			if (ag->atoms[i1].chrg!=0) {
				for(j=0; j<n2; j++) {
					double s,s2,brij,expo,fexp,rij2,rij,cij,fac2,sw,dsw,fac3,fac4,dij,dji,facc2,fac5,fx,fy,fz;
					i2=p[j];
					if ((ag->atoms[i2].chrg!=0) && (darr[ij]>0)){
						dx=x1-ag->atoms[i2].X;
						dy=y1-ag->atoms[i2].Y;
						dz=z1-ag->atoms[i2].Z;
						s     = darr[ij];
						s2    = s*s;
						brij  = rborn[i1]*rborn[i2];
						expo  = s2/(4.0*brij);
						fexp  = exp(-expo);
						rij2  = s2 + brij*fexp;
						rij   = sqrt(rij2);
						cij  = ag->atoms[i1].chrg*ag->atoms[i2].chrg;
						fac2 = fac1*cij/rij;
						sw    = swarr[ij];
						dsw   = dswarr[ij];
						etotal+= fac2*sw;
						fac3 = fac2/rij2;
						fac4 = 0.5*fac3*(1+expo)*fexp;
						fac4 = fac4*sw;
						dij   = fac4*rborn[i2]*dbrdes[i1];
						dji   = fac4*rborn[i1]*dbrdes[i2];
						diarr[i1]+=dij;
						diarr[i2]+=dji;
						//full field force (incl. coulomb) for non-bonded pairs
						facc2= facc1*cij/s;
						fac5 = fac3*(0.25*fexp-1.0)-facc2/s2;
						ecoul +=facc2*sw;
						fac5 = sw*fac5 + (fac2+facc2)*dsw;
						fx   = fac5*dx;
						fy   = fac5*dy;
						fz   = fac5*dz;
						ag->atoms[i1].GX+=fx;
						ag->atoms[i1].GY+=fy;
						ag->atoms[i1].GZ+=fz;
						ag->atoms[i2].GX-=fx;
						ag->atoms[i2].GY-=fy;
						ag->atoms[i2].GZ-=fz;
					}
					ij++;
				}
			} else {
				ij+=n2;
			}
		}
			//1-2-3-4 second loop
		for(i=0; i<n0123; i++) {
			double s,s2,brij,expo,fexp,rij2,rij,cij,fac2,sw,dsw,fac3,fac4,dij,dji,facc2,fac5,fx,fy,fz;
			i1=list0123[2*i];
			i2=list0123[2*i+1];
			x1=ag->atoms[i1].X;
			y1=ag->atoms[i1].Y;
			z1=ag->atoms[i1].Z;
			dx=x1-ag->atoms[i2].X;
			dy=y1-ag->atoms[i2].Y;
			dz=z1-ag->atoms[i2].Z;
			if ((darr[ij]>0)&& (ag->atoms[i1].chrg!=0)&& (ag->atoms[i2].chrg!=0)){
				s     = darr[ij];
				s2    = s*s;
				brij  = rborn[i1]*rborn[i2];
				expo  = s2/(4.0*brij);
				fexp  = exp(-expo);
				rij2  = s2 + brij*fexp;
				rij   = sqrt(rij2);
				cij  = ag->atoms[i1].chrg*ag->atoms[i2].chrg;
				fac2 = fac1*cij/rij;
				sw    = swarr[ij];
				dsw   = dswarr[ij];
				etotal+= fac2*sw;
				fac3 = fac2/rij2;
				fac4 = 0.5*fac3*(1+expo)*fexp;
				fac4 = fac4*sw;
				dij   = fac4*rborn[i2]*dbrdes[i1];
				dji   = fac4*rborn[i1]*dbrdes[i2];
				diarr[i1]+=dij;
				diarr[i2]+=dji;
				//1-4 switching
				if (i<ags->nf03){
					facc2= ac_s->efac*facc1*cij/s;
					fac5 = fac3*(0.25*fexp-1.0)-facc2/s2;
					ecoul +=facc2*sw;
					fac5 = sw*fac5 + (fac2+facc2)*dsw;
				}else{
					//1-3 1-2 interactions
					fac5 = fac3*(0.25*fexp-1.0);
					fac5 = sw*fac5 + (fac2)*dsw;
				}

				fx   = fac5*dx;
				fy   = fac5*dy;
				fz   = fac5*dz;
				ag->atoms[i1].GX+=fx;
				ag->atoms[i1].GY+=fy;
				ag->atoms[i1].GZ+=fz;
				ag->atoms[i2].GX-=fx;
				ag->atoms[i2].GY-=fy;
				ag->atoms[i2].GZ-=fz;
			}
			ij++;
		}
	}
	*en=ecoul+etotal+ehydr;

	for (i=0;i<ag->natoms;i++) {
		double fdiarr=-factor_E*diarr[i];
		ag->atoms[i].GX+=fdiarr*xf[i];
		ag->atoms[i].GY+=fdiarr*yf[i];
		ag->atoms[i].GZ+=fdiarr*zf[i];
	}
	ij=0;
	for(i=0; i<ags->nblst->nfat; i++) {
		i1=ags->nblst->ifat[i];
		n2=ags->nblst->nsat[i];
		p=ags->nblst->isat[i];
		for(j=0; j<n2; j++) {
			if (darr[ij]>0) {
				double fdiarr1;
				double fdiarr2;
				i2=p[j];
				fdiarr1=-factor_E*diarr[i1];
				fdiarr2=-factor_E*diarr[i2];
				ag->atoms[i2].GX+=fdiarr1*xsf[ij];
				ag->atoms[i2].GY+=fdiarr1*ysf[ij];
				ag->atoms[i2].GZ+=fdiarr1*zsf[ij];
				ag->atoms[i1].GX+=fdiarr2*xsf[ij+nbsize];
				ag->atoms[i1].GY+=fdiarr2*ysf[ij+nbsize];
				ag->atoms[i1].GZ+=fdiarr2*zsf[ij+nbsize];
			}
			ij++;
		}
	}
	//Loop through 1-2-3-4 list
	for(i=0; i<n0123; i++) {
		i1=list0123[2*i];
		i2=list0123[2*i+1];
		if (darr[ij]>0){
			double fdiarr1=-factor_E*diarr[i1];
			double fdiarr2=-factor_E*diarr[i2];
			ag->atoms[i2].GX+=fdiarr1*xsf[ij];
			ag->atoms[i2].GY+=fdiarr1*ysf[ij];
			ag->atoms[i2].GZ+=fdiarr1*zsf[ij];
			ag->atoms[i1].GX+=fdiarr2*xsf[ij+nbsize];
			ag->atoms[i1].GY+=fdiarr2*ysf[ij+nbsize];
			ag->atoms[i1].GZ+=fdiarr2*zsf[ij+nbsize];
		}
		ij++;
	}
	for (i=0;i<ag->natoms;i++) {
		ag->atoms[i].GX*=-1;
		ag->atoms[i].GY*=-1;
		ag->atoms[i].GZ*=-1;
	}


}

void aceeng(struct atomgrp* ag,double* en,struct acesetup* ac_s,struct agsetup* ags)
{
	aceeng_internal(ag, en, ac_s, ags, ACE_ALL);
}
void aceeng_nonpolar(struct atomgrp* ag,double* en,struct acesetup* ac_s,struct agsetup* ags)
{
	aceeng_internal(ag, en, ac_s, ags, ACE_NONPOLAR);
}
void aceeng_polar(struct atomgrp* ag,double* en,struct acesetup* ac_s,struct agsetup* ags)
{
	aceeng_internal(ag, en, ac_s, ags, ACE_POLAR);
}
