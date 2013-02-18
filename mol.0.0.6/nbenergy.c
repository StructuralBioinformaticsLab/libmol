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
#include <math.h>
#include <errno.h>

#include _MOL_INCLUDE_

#define MAXPAIR 36
#define CCELEC 332.0716


void test_nbgrads(struct atomgrp *ag, double d, struct nblist *nblst,
                  int n03, int* list03)
{
int n=ag->natoms, i;
        double en, en1, t, f=0.5;
	double *fs=_mol_malloc(3*n*sizeof(double));
//en0
        en=0;
        vdwengs03(f, nblst->nbcof, ag, &en, n03, list03);
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
                vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
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
                vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
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
                vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
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
        vdwengs03(f, nblst->nbcof, ag, &en, n03, list03);
//        vdweng03(f, ag, &en, n03, list03);
//        eleng(f, ag, 1.0, &en, n03, list03);
//        elengs03(f, nblst->nbcof, ag, 1.0, &en, n03, list03);
//        vdweng(ag, &en, nblst);
//        eleng(ag, 1.0, &en, nblst);
        for(i=0; i<n; i++)
        {
                printf("%d calculated: %lf %lf %lf\n",
                i,ag->atoms[i].GX,ag->atoms[i].GY,ag->atoms[i].GZ);
                printf("%d numerical : %lf %lf %lf\n",
                i,fs[3*i],fs[3*i+1],fs[3*i+2]);
        }
        free(fs);
}

void vdweng03(double f, struct atomgrp *ag, double* ven, int n03, int* list03)
{
        int i, i1, i2;
        struct atom *a1, *a2;
        double ei, ej, eij, ri, rj, rij;
        double dx, dy, dz, d2;
        double Rd6, Rd12, dven, g;
        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           ei=f*(a1->eps03);
           ri=a1->rminh03;

           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           ej=a2->eps03;
           rj=a2->rminh03;

           eij=ei*ej;
           rij=ri+rj;
           rij*=rij;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;
           Rd6=rij/d2;
           Rd6=Rd6*Rd6*Rd6;
           Rd12=Rd6*Rd6;
           (*ven)+=eij*(Rd12-2*Rd6);
           dven=-eij*12*(-Rd12+Rd6)/d2;
           g=dven*dx;
           (a1->GX)+=g;
           (a2->GX)-=g;
           g=dven*dy;
           (a1->GY)+=g;
           (a2->GY)-=g;
           g=dven*dz;
           (a1->GZ)+=g;
           (a2->GZ)-=g;
        }
}

void vdwengs03(double f, double rc, struct atomgrp *ag, double* ven,
               int n03, int* list03)
{
        int i, i1, i2;
        double ei, ri, ej, rj, dx, dy, dz, eij;
        double d2, Rd12, Rd6, Rr12, Rr6, dr6;
        double rij, dven, g;
        struct atom *a1, *a2;
        double rc2=rc*rc;

        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           i2=list03[2*i+1];
           if(i1 == i2)continue;
           a1=&(ag->atoms[i1]);
           ei=f*(a1->eps03);
           ri=a1->rminh03;

           a2=&(ag->atoms[i2]);
           ej=a2->eps03;
           rj=a2->rminh03;

           eij=ei*ej;
           rij=ri+rj;
           rij*=rij;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;

           Rd6=rij/d2;
           Rd6=Rd6*Rd6*Rd6;
           Rd12=Rd6*Rd6;
           Rr6=rij/rc2;
           Rr6=Rr6*Rr6*Rr6;
           Rr12=Rr6*Rr6;
           dr6=d2/rc2;
           dr6=dr6*dr6*dr6;

           (*ven)+=eij*(Rd12-2*Rd6+Rr6*(4.0-2*dr6)+Rr12*(2*dr6-3.0));
           dven=-eij*12*(-Rd12+Rd6+dr6*(Rr12-Rr6))/d2;
           g=dven*dx;
           (a1->GX)+=g;
           (a2->GX)-=g;
           g=dven*dy;
           (a1->GY)+=g;
           (a2->GY)-=g;
           g=dven*dz;
           (a1->GZ)+=g;
           (a2->GZ)-=g;
	}
}

void vdweng(struct atomgrp *ag, double* ven, struct nblist *nblst)
{
        int i, i1, n2, *p, j, i2;
        double ei, ri, ej, rj, x1, y1, z1, dx, dy, dz, eij;
        double d2, Rd12, Rd6, Rr12, Rr6, dr6;
        double rij, dven, g;
        struct atom *a1, *a2;
        double rc=nblst->nbcof;
        double rc2=rc*rc;

        for(i=0; i<nblst->nfat; i++)
        {
           i1=nblst->ifat[i];
           a1=&(ag->atoms[i1]);
           ei=a1->eps;
           ri=a1->rminh;
           x1=ag->atoms[i1].X;
           y1=ag->atoms[i1].Y;
           z1=ag->atoms[i1].Z;
           n2=nblst->nsat[i];
           p=nblst->isat[i];
           for(j=0; j<n2; j++)
           {
              i2=p[j];
              a2=&(ag->atoms[i2]);
              rj=a2->rminh;
              ej=a2->eps;
              dx=x1-ag->atoms[i2].X;
              dy=y1-ag->atoms[i2].Y;
              dz=z1-ag->atoms[i2].Z;
              d2=dx*dx+dy*dy+dz*dz;
              if(d2<rc2)
              {
                 eij=ei*ej;
                 rij=ri+rj;
                 rij*=rij;

                 Rd6=rij/d2;
                 Rd6=Rd6*Rd6*Rd6;
                 Rd12=Rd6*Rd6;
                 Rr6=rij/rc2;
                 Rr6=Rr6*Rr6*Rr6;
                 Rr12=Rr6*Rr6;
                 dr6=d2/rc2;
                 dr6=dr6*dr6*dr6;


                 (*ven)+=eij*(Rd12-2*Rd6+Rr6*(4.0-2*dr6)+Rr12*(2*dr6-3.0));
                 dven=-eij*12*(-Rd12+Rd6+dr6*(Rr12-Rr6))/d2;
                 g=dven*dx;
                 (a1->GX)+=g;
                 (a2->GX)-=g;
                 g=dven*dy;
                 (a1->GY)+=g;
                 (a2->GY)-=g;
                 g=dven*dz;
                 (a1->GZ)+=g;
                 (a2->GZ)-=g;
              }
           }
        }
}



void eleng03(double f, struct atomgrp *ag, double eps, double* een,
             int n03, int* list03)
{
        int i, i1, i2;
        struct atom *a1, *a2;
        double chi, chj, chij, ee, deen, g;
        double dx, dy, dz, d, d2;
        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           chi=a1->chrg;

           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           chj=a2->chrg;

           chij=f*CCELEC*chi*chj/eps;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;
           d=sqrt(d2);
           ee=chij/d;
           (*een)+=ee;
           deen=ee/d2;
           g=deen*dx;
           (a1->GX)+=g;
           (a2->GX)-=g;
           g=deen*dy;
           (a1->GY)+=g;
           (a2->GY)-=g;
           g=deen*dz;
           (a1->GZ)+=g;
           (a2->GZ)-=g;
        }
}

void elengs03(double f, double rc, struct atomgrp *ag, double eps, double* een,
              int n03, int* list03)
{
        int i, i1, i2;
        double dx, dy, dz;
        double d2, d1, g, ch1, ch2;
        double esh, desh;
        struct atom *a1, *a2;

        double pf=f*CCELEC/eps;
        double rc2=1.0/(rc*rc);
        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           ch1=pf*(a1->chrg);

           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           ch2=a2->chrg;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;
           d1=sqrt(d2);
           if(d1<rc)
           {
              esh=1.0-d1/rc;
              esh*=esh/d1;
              desh=(1.0/d2-rc2)/d1;
              ch2=ch2*ch1;
              (*een)+=ch2*esh;
              g=ch2*desh*dx;
              (a1->GX)+=g;
              (a2->GX)-=g;
              g=ch2*desh*dy;
              (a1->GY)+=g;
              (a2->GY)-=g;
              g=ch2*desh*dz;
              (a1->GZ)+=g;
              (a2->GZ)-=g;
           }
        }
}


void eleng(struct atomgrp *ag, double eps, double* een, struct nblist *nblst)
{
        int i, j, i1, n2, i2, *p;
        double dx, dy, dz;
        double d2, d1, g, ch1, ch2, x1, y1, z1;
        double esh, desh;
        struct atom *a1, *a2;

        double pf=CCELEC/eps;
        double rc=nblst->nbcof;
        double rc2=1.0/(rc*rc);
	for(i=0; i<nblst->nfat; i++)
	{
	   i1=nblst->ifat[i];
           a1=&(ag->atoms[i1]);
           ch1=pf*(a1->chrg);
           x1=ag->atoms[i1].X;
           y1=ag->atoms[i1].Y;
           z1=ag->atoms[i1].Z;
           n2=nblst->nsat[i];
           p=nblst->isat[i];
           for(j=0; j<n2; j++)
	   {
	      i2=p[j];
              a2=&(ag->atoms[i2]);
              ch2=a2->chrg;
              dx=x1-ag->atoms[i2].X;
              dy=y1-ag->atoms[i2].Y;
              dz=z1-ag->atoms[i2].Z;
              d2=dx*dx+dy*dy+dz*dz;
              d1=sqrt(d2);
              if(d1<rc)
              {
                 esh=1.0-d1/rc;
                 esh*=esh/d1;
                 desh=(1.0/d2-rc2)/d1;
                 ch2=ch2*ch1;
                 (*een)+=ch2*esh;
                 g=ch2*desh*dx;
                 (a1->GX)+=g;
                 (a2->GX)-=g;
                 g=ch2*desh*dy;
                 (a1->GY)+=g;
                 (a2->GY)-=g;
                 g=ch2*desh*dz;
                 (a1->GZ)+=g;
                 (a2->GZ)-=g;
              }
	   }
	}
}

void destroy_nblist(struct nblist *nblst)
{
        for(int i=0; i<nblst->nfat; i++)
        {
           if(nblst->nsat[i]>0)
              free(nblst->isat[i]);
        }
        free(nblst->isat);
        free(nblst->nsat);
        free(nblst->ifat);
        free(nblst->crds);
}
void free_nblist(struct nblist *nblst)
{
	destroy_nblist(nblst);
	free(nblst);
}

//! Generate a nonbonded list nblst from cluster and cube subdivisions and exclusion list arrays.
/*! To generate nblist
    1. Loop over the filled cubes
    2. Loop over linked clusters in the first cube
    3. Loop over neighboring cubes including first cube.
    4. Loop through the linked list clusters of the second cube.
    5. Go only for pairs of clusters which are closer than an nbcut plus lagest distance from any atom to the
       geometry center of cluster 1 plus the largest distance of the cluster 2
    6. For the pair of clusters icl1 and icl2 go through all pairs of their atoms
    7. Write the pair to the list if the distance between clusters is less than nbcut minus largest dist of
       cluster one minus largest dist of cluster two.
    8. Check the atomic pair distance  if the condition in the step 7 is not true
    9. Write the pair if distance less than nonbonded cutoff.

    The structure of the nblist is:
    number of first atoms,
    array of first atom indices,
    array of numbers of second atoms for each first one,
    array of pointers to arrays of all second atoms,
    arrays of all second atoms. */

void gen_nblist(struct atomgrp *ag, struct cubeset *cust, struct clusterset *clst,
                     int *excl_list, int **pd1, int **pd2, int ndm, struct nblist *nblst)
{
        int i, j, ic1, ic2, icl1, icl2;
        int k1, k2, ak1, ak2, ka1, ka2, ka3, ip;
        double mdc1, mdc2, dna, dda, dx, dy, dz, d;
        double x1, x2, y1, y2, z1, z2, nbcut;
        double dd, dista;

        int natoms=ag->natoms;

/* prepare nonbond list arrays. */
        int *p;
        for(i=0; i<nblst->nfat; i++)
        {
           if((nblst->nsat[i])>0)
               free(nblst->isat[i]);
        }
        nblst->nfat=0;
        nblst->ifat=_mol_realloc(nblst->ifat,natoms*sizeof(int));
        nblst->nsat=_mol_realloc(nblst->nsat,natoms*sizeof(int));
        nblst->isat=_mol_realloc(nblst->isat,natoms*sizeof(int*));
//ifat is temporal
        int *ifat=_mol_malloc(natoms*sizeof(int));
        for(i=0; i<natoms; i++)
        {
           nblst->nsat[i]=0;
           nblst->ifat[i]=-1;
           ifat[i]=-1;
        }

        nbcut=nblst->nbcut;
        double nbcuts=nbcut*nbcut;

/* loop over the filled cubes. */
	for (i=0; i<cust->nfcubes; i++)
	{
           ic1=cust->ifcubes[i];
           icl1=cust->cubes[ic1].hstincube;
/* loop over all clusters in cube1. */
	   while(icl1>=0)
           {
              mdc1=clst->clusters[icl1].mdc;
              x1=clst->clusters[icl1].gcent[0];
              y1=clst->clusters[icl1].gcent[1];
              z1=clst->clusters[icl1].gcent[2];
/* loop over neighboring cubes including cube1. */
              for(j=-1; j<cust->cubes[ic1].nncubes; j++)
              {
                 if(j==-1)                                  // cube-self part
                    icl2=clst->clusters[icl1].nextincube;
                 else                                       // cube-neighbor cube part
                 {
                    ic2=cust->cubes[ic1].icubes[j];
                    icl2=cust->cubes[ic2].hstincube;
                 }
/* second loop over clusters. */
                 while(icl2>=0)
                 {
                    mdc2=clst->clusters[icl2].mdc;
                    x2=clst->clusters[icl2].gcent[0];
                    y2=clst->clusters[icl2].gcent[1];
                    z2=clst->clusters[icl2].gcent[2];
                    dx=x1-x2;
                    dx*=dx;
                    dy=y1-y2;
                    dy*=dy;
                    dz=z1-z2;
                    dz*=dz;
                    d=dx+dy+dz;
                    dna=nbcut-mdc1-mdc2;
                    dna*=dna;
                    dda=nbcut+mdc1+mdc2;
                    dda*=dda;
                    if(d<=dda)
                    {
/* for clusters icl1 and icl2 go through all pairs of their atoms. */
                       for(k1=0; k1<clst->clusters[icl1].natoms; k1++)
                       {
                          ak1=clst->clusters[icl1].iatom[k1];
                          for(k2=0; k2<clst->clusters[icl2].natoms; k2++)
                          {
                             ak2=clst->clusters[icl2].iatom[k2];
                             if(ag->atoms[ak1].fixed+ag->atoms[ak2].fixed==2)continue;
                             if(ak1<ak2)
                             {
                                ka1=ak1;
                                ka2=ak2;
                             }
                             else
                             {
                                ka1=ak2;
                                ka2=ak1;
                             }
/* check the exclusion table. */
                             ka3=exta(ka1, ka2, excl_list, pd1, pd2, ndm);
                             if(ka3>0)continue;
/* check distances of distant atom pairs (d>dna). */
                             if(d>dna)
                             {
                                dista=0.0;
                                dd=ag->atoms[ka1].X-ag->atoms[ka2].X;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Y-ag->atoms[ka2].Y;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Z-ag->atoms[ka2].Z;
                                dista+=dd*dd;
                                if(dista>nbcuts)continue;
                             }
/* write a successful atom pair to a nonbond list structure. */
                             if(ifat[ka1]==-1)
                             {
                                nblst->ifat[nblst->nfat]=ka1;
                                ifat[ka1]=nblst->nfat;
                                nblst->isat[nblst->nfat]=_mol_malloc((natoms-ka1-1)*sizeof(int));
                                (nblst->nfat)++;
                             }
                             ip=ifat[ka1];
                             p=nblst->isat[ip];
                             p[nblst->nsat[ip]]=ka2;
                             (nblst->nsat[ip])++;
                          }
                       }
                    }
                    icl2=clst->clusters[icl2].nextincube;
                 }
              }
              icl1=clst->clusters[icl1].nextincube;
           }
	}
/* cleanup, total pair calculation. */
        int ntotp=0;
        free(ifat);
        if(nblst->nfat > 0)  {
            nblst->ifat=_mol_realloc(nblst->ifat,(nblst->nfat)*sizeof(int));
            nblst->nsat=_mol_realloc(nblst->nsat,(nblst->nfat)*sizeof(int));
            nblst->isat=_mol_realloc(nblst->isat,(nblst->nfat)*sizeof(int*));
        }

        for(i=0; i<nblst->nfat; i++)
        {
           ntotp+=nblst->nsat[i];
           nblst->isat[i]=_mol_realloc(nblst->isat[i],(nblst->nsat[i])*sizeof(int));
        }
        nblst->npairs=ntotp;
}



void gen_nblist_new(struct atomgrp *ag, struct cubeset *cust, struct clusterset *clst,
                int *excl_list, int **pd1, int **pd2, int ndm, struct nblist *nblst)
{
        int i, j, ic1, ic2, icl1, icl2;
        int k1, k2, ak1, ak2, ka1, ka2, ka3, ip;
        double mdc1, mdc2, dna, dda, dx, dy, dz, d;
        double x1, x2, y1, y2, z1, z2, nbcut;
        double dd, dista;

        int natoms=ag->natoms;

//prepare nonbond list arrays
        int *p;
        for(i=0; i<nblst->nfat; i++)
        {
           if((nblst->nsat[i])>0)
               free(nblst->isat[i]);
        }
        nblst->nfat=0;
        nblst->ifat=_mol_realloc(nblst->ifat,natoms*sizeof(int));
        nblst->nsat=_mol_realloc(nblst->nsat,natoms*sizeof(int));
        nblst->isat=_mol_realloc(nblst->isat,natoms*sizeof(int*));
//ifat is temporal
        int *ifat=_mol_malloc(natoms*sizeof(int));
        for(i=0; i<natoms; i++)
        {
           nblst->nsat[i]=0;
           nblst->ifat[i]=-1;
           ifat[i]=-1;
        }

        nbcut=nblst->nbcut;
        double nbcuts=nbcut*nbcut;

//loop over the filled cubes
	for (i=0; i<cust->nfcubes; i++)
	{
           ic1=cust->ifcubes[i];
           icl1=cust->cubes[ic1].hstincube;
//loop over all clusters in cube1
	   while(icl1>=0)
           {
              mdc1=clst->clusters[icl1].mdc;
              x1=clst->clusters[icl1].gcent[0];
              y1=clst->clusters[icl1].gcent[1];
              z1=clst->clusters[icl1].gcent[2];
//loop over neighboring cubes including cube1
              for(j=-1; j<cust->cubes[ic1].nncubes; j++)
              {
                 if(j==-1)                                  // cube-self part
                    icl2=clst->clusters[icl1].nextincube;
                 else                                       // cube-neighbor cube part
                 {
                    ic2=cust->cubes[ic1].icubes[j];
                    icl2=cust->cubes[ic2].hstincube;
                 }
//second loop over clusters
                 while(icl2>=0)
                 {
                    mdc2=clst->clusters[icl2].mdc;
                    x2=clst->clusters[icl2].gcent[0];
                    y2=clst->clusters[icl2].gcent[1];
                    z2=clst->clusters[icl2].gcent[2];
                    dx=x1-x2;
                    dx*=dx;
                    dy=y1-y2;
                    dy*=dy;
                    dz=z1-z2;
                    dz*=dz;
                    d=dx+dy+dz;
                    dna=nbcut-mdc1-mdc2;
                    dna*=dna;
                    dda=nbcut+mdc1+mdc2;
                    dda*=dda;
                    if(d<=dda)
                    {
//for clusters icl1 and icl2 go through all pairs of their atoms
                       for(k1=0; k1<clst->clusters[icl1].natoms; k1++)
                       {
                          ak1=clst->clusters[icl1].iatom[k1];
                          for(k2=0; k2<clst->clusters[icl2].natoms; k2++)
                          {
                             ak2=clst->clusters[icl2].iatom[k2];
                             if(ag->atoms[ak1].fixed+ag->atoms[ak2].fixed==2)continue;
                             if(ak1<ak2)
                             {
                                ka1=ak1;
                                ka2=ak2;
                             }
                             else
                             {
                                ka1=ak2;
                                ka2=ak1;
                             }
//check the exclusion table
                             ka3=exta(ka1, ka2, excl_list, pd1, pd2, ndm);
                             if(ka3>0)continue;
//check distances of distant atom pairs (d>dna)
                             if(d>dna)
                             {
                                dista=0.0;
                                dd=ag->atoms[ka1].X-ag->atoms[ka2].X;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Y-ag->atoms[ka2].Y;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Z-ag->atoms[ka2].Z;
                                dista+=dd*dd;
                                if(dista>nbcuts)continue;
                             }
//write a successful atom pair to a nonbond list structure
                             if(ifat[ka1]==-1)
                             {
//                                nblst->ifat[nblst->nfat]=ka1;
                                ifat[ka1]=nblst->nfat;
//                                nblst->isat[nblst->nfat]=_mol_malloc((natoms-ka1-1)*sizeof(int));
                                (nblst->nfat)++;
                             }
                             ip=ifat[ka1];
//                             p=nblst->isat[ip];
//                             p[nblst->nsat[ip]]=ka2;
                             (nblst->nsat[ip])++;
                          }
                       }
                    }
                    icl2=clst->clusters[icl2].nextincube;
                 }
              }
              icl1=clst->clusters[icl1].nextincube;
           }
	}

        for(i=0; i<natoms; i++)
        {
           ip=ifat[i];
           if ( ip > -1 )
             {
               if ( nblst->nsat[ip] > 0 ) nblst->isat[ip]=_mol_malloc(nblst->nsat[ip]*sizeof(int));
               nblst->nsat[ip]=0;
             }  
           nblst->ifat[i]=-1;           
           ifat[i]=-1;
        }
        nblst->nfat=0;	
	
//loop over the filled cubes
	for (i=0; i<cust->nfcubes; i++)
	{
           ic1=cust->ifcubes[i];
           icl1=cust->cubes[ic1].hstincube;
//loop over all clusters in cube1
	   while(icl1>=0)
           {
              mdc1=clst->clusters[icl1].mdc;
              x1=clst->clusters[icl1].gcent[0];
              y1=clst->clusters[icl1].gcent[1];
              z1=clst->clusters[icl1].gcent[2];
//loop over neighboring cubes including cube1
              for(j=-1; j<cust->cubes[ic1].nncubes; j++)
              {
                 if(j==-1)                                  // cube-self part
                    icl2=clst->clusters[icl1].nextincube;
                 else                                       // cube-neighbor cube part
                 {
                    ic2=cust->cubes[ic1].icubes[j];
                    icl2=cust->cubes[ic2].hstincube;
                 }
//second loop over clusters
                 while(icl2>=0)
                 {
                    mdc2=clst->clusters[icl2].mdc;
                    x2=clst->clusters[icl2].gcent[0];
                    y2=clst->clusters[icl2].gcent[1];
                    z2=clst->clusters[icl2].gcent[2];
                    dx=x1-x2;
                    dx*=dx;
                    dy=y1-y2;
                    dy*=dy;
                    dz=z1-z2;
                    dz*=dz;
                    d=dx+dy+dz;
                    dna=nbcut-mdc1-mdc2;
                    dna*=dna;
                    dda=nbcut+mdc1+mdc2;
                    dda*=dda;
                    if(d<=dda)
                    {
//for clusters icl1 and icl2 go through all pairs of their atoms
                       for(k1=0; k1<clst->clusters[icl1].natoms; k1++)
                       {
                          ak1=clst->clusters[icl1].iatom[k1];
                          for(k2=0; k2<clst->clusters[icl2].natoms; k2++)
                          {
                             ak2=clst->clusters[icl2].iatom[k2];
                             if(ag->atoms[ak1].fixed+ag->atoms[ak2].fixed==2)continue;
                             if(ak1<ak2)
                             {
                                ka1=ak1;
                                ka2=ak2;
                             }
                             else
                             {
                                ka1=ak2;
                                ka2=ak1;
                             }
//check the exclusion table
                             ka3=exta(ka1, ka2, excl_list, pd1, pd2, ndm);
                             if(ka3>0)continue;
//check distances of distant atom pairs (d>dna)
                             if(d>dna)
                             {
                                dista=0.0;
                                dd=ag->atoms[ka1].X-ag->atoms[ka2].X;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Y-ag->atoms[ka2].Y;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Z-ag->atoms[ka2].Z;
                                dista+=dd*dd;
                                if(dista>nbcuts)continue;
                             }
//write a successful atom pair to a nonbond list structure
                             if(ifat[ka1]==-1)
                             {
                                nblst->ifat[nblst->nfat]=ka1;
                                ifat[ka1]=nblst->nfat;
//                                nblst->isat[nblst->nfat]=_mol_malloc((natoms-ka1-1)*sizeof(int));
                                (nblst->nfat)++;
                             }
                             ip=ifat[ka1];
                             p=nblst->isat[ip];
                             p[nblst->nsat[ip]]=ka2;
                             (nblst->nsat[ip])++;
                          }
                       }
                    }
                    icl2=clst->clusters[icl2].nextincube;
                 }
              }
              icl1=clst->clusters[icl1].nextincube;
           }
	}	
//cleanup, total pair calculation
        int ntotp=0;
        free(ifat);
        if(nblst->nfat > 0)  {
            nblst->ifat=_mol_realloc(nblst->ifat,(nblst->nfat)*sizeof(int));
            nblst->nsat=_mol_realloc(nblst->nsat,(nblst->nfat)*sizeof(int));
            nblst->isat=_mol_realloc(nblst->isat,(nblst->nfat)*sizeof(int*));
        }

        for(i=0; i<nblst->nfat; i++)
        {
           ntotp+=nblst->nsat[i];
//           nblst->isat[i]=_mol_realloc(nblst->isat[i],(nblst->nsat[i])*sizeof(int));
        }
        nblst->npairs=ntotp;
}


void free_cubeset(struct cubeset *cust)
{
        int i,ic;

        for (i=0; i<(cust->nfcubes); i++)
        {
             ic=cust->ifcubes[i];
             if (cust->cubes[ic].nncubes>0)
                   free(cust->cubes[ic].icubes);
        }
        free(cust->cubes);
        free(cust->ifcubes);
}

//! Generate a set of cubes cust based on the clusterset clst.
/*! Called in update_nblst. */
        /*! 0. Sum of nblist cutoff and margin size gives the cubesize.
            1. Calculate coordinate span of cluster centers
            2. Number of cubes in each direction is span over cubesize plus 1
            3. For each cluster find a cube it belongs to.
            4. Each cube knows last cluster contained in it.
            5. Each cluster knows another cluster belongin to the same cube (linked list)
            6. Each cube knows its occupied neighbouring cubes. */

void gen_cubeset(double nbcut, struct clusterset *clst, struct cubeset *cust)
{
	double xmin, xmax, ymin, ymax, zmin, zmax, cubel;
        int i, nx, ny, nz, nxy, n, j, k, l;
        int    ix, iy, iz, ic, ic1;
        int    ix1, ix2, iy1, iy2, iz1, iz2;
        int    nnbr, temp[27];
        /* sum of nblist cutoff and margin size gives the cubesize cubel. */
        cubel=nbcut+clst->marg;
        cust->cubel=cubel;
        /* calculate coordinate span of cluster centers. */
        xmin=clst->clusters[0].gcent[0];
        xmax=xmin;
        ymin=clst->clusters[0].gcent[1];
        ymax=ymin;
        zmin=clst->clusters[0].gcent[2];
        zmax=zmin;
        for(i=1; i<clst->nclusters; i++)
	{
 	   if(xmin>clst->clusters[i].gcent[0])xmin=clst->clusters[i].gcent[0];
           if(xmax<clst->clusters[i].gcent[0])xmax=clst->clusters[i].gcent[0];
           if(ymin>clst->clusters[i].gcent[1])ymin=clst->clusters[i].gcent[1];
           if(ymax<clst->clusters[i].gcent[1])ymax=clst->clusters[i].gcent[1];
           if(zmin>clst->clusters[i].gcent[2])zmin=clst->clusters[i].gcent[2];
           if(zmax<clst->clusters[i].gcent[2])zmax=clst->clusters[i].gcent[2];
	}
        /* number of cubes in each direction is span over cubesize plus 1. */
        nx=1+(int)((xmax-xmin)/cubel);
        ny=1+(int)((ymax-ymin)/cubel);
        nz=1+(int)((zmax-zmin)/cubel);
        nxy=nx*ny;
        n=nxy*nz;
        cust->ncubes=n;
        cust->nfcubes=0;
        cust->cubes=_mol_malloc(n*sizeof(struct cube));
        cust->ifcubes=_mol_malloc(n*sizeof(int));
        for(i=0; i<n; i++)
          {
	   cust->cubes[i].hstincube=-1;
	   cust->cubes[i].nncubes=0;
	  }
        /* for each cluster find a cube it belongs. */
        for(i=0; i<clst->nclusters; i++)
        {
	   ix=(int)((clst->clusters[i].gcent[0]-xmin)/cubel);
           iy=(int)((clst->clusters[i].gcent[1]-ymin)/cubel);
           iz=(int)((clst->clusters[i].gcent[2]-zmin)/cubel);
           ic=ix+iy*nx+iz*nxy;
           if(cust->cubes[ic].hstincube==-1)
           {
              cust->cubes[ic].ix=ix;
              cust->cubes[ic].iy=iy;
              cust->cubes[ic].iz=iz;
              cust->ifcubes[(cust->nfcubes)++]=ic;
           }
           /* each cluster knows another cluster belonging to the same cube (linked list). */
	   clst->clusters[i].nextincube=cust->cubes[ic].hstincube;
           /* each cube knows last cluster contained in it. */ 
	   cust->cubes[ic].hstincube=i;
        }
        /* loop over occupied cubes. */
        for(i=0; i<cust->nfcubes; i++)
        {
	   ic=cust->ifcubes[i];
           if(cust->cubes[ic].hstincube==-1)continue;
           ix=cust->cubes[ic].ix;
              ix1=ix-1;
              if(ix1<0)ix1=0;
              ix2=ix+1;
              if(ix2>nx-1)ix2=nx-1;
           iy=cust->cubes[ic].iy;
              iy1=iy-1;
              if(iy1<0)iy1=0;
              iy2=iy+1;
              if(iy2>ny-1)iy2=ny-1;
           iz=cust->cubes[ic].iz;
              iz1=iz-1;
              if(iz1<0)iz1=0;
              iz2=iz+1;
              if(iz2>nz-1)iz2=nz-1;
           nnbr=0;
           for(j=iz1; j<=iz2; j++)
           {
              for(k=iy1; k<=iy2; k++)
              {
                 for(l=ix1; l<=ix2; l++)
                 {
                    ic1=l+k*nx+j*nxy;
                    if(ic1<=ic)continue;
                    if(cust->cubes[ic1].hstincube==-1)continue;
                    temp[nnbr++]=ic1;
                 }
              }
           }
           cust->cubes[ic].nncubes=nnbr;
           /* each cube knows its occupied neighbouring cubes. */
           if(nnbr>0)
           {
              cust->cubes[ic].icubes=_mol_malloc(nnbr*sizeof(int));
              for(j=0; j<nnbr; j++)cust->cubes[ic].icubes[j]=temp[j];
           }
	}
}

//! Access margin size of future cubes.
/*! Calculate geometrical center of each cluster and margine size which is a sum of 
    2 longest distances from
    the geometrical center of any cluster to one of its member atoms. */
void findmarg(struct atomgrp *ag, struct clusterset *clst)
{
	int i, j, k, a;
        double dist1=0, dist2=0;
        double dx, dy, dz, d, md, gc[3];
        for(i=0; i<clst->nclusters; i++)
	{
           md=0.0;
           for(j=0; j<3; j++)gc[j]=0.0;
           a=clst->clusters[i].natoms;
           for(j=0; j<a; j++)
           {
              k=clst->clusters[i].iatom[j];
              gc[0]+=(ag->atoms[k].X)/a;  /* geometrical center of a cluster */
              gc[1]+=(ag->atoms[k].Y)/a;
              gc[2]+=(ag->atoms[k].Z)/a;
           }
           for(j=0; j<a; j++)
           {
              k=clst->clusters[i].iatom[j];
              dx=ag->atoms[k].X-gc[0];
              dx*=dx;
              dy=ag->atoms[k].Y-gc[1];
              dy*=dy;
              dz=ag->atoms[k].Z-gc[2];
              dz*=dz;
              d=sqrt(dx+dy+dz);
              if(d>md)md=d;    /* maximum distance for a given cluster */
              if(d>dist1)
              {
                 dist2=dist1; /* dist1 and dist2 are two longest distances among the set of clusters. */
                 dist1=d;
              }
              else if(d>dist2)dist2=d;
           }
           clst->clusters[i].mdc=md; /* save maximum distance. */
           for(j=0; j<3; j++)clst->clusters[i].gcent[j]=gc[j]; /*! save cluster center. */
	}
        clst->marg=dist1+dist2; /* save a margine distance for a cluster set. */
}

void free_clset(struct clusterset *clst)
{
        int i;
        for(i=0; i<clst->nclusters; i++)
            free(clst->clusters[i].iatom);
        free(clst->clusters);
}

void gen_clset(int natoms, struct clusterset *clst,
               int nclust, int *clust)
{
	clst->nclusters=nclust;
        clst->clusters=_mol_malloc(nclust*sizeof(struct cluster));
        int i, j;
        for(i=0; i<nclust; i++)
            clst->clusters[i].natoms=0;
        for(i=0; i<natoms; i++)
	    (clst->clusters[clust[i]].natoms)++;
        for(i=0; i<nclust; i++)
	{
	    clst->clusters[i].iatom=
            _mol_malloc((clst->clusters[i].natoms)*sizeof(int));
            clst->clusters[i].natoms=0;
	}
        for(i=0; i<natoms; i++)
        {
            j=clst->clusters[clust[i]].natoms;
            clst->clusters[clust[i]].iatom[j++]=i;
            clst->clusters[clust[i]].natoms=j;
        }
}

//! Assign atoms to clusters based on three bonds connectivity.
/*! 1-4 clustering was performed the following way.
            1. All atoms were unassigned (flag -1) to any cluster.
            2. Loop through all atoms
            3. If atom is unassigned find all anassigned atoms connected to it through
               1-3 connection and save a list of all such connections (bonds)
            4. loop through the list of saved connections and for each calculate the
               number of unassigned atoms for which this connection is 2-3 connection.
               Connection with the maximum such number is the 2-3 connection of
               a new cluster.
            5. No restriction on carbon being 1 or 4 atom in the cluster was made
               (differ from the original published algorithm). */
void clust14(struct atomgrp* ag, int *nclust, int *clust)
{

	*nclust=0;
	int npair;
	int ma, mapair, ipair;
	int *pairs=_mol_malloc(MAXPAIR*2*sizeof(int));
	int i, j, k, l, m, n=ag->natoms;
        struct atombond *b;
	struct atom *a0, *a1, *a2;
        /*All atoms were unassigned (flag -1) to any cluster. */
	for(i=0; i<n; i++)
	{
		clust[i]=-1;
	}
        /*Loop through all atoms. */
        for(i=0; i<n; i++)
        {
	   if(clust[i]==-1)
	   {
		npair=0;
		a0=&(ag->atoms[i]);
                /* If atom is unassigned find all anassigned atoms connected to it through
                    1-3 connection and save a list of all such connections (bonds) in an array
                    pairs. */
		for(j=0; j<(a0->nbonds); j++)
		{
	  	   b=a0->bonds[j];
		   a1=b->a1;
		   if(a1==a0)a1=b->a0;
		   k=(a1->ingrp);
		   if(clust[k]==-1)
		   {
			addpair(i, k, &npair, pairs);
			for(l=0; l<(a1->nbonds); l++)
			{
			   b=a1->bonds[l];
			   a2=b->a1;
			   if(a2==a1)a2=b->a0;
			   if(a2==a0)continue;
                           m=(a2->ingrp);
                           if(clust[m]==-1)
				addpair(k, m, &npair, pairs);
			}
		   }
		}
		if(npair==0)clust[i]=*nclust;
		else
		{
		   mapair=2;
		   ipair=0;
                   /* loop through the list of saved connections and for each calculate the
                       number of unassigned atoms for which this connection is 2-3 connection.
                       Connection with the maximum such number is the 2-3 connection of
                       a new cluster. */
		   for(l=0; l<npair; l++)
		   {
                        /* return the number of atoms in 1-4 cluster for a given 2-3 connection. */
			ma=natpair(ag, &pairs[2*l], clust);
			if(ma>mapair)
			{
			   ipair=l;
			   mapair=ma;
			}
		   }
                   /* Assign cluster number to the atom indices in the clust array. */
		   addclust(ag, &pairs[2*ipair], clust, *nclust);
		}
		(*nclust)++;
	   }
        }
	free(pairs);
}

//! Increase the list of atom pairs
/*! used in clust14 function. */ 
void addpair(int i, int k, int *npair, int* pairs)
{
   int l=*npair;
   if(l+1>MAXPAIR)
   {
       printf("Increase MAXPAIR in clust14\n");
       exit (EXIT_FAILURE);
   }
   pairs[2*l]=i;
   pairs[2*l+1]=k;
   (*npair)+=1;
}


//! Return the number of atoms in 1-4 cluster for a given 2-3 connection.
/* used in clust14 function. */
int natpair(struct atomgrp* ag, int *pair, int *clust)
{
	int n=0, j, k;
	struct atom *a1=&(ag->atoms[pair[0]]);
	struct atom *a2=&(ag->atoms[pair[1]]);
	struct atom *a;
	struct atombond *b;

        for(j=0; j<(a1->nbonds); j++)
        {
           b=a1->bonds[j];
           a=b->a1;
           if(a==a1)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)n++;
	}
        for(j=0; j<(a2->nbonds); j++)
        {
           b=a2->bonds[j];
           a=b->a1;
           if(a==a2)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)n++;
        }
	return n--;
}

//! Assign cluster number to the atom indices in the clust array.
/* used in clust14 function. */

void addclust(struct atomgrp* ag, int *pair, int *clust, int nclust)
{
        int j, k;
        struct atom *a1=&(ag->atoms[pair[0]]);
        struct atom *a2=&(ag->atoms[pair[1]]);
        struct atom *a;
        struct atombond *b;

        for(j=0; j<(a1->nbonds); j++)
        {
           b=a1->bonds[j];
           a=b->a1;
           if(a==a1)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)clust[k]=nclust;
        }
        for(j=0; j<(a2->nbonds); j++)
        {
           b=a2->bonds[j];
           a=b->a1;
           if(a==a2)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)clust[k]=nclust;
        }
}

//! Calculate a one dimensional bonded list list01 for the atomgroup ag.
/*! na01[i] - number of bonded atoms for the atom i
    list01 - sequential blocks of numbers of bonded atoms
    pna01[i] - pointer to the block of bonded atoms in list01 for the atom i */
void comp_list01(struct atomgrp* ag, int *list01, int *na01, int **pna01)
{
	int i, j, k, nb;
	struct atom *a, *a0;
	struct atombond *b;
	for(i=0; i<ag->natoms; i++)
	{
		na01[i]=0;
		a0=&(ag->atoms[i]);
		nb=a0->nbonds;
		pna01[i]=list01;
		for(j=0; j<nb; j++)
		{
			b=a0->bonds[j];
			a=b->a1;
			if(a==a0)a=b->a0;
			k=(a->ingrp);
			*list01=k;
			na01[i]++;
			list01++;
		}
		if(na01[i]==0)pna01[i]=NULL;
	}
}

//! Calculate a total number of bonded triplets of atoms n02 and quadruplets n03.
/*! Based on number of connections na01 and an array of pointers pna01 to the list01. */
void comp_n23(int natoms, int *na01, int **pna01, int *n02, int *n03)
{
	int i, j, k, l, m, *p;
	*n02=0;
	*n03=0;
	for(i=0; i<natoms; i++)
        {
	   j=na01[i];
	   if(j>1)
	   {
		(*n02)+=j*(j-1)/2;
		p=pna01[i];
		for(k=0; k<j; k++)
		{
		   l=p[k];
		   if(l>i)
		   {
			m=na01[l];
			if(m>1)(*n03)+=(j-1)*(m-1);
		   }
		}
	   }
	}
}

int trim_comp(const void *s1, const void *s2)
{
     int *v1, *v2;
     v1=(int *)s1;
     v2=(int *)s2;

     int i11=*v1;
     int i21=*v2;
     if(i11<i21)
        return -1;
     if(i11>i21)
        return 1;

     int i12=*(v1+1);
     int i22=*(v2+1);
     if(i12<i22)
        return -1;
     if(i12>i22)
        return 1;
     return 0;
}

//! Calculate a list of atoms connected by two bonds list02.
/*! Based on number of connections na01 and an array of pointers pna01 to the list01. 
    na02[i] number of atom pairs connected by two bonds with the smallest atom index i
    na02[i]=0 for some i
    pna02[i] pointer to the block in list02 associated with the smallest atom index i
    pna02[i]=NULL for some i.*/
void comp_list02(int natoms, int *na01, int **pna01,
                 int *na02, int **pna02, int *n02, int *list02)
{
	int i,j,k,l,m,n, *p;
        int *list=list02;

	for(i=0; i<natoms; i++)
        {
           j=na01[i];
           if(j>1)
           {
		p=pna01[i];
		for(k=0; k<j-1; k++)
		{
			n=p[k];
			for(l=k+1; l<j; l++)
			{
			   m=p[l];
			   if(m>n)
			   {
				*list02=n;
				list02++;
				*list02=m;
				list02++;
			   }
			   else
                           {
                                *list02=m;
                                list02++;
                                *list02=n;
                                list02++;
                           }
			}
		}
	   }
	}

        for(i=0; i<*n02; i++)
        {
           j=list[2*i];
           k=list[2*i+1];
           l=na01[j];
           if(l>0)
           {
              p=pna01[j];
              for(m=0; m<l; m++)
              {
                 n=p[m];
                 if(n==k)
                 {
                    list[2*i]=natoms;
                    list[2*i+1]=natoms;
                    break;
                 }
              }
           }
        }
        qsort(list, *n02, 2*sizeof(list[0]), trim_comp);
        j=(*n02)-1;
        for(i=j; i>0; i--)
        {
           if(list[2*i]==natoms)
               (*n02)--;
           else
               break;
        }
        for(i=0; i<natoms; i++)
        {
                na02[i]=0;
                pna02[i]=NULL;
        }
        for(i=0; i<*n02; i++)
        {
                j=list[2*i];
                (na02[j])++;
                if(i==0||j>list[2*i-2])
                   pna02[j]=&(list[2*i]);
        }
}

//! Calculate list of atom pairs connected by three bonds list03.
/*! It is redundand so far. */
void comp_list03(int natoms, int *na01, int **pna01, int *list03)
{
	int i,j,k,l,m,n,q,r,s, *p, *o;
	for(i=0; i<natoms; i++)
	{
	   j=na01[i];
	   if(j>1)
	   {
		p=pna01[i];
		for(k=0; k<j; k++)
                {
		   l=p[k];
		   if(l>i)
		   {
			m=na01[l];
			if(m>1)
			{
			   o=pna01[l];
			   for(n=0; n<m; n++)
			   {
			      q=o[n];
			      if(q==i)continue;
			      for(r=0; r<j; r++)
			      {
			         if(r==k)continue;
			         s=p[r];
	                         if(q>s)
                           	 {
                                    *list03=s;
                                    list03++;
                                    *list03=q;
                                    list03++;
                           	 }
                                 else
                                 {
                                    *list03=q;
                                    list03++;
                                    *list03=s;
                                    list03++;
                                 }
			      }
			   }
			}
		   }
		}
	   }
	}
}

//! Create list listf03 from list03 by removing pairs of "fixed" atoms
/*! nf03 is a number of nonfixed pairs */
void fix_list03(struct atomgrp *ag, int n03, int* list03, int *nf03, int *listf03)
{
    int i, i1, i2;
    int n=0;
    struct atom *a1, *a2;
    for(i=0; i<n03; i++)
    {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           if(a2->fixed == 1 && a1->fixed == 1)continue;
           listf03[2*n]=i1;
           listf03[2*n+1]=i2;
           n++;
    }
    *nf03=n;
}

//! Remove redundancies from list03.
/*! Reset the number n03 */
void trim_list03(int natoms,
                 int *na01, int **pna01,
                 int *na02, int **pna02,
                 int *n03, int *list03)
{
	qsort(list03, *n03, 2*sizeof(list03[0]), trim_comp);
        int i, j, k, l, m, n, o, *p;
        int *list=_mol_malloc(*n03*sizeof(int));
        for(i=1; i<*n03; i++)
        {
           j=2*i;
           list[i]=0;
           k=list03[j];
           l=list03[j+1];
// redundant pairs
           if(k==list03[j-2] && l==list03[j-1])
           {
              list[i]=1;
              continue;
           }
// directly bonded pairs
           m=na01[k];
           if(m>0)
           {
              p=pna01[k];
              for(n=0; n<m; n++)
              {
                 o=p[n];
                 if(l==o)
                 {
                    list[i]=1;
                    continue;
                 }
              }
           }
// bonded through an angle pairs
           m=na02[k];
           if(m>0)
           {
              p=pna02[k];
              for(n=0; n<m; n++)
              {
                 o=p[2*n+1];
                 if(l==o)
                 {
                    list[i]=1;
                    continue;
                 }
              }
           }
        }
        for(i=1; i<*n03; i++)
        {
           if(list[i]==1)
           {
               j=2*i;
               list03[j]=natoms;
               list03[j+1]=natoms;
           }
        }
        qsort(list03, *n03, 2*sizeof(list03[0]), trim_comp);
        j=*n03-1;
        for(i=j; i>0; i--)
        {
           if(list03[2*i]==natoms)(*n03)--;
           else break;
        }
        free(list);
//        for(i=0; i<*n03; i++)
//           printf("list03 %d, %d, %d, %d\n", i, *n03, list03[2*i], list03[2*i+1]);
//        exit(0);
}

//! Evaluate dimensions of the exclusion list.
/*! The dimensions of the exclusion
    list relate to how close the connected atoms are
    in terms of their numbers. */
void excl_dims(int natoms, int *na01, int **pna01,
               int n02, int *list02, int n03, int *list03,
               int *nd1, int *nd2, int *ndm, int *atmind)
{
	int i,  j,  k, l, *p;
	for(i=0; i<2*natoms; i++)atmind[i]=0;
	*ndm=20;
	*nd2=0;
	*nd1=0;
	for(i=0; i<natoms; i++)
	{
		p=pna01[i];
		for(k=0; k<na01[i]; k++)
		{
		   j=p[k];
		   l=j-i;
		   if(l>*ndm)*ndm=l;
		   l--;
		   if(l>19 && atmind[2*i+1]==0)
		   {
			atmind[2*i+1]=1;
			(*nd2)++;
		   }
		   else if(l>9 && atmind[2*i]==0)
		   {
			atmind[2*i]=1;
			(*nd1)++;
		   }
		}
	}
	for(i=0; i<n02; i++)
	{
		j=list02[2*i];
		k=list02[2*i+1];
		l=k-j;
		if(l>*ndm)*ndm=l;
		l--;
                if(l>19 && atmind[2*j+1]==0)
                {
                     atmind[2*j+1]=1;
                     (*nd2)++;
                }
                else if(l>9 && atmind[2*j]==0)
                {
                     atmind[2*j]=1;
                     (*nd1)++;
                }
	}
        for(i=0; i<n03; i++)
        {
                j=list03[2*i];
                k=list03[2*i+1];
                l=k-j;
                if(l>*ndm)
                {
                   *ndm=l;
//                   printf("ndm %d %d %d\n",j,k,*ndm);
                }
                l--;
                if(l>19 && atmind[2*j+1]==0)
                {
                     atmind[2*j+1]=1;
                     (*nd2)++;
                }
                else if(l>9 && atmind[2*j]==0)
                {
                     atmind[2*j]=1;
                     (*nd1)++;
                }
        }
}

//! Evaluate if a pair of atoms should be excluded from nonbonded list.
/*! Based on connectivity data in exclusion list excl_list/pd1/pd2. */
inline int exta(int a1, int a2, int *excl_list, int **pd1, int **pd2, int ndm)
{
	int i=a2-a1-1;
        if(i<0)
        {
           printf("exta: ERROR a1<=a2\n");
           exit (EXIT_FAILURE);
        }
        if(i>=ndm)
           return 0;
        if(i<10)
           return excl_list[10*a1+i];
        int *p;
        if(i<20)
        {
           p=pd1[a1];
           return p[i-10];
        }
        p=pd2[a1];
        return p[i-20];
}

/*
#ifdef _BGL_

// assuming ndm <= 20
int exta2( int a1, int a2, int *excl_list, int **pd1, int ndm )
{
    int i = a2 - a1 - 1;

    if ( i >= ndm ) return 0;

    if ( i < 10 ) return excl_list[ ( a1 << 3 ) + ( a1 << 1 ) + i ];
    else return pd1[ a1 ][ i - 10 ];
}

#endif
*/

//! Create an exclusion list (excl_list).
/*! Exclusion list is a compact array of numbers coding connected through 1,
    2, or 3 bonds pairs of atoms. It is used later by the function exta() to
    quickly evaluate if a pair is excluded. The compactness of the exclusion
    list relates to the convention, that connected atoms are generally close
    to each other by the numbering. */
void excl_tab(int natoms, int *na01, int **pna01,
               int n02, int *list02, int n03, int *list03,
               int nd1, int nd2, int ndm, int *atmind,
               int *excl_list, int **pd1, int **pd2)
{
        int i,  j,  k, l, *p, *p1;
        int md1=-1, md2=-1;
// zero the exclusion list
        i=(natoms+nd1+1)*10+(nd2+1)*(ndm-20);
        for(j=0; j<i; j++)excl_list[j]=0;
// populate atom pointers to the exclusion list
        for(i=0; i<natoms; i++)
        {
           if(atmind[2*i]==1)
              pd1[i]=&(excl_list[10*(++md1+natoms)]);
           else
              pd1[i]=&(excl_list[10*(nd1+natoms)]);
           if(atmind[2*i+1]==1)
              pd2[i]=&(excl_list[10*(natoms+nd1+1)+(ndm-20)*(++md2)]);
           else
              pd2[i]=&(excl_list[10*(natoms+nd1+1)+(ndm-20)*nd2]);
        }
        for(i=0; i<natoms; i++)
        {
                p=pna01[i];
                for(k=0; k<na01[i]; k++)
                {
                   j=p[k];
                   l=j-i-1;
                   if(l>19)
                   {
                      p1=pd2[i];
                      p1[l-20]=1;
                   }
                   else if(l>9)
                   {
                      p1=pd1[i];
                      p1[l-10]=1;
                   }
                   else if(l>=0)
                   {
                      excl_list[10*i+l]=1;
                   }
                }
        }
        for(i=0; i<n02; i++)
        {
                j=list02[2*i];
                k=list02[2*i+1];
                l=k-j-1;
                if(l>19)
                {
                   p1=pd2[j];
                   p1[l-20]=2;
                }
                else if(l>9)
                {
                   p1=pd1[j];
                   p1[l-10]=2;
                }
                else if(l>=0)
                {
                   excl_list[10*j+l]=2;
                }
        }
        for(i=0; i<n03; i++)
        {
                j=list03[2*i];
                k=list03[2*i+1];
                l=k-j-1;
                if(l>19)
                {
                   p1=pd2[j];
                   p1[l-20]=3;
                }
                else if(l>9)
                {
                   p1=pd1[j];
                   p1[l-10]=3;
                }
                else if(l>=0)
                {
                   excl_list[10*j+l]=3;
                }
        }
}

//free nblist arrays in agsetup
void destroy_agsetup(struct agsetup* ags)
{
    free(ags->excl_list);
    free(ags->pd1);
    free(ags->pd2);
    free(ags->listf03);
    free(ags->list03);
    free(ags->list02);
    free_nblist(ags->nblst);
    free_clset(ags->clst);
    free(ags->clst);
}
void free_agsetup(struct agsetup* ags)
{
    destroy_agsetup(ags);
    free(ags);
}


//Setup nblist calculation
//! Initialize a nonbonded list data structure ags for the atomgroup ag.
/*! Nonbonded list generation involves initialization step performed once based on
    the atom connectivity in this function and also update step covered by update_nblst 
    function using spatial positions of atoms. The algorithm is described at:
    Petrella et all, "An Improved Method For Nonbonded List Generation: Rapid Determination
    of Near-Neighbour Pairs" J. Comp. Chem. 2002, v.24 pp. 222-231.
    The algorithm involves assignement of atoms to clusters connected by at most 3 bonds and 
    distribution of the clusters between cubes of nonbonded cutoff length with some padding.
    The speedup of the list generation is achieved by searching only neighbouring cubes and
    by using mostly cluster-cluster distances instead of atom-atom ones.
    */
void init_nblst(struct atomgrp* ag, struct agsetup* ags)
{
        int *list01=_mol_malloc(2*(ag->nbonds)*sizeof(int)); 
        int *na01=_mol_malloc((ag->natoms)*sizeof(int));
	int **pna01=_mol_malloc((ag->natoms)*sizeof(int*));
        int *na02=_mol_malloc((ag->natoms)*sizeof(int));
        int **pna02=_mol_malloc((ag->natoms)*sizeof(int*));

        // Calculate a one dimensional bonded list list01 for the atomgroup ag.
        /* na01[i] - number of bonded atoms for the atom i
            list01 - sequential blocks of numbers of bonded atoms
            pna01[i] - pointer to the block of bonded atoms in list01 for the atom i. */
        comp_list01(ag, list01, na01, pna01);
        int i,j;

        int n02, n03,nf03;
        // Calculate a total number of bonded triplets of atoms n02 and quadruplets n03.
        /* Based on number of connections na01 and an array of pointers pna01 to the list01. */
	comp_n23(ag->natoms, na01, pna01, &n02, &n03);
        int *list02=_mol_malloc(2*n02*sizeof(int));
        int *list03=_mol_malloc(2*n03*sizeof(int));
	int *listf03=_mol_malloc(2*n03*sizeof(int));
        // Calculate a list of atoms connected by two bonds list02.
        /* Based on number of connections na01 and an array of pointers pna01 to the list01.
            na02[i] number of atom pairs connected by two bonds with the smallest atom index i
            na02[i]=0 for some i
            pna02[i] pointer to the block in list02 associated with the smallest atom index i
            pna02[i]=NULL for some i.*/
	comp_list02(ag->natoms, na01, pna01, na02, pna02, &n02, list02);
        // Calculate list of atom pairs connected by three bonds list03.
        /* It is redundand so far. */
        comp_list03(ag->natoms, na01, pna01, list03);
        //! Remove redundancies from list03.
        /*! Reset the number n03. */
        trim_list03(ag->natoms, na01, pna01, na02, pna02, &n03, list03);
        // Create list listf03 from list03 by removing pairs of "fixed" atoms
        /* nf03 is a number of nonfixed pairs. */
	fix_list03(ag,n03,list03,&nf03,listf03);
	if(nf03>0)
           listf03=_mol_realloc(listf03,2*nf03*sizeof(int));
        else
           listf03=_mol_realloc(listf03,2*sizeof(int));
        int nd1,nd2,ndm;
	int *atmind=_mol_malloc(2*(ag->natoms)*sizeof(int));
        // Evaluate dimensions of the exclusion list.
        /* The dimensions of the exclusion
            list relate to how close the connected atoms are
            in terms of their numbers. */
	excl_dims(ag->natoms, na01, pna01,
                       n02, list02, n03, list03,
                       &nd1, &nd2, &ndm, atmind);

	i=10*((ag->natoms)+nd1+1)+(nd2+1)*(ndm-20);
	int *excl_list=_mol_malloc(i*sizeof(int));
        int **pd1=(int**)_mol_malloc((ag->natoms)*sizeof(int*));
        int **pd2=(int**)_mol_malloc((ag->natoms)*sizeof(int*));
        // Create an exclusion list (excl_list).
        /* Exclusion list is a compact array of numbers coding connected through 1,
            2, or 3 bonds pairs of atoms. It is used later by the function exta() to
            quickly evaluate if a pair is excluded. The compactness of the exclusion
            list relates to the convention, that connected atoms are generally close
            to each other by the numbering. */
        excl_tab(ag->natoms, na01, pna01,
               n02, list02, n03, list03,
               nd1, nd2, ndm, atmind,
               excl_list, pd1, pd2);
        // Pack connectivity data into nonbonded list ags structure.
        /* Assemble all lists pointers and dimensions. */ 
	ags->list02=list02;
	ags->list03=list03;
	ags->listf03=listf03;
	ags->n02=n02;
	ags->n03=n03;
	ags->nf03=nf03;
	ags->ndm=ndm;
	ags->excl_list=excl_list;
	ags->pd1=pd1;
	ags->pd2=pd2;
        
        // Free temporal storage and data.
        /* 01 and 02 lists and temporal atmind go. */
	free(pna01);
        free(na01);
        free(list01);

        free(pna02);
        free(na02);
        free(atmind);
//***************Generating cluster------------------------
	int *clust=_mol_malloc((ag->natoms)*sizeof(int));
        int nclust;
        // Assign atoms to clusters based on three bonds connectivity.
        /* 1-4 clustering was performed the following way.
            1. All atoms were unassigned (flag -1) to any cluster.
            2. Loop through all atoms
            3. If atom is unassigned find all anassigned atoms connected to it through
               1-3 connection and save a list of all such connections (bonds)
            4. loop through the list of saved connections and for each calculate the 
               number of unassigned atoms for which this connection is 2-3 connection.
               Connection with the maximum such number is the 2-3 connection of
               a new cluster.
            5. No restriction on carbon being 1 or 4 atom in the cluster was made 
               (differ from the original published algorithm). */

        
        clust14(ag, &nclust, clust);

        struct clusterset *clst= _mol_malloc(sizeof(struct clusterset));

        // Convert an array clust to a structure clusterset clst.
        gen_clset(ag->natoms, clst, nclust, clust);
	ags->clst=clst;
//**********************free****************************************
        free(clust);

        struct nblist *nblst=_mol_malloc(sizeof(struct nblist));
        double nbcut=13.0;
        double nbcof=12.0;
        nblst->nbcut=nbcut;
        nblst->nbcof=nbcof;
        nblst->crds=_mol_malloc(3*(ag->natoms)*sizeof(float));
        nblst->nfat=0;
        nblst->ifat=_mol_malloc((ag->natoms)*sizeof(int));
        nblst->nsat=_mol_malloc((ag->natoms)*sizeof(int));
        for(j=0; j<(ag->natoms); j++)
                            nblst->nsat[j]=0;
        nblst->isat=_mol_malloc((ag->natoms)*sizeof(int*));
	ags->nblst=nblst;
}


//! Spatial part of the nblist generation altorithm
/*! See comments to init_nblst function. This part is performed
    periodically when atom coordinates change sufficiently. 
    The steps of the function are (1) Save snapshot of coordinates at the update stage.
    (2) Calculate the margin size of future cubes. (3) Generate a set of cubes cust based 
    on the clusterset ags->clst (4) Generate a nonbonded list nblst from cluster and cube 
    subdivisions and exclusion list arrays. */
void update_nblst(struct atomgrp* ag, struct agsetup* ags)
{
	int i;
        /* save snapshot of coordinates at the update stage. */
        for(i=0; i<ag->natoms; i++)
        {
          ags->nblst->crds[3*i]=ag->atoms[i].X;
	  ags->nblst->crds[3*i+1]=ag->atoms[i].Y;
	  ags->nblst->crds[3*i+2]=ag->atoms[i].Z;
        }
        // Access margin size of future cubes.
        /* sum of nblist cutoff and margin size gives the cubesize. */
        findmarg(ag, ags->clst);

        struct cubeset *cust = _mol_malloc(sizeof(struct cubeset));
        // Generate a set of cubes cust based on the clusterset ags->clst.
        /* 0. sum of nblist cutoff and margin size gives the cubesize. 
            1. calculate coordinate span of cluster centers
            2. number of cubes in each direction is span over cubesize plus 1
            3. for each cluster find a cube it belongs to.
            4. each cube knows last cluster contained in it.
            5. each cluster knows another cluster belongin to the same cube (linked list)
            6. each cube knows its occupied neighbouring cubes. */
        gen_cubeset(ags->nblst->nbcut, ags->clst, cust);
        // Generate a nonbonded list nblst from cluster and cube subdivisions and exclusion list arrays.
        /* To generate nblist
            1. loop over the filled cubes
            2. loop over linked clusters in the first cube
            3. loop over neighboring cubes including first cube.
            4. loop through the linked list clusters of the second cube.
            5. go only for pairs of clusters which are closer than an nbcut plus lagest distance from any atom to the
               geometry center of cluster 1 plus the largest distance of the cluster 2
            6. for the pair of clusters icl1 and icl2 go through all pairs of their atoms
            7. write the pair to the list if the distance between clusters is less than nbcut minus largest dist of 
               cluster one minus largest dist of cluster two.
            8. check the atomic pair distance  if the condition in the step 7 is not true
            9. write the pair if distance less than nonbonded cutoff.

            The structure of the nblist is:
            number of first atoms
            array of first atom indices
            array of numbers of second atoms for each first one
            array of pointers to arrays of all second atoms
            arrays of all second atoms. */
        gen_nblist(ag, cust, ags->clst, ags->excl_list, ags->pd1, ags->pd2, ags->ndm, ags->nblst);
	free_cubeset(cust);
        free(cust);

}
//Checks wether cluster needs update
int check_clusterupdate(struct atomgrp* ag,struct agsetup* ags){
 int i=0,j;
 short flag=0;
 double delta0,delta;
 delta0=fabs(ags->nblst->nbcut-ags->nblst->nbcof)/2.0;
 delta0=delta0*delta0;
    for (i=0;i<ag->nactives;i++){
        j=ag->activelist[i];
        delta=(ags->nblst->crds[3*j]-ag->atoms[j].X)*(ags->nblst->crds[3*j]-ag->atoms[j].X)+
	   (ags->nblst->crds[3*j+1]-ag->atoms[j].Y)*(ags->nblst->crds[3*j+1]-ag->atoms[j].Y)+
	   (ags->nblst->crds[3*j+2]-ag->atoms[j].Z)*(ags->nblst->crds[3*j+2]-ag->atoms[j].Z);
        delta=delta*delta;
	if (delta>delta0){
	    flag=1;
	    break;
	}
	}
    if (flag){
	update_nblst(ag,ags);
	return 1;
	    } else {return 0;}
}

void give_012(struct atomgrp* ag, struct agsetup* ags,
              int* na012, int** la012, int* nf012, int** lf012)
{
     struct atom *a0, *a1;
     int ia0, ia1;
     int n01=ag->nbonds;
     int n02=ags->n02;
     int n=n01+n02;
     int na=0;
     int nf=0;

     if(n>0)
     {
        *la012=_mol_malloc(2*n*sizeof(int));
        *lf012=_mol_malloc(2*n*sizeof(int));
        int i;
        for(i=0; i<n01; i++)
        {
           a0=ag->bonds[i].a0;
           a1=ag->bonds[i].a1;
           if(a0->fixed == 0 || a1->fixed == 0)
           {
              *la012[2*na]=a0->ingrp;
              *la012[2*na+1]=a1->ingrp;
              na++;
           }
           else
           {
              *lf012[2*nf]=a0->ingrp;
              *lf012[2*nf+1]=a1->ingrp;
              nf++;
           }
        }
        for(i=0; i<n02; i++)
        {
           ia0=ags->list02[2*i];
           ia1=ags->list02[2*i+1];
           if(ag->atoms[ia0].fixed == 0 || ag->atoms[ia1].fixed == 0)
           {
              *la012[2*na]=ia0;
              *la012[2*na+1]=ia1;
              na++;
           }
           else
           {
              *lf012[2*nf]=ia0;
              *lf012[2*nf+1]=ia1;
              nf++;
           }
        }
        if(na>0)*la012=_mol_realloc(*la012, 2*na*sizeof(int));
        else free(*la012);
        if(nf>0)*lf012=_mol_realloc(*lf012, 2*nf*sizeof(int));
        else free(*lf012);
     }
     *na012=na;
     *nf012=nf;
}

void springeng(struct springset *sprst, double* spren)
{
     int i,j,nat,lj;
     double xtot,ytot,ztot,fk;
     struct atomgrp *ag;
     for (i=0; i<sprst->nsprings; i++)
     {
         nat=sprst->springs[i].naspr;
         if(nat>0)
         {
            xtot=0.0;
            ytot=0.0;
            ztot=0.0;
            ag =sprst->springs[i].agp;
            for (j=0; j<nat; j++)
            {
               lj=sprst->springs[i].laspr[j];
               xtot+=ag->atoms[lj].X;
               ytot+=ag->atoms[lj].Y;
               ztot+=ag->atoms[lj].Z;
            }
            xtot=xtot/nat-sprst->springs[i].X0;
            ytot=ytot/nat-sprst->springs[i].Y0;
            ztot=ztot/nat-sprst->springs[i].Z0;
            fk=sprst->springs[i].fkspr;
            *spren+=fk*(xtot*xtot+ytot*ytot+ztot*ztot);
            fk=2*fk/nat;
            for (j=0; j<nat; j++)
            {
                lj=sprst->springs[i].laspr[j];
                ag->atoms[lj].GX -= xtot*fk;
                ag->atoms[lj].GY -= ytot*fk;
                ag->atoms[lj].GZ -= ztot*fk;
            }
         }
     }
}

void modify_vdw_save(struct atomgrp* ag, int nmod, int* mod, double* enew, double* rnew,
		                                             double* eold, double* rold)
{
   int i;
   for(i=0; i<nmod; i++)
   {
      eold[i]=ag->atoms[mod[i]].eps;
      if(enew[i]<=0)ag->atoms[mod[i]].eps=sqrt(-enew[i]);
      rold[i]=ag->atoms[mod[i]].rminh;
      ag->atoms[mod[i]].rminh=rnew[i];
   }

}

void modify_vdw_back(struct atomgrp* ag, int nmod, int* mod, double* eold, double* rold)
{
   int i;
   for(i=0; i<nmod; i++)
   {
      ag->atoms[mod[i]].eps=eold[i];
      ag->atoms[mod[i]].rminh=rold[i];
   }

}

void get_mod_vdw_all(double lambe, double lambr, struct atomgrp* ag, int *nmod, int **mod, double **modeps, double **modrminh)
{
   int i;
   *nmod=ag->natoms;
   *mod=mymalloc(*nmod*sizeof(int));
   *modeps=mymalloc(*nmod*sizeof(double));
   *modrminh=mymalloc(*nmod*sizeof(double));
   for(i=0; i<*nmod; i++)
   {
      (*mod)[i]=i;
      (*modeps)[i]=lambe*(ag->atoms[i].eps);
      (*modrminh)[i]=lambr*(ag->atoms[i].rminh);
   }
}



/* ####################################################################### */
/* #                       START: ADDED BY REZAUL                        # */
/* ####################################################################### */


//#define DEBUG_VDW

/**
   This function computes the van der Waals interaction energy between the 
   atoms of the two octree nodes specified in 'OCTREE_PARAMS *octpar' using 
   the distance cutoff value and other parameters supplied through the same 
   structure. The structure contains two octree pointers 'octpar->octree_static' 
   and 'octpar->octree_moving', but this function assumes that either both point 
   to the same octree, or 'octpar->octree_moving' includes only the nonfixed
   atoms of 'octpar->octree_static'. Interaction energy
   is computed between the nonfixed atoms of the octree node with
   index 'octpar->node_moving' and all atoms of the node with
   index 'octpar->node_static'. The energy value is returned 
   in 'double *energy'.
*/
void vdweng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy )
{
    OCTREE *octree_static = octpar->octree_static;
    OCTREE *octree_moving = octpar->octree_moving;
    double *trans_mat = octpar->trans;

    OCTREE_PARAMS *prms = ( OCTREE_PARAMS * ) octpar->proc_func_params;
    struct agsetup *ags = prms->ags;

    OCTREE_NODE *snode = &( octree_static->nodes[ octpar->node_static ] );  /* the static node */
    OCTREE_NODE *mnode = &( octree_moving->nodes[ octpar->node_moving ] );  /* the moving node */

    double rc2 = octpar->dist_cutoff * octpar->dist_cutoff;
    double irc2 = 1 / rc2;

    int nf = snode->nfixed;

    double en = 0;

#ifdef DEBUG_VDW
    int m = 0, n = 0;
    double minD2 = 100000000000.0;
    double maxD2 = 0;
#endif

    /* if no transform matrix is given or the number of nonfixed atoms 
       in the moving node is smaller than the total number of atoms in the 
       static node, then iterate through the nonfixed atoms in the moving 
       node in the outer loop */
    if ( ( trans_mat != NULL ) || ( mnode->n - mnode->nfixed <= snode->n ) )
      {        
        /* iterate through the nofixed atoms of the moving node */
        for ( int i = mnode->nfixed; i < mnode->n; i++ )
          {
            int ai = mnode->indices[ i ];
            mol_atom *atom_i = &( octree_moving->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
      
            /* transform the atomic positions if a transformation matrix is given */ 
            if ( trans_mat != NULL ) transform_point( x, y, z, trans_mat, &x, &y, &z );    // defined in octree.h
      
            /* calculate the smallest distance from the atom center < x, y, z >
               to the static node bounding box */
            double d2 = min_pt2bx_dist2( snode->lx, snode->ly, snode->lz, snode->dim, x, y, z );

#ifdef DEBUG_VDW                        
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif
      
            /* if not definitely outside the distance cutoff */
            if ( d2 < rc2 )
              {
#ifdef DEBUG_VDW              
                m++;
#endif      
                double e_i = atom_i->eps, r_i = atom_i->rminh;
      
                double gx = 0, gy = 0, gz = 0;
              
                /* iterate through all atoms of the static node */        
                for ( int j = 0; j < snode->n; j++ )
                  {
                    int aj = snode->indices[ j ];
      
                    if ( ( j >= nf ) && ( aj <= ai ) ) continue; /* avoing double-counting nonfixed-nonfixed interactions */
      
                    mol_atom *atom_j = &( octree_static->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;
      
                    /* squared distance between atom_i (of moving node) and atom_j (of static node) */
                    d2 = dx * dx + dy * dy + dz * dz;
            
                    /* ignore if outside distance cutoff */        
                    if ( d2 >= rc2 ) continue;
      
                    /* if too close check for bonds */
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;  /* ignore atom pairs within three bond distance */
                      }

#ifdef DEBUG_VDW      
                    n++;
#endif      
                    /* compute 6-12 vdW interaction */
                    
                    double id2 = 1 / d2;
      
                    double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
                    r_ij *= r_ij;
      
                    double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;
      
                    Rd6 = Rd6 * Rd6 * Rd6;
                    Rr6 = Rr6 * Rr6 * Rr6;
                    dr6 = dr6 * dr6 * dr6;
      
                    double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;
      
                    /* compute energy */
                    en += e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );

                    /* compute gradients */      
                    double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;

                    double g = dven * dx;
      
                    gx += g;
                    ( atom_j->GX ) -= g;
      
                    g = dven * dy;
      
                    gy += g;
                    ( atom_j->GY ) -= g;
      
                    g = dven * dz;
      
                    gz += g;
                    ( atom_j->GZ ) -= g;
                  }
      
                /* update gradients */
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }
    else
      {
        /* iterate through all atoms of the static node */
        for ( int i = 0; i < snode->n; i++ )
          {
            int ai = snode->indices[ i ];
            mol_atom *atom_i = &( octree_static->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
      
            /* calculate the smallest distance from the atom center < x, y, z >
               to the moving node bounding box */
            double d2 = min_pt2bx_dist2( mnode->lx, mnode->ly, mnode->lz, mnode->dim, x, y, z );
      
#ifdef DEBUG_VDW            
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif

            /* if not definitely outside the distance cutoff */      
            if ( d2 < rc2 )
              {
#ifdef DEBUG_VDW              
                m++;
#endif                
      
                double e_i = atom_i->eps, r_i = atom_i->rminh;
      
                double gx = 0, gy = 0, gz = 0;
      
                /* iterate through the nonfixed atoms of the moving node */ 
                for ( int j = mnode->nfixed; j < mnode->n; j++ )
                  {
                    int aj = mnode->indices[ j ];
      
                    if ( ( i >= nf ) && ( aj >= ai ) ) continue;  /* avoing double-counting nonfixed-nonfixed interactions */
      
                    mol_atom *atom_j = &( octree_moving->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;

                    /* squared distance between atom_i (of moving node) and atom_j (of static node) */      
                    d2 = dx * dx + dy * dy + dz * dz;

                    /* ignore if outside distance cutoff */                            
                    if ( d2 >= rc2 ) continue;

                    /* if too close check for bonds */      
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;  /* ignore atom pairs within three bond distance */
                      }

#ifdef DEBUG_VDW      
                    n++;
#endif                    

                    /* compute 6-12 vdW interaction */      
                    
                    double id2 = 1 / d2;
      
                    double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
                    r_ij *= r_ij;
      
                    double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;
      
                    Rd6 = Rd6 * Rd6 * Rd6;
                    Rr6 = Rr6 * Rr6 * Rr6;
                    dr6 = dr6 * dr6 * dr6;
      
                    double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;
      
                    /* compute energy */      
                    en += e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );
      
                    /* compute gradients */            
                    double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;
      
                    double g = dven * dx;
      
                    gx += g;
                    ( atom_j->GX ) -= g;
      
                    g = dven * dy;
      
                    gy += g;
                    ( atom_j->GY ) -= g;
      
                    g = dven * dz;
      
                    gz += g;
                    ( atom_j->GZ ) -= g;
                  }
      
                /* update gradients */      
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }  

#ifdef DEBUG_VDW
    int d = within_distance_cutoff( snode->lx, snode->ly, snode->lz, snode->dim,
                                    mnode->lx, mnode->ly, mnode->lz, mnode->dim,
                                    octpar->dist_cutoff );                                                           
                                         
    printf( "< %2d %2d %9lf > < %2d %2d %9lf > %2d %3d %15.12lf %9lf %9lf %d ( %9lf, %9lf, %9lf | %9lf %9lf %9lf )\n", 
            mnode->nfixed, mnode->n, mnode->dim, snode->nfixed, snode->n, snode->dim, m, n, en, sqrt( minD2 ), sqrt( maxD2 ), 
            d, snode->lx, snode->ly, snode->lz, mnode->lx, mnode->ly, mnode->lz );
#endif

    *energy = en;
}


//#define DEBUG_ELENG

/**
   This function computes the electrostatic interaction energy between the 
   atoms of the two octree nodes specified in 'OCTREE_PARAMS *octpar' using 
   the distance cutoff value and other parameters supplied through the same 
   structure. The structure contains two octree pointers 'octpar->octree_static' 
   and 'octpar->octree_moving', but this function assumes that either both point 
   to the same octree, or 'octpar->octree_moving' includes only the nonfixed
   atoms of 'octpar->octree_static'. Interaction energy
   is computed between the nonfixed atoms of the octree node with
   index 'octpar->node_moving' and all atoms of the node with
   index 'octpar->node_static'. The energy value is returned 
   in 'double *energy'.
*/
void eleng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy )
{
    OCTREE *octree_static = octpar->octree_static;
    OCTREE *octree_moving = octpar->octree_moving;
    double *trans_mat = octpar->trans;

    OCTREE_PARAMS *prms = ( OCTREE_PARAMS * ) octpar->proc_func_params;
    struct agsetup *ags = prms->ags;

    OCTREE_NODE *snode = &( octree_static->nodes[ octpar->node_static ] );  /* the static node */
    OCTREE_NODE *mnode = &( octree_moving->nodes[ octpar->node_moving ] );  /* the moving node */

    double pf = CCELEC / prms->eps;
    double rc = octpar->dist_cutoff;
    double rc2 = rc * rc;
    double irc2 = 1 / rc2;

    int nf = snode->nfixed;

    double en = 0;

#ifdef DEBUG_ELENG
    int m = 0, n = 0;
    double minD2 = 100000000000.0;
    double maxD2 = 0;
#endif

    /* if no transform matrix is given or the number of nonfixed atoms 
       in the moving node is smaller than the total number of atoms in the 
       static node, then iterate through the nonfixed atoms in the moving 
       node in the outer loop */
    if ( ( trans_mat != NULL ) || ( mnode->n - mnode->nfixed <= snode->n ) )
      {
        /* iterate through the nofixed atoms of the moving node */      
        for ( int i = mnode->nfixed; i < mnode->n; i++ )
          {
            int ai = mnode->indices[ i ];
            mol_atom *atom_i = &( octree_moving->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;

            /* transform the atomic positions if a transformation matrix is given */       
            if ( trans_mat != NULL ) transform_point( x, y, z, trans_mat, &x, &y, &z );    // defined in octree.h
      
            /* calculate the smallest distance from the atom center < x, y, z >
               to the static node bounding box */
            double d2 = min_pt2bx_dist2( snode->lx, snode->ly, snode->lz, snode->dim, x, y, z );

#ifdef DEBUG_ELENG                        
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif

            /* if not definitely outside the distance cutoff */      
            if ( d2 < rc2 )
              {
#ifdef DEBUG_ELENG              
                m++;
#endif      
                double c_i = pf * atom_i->chrg;
      
                double gx = 0, gy = 0, gz = 0;

                /* iterate through all atoms of the static node */                              
                for ( int j = 0; j < snode->n; j++ )
                  {
                    int aj = snode->indices[ j ];
      
                    if ( ( j >= nf ) && ( aj <= ai ) ) continue;   /* avoing double-counting nonfixed-nonfixed interactions */
      
                    mol_atom *atom_j = &( octree_static->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;

                    /* squared distance between atom_i (of moving node) and atom_j (of static node) */      
                    d2 = dx * dx + dy * dy + dz * dz;

                    /* ignore if outside distance cutoff */                            
                    if ( d2 >= rc2 ) continue;

                    /* if too close check for bonds */      
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;  /* ignore atom pairs within three bond distance */
                      }

#ifdef DEBUG_ELENG      
                    n++;
#endif      
                    /* compute electrostatic interaction */
                    double d1 = sqrt( d2 );
                    double id2 = 1 / d2;
                    
                    double c_ij = c_i * atom_j->chrg;
                    
                    double esh = 1.0 - d1 / rc, desh = ( id2 - irc2 ) / d1;
                    
                    esh *= ( esh / d1 );
                    
                    /* compute energy */      
                    en += c_ij * esh;

                    /* compute gradients */                          
                    double g = c_ij * desh * dx;
                    
                    gx -= g;
                    ( atom_j->GX ) += g;
                    
                    g = c_ij * desh * dy;
                    
                    gy -= g;
                    ( atom_j->GY ) += g;
                    
                    g = c_ij * desh * dz;
                    
                    gz -= g;
                    ( atom_j->GZ ) += g;
                  }
      
                /* update gradients */
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }
    else
      {
        /* iterate through all atoms of the static node */      
        for ( int i = 0; i < snode->n; i++ )
          {
            int ai = snode->indices[ i ];
            mol_atom *atom_i = &( octree_static->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;

            /* calculate the smallest distance from the atom center < x, y, z >
               to the moving node bounding box */      
            double d2 = min_pt2bx_dist2( mnode->lx, mnode->ly, mnode->lz, mnode->dim, x, y, z );
      
#ifdef DEBUG_ELENG            
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif
      
            /* if not definitely outside the distance cutoff */            
            if ( d2 < rc2 )
              {
#ifdef DEBUG_ELENG              
                m++;
#endif                      
                double c_i = pf * atom_i->chrg;
      
                double gx = 0, gy = 0, gz = 0;

                /* iterate through the nonfixed atoms of the moving node */       
                for ( int j = mnode->nfixed; j < mnode->n; j++ )
                  {
                    int aj = mnode->indices[ j ];
      
                    if ( ( i >= nf ) && ( aj >= ai ) ) continue;  /* avoing double-counting nonfixed-nonfixed interactions */
      
                    mol_atom *atom_j = &( octree_moving->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;
      
                   /* squared distance between atom_i (of moving node) and atom_j (of static node) */      
                    d2 = dx * dx + dy * dy + dz * dz;
                    
                    /* ignore if outside distance cutoff */                                                
                    if ( d2 >= rc2 ) continue;

                    /* if too close check for bonds */            
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;  /* ignore atom pairs within three bond distance */
                      }

#ifdef DEBUG_ELENG      
                    n++;
#endif                    
                    /* compute electrostatic interaction */
                    double d1 = sqrt( d2 );
                    double id2 = 1 / d2;
                    
                    double c_ij = c_i * atom_j->chrg;
                    
                    double esh = 1.0 - d1 / rc, desh = ( id2 - irc2 ) / d1;
                    
                    esh *= ( esh / d1 );

                    /* compute energy */                          
                    en += c_ij * esh;

                    /* compute gradients */                                
                    double g = c_ij * desh * dx;
                    
                    gx -= g;
                    ( atom_j->GX ) += g;
                    
                    g = c_ij * desh * dy;
                    
                    gy -= g;
                    ( atom_j->GY ) += g;
                    
                    g = c_ij * desh * dz;
                    
                    gz -= g;
                    ( atom_j->GZ ) += g;                          
                  }

                /* update gradients */            
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }  

#ifdef DEBUG_ELENG
    int d = within_distance_cutoff( snode->lx, snode->ly, snode->lz, snode->dim,
                                    mnode->lx, mnode->ly, mnode->lz, mnode->dim,
                                    octpar->dist_cutoff );                                                           
                                         
    printf( "< %2d %2d %9lf > < %2d %2d %9lf > %2d %3d %15.12lf %9lf %9lf %d ( %9lf, %9lf, %9lf | %9lf %9lf %9lf )\n", 
            mnode->nfixed, mnode->n, mnode->dim, snode->nfixed, snode->n, snode->dim, m, n, en, sqrt( minD2 ), sqrt( maxD2 ), 
            d, snode->lx, snode->ly, snode->lz, mnode->lx, mnode->ly, mnode->lz );
#endif

    *energy = en;
}


/* ####################################################################### */
/* #                        END: ADDED BY REZAUL                         # */
/* ####################################################################### */
