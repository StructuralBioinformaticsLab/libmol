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
#include <math.h>

#include _MOL_INCLUDE_

/* comparison function for qsort */
//int accs_comp(const void *s1, const void *s2)
int accs_comp(const void *s1, const void *s2)
{
     const float f1=*(const float *)s1;
     const float f2=*(const float *)s2;
     if(f1<f2)return -1;
     if(f1==f2)return 0;
     return 1;
}

/* comparison function for qsort double version*/
//int accs_comp(const void *s1, const void *s2)
int accs_comp1(const void *s1, const void *s2)
{
     const double f1=*(const double *)s1;
     const double f2=*(const double *)s2;
     if(f1<f2)return -1;
     if(f1==f2)return 0;
     return 1;
}

void mark_sasa (struct atomgrp* ag, int* sasas)
{
	int nsasas = numbersasas (sasas);
	if (ag->natoms != nsasas)
	{
		fprintf (stderr, "error: ag->natoms (%d) != nsasas (%d)\n", ag->natoms, nsasas);
		exit (EXIT_FAILURE);
	}

	int atomi;
	for (atomi = 0; atomi < ag->natoms; atomi++)
	{
		ag->atoms[atomi].sa = sasas[atomi];
	}
}

int* read_sasa (const char* path)
{
	FILE* fp = myfopen (path, "r");

	int nsasas = 100; // just a guess, realloc as necessary
	int* sasas = (int*) _mol_malloc (sizeof (int) * nsasas);

	char* line = NULL;
	size_t len = 0;


	int i = 0;
	while (getline (&line, &len, fp) != -1)
	{
		if (i+1 > nsasas)
		{
			nsasas *= 2;
			sasas = (int*) _mol_realloc (sasas, sizeof (int) * nsasas);
		}

		if (sscanf (line, "%d", &sasas[i]) < 1)
		{
			fprintf (stderr, "error: in file %s line %s: incorrect sasa line\n", path, line);
			exit (EXIT_FAILURE);
		}
		if (sasas[i] != 0 && sasas[i] != 1)
		{
			fprintf (stderr, "error: in file %s line %s integer should be 0 or 1\n", path, line);
			exit (EXIT_FAILURE);
		}

		i++;
	}
	if (line)
		free (line);
	myfclose (fp);

	// final realloc of the arrays to make them tight
	nsasas = i;
	sasas = (int*) _mol_realloc (sasas, sizeof (int) * (nsasas+1)); // one extra for -1
	sasas[nsasas] = -1;

	return sasas;
}

int numbersasas (int* sasas)
{
	int i = 0;
	while (sasas[i] != -1)
	{
		i++;
	}
	return i;
}

//void msur (struct atomgrp* ag, struct prm* prm)
void msur (struct atomgrp* ag, struct prm* prm, float msur_k)
{
	if (msur_k < 0)
		return;

	float r_solv=1.4*msur_k;

	struct prm* rkprm = copy_prm (prm);
	modify_prms_radii (rkprm, msur_k);
	float sthresh=0.0;
	short cont_acc=1;
	short rpr=1;
	baccs(ag, rkprm, r_solv, cont_acc, rpr, sthresh);

	free_prm (rkprm); // (free_prms needs to be fixed)
}

/* convert float atomic surface arrays into 0/1 representation (sa in the structure atom) */
/* atoms with zero radii in prm structure are excluded from calculation */
/* r_solv -solvent radius */
/* cont_acc =1 - contact surface
			=0 - accessible surface */
/* rpr      =1 - use radii from the parameter structure prm
			=0 - use preset radii */
/* sthresh   if as > sthresh sa=1
	     if as<= sthresh sa=0   */
void baccs(struct atomgrp* ag, struct prm* prm, 
		float r_solv, short cont_acc, short rpr, float sthresh)
{
	int n_at=ag->natoms;
	int i;
	float* as=_mol_malloc(n_at*sizeof(float));
	accs(ag, prm, r_solv, cont_acc, rpr, as);

	for(i=0; i<n_at; i++)ag->atoms[i].sa=as[i]>sthresh?1:0;
	free(as);
}

/*This function does the same thing as baccs except it call accs1 (the traditional
surface area accesibility algorith) instead of accs (in-house protein-protein 
concerned algorithm).*/
void baccs1(struct atomgrp* ag, int n_at, int* restat,
		double r_solv, short cont_acc, float sthresh)
{
	int i;
	double* as=_mol_malloc(n_at*sizeof(double));
	
	accs1(ag, n_at, restat, r_solv, cont_acc, as);

	for(i=0; i<n_at; i++)ag->atoms[i].sa=as[i]>sthresh?1:0;
	free(as);
}

/* atoms with zero radii in prm structure are excluded from calculation */
/* r_solv -solvent radius */
/* cont_acc =1 - contact surface
			=0 - accessible surface */
/* rpr      =1 - use radii from the parameter structure prm
			=0 - use preset radii */				
/* as - atomic surface (output) */
void accs (struct atomgrp* ag, struct prm* prm, float r_solv, short cont_acc, short rpr, float* as)
{
	const int NAC=5000;  /* max number of atoms in a cube */

/* radii in preset mode */
	int at;
	char an;
	const float R_C=1.9;
	const float R_N=1.7;
	const float R_O=1.4;
	const float R_S=1.8;
	const float R_H=0.8;
	const float R_ELSE=1.9;

/* integration increment */
	const float P=0.01;

	int i, sint=sizeof(int), sflo=sizeof(float);
	int n_at1=ag->natoms;
	float  xmin, xmax, ymin, ymax, zmin, zmax;
	float rmax=0;

	float pi=acos(-1.0);
	float pix2=2.0*pi;
	float ri, xi, yi, zi;

/* eliminate atoms with zero radii */
	int n_at=0;
	i=n_at1*sint;
	int* restat=_mol_malloc(i);
	for(i=0; i<n_at1; i++)
	{
		as[i]=0.0;
		at=ag->atoms[i].atom_typen;
		ri=prm->atoms[at].r;
		if(ri>0.0)restat[n_at++]=i;
	}

/* initialize boundary constants */
	xmin=ag->atoms[restat[0]].X;
	xmax=xmin;
	ymin=ag->atoms[restat[0]].Y;
	ymax=ymin;
	zmin=ag->atoms[restat[0]].Z;
	zmax=zmin;

/* allocate general atom related arrays */
	i=n_at*sflo;
	float* x=_mol_malloc(i);
	float* y=_mol_malloc(i);
	float* z=_mol_malloc(i);
	float* r=_mol_malloc(i);
	float* r2=_mol_malloc(i);

/* allocate arrays for neighbouring atoms */
	float* dx=_mol_malloc(i);
	float* dy=_mol_malloc(i);
	float* d=_mol_malloc(i);
	float* dsq=_mol_malloc(i);
	float* arcif=_mol_malloc(2*2*i);
	i=n_at*sint;
	int* inov=_mol_malloc(i); 

/* calculate sizes and dimensions*/
	for(i=0; i<n_at; i++)
	{
		at=ag->atoms[restat[i]].atom_typen;
		if(rpr)
		{
			ri=prm->atoms[at].r;
		}
		else
		{
			an=*(prm->atoms[at].typemin);
			if(an=='C')ri=R_C;
			else if(an=='N')ri=R_N;
			else if(an=='O')ri=R_O;
			else if(an=='S')ri=R_S;
			else if(an=='H')ri=R_H;
/*			{
				fprintf (stderr, "error: accs.c> hydrogens are not allowed with rpr=0\n");
				exit(EXIT_FAILURE);
			} */
			else ri=R_ELSE;
		}
		ri=ri+r_solv;
		r[i]=ri;
		r2[i]=ri*ri;
		if(ri>rmax)rmax=ri;
		x[i]=ag->atoms[restat[i]].X;
		if(xmin>x[i])xmin=x[i];
		if(xmax<x[i])xmax=x[i];
		y[i]=ag->atoms[restat[i]].Y;
		if(ymin>y[i])ymin=y[i];
		if(ymax<y[i])ymax=y[i];
		z[i]=ag->atoms[restat[i]].Z;
		if(zmin>z[i])zmin=z[i];
		if(zmax<z[i])zmax=z[i];
	}
	float dmax=rmax*2.0;

	i=(xmax-xmin)/dmax+1;
	int idim=i<3?3:i;

	i=(ymax-ymin)/dmax+1;
	int jidim=i<3?3:i;
	jidim*=idim;

	i=(zmax-zmin)/dmax+1;
	int kjidim=i<3?3:i;
	kjidim*=jidim;				/* total number of cubes */
	

/* map atoms to adjacent cubes */
/* allocate cubical arrays */

	i=kjidim*sint;
	int* itab=_mol_malloc(i);		/* number of atoms in each cube */
	for(i=0; i<kjidim; i++)itab[i]=0;

	i=NAC*kjidim*sint;
	int* natm=_mol_malloc(i);		/* atom index in each cube */

	i=n_at*sint;
	int* cube=_mol_malloc(i);		/* cube number for each atom */

	int j, k, l, m, n, kji;

	for(l=0; l<n_at; l++)
	{
		i=(x[l]-xmin)/dmax;
		j=(y[l]-ymin)/dmax;
		k=(z[l]-zmin)/dmax;
		kji=k*jidim+j*idim+i;	/* cube number */
		n=itab[kji]+1;	
		if(n>NAC)
		{
			_mol_error("number of atoms in a cube %d is above the maximum  NAC= %d\n", n, NAC);
			exit(EXIT_FAILURE);
		}
		itab[kji]=n;
		natm[kji*NAC+n-1]=l;
		cube[l]=kji;
	}
	

	int ir, io, in, mkji, nm, nzp, karc;
	float area, xr, yr, zr, rr, rrx2, rr2, b, zres, zgrid;
	float rsec2r, rsecr, rsec2n, rsecn;
	float calpha, alpha, beta, ti, tf, arcsum, parea, t, tt;


/* main loop over atoms */
	zi=1.0/P+0.5;
	nzp=zi;		 /* number of z planes */
	
	for (ir=0; ir<n_at; ir++)
	{
		kji=cube[ir];
		io=0;					/* number of neighbouring atoms */
		area=0.0;
		xr=x[ir];
		yr=y[ir];
		zr=z[ir];
		rr=r[ir];
		rrx2=rr*2;
		rr2=r2[ir];
/* loops over neighbouring cubes */
		for (k=-1; k<2; k++)
		{
			for(j=-1; j<2; j++)
			{
				for(i=-1; i<2; i++)
				{
					mkji=kji+k*jidim+j*idim+i;
					if(mkji<0)continue;
					if(mkji>=kjidim)goto esc_cubes;
					nm=itab[mkji];
					if(nm<1)continue;
					for(m=0; m<nm; m++)
					{
						in=natm[mkji*NAC+m];
						if(in!=ir)
						{
							xi=xr-x[in];
							yi=yr-y[in];
							dx[io]=xi;
							dy[io]=yi;
							ri=xi*xi+yi*yi;
							dsq[io]=ri;
							d[io]=sqrtf(ri);
							inov[io]=in;
							io++;
						}
					}
				}
			}
		}
		
		esc_cubes:
		if(io!=0)
		{
			zres=rrx2/nzp;	/* separation between planes */
			zgrid=z[ir]-rr-zres/2.0;	/* z level */

			for(i=0; i<nzp; i++)
			{
				zgrid+=zres;
/* radius of the circle intersection with a z-plane */
				zi=zgrid-zr;
				rsec2r=rr2-zi*zi;
				rsecr=sqrtf(rsec2r);

				karc=0;
				for(j=0; j<io; j++)
				{
					in=inov[j];
/* radius of the circle for a neighbour */
					zi=zgrid-z[in];
					rsec2n=r2[in]-zi*zi;
					if(rsec2n<=0.0)continue;
					rsecn=sqrtf(rsec2n);
/* are they close? */
					if(d[j]>=rsecr+rsecn)continue;
/* do they intersect? */
					b=rsecr-rsecn;
					if(b<=0.0)
					{
						if(d[j]<=-b)goto next_plane;
					}
					else
					{
						if(d[j]<=b)continue;
					}
					calpha=(dsq[j]+rsec2r-rsec2n)/(2.0*d[j]*rsecr);
					if(calpha>=1.0)continue;
/* yes, they do */
					alpha=acosf(calpha);
					beta=atan2f(dy[j],dx[j])+pi;
					ti=beta-alpha;
					tf=beta+alpha;
					if(ti<0.0)ti+=pix2;
					if(tf>pix2)tf-=pix2;
					arcif[karc]=ti;
					if(tf<ti)
					{
						arcif[karc+1]=pix2;
						karc+=2;
						arcif[karc]=0.0;
					}
					arcif[karc+1]=tf;
					karc+=2;
				}
/* find the atom accessible surface increment in z-plane */

				karc/=2;
				if(karc==0) 
					arcsum=pix2;
				else
				{
					qsort(arcif, karc, 2*sizeof(arcif[0]), accs_comp);
					arcsum=arcif[0];
					t=arcif[1];
					if(karc>1)
					{
						for(k=2; k<karc*2; k+=2)
						{
							if(t<arcif[k])arcsum+=(arcif[k]-t);
							tt=arcif[k+1];
							if(tt>t)t=tt;
						}
					}
					arcsum+=(pix2-t);
				}
				parea=arcsum*zres;
				area+=parea;
				next_plane:;
			}
		}
		else
		{
			area=pix2*rrx2;
		}

		ri=rr-r_solv;
		if(cont_acc)
			b=area*ri*ri/rr;
		else
			b=area*rr;
		as[restat[ir]]=b;
	}
/* free all */
	free(x);
	free(y);
	free(z);
	free(r);
	free(r2);
	free(dx);
	free(dy);
	free(d);
	free(dsq);
	free(arcif);
	free(inov);
	free(itab);
	free(natm);
	free(cube);
	free(restat);
}

/*This function is the unaltered, traditional surface area accessible algorithm.
Ryan Brenke altered accs to better account for protein-protein docking consideration.
Dmitri Beglov apparently restored accs to its original form.*/
void accs1 (struct atomgrp* ag, int n_at, int* restat, double r_solv, short cont_acc, double* as)
{
	const int NAC=800;  /* max number of atoms in a cube */

/* integration increment */
	const double P=0.01;

	int i, sint=sizeof(int), sdou=sizeof(double);
	int n_at1=ag->natoms;
	double  xmin, xmax, ymin, ymax, zmin, zmax;
	double rmax=0;

	double pi=acos(-1.0);
	double pix2=2.0*pi;
	double ri, xi, yi, zi;

	for(i=0; i<n_at1; i++)as[i]=0.0;

/* initialize boundary constants */
	xmin=ag->atoms[restat[0]].X;
	xmax=xmin;
	ymin=ag->atoms[restat[0]].Y;
	ymax=ymin;
	zmin=ag->atoms[restat[0]].Z;
	zmax=zmin;

/* allocate general atom related arrays */
	i=n_at*sdou;
	double* x=_mol_malloc(i);
	double* y=_mol_malloc(i);
	double* z=_mol_malloc(i);
	double* r=_mol_malloc(i);
	double* r2=_mol_malloc(i);

/* allocate arrays for neighbouring atoms */
	double* dx=_mol_malloc(i);
	double* dy=_mol_malloc(i);
	double* d=_mol_malloc(i);
	double* dsq=_mol_malloc(i);
	double* arcif=_mol_malloc(2*2*i);
	i=n_at*sint;
	int* inov=_mol_malloc(i);

/* calculate sizes and dimensions*/
	for(i=0; i<n_at; i++)
	{
		ri=ag->atoms[restat[i]].rminh;
		ri=ri+r_solv;
		r[i]=ri;
		r2[i]=ri*ri;
		if(ri>rmax)rmax=ri;
		x[i]=ag->atoms[restat[i]].X;
		if(xmin>x[i])xmin=x[i];
		if(xmax<x[i])xmax=x[i];
		y[i]=ag->atoms[restat[i]].Y;
		if(ymin>y[i])ymin=y[i];
		if(ymax<y[i])ymax=y[i];
		z[i]=ag->atoms[restat[i]].Z;
		if(zmin>z[i])zmin=z[i];
		if(zmax<z[i])zmax=z[i];
	}
	double dmax=rmax*2;

	i=(xmax-xmin)/dmax+1;
	int idim=i<3?3:i;

	i=(ymax-ymin)/dmax+1;
	int jidim=i<3?3:i;
	jidim*=idim;

	i=(zmax-zmin)/dmax+1;
	int kjidim=i<3?3:i;
	kjidim*=jidim;			  /* total number of cubes */


/* map atoms to adjacent cubes */
/* allocate cubical arrays */

	i=kjidim*sint;
	int* itab=_mol_malloc(i);	  /* number of atoms in each cube */
	for(i=0; i<kjidim; i++)itab[i]=0;

	i=NAC*kjidim*sint;
	int* natm=_mol_malloc(i);	  /* atom index in each cube */

	i=n_at*sint;
	int* cube=_mol_malloc(i);	  /* cube number for each atom */

	int j, k, l, m, n, kji;


	for(l=0; l<n_at; l++)
	{
		i=(x[l]-xmin)/dmax;
		j=(y[l]-ymin)/dmax;
		k=(z[l]-zmin)/dmax;
		kji=k*jidim+j*idim+i;   /* cube number */
		n=itab[kji]+1;
		if(n>NAC)
		{
			_mol_error("number of atoms in a cube %d is above the maximum  NAC= %d\n", n, NAC);
			exit(EXIT_FAILURE);
		}
		itab[kji]=n;
		natm[kji*NAC+n-1]=l;
		cube[l]=kji;
	}


	int ir, io, in, mkji, nm, nzp, karc;
	double area, xr, yr, zr, rr, rrx2, rr2, b, zres, zgrid;
	double rsec2r, rsecr, rsec2n, rsecn;
	double calpha, alpha, beta, ti, tf, arcsum, parea, t, tt;


/* main loop over atoms */
	zi=1.0/P+0.5;
	nzp=zi;		 /* number of z planes */

	for (ir=0; ir<n_at; ir++)
	{
		kji=cube[ir];
		io=0;				   /* number of neighbouring atoms */
		area=0.0;
		xr=x[ir];
		yr=y[ir];
		zr=z[ir];
		rr=r[ir];
		rrx2=rr*2;
		rr2=r2[ir];
/* loops over neighbouring cubes */
		for (k=-1; k<2; k++)
		{
			for(j=-1; j<2; j++)
			{
				for(i=-1; i<2; i++)
				{
					mkji=kji+k*jidim+j*idim+i;
					if(mkji<0)continue;
					if(mkji>=kjidim)goto esc_cubes;
					nm=itab[mkji];
					if(nm<1)continue;
					for(m=0; m<nm; m++)
					{
						in=natm[mkji*NAC+m];
						if(in!=ir)
						{
							xi=xr-x[in];
							yi=yr-y[in];
							dx[io]=xi;
							dy[io]=yi;
							ri=xi*xi+yi*yi;
							dsq[io]=ri;
							d[io]=sqrtf(ri);
							inov[io]=in;
							io++;
						}
					}
				}
			}
		}

		esc_cubes:
		if(io!=0)
		{
			zres=rrx2/nzp;  /* separation between planes */
			zgrid=z[ir]-rr-zres/2.0;	/* z level */

			for(i=0; i<nzp; i++)
			{
				zgrid+=zres;
/* radius of the circle intersection with a z-plane */
				zi=zgrid-zr;
				rsec2r=rr2-zi*zi;
				rsecr=sqrtf(rsec2r);

				karc=0;
				for(j=0; j<io; j++)
				{
					in=inov[j];
/* radius of the circle for a neighbour */
					zi=zgrid-z[in];
					rsec2n=r2[in]-zi*zi;
					if(rsec2n<=0.0)continue;
					rsecn=sqrtf(rsec2n);
/* are they close? */
					if(d[j]>=rsecr+rsecn)continue;
/* do they intersect? */
					b=rsecr-rsecn;
					if(b<=0.0)
					{
						if(d[j]<=-b)goto next_plane;
					}
					else
					{
						if(d[j]<=b)continue;
					}
/* yes, they do */
					calpha=(dsq[j]+rsec2r-rsec2n)/(2.*d[j]*rsecr);
					if(calpha >= 1.0)continue;

					alpha=acosf(calpha);
					beta=atan2f(dy[j],dx[j])+pi;
					ti=beta-alpha;
					tf=beta+alpha;
					if(ti<0.0)ti+=pix2;
					if(tf>pix2)tf-=pix2;
					arcif[karc]=ti;
					if(tf<ti)
					{
						arcif[karc+1]=pix2;
						karc+=2;
						arcif[karc]=0.0;
					}
					arcif[karc+1]=tf;
					karc+=2;
				}
/* find the atom accessible surface increment in z-plane */

				karc/=2;
				if(karc==0)
					arcsum=pix2;
				else
				{
					qsort(arcif, karc, 2*sizeof(arcif[0]), accs_comp1);
					arcsum=arcif[0];
					t=arcif[1];
					if(karc>1)
					{
						for(k=2; k<karc*2; k+=2)
						{
							if(t<arcif[k])arcsum+=(arcif[k]-t);
							tt=arcif[k+1];
							if(tt>t)t=tt;
						}
					}
					arcsum+=(pix2-t);
				}
				parea=arcsum*zres;
				area+=parea;
				next_plane:;
			}
		}
		else
		{
			area=pix2*rrx2;
		}

		ri=rr-r_solv;
		if(cont_acc)
			b=area*ri*ri/rr;
		else
			b=area*rr;
		as[restat[ir]]=b;
	}
/* free all */
	free(x);
	free(y);
	free(z);
	free(r);
	free(r2);
	free(dx);
	free(dy);
	free(d);
	free(dsq);
	free(arcif);
	free(inov);
	free(itab);
	free(natm);
	free(cube);
}

