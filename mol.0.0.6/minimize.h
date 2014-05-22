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
#ifndef _MOL_MINIMIZE_H_
#define _MOL_MINIMIZE_H_

/** \file minimize.h
        This file contains  functions
        for local structure minimization
*/

#include "lbfgs.h"

typedef enum
{
    MOL_LBFGS,
    MOL_CONJUGATE_GRADIENTS,
    MOL_POWELL,
    MOL_RIGID,
} mol_min_method;


//Wrapper for atom group minimizers;
void minimize_ag(mol_min_method min_type, unsigned int maxIt, double tol, struct atomgrp* ag,void* minprms, void (*egfun)(int , double* , void* , double* , double*));

// Powell/dirPowell below

int progress(void *instance,
const lbfgsfloatval_t *x,
const lbfgsfloatval_t *g,
const lbfgsfloatval_t fx,
const lbfgsfloatval_t xnorm,
const lbfgsfloatval_t gnorm,
const lbfgsfloatval_t step,
int n,int k,int ls);

void bracket(double* orig, double* dir, double step,
             int ndim, void* prms,
             void (*egfun)(int , double* , void* , double* , double*),
             double *fb,  double *la, double *lb, double *lc);

void brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, unsigned int maxtimes,
           double*  min, double* fmim);

void dirbrent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, int maxtimes,
           double*  min, double* fmim, double*  grad);

void powell(double* orig, double* directions, unsigned int maxIt, double tol,
            int ndim, void* prms,
            void (*egfun)(int , double* , void* , double*, double* ),
            double* min, double* fmim);

void dirMin(double* orig, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim);

void limin(double* orig, double* dir, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim);

void old_brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double*, double* ),
           double tol, int maxtimes,
           double*  min, double* fmim);
// Powell/dirPowell above

#endif
