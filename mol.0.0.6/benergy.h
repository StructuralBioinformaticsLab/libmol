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
#ifndef _MOL_BENERGY_H_
#define _MOL_BENERGY_H_

/** \file benergy.h
        This file contains  functions
        for calculating bonded energies and forces
	(bonds, angles, dihedrals, impropers)
*/

/**
  find the bond energy and forces
*/
void beng(struct atomgrp *ag, double* en);

/**
  find the angle energy and forces
*/
void aeng(struct atomgrp *ag, double* en);

/**
  find the improper energy and forces
*/
void ieng(struct atomgrp *ag, double* en);

/**
  find the dihedral energy and forces
*/
void teng(struct atomgrp *ag, double* en);

void zero_grads(struct atomgrp *ag);

void check_b_grads(struct atomgrp *ag, double d,
                void (*efun)(struct atomgrp *, double*));
void check_speng_grads(int nstart, int nend,
                struct atomgrp *ag, double d, double stens,
                double* hx0, double* hy0, double* hz0, 
                int nx, int ny, int nz,
                double dcx, double dcy, double dcz,
                double cx, double cy, double cz, double w,
                void (*efun)(double, double, struct atomgrp *, double*,
                double*, double*, double*, 
                int,int,int,
                double,double,double,
                double,double,double,double));


#endif
