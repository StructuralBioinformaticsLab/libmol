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
#ifndef _MOL_ENERGY_H_
#define _MOL_ENERGY_H_

/** \file energy.h
	This file contains structures and functions
	that are used for calculating offgrid energies
*/

/**
  find the pairwise potential on the two structures and return it
*/
float pairwise_potential_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm, int only_sab);

/**
  find the coulombic energy the two structures and return it
*/
float coulombic_elec_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm);

/**
  numerical module energy
*/
float nummod_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm);

void test_energy_grads(struct atomgrp *ag, void *minprms, void (*egfun)(int , double* , void* , double* , double*), double delta);

#endif
