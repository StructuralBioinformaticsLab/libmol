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
#ifndef _MOL_SASA_H_
#define _MOL_SASA_H_

/*
   wrapper for baccs: fill sa term in ag
   using a common set of parameters
*/
//void msur (struct atomgrp* ag, struct prm* prm);
// msur_k == -1 => don't msur
void msur (struct atomgrp* ag, struct prm* prm, float msur_k);

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
                float r_solv, short cont_acc, short rpr, float sthresh);
void baccs1(struct atomgrp* ag, int n_at, int* restat,
                double r_solv, short cont_acc, float sthresh);

/* atoms with zero radii in prm structure are excluded from calculation */
/* r_solv -solvent radius */
/* cont_acc =1 - contact surface
			=0 - accessible surface */
/* rpr      =1 - use radii from the parameter structure prm
			=0 - use preset radii */				
/* as - atomic surface (output) */
void accs (struct atomgrp* ag, struct prm* prm, float r_solv, short cont_acc, short rpr, float* as);
void accs1 (struct atomgrp* ag, int n_at, int* restat, double r_solv, short cont_acc, double* as);
void mark_sasa (struct atomgrp* ag, int* sasas);

int* read_sasa (const char* path);

int numbersasas (int* sasas);

#endif
