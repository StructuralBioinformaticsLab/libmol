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
#ifndef _MOL_GBSA_H_
#define _MOL_GBSA_H_
struct acesetup {
    int ntypes;
    int nbsize;//Size of nblist
    double efac;
    int* list0123;
    double* eself;
    double* rborn;
    double* swarr;
    double* dswarr;
    double* darr;
    //d(rb)/d(eself)
    double* dbrdes;
    double* xsf;
    double* ysf;
    double* zsf;
    double* xf;
    double* yf;
    double* zf;
    double* diarr;
    double  *lwace,*rsolv,*vsolv,*s2ace,*uace,*wace,*hydr;
    int n0123;
   
};
//Initialize ace data types
void ace_ini(struct atomgrp* ag,struct acesetup* ac_s);
//Update ace lists once fixedlist was updated
void ace_fixedupdate(struct atomgrp* ag,struct agsetup* ags ,struct acesetup* ac_s);
//Update nblst once nblist is updated
void ace_updatenblst(const struct agsetup* const restrict ags, struct acesetup* const restrict ac_s);
//Calculate ace energy and gradients
void aceeng(struct atomgrp* ag,double *en,struct acesetup* ac_s,struct agsetup* ags);
void aceeng_nonpolar(struct atomgrp* ag,double* en,struct acesetup* ac_s,struct agsetup* ags);
void aceeng_polar(struct atomgrp* ag,double* en,struct acesetup* ac_s,struct agsetup* ags);
//Free ace data
void destroy_acesetup(struct acesetup* ac_s);
void free_acesetup(struct acesetup* ac_s);
//Test ace gradients
void test_acegrads(struct atomgrp *ag,struct agsetup* ags, struct acesetup* acs,double d);
void test_acegrads_nonpolar(struct atomgrp *ag,struct agsetup* ags, struct acesetup* acs,double d);
#endif
