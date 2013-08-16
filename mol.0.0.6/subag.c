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
#ifdef _WIN32
#include "../mol.0.0.6.h"
#else
#include _MOL_INCLUDE_
#endif

void setup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist,
		 int *list)
{
	int i;
	for (i = 0; i < nlist; i++) {
		int j = list[i];
		ag->atoms[j].X = new_actives->atoms[i].X;
		ag->atoms[j].Y = new_actives->atoms[i].Y;
		ag->atoms[j].Z = new_actives->atoms[i].Z;
	}
}

void unsetup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist,
		   int *list)
{
	int i;
	for (i = 0; i < nlist; i++) {
		int j = list[i];
		new_actives->atoms[i].X = ag->atoms[j].X;
		new_actives->atoms[i].Y = ag->atoms[j].Y;
		new_actives->atoms[i].Z = ag->atoms[j].Z;
	}
}
