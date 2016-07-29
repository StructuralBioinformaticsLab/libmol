/*
Copyright (c) 2010-2012, Structural Bioinformatics Laboratory, Boston University
Copyright (c) 2013, Acpharis Inc
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
#include <math.h>

#include _MOL_INCLUDE_

void get_detransform(const struct atomgrp *ag, const int residue_index,
		     struct mol_matrix3f *rot_to_bb,
		     struct mol_vector3f *trans_to_bb)
{
	struct mol_vector3f e1, e2, e3;
	float d;

	struct atom *nl = &(ag->atoms[ag->iares[residue_index]]);
	struct atom *cal = &(ag->atoms[ag->iares[residue_index] + 2]);
	struct atom *cl = &(ag->atoms[ag->iares[residue_index + 1] - 2]);	//second to last atom in res 

	e1.X = nl->X - cal->X;
	e1.Y = nl->Y - cal->Y;
	e1.Z = nl->Z - cal->Z;

	//normalize e1
	d = sqrtf((e1.X * e1.X) + (e1.Y * e1.Y) + (e1.Z * e1.Z));
	e1.X /= d;
	e1.Y /= d;
	e1.Z /= d;

	e2.X = cl->X - cal->X;
	e2.Y = cl->Y - cal->Y;
	e2.Z = cl->Z - cal->Z;

	//Orthogonalize:
	//dot product to get component of e1 along e2
	d = (e1.X * e2.X) + (e1.Y * e2.Y) + (e1.Z * e2.Z);
	printf("d = %f\n", d);
	//subtract component
	e2.X -= e1.X * d;
	e2.Y -= e1.Y * d;
	e2.Z -= e1.Z * d;

	//normalize e2
	d = sqrtf((e2.X * e2.X) + (e2.Y * e2.Y) + (e2.Z * e2.Z));
	e2.X /= d;
	e2.Y /= d;
	e2.Z /= d;

	//generate third orthonormal component as cross product
	e3.X = (e1.Y * e2.Z) - (e1.Z * e2.Y);
	e3.Y = (e1.Z * e2.X) - (e1.X * e2.Z);
	e3.Z = (e1.X * e2.Y) - (e1.Y * e2.X);

	//assign rotation matrix and translation vector
	rot_to_bb->m11 = e1.X;
	rot_to_bb->m12 = e2.X;
	rot_to_bb->m13 = e3.X;
	rot_to_bb->m21 = e1.Y;
	rot_to_bb->m22 = e2.Y;
	rot_to_bb->m23 = e3.Y;
	rot_to_bb->m31 = e1.Z;
	rot_to_bb->m32 = e2.Z;
	rot_to_bb->m33 = e3.Z;

	trans_to_bb->X = cal->X;
	trans_to_bb->Y = cal->Y;
	trans_to_bb->Z = cal->Z;
}
