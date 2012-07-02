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
#include _MOL_INCLUDE_

float* moment_of_inertia (struct atomgrp* ag)
{
	struct tvector* com = center_of_mass (ag);

	float sq_x_com = powf (com->X, 2.0);
	float sq_y_com = powf (com->Y, 2.0);
	float sq_z_com = powf (com->Z, 2.0);

	float sum_sq_x = 0;
	float sum_sq_y = 0;
	float sum_sq_z = 0;

	float sum_mult_xy = 0;
	float sum_mult_xz = 0;
	float sum_mult_yz = 0;

	int i;
	for (i = 0; i < ag->natoms; i++)
	{
		sum_sq_x += powf (ag->atoms[i].X, 2.0);
		sum_sq_y += powf (ag->atoms[i].Y, 2.0);
		sum_sq_z += powf (ag->atoms[i].Z, 2.0);

		sum_mult_xy += ag->atoms[i].X * ag->atoms[i].Y;
		sum_mult_xz += ag->atoms[i].X * ag->atoms[i].Z;
		sum_mult_yz += ag->atoms[i].Y * ag->atoms[i].Z;
	}

	float* moi_matrix = (float*) _mol_malloc (sizeof (float) * 9);

	moi_matrix[0] = sum_sq_y + sum_sq_z - (sq_y_com + sq_z_com) * ag->natoms;
	moi_matrix[1] = -sum_mult_xy + com->X * com->Y * ag->natoms;
	moi_matrix[2] = -sum_mult_xz + com->X * com->Z * ag->natoms;

	moi_matrix[3] = -sum_mult_xy + com->X * com->Y * ag->natoms;
	moi_matrix[4] = sum_sq_x + sum_sq_z - (sq_x_com + sq_z_com) * ag->natoms;
	moi_matrix[5] = -sum_mult_yz + com->Y * com->Z * ag->natoms;

	moi_matrix[6] = -sum_mult_xz + com->X * com->Z * ag->natoms;
	moi_matrix[7] = -sum_mult_yz + com->Y * com->Z * ag->natoms;
	moi_matrix[8] = sum_sq_x + sum_sq_y - (sq_x_com + sq_y_com) * ag->natoms;

	free (com);

	return moi_matrix;
}

void print_moment_of_inertia_matrix (float* moimat)
{
	printf ("%.3f\t%.3f\t%.3f\n", moimat[0], moimat[1], moimat[2]);
	printf ("%.3f\t%.3f\t%.3f\n", moimat[3], moimat[4], moimat[5]);
	printf ("%.3f\t%.3f\t%.3f\n", moimat[6], moimat[7], moimat[8]);
}
