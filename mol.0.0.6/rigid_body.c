/*
Copyright (c) 2009-2014, Structural Bioinformatics Laboratory, Boston University
Copyright (c) 2014, Acpharis Inc
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
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#include <math.h>
#include <string.h>

#include _MOL_INCLUDE_

#define Epsilon 1e-7

static double Norm(double v[3])
{
	return (sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
}

//---------------------------------------------
//Exp_To_Skew given an exponential map v will construt its corresponding
//Skew matrix and save it in Z

static void Exp_To_Skew(double v[3], double Z[3][3])
{
	Z[0][0] =   0.0;
	Z[0][1] = -v[2];
	Z[0][2] =  v[1];

	Z[1][0] =  v[2];
	Z[1][1] =   0.0;
	Z[1][2] = -v[0];

	Z[2][0] = -v[1];
	Z[2][1] =  v[0];
	Z[2][2] =   0.0;
}

static void Exp_To_Skew_Square(double v[3], double Z[3][3])
{
	Z[0][0] = -(v[2]*v[2] + v[1]*v[1]);
	Z[1][1] = -(v[2]*v[2] + v[0]*v[0]);
	Z[2][2] = -(v[1]*v[1] + v[0]*v[0]);
	Z[1][0] = Z[0][1] = v[1]*v[0];
	Z[2][0] = Z[0][2] = v[2]*v[0];
	Z[2][1] = Z[1][2] = v[2]*v[1];
}

//-------------------------------------------
//ConstMult given matrix X and double c will multiply the whole matrix by constant
static void ConstMult(double X[3][3], double c)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			X[i][j] *= c;
		}
	}
}

//-----------------------------------------
// Sum given two matrix X and Y will calculate sum of them 

static void Sum(double X[3][3], double Y[3][3], double Z[3][3])
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Z[i][j] = X[i][j] + Y[i][j];
		}
	}
}

//------------------------------------------------
//Exp_To_R given an exponential map v use Rodrigus formula and
// will generate the corresponding rotation matrix

static void Exp_To_R(double v[3], double R[3][3])
{
	double X[3][3];
	double Y[3][3];

	double theta = Norm(v);

	if (theta > Epsilon) {
		Exp_To_Skew(v, X);
		ConstMult(X, sin(theta) / theta);

		Exp_To_Skew_Square(v, Y);
		ConstMult(Y, (1 - cos(theta)) / (theta * theta));
	} else {
		Exp_To_Skew(v, X);
		ConstMult(X, 1 - (theta * theta / 6));

		Exp_To_Skew_Square(v, Y);
		ConstMult(Y, 0.5 - (theta * theta / 24));
	}

	Sum(X, Y, R);

	for (int i = 0; i < 3; i++) {
		R[i][i] += 1.0;
	}
}

//-------------------------------------------------
//Partial_Skew will compute the partial derivative of the skew matrix
// corresponding to the exponential map v wrt the component of v

static void Partial_Skew(double dZdv[3][3][3])
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				dZdv[i][j][k] = 0;
			}
		}
	}

	// derivative wrt v[0]
	dZdv[1][2][0] = -1.0;
	dZdv[2][1][0] = 1.0;

	//derivative wrt v[1]
	dZdv[0][2][1] = 1.0;
	dZdv[2][0][1] = -1.0;

	//derivative wrt v[2]
	dZdv[0][1][2] = -1.0;
	dZdv[1][0][2] = 1.0;

}

 //-----------------------------------------------
 //Partial_SquareSkew will compute the partial derivative of the square of the
 //skew matrix corresponding to the exponential map v wrt the component of v

static void Partial_SquareSkew(double v[3], double dZdv[3][3][3])
{

	//derivative wrt v[0]
	dZdv[0][0][0] = 0.0;
	dZdv[0][1][0] = v[1];
	dZdv[0][2][0] = v[2];
	dZdv[1][0][0] = v[1];
	dZdv[1][1][0] = -2.0 * v[0];
	dZdv[1][2][0] = 0.0;
	dZdv[2][0][0] = v[2];
	dZdv[2][1][0] = 0.0;
	dZdv[2][2][0] = -2.0 * v[0];

	//derivative wrt v[1]
	dZdv[0][0][1] = -2.0 * v[1];
	dZdv[0][1][1] = v[0];
	dZdv[0][2][1] = 0.0;
	dZdv[1][0][1] = v[0];
	dZdv[1][1][1] = 0.0;
	dZdv[1][2][1] = v[2];
	dZdv[2][0][1] = 0.0;
	dZdv[2][1][1] = v[2];
	dZdv[2][2][1] = -2.0 * v[1];

	//derivative wrt v[2]
	dZdv[0][0][2] = -2.0 * v[2];
	dZdv[0][1][2] = 0.0;
	dZdv[0][2][2] = v[0];
	dZdv[1][0][2] = 0.0;
	dZdv[1][1][2] = -2.0 * v[2];
	dZdv[1][2][2] = v[1];
	dZdv[2][0][2] = v[0];
	dZdv[2][1][2] = v[1];
	dZdv[2][2][2] = 0.0;

}

static void ConstMult3(double X[3][3][3], double c)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				X[i][j][k] *= c;
			}
		}
	}
}

//------------------------------------------------
//Multiply each element of given vector v by a constanct c

static void Scale(double v[3], double c)
{
	for (int i = 0; i < 3; i++) {
		v[i] *= c;
	}
}

//---------------------------------------------------
//As we parametrize the space of the rotation as a ball in R3 with
// radius 2*PI if the point is out of the ball we map it to a point
// with the same rotation matrix

static void Check_Inside(double v[3])
{

	double theta = Norm(v);

	if (theta > 2.0 * M_PI) {

		double deg = fmod(theta, 2.0 * M_PI);
		Scale(v, deg / theta);
	}
}

//-------------------------------------------------
//Partial_R_Exp will calculate the partial derivative of the rotation matrix corresponding
// to the exponenial map v wrt all components of v

static void Partial_R_Exp(double v[3], double dRdv[3][3][3])
{

	Check_Inside(v);

	double theta = Norm(v);

	// cout << theta << endl;

	if (theta >= Epsilon) {
		double Dz1[3][3][3];

		Partial_Skew(Dz1);

		double c = sin(theta) / theta;

		ConstMult3(Dz1, c);

		double Dz2[3][3][3];

		Partial_SquareSkew(v, Dz2);
		c = (1 - cos(theta)) / (theta * theta);
		ConstMult3(Dz2, c);

		for (int i = 0; i < 3; i++) {
			double Z1[3][3];

			//int j;
			Exp_To_Skew(v, Z1);
			double Add =
			    (v[i] / (theta * theta)) * (cos(theta) -
							sin(theta) / theta);
			ConstMult(Z1, Add);

			double Z2[3][3];

			Add =
			    (v[i] / (theta * theta * theta)) * (sin(theta) -
								2 * (1 -
								     cos(theta))
								/ theta);
			Exp_To_Skew_Square(v, Z2);
			ConstMult(Z2, Add);

			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					dRdv[j][k][i] =
					    Z1[j][k] + Z2[j][k] + Dz1[j][k][i] +
					    Dz2[j][k][i];
				}
			}
		}
	}

	else
	{
		double dXdv[3][3][3];

		Partial_Skew(dXdv);
		ConstMult3(dXdv, 1 - theta * theta / 6);

		double dYdv[3][3][3];

		Partial_SquareSkew(v, dYdv);
		ConstMult3(dYdv, 0.5 - theta * theta / 24);

		for (int i = 0; i < 3; i++) {

			double Z[3][3];

			Exp_To_Skew(v, Z);
			ConstMult(Z, -1.0 / 3.0 * v[i]);

			double S[3][3];

			Exp_To_Skew_Square(v, S);
			ConstMult(S, -v[i] / 12.0);

			double ZS[3][3];

			Sum(Z, S, ZS);

			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					dRdv[j][k][i] =
					    dXdv[j][k][i] + dYdv[j][k][i] +
					    Z[j][k] + S[j][k];
		}
	}
}

static void Mult_Partial(double v[3], double pos[3], double grad[3], double result[3])
{

	double dvdR[3][3][3];
	int i, j, k;

	Partial_R_Exp(v, dvdR);

	double temp[3];

	for (i = 0; i < 3; i++) {
		result[i] = 0;

		for (j = 0; j < 3; j++)
			temp[j] = 0;

		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				temp[j] += dvdR[j][k][i] * pos[k];

		for (j = 0; j < 3; j++)
			result[i] += temp[j] * grad[j];
	}
}

void mol_rigidbody_grad(double *grad, struct atomgrp *ag, double *inp, double *origin)
{
	memset(grad, 0, 6 * sizeof(double));

	for (int j = 0; j < ag->nactives; j++) {
		int i = ag->activelist[j];

		double gr[3];
		gr[0] = -ag->atoms[i].GX;
		gr[1] = -ag->atoms[i].GY;
		gr[2] = -ag->atoms[i].GZ;

		double pos[3];
		pos[0] = origin[3 * j];
		pos[1] = origin[3 * j + 1];
		pos[2] = origin[3 * j + 2];

		double result[3];
		Mult_Partial(inp, pos, gr, result);

		grad[0] += result[0];
		grad[1] += result[1];
		grad[2] += result[2];
		grad[3] += gr[0];
		grad[4] += gr[1];
		grad[5] += gr[2];
	}
}

//move the ligand atoms according to the given transformation 
//change has 6 members, the first 3 is the exponential map for the rotation and the other is for the translation

void rigidbody2ag(double *change, struct atomgrp *ag, struct rigidbody *rigidbody)
{
	double rotate[3][3];
	Exp_To_R(change, rotate);


	for (int j = 0; j < ag->nactives; j++) {
		double pos[3];
		double newpos[3];
		pos[0] = rigidbody->origin[3 * j];
		pos[1] = rigidbody->origin[3 * j + 1];
		pos[2] = rigidbody->origin[3 * j + 2];

		for (int i = 0; i < 3; i++)
			newpos[i] = 0;

		for (int i = 0; i < 3; i++)
			for (int jj = 0; jj < 3; jj++)
				newpos[i] += rotate[i][jj] * pos[jj];

		for (int i = 0; i < 3; i++)
			newpos[i] += change[i + 3];

		int i = ag->activelist[j];

		ag->atoms[i].X = newpos[0] + rigidbody->center[0];
		ag->atoms[i].Y = newpos[1] + rigidbody->center[1];
		ag->atoms[i].Z = newpos[2] + rigidbody->center[2];

	}
}

void ag2rigidbody( struct rigidbody *rigidbody, struct atomgrp *ag)
{
	double X, Y, Z;
	X = Y = Z = 0.0;
	for (int j = 0; j < ag->nactives; j++) {
		int i = ag->activelist[j];
		X += ag->atoms[i].X;
		Y += ag->atoms[i].Y;
		Z += ag->atoms[i].Z;
	}
	X /= ag->nactives;
	Y /= ag->nactives;
	Z /= ag->nactives;
	rigidbody->center[0] = X;
	rigidbody->center[1] = Y;
	rigidbody->center[2] = Z;

	for (int j = 0; j < ag->nactives; j++) {
		int i = ag->activelist[j];
		rigidbody->origin[(3*j)+0] = ag->atoms[i].X - X;
		rigidbody->origin[(3*j)+1] = ag->atoms[i].Y - Y;
		rigidbody->origin[(3*j)+2] = ag->atoms[i].Z - Z;
	}
}
