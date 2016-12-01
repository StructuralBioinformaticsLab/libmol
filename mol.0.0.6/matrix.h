/*
Copyright (c) 2009-2012, Structural Bioinformatics Laboratory, Boston University
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
#ifndef _MOL_MATRIX_H_
#define _MOL_MATRIX_H_

/** \file matrix.h
	structures and functions for matrix operations
*/

struct matrix2di // 2d matrix of integers
{
	int** vals; // pointer to array of rows
	int ni, nj; // number of rows, number of cols
};

struct matrix2df // 2d matrix of floats
{
	float** vals; // pointer to array of rows
	int ni, nj; // number of rows, number of cols
};

struct mol_matrix3f {
	float m11; float m12; float m13;
	float m21; float m22; float m23;
	float m31; float m32; float m33;
};

struct mol_matrix3 {
	double m11; double m12; double m13;
	double m21; double m22; double m23;
	double m31; double m32; double m33;
};

struct mol_matrix3l {
	long double m11; long double m12; long double m13;
	long double m21; long double m22; long double m23;
	long double m31; long double m32; long double m33;
};

/**
  creates a matrix2di of size ni x nj,
  with vals (vals[0]) allocated as a contiguous block of memory
*/
struct matrix2di* matrix2di_create (int ni, int nj);
struct matrix2df* matrix2df_create (int ni, int nj);

/**
  frees the matrix A and all of its allocated memory
*/
void matrix2di_destroy (struct matrix2di* A);
void matrix2df_destroy (struct matrix2df* A);

/**
  initializes the vals of A to initval
*/
void matrix2di_init (struct matrix2di* A, int initval);
void matrix2df_init (struct matrix2df* A, float initval);

/**
  print the matrix A to stdout
*/
void matrix2di_print (struct matrix2di* A);
void matrix2df_print (struct matrix2df* A);

/**
  fprint the matrix A to the file whose name is the string pointed to by path
*/
void matrix2di_fprint (struct matrix2di* A, const char* path);
void matrix2df_fprint (struct matrix2df* A, const char* path);

/**
  reads the matrix in the file at path
*/
struct matrix2di* matrix2di_read (const char* path);
struct matrix2df* matrix2df_read (const char* path);

/**
  compares ni of A and B, and nj of A and B
  \return 0 if the sizes are equal, 1 if the sizes differ
*/
int matrix2di_size_diff (struct matrix2di* A, struct matrix2di* B);
int matrix2df_size_diff (struct matrix2df* A, struct matrix2df* B);

/**
  perform operation on matrices A and B and store the result in C
*/
void matrix2di_add (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_sub (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_pairwise_mult (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_pairwise_div (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_pairwise_arith (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C, int op);

void matrix2df_add (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_sub (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_pairwise_mult (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_pairwise_div (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_pairwise_arith (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C, int op);

void matrix2df_log (struct matrix2df* A, struct matrix2df* B);
float matrix2df_element_sum (struct matrix2df* A);
void matrix2df_scalar_mult (struct matrix2df* A, float val, struct matrix2df* B);

void print_mol_matrix3f (struct mol_matrix3f* rmatrix);

#endif
