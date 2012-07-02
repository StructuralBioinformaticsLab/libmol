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
#ifndef _MOL_MEM_H_
#define _MOL_MEM_H_

#ifdef _DEBUG_
#define _mol_malloc(size) malloc(size)
#define _mol_realloc(ptr, size) realloc(ptr, size)
#define _mol_calloc(nmemb, size) calloc(nmemb, size)
#endif

/**
	This is a wrapper to malloc that will print
	an error message if malloc returns NULL.
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void*
_mol_malloc (size_t size);
#endif

/**
	This is a wrapper to calloc that will print
	an error message if malloc returns NULL.
	\param nmemb size of array
	\param size size of each element in the array
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void*
_mol_calloc (size_t nmemb, size_t size);
#endif

/**
	This is a wrapper to realloc that will print
	an error message accordingly.
	\param ptr pointer to the original memory
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void*
_mol_realloc (void* ptr, size_t size);
#endif


#endif
