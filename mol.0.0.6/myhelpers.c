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
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#ifdef _WIN32
#include "../mol.0.0.6.h"
#else
#include _MOL_INCLUDE_
#endif

#ifndef _DEBUG_
void* mymalloc (size_t size)
{
//	_PRINT_DEPRECATED_
//	fprintf (stderr, "\t(please write your own malloc wrapper)\n");

	void* v = (void*) malloc (size);
	if (v == NULL)
	{
		exit (EXIT_FAILURE);
	}
	return v;
}

void* mycalloc (size_t nmemb, size_t size)
{
//	_PRINT_DEPRECATED_
//	fprintf (stderr, "\t(please write your own calloc wrapper)\n");

	void* v = (void*) calloc (nmemb, size);
	if (v == NULL)
	{
		perror ("calloc"), exit (EXIT_FAILURE);
	}
	return v;
}

void* myrealloc (void* ptr, size_t size)
{
//	_PRINT_DEPRECATED_
//	fprintf (stderr, "\t(please write your own realloc wrapper)\n");

	void* v;
	if (size < 1)
	{
		fprintf (stderr, "warning: _mol_realloc called with size 0, no realloc will occur\n");
		return ptr;
		
	}
	v = (void*) realloc (ptr, size);
	if (v == NULL && ptr != NULL)
	{
		perror ("realloc"), exit (EXIT_FAILURE);
	}
	return v;
}
#endif /* _DEBUG_ */


#if (defined _DARWIN_  && defined _DARWIN_SNOW_LEOPARD_) || defined(_WIN32)

int
getline2 (char **lineptr, size_t *n, FILE *stream)
{
 char* read=(char*)malloc(10000*sizeof(char));
 char* result= fgets(read, 10000, stream);
 *lineptr=read;
 *n=strlen(read);
 if (result==NULL) return -1;
 return 0;
}

#endif /* _DARWIN_ && _DARWIN_SNOW_LEOPARD_ */


void myexit (int status)
{
	exit (status);
}

FILE* myfopen (const char* path, const char* mode)
{
	FILE* fp;
    if ( strncmp(path, "-", 2) == 0 ) {
        if (strncmp(mode, "r", 1) == 0 ) {
            return stdin;
        } else if (strncmp(mode, "w", 1) == 0) {
            return stdout;
        } else {
            fprintf(stderr, "case in myfopen not handled\n");
            exit( EXIT_FAILURE );
        }
    }
	fp = fopen (path, mode);
	if (fp == NULL) // file could not be opened
	{
		fprintf (stderr, "fopen error on %s:\n", path);
		perror ("fopen"), exit (EXIT_FAILURE);
	}

	return fp;
}

void myfclose (FILE* fp)
{
	int retval;
    if (fp == stdin || fp == stdout) {
        rewind(fp); //Ryan thinks this is a hack, David thinks it's beautiful
        return;     //if we are using a pdb from stdin, we often need to reuse
    }               //it, so we just rewind to the beginning

	retval = fclose (fp);

	if (retval != 0) // file could not be closed
	{
		fprintf (stderr, "fclose error:\n");
		perror ("fclose"), exit (EXIT_FAILURE);
	}
}

int iswhiteline (char* line)
{
	int i;
	for (i = 0; line[i] != '\0'; i++)
	{
		int ch = line[i];
		if (ch != ' ' && ch != '\t' && ch != '\n')
		{
			return 0;
		}
	}
	return 1;
}

char* rstrip (char* string)
{
    size_t length = strlen(string);

    char* end = string + length - 1;

    while (end >= string && isspace(*end))
        end--;
    *(end+1) = '\0';

    return string;
}

/**
	Prints a formatted error message.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]	
*/
void print_error( char *format, ... )
{
   char eMsg[ 500 ];
   va_list args;
   
   va_start( args, format );
   
   vsprintf( eMsg, format, args );
   
   va_end( args );
   
   printf( "\nError: %s\n\n", eMsg );   
   
   fflush( stdout );   
}
