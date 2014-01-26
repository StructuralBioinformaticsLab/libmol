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
#include <string.h>

#ifdef _WIN32
#define _MOL_VERSION_ "0.0.6"
#define _GIT_VERSION_ "Windows Unknown"
#endif
#include _MOL_INCLUDE_

void
mol_version (char** lineptr, size_t *n)
{
    size_t length = strlen(_MOL_VERSION_);

    if (*lineptr == NULL) {
        *lineptr = _mol_malloc( sizeof(char) * length );
        *n = length;
    } else if( *n < length) {
        *lineptr = _mol_realloc(*lineptr, sizeof(char) * length );
        *n = length;
    }

	sprintf(*lineptr, "%s", _MOL_VERSION_);
}

void
mol_git_version (char** lineptr, size_t *n)
{
    size_t length = strlen(_GIT_VERSION_);

    if (*lineptr == NULL) {
        *lineptr = _mol_malloc( sizeof(char) * length );
        *n = length;
    } else if( *n < length) {
        *lineptr = _mol_realloc(*lineptr, sizeof(char) * length );
        *n = length;
    }

	sprintf(*lineptr, "%s", _GIT_VERSION_);
}
