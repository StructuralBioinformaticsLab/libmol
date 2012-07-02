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
#ifndef _MOL_SUBAG_H_
#define _MOL_SUBAG_H_

/**
 * Replaces atoms in a larger atomgrp with the atoms
 * from another atomgrp.  The first item in list should
 * correspond to the first atom in new_actives, and so on.
 * @param[out] ag          The atomgrp that has its atoms replaced
 * @param[in]  new_actives The source of atoms to be placed
 * @param[in]  nlist       The number of atoms to be placed
 * @param[in]  list        The indices of the replaced atoms in ag
 */
void setup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist, int *list);

/**
 * The inverse of seutp_subag.  Takes atoms from a larger atomgrp
 * and places them in another atomgrp.  The first item in list should
 * correspond to the first atom in new_actives, and so on.
 * @param[in]  ag          The source of atoms to be placed
 * @param[out] new_actives The atomgrp that hs its atoms replaced
 * @param[in]  nlist       The number of atoms to be placed
 * @param[in]  list        The indices of the atoms in ag to be copied to new_actives
 */
void unsetup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist, int *list);

#endif
