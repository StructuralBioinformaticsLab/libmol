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
#ifndef _MOL_ROTAMER_H_
#define _MOL_ROTAMER_H_

/**
 * Gets the rotation matrix and translation vector for moving
 * a rotamer from the canonical position into the
 * backbone of residue in atomgrp
 * \param[in]  ag            atomgrp containing residue
 * \param[in]  residue_index index inside iares of residue
 * \param[out] rot_to_bb     matrix to rotate residue to
 *                           backbone
 * \param[out] trans_to_bb   vector to translate residue to
 *                           backbone
 */
void get_detransform(const struct atomgrp *ag, const int residue_index,
		     struct mol_matrix3f *rot_to_bb, struct mol_vector3f *trans_to_bb);

#endif
