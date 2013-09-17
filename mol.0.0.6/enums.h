/*
Copyright (c) 2009-2012, Structural Bioinformatics Laboratory, Boston University
Copyright (c) 2013, Acpharis
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
#ifndef _MOL_ENUMS_H_
#define _MOL_ENUMS_H_

enum mol_yeti { MOL_YETI_NONE, MOL_YETI_CARBONYL, MOL_YETI_HYDROXYL,
                MOL_YETI_SULFONAMIDE, MOL_YETI_N5_AROMATIC,
                MOL_YETI_N6_AROMATIC };

enum HBondProp { UNKNOWN_HPROP = 0, HBOND_DONOR = 1, HBOND_ACCEPTOR = 2,
                 DONATABLE_HYDROGEN = 4, ROTATABLE_HYDROGEN = 8,
                 FIXED_HYDROGEN = 16 };

enum Hybridization_State { UNKNOWN_HYBRID = 0, SP1_HYBRID, SP2_HYBRID,
                           SP3_HYBRID, RING_HYBRID };
enum HB_Donor_Type { hbdon_NO = 0, hbdon_BB, hbdon_SC, hbdon_SM, hbdon_MAX };
enum HB_Acceptor_Type { hbacc_NO = 0, hbacc_BB, hbacc_SP2,  hbacc_SP3,
                        hbacc_RING, hbacc_SP2_SM,  hbacc_SP3_SM, hbacc_RING_SM,
                        hbacc_MAX };

enum HB_Type_Evaluation { hbe_NONE = 0, hbe_BB, hbe_BBTURN, hbe_BBHELIX,
                          hbe_BBOTHER, hbe_SP2B, hbe_SP3B, hbe_RINGB,
                          hbe_BSC, hbe_SP2SC, hbe_SP3SC, hbe_RINGSC,
                          /* Small Molecule Section -> */
                          hbe_BSM, hbe_SP2SCSM, hbe_SP3SCSM, hbe_RINGSCSM,
                          hbe_SP2SMSC, hbe_SP3SMSC, hbe_RINGSMSC, hbe_SP2SMB,
                          hbe_SP3SMB, hbe_RINGSMB,hbe_SP2SM, hbe_SP3SM,
                          hbe_RINGSM, hbe_MAX };
enum HB_Weight_Type { hbw_NONE = -1,  hbw_TOTAL, hbw_SR_BB, hbw_LR_BB,
                       hbw_BB_SC, hbw_SC, hbw_BB_SM, hbw_SC_SM, hbw_SM, hbw_H2O};

#endif
