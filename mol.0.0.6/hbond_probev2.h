#ifndef _MOL_HBONDPV2_H_
#define _MOL_HBONDPV2_H_
void hbondeng_probe_weighted(struct atomgrp *ag, double *energy, struct nblist *nblst);
void hbondeng_split_probe_weighted(struct atomgrp *ag, double *energy, struct nblist *nblst, int atom_split);
#endif
