#ifndef _MOL_MOL2_H_
#define _MOL_MOL2_H_

/**
  routine for reading hybridization states from CHIMERA MOL2 file 
*/
enum Hybridization_State mol_hydridization_from_sybyl(const char *sybyl_type);

int read_hybridization_states_from_mol2( const char* mol2file, struct atomgrp* ag );

#endif
