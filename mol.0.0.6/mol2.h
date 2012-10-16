#ifndef _MOL_MOL2_H_
#define _MOL_MOL2_H_

/**
  routine for reading hybridization states from CHIMERA MOL2 file 
*/

int read_hybridization_states_from_mol2( const char* mol2file, const char* pdbfile, struct atomgrp* ag );

#endif
