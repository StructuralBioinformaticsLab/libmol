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
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include _MOL_INCLUDE_

#ifndef linesize
   #undef linesize
#endif   

#define linesize 128

enum Hybridization_State mol_hydridization_from_sybyl(const char *sybyl_type)
{
	if (sybyl_type == NULL)
		return UNKNOWN_HYBRID;

	/* Carbons */
	if (strcmp("C.3", sybyl_type) == 0) /* carbon sp3 */
		return SP3_HYBRID;
	if (strcmp("C.2", sybyl_type) == 0) /* carbon sp2 */
		return SP2_HYBRID;
	if (strcmp("C.cat", sybyl_type) == 0) /* carbocation (C+) used only in a guadinium group */
		return SP2_HYBRID;
	if (strcmp("C.1", sybyl_type) == 0) /* carbon sp */
		return SP1_HYBRID;
	if (strcmp("C.ar", sybyl_type) == 0) /* carbon aromatic */
		return RING_HYBRID;

	/* Nitrogens */
	if (strcmp("N.3", sybyl_type) == 0) /* nitrogen sp3 */
		return SP3_HYBRID;
	if (strcmp("N.4", sybyl_type) == 0) /* nitrogen sp3 positively charged */
		return SP3_HYBRID;
	if (strcmp("N.2", sybyl_type) == 0) /* nitrogen sp2 */
		return SP2_HYBRID;
	if (strcmp("N.am", sybyl_type) == 0) /* nitrogen amide */
		return SP2_HYBRID;
	if (strcmp("N.pl3", sybyl_type) == 0) /* nitrogen trigonal planar */
		return SP2_HYBRID;
	if (strcmp("N.1", sybyl_type) == 0) /* nitrogen sp */
		return SP1_HYBRID;
	if (strcmp("N.ar", sybyl_type) == 0) /* nitrogen aromatic */
		return RING_HYBRID;

	/* Oxygens */
	if (strcmp("O.3", sybyl_type) == 0) /* oxygen sp3 */
		return SP3_HYBRID;
	if (strcmp("O.spc", sybyl_type) == 0) /* oxygen in single point charge (SPC) water model */
		return SP3_HYBRID;
	if (strcmp("O.t3p", sybyl_type) == 0) /* oxygen in Transferable Intermolecular Potential (TIP3P) water model */
		return SP3_HYBRID;
	if (strcmp("O.2", sybyl_type) == 0) /* oxygen sp2 */
		return SP2_HYBRID;
	if (strcmp("O.co2", sybyl_type) == 0) /* oxygen in carboxylate and phosphate groups */
		return SP2_HYBRID;

	/* Sulfurs */
	if (strcmp("S.3", sybyl_type) == 0) /* sulfur sp3 */
		return SP3_HYBRID;
	if (strcmp("S.O", sybyl_type) == 0) /* sulfoxide sulfur */
		return SP3_HYBRID;
	if (strcmp("S.O2", sybyl_type) == 0) /* sulfone sulfur */
		return SP3_HYBRID;
	if (strcmp("S.2", sybyl_type) == 0) /* sulfur sp2 */
		return SP2_HYBRID;

	/* Phosphorous */
	if (strcmp("P.3", sybyl_type) == 0) /* phosphorous sp3 */
		return SP3_HYBRID;

	/* everything else */
	return UNKNOWN_HYBRID;
}


int read_hybridization_states_from_mol2( const char* mol2file, struct atomgrp* ag )
{
   char buffer[ linesize + 1 ];
   int mol_info_found = 0;

   FILE* fp = myfopen( mol2file, "r" );
   
   if ( fp == NULL )
     {
       print_error( "Failed to open %s!", mol2file );
       return 0;
     }
     
   while ( fgets( buffer, linesize, fp ) != NULL )
     {
       if ( strstr( buffer, "<TRIPOS>MOLECULE" ) != NULL )     
         {
           int nat;
           if ( fgets( buffer, linesize, fp ) == NULL )
             {
               print_error( "Failed to read %s!", mol2file );
               return 0;
             }
           
           if ( fgets( buffer, linesize, fp ) == NULL )
             {
               print_error( "Failed to read %s!", mol2file );
               return 0;
             } 
             
           nat = atoi( buffer );                         
           
           if ( nat != ag->natoms )
             {
               print_error( "%s has %d atoms instead of %d atoms!", mol2file, nat, ag->natoms );
               return 0;
             }       
             
           mol_info_found = 1;       
         }

       if ( strstr( buffer, "<TRIPOS>ATOM" ) != NULL )     
         {
           int i;
           if ( !mol_info_found )
             {
               print_error( "MOLECULE Record missing in %s!", mol2file );
               return 0;
             }
           
           for ( i = 0; i < ag->natoms; i++ )
             {
               int atom_id;
               char atom_name[ 20 ], atom_type[ 20 ];
               double x, y, z;
               
               if ( fgets( buffer, linesize, fp ) == NULL )
                 {
                   print_error( "Failed to read %s!", mol2file );
                   return 0;
                 }
                 
               if ( sscanf( buffer, "%d %s %lf %lf %lf %s", &atom_id, atom_name, &x, &y, &z, atom_type ) != 6 )
                 {
                   print_error( "Failed to read ATOM record from %s!", mol2file );
                   return 0;                 
                 }  

	       /* artem's modifications
               if ( ( x != ag->atoms[ i ].X ) || ( y != ag->atoms[ i ].Y ) || ( z != ag->atoms[ i ].Z ) )  
                 {
                   print_error( "Coordinates of atom %d in %s do not match with those in %s!", atom_id, mol2file, pdbfile );
                   return 0;                                  
                 }
	       */
                 
               ag->atoms[ i ].hybridization = UNKNOWN_HYBRID;
                              
               if ( atom_type[ 1 ] == '.' )
                 {
                   switch ( atom_type[ 2 ] )
                     {
                       case '4' : 
                       case '3' : ag->atoms[ i ].hybridization = SP3_HYBRID; break;
                       
                       case '2' : ag->atoms[ i ].hybridization = SP2_HYBRID; break;
                       
                       case '1' : ag->atoms[ i ].hybridization = SP1_HYBRID; break;
                       
                       case 'a' : if ( atom_type[ 3 ] == 'r' ) ag->atoms[ i ].hybridization = RING_HYBRID; 
                                  break;
                                  
                       case 'c' : if ( ( atom_type[ 3 ] == 'o' ) && ( atom_type[ 4 ] == '2' ) ) ag->atoms[ i ].hybridization = SP2_HYBRID; 
                                  break;                                  
                     } 
                 }
                 
//               printf( "%d %s ( %d, %d )\n", i + 1, atom_type, ag->atoms[ i ].hprop, ag->atoms[ i ].hybridization );                 
             }                
         }
     }
   
   myfclose( fp );
   
   return 1;
}
