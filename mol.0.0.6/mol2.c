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

void extract_filename( const char *str, char *fname )
{
   int l = 0;
   
   while ( str[ l ] ) l++;
      
   while ( ( l > 0 ) && isspace( str[ l - 1 ] ) ) l--;
   
   int k = l - 1;
   
   while ( k >= 0 )
     {
       if ( str[ k ] == '/' ) break;
       k--;
     }
     
   int i = 0;
   
   for ( int j = k + 1; j < l; j++ )
     fname[ i++ ] = str[ j ];
     
   fname[ i ] = 0;    
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
           if ( fgets( buffer, linesize, fp ) == NULL )
             {
               print_error( "Failed to read %s!", mol2file );
               return 0;
             }
           
	   /* artem modifications
           char fname1[ linesize + 1 ], fname2[ linesize + 1 ];
           
           extract_filename( pdbfile, fname1 );
           extract_filename( buffer, fname2 );           
           
           if ( strcmp( fname1, fname2 ) )
             {
               print_error( "%s was generated from %s, not %s!", mol2file, fname2, fname1 );
               return 0;
             }
	   */
             
           if ( fgets( buffer, linesize, fp ) == NULL )
             {
               print_error( "Failed to read %s!", mol2file );
               return 0;
             } 
             
           int nat = atoi( buffer );                         
           
           if ( nat != ag->natoms )
             {
               print_error( "%s has %d atoms instead of %d atoms!", mol2file, nat, ag->natoms );
               return 0;
             }       
             
           mol_info_found = 1;       
         }

       if ( strstr( buffer, "<TRIPOS>ATOM" ) != NULL )     
         {
           if ( !mol_info_found )
             {
               print_error( "MOLECULE Record missing in %s!", mol2file );
               return 0;
             }
           
           for ( int i = 0; i < ag->natoms; i++ )
             {
               if ( fgets( buffer, linesize, fp ) == NULL )
                 {
                   print_error( "Failed to read %s!", mol2file );
                   return 0;
                 }
                 
               int atom_id;
               char atom_name[ 20 ], atom_type[ 20 ];
               double x, y, z;
               
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
