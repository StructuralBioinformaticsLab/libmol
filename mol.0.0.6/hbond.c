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
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#ifdef _WIN32
#include "mol_stdbool.h"
#else
#include <stdbool.h>
#endif

#ifdef _WIN32
#include "../mol.0.0.6.h"
#else
#include _MOL_INCLUDE_
#endif

// POLY_OUT is the polynomial value outside of clip interval [xmin xmax]
#define POLY_OUT 0.0

#define poly2generator(name, xmin, xmax, c1, c0) static void name(FLOAT x, FLOAT *val, FLOAT *der) { if ( (x<=xmin) || (x>=xmax) ) {(*val) = POLY_OUT; (*der) = 0.0; return; } (*val) = (*der) = c1; (*val) = (*val)*x+c0; }

#define poly3generator(name, xmin, xmax, c2, c1, c0) static void name(FLOAT x, FLOAT *val, FLOAT *der) { if ( (x<=xmin) || (x>=xmax) ) {(*val) = POLY_OUT; (*der) = 0.0; return; } (*val) = (*der) = c2; (*val) = (*val)*x+c1; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c0; }

#define poly4generator(name, xmin, xmax, c3, c2, c1, c0) static void name(FLOAT x, FLOAT *val, FLOAT *der) { if ( (x<=xmin) || (x>=xmax) ) {(*val) = POLY_OUT; (*der) = 0.0; return; } (*val) = (*der) = c3; (*val) = (*val)*x+c2; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c1; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c0; }

#define poly5generator(name, xmin, xmax, c4, c3, c2, c1, c0) static void name(FLOAT x, FLOAT *val, FLOAT *der) { if ( (x<=xmin) || (x>=xmax) ) {(*val) = POLY_OUT; (*der) = 0.0; return; } (*val) = (*der) = c4; (*val) = (*val)*x+c3; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c2; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c1; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c0; }

#define poly6generator(name, xmin, xmax, c5, c4, c3, c2, c1, c0) static void name(FLOAT x, FLOAT *val, FLOAT *der) { if ( (x<=xmin) || (x>=xmax) ) {(*val) = POLY_OUT; (*der) = 0.0; return; } (*val) = (*der) = c5; (*val) = (*val)*x+c4; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c3; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c2; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c1; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c0; }

#define poly7generator(name, xmin, xmax, c6, c5, c4, c3, c2, c1, c0) static void name(FLOAT x, FLOAT *val, FLOAT *der) { if ( (x<=xmin) || (x>=xmax) ) {(*val) = POLY_OUT; (*der) = 0.0; return; } (*val) = (*der) = c6; (*val) = (*val)*x+c5; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c4; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c3; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c2; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c1; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c0; }

#define poly8generator(name, xmin, xmax, c7, c6, c5, c4, c3, c2, c1, c0) static void name(FLOAT x, FLOAT *val, FLOAT *der) { if ( (x<=xmin) || (x>=xmax) ) {(*val) = POLY_OUT; (*der) = 0.0; return; } (*val) = (*der) = c7; (*val) = (*val)*x+c6; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c5; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c4; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c3; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c2; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c1; (*der) = (*der)*x+(*val); (*val) = (*val)*x+c0; }

// Defining polynomials for hbond score calculation

#ifdef NON_POSITIVE_POLY

poly8generator(AH_BBHelix, 1.78218633, 2.6757, 12.93768086, -221.0155722,
	       1604.391304, -6409.335773, 15200.86425, -21375.00216,
	       16475.98811, -5361.55644)
poly8generator(AH_BBOther, 1.6971523, 2.679339, 13.58980244, -224.0452428,
	       1568.933094, -6044.257847, 13820.1498, -18730.96076, 13912.92238,
	       -4361.995425)
poly5generator(AH_SP2, 1.6941, 2.5, 10.98727738, -100.2401419, 340.9733405,
	       -511.6111233, 285.0061262)
poly5generator(AH_SP3, 1.755, 2.521385, 7.011735538, -68.99968829, 251.820931,
	       -403.3593133, 238.7378958)
poly8generator(xD_BBHelix, 0.3746, 1.04, 223.5268153, -757.7254095, 1019.593508,
	       -689.2232431, 240.1436064, -37.84119583, 0.85868904, 0.278181985)
poly8generator(xD_BBOther, 0.76, 1.09, 111.9877946, -380.3066184, 514.7650204,
	       -352.4092342, 124.6219703, -19.94401946, 0.149314979,
	       0.635771774)
poly3generator(xD_SP2short, 0.7071, 1.01, -0.562582503, -0.746682668,
	       0.809265171)
poly3generator(xD_SP2long, 0.0, 1.01, 0.094962885, -0.254313172, 0.0)
poly3generator(xD_SP3short, 0.61566, 1.01, -0.100140144, -1.139139041,
	       0.739279186)
poly3generator(xD_SP3long, 0.0, 1.01, 0.089380221, -0.207503776, 0.0)
poly8generator(xH_BBHelix, 0.156, 1.03, 54.80664331, -196.8196655, 295.9418886,
	       -232.105602, 96.99124565, -20.60918361, 1.573169816, 0.000745458)
poly8generator(xH_BBOther, 0.61566, 1.07, 43.94483847, -144.3836033,
	       193.5865176, -132.4469355, 47.28137288, -8.945888012, -0.227035135,
	       0.791902995)
poly3generator(xH_SP2short, 0.0, 1.08, 1.720984644, -1.855254573, 0.0)
poly3generator(xH_SP2long, 0.0, 1.01, 0.439598249, -0.444673076, 0.0)
poly3generator(xH_SP3, 0.0, 1.06, 1.761487842, -1.876959406, 0.0)
poly7generator(xH_Ring, 0.7608, 1.089, 37.744316, -117.731674, 143.0759275,
	       -86.2258835, 26.7448175, -4.4699705, 0.6458455)
#else

poly8generator(AH_BBHelix, MIN_AH, 2.8, 12.93768086, -221.0155722, 1604.391304,
	       -6409.335773, 15200.86425, -21375.00216, 16475.98811,
	       -5361.55644)
poly8generator(AH_BBOther, MIN_AH, 2.745, 13.58980244, -224.0452428,
	       1568.933094, -6044.257847, 13820.1498, -18730.96076, 13912.92238,
	       -4361.995425)
poly5generator(AH_SP2, MIN_AH, 2.5, 10.98727738, -100.2401419, 340.9733405,
	       -511.6111233, 285.0061262)
poly5generator(AH_SP3, MIN_AH, 2.5, 7.011735538, -68.99968829, 251.820931,
	       -403.3593133, 238.7378958)
poly8generator(xD_BBHelix, MIN_xD, MAX_xD, 223.5268153, -757.7254095,
	       1019.593508, -689.2232431, 240.1436064, -37.84119583, 0.85868904,
	       0.278181985)
poly8generator(xD_BBOther, MIN_xD, MAX_xD, 111.9877946, -380.3066184,
	       514.7650204, -352.4092342, 124.6219703, -19.94401946, 0.149314979,
	       0.635771774)
poly3generator(xD_SP2short, MIN_xD, MAX_xD, -0.562582503, -0.746682668,
	       0.809265171)
poly3generator(xD_SP2long, MIN_xD, MAX_xD, 0.094962885, -0.254313172, 0.0)
poly3generator(xD_SP3short, MIN_xD, MAX_xD, -0.100140144, -1.139139041,
	       0.739279186)
poly3generator(xD_SP3long, MIN_xD, MAX_xD, 0.089380221, -0.207503776, 0.0)
poly8generator(xH_BBHelix, MIN_xH, MAX_xH, 54.80664331, -196.8196655,
	       295.9418886, -232.105602, 96.99124565, -20.60918361, 1.573169816,
	       0.000745458)
poly8generator(xH_BBOther, MIN_xH, MAX_xH, 43.94483847, -144.3836033,
	       193.5865176, -132.4469355, 47.28137288, -8.945888012, -0.227035135,
	       0.791902995)
poly3generator(xH_SP2short, MIN_xH, MAX_xH, 1.720984644, -1.855254573, 0.0)
poly3generator(xH_SP2long, MIN_xH, MAX_xH, 0.439598249, -0.444673076, 0.0)
poly3generator(xH_SP3, MIN_xH, MAX_xH, 1.761487842, -1.876959406, 0.0)
poly7generator(xH_Ring, MIN_xH, MAX_xH, 37.744316, -117.731674, 143.0759275,
	       -86.2258835, 26.7448175, -4.4699705, 0.6458455)
#endif
#if FADING_FUNCTION == LINEAR_FADE
#define create_fade_interval( name, min0, fmin, fmax, max0 ) \
static void name##_value_deriv(FLOAT x, FLOAT *val, FLOAT *der) { \
     const FLOAT name##_dfade_min = ( fmin == min0 ) ? 0.0 : ( 1.0 / ( fmin - min0 ) ), name##_dfade_max = ( max0 == fmax ) ? 0.0 : ( 1.0 / ( max0 - fmax ) ); \
     if ( x < fmax ) { if ( x >= fmin ) { *val = 1.0; *der = 0.0; return; } \
                       if ( x <= min0 ) { *val = *der = 0.0; return; } \
                       *der = name##_dfade_min; *val = ( x - min0 ) * name##_dfade_min; \
                     } \
     else { if ( x >= max0 ) { *val = *der = 0.0; return; } \
 	    *der = -name##_dfade_max; *val = ( max0 - x ) * name##_dfade_max; } \
}
// For adjusting xD and xH
create_fade_interval(fade_rBB, MIN_AH, MIN_AH, INTPOL_EDGE_DIST, MAX_AH)
create_fade_interval(fade_rshort, MIN_AH, MIN_AH, INTPOL_MIN, INTPOL_MAX)
create_fade_interval(fade_rlong, INTPOL_MIN, INTPOL_MAX, INTPOL_MAX, MAX_AH)
// fading theta(xD)
create_fade_interval(fade_xD, MIN_xD, INTPOL_EDGE_ANGL, 1, 1)
// fading r,xD
create_fade_interval(fade_xH, MIN_xH, INTPOL_EDGE_ANGL, 1, 1)
#endif //LINEAR_FADE
#ifdef USE_LONG_DOUBLE
#define create_log_fade_interval( name, a, b, c, d, e, f ) \
  static void name##_log_value_deriv( FLOAT x, FLOAT *val, FLOAT *der ) { \
       FLOAT lnx = logl( x ); \
       *val = ( ( ( ( f * lnx + e ) * lnx + d ) * lnx + c ) * lnx + b ) * lnx + a; \
       *der = ( ( ( ( 5 * f * lnx + 4 * e ) * lnx + 3 * d ) * lnx + 2 * c ) * lnx + b ) / x; \
  }
#else
#define create_log_fade_interval( name, a, b, c, d, e, f ) \
  static void name##_log_value_deriv( FLOAT x, FLOAT *val, FLOAT *der ) { \
       FLOAT lnx = log( x ); \
       *val = ( ( ( ( f * lnx + e ) * lnx + d ) * lnx + c ) * lnx + b ) * lnx + a; \
       *der = ( ( ( ( 5 * f * lnx + 4 * e ) * lnx + 3 * d ) * lnx + 2 * c ) * lnx + b ) / x; \
  }
#endif
// used to adjust xD,xH
#if FADING_FUNCTION == LOG_FADE
create_log_fade_interval(fade_rBB, -37.3180547817287, 285.741141861363,
			 -834.307606917721, 1189.6798830102, -826.409671561707,
			 222.919001825369)
create_log_fade_interval(fade_rshort, 5.17967058862565, -73.0669137757406,
			 315.57016260663, -559.689100793323, 438.586611628675,
			 -126.548509998553)
create_log_fade_interval(fade_rlong, -14.0813542571684, 130.05181593123,
			 -453.153559518528, 740.025843886255, -562.47374254511,
			 160.168566150279)
create_log_fade_interval(fade_xD, 0.133459218358311, 13.0816902916604, -57.2491258523768, 86.8476548565786, -3.78050184877689, -58.7848249871094)	//      theta(xD) should fade r,xH sooner
create_log_fade_interval(fade_xH, 0.133459218358311, 13.0816902916604, -57.2491258523768, 86.8476548565786, -3.78050184877689, -58.7848249871094)	// fading r and xD
#endif //LOG_FADE
// The following cubic splines were generated using DataFit 9
// fade_interval: 1.4, 1.4, 2.1, 3.0
int fade_rBB_n = 44;
FLOAT fade_rBB_x[] = { 0, 0.100, 0.500, 0.750, 1.000, 1.200,
                          1.250, 1.280, 1.300, 1.320, 1.350,
                          1.370, 1.390, 1.393, 1.395, 1.398,
                          1.400, 1.405, 1.410, 1.415, 1.420,
                          1.425, 1.430, 1.440, 1.460, 1.480,
                          1.500, 1.510, 1.530, 1.550, 1.600,
                          1.700, 1.800, 1.900, 2.000, 2.100,
                          2.325, 2.550, 2.775, 3.000, 3.100,
                          3.200, 3.300, 3.400, 3.500
};

FLOAT fade_rBB_y[] = { 0, 0.00, 0.00, 0.00, 0.00, 0.00,
                          0.00, 0.00, 0.00, 0.00, 0.00,
                          0.00, 0.00, 0.25, 0.50, 0.75,
                          1.00, 1.00, 1.00, 1.00, 1.00,
                          1.00, 1.00, 1.00, 1.00, 1.00,
                          1.00, 1.00, 1.00, 1.00, 1.00,
                          1.00, 1.00, 1.00, 1.00, 1.00,
                          0.75, 0.50, 0.25, 0.00, 0.00,
                          0.00, 0.00, 0.00, 0.00
};

FLOAT fade_rBB_y2[] =
    { 0, 0.000000000000, 0.002691890564, -0.013997830935, 0.053299433177,
         -0.222350160627, 2.010303873564, -10.351037057964, 48.739729479473,
         -184.607880859929, 582.866449880115, -2637.420428110680,
         9966.815262562670, 31424.968828103199, -47075.067034355103,
         52633.577562452803, -67555.287260737794, 18101.373305087502,
         -4850.205959610880, 1299.450533356390, -347.596173814588,
         90.934161901984, -16.140473793343, 2.954340429037, -0.792784390439,
         0.216797132720, -0.074404140441, 0.012830577206, -0.001289661396,
         -0.007671931621, 0.021997273096, -0.062155853479, 0.226626140818,
         -0.844348709794, 3.150768698359, -11.758726083643, 2.939681857178,
         -0.000001345070, -2.939676476899, 11.758707252664, -3.150658402628,
         0.843926357847, -0.225047028759, 0.056261757190, 0.000000000000
};

// fade_interval: 1.4, 1.4, 1.9, 2.3
int fade_rshort_n = 49;
FLOAT fade_rshort_x[] = { 0, 0.100, 0.500, 0.750, 1.000, 1.200,
                             1.250, 1.280, 1.300, 1.320, 1.350,
                             1.360, 1.370, 1.380, 1.385, 1.390,
                             1.393, 1.395, 1.398, 1.400, 1.405,
                             1.410, 1.415, 1.420, 1.425, 1.430,
                             1.440, 1.460, 1.480, 1.500, 1.510,
                             1.530, 1.550, 1.600, 1.700, 1.800,
                             1.850, 1.900, 2.000, 2.100, 2.200,
                             2.300, 2.350, 2.400, 2.500, 2.600,
                             2.700, 2.800, 2.900, 3.000
};

FLOAT fade_rshort_y[] = { 0, 0.00, 0.00, 0.00, 0.00, 0.00,
                             0.00, 0.00, 0.00, 0.00, 0.00,
                             0.00, 0.00, 0.00, 0.00, 0.00,
                             0.25, 0.50, 0.75, 1.00, 1.00,
                             1.00, 1.00, 1.00, 1.00, 1.00,
                             1.00, 1.00, 1.00, 1.00, 1.00,
                             1.00, 1.00, 1.00, 1.00, 1.00,
                             1.00, 1.00, 0.75, 0.50, 0.25,
                             0.00, 0.00, 0.00, 0.00, 0.00,
                             0.00, 0.00, 0.00, 0.00
};

FLOAT fade_rshort_y2[] =
    { 0, 0.000000000000, -0.000062723877, 0.000326164163, -0.001241932774,
         0.005180992280, -0.046842191700, 0.241190035268, -1.135686888788,
         4.301557519885, -13.581400473759, 95.746531230416, -369.404724447903,
         1381.872366561210, -7552.424750471250, 28827.826635324400,
         25505.632529047500, -45769.902598219698, 52229.253641372103,
         -67491.414309537300, 18084.258610158398, -4845.620131095010,
         1298.221914221990, -347.267525792860, 90.848188949475,
         -16.125230005034, 2.951595540365, -0.792171618577, 0.217090933945,
         -0.076192117203, 0.022970835326, -0.030816447377, 0.100294954181,
         -0.268499292757, 0.755350401180, -2.752902311963, 15.006713069415,
         -57.273949965699, 14.318493362390, -0.000023483861, -14.318399426946,
         57.273621191643, -15.004928295967, 2.746091992223, -0.735811828685,
         0.197155322519, -0.052809461389, 0.014082523037, -0.003520630759,
         0.000000000000
};

// fade_interval: 1.9, 2.3, 2.3, 3.0
int fade_rlong_n = 23;
FLOAT fade_rlong_x[] = { 0, 1.500, 1.600, 1.700, 1.800, 1.850,
                            1.900, 2.000, 2.100, 2.200, 2.300,
                            2.400, 2.500, 2.600, 2.700, 2.800,
                            2.900, 3.000, 3.050, 3.100, 3.200,
                            3.300, 3.400, 3.500
};

FLOAT fade_rlong_y[] = { 0, 0.00, 0.00, 0.00, 0.00, 0.00,
                            0.00, 0.25, 0.50, 0.75, 1.00,
                            0.85, 0.71, 0.56, 0.42, 0.28,
                            0.14, 0.00, 0.00, 0.00, 0.00,
                            0.00, 0.00, 0.00
};

FLOAT fade_rlong_y2[] =
    { 0, 0.000000000000, 0.182712031224, -0.730848124897, 2.740680468365,
         -14.982386560393, 57.188865773208, -14.075404039428, -0.887249615495,
         17.624402501407, -69.610360390132, 20.817039059120, -7.657795846350,
         3.814144326278, -1.598781458762, 2.580981508772, -8.725144576323,
         32.319596796522, -8.467291626486, 1.549569709422, -0.415063315024,
         0.110683550673, -0.027670887668, 0.000000000000
};

// fade_interval: 0, 0.05, 1, 1
int fade_xD_n = 37;
FLOAT fade_xD_x[] = { 0, -0.100, 0.000, 0.010, 0.020, 0.030,
                          0.040, 0.050, 0.060, 0.070, 0.100,
                          0.115, 0.130, 0.150, 0.200, 0.300,
                          0.400, 0.500, 0.600, 0.700, 0.800,
                          0.850, 0.870, 0.885, 0.892, 0.900,
                          0.930, 0.950, 0.960, 0.970, 0.980,
                          0.990, 0.992, 0.994, 0.996, 0.998,
                          1.000, 1.100
};

FLOAT fade_xD_y[] = { 0, 0.00, 0.00, 0.20, 0.40, 0.60,
                         0.80, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 0.80, 0.60, 0.40, 0.20,
                         0.00, 0.00
};

FLOAT fade_xD_y2[] =
    { 0, 0.000000000000, 552.940610336330, -164.693427399256, 105.833099260692,
         -258.638969643515, 928.722779313369, -3456.252147609970,
         896.285811126504, -128.891096896051, 44.947654680636, -11.903734291710,
         2.667282486206, -0.407687982937, 0.074613357741, -0.019996081755,
         0.005370969278, -0.001487795357, 0.000580212150, -0.000833053243,
         0.002752000823, -0.014845898451, 0.097041287097, -0.433064808519,
         2.514176038338, -9.049228436313, 22.254265095103, -97.697482821044,
         541.676366736059, -2069.007984123190, 7734.355569756700,
         -28868.414294903600, 7749.193690060420, -2128.360465338230,
         764.248171292678, -928.632219832562, 2950.280708037570, 0.000000000000
};

// fade_interval: 0, 0.05, 1, 1
int fade_xH_n = 37;
FLOAT fade_xH_x[] = { 0, -0.100, 0.000, 0.010, 0.020, 0.030,
                          0.040, 0.050, 0.060, 0.070, 0.100,
                          0.115, 0.130, 0.150, 0.200, 0.300,
                          0.400, 0.500, 0.600, 0.700, 0.800,
                          0.850, 0.870, 0.885, 0.892, 0.900,
                          0.930, 0.950, 0.960, 0.970, 0.980,
                          0.990, 0.992, 0.994, 0.996, 0.998,
                          1.000, 1.100
};

FLOAT fade_xH_y[] = { 0, 0.00, 0.00, 0.20, 0.40, 0.60,
                         0.80, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 1.00, 1.00, 1.00, 1.00,
                         1.00, 0.80, 0.60, 0.40, 0.20,
                         0.00, 0.00
};

FLOAT fade_xH_y2[] =
    { 0, 0.000000000000, 552.940610336330, -164.693427399256, 105.833099260692,
         -258.638969643515, 928.722779313369, -3456.252147609970,
         896.285811126504, -128.891096896051, 44.947654680636, -11.903734291710,
         2.667282486206, -0.407687982937, 0.074613357741, -0.019996081755,
         0.005370969278, -0.001487795357, 0.000580212150, -0.000833053243,
         0.002752000823, -0.014845898451, 0.097041287097, -0.433064808519,
         2.514176038338, -9.049228436313, 22.254265095103, -97.697482821044,
         541.676366736059, -2069.007984123190, 7734.355569756700,
         -28868.414294903600, 7749.193690060420, -2128.360465338230,
         764.248171292678, -928.632219832562, 2950.280708037570, 0.000000000000
};

int hbeng_HOH_n = 11;
FLOAT hbeng_HOH_x[] = { 0, 2.500, 2.550, 2.600, 2.800, 2.850,
                           2.900, 3.200, 3.300, 3.400, 3.600,
                           3.700
};

FLOAT hbeng_HOH_y[] = { 0, 0.000, 0.000, 0.000, -0.950, -0.975,
                          -1.000, -0.150, -0.050, -0.010, 0.000,
                          0.000
};

FLOAT hbeng_HOH_y2[] =
    { 0, 0.000000000000, 24.162670008471, -96.6506800338837, 93.0860325825914,
         -34.2576056903785, 43.9443901789217, -30.1606428024209,
         -0.548028117397894, -3.64724472798789, 0.71574824266263, 0.000000000000
};

#if FADING_FUNCTION == BSPLINE_FADE
FLOAT eval_spline(FLOAT * x, FLOAT * y, FLOAT * y2, int n, FLOAT v);

#define create_bspline_fade_interval( name ) \
static void name##_bspline_value_deriv( FLOAT x, FLOAT *val, FLOAT *der ) { \
     if ( ( x < name##_x[ 2 ] ) || ( x > name##_x[ name##_n - 1 ] ) ) { *val = *der = 0; } \
     else { \
            FLOAT h = 0.00001; \
            *val = eval_spline( name##_x, name##_y, name##_y2, name##_n, x ); \
            *der = ( eval_spline( name##_x, name##_y, name##_y2, name##_n, x + h ) - ( *val ) ) / h; \
          } \
}

create_bspline_fade_interval(fade_rBB)	// used to adjust xD, xH
create_bspline_fade_interval(fade_rshort)
create_bspline_fade_interval(fade_rlong)
// xD=theta should fade r,xH sooner!
create_bspline_fade_interval(fade_xD)
// fades r,xD
create_bspline_fade_interval(fade_xH)
// water-mediated hbond energy
create_bspline_fade_interval(hbeng_HOH)
// generated using DataFit 9
//
FLOAT eval_spline(FLOAT * x, FLOAT * y, FLOAT * y2, int n, FLOAT v)
{
	int l = 1, h = n;
	FLOAT d;
	FLOAT a;
	FLOAT b;

	while (h - l > 1) {
		int k = (h + l) >> 1;
		if (x[k] > v)
			h = k;
		else
			l = k;
	}

	d = x[h] - x[l];
	a = (x[h] - v) / d;
	b = (v - x[l]) / d;

	return (a * y[l] + b * y[h] +
		(a * (a * a - 1) * y2[l] +
		 b * (b * b - 1) * y2[h]) * (d * d) / 6.0);
}
#endif //BSPLINE_FADE

double dot_product(struct dvector *va, struct dvector *vb)
{
	return ((va->X * vb->X) + (va->Y * vb->Y) + (va->Z * vb->Z));
}

bool is_hydrogen(mol_atom * atom, struct prm * prm)
{
	uint_fast8_t i = 0;
	while (isdigit((prm->atoms[atom->atom_typen].typemin[i])))
		i++;

	if (prm->atoms[atom->atom_typen].typemin[i] == 'H')
		return true;

	return false;
}

// returns the non-self index of the atom bonded to atom
// atomi at atomi's bond index bondi
int bonded_atom_index(mol_atom_group * ag, int atomi, int bondi)
{
	int bonded_atomi;
	mol_atom *atom;
	int bondsi;

	atom = &ag->atoms[atomi];
	assert(bondi < atom->nbondis);
	bondsi = atom->bondis[bondi];
	assert(bondsi < ag->nbonds);

	bonded_atomi = ag->bonds[bondsi].ai;
	if (bonded_atomi == atomi)
		bonded_atomi = ag->bonds[bondsi].aj;

	return bonded_atomi;
}

//#define DEBUG_HBOND_DONOR 1
#define DEBUG_HBOND_DONOR 0
void mark_hbond_donors(struct atomgrp *ag, struct prm *prm)
{
	int atomi;

	for (atomi = 0; atomi < ag->natoms; atomi++) {
		int j;
		int donor_atomi;
		int heavy_donor_base_atomi = 0;
		mol_atom *donor_atom;
		mol_atom *heavy_donor_base_atom = NULL;
		const int hydrogen_atomi = atomi;
		mol_atom *hydrogen_atom;
		int heavy_donor_base_count;

		if (!is_hydrogen(&ag->atoms[atomi], prm)) 
			continue;

		if (DEBUG_HBOND_DONOR)
			printf("atom index: %d\n", atomi);

		hydrogen_atom = &(ag->atoms[hydrogen_atomi]);
		assert(hydrogen_atom->nbondis == 1);

		if (DEBUG_HBOND_DONOR) {
			printf("\t----------------\n");
			printf("\t(hydrogen)\n");
			printf("\tatom type name prefix: %s\n",
			       prm->atoms[hydrogen_atom->atom_typen].typemaj);
			printf("\tatom type name suffix: %s\n",
			       prm->atoms[hydrogen_atom->atom_typen].typemin);
		}

		donor_atomi = bonded_atom_index(ag, hydrogen_atomi, 0);
		assert(donor_atomi < ag->natoms);
		hydrogen_atom->base = donor_atomi;
		donor_atom = &(ag->atoms[donor_atomi]);
		donor_atom->hprop |= HBOND_DONOR;

		if (DEBUG_HBOND_DONOR) {
			printf("\tdonor atom: \n");
			printf("\tatom index: %d\n", donor_atomi);
			printf("\t\tatom type name prefix: %s\n",
			       prm->atoms[donor_atom->atom_typen].typemaj);
			printf("\t\tatom type name suffix: %s\n",
			       prm->atoms[donor_atom->atom_typen].typemin);
			printf("\tdonor bonds: \n");
		}

		heavy_donor_base_count = 0;	// number of non-H donor bases
		for (j = 0; j < donor_atom->nbondis; j++) {
			int donor_base_atomi;
			mol_atom *donor_base_atom;

			donor_base_atomi = bonded_atom_index(ag, donor_atomi, j);
			assert(donor_base_atomi < ag->natoms);
			donor_base_atom = &(ag->atoms[donor_base_atomi]);

			if (is_hydrogen(donor_base_atom, prm))
				continue;

			if (DEBUG_HBOND_DONOR)
				printf("\t\t(heavy atom)\n");
			heavy_donor_base_count++;
			// take the first one as the donor base
			if (heavy_donor_base_count == 1) {
				heavy_donor_base_atomi = donor_base_atomi;
				heavy_donor_base_atom = donor_base_atom;
			}
			if (DEBUG_HBOND_DONOR)
				printf("\t\tbond %d:\n", j);
		}

		if (heavy_donor_base_count == 1)	// rotatable hydrogen
			hydrogen_atom->hprop |= ROTATABLE_HYDROGEN;
		else
			hydrogen_atom->hprop |= FIXED_HYDROGEN;

		if (DEBUG_HBOND_DONOR) {
			printf("\t\tthe chosen one:\n");
			printf("\t\tatom index: %d\n",
			       heavy_donor_base_atomi);
			printf("\t\t\tatom type name prefix: %s\n",
			       prm->
			       atoms[heavy_donor_base_atom->atom_typen].
			       typemaj);
			printf("\t\t\tatom type name suffix: %s\n",
			       prm->
			       atoms[heavy_donor_base_atom->atom_typen].
			       typemaj);
		}

		if (DEBUG_HBOND_DONOR) {
			printf("\n");
			printf("\t----------------\n");
		}
	}
}

//#define DEBUG_HBOND_ACCEPTOR 1
#define DEBUG_HBOND_ACCEPTOR 0
void mark_hbond_acceptors(struct atomgrp *ag, struct prm *prm)
{
	int i;
	mol_atom *atom;

	if (DEBUG_HBOND_ACCEPTOR)
		printf("acceptor\n");

	for (i = 0; i < ag->natoms; i++) {
		atom = &(ag->atoms[i]);
		if (DEBUG_HBOND_ACCEPTOR)
			printf("atom index: %d\n", i);

		if (prm->atoms[atom->atom_typen].typemin[0] == 'O'
		    || prm->atoms[atom->atom_typen].typemin[0] == 'N') {
			if (DEBUG_HBOND_ACCEPTOR) {
				printf("\t----------------\n");
				printf("\t(oxygen or nitrogen)\n");
				printf("\tatom type name prefix: %s\n",
				       prm->atoms[atom->atom_typen].typemaj);
				printf("\tatom type name suffix: %s\n",
				       prm->atoms[atom->atom_typen].typemin);
			}

			assert(atom->nbondis >= 1);

			if (atom->nbondis < 3)	// acceptor
			{
				int j;
				int heavy_acceptor_base_count = 0;

				if (DEBUG_HBOND_ACCEPTOR)
					printf("\t(fewer than three bonds)\n");

				atom->hprop |= HBOND_ACCEPTOR;
				atom->base = atom->base2 = -1;

				for (j = 0; j < atom->nbondis; j++) {
					int acceptor_base_atomi;
					mol_atom *acceptor_base_atom;

					acceptor_base_atomi =
					    bonded_atom_index(ag, i, j);
					assert(acceptor_base_atomi <
					       ag->natoms);
					acceptor_base_atom =
					    &(ag->atoms[acceptor_base_atomi]);

					if (!is_hydrogen
					    (acceptor_base_atom, prm)) {
						if (DEBUG_HBOND_ACCEPTOR)
							printf
							    ("\t\t(heavy atom)\n");
						heavy_acceptor_base_count++;
						if (heavy_acceptor_base_count == 1)	// take the first one as the acceptor base
							atom->base =
							    acceptor_base_atomi;
						else
							atom->base2 =
							    acceptor_base_atomi;
					}
					if (DEBUG_HBOND_ACCEPTOR)
						printf("\t\tbond %d:\n", j);
				}

				if (DEBUG_HBOND_ACCEPTOR) {
					printf("\n");
					printf("\t----------------\n");
				}
			}
		}

		if (DEBUG_HBOND_ACCEPTOR) {
			printf("\tatom type name prefix: %s\n",
			       prm->atoms[atom->atom_typen].typemaj);
			printf("\tatom type name suffix: %s\n",
			       prm->atoms[atom->atom_typen].typemin);
		}
	}
}

void fix_acceptor_bases(struct atomgrp *ag, struct prm *prm)
{
	int i;
	for (i = 0; i < ag->natoms; i++) {
		int j;
		mol_atom *atom = &(ag->atoms[i]);

		if (!(atom->hprop & HBOND_ACCEPTOR))
			continue;

		atom->base2 = -1;

		for (j = 0; j < atom->nbondis; j++) {
			int base_i = bonded_atom_index(ag, i, j);
			mol_atom *base = &(ag->atoms[base_i]);

			if (is_hydrogen(base, prm))
				continue;

			if (atom->base == -1)
				atom->base = base_i;
			else if ((atom->base2 == -1) && (atom->base != base_i))
				atom->base2 = base_i;
		}

		if (atom->base2 == -1) {
			for (j = 0; j < atom->nbondis; j++) {
				int base_i = bonded_atom_index(ag, i, j);

				if (atom->base != base_i) {
					atom->base2 = base_i;
					break;
				}
			}
		}
	}
}

static enum HB_Weight_Type get_hbond_weight_type(int hbe_type)
{
	switch (hbe_type) {
	case hbe_BBTURN:
	case hbe_BBHELIX:
		return hbw_SR_BB;
	case hbe_BBOTHER:
		return hbw_LR_BB;
	case hbe_SP2B:
	case hbe_SP3B:
	case hbe_RINGB:
	case hbe_BSC:
		return hbw_BB_SC;
	case hbe_SP2SC:
	case hbe_SP3SC:
	case hbe_RINGSC:
		return hbw_SC;
	default:
		return hbw_NONE;
	}
	return hbw_NONE;
}

static int classify_BB_by_separation(mol_atom * don, mol_atom * acc)
{
	int hbe;
	switch (don->comb_res_seq - acc->comb_res_seq) {
	case 0:
		hbe = hbe_NONE;
		break;
	case -2:
	case -1:
	case 1:
	case 2:
	case -3:
	case 3:
		hbe = hbe_BBTURN;
		break;
	case -4:
	case 4:
		hbe = hbe_BBHELIX;
		break;
	default:
		hbe = hbe_BBOTHER;
		break;
	}
	return hbe;
}

static int get_donor_chem_type(mol_atom * don)
{
	if (don->backbone)
		return hbdon_BB;
	else
		return hbdon_SC;
}

static int get_acceptor_chem_type(mol_atom * acc)
{
	if (acc->backbone)
		return hbacc_BB;

	switch (acc->hybridization) {
	case SP2_HYBRID:
		return hbacc_SP2;
	case SP3_HYBRID:
		return hbacc_SP3;
	case RING_HYBRID:
		return hbacc_RING;
	default:
		break;
	}

	return hbacc_NO;
}

static int get_hbe_type(mol_atom * atoms, mol_atom * hydro, mol_atom * acc)
{
	mol_atom *don;

	int don_type;
	int acc_type;
	if ((hydro->base < 0) && (!(acc->hprop & HBOND_ACCEPTOR)))
		return hbe_NONE;

	don = &(atoms[hydro->base]);

	don_type = get_donor_chem_type(don);
	acc_type = get_acceptor_chem_type(acc);

	if ((don_type == hbdon_BB) && (acc_type == hbacc_BB))
		return classify_BB_by_separation(don, acc);

	if (acc_type == hbacc_BB)
		return hbe_BSC;

	if (don_type == hbdon_BB) {
		switch (acc_type) {
		case hbacc_SP2:
			return hbe_SP2B;
		case hbacc_SP3:
			return hbe_SP3B;
		case hbacc_RING:
			return hbe_RINGB;
		default:
			return hbe_NONE;
		}
	} else {
		switch (acc_type) {
		case hbacc_SP2:
			return hbe_SP2SC;
		case hbacc_SP3:
			return hbe_SP3SC;
		case hbacc_RING:
			return hbe_RINGSC;
		default:
			return hbe_NONE;
		}
	}
}

double *alloc_categorized_hbondeng(void)
{
	return (double *)_mol_malloc((hbw_SC + 1) * sizeof(double));
}

void init_categorized_hbondeng(double *engcat)
{
	int i;
	if (engcat == NULL)
		return;

	for (i = 0; i <= hbw_SC; i++)
		engcat[i] = 0.0;
}

void print_categorized_hbondeng(double *engcat)
{
	if (engcat != NULL)
		printf
		    ("\t backbone-backbone ( short-range ) = %lf\n\t backbone-backbone ( long-range ) = %lf\n\t backbone-sidechain = %lf\n\t sidechain-sidechain = %lf\n",
		     engcat[hbw_SR_BB], engcat[hbw_LR_BB], engcat[hbw_BB_SC],
		     engcat[hbw_SC]);
}

void get_categorized_hbondeng(double *engcat, double *bb_bb_sr,
			      double *bb_bb_lr, double *bb_sc, double *sc_sc)
{
	*bb_bb_sr = engcat[hbw_SR_BB];
	*bb_bb_lr = engcat[hbw_LR_BB];
	*bb_sc = engcat[hbw_BB_SC];
	*sc_sc = engcat[hbw_SC];
}

void set_categorized_hbondeng(double *engcat, double bb_bb_sr, double bb_bb_lr,
			      double bb_sc, double sc_sc)
{
	engcat[hbw_SR_BB] = bb_bb_sr;
	engcat[hbw_LR_BB] = bb_bb_lr;
	engcat[hbw_BB_SC] = bb_sc;
	engcat[hbw_SC] = sc_sc;
}

/* create a unit vector pointing from the hydrogen towards the donor */
/* returns 0 on failure, 1 otherwise */

int create_donor_orientation_unit_vector(mol_atom * atoms, mol_atom * hydro,
					 struct dvector *HDunit,
					 FLOAT * invHDdis)
{
	mol_atom *donor;
	FLOAT HDdis;
	if (hydro->base < 0)
		return 0;

	donor = &(atoms[hydro->base]);

	HDunit->X = donor->X - hydro->X;
	HDunit->Y = donor->Y - hydro->Y;
	HDunit->Z = donor->Z - hydro->Z;

	HDdis =
	    HDunit->X * HDunit->X + HDunit->Y * HDunit->Y +
	    HDunit->Z * HDunit->Z;

	if (HDdis <= 0) {
		print_error
		    ("Failed to normalize the hydrogen to donor vector!");
		return 0;
	}

#ifdef USE_LONG_DOUBLE
	HDdis = sqrtl(HDdis);
#else
	HDdis = sqrt(HDdis);
#endif
	*invHDdis = 1 / HDdis;

	HDunit->X *= (*invHDdis);
	HDunit->Y *= (*invHDdis);
	HDunit->Z *= (*invHDdis);

	return 1;
}

// To construct the unit vector from acceptor base to the acceptor atom:

int create_base_to_acceptor_unit_vector(int hbe_type, mol_atom * atoms,
					mol_atom * acc, struct dvector *B,
					struct dvector *BAunit,
					FLOAT * invBAdis)
{
	mol_atom *base, *base2;
	FLOAT BAdis;

	if (acc->base < 0)
		return 0;
	else
		base = &(atoms[acc->base]);

	if ((hbe_type == hbe_RINGSC) || (hbe_type == hbe_RINGB)
	    || (hbe_type == hbe_SP3SC) || (hbe_type == hbe_SP3B)) {
		if (acc->base2 < 0)
			return 0;
		else
			base2 = &(atoms[acc->base2]);
	}

	switch (hbe_type) {
	case hbe_RINGSC:
	case hbe_RINGB: // RING: B is the midpoint between base and base2
		B->X = (base->X + base2->X) * 0.5;
		B->Y = (base->Y + base2->Y) * 0.5;
		B->Z = (base->Z + base2->Z) * 0.5;
		break;

	case hbe_SP3SC:
	case hbe_SP3B: // SP3: B is base2
		B->X = base2->X;
		B->Y = base2->Y;
		B->Z = base2->Z;
		break;

	default: // SP2: B is base
		B->X = base->X;
		B->Y = base->Y;
		B->Z = base->Z;
		break;
	}

	BAunit->X = acc->X - B->X;
	BAunit->Y = acc->Y - B->Y;
	BAunit->Z = acc->Z - B->Z;

	BAdis =
	    _mol_sq(BAunit->X) + _mol_sq(BAunit->Y) + _mol_sq(BAunit->Z);

	if (BAdis <= 0) {
		print_error("Failed to normalize the base to acceptor vector!");
		return 0;
	}
#ifdef USE_LONG_DOUBLE
	BAdis = sqrtl(BAdis);
#else
	BAdis = sqrt(BAdis);
#endif
	*invBAdis = 1 / BAdis;

	BAunit->X *= (*invBAdis);
	BAunit->Y *= (*invBAdis);
	BAunit->Z *= (*invBAdis);

	return 1;
}

static int hbe_is_BB_type(int hbe)
{
	return ((hbe >= hbe_BB) && (hbe <= hbe_BBOTHER));
}

int hbond_energy_computation(int hbe,
			     FLOAT AHdis,
			     FLOAT xD,
			     FLOAT xH,
			     FLOAT * energy,
			     FLOAT * dE_dr, FLOAT * dE_dxD, FLOAT * dE_dxH)
{
	FLOAT dAHdis = AHdis;
	FLOAT dxD = xD;
	FLOAT dxH = xH;

	FLOAT FSr = 0.0, FLr = 0.0, FxD = 0.0, FxH = 0.0; // fading intervals values
	FLOAT dFSr = 0.0, dFLr = 0.0, dFxD = 0.0, dFxH = 0.0; // fading intervals derivatives
	FLOAT Pr = 0.0, PSxD = 0.0, PSxH = 0.0, PLxD = 0.0, PLxH = 0.0; // polynomials values
	FLOAT dPr = 0.0, dPSxD = 0.0, dPSxH = 0.0, dPLxD = 0.0, dPLxH = 0.0; // polynomials derivatives

	if (energy != NULL)
		*energy = HB_ENG_MAX + 1.0f;
	if (dE_dr != NULL)
		*dE_dr = 0.0;
	if (dE_dxD != NULL)
		*dE_dxD = 0.0;
	if (dE_dxH != NULL)
		*dE_dxH = 0.0;

#ifdef USE_LONG_DOUBLE
	if ((fabsl(xD) > 1.0) || (fabsl(xH) > 1.0)) {
		print_error
		    ("Invalid angle value in hbond_energy_computation: xH = %Lf, xD = %Lf\n",
		     xH, xD);
		return 0;
	}
#else
	if ((fabs(xD) > 1.0) || (fabs(xH) > 1.0)) {
		print_error
		    ("Invalid angle value in hbond_energy_computation: xH = %lf, xD = %lf\n",
		     xH, xD);
		return 0;
	}
#endif

	if ((AHdis > MAX_AH) || (AHdis < MIN_AH) || (xH < MIN_xH)
	    || (xD < MIN_xD) || (xH > MAX_xH) || (xD > MAX_xD))
		return 0;

// fade functions

#if FADING_FUNCTION == LINEAR_FADE
	fade_xD_value_deriv(xD, &FxD, &dFxD);
	fade_xH_value_deriv(xH, &FxH, &dFxH);

	if (hbe_is_BB_type(hbe))
		fade_rBB_value_deriv(AHdis, &FSr, &dFSr);
	else {
		fade_rshort_value_deriv(AHdis, &FSr, &dFSr);
		fade_rlong_value_deriv(AHdis, &FLr, &dFLr);
	}
#elif FADING_FUNCTION == LOG_FADE
	fade_xD_log_value_deriv(xD + 1, &FxD, &dFxD);
	FxD -= 0.1;
	fade_xH_log_value_deriv(xH + 1, &FxH, &dFxH);
	FxH -= 0.1;

	if (hbe_is_BB_type(hbe))
		fade_rBB_log_value_deriv(AHdis, &FSr, &dFSr);
	else {
		fade_rshort_log_value_deriv(AHdis, &FSr, &dFSr);
		fade_rlong_log_value_deriv(AHdis, &FLr, &dFLr);
	}
#elif FADING_FUNCTION == BSPLINE_FADE
	fade_xD_bspline_value_deriv(xD /* + 1 */ , &FxD, &dFxD);	//FxD -= 0.1;
	fade_xH_bspline_value_deriv(xH /* + 1 */ , &FxH, &dFxH);	//FxH -= 0.1;

	if (hbe_is_BB_type(hbe))
		fade_rBB_bspline_value_deriv(AHdis, &FSr, &dFSr);
	else {
		fade_rshort_bspline_value_deriv(AHdis, &FSr, &dFSr);
		fade_rlong_bspline_value_deriv(AHdis, &FLr, &dFLr);
	}
#endif
	switch (hbe) {
	case hbe_NONE:
		break;

	case hbe_BBHELIX:	// hbstat = 1, polynomials 1,5,11
		AH_BBHelix(dAHdis, &Pr, &dPr);
		xD_BBHelix(dxD, &PSxD, &dPSxD);
		xH_BBHelix(dxH, &PSxH, &dPSxH);
		break;

	case hbe_BBOTHER:	// hbstat = 2, polynomials 2,6,11
		AH_BBOther(dAHdis, &Pr, &dPr);
		xD_BBOther(dxD, &PSxD, &dPSxD);
		xH_BBOther(dxH, &PSxH, &dPSxH);
		break;

	case hbe_BB:		// hbstat = 3, polynomials 3,7,13,8,14
	case hbe_BBTURN:
	case hbe_BSC:
	case hbe_SP2B:
	case hbe_SP2SC:
		AH_SP2(dAHdis, &Pr, &dPr);
		xD_SP2short(dxD, &PSxD, &dPSxD);
		xH_SP2short(dxH, &PSxH, &dPSxH);
		xD_SP2long(dxD, &PLxD, &dPLxD);
		xH_SP2long(dxH, &PLxH, &dPLxH);
		break;

	case hbe_SP3B:		// hbstat = 4, polynomials 4,9,15,10,15
	case hbe_SP3SC:
		AH_SP3(dAHdis, &Pr, &dPr);
		xD_SP3short(dxD, &PSxD, &dPSxD);
		xH_SP3(dxH, &PSxH, &dPSxH);
		xD_SP3long(dxD, &PLxD, &dPLxD);
		PLxH = PSxH;
		dPLxH = dPSxH;
		break;

	case hbe_RINGB:	// hbstat = 5, polynomials 3,7,16,8,16
	case hbe_RINGSC:
		AH_SP2(dAHdis, &Pr, &dPr);
		xD_SP2short(dxD, &PSxD, &dPSxD);
		xH_Ring(dxH, &PSxH, &dPSxH);
		xD_SP2long(dxD, &PLxD, &dPLxD);
		PLxH = PSxH;
		dPLxH = dPSxH;
		break;
	}

	if (energy != NULL)
		*energy =
		    Pr * FxD * FxH + FSr * (PSxD * FxH + FxD * PSxH) +
		    FLr * (PLxD * FxH + FxD * PLxH);
	if (dE_dr != NULL)
		*dE_dr =
		    dPr * FxD * FxH + dFSr * (PSxD * FxH + FxD * PSxH) +
		    dFLr * (PLxD * FxH + FxD * PLxH);
	if (dE_dxD != NULL)
		*dE_dxD =
		    dFxD * (Pr * FxH + FLr * PLxH + FSr * PSxH) +
		    FxH * (FSr * dPSxD + FLr * dPLxD);
	if (dE_dxH != NULL)
		*dE_dxH =
		    dFxH * (Pr * FxD + FLr * PLxD + FSr * PSxD) +
		    FxD * (FSr * dPSxH + FLr * dPLxH);

	return 1;
}

int water_mediated_hbond_compute_energy_and_gradient(int hbe, mol_atom * atoms_hydro, // atom group containing the hydrogen atom
						     int hydro_id, // index of the hydrogen atom in atoms_hydro
						     mol_atom * atoms_acc, // atom group containing the acceptor atom
						     int acc_id, // index of the acceptor atom in atoms_acc
						     mol_atom * atoms_wox, // atom group containing the water oxygen atom
						     int wox_id, // index of the water oxygen atom in atoms_wox
						     double *energy,
						     int comp_grad)
{
	struct dvector WHunit, HDunit, BAunit, AWunit, WDunit, B;
	FLOAT invHDdis, invBAdis, invAWdis, invWHdis, d_wp;

	mol_atom *hydro = &(atoms_hydro[hydro_id]);
	mol_atom *acc = &(atoms_acc[acc_id]);
	mol_atom *don = &(atoms_hydro[hydro->base]);
	mol_atom *wox = &(atoms_wox[wox_id]);
	FLOAT AWdis;
	FLOAT WHdis;
	FLOAT WDdis;
	FLOAT cosPsi;
	FLOAT cosOmg;
	FLOAT cosTht;
	FLOAT dE_dr;
	FLOAT u;

	if (energy != NULL)
		*energy = HB_ENG_MAX + 1.0f;

	if (!create_donor_orientation_unit_vector
	    (atoms_hydro, hydro, &HDunit, &invHDdis))
		return 0;

//        printf( "invHDdis = %lf\n", invHDdis );

//      if ( ( invHDdis < 0.8 ) || ( invHDdis > 1.25 ) ) return 1;

	if (!create_base_to_acceptor_unit_vector
	    (hbe, atoms_acc, acc, &B, &BAunit, &invBAdis))
		return 0;

	AWunit.X = wox->X - acc->X;
	AWunit.Y = wox->Y - acc->Y;
	AWunit.Z = wox->Z - acc->Z;

	WHunit.X = hydro->X - wox->X;
	WHunit.Y = hydro->Y - wox->Y;
	WHunit.Z = hydro->Z - wox->Z;

	WDunit.X = don->X - wox->X;
	WDunit.Y = don->Y - wox->Y;
	WDunit.Z = don->Z - wox->Z;

#ifdef USE_LONG_DOUBLE
	AWdis =
	    sqrtl(AWunit.X * AWunit.X + AWunit.Y * AWunit.Y +
		  AWunit.Z * AWunit.Z);
	WHdis =
	    sqrtl(WHunit.X * WHunit.X + WHunit.Y * WHunit.Y +
		  WHunit.Z * WHunit.Z);
	WDdis =
	    sqrtl(WDunit.X * WDunit.X + WDunit.Y * WDunit.Y +
		  WDunit.Z * WDunit.Z);
#else
	AWdis =
	    sqrt(AWunit.X * AWunit.X + AWunit.Y * AWunit.Y +
		 AWunit.Z * AWunit.Z);
	WHdis =
	    sqrt(WHunit.X * WHunit.X + WHunit.Y * WHunit.Y +
		 WHunit.Z * WHunit.Z);
	WDdis =
	    sqrt(WDunit.X * WDunit.X + WDunit.Y * WDunit.Y +
		 WDunit.Z * WDunit.Z);
#endif

	if ((AWdis - WDdis > D_WP_DIFFERENCE_TOLERANCE)
	    || (AWdis - WDdis < -D_WP_DIFFERENCE_TOLERANCE))
		return 0;

	if (AWdis <= 0) {
		print_error("Overlapping acceptor and water-oxygen atoms!");
		return 0;
	}

	if (WHdis <= 0) {
		print_error("Overlapping hydrogen and water-oxygen atoms!");
		return 0;
	}

	invAWdis = 1 / AWdis;

	AWunit.X *= invAWdis;
	AWunit.Y *= invAWdis;
	AWunit.Z *= invAWdis;

	invWHdis = 1 / WHdis;

	WHunit.X *= invWHdis;
	WHunit.Y *= invWHdis;
	WHunit.Z *= invWHdis;

	cosPsi = -dot_product(&BAunit, &AWunit);
	cosOmg = -dot_product(&AWunit, &WHunit);
	cosTht = -dot_product(&WHunit, &HDunit);

	if ((cosPsi > 0.0) || (cosPsi < -0.939693))
		return 0; // Psi (in degrees) is not in the allowed range [ 90, 160 ]
	if ((cosOmg > 0.173648) || (cosOmg < -0.766044))
		return 0; // Omega (in degrees) is not in the allowed range [ 80, 140 ]
	if ((cosTht > -0.5) || (cosTht < -1.0))
		return 0; // Theta (in degrees) is not in the allowed range [ 120, 180 ]

	d_wp = (AWdis + WDdis) * 0.5;

	if ((d_wp < 2.6) || (d_wp > 3.6))
		return 0;

	hbeng_HOH_bspline_value_deriv(d_wp, energy, &dE_dr);

	if (!comp_grad)
		return 1;

	// RAC: not sure about the gradients computed below...

	u = dE_dr * AWunit.X;
	wox->GX += -u;
	acc->GX += u;

	u = dE_dr * AWunit.Y;
	wox->GY += -u;
	acc->GY += u;

	u = dE_dr * AWunit.Z;
	wox->GZ += -u;
	acc->GZ += u;

	u = dE_dr * WDunit.X;
	wox->GX += u;
	don->GX += -u;

	u = dE_dr * WDunit.Y;
	wox->GY += u;
	don->GY += -u;

	u = dE_dr * WDunit.Z;
	wox->GZ += u;
	don->GZ += -u;

	return 1;
}

int hbond_energy_and_gradient_computation(int hbe, mol_atom * atoms_hydro,	// atom group containing the hydrogen atom
					  int hydro_id,	// index of the hydrogen atom in atoms_hydro
					  mol_atom * atoms_acc,	// atom group containing the acceptor atom
					  int acc_id,	// index of the acceptor atom in atoms_acc
					  double *energy, int comp_grad)
{
	struct dvector AHunit, HDunit, BAunit, B;
	FLOAT invAHdis, invHDdis, invBAdis;
	FLOAT AHdis;
	FLOAT xD;
	FLOAT xH;
	FLOAT dE_dr, dE_dxD, dE_dxH;
	FLOAT en;
	FLOAT u;
	FLOAT v;
	int en_comp;
	mol_atom *base;
	mol_atom *base2;

	mol_atom *hydro = &(atoms_hydro[hydro_id]);
	mol_atom *acc = &(atoms_acc[acc_id]);
	mol_atom *don = &(atoms_hydro[hydro->base]);

	if (!create_donor_orientation_unit_vector
	    (atoms_hydro, hydro, &HDunit, &invHDdis))
		return 0;

	if ((invHDdis < 0.8) || (invHDdis > 1.25)) {
		*energy = 0;
		return 1;
	}

	if (!create_base_to_acceptor_unit_vector
	    (hbe, atoms_acc, acc, &B, &BAunit, &invBAdis))
		return 0;

	AHunit.X = hydro->X - acc->X;
	AHunit.Y = hydro->Y - acc->Y;
	AHunit.Z = hydro->Z - acc->Z;

#ifdef USE_LONG_DOUBLE
	AHdis =
	    sqrtl(AHunit.X * AHunit.X + AHunit.Y * AHunit.Y +
		  AHunit.Z * AHunit.Z);
#else
	AHdis =
	    sqrt(AHunit.X * AHunit.X + AHunit.Y * AHunit.Y +
		 AHunit.Z * AHunit.Z);
#endif

	if (AHdis <= 0) {
		print_error("Overlapping acceptor and hydrogen!");
		return 0;
	}

	invAHdis = 1 / AHdis;

	AHunit.X *= invAHdis;
	AHunit.Y *= invAHdis;
	AHunit.Z *= invAHdis;

	xD = dot_product(&AHunit, &HDunit);
	xH = dot_product(&BAunit, &AHunit);

	en_comp =
	    hbond_energy_computation(hbe, AHdis, xD, xH, &en, &dE_dr, &dE_dxD,
				     &dE_dxH);

	if (!en_comp)
		return 0;

	*energy = en;

	if (*energy >= HB_ENG_MAX) {
		*energy = 0;
		return 1;
	}

	if (!comp_grad)
		return 1;

	base = &(atoms_acc[acc->base]);
	base2 = NULL;

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB) || (hbe == hbe_SP3SC)
	    || (hbe == hbe_SP3B)) {
		if (acc->base2 < 0) {
			print_error("Acceptor base 2 missing!");
			return 0;
		}
		base2 = &(atoms_acc[acc->base2]);
	}

	u = dE_dr * AHunit.X;
	hydro->GX += -u;
	acc->GX += u;

	u = dE_dr * AHunit.Y;
	hydro->GY += -u;
	acc->GY += u;

	u = dE_dr * AHunit.Z;
	hydro->GZ += -u;
	acc->GZ += u;

	u = -dE_dxD * invAHdis * (xD * AHunit.X - HDunit.X);
	v = -dE_dxD * invHDdis * (AHunit.X - xD * HDunit.X);
	hydro->GX += -(u + v);
	acc->GX += u;
	don->GX += v;

	u = -dE_dxD * invAHdis * (xD * AHunit.Y - HDunit.Y);
	v = -dE_dxD * invHDdis * (AHunit.Y - xD * HDunit.Y);
	hydro->GY += -(u + v);
	acc->GY += u;
	don->GY += v;

	u = -dE_dxD * invAHdis * (xD * AHunit.Z - HDunit.Z);
	v = -dE_dxD * invHDdis * (AHunit.Z - xD * HDunit.Z);
	hydro->GZ += -(u + v);
	acc->GZ += u;
	don->GZ += v;

	u = -dE_dxH * invAHdis * (BAunit.X - xH * AHunit.X);
	v = -dE_dxH * invBAdis * (xH * BAunit.X - AHunit.X);
	hydro->GX += u;
	acc->GX += -(u + v);

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB)) {
		base->GX += 0.5 * v;
		base2->GX += 0.5 * v;
	} else if ((hbe == hbe_SP3SC) || (hbe == hbe_SP3B))
		base2->GX += v;
	else
		base->GX += v;

	u = -dE_dxH * invAHdis * (BAunit.Y - xH * AHunit.Y);
	v = -dE_dxH * invBAdis * (xH * BAunit.Y - AHunit.Y);
	hydro->GY += u;
	acc->GY += -(u + v);

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB)) {
		base->GY += 0.5 * v;
		base2->GY += 0.5 * v;
	} else if ((hbe == hbe_SP3SC) || (hbe == hbe_SP3B))
		base2->GY += v;
	else
		base->GY += v;

	u = -dE_dxH * invAHdis * (BAunit.Z - xH * AHunit.Z);
	v = -dE_dxH * invBAdis * (xH * BAunit.Z - AHunit.Z);
	hydro->GZ += u;
	acc->GZ += -(u + v);

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB)) {
		base->GZ += 0.5 * v;
		base2->GZ += 0.5 * v;
	} else if ((hbe == hbe_SP3SC) || (hbe == hbe_SP3B))
		base2->GZ += v;
	else
		base->GZ += v;

	return 1;
}

static double get_pairwise_hbondeng(mol_atom * atoms_hydro, int hydro_id,
				    mol_atom * atoms_acc, int acc_id,
				    struct agsetup *ags, double *engcat,
				    double rc2, int comp_grad)
{
	mol_atom *hydro = &(atoms_hydro[hydro_id]);
	mol_atom *acc = &(atoms_acc[acc_id]);

	double dx = hydro->X - acc->X;
	double dy = hydro->Y - acc->Y;
	double dz = hydro->Z - acc->Z;

	double d2 = dx * dx + dy * dy + dz * dz;
	int hbe;
	int ka1, ka2, ka3;
	double en;

	if (d2 > rc2)
		return 0;

	hbe = get_hbe_type(atoms_hydro, hydro, acc);

	if (hbe == hbe_NONE)
		return 0;

	if (hydro_id < acc_id) {
		ka1 = hydro_id;
		ka2 = acc_id;
	} else {
		ka1 = acc_id;
		ka2 = hydro_id;
	}

	ka3 = exta(ka1, ka2, ags->excl_list, ags->pd1, ags->pd2, ags->ndm);

	if (ka3 > 0)
		return 0;

	en = 0;

	if (hbond_energy_and_gradient_computation
	    (hbe, atoms_hydro, hydro_id, atoms_acc, acc_id, &en, comp_grad)) {
		if (engcat != NULL) {
			enum HB_Weight_Type hbw = get_hbond_weight_type(hbe);
			engcat[hbw] += en;
		}

		return en;
	}

	return 0;
}

double get_pairwise_hbondeng_nblist(mol_atom * atoms_hydro, int hydro_id,
				    mol_atom * atoms_acc, int acc_id,
				    double *engcat, double rc2, int comp_grad)
{
	mol_atom *hydro = &(atoms_hydro[hydro_id]);
	mol_atom *acc = &(atoms_acc[acc_id]);

	double dx = hydro->X - acc->X;
	double dy = hydro->Y - acc->Y;
	double dz = hydro->Z - acc->Z;

	double d2 = dx * dx + dy * dy + dz * dz;

	int hbe;
	double en;

	if (d2 > rc2)
		return 0;

	hbe = get_hbe_type(atoms_hydro, hydro, acc);

	if (hbe == hbe_NONE)
		return 0;

	en = 0;

	if (hbond_energy_and_gradient_computation
	    (hbe, atoms_hydro, hydro_id, atoms_acc, acc_id, &en, comp_grad)) {
		if (engcat != NULL) {
			enum HB_Weight_Type hbw = get_hbond_weight_type(hbe);
			engcat[hbw] += en;
		}

		return en;
	}

	return 0;
}

static double get_water_mediated_pairwise_hbondeng(mol_atom * atoms_hydro,
						   int hydro_id,
						   mol_atom * atoms_acc,
						   int acc_id,
						   mol_atom * atoms_wox,
						   int wox_id, int comp_grad)
{
	mol_atom *hydro = &(atoms_hydro[hydro_id]);
	mol_atom *acc = &(atoms_acc[acc_id]);

	int hbe = get_hbe_type(atoms_hydro, hydro, acc);

	double en = 0;

	water_mediated_hbond_compute_energy_and_gradient(hbe, atoms_hydro,
							 hydro_id, atoms_acc,
							 acc_id, atoms_wox,
							 wox_id, &en,
							 comp_grad);

	return en;
}

void hbondeng_octree_single_mol(OCTREE_PARAMS * octpar, double *energy)
{
	OCTREE *octree_static = octpar->octree_static;
	OCTREE *octree_moving = octpar->octree_moving;
	double dist_cutoff = octpar->dist_cutoff;
	double *trans_mat = octpar->trans;

	OCTREE_PARAMS *prms = (OCTREE_PARAMS *) octpar->proc_func_params;
	struct agsetup *ags = prms->ags;

	OCTREE_NODE *snode = &(octree_static->nodes[octpar->node_static]);
	OCTREE_NODE *mnode = &(octree_moving->nodes[octpar->node_moving]);

	double rc = dist_cutoff;
	double rc2 = rc * rc;

	int nf = snode->nfixed;

	double *engcat = prms->engcat;

	*energy = 0;

	if ((trans_mat != NULL) || (mnode->n - mnode->nfixed <= snode->n)) {
		int i;
		for (i = mnode->nfixed; i < mnode->n; i++) {
			int ai = mnode->indices[i];
			mol_atom *atom_i = &(octree_moving->atoms[ai]);
			double x, y, z;
			double d2;

			if (!(atom_i->hprop & DONATABLE_HYDROGEN)
			    && !(atom_i->hprop & HBOND_ACCEPTOR))
				continue;

			x = atom_i->X, y = atom_i->Y, z = atom_i->Z;

			if (trans_mat != NULL)
				transform_point(x, y, z, trans_mat, &x, &y, &z);	// defined in octree.h

			d2 =
			    min_pt2bx_dist2(snode->lx, snode->ly, snode->lz,
					    snode->dim, x, y, z);

			if (rc2 < d2)
				continue;

			if (atom_i->hprop & DONATABLE_HYDROGEN) {
				int j;
				for (j = 0; j < snode->n; j++) {
					int aj = snode->indices[j];
					mol_atom *atom_j;

					if ((j >= nf) && (aj <= ai))
						continue;

					atom_j =
					    &(octree_static->atoms[aj]);

					if (!(atom_j->hprop & HBOND_ACCEPTOR))
						continue;

					(*energy) += get_pairwise_hbondeng
					    (octree_moving->atoms, ai,
					     octree_static->atoms, aj,
					     ags, engcat, rc2, 1);
				}
			} else {
				int j;
				for (j = 0; j < snode->n; j++) {
					int aj = snode->indices[j];
					mol_atom *atom_j;

					if ((j >= nf) && (aj <= ai))
						continue;

					atom_j =
					    &(octree_static->atoms[aj]);

					if (!(atom_j->hprop &
					     DONATABLE_HYDROGEN))
						continue;

					(*energy) += get_pairwise_hbondeng
					    (octree_static->atoms, aj,
					     octree_moving->atoms, ai,
					     ags, engcat, rc2, 1);
				}
			}
		}
	} else {
		int i;
		for (i = 0; i < snode->n; i++) {
			int ai = snode->indices[i];
			mol_atom *atom_i = &(octree_static->atoms[ai]);
			double x, y, z;
			double d2;

			if (!(atom_i->hprop & DONATABLE_HYDROGEN)
			    && !(atom_i->hprop & HBOND_ACCEPTOR))
				continue;

			x = atom_i->X, y = atom_i->Y, z = atom_i->Z;

			d2 =
			    min_pt2bx_dist2(mnode->lx, mnode->ly, mnode->lz,
					    mnode->dim, x, y, z);

			if (rc2 < d2)
				continue;

			if (atom_i->hprop & DONATABLE_HYDROGEN) {
				int j;
				for (j = mnode->nfixed; j < mnode->n; j++) {
					int aj = mnode->indices[j];
					mol_atom *atom_j;

					if ((i >= nf) && (aj >= ai))
						continue;

					atom_j =
					    &(octree_moving->atoms[aj]);

					if (!(atom_j->hprop & HBOND_ACCEPTOR))
						continue;

					(*energy) += get_pairwise_hbondeng
					    (octree_static->atoms, ai,
					     octree_moving->atoms, aj,
					     ags, engcat, rc2, 1);
				}
			} else {
				int j;
				for (j = mnode->nfixed; j < mnode->n; j++) {
					int aj = mnode->indices[j];
					mol_atom *atom_j;

					if ((i >= nf) && (aj >= ai))
						continue;

					atom_j =
					    &(octree_moving->atoms[aj]);

					if (!(atom_j->hprop &
					     DONATABLE_HYDROGEN))
						continue;

					(*energy) += get_pairwise_hbondeng
					    (octree_moving->atoms, aj,
					     octree_static->atoms, ai,
					     ags, engcat, rc2, 1);
				}
			}
		}
	}
}

void hbondeng(struct atomgrp *ag, double *energy, struct nblist *nblst)
{
	double rc = nblst->nbcof;
	double rc2 = rc * rc;
	int i;

	(*energy) = 0;

	for (i = 0; i < nblst->nfat; i++) {
		int ai = nblst->ifat[i];
		mol_atom *atom_i = &(ag->atoms[ai]);
		int n2;
		int *p;
		int j;

		if (!(atom_i->hprop & DONATABLE_HYDROGEN)
		    && !(atom_i->hprop & HBOND_ACCEPTOR))
			continue;

		n2 = nblst->nsat[i];
		p = nblst->isat[i];

		for (j = 0; j < n2; j++) {
			int aj = p[j];
			mol_atom *atom_j = &(ag->atoms[aj]);

			if (atom_i->hprop & DONATABLE_HYDROGEN) {
				if (!(atom_j->hprop & HBOND_ACCEPTOR))
					continue;

				(*energy) +=
				    get_pairwise_hbondeng_nblist(ag->atoms, ai,
								 ag->atoms, aj,
								 NULL, rc2, 1);
			} else {
				if (!(atom_j->hprop & DONATABLE_HYDROGEN))
					continue;

				(*energy) +=
				    get_pairwise_hbondeng_nblist(ag->atoms, aj,
								 ag->atoms, ai,
								 NULL, rc2, 1);
			}
		}
	}
}

void hbondengcat(struct atomgrp *ag, double *energy, struct nblist *nblst)
{
	double rc = nblst->nbcof;
	double rc2 = rc * rc;
	int i;
	for (i = 0; i < nblst->nfat; i++) {
		int ai = nblst->ifat[i];
		mol_atom *atom_i = &(ag->atoms[ai]);
		int n2;
		int *p;
		int j;

		if (!(atom_i->hprop & DONATABLE_HYDROGEN)
		    && !(atom_i->hprop & HBOND_ACCEPTOR))
			continue;

		n2 = nblst->nsat[i];
		p = nblst->isat[i];

		for (j = 0; j < n2; j++) {
			int aj = p[j];
			mol_atom *atom_j = &(ag->atoms[aj]);

			if (atom_i->hprop & DONATABLE_HYDROGEN) {
				if (!(atom_j->hprop & HBOND_ACCEPTOR))
					continue;

				get_pairwise_hbondeng_nblist(ag->atoms, ai,
							     ag->atoms, aj,
							     energy, rc2, 1);
			} else {
				if (!(atom_j->hprop & DONATABLE_HYDROGEN))
					continue;

				get_pairwise_hbondeng_nblist(ag->atoms, aj,
							     ag->atoms, ai,
							     energy, rc2, 1);
			}
		}
	}
}

static int hb_bonded(mol_atom * ai, mol_atom * aj)
{
	int i;
	if (ai->nbonds > aj->nbonds) {
		mol_atom *t = ai;
		ai = aj;
		aj = t;
	}

	for (i = 0; i < ai->nbonds; i++) {
		mol_bond *b = ai->bonds[i];

		if ((b->a0 == aj) || (b->a1 == aj))
			return 1;
	}

	return 0;
}

void hbondeng_all(struct atomgrp *ag, double *energy, struct nblist *nblst)
{
	double rc = nblst->nbcof;
	double rc2 = rc * rc;
	int i;
	double en[hbw_SC + 1];


	(*energy) = 0;

	for (i = 0; i < hbw_SC; i++)
		en[i] = 0;

	printf("ag->natoms = %d\n", ag->natoms);

	for (i = 0; i < ag->natoms; i++) {
		int ai = i;
		mol_atom *atom_i = &(ag->atoms[ai]);
		int j;

		if (!(atom_i->hprop & DONATABLE_HYDROGEN))
			continue;

		for (j = 0; j < ag->natoms; j++) {
			int aj = j;
			mol_atom *atom_j = &(ag->atoms[aj]);
			int hbe;
			enum HB_Weight_Type hbw;

			if ((atom_i->comb_res_seq ==
			     atom_j->
			     comb_res_seq) /*( i == j ) */ ||hb_bonded(atom_i,
								       atom_j))
				continue;

			if (!(atom_j->hprop & HBOND_ACCEPTOR))
				continue;

			en[0] =
			    get_pairwise_hbondeng_nblist(ag->atoms, ai,
							 ag->atoms, aj, NULL,
							 rc2, 1);

			hbe = get_hbe_type(ag->atoms, atom_i, atom_j);
			hbw = get_hbond_weight_type(hbe);

			if (en[0] >= 0)
				continue;

			(*energy) += en[0];

			if (hbw == hbw_NONE)
				continue;

			en[hbw] += en[0];
			printf("%d %d %s  %d %d %s  %d %d   %lf\n", ai + 1,
			       atom_i->res_seq, atom_i->name, aj + 1,
			       atom_j->res_seq, atom_j->name, hbe, hbw, en[0]);
		}
	}

	printf("total = %lf, SR_BB = %lf, LR_BB = %lf, BB_SC = %lf, SC = %lf\n",
	       (*energy), en[hbw_SR_BB], en[hbw_LR_BB], en[hbw_BB_SC],
	       en[hbw_SC]);
}

void hbondeng_bbexc(struct atomgrp *ag, double *energy, struct nblist *nblst)
{

	double rc = nblst->nbcof;
	double rc2 = rc * rc;
	int i;
	int ai;
	double en[hbw_SC + 1];
	int *blacklist;

	(*energy) = 0;

	for (i = 0; i < hbw_SC; i++)
		en[i] = 0;

	printf("ag->natoms = %d\n", ag->natoms);

	blacklist = _mol_malloc(sizeof(int) * ag->natoms);
	for (i = 0; i < ag->natoms; i++)
		blacklist[i] = 0;

	for (ai = 0; ai < ag->natoms; ai++) {
		int aj;
		mol_atom *atom_i = &(ag->atoms[ai]);

		if (!(atom_i->hprop & DONATABLE_HYDROGEN))
			continue;

		for (aj = 0; aj < ag->natoms; aj++) {
			int hbe;
			mol_atom *atom_j = &(ag->atoms[aj]);

			if ((atom_i->comb_res_seq == atom_j->comb_res_seq) 
			    ||hb_bonded(atom_i, atom_j))
				continue;

			if (!(atom_j->hprop & HBOND_ACCEPTOR))
				continue;

			en[0] =
			    get_pairwise_hbondeng_nblist(ag->atoms, ai,
							 ag->atoms, aj, NULL,
							 rc2, 1);

			hbe = get_hbe_type(ag->atoms, atom_i, atom_j);

			if ((en[0] < 0)
			    && ((hbe == 1) || (hbe == 2) || (hbe == 3)
				|| (hbe == 4))) {
				blacklist[ai] = 1;
				blacklist[aj] = 1;
			}
		}
	}

	for (ai = 0; ai < ag->natoms; ai++) {
		int aj;
		mol_atom *atom_i;
		if (blacklist[ai] == 1)
			continue;

		atom_i = &(ag->atoms[ai]);

		if (!(atom_i->hprop & DONATABLE_HYDROGEN))
			continue;

		for (aj = 0; aj < ag->natoms; aj++) {
			mol_atom *atom_j;
			int hbe;
			enum HB_Weight_Type hbw;

			if (blacklist[aj] == 1)
				continue;

			atom_j = &(ag->atoms[aj]);

			if ((atom_i->comb_res_seq == atom_j-> comb_res_seq) 
			    ||hb_bonded(atom_i, atom_j))
				continue;

			if (!(atom_j->hprop & HBOND_ACCEPTOR))
				continue;

			en[0] =
			    get_pairwise_hbondeng_nblist(ag->atoms, ai,
							 ag->atoms, aj, NULL,
							 rc2, 1);

			hbe = get_hbe_type(ag->atoms, atom_i, atom_j);
			hbw = get_hbond_weight_type(hbe);

			if (en[0] >= 0)
				continue;

			(*energy) += en[0];

			if (hbw == hbw_NONE)
				continue;

			en[hbw] += en[0];
		}
	}

	printf("total = %lf, SR_BB = %lf, LR_BB = %lf, BB_SC = %lf, SC = %lf\n",
	       (*energy), en[hbw_SR_BB], en[hbw_LR_BB], en[hbw_BB_SC],
	       en[hbw_SC]);

	// To display the excluded atoms
/*
    for( int i = 0; i < ag->natoms; i++)
       if(blacklist[i]==1){
	   int ai = i;
	   mol_atom *atom_i = &( ag->atoms[ ai ] );
	   printf("%d %d %s\n", ai + 1, atom_i->res_seq, atom_i->name);
	}
*/
}

void water_mediated_hbondeng(struct atomgrp *ag, double *energy)
{
	struct prm *prm = (struct prm *)ag->prm;
	int ak;

	(*energy) = 0;

	for (ak = 0; ak < ag->natoms; ak++)	// Water Oxygen
	{
		mol_atom *atom_k = &(ag->atoms[ak]);
		int nb;
		int ai;

		if ((ag->res_type[atom_k->res_num] != HOH)
		    || (prm->atoms[atom_k->atom_typen].typemin[0] != 'O'))
			continue;

		nb = 0;

		for (ai = 0; ai < ag->natoms; ai++) {
			mol_atom *atom_i = &(ag->atoms[ai]);

			double dx = atom_k->X - atom_i->X;
			double dy = atom_k->Y - atom_i->Y;
			double dz = atom_k->Z - atom_i->Z;

			double d2 = dx * dx + dy * dy + dz * dz;

			if (d2 > 64)
				continue;

			nb++;

			if (nb == 10)
				break;
		}

		if (nb < 10)
			continue;

		for (ai = 0; ai < ag->natoms; ai++)	// Donated Hydrogen which is NOT H of a Water molecule
		{
			mol_atom *atom_i = &(ag->atoms[ai]);
			int aj;

			if (!(atom_i->hprop & DONATABLE_HYDROGEN)
			    || (ag->res_type[atom_i->res_num] == HOH))
				continue;

			for (aj = 0; aj < ag->natoms; aj++)	// Acceptor that is not the same as Atom_i and has no hbond with Atom_i
			{
				mol_atom *atom_j;
				double en;
				int al;
				mol_atom *atom_l;
				if (ak == aj)
					continue;

				atom_j = &(ag->atoms[aj]);

				if ((atom_i->comb_res_seq ==
				     atom_j->comb_res_seq)
				    || hb_bonded(atom_i, atom_j))
					continue;

				if (!(atom_j->hprop & HBOND_ACCEPTOR))
					continue;

				en =
				    get_water_mediated_pairwise_hbondeng(ag->
									 atoms,
									 ai,
									 ag->
									 atoms,
									 aj,
									 ag->
									 atoms,
									 ak, 1);

				if (en >= 0)
					continue;

				(*energy) += en;

				al = atom_i->base;
				atom_l = &(ag->atoms[al]);	//   Base of the Hydrogen atom_i

				printf
				    ("Water Oxygen: %s ( %d %d ), Acceptor: %s ( %d %d ), Donor: %s ( %d %d ), Hydrogen: %s ( %d %d ), Energy: %lf\n",
				     atom_k->name, ak + 1, atom_k->res_seq,
				     atom_j->name, aj + 1, atom_j->res_seq,
				     atom_l->name, al + 1, atom_l->res_seq,
				     atom_i->name, ai + 1, atom_i->res_seq, en);
			}
		}
	}

//  printf( "total water-mediated hbond energy = %lf\n", ( *energy ) );
}

int residual_acceptor_valency(mol_atom * atom, struct prm *prm)
{
	switch (prm->atoms[atom->atom_typen].typemin[0]) {
	case 'O':
		return 2; //( 3 - atom->nbondis );
	case 'N':
		return 1; //( 4 - atom->nbondis );
	default:
		return 0;
	}
}

int init_flow_struct(FLOW_STRUCT * flow_struct)
{
	flow_struct->max_n_fedge = 0;
	flow_struct->max_n = 0;

	flow_struct->flow_edge = NULL;
	flow_struct->flow_map = NULL;
	flow_struct->cap = NULL;
	flow_struct->cost = NULL;
	flow_struct->deg = NULL;
	flow_struct->par = NULL;
	flow_struct->q = NULL;
	flow_struct->inq = NULL;
	flow_struct->pi = NULL;
	flow_struct->d = NULL;
	flow_struct->adj = NULL;
	flow_struct->fnet = NULL;

	flow_struct->max_n_fedge = 100;
	flow_struct->max_n = 50;

	flow_struct->flow_edge =
	    (FLOW_EDGE *) _mol_calloc(flow_struct->max_n_fedge,
				      sizeof(FLOW_EDGE));
	flow_struct->flow_map =
	    (FLOW_MAP *) _mol_calloc(flow_struct->max_n, sizeof(FLOW_MAP));
	flow_struct->cap =
	    (int *)_mol_calloc(flow_struct->max_n * flow_struct->max_n,
			       sizeof(int));
	flow_struct->cost =
	    (FLOAT *) _mol_calloc(flow_struct->max_n * flow_struct->max_n,
				  sizeof(FLOAT));
	flow_struct->deg = (int *)_mol_calloc(flow_struct->max_n, sizeof(int));
	flow_struct->par = (int *)_mol_calloc(flow_struct->max_n, sizeof(int));
	flow_struct->q = (int *)_mol_calloc(flow_struct->max_n, sizeof(int));
	flow_struct->inq = (int *)_mol_calloc(flow_struct->max_n, sizeof(int));
	flow_struct->pi =
	    (FLOAT *) _mol_calloc(flow_struct->max_n, sizeof(FLOAT));
	flow_struct->d =
	    (FLOAT *) _mol_calloc(flow_struct->max_n, sizeof(FLOAT));
	flow_struct->adj =
	    (int *)_mol_calloc(flow_struct->max_n * flow_struct->max_n,
			       sizeof(int));
	flow_struct->fnet =
	    (int *)_mol_calloc(flow_struct->max_n * flow_struct->max_n,
			       sizeof(int));

	if ((flow_struct->flow_edge == NULL) || (flow_struct->flow_map == NULL)
	    || (flow_struct->cap == NULL) || (flow_struct->cost == NULL)
	    || (flow_struct->deg == NULL)
	    || (flow_struct->par == NULL) || (flow_struct->q == NULL)
	    || (flow_struct->inq == NULL)
	    || (flow_struct->pi == NULL) || (flow_struct->d == NULL)
	    || (flow_struct->adj == NULL) || (flow_struct->fnet == NULL)) {
		print_error
		    ("Failed to allocate memory for hbond flow network!");
		free_flow_struct(flow_struct);
		return 0;
	}

	return 1;
}

void free_flow_struct(FLOW_STRUCT * flow_struct)
{
	flow_struct->max_n_fedge = 0;
	flow_struct->max_n = 0;

	freeMem(flow_struct->flow_edge);
	freeMem(flow_struct->flow_map);
	freeMem(flow_struct->cap);
	freeMem(flow_struct->cost);
	freeMem(flow_struct->deg);
	freeMem(flow_struct->par);
	freeMem(flow_struct->q);
	freeMem(flow_struct->inq);
	freeMem(flow_struct->pi);
	freeMem(flow_struct->d);
	freeMem(flow_struct->adj);
	freeMem(flow_struct->fnet);
}

int reinit_flow_struct(FLOW_STRUCT * flow_struct, int max_n_fedge, int max_n)
{
	if (max_n_fedge > flow_struct->max_n_fedge) {
		flow_struct->max_n_fedge = max_n_fedge;
		flow_struct->flow_edge =
		    (FLOW_EDGE *) _mol_realloc(flow_struct->flow_edge,
					       flow_struct->max_n_fedge *
					       sizeof(FLOW_EDGE));
	}

	if (max_n > flow_struct->max_n) {
		flow_struct->max_n = max_n;

		flow_struct->flow_map =
		    (FLOW_MAP *) _mol_realloc(flow_struct->flow_map,
					      flow_struct->max_n *
					      sizeof(FLOW_MAP));
		flow_struct->cap =
		    (int *)_mol_realloc(flow_struct->cap,
					flow_struct->max_n *
					flow_struct->max_n * sizeof(int));
		flow_struct->cost =
		    (FLOAT *) _mol_realloc(flow_struct->cost,
					   flow_struct->max_n *
					   flow_struct->max_n * sizeof(FLOAT));
		flow_struct->deg =
		    (int *)_mol_realloc(flow_struct->deg,
					flow_struct->max_n * sizeof(int));
		flow_struct->par =
		    (int *)_mol_realloc(flow_struct->par,
					flow_struct->max_n * sizeof(int));
		flow_struct->q =
		    (int *)_mol_realloc(flow_struct->q,
					flow_struct->max_n * sizeof(int));
		flow_struct->inq =
		    (int *)_mol_realloc(flow_struct->inq,
					flow_struct->max_n * sizeof(int));
		flow_struct->pi =
		    (FLOAT *) _mol_realloc(flow_struct->pi,
					   flow_struct->max_n * sizeof(FLOAT));
		flow_struct->d =
		    (FLOAT *) _mol_realloc(flow_struct->d,
					   flow_struct->max_n * sizeof(FLOAT));
		flow_struct->adj =
		    (int *)_mol_realloc(flow_struct->adj,
					flow_struct->max_n *
					flow_struct->max_n * sizeof(int));
		flow_struct->fnet =
		    (int *)_mol_realloc(flow_struct->fnet,
					flow_struct->max_n *
					flow_struct->max_n * sizeof(int));
	}

	if ((flow_struct->flow_edge == NULL) || (flow_struct->flow_map == NULL)
	    || (flow_struct->cap == NULL) || (flow_struct->cost == NULL)
	    || (flow_struct->deg == NULL)
	    || (flow_struct->par == NULL) || (flow_struct->q == NULL)
	    || (flow_struct->inq == NULL)
	    || (flow_struct->pi == NULL) || (flow_struct->d == NULL)
	    || (flow_struct->adj == NULL) || (flow_struct->fnet == NULL)) {
		print_error
		    ("Failed to reallocate memory for hbond flow network!");
		free_flow_struct(flow_struct);
		return 0;
	}

	return 1;
}

int dijkstra(int n, int *deg, int *adj, int *cap, int *fnet, FLOAT * cost,
	     int *q, int *inq, int *par, FLOAT * pi, FLOAT * d, FLOAT inf)
{
	int i;
	int qs;
	for (i = 0; i < n; i++) {
		d[i] = inf;
		par[i] = inq[i] = -1;
	}

	d[0] = 0;
	qs = 1;

	q[0] = inq[0] = 0;
	par[0] = n;

	while (qs) {
		int u = q[0];
		int j, k, v;
		inq[u] = -1;

		q[0] = q[--qs];
		if (qs)
			inq[q[0]] = 0;

		for (i = 0, j = 2 * i + 1; j < qs; i = j, j = 2 * i + 1) {
			int t;
			if ((j + 1 < qs) && (d[q[j + 1]] < d[q[j]]))
				j++;
			if (d[q[j]] >= d[q[i]])
				break;

			t = q[i];
			q[i] = q[j];
			q[j] = t;
			t = inq[q[i]];
			inq[q[i]] = inq[q[j]];
			inq[q[j]] = t;
		}

		for (k = 0, v = ARY(adj, n, u, k); k < deg[u];
		     v = ARY(adj, n, u, ++k)) {
			if (ARY(fnet, n, v, u)
			    && (d[v] >
				d[u] + pi[u] - pi[v] - ARY(cost, n, v, u))) {
				d[v] =
				    d[u] + pi[u] - pi[v] - ARY(cost, n, v, u);
				par[v] = u;
			}

			if ((ARY(fnet, n, u, v) < ARY(cap, n, u, v))
			    && (d[v] >
				d[u] + pi[u] - pi[v] + ARY(cost, n, u, v))) {
				d[v] =
				    d[u] + pi[u] - pi[v] + ARY(cost, n, u, v);
				par[v] = u;
			}

			if (par[v] == u) {
				if (inq[v] < 0) {
					q[qs] = v;
					inq[v] = qs++;
				}

				for (i = inq[v], j = (i >> 1);
				     d[q[i]] < d[q[j]]; i = j, j = (i >> 1)) {
					int t;

					t = q[i];
					q[i] = q[j];
					q[j] = t;
					t = inq[q[i]];
					inq[q[i]] = inq[q[j]];
					inq[q[j]] = t;
				}
			}
		}
	}

	for (i = 0; i < n; i++)
		if (pi[i] < inf)
			pi[i] += d[i];

	return (par[n - 1] >= 0);
}

int min_cost_max_flow(int n, FLOW_STRUCT * fs, FLOAT inf, FLOAT * fcost)
{
	int *cap = fs->cap;
	FLOAT *cost = fs->cost;
	int *deg = fs->deg;
	int *par = fs->par;
	int *q = fs->q;
	int *inq = fs->inq;
	FLOAT *pi = fs->pi;
	FLOAT *d = fs->d;
	int *adj = fs->adj;
	int *fnet = fs->fnet;
	int i;
	int flow;

	for (i = 0; i < n; i++) {
		deg[i] = 0;
		pi[i] = 0;
	}

	for (i = 0; i < n * n; i++)
		fnet[i] = 0;

	for (i = 0; i < n; i++) {
		int j;
		for (j = 0; j < n; j++) {
			if (ARY(cap, n, i, j) || ARY(cap, n, j, i)) {
				ARY(adj, n, i, deg[i]) = j;
				deg[i]++;
			}
		}
	}

	flow = 0;
	*fcost = 0;

	while (dijkstra(n, deg, adj, cap, fnet, cost, q, inq, par, pi, d, inf)) {
		int bot = INT_MAX;
		int u, v;

		for (v = n - 1, u = par[v]; v != 0; u = par[v]) {
			int f;

			if (ARY(fnet, n, v, u))
				f = ARY(fnet, n, v, u);
			else
				f = ARY(cap, n, u, v) - ARY(fnet, n, u, v);

			if (f < bot)
				bot = f;

			v = u;
		}

		for (v = n - 1, u = par[v]; v != 0; u = par[v]) {
			if (ARY(fnet, n, v, u)) {
				ARY(fnet, n, v, u) -= bot;
				(*fcost) -= bot * ARY(cost, n, v, u);
			} else {
				ARY(fnet, n, u, v) += bot;
				(*fcost) += bot * ARY(cost, n, u, v);
			}

			v = u;
		}

		flow += bot;
	}

	for (i = 0; i < n * n; i++)
		cap[i] = fnet[i];

	return flow;
}

void flow_hbondeng(struct atomgrp *ag, double *energy, struct nblist *nblst)
{
	int n_fedge = 0;

	FLOW_STRUCT *fs = (FLOW_STRUCT *) (ag->flow_struct);
	FLOW_EDGE *flow_edge = fs->flow_edge;

	double rc = nblst->nbcof;
	double rc2 = rc * rc;
	double min_en = 0;
	int i, j, k;
	int n_hydro;
	int n_acc;
	int n;
	int l;
	int *cap;
	FLOAT *cost;
	FLOAT inf;
	FLOAT fcost;
	int flow;
	FLOW_MAP *flow_map;

	*energy = 0;

	for (i = 0; i < nblst->nfat; i++) {
		int ai = nblst->ifat[i];
		mol_atom *atom_i = &(ag->atoms[ai]);
		int n2;
		int *p;

		if (!(atom_i->hprop & DONATABLE_HYDROGEN)
		    && !(atom_i->hprop & HBOND_ACCEPTOR))
			continue;

		n2 = nblst->nsat[i];
		p = nblst->isat[i];

		for (j = 0; j < n2; j++) {
			int aj = p[j];
			mol_atom *atom_j = &(ag->atoms[aj]);

			double en = 0;

			if (atom_i->hprop & DONATABLE_HYDROGEN) {
				if (!(atom_j->hprop & HBOND_ACCEPTOR))
					continue;

				en = get_pairwise_hbondeng_nblist(ag->atoms, ai,
								  ag->atoms, aj,
								  NULL, rc2, 0);
			} else {
				if (!(atom_j->hprop & DONATABLE_HYDROGEN))
					continue;

				en = get_pairwise_hbondeng_nblist(ag->atoms, aj,
								  ag->atoms, ai,
								  NULL, rc2, 0);
			}

			if (en < 0) {
				if (n_fedge >= fs->max_n_fedge) {
					if (!reinit_flow_struct
					    (fs, (n_fedge << 1), fs->max_n))
						return;
					flow_edge = fs->flow_edge;
				}

				if (atom_i->hprop & DONATABLE_HYDROGEN) {
					flow_edge[n_fedge].hydro_id = ai;
					flow_edge[n_fedge].acc_id = aj;
				} else {
					flow_edge[n_fedge].hydro_id = aj;
					flow_edge[n_fedge].acc_id = ai;
				}

				flow_edge[n_fedge].en = en;

				n_fedge++;

				if (en < min_en)
					min_en = en;
			}
		}
	}

	if (n_fedge == 0)
		return;

	for (i = 1; i < n_fedge; i++) {
		FLOW_EDGE e = flow_edge[i];
		j = i - 1;

		while ((j >= 0) && (e.hydro_id < flow_edge[j].hydro_id)) {
			flow_edge[j + 1] = flow_edge[j];
			j--;
		}

		flow_edge[j + 1] = e;
	}

	n_hydro = 1;

	for (i = 1; i < n_fedge; i++) {
		if (flow_edge[i].hydro_id > flow_edge[i - 1].hydro_id)
			n_hydro++;
	}

	if (fs->max_n < n_hydro + 1) {
		if (!reinit_flow_struct
		    (fs, fs->max_n_fedge, ((n_hydro + 1) << 1)))
			return;
	}

	flow_map = fs->flow_map;

	k = 0;

	for (i = 0; i < n_fedge; i++) {
		if (!i || (flow_edge[i].hydro_id > flow_map[k].id)) {
			flow_map[++k].cap = 1;
			flow_map[k].id = flow_edge[i].hydro_id;
		}

		flow_edge[i].hydro_id = k;
	}

	for (i = 1; i < n_fedge; i++) {
		FLOW_EDGE e = flow_edge[i];
		j = i - 1;

		while ((j >= 0) && (e.acc_id < flow_edge[j].acc_id)) {
			flow_edge[j + 1] = flow_edge[j];
			j--;
		}

		flow_edge[j + 1] = e;
	}

	n_acc = 1;

	for (i = 1; i < n_fedge; i++) {
		if (flow_edge[i].acc_id > flow_edge[i - 1].acc_id)
			n_acc++;
	}

	if (fs->max_n < n_hydro + n_acc + 2) {
		if (!reinit_flow_struct
		    (fs, fs->max_n_fedge, ((n_hydro + n_acc + 2) << 1)))
			return;
		flow_map = fs->flow_map;
	}

	for (i = 0; i < n_fedge; i++) {
		if (!i || (flow_edge[i].acc_id > flow_map[k].id)) {
			flow_map[++k].cap =
			    residual_acceptor_valency(&
						      (ag->
						       atoms[flow_edge[i].
							     acc_id]), ag->prm);
			flow_map[k].id = flow_edge[i].acc_id;
		}

		flow_edge[i].acc_id = k;
	}

	n = n_hydro + n_acc + 2;

	printf("n_fedge = %d, n_hydro = %d, n_acc = %d, n = %d\n", n_fedge,
	       n_hydro, n_acc, n);
	fflush(stdout);

	if ((n_fedge == n_hydro) && (n_fedge == n_acc)) {
		for (l = 0; l < n_fedge; l++) {
			int ai, aj;
			double en;
			i = flow_edge[l].hydro_id;
			j = flow_edge[l].acc_id;

			ai = flow_map[i].id;
			aj = flow_map[j].id;

			en =
			    get_pairwise_hbondeng_nblist(ag->atoms, ai,
							 ag->atoms, aj, NULL,
							 rc2, 1);

			(*energy) += en;
		}

		return;
	}

	cap = fs->cap;
	cost = fs->cost;

	for (i = 0; i < n * n; i++) {
		cap[i] = 0;
		cost[i] = 0;
	}

	for (i = 1; i <= n_hydro; i++)
		ARY(cap, n, 0, i) = 1;

	for (i = n_hydro + 1; i < n - 1; i++)
		ARY(cap, n, i, n - 1) = flow_map[i].cap;

	inf = 1;

	for (l = 0; l < n_fedge; l++) {
		i = flow_edge[l].hydro_id;
		j = flow_edge[l].acc_id;

		ARY(cap, n, i, j) = 1;
		ARY(cost, n, i, j) = flow_edge[l].en - (min_en - 1);

		inf += ARY(cost, n, i, j);
	}

	inf *= 5.0;

	flow = min_cost_max_flow(n, fs, inf, &fcost);

	if (flow > 0) {
		for (l = 0; l < n_fedge; l++) {
			int ai, aj, f;
			double en;
			i = flow_edge[l].hydro_id;
			j = flow_edge[l].acc_id;

			ai = flow_map[i].id;
			aj = flow_map[j].id;

			f = ARY(cap, n, i, j) - ARY(cap, n, j, i);

			if (f < 1)
				continue;

			en =
			    get_pairwise_hbondeng_nblist(ag->atoms, ai,
							 ag->atoms, aj, NULL,
							 rc2, 1);

			(*energy) += en;
		}
	}

	printf("flow = %d, fcost = %lf, *energy = %lf\n", flow,
	       (double)(fcost + (flow * (min_en - 1))), (double)(*energy));
	fflush(stdout);
}

/* --------------------------------------------------------------------------

Hydrogen Bonding Energy Functions in presence of Small Molecules

Receptor - Small Molecule Complex

By: Mohammad Moghadasi 		mohamad@bu.edu

March 2013

-------------------------------------------------------------------------- */
static double get_pairwise_hbondeng_nblist_smallmol(mol_atom * atoms_hydro,
						    int hydro_id,
						    mol_atom * atoms_acc,
						    int acc_id, double *engcat,
						    double rc2, int comp_grad,
						    struct rig_forest *prot);
static int get_hbe_type_smallmol(mol_atom * atoms, mol_atom * hydro,
				 mol_atom * acc, struct rig_forest *prot);
static int get_donor_chem_type_smallmol(mol_atom * don,
					struct rig_forest *prot);
static int get_acceptor_chem_type_smallmol(mol_atom * acc,
					   struct rig_forest *prot);
static int get_hbe_type_smallmol(mol_atom * atoms, mol_atom * hydro,
				 mol_atom * acc, struct rig_forest *prot);
int hbond_energy_and_gradient_computation_smallmol(int hbe, mol_atom * atoms_hydro,	// atom group containing the hydrogen atom
						   int hydro_id,	// index of the hydrogen atom in atoms_hydro
						   mol_atom * atoms_acc,	// atom group containing the acceptor atom
						   int acc_id,	// index of the acceptor atom in atoms_acc
						   double *energy,
						   int comp_grad);
static int get_hbond_weight_type_smallmol(int hbe_type);	// DONE
int hbond_energy_computation_smallmol(int hbe,	// DONE
				      FLOAT AHdis,
				      FLOAT xD,
				      FLOAT xH,
				      FLOAT * energy,
				      FLOAT * dE_dr,
				      FLOAT * dE_dxD, FLOAT * dE_dxH);
int create_base_to_acceptor_unit_vector_smallmol(int hbe_type, mol_atom * atoms, mol_atom * acc, struct dvector *B, struct dvector *BAunit, FLOAT * invBAdis);	// DONE

void hbondengcat_smallmol(struct atomgrp *ag, double *energy,
			  struct nblist *nblst, struct rig_forest *prot)
{

	/*
	   // Display Small Molecule Atoms
	   printf("Number of Trees = %d\n", prot->numt);
	   for(int i = 0 ; i < prot->numt; i++){
	   printf("Size of Tree[%d] = %d\n",i, prot->tr[i]->total_atm_size);
	   for(int j = 0 ; j < prot->tr[i]->total_atm_size; j++){
	   printf("Atom Index[%d] = %d\n",j, prot->tr[i]->tree_atoms[j]);
	   }
	   }
	 */

	double rc = nblst->nbcof;
	double rc2 = rc * rc;
	int i;

	for (i = 0; i < nblst->nfat; i++) {
		int ai = nblst->ifat[i];
		mol_atom *atom_i = &(ag->atoms[ai]);
		int j;
		int n2;
		int *p;

		if (!(atom_i->hprop & DONATABLE_HYDROGEN)
		    && !(atom_i->hprop & HBOND_ACCEPTOR))
			continue;

		n2 = nblst->nsat[i];
		p = nblst->isat[i];

		for (j = 0; j < n2; j++) {
			int aj = p[j];
			mol_atom *atom_j = &(ag->atoms[aj]);

			if (atom_i->hprop & DONATABLE_HYDROGEN) {
				if (!(atom_j->hprop & HBOND_ACCEPTOR))
					continue;

				get_pairwise_hbondeng_nblist_smallmol(ag->atoms,
								      ai,
								      ag->atoms,
								      aj,
								      energy,
								      rc2, 1,
								      prot);
			} else {
				if (!(atom_j->hprop & DONATABLE_HYDROGEN))
					continue;

				get_pairwise_hbondeng_nblist_smallmol(ag->atoms,
								      aj,
								      ag->atoms,
								      ai,
								      energy,
								      rc2, 1,
								      prot);
			}
		}
	}

	energy[hbw_TOTAL] =
	    energy[hbw_SR_BB] + energy[hbw_LR_BB] + energy[hbw_BB_SC] +
	    energy[hbw_SC] + energy[hbw_BB_SM] + energy[hbw_SC_SM] +
	    energy[hbw_SM];
}

static double get_pairwise_hbondeng_nblist_smallmol(mol_atom * atoms_hydro,
						    int hydro_id,
						    mol_atom * atoms_acc,
						    int acc_id, double *engcat,
						    double rc2, int comp_grad,
						    struct rig_forest *prot)
{
	mol_atom *hydro = &(atoms_hydro[hydro_id]);
	mol_atom *acc = &(atoms_acc[acc_id]);

	double dx = hydro->X - acc->X;
	double dy = hydro->Y - acc->Y;
	double dz = hydro->Z - acc->Z;

	double d2 = dx * dx + dy * dy + dz * dz;
	int hbe;
	double en = 0;

	if (d2 > rc2)
		return 0;

	hbe = get_hbe_type_smallmol(atoms_hydro, hydro, acc, prot);
/*
    if(hbe >= 12)	// Only SM-involved
	printf("acc-id: %d  don-id: %d     hbe-type: %d\n", acc->ingrp, hydro->ingrp,hbe);
*/

	if (hbe == hbe_NONE)
		return 0;

	if (hbond_energy_and_gradient_computation_smallmol
	    (hbe, atoms_hydro, hydro_id, atoms_acc, acc_id, &en, comp_grad)) {
		if (engcat != NULL) {
			int hbw = get_hbond_weight_type_smallmol(hbe);
/*
 	    if(hbe >= 12)       // Only SM-involved	TO TEST
                printf("acc-id: %d  don-id: %d     hbe-type: %d	en: %f  +++++++++++++++++ \n", acc->ingrp, hydro->ingrp,hbe,en);
*/
			engcat[hbw] += en;
		}

		return en;
	}
	return 0;
}

static int get_hbe_type_smallmol(mol_atom * atoms, mol_atom * hydro,
				 mol_atom * acc, struct rig_forest *prot)
{
	mol_atom *don;
	int don_type;
	int acc_type;

	if ((hydro->base < 0) && (!(acc->hprop & HBOND_ACCEPTOR)))
		return hbe_NONE;

	don = &(atoms[hydro->base]);

	don_type = get_donor_chem_type_smallmol(don, prot);
	acc_type = get_acceptor_chem_type_smallmol(acc, prot);

	if ((don_type == hbdon_BB) && (acc_type == hbacc_BB)) { // acc BB  don BB
		return classify_BB_by_separation(don, acc);

	} else if (acc_type == hbacc_BB) { // acc BB  don S*
		switch (don_type) {
		case hbdon_SC:
			return hbe_BSC; // acc BB  don SC
		case hbdon_SM:
			return hbe_BSM; // acc BB  don SM
		} 
	} else if (don_type == hbdon_BB) { // acc S*  don BB
		switch (acc_type) {
		case hbacc_SP2:
			return hbe_SP2B; // acc SC-SP2   don BB
		case hbacc_SP3:
			return hbe_SP3B; // acc SC-SP3   don BB
		case hbacc_RING:
			return hbe_RINGB; // acc SC-RING  don BB

		case hbacc_SP2_SM:
			return hbe_SP2SMB; // acc SM-SP2   don BB
		case hbacc_SP3_SM:
			return hbe_SP3SMB; // acc SM-SP3   don BB
		case hbacc_RING_SM:
			return hbe_RINGSMB; // acc SM-RING  don BB
		}
	} else if (don_type == hbdon_SC) { // acc SM  don SC
		switch (acc_type) {
		case hbacc_SP2:
			return hbe_SP2SC; // acc SC-SP2   don SC
		case hbacc_SP3:
			return hbe_SP3SC; // acc SC-SP3   don SC
		case hbacc_RING:
			return hbe_RINGSC; // acc SC-RING  don SC

		case hbacc_SP2_SM:
			return hbe_SP2SMSC; // acc SM-SP2   don SC
		case hbacc_SP3_SM:
			return hbe_SP3SMSC; // acc SM-SP3   don SC
		case hbacc_RING_SM:
			return hbe_RINGSMSC; // acc SM-RING  don SC
		}
	} else if (don_type == hbdon_SM) { // acc SM  don SC
		switch (acc_type) {
		case hbacc_SP2:
			return hbe_SP2SCSM; // acc SC-SP2   don SM
		case hbacc_SP3:
			return hbe_SP3SCSM; // acc SC-SP3   don SM
		case hbacc_RING:
			return hbe_RINGSCSM; // acc SC-RING  don SM

		case hbacc_SP2_SM:
			return hbe_SP2SM; // acc SM-SP2   don SM
		case hbacc_SP3_SM:
			return hbe_SP3SM; // acc SM-SP3   don SM
		case hbacc_RING_SM:
			return hbe_RINGSM; // acc SM-RING  don SM
		}
	}
	return hbe_NONE;
}

static int get_donor_chem_type_smallmol(mol_atom * don, struct rig_forest *prot)
{
	int i;
	int j;
	if (don->backbone)
		return hbdon_BB;

	for (i = 0; i < prot->numt; i++) { // number of trees (small molecules)
		for (j = 0; j < prot->tr[i]->total_atm_size; j++) // number of leaves (atoms) of each tree
			if (don->ingrp == prot->tr[i]->tree_atoms[j]) {
				return hbdon_SM;
			}
	}
	return hbdon_SC;
}

static int get_acceptor_chem_type_smallmol(mol_atom * acc,
					   struct rig_forest *prot)
{
	if (acc->backbone)
		return hbacc_BB;
	else {
		int smallmol_flag = 0;
		int i, j;
		for (i = 0; i < prot->numt; i++) {	// number of trees (small molecules)
			for (j = 0; j < prot->tr[i]->total_atm_size; j++)	// number of leaves (atoms) of each tree
				if (acc->ingrp == prot->tr[i]->tree_atoms[j]) {
					smallmol_flag = 1;
					break;
				}
			if (smallmol_flag)
				break;
		}
		if (smallmol_flag) { // small molecule
			switch (acc->hybridization) {
			case SP2_HYBRID:
				return hbacc_SP2_SM;
			case SP3_HYBRID:
				return hbacc_SP3_SM;
			case RING_HYBRID:
				return hbacc_RING_SM;
			default:
				break;
			}
		} else { // side chain
			switch (acc->hybridization) {
			case SP2_HYBRID:
				return hbacc_SP2;
			case SP3_HYBRID:
				return hbacc_SP3;
			case RING_HYBRID:
				return hbacc_RING;
			default:
				break;
			}
		}
	}
	return hbacc_NO;
}

int hbond_energy_and_gradient_computation_smallmol(int hbe, mol_atom * atoms_hydro,	// atom group containing the hydrogen atom
						   int hydro_id,	// index of the hydrogen atom in atoms_hydro
						   mol_atom * atoms_acc,	// atom group containing the acceptor atom
						   int acc_id,	// index of the acceptor atom in atoms_acc
						   double *energy,
						   int comp_grad)
{
	struct dvector AHunit, HDunit, BAunit, B;
	FLOAT invAHdis, invHDdis, invBAdis;

	mol_atom *hydro = &(atoms_hydro[hydro_id]);
	mol_atom *acc = &(atoms_acc[acc_id]);
	mol_atom *don = &(atoms_hydro[hydro->base]);
	FLOAT AHdis;
	FLOAT dE_dr, dE_dxD, dE_dxH;
	FLOAT en;
	FLOAT xD, xH;
	int en_comp;
	mol_atom *base;
	mol_atom *base2;
	FLOAT u;
	FLOAT v;

	if (!create_donor_orientation_unit_vector
	    (atoms_hydro, hydro, &HDunit, &invHDdis))
		return 0;

	if ((invHDdis < 0.8) || (invHDdis > 1.25)) {
		*energy = 0;
		return 1;
	}

	if (!create_base_to_acceptor_unit_vector_smallmol
	    (hbe, atoms_acc, acc, &B, &BAunit, &invBAdis))
		return 0;

	AHunit.X = hydro->X - acc->X;
	AHunit.Y = hydro->Y - acc->Y;
	AHunit.Z = hydro->Z - acc->Z;

#ifdef USE_LONG_DOUBLE
	AHdis =
	    sqrtl(AHunit.X * AHunit.X + AHunit.Y * AHunit.Y +
		  AHunit.Z * AHunit.Z);
#else
	AHdis =
	    sqrt(AHunit.X * AHunit.X + AHunit.Y * AHunit.Y +
		 AHunit.Z * AHunit.Z);
#endif

	if (AHdis <= 0) {
		print_error("Overlapping acceptor and hydrogen!");
		return 0;
	}

	invAHdis = 1 / AHdis;

	AHunit.X *= invAHdis;
	AHunit.Y *= invAHdis;
	AHunit.Z *= invAHdis;

	xD = dot_product(&AHunit, &HDunit), xH = dot_product(&BAunit, &AHunit);

	en_comp =
	    hbond_energy_computation_smallmol(hbe, AHdis, xD, xH, &en, &dE_dr,
					      &dE_dxD, &dE_dxH);

	if (!en_comp)
		return 0;

	*energy = en;

	if (*energy >= HB_ENG_MAX) {
		*energy = 0;
		return 1;
	}

	if (!comp_grad)
		return 1;

	base = &(atoms_acc[acc->base]);
	base2 = NULL;

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB) || (hbe == hbe_RINGSCSM)
	    || (hbe == hbe_RINGSMSC) || (hbe == hbe_RINGSMB)
	    || (hbe == hbe_RINGSM) || (hbe == hbe_SP3SC) || (hbe == hbe_SP3B)
	    || (hbe == hbe_SP3SCSM) || (hbe == hbe_SP3SMSC)
	    || (hbe == hbe_SP3SMB) || (hbe == hbe_SP3SM)) {
		if (acc->base2 < 0) {
			print_error("Acceptor base 2 missing!");
			return 0;
		}
		base2 = &(atoms_acc[acc->base2]);
	}

	u = dE_dr * AHunit.X;
	hydro->GX += -u;
	acc->GX += u;

	u = dE_dr * AHunit.Y;
	hydro->GY += -u;
	acc->GY += u;

	u = dE_dr * AHunit.Z;
	hydro->GZ += -u;
	acc->GZ += u;

	u = -dE_dxD * invAHdis * (xD * AHunit.X - HDunit.X);
	v = -dE_dxD * invHDdis * (AHunit.X - xD * HDunit.X);
	hydro->GX += -(u + v);
	acc->GX += u;
	don->GX += v;

	u = -dE_dxD * invAHdis * (xD * AHunit.Y - HDunit.Y);
	v = -dE_dxD * invHDdis * (AHunit.Y - xD * HDunit.Y);
	hydro->GY += -(u + v);
	acc->GY += u;
	don->GY += v;

	u = -dE_dxD * invAHdis * (xD * AHunit.Z - HDunit.Z);
	v = -dE_dxD * invHDdis * (AHunit.Z - xD * HDunit.Z);
	hydro->GZ += -(u + v);
	acc->GZ += u;
	don->GZ += v;

	u = -dE_dxH * invAHdis * (BAunit.X - xH * AHunit.X);
	v = -dE_dxH * invBAdis * (xH * BAunit.X - AHunit.X);
	hydro->GX += u;
	acc->GX += -(u + v);

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB) || (hbe == hbe_RINGSCSM)
	    || (hbe == hbe_RINGSMSC) || (hbe == hbe_RINGSMB)
	    || (hbe == hbe_RINGSM)) {
		base->GX += 0.5 * v;
		base2->GX += 0.5 * v;
	} else if ((hbe == hbe_SP3SC) || (hbe == hbe_SP3B)
		   || (hbe == hbe_SP3SCSM) || (hbe == hbe_SP3SMSC)
		   || (hbe == hbe_SP3SMB) || (hbe == hbe_SP3SM))
		base2->GX += v;
	else
		base->GX += v;

	u = -dE_dxH * invAHdis * (BAunit.Y - xH * AHunit.Y);
	v = -dE_dxH * invBAdis * (xH * BAunit.Y - AHunit.Y);
	hydro->GY += u;
	acc->GY += -(u + v);

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB) || (hbe == hbe_RINGSCSM)
	    || (hbe == hbe_RINGSMSC) || (hbe == hbe_RINGSMB)
	    || (hbe == hbe_RINGSM)) {
		base->GY += 0.5 * v;
		base2->GY += 0.5 * v;
	} else if ((hbe == hbe_SP3SC) || (hbe == hbe_SP3B)
		   || (hbe == hbe_SP3SCSM) || (hbe == hbe_SP3SMSC)
		   || (hbe == hbe_SP3SMB) || (hbe == hbe_SP3SM))
		base2->GY += v;
	else
		base->GY += v;

	u = -dE_dxH * invAHdis * (BAunit.Z - xH * AHunit.Z);
	v = -dE_dxH * invBAdis * (xH * BAunit.Z - AHunit.Z);
	hydro->GZ += u;
	acc->GZ += -(u + v);

	if ((hbe == hbe_RINGSC) || (hbe == hbe_RINGB) || (hbe == hbe_RINGSCSM)
	    || (hbe == hbe_RINGSMSC) || (hbe == hbe_RINGSMB)
	    || (hbe == hbe_RINGSM)) {
		base->GZ += 0.5 * v;
		base2->GZ += 0.5 * v;
	} else if ((hbe == hbe_SP3SC) || (hbe == hbe_SP3B)
		   || (hbe == hbe_SP3SCSM) || (hbe == hbe_SP3SMSC)
		   || (hbe == hbe_SP3SMB) || (hbe == hbe_SP3SM))
		base2->GZ += v;
	else
		base->GZ += v;

	return 1;
}

// To construct the unit vector from acceptor base to the acceptor atom:

int create_base_to_acceptor_unit_vector_smallmol(int hbe_type, mol_atom * atoms, mol_atom * acc, struct dvector *B, struct dvector *BAunit, FLOAT * invBAdis)	// DONE
{
	mol_atom *base, *base2;
	FLOAT BAdis;

	if (acc->base < 0)
		return 0;
	else
		base = &(atoms[acc->base]);

	if ((hbe_type == hbe_RINGSC) || (hbe_type == hbe_RINGB)
	    || (hbe_type == hbe_RINGSCSM) || (hbe_type == hbe_RINGSMSC)
	    || (hbe_type == hbe_RINGSMB) || (hbe_type == hbe_RINGSM)
	    || (hbe_type == hbe_SP3SC) || (hbe_type == hbe_SP3B)
	    || (hbe_type == hbe_SP3SCSM) || (hbe_type == hbe_SP3SMSC)
	    || (hbe_type == hbe_SP3SMB) || (hbe_type == hbe_SP3SM)) {
		if (acc->base2 < 0)
			return 0;
		else
			base2 = &(atoms[acc->base2]);
	}

	switch (hbe_type) {
	case hbe_RINGSC:
	case hbe_RINGSM:
	case hbe_RINGSMSC:
	case hbe_RINGSCSM:
	case hbe_RINGB:
	case hbe_RINGSMB:	// RING: B is the midpoint between base and base2
		B->X = (base->X + base2->X) * 0.5;
		B->Y = (base->Y + base2->Y) * 0.5;
		B->Z = (base->Z + base2->Z) * 0.5;
		break;

	case hbe_SP3SC:
	case hbe_SP3SM:
	case hbe_SP3SCSM:
	case hbe_SP3SMSC:
	case hbe_SP3B:
	case hbe_SP3SMB:	// SP3: B is base2
		B->X = base2->X;
		B->Y = base2->Y;
		B->Z = base2->Z;
		break;

	default:		// SP2: B is base
		B->X = base->X;
		B->Y = base->Y;
		B->Z = base->Z;
		break;
	}

	BAunit->X = acc->X - B->X;
	BAunit->Y = acc->Y - B->Y;
	BAunit->Z = acc->Z - B->Z;

	BAdis =
	    BAunit->X * BAunit->X + BAunit->Y * BAunit->Y +
	    BAunit->Z * BAunit->Z;

	if (BAdis <= 0) {
		print_error("Failed to normalize the base to acceptor vector!");
		return 0;
	}
#ifdef USE_LONG_DOUBLE
	BAdis = sqrtl(BAdis);
#else
	BAdis = sqrt(BAdis);
#endif
	*invBAdis = 1 / BAdis;

	BAunit->X *= (*invBAdis);
	BAunit->Y *= (*invBAdis);
	BAunit->Z *= (*invBAdis);

	return 1;
}

int hbond_energy_computation_smallmol(int hbe,	// DONE
				      FLOAT AHdis,
				      FLOAT xD,
				      FLOAT xH,
				      FLOAT * energy,
				      FLOAT * dE_dr,
				      FLOAT * dE_dxD, FLOAT * dE_dxH)
{
	FLOAT dAHdis;
	FLOAT dxD;
	FLOAT dxH;

	FLOAT FSr = 0.0, FLr = 0.0, FxD = 0.0, FxH = 0.0;	// fading intervals values
	FLOAT dFSr = 0.0, dFLr = 0.0, dFxD = 0.0, dFxH = 0.0;	// fading intervals derivatives
	FLOAT Pr = 0.0, PSxD = 0.0, PSxH = 0.0, PLxD = 0.0, PLxH = 0.0;	// polynomials values
	FLOAT dPr = 0.0, dPSxD = 0.0, dPSxH = 0.0, dPLxD = 0.0, dPLxH = 0.0;	// polynomials derivatives
	if (energy != NULL)
		*energy = HB_ENG_MAX + 1.0f;
	if (dE_dr != NULL)
		*dE_dr = 0.0;
	if (dE_dxD != NULL)
		*dE_dxD = 0.0;
	if (dE_dxH != NULL)
		*dE_dxH = 0.0;

#ifdef USE_LONG_DOUBLE
	if ((fabsl(xD) > 1.0) || (fabsl(xH) > 1.0)) {
		print_error
		    ("Invalid angle value in hbond_energy_computation: xH = %Lf, xD = %Lf\n",
		     xH, xD);
		return 0;
	}
#else
	if ((fabs(xD) > 1.0) || (fabs(xH) > 1.0)) {
		print_error
		    ("Invalid angle value in hbond_energy_computation: xH = %lf, xD = %lf\n",
		     xH, xD);
		return 0;
	}
#endif

	if ((AHdis > MAX_AH) || (AHdis < MIN_AH) || (xH < MIN_xH)
	    || (xD < MIN_xD) || (xH > MAX_xH) || (xD > MAX_xD))
		return 0;

	dAHdis = AHdis;
	dxD = xD;
	dxH = xH;

// fade functions

#if FADING_FUNCTION == LINEAR_FADE
	fade_xD_value_deriv(xD, &FxD, &dFxD);
	fade_xH_value_deriv(xH, &FxH, &dFxH);

	if (hbe_is_BB_type(hbe))
		fade_rBB_value_deriv(AHdis, &FSr, &dFSr);
	else {
		fade_rshort_value_deriv(AHdis, &FSr, &dFSr);
		fade_rlong_value_deriv(AHdis, &FLr, &dFLr);
	}
#elif FADING_FUNCTION == LOG_FADE
	fade_xD_log_value_deriv(xD + 1, &FxD, &dFxD);
	FxD -= 0.1;
	fade_xH_log_value_deriv(xH + 1, &FxH, &dFxH);
	FxH -= 0.1;

	if (hbe_is_BB_type(hbe))
		fade_rBB_log_value_deriv(AHdis, &FSr, &dFSr);
	else {
		fade_rshort_log_value_deriv(AHdis, &FSr, &dFSr);
		fade_rlong_log_value_deriv(AHdis, &FLr, &dFLr);
	}
#elif FADING_FUNCTION == BSPLINE_FADE
	fade_xD_bspline_value_deriv(xD /* + 1 */ , &FxD, &dFxD);	//FxD -= 0.1;
	fade_xH_bspline_value_deriv(xH /* + 1 */ , &FxH, &dFxH);	//FxH -= 0.1;

	if (hbe_is_BB_type(hbe))
		fade_rBB_bspline_value_deriv(AHdis, &FSr, &dFSr);
	else {
		fade_rshort_bspline_value_deriv(AHdis, &FSr, &dFSr);
		fade_rlong_bspline_value_deriv(AHdis, &FLr, &dFLr);
	}
#endif
	switch (hbe) {
	case hbe_NONE:
		break;

	case hbe_BBHELIX:	// hbstat = 1, polynomials 1,5,11
		AH_BBHelix(dAHdis, &Pr, &dPr);
		xD_BBHelix(dxD, &PSxD, &dPSxD);
		xH_BBHelix(dxH, &PSxH, &dPSxH);
		break;

	case hbe_BBOTHER:	// hbstat = 2, polynomials 2,6,11
		AH_BBOther(dAHdis, &Pr, &dPr);
		xD_BBOther(dxD, &PSxD, &dPSxD);
		xH_BBOther(dxH, &PSxH, &dPSxH);
		break;

	case hbe_BB:		// hbstat = 3, polynomials 3,7,13,8,14
	case hbe_BBTURN:
	case hbe_BSC:
	case hbe_SP2B:
	case hbe_SP2SC:
	case hbe_BSM:
	case hbe_SP2SMB:
	case hbe_SP2SMSC:
	case hbe_SP2SM:
	case hbe_SP2SCSM:
		AH_SP2(dAHdis, &Pr, &dPr);
		xD_SP2short(dxD, &PSxD, &dPSxD);
		xH_SP2short(dxH, &PSxH, &dPSxH);
		xD_SP2long(dxD, &PLxD, &dPLxD);
		xH_SP2long(dxH, &PLxH, &dPLxH);
		break;

	case hbe_SP3B:		// hbstat = 4, polynomials 4,9,15,10,15
	case hbe_SP3SMB:	// hbstat = 4, polynomials 4,9,15,10,15
	case hbe_SP3SC:
	case hbe_SP3SCSM:
	case hbe_SP3SMSC:
	case hbe_SP3SM:
		AH_SP3(dAHdis, &Pr, &dPr);
		xD_SP3short(dxD, &PSxD, &dPSxD);
		xH_SP3(dxH, &PSxH, &dPSxH);
		xD_SP3long(dxD, &PLxD, &dPLxD);
		PLxH = PSxH;
		dPLxH = dPSxH;
		break;

	case hbe_RINGB:	// hbstat = 5, polynomials 3,7,16,8,16
	case hbe_RINGSMB:	// hbstat = 5, polynomials 3,7,16,8,16
	case hbe_RINGSC:
	case hbe_RINGSCSM:
	case hbe_RINGSMSC:
	case hbe_RINGSM:
		AH_SP2(dAHdis, &Pr, &dPr);
		xD_SP2short(dxD, &PSxD, &dPSxD);
		xH_Ring(dxH, &PSxH, &dPSxH);
		xD_SP2long(dxD, &PLxD, &dPLxD);
		PLxH = PSxH;
		dPLxH = dPSxH;
		break;
	}

	if (energy != NULL)
		*energy =
		    Pr * FxD * FxH + FSr * (PSxD * FxH + FxD * PSxH) +
		    FLr * (PLxD * FxH + FxD * PLxH);
	if (dE_dr != NULL)
		*dE_dr =
		    dPr * FxD * FxH + dFSr * (PSxD * FxH + FxD * PSxH) +
		    dFLr * (PLxD * FxH + FxD * PLxH);
	if (dE_dxD != NULL)
		*dE_dxD =
		    dFxD * (Pr * FxH + FLr * PLxH + FSr * PSxH) +
		    FxH * (FSr * dPSxD + FLr * dPLxD);
	if (dE_dxH != NULL)
		*dE_dxH =
		    dFxH * (Pr * FxD + FLr * PLxD + FSr * PSxD) +
		    FxD * (FSr * dPSxH + FLr * dPLxH);

	return 1;
}

static int get_hbond_weight_type_smallmol(int hbe_type)	// DONE
{
	switch (hbe_type) {
	case hbe_BBTURN:
	case hbe_BBHELIX:
		return hbw_SR_BB;
	case hbe_BBOTHER:
		return hbw_LR_BB;
	case hbe_SP2B:
	case hbe_SP3B:
	case hbe_RINGB:
	case hbe_BSC:
		return hbw_BB_SC;
	case hbe_SP2SMB:
	case hbe_SP3SMB:
	case hbe_RINGSMB:
	case hbe_BSM:
		return hbw_BB_SM;
	case hbe_SP2SC:
	case hbe_SP3SC:
	case hbe_RINGSC:
		return hbw_SC;
	case hbe_SP2SM:
	case hbe_SP3SM:
	case hbe_RINGSM:
		return hbw_SM;
	case hbe_SP2SCSM:
	case hbe_SP2SMSC:
	case hbe_SP3SCSM:
	case hbe_SP3SMSC:
	case hbe_RINGSCSM:
	case hbe_RINGSMSC:
		return hbw_SC_SM;
	default:
		return hbw_NONE;
	}
	return hbw_NONE;
}

void get_categorized_hbondeng_smallmol(double *engcat, double *bb_bb_sr,
				       double *bb_bb_lr, double *bb_sc,
				       double *sc_sc, double *bb_sm,
				       double *sc_sm, double *sm_sm)
{
	*bb_bb_sr = engcat[hbw_SR_BB];
	*bb_bb_lr = engcat[hbw_LR_BB];
	*bb_sc = engcat[hbw_BB_SC];
	*sc_sc = engcat[hbw_SC];
	*bb_sm = engcat[hbw_BB_SM];
	*sc_sm = engcat[hbw_SC_SM];
	*sm_sm = engcat[hbw_SM];

}

void set_categorized_hbondeng_smallmol(double *engcat, double bb_bb_sr,
				       double bb_bb_lr, double bb_sc,
				       double sc_sc, double bb_sm, double sc_sm,
				       double sm_sm)
{
	engcat[hbw_SR_BB] = bb_bb_sr;
	engcat[hbw_LR_BB] = bb_bb_lr;
	engcat[hbw_BB_SC] = bb_sc;
	engcat[hbw_SC] = sc_sc;
	engcat[hbw_BB_SM] = bb_sm;
	engcat[hbw_SC_SM] = sc_sm;
	engcat[hbw_SM] = sm_sm;
}

void init_categorized_hbondeng_smallmol(double *engcat)
{
	int i;
	if (engcat == NULL)
		return;

	for (i = 0; i <= hbw_SM; i++)
		engcat[i] = 0.0;
}

double *alloc_categorized_hbondeng_smallmol(void)
{
	return (double *)_mol_malloc((hbw_SM + 1) * sizeof(double));
}
