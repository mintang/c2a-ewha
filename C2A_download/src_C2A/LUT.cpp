/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             E. Larsen
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "LUT.h"
#include "MatVec.h"
#include "RectDist.h"
#include "OBB_Disjoint.h"
#include "PQP.h"



static const int north[3] = {112, 480, 1984};
static const int south[3] = {113, 481, 1985};
static const int bases[3] = {0, 4, 12};
static const PQP_REAL aincr[3] = { 0.125*M_PI, 0.0625*M_PI, 0.03125*M_PI };
static const PQP_REAL tan_xpi_y[28] = {
    tan( 0.0625*M_PI ), tan( 0.1875*M_PI ), tan( 0.3125*M_PI ), tan( 0.4375*M_PI ),

    tan( 0.03125*M_PI ), tan( 0.09375*M_PI ), tan( 0.15625*M_PI ),
    tan( 0.21875*M_PI ), tan( 0.28125*M_PI ), tan( 0.34375*M_PI ),
    tan( 0.40625*M_PI ), tan( 0.46875*M_PI ),

    tan( 0.015625*M_PI ), tan( 0.046875*M_PI ), tan( 0.078125*M_PI ),
    tan( 0.109375*M_PI ), tan( 0.140625*M_PI ), tan( 0.171875*M_PI ),
    tan( 0.203125*M_PI ), tan( 0.234375*M_PI ), tan( 0.265625*M_PI ),
    tan( 0.296875*M_PI ), tan( 0.328125*M_PI ), tan( 0.359375*M_PI ),
    tan( 0.390625*M_PI ), tan( 0.421875*M_PI ), tan( 0.453125*M_PI ),
    tan( 0.484375*M_PI ) };
static const PQP_REAL tan_sq_xpi_y[28] = {
  tan_xpi_y[0] * tan_xpi_y[0], tan_xpi_y[1] * tan_xpi_y[1],
  tan_xpi_y[2] * tan_xpi_y[2], tan_xpi_y[3] * tan_xpi_y[3],

  tan_xpi_y[4] * tan_xpi_y[4], tan_xpi_y[5] * tan_xpi_y[5],
  tan_xpi_y[6] * tan_xpi_y[6], tan_xpi_y[7] * tan_xpi_y[7],
  tan_xpi_y[8] * tan_xpi_y[8], tan_xpi_y[9] * tan_xpi_y[9],
  tan_xpi_y[10] * tan_xpi_y[10], tan_xpi_y[11] * tan_xpi_y[11],

  tan_xpi_y[12] * tan_xpi_y[12], tan_xpi_y[13] * tan_xpi_y[13],
  tan_xpi_y[14] * tan_xpi_y[14], tan_xpi_y[15] * tan_xpi_y[15],
  tan_xpi_y[16] * tan_xpi_y[16], tan_xpi_y[17] * tan_xpi_y[17],
  tan_xpi_y[18] * tan_xpi_y[18], tan_xpi_y[19] * tan_xpi_y[19],
  tan_xpi_y[20] * tan_xpi_y[20], tan_xpi_y[21] * tan_xpi_y[21],
  tan_xpi_y[22] * tan_xpi_y[22], tan_xpi_y[23] * tan_xpi_y[23],
  tan_xpi_y[24] * tan_xpi_y[24], tan_xpi_y[25] * tan_xpi_y[25],
  tan_xpi_y[26] * tan_xpi_y[26], tan_xpi_y[27] * tan_xpi_y[27] };

PQP_VECTOR3** Compute_Sphere_Pts( )
{
    PQP_VECTOR3** result = new PQP_VECTOR3*[3];
    int i, j, idx, k;

    for( idx = 0; idx < 3; idx++ ) {
        PQP_VECTOR3* sub_result = new PQP_VECTOR3[south[idx]+1];
        // Sweep from north to south
        for( i = (4<<idx)-1; i >= 1-(4<<idx); i-- ) {
            const PQP_REAL al = (PQP_REAL)i * aincr[idx];
            const PQP_REAL cos_al = cos( al );
            const PQP_REAL sin_al = sin( al );
            for( j = 0; j < (16<<idx); j++ ) {
                const PQP_REAL az = (PQP_REAL)j * aincr[idx];
								k = (((4<<idx)-1-i)<<(4+idx))+j;
								sub_result[ k ][0]= cos_al * cos(az);
								sub_result[ k ][1]= cos_al * sin(az);
								sub_result[ k ][2]= sin_al;
            }
        }

        sub_result[north[idx]][0]=0.0;
				sub_result[north[idx]][1]=0.0;
				sub_result[north[idx]][2]=1.0;
				sub_result[south[idx]][0]=0.0;
				sub_result[south[idx]][1]=0.0;
				sub_result[south[idx]][2]=0.0;

        result[idx] = sub_result;
    }

    return result;
}

PQP_VECTOR3** sphere_pts = Compute_Sphere_Pts();


void LUT::Create( PQP_Model* p )
{
    // Compute the lookup table.  Find object points that are nearest to points
    // on the sphere of a radius a factor bigger than the object radius.
    int i;
    PQP_REAL lut_radius = 1.0;

    type = LUT_22_5;

		lut  = new PQP_REAL[ south[(int)type]+1 ];

		PQP_REAL c[3];
    for( i = 0; i < north[(int)type]; i++ ) {
				c[0]=sphere_pts[(int)type][i][0] + p->com[0];
				c[1]=sphere_pts[(int)type][i][1] + p->com[1];
				c[2]=sphere_pts[(int)type][i][2] + p->com[2];
        lut[i] = p->Max_Cross_Product( c );
    }

		c[0]= p->com[0];
		c[1]= p->com[1];
		c[2]= lut_radius + p->com[2];
    lut[north[(int)type]] = p->Max_Cross_Product( c );

		c[0]= p->com[0];
		c[1]= p->com[1];
		c[2]= -lut_radius + p->com[2];
    lut[south[(int)type]] = p->Max_Cross_Product( c );
}


PQP_REAL LUT::Lookup_Internal( const PQP_REAL dir[3] )
{
    const PQP_REAL xy_sq = dir[0] * dir[0] + dir[1] * dir[1];
    const PQP_REAL z_sq = dir[2] * dir[2];
    const int idx = (int)type;
    const int tan_len_shift = 4 + idx;
    const int tan_len = 4 << idx;
    const PQP_REAL* tan_sq_xpi_y_base = tan_sq_xpi_y+bases[idx];
    const PQP_REAL* tan_xpi_y_base = tan_xpi_y+bases[idx];
    int i = (2 << idx) - 1;
    int j = i;
    int stride = 1 << idx;

    while( stride != 0 ) {
        if( z_sq < tan_sq_xpi_y_base[i] * xy_sq ) {
            i -= stride;
        } else {
            i += stride;
        }
        if( fabs(dir[1]) < tan_xpi_y_base[j] * fabs(dir[0]) ) {
            j -= stride;
        } else {
            j += stride;
        }
        stride >>= 1;
    }

    if( z_sq > tan_sq_xpi_y_base[i] * xy_sq ) {
        i++;
    }

    if( i == tan_len - 1 && z_sq > tan_sq_xpi_y_base[i] * xy_sq ) {
        return lut[(dir[2] > 0.0 ? north[idx] : south[idx])];
    } else {
        i = (dir[2] > 0.0 ? (tan_len - i - 1) : (tan_len + i - 1));
        i <<= tan_len_shift;
    }

    if( fabs(dir[1]) > tan_xpi_y_base[j] * fabs(dir[0]) ) {
        j++;
    }

    if( j == tan_len - 1 && fabs(dir[1]) > tan_xpi_y_base[j] * fabs(dir[0])
    ) {
        return lut[(dir[1] > 0.0 ? i+tan_len : i+(3 << (tan_len_shift-2)) )];
    } else if( j == 0 ) {
        return lut[(dir[0] > 0.0 ? i : i+(tan_len<<1))];
    } else {
        return lut[(dir[1] > 0.0 ?
                   (dir[0] > 0.0 ? i+j : i+(tan_len<<1)-j) :
                   (dir[0] > 0.0 ? i+(tan_len<<2)-j : i+(tan_len<<1)+j) )];
    }
}

