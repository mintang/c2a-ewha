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

#ifndef PQP_BV_H
#define PQP_BV_H

#include <math.h>
#include "PQP_Tri.h"
#include "PQP_Compile.h"

class BV
{
public:

  PQP_REAL R[3][3];     // orientation of RSS & OBB
  PQP_REAL R_abs[3][3]; // abs orientation of RSS & OBB: Liangjun
  PQP_REAL R_loc[3][3];
  PQP_REAL trilength;

#if PQP_BV_TYPE & RSS_TYPE
  PQP_REAL Tr[3];       // position of rectangle
  PQP_REAL Tr_abs[3];   // abs position of rectangle: Liangjun
  PQP_REAL l[2];        // side lengths of rectangle
  PQP_REAL r;           // radius of sphere summed with rectangle to form RSS
  PQP_REAL Corner[3][3];// the four corners of rectangle 
  PQP_REAL Length_C1;// the length of C1
  PQP_REAL Length_C2;// the length of C2
  
  int deep;
 

  PQP_REAL max_Length;
#endif

#if PQP_BV_TYPE & OBB_TYPE
  PQP_REAL To[3];       // position of obb
  PQP_REAL To_abs[3];   // abs position of obb
  PQP_REAL d[3];        // (half) dimensions of obb
#endif

  int first_child;      // positive value is index of first_child bv
                        // negative value is -(index + 1) of triangle
  
  int parentID; //the id of parent BV
  PQP_REAL com[3]; // the center of Mass, 
  PQP_REAL comRoot[3]; // the center of model that BV belongs to, 
  PQP_REAL angularRadius; // the angular radius of the BV to comRoot, 
  PQP_REAL BV_Size;
  void ComputeCenterOfMassBV(PQP_Tri *tris, int num_tris);
  void ComputeAngularRadius(PQP_REAL O[3][3], PQP_Tri *tris, int num_tris);
  

  BV();
  ~BV();
  int      Leaf()    { return first_child < 0; }
  PQP_REAL GetSize(); 
  void     FitToTris(PQP_REAL O[3][3], PQP_Tri *tris, int num_tris, PQP_REAL comR[3]);
  void     FitToTris_Corner(PQP_REAL O[3][3], PQP_Tri *tris, int num_tris, PQP_REAL comR[3]);
};

inline
PQP_REAL 
BV::GetSize()
{
#if PQP_BV_TYPE & RSS_TYPE
  return (sqrt(l[0]*l[0] + l[1]*l[1]) + 2*r);
#else
  return (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
#endif
}

int
BV_Overlap(PQP_REAL R[3][3], PQP_REAL T[3], BV *b1, BV *b2);

#if PQP_BV_TYPE & RSS_TYPE
PQP_REAL
BV_Distance(PQP_REAL R[3][3], PQP_REAL T[3], BV *b1, BV *b2, PQP_REAL S[3]);

// The closest points between the two rectangles not the two RSS!!!
PQP_REAL
BV_Distance(PQP_REAL R[3][3], PQP_REAL T[3], BV *b1, BV *b2,
			PQP_REAL P[3], PQP_REAL Q[3],
			bool &bValidPQ, PQP_REAL S[3]);
#endif

#endif



