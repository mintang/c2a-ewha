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

#ifndef PQP_TRIDIST_H
#define PQP_TRIDIST_H

#include "PQP_Compile.h"

// TriDist()
//
// computes the closest points on two triangles, and returns the 
// distance between them.
// 
// s and t are the triangles, stored tri[point][dimension].
//
// If the triangles are disjoint, p and q give the closest points of 
// s and t respectively. However, if the triangles overlap, p and q 
// are basically a random pair of points from the triangles, not 
// coincident points on the intersection of the triangles, as might 
// be expected.

PQP_REAL 
TriDist(PQP_REAL p[3], PQP_REAL q[3], 
        const PQP_REAL s[3][3], const PQP_REAL t[3][3]);
//TriDist with contact feature returned
PQP_REAL 
TriDist(PQP_REAL p[3], PQP_REAL q[3], 
        const PQP_REAL s[3][3], const PQP_REAL t[3][3], ContactFeature &f1, ContactFeature &f2, bool &bCollided);

PQP_REAL 
TriDist_Direction(PQP_REAL p[3], PQP_REAL q[3], 
        PQP_REAL s[3][3], PQP_REAL t[3][3], PQP_REAL Direct[3], bool *bExisted);

#endif