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

//--------------------------------------------------------------------------
// File:   TriDist.cpp
// Author: Eric Larsen
// Description:
// contains SegPoints() for finding closest points on a pair of line
// segments and TriDist() for finding closest points on a pair of triangles
//--------------------------------------------------------------------------

#include "MatVec.h"
#ifdef _WIN32
#include <float.h>
#define isnan _isnan
#endif

#define Min_Value 1e-30

//--------------------------------------------------------------------------
// SegPoints() 
//
// Returns closest points between an segment pair.
// Implemented from an algorithm described in
//
// Vladimir J. Lumelsky,
// On fast computation of distance between line segments.
// In Information Processing Letters, no. 21, pages 55-61, 1985.   
//--------------------------------------------------------------------------

void
SegPoints(PQP_REAL VEC[3], 
	  PQP_REAL X[3], PQP_REAL Y[3],             // closest points
          const PQP_REAL P[3], const PQP_REAL A[3], // seg 1 origin, vector
          const PQP_REAL Q[3], const PQP_REAL B[3]) // seg 2 origin, vector
{
  PQP_REAL T[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;
  PQP_REAL TMP[3];

  VmV(T,Q,P);
  A_dot_A = VdotV(A,A);
  B_dot_B = VdotV(B,B);
  A_dot_B = VdotV(A,B);
  A_dot_T = VdotV(A,T);
  B_dot_T = VdotV(B,T);

  // t parameterizes ray P,A 
  // u parameterizes ray Q,B 

  PQP_REAL t,u;

  // compute t for the closest point on ray P,A to
  // ray Q,B

  PQP_REAL denom = A_dot_A*B_dot_B - A_dot_B*A_dot_B;

  t = (A_dot_T*B_dot_B - B_dot_T*A_dot_B) / denom;

  // clamp result so t is on the segment P,A

  if ((t < 0) || isnan(t)) t = 0; else if (t > 1) t = 1;

  // find u for point on ray Q,B closest to point at t

  u = (t*A_dot_B - B_dot_T) / B_dot_B;

  // if u is on segment Q,B, t and u correspond to 
  // closest points, otherwise, clamp u, recompute and
  // clamp t 

  if ((u <= 0) || isnan(u)) {

    VcV(Y, Q);

    t = A_dot_T / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Q, P);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Q, X);
    }
    else {
      VpVxS(X, P, A, t);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  }
  else if (u >= 1) {

    VpV(Y, Q, B);

    t = (A_dot_B + A_dot_T) / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Y, P);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Y, X);
    }
    else {
      VpVxS(X, P, A, t);
      VmV(T, Y, P);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  }
  else {

    VpVxS(Y, Q, B, u);

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(T, Q, X);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    }
    else {
      VpVxS(X, P, A, t);
      VcrossV(VEC, A, B);
      if (VdotV(VEC, T) < 0) {
        VxS(VEC, VEC, -1);
      }
    }
  }
}

//--------------------------------------------------------------------------
// TriDist() 
//
// Computes the closest points on two triangles, and returns the 
// distance between them.
// 
// S and T are the triangles, stored tri[point][dimension].
//
// If the triangles are disjoint, P and Q give the closest points of 
// S and T respectively. However, if the triangles overlap, P and Q 
// are basically a random pair of points from the triangles, not 
// coincident points on the intersection of the triangles, as might 
// be expected.
//--------------------------------------------------------------------------

PQP_REAL 
TriDist(PQP_REAL P[3], PQP_REAL Q[3],
        const PQP_REAL S[3][3], const PQP_REAL T[3][3])  
{
  // Compute vectors along the 6 sides

  PQP_REAL Sv[3][3], Tv[3][3];
  PQP_REAL VEC[3];

  VmV(Sv[0],S[1],S[0]);
  VmV(Sv[1],S[2],S[1]);
  VmV(Sv[2],S[0],S[2]);

  VmV(Tv[0],T[1],T[0]);
  VmV(Tv[1],T[2],T[1]);
  VmV(Tv[2],T[0],T[2]);

  // For each edge pair, the vector connecting the closest points 
  // of the edges defines a slab (parallel planes at head and tail
  // enclose the slab). If we can show that the off-edge vertex of 
  // each triangle is outside of the slab, then the closest points
  // of the edges are the closest points for the triangles.
  // Even if these tests fail, it may be helpful to know the closest
  // points found, and whether the triangles were shown disjoint

  PQP_REAL V[3];
  PQP_REAL Z[3];
  PQP_REAL minP[3], minQ[3], mindd;
  int shown_disjoint = 0;

  mindd = VdistV2(S[0],T[0]) + 1;  // Set first minimum safely high

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      // Find closest points on edges i & j, plus the 
      // vector (and distance squared) between these points

      SegPoints(VEC,P,Q,S[i],Sv[i],T[j],Tv[j]);
      
      VmV(V,Q,P);
      PQP_REAL dd = VdotV(V,V);

      // Verify this closest point pair only if the distance 
      // squared is less than the minimum found thus far.

      if (dd <= mindd)
      {
        VcV(minP,P);
        VcV(minQ,Q);
        mindd = dd;

        VmV(Z,S[(i+2)%3],P);
        PQP_REAL a = VdotV(Z,VEC);
        VmV(Z,T[(j+2)%3],Q);
        PQP_REAL b = VdotV(Z,VEC);

        if ((a <= 0) && (b >= 0)) return sqrt(dd);

        PQP_REAL p = VdotV(V, VEC);

        if (a < 0) a = 0;
        if (b > 0) b = 0;
        if ((p - a + b) > 0) shown_disjoint = 1;	
      }
    }
  }

  // No edge pairs contained the closest points.  
  // either:
  // 1. one of the closest points is a vertex, and the
  //    other point is interior to a face.
  // 2. the triangles are overlapping.
  // 3. an edge of one triangle is parallel to the other's face. If
  //    cases 1 and 2 are not true, then the closest points from the 9
  //    edge pairs checks above can be taken as closest points for the
  //    triangles.
  // 4. possibly, the triangles were degenerate.  When the 
  //    triangle points are nearly colinear or coincident, one 
  //    of above tests might fail even though the edges tested
  //    contain the closest points.

  // First check for case 1

  PQP_REAL Sn[3], Snl;       
  VcrossV(Sn,Sv[0],Sv[1]); // Compute normal to S triangle
  Snl = VdotV(Sn,Sn);      // Compute square of length of normal
  
  // If cross product is long enough,

  if (Snl > 1e-15)  
  {
    // Get projection lengths of T points

    PQP_REAL Tp[3]; 

    VmV(V,S[0],T[0]);
    Tp[0] = VdotV(V,Sn);

    VmV(V,S[0],T[1]);
    Tp[1] = VdotV(V,Sn);

    VmV(V,S[0],T[2]);
    Tp[2] = VdotV(V,Sn);

    // If Sn is a separating direction,
    // find point with smallest projection

    int point = -1;
    if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0))
    {
      if (Tp[0] < Tp[1]) point = 0; else point = 1;
      if (Tp[2] < Tp[point]) point = 2;
    }
    else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0))
    {
      if (Tp[0] > Tp[1]) point = 0; else point = 1;
      if (Tp[2] > Tp[point]) point = 2;
    }

    // If Sn is a separating direction, 

    if (point >= 0) 
    {
      shown_disjoint = 1;

      // Test whether the point found, when projected onto the 
      // other triangle, lies within the face.
    
      VmV(V,T[point],S[0]);
      VcrossV(Z,Sn,Sv[0]);
      if (VdotV(V,Z) > 0)
      {
        VmV(V,T[point],S[1]);
        VcrossV(Z,Sn,Sv[1]);
        if (VdotV(V,Z) > 0)
        {
          VmV(V,T[point],S[2]);
          VcrossV(Z,Sn,Sv[2]);
          if (VdotV(V,Z) > 0)
          {
            // T[point] passed the test - it's a closest point for 
            // the T triangle; the other point is on the face of S

            VpVxS(P,T[point],Sn,Tp[point]/Snl);
            VcV(Q,T[point]);
            return sqrt(VdistV2(P,Q));
          }
        }
      }
    }
  }

  PQP_REAL Tn[3], Tnl;       
  VcrossV(Tn,Tv[0],Tv[1]); 
  Tnl = VdotV(Tn,Tn);      
  
  if (Tnl > 1e-15)  
  {
    PQP_REAL Sp[3]; 

    VmV(V,T[0],S[0]);
    Sp[0] = VdotV(V,Tn);

    VmV(V,T[0],S[1]);
    Sp[1] = VdotV(V,Tn);

    VmV(V,T[0],S[2]);
    Sp[2] = VdotV(V,Tn);

    int point = -1;
    if ((Sp[0] > 0) && (Sp[1] > 0) && (Sp[2] > 0))
    {
      if (Sp[0] < Sp[1]) point = 0; else point = 1;
      if (Sp[2] < Sp[point]) point = 2;
    }
    else if ((Sp[0] < 0) && (Sp[1] < 0) && (Sp[2] < 0))
    {
      if (Sp[0] > Sp[1]) point = 0; else point = 1;
      if (Sp[2] > Sp[point]) point = 2;
    }

    if (point >= 0) 
    { 
      shown_disjoint = 1;

      VmV(V,S[point],T[0]);
      VcrossV(Z,Tn,Tv[0]);
      if (VdotV(V,Z) > 0)
      {
        VmV(V,S[point],T[1]);
        VcrossV(Z,Tn,Tv[1]);
        if (VdotV(V,Z) > 0)
        {
          VmV(V,S[point],T[2]);
          VcrossV(Z,Tn,Tv[2]);
          if (VdotV(V,Z) > 0)
          {
            VcV(P,S[point]);
            VpVxS(Q,S[point],Tn,Sp[point]/Tnl);
            return sqrt(VdistV2(P,Q));
          }
        }
      }
    }
  }

  // Case 1 can't be shown.
  // If one of these tests showed the triangles disjoint,
  // we assume case 3 or 4, otherwise we conclude case 2, 
  // that the triangles overlap.
  
  if (shown_disjoint)
  {
    VcV(P,minP);
    VcV(Q,minQ);
    return sqrt(mindd);
  }
  else return 0;
}


#include"PQP_Internal.h"

PQP_REAL 
TriDist(PQP_REAL P[3], PQP_REAL Q[3],
        const PQP_REAL S[3][3], const PQP_REAL T[3][3], ContactFeature &f1, ContactFeature &f2, bool &bCollided)  
{
  // Compute vectors along the 6 sides

  PQP_REAL Sv[3][3], Tv[3][3];
  PQP_REAL VEC[3];

  VmV(Sv[0],S[1],S[0]);
  VmV(Sv[1],S[2],S[1]);
  VmV(Sv[2],S[0],S[2]);

  VmV(Tv[0],T[1],T[0]);
  VmV(Tv[1],T[2],T[1]);
  VmV(Tv[2],T[0],T[2]);

  // For each edge pair, the vector connecting the closest points 
  // of the edges defines a slab (parallel planes at head and tail
  // enclose the slab). If we can show that the off-edge vertex of 
  // each triangle is outside of the slab, then the closest points
  // of the edges are the closest points for the triangles.
  // Even if these tests fail, it may be helpful to know the closest
  // points found, and whether the triangles were shown disjoint

  PQP_REAL V[3];
  PQP_REAL Z[3];
  PQP_REAL minP[3], minQ[3], mindd;
  int shown_disjoint = 0;

  mindd = VdistV2(S[0],T[0]) + 1;  // Set first minimum safely high

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      // Find closest points on edges i & j, plus the 
      // vector (and distance squared) between these points

      SegPoints(VEC,P,Q,S[i],Sv[i],T[j],Tv[j]);
      
      VmV(V,Q,P);
      PQP_REAL dd = VdotV(V,V);

      // Verify this closest point pair only if the distance 
      // squared is less than the minimum found thus far.

      if (dd <= mindd)
      {
        VcV(minP,P);
        VcV(minQ,Q);
        mindd = dd;

		f1.type = 1;
		f1.fid[0] = i;
		f2.type = 1;
		f2.fid[0] = j;

        VmV(Z,S[(i+2)%3],P);
        PQP_REAL a = VdotV(Z,VEC);
        VmV(Z,T[(j+2)%3],Q);
        PQP_REAL b = VdotV(Z,VEC);

        if ((a <= 0) && (b >= 0)) return sqrt(dd);

        PQP_REAL p = VdotV(V, VEC);

        if (a < 0) a = 0;
        if (b > 0) b = 0;
        if ((p - a + b) > 0) shown_disjoint = 1;	
      }
    }
  }

  // No edge pairs contained the closest points.  
  // either:
  // 1. one of the closest points is a vertex, and the
  //    other point is interior to a face.
  // 2. the triangles are overlapping.
  // 3. an edge of one triangle is parallel to the other's face. If
  //    cases 1 and 2 are not true, then the closest points from the 9
  //    edge pairs checks above can be taken as closest points for the
  //    triangles.
  // 4. possibly, the triangles were degenerate.  When the 
  //    triangle points are nearly colinear or coincident, one 
  //    of above tests might fail even though the edges tested
  //    contain the closest points.

  // First check for case 1

  PQP_REAL Sn[3], Snl;       
  VcrossV(Sn,Sv[0],Sv[1]); // Compute normal to S triangle
  Snl = VdotV(Sn,Sn);      // Compute square of length of normal
  
  // If cross product is long enough,

  if (Snl > 1e-15)  
  {
    // Get projection lengths of T points

    PQP_REAL Tp[3]; 

    VmV(V,S[0],T[0]);
    Tp[0] = VdotV(V,Sn);

    VmV(V,S[0],T[1]);
    Tp[1] = VdotV(V,Sn);

    VmV(V,S[0],T[2]);
    Tp[2] = VdotV(V,Sn);

    // If Sn is a separating direction,
    // find point with smallest projection

    int point = -1;
    if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0))
    {
      if (Tp[0] < Tp[1]) point = 0; else point = 1;
      if (Tp[2] < Tp[point]) point = 2;
    }
    else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0))
    {
      if (Tp[0] > Tp[1]) point = 0; else point = 1;
      if (Tp[2] > Tp[point]) point = 2;
    }

    // If Sn is a separating direction, 

    if (point >= 0) 
    {
      shown_disjoint = 1;

      // Test whether the point found, when projected onto the 
      // other triangle, lies within the face.
    
      VmV(V,T[point],S[0]);
      VcrossV(Z,Sn,Sv[0]);
      if (VdotV(V,Z) > 0)
      {
        VmV(V,T[point],S[1]);
        VcrossV(Z,Sn,Sv[1]);
        if (VdotV(V,Z) > 0)
        {
          VmV(V,T[point],S[2]);
          VcrossV(Z,Sn,Sv[2]);
          if (VdotV(V,Z) > 0)
          {
            // T[point] passed the test - it's a closest point for 
            // the T triangle; the other point is on the face of S

			  f1.type = 2;
			  f1.fid[0] = -1;
			  f2.type = 0;
			  f2.fid[0] = point;
			  
            VpVxS(P,T[point],Sn,Tp[point]/Snl);
            VcV(Q,T[point]);
            return sqrt(VdistV2(P,Q));
          }
        }
      }
    }
  }

  PQP_REAL Tn[3], Tnl;       
  VcrossV(Tn,Tv[0],Tv[1]); 
  Tnl = VdotV(Tn,Tn);      
  
  if (Tnl > 1e-15)  
  {
    PQP_REAL Sp[3]; 

    VmV(V,T[0],S[0]);
    Sp[0] = VdotV(V,Tn);

    VmV(V,T[0],S[1]);
    Sp[1] = VdotV(V,Tn);

    VmV(V,T[0],S[2]);
    Sp[2] = VdotV(V,Tn);

    int point = -1;
    if ((Sp[0] > 0) && (Sp[1] > 0) && (Sp[2] > 0))
    {
      if (Sp[0] < Sp[1]) point = 0; else point = 1;
      if (Sp[2] < Sp[point]) point = 2;
    }
    else if ((Sp[0] < 0) && (Sp[1] < 0) && (Sp[2] < 0))
    {
      if (Sp[0] > Sp[1]) point = 0; else point = 1;
      if (Sp[2] > Sp[point]) point = 2;
    }

    if (point >= 0) 
    { 
      shown_disjoint = 1;

      VmV(V,S[point],T[0]);
      VcrossV(Z,Tn,Tv[0]);
      if (VdotV(V,Z) > 0)
      {
        VmV(V,S[point],T[1]);
        VcrossV(Z,Tn,Tv[1]);
        if (VdotV(V,Z) > 0)
        {
          VmV(V,S[point],T[2]);
          VcrossV(Z,Tn,Tv[2]);
          if (VdotV(V,Z) > 0)
          {
			  f1.type = 0;
			  f1.fid[0] = point;
			  f2.type = 2;
			  f2.fid[0] = -1;
			  
            VcV(P,S[point]);
            VpVxS(Q,S[point],Tn,Sp[point]/Tnl);
            return sqrt(VdistV2(P,Q));
          }
        }
      }
    }
  }

  // Case 1 can't be shown.
  // If one of these tests showed the triangles disjoint,
  // we assume case 3 or 4, otherwise we conclude case 2, 
  // that the triangles overlap.
  
  if (shown_disjoint)
  {
    VcV(P,minP);
    VcV(Q,minQ);
    return sqrt(mindd);
  }
  else
  {
	  bCollided = true;
	   return 0;
  }
  
}
PQP_REAL Dis_Point_Plane(PQP_REAL point[3], PQP_REAL Norm[3], PQP_REAL point_plane[3])
{
	PQP_REAL temp[3];
    VmV(temp, point_plane, point);
	return fabs(VdotV(temp,Norm));
	
}

/*
* Ordinary inside-triangle test for p. The triangle normal is computed from the vertices.
*/
static inline bool _insideTriangle(PQP_REAL a[3], PQP_REAL  b[3], PQP_REAL c[3], PQP_REAL p[3])
{
	PQP_REAL n[3], da[3], db[3], dc[3], ba[3], ca[3],temp[3];
	double wa, wb, wc;

	VmV(ba,b,a);
	VmV(ca,c,a);
	VcrossV(n,ba,ca);

	VmV(da,a,p);
	VmV(db,b,p);
	VmV(dc,c,p);
    
	VcrossV(temp,db,dc);
	if ((wa = VdotV(temp,n))<1e-30) return false;

	VcrossV(temp,dc,da);
	if ((wb = VdotV(temp,n))<1e-30) return false;

	VcrossV(temp,da,db);
	if ((wc = VdotV(temp,n))<1e-30) return false;


	return true;

	
}


/****************************************************************************************************/
/* intersection with ray and plane
/* the direction is start from endpoint
/****************************************************************************************************/
bool RayPlaneIntersection(PQP_REAL Inter_point[3], PQP_REAL Endpoint[3], PQP_REAL Direction[3], PQP_REAL Norm[3], PQP_REAL D)
{
	PQP_REAL temp;
	temp=VdotV(Direction, Norm);
	if (temp<=0.0f)
	{
		return false;//intersection does not exist
	}

	PQP_REAL t=-( VdotV(Endpoint,Norm) + D)/temp;
	VpVxS(Inter_point,Endpoint,Direction,t);
	return true;
  
}

//******************************************************************************//
//copy from http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm     //
//******************************************************************************//
inline PQP_REAL Disline_line(PQP_REAL L1[2][3],PQP_REAL L2[2][3],PQP_REAL p[3], PQP_REAL q[3])
{
	PQP_REAL u[3], v[3], w0[3];

	VmV(u,L1[1],L1[0]);

	VmV(v,L2[1],L2[0]);

	VmV(w0,L1[0],L2[0]);

	PQP_REAL a,b,c,d,e,D,sc,tc,sN,tN;
	PQP_REAL sD,tD;

	a=VdotV(u,u); b=VdotV(v,u); c=VdotV(v,v); d=VdotV(w0,u);e=VdotV(v,w0);

	D=a*c-b*b;
	sD=tD=D;
	if (D< 1e-30 )//0.0
	{
		sN=0.0;
		sD=1.0;
		tN=e;
		tD=c;
	}
	else 
	{                // get the closest points on the infinite lines
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0.0)
			sN = 0.0;
		else if (-d > a)
			sN = sD;
		else {
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
		tN = tD;
		// recompute sc for this edge
		if ((-d + b) < 0.0)
			sN = 0;
		else if ((-d + b) > a)
			sN = sD;
		else {
			sN = (-d + b);
			sD = a;
		}
	}
	// finally do the division to get sc and tc
	sc = (abs(sN) < 1e-30 ? 0.0 : sN / sD);
	tc = (abs(tN) < 1e-30 ? 0.0 : tN / tD);//

	// get the difference of the two closest points
	PQP_REAL dp[3];
	VpVxS(p,L1[0],u,sc);
	VpVxS(q,L2[0],v,tc);
	VmV(dp,p,q);
	return Vlength(dp);

}
bool InterSegments(PQP_REAL a[3],PQP_REAL d_a[3],PQP_REAL b[3],PQP_REAL d_b[3],PQP_REAL p[3],PQP_REAL q[3])
{
	PQP_REAL temp[3];
	VcrossV(temp,d_a,d_b);
	if (Vlength(temp)<1e-30)
	{
		return false;
	}
	
	PQP_REAL L1[2][3];PQP_REAL L2[2][3];
	VcV(L1[0],a);
	VpV(L1[1],a,d_a);

	VcV(L2[0],b);
	VpV(L2[1],b,d_b);

	if(Disline_line(L1,L2,p,q)<1e-30) return true;
	return false;



}
void IntersectionTriangles_CoPlane(PQP_REAL tri1[3][3], PQP_REAL tri2[3][3],PQP_REAL InterPoints[6][3],int* number)
{
	PQP_REAL p[3], q[3], a[3], b[3], d_a[3], d_b[3];
	bool b_inter;
    number[0]=0;

	VcV(a,tri2[0]);
	VmV(d_a,tri2[1],tri2[0]);

	VcV(b,tri1[0]);
	VmV(d_b,tri1[1],tri1[0]);

	
	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}

	VcV(b,tri1[1]);
	VmV(d_b,tri1[2],tri1[1]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}

	VcV(b,tri1[2]);
	VmV(d_b,tri1[0],tri1[2]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}
	///////////////////

	VcV(a,tri2[1]);
	VmV(d_a,tri2[2],tri2[1]);

	VcV(b,tri1[0]);
	VmV(d_b,tri1[1],tri1[0]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}

	VcV(b,tri1[1]);
	VmV(d_b,tri1[2],tri1[1]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}

	VcV(b,tri1[2]);
	VmV(d_b,tri1[0],tri1[2]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}

	////////////////////////////////////

	VcV(a,tri2[2]);
	VmV(d_a,tri2[0],tri2[2]);

	VcV(b,tri1[0]);
	VmV(d_b,tri1[1],tri1[0]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}

	VcV(b,tri1[1]);
	VmV(d_b,tri1[2],tri1[1]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}

	VcV(b,tri1[2]);
	VmV(d_b,tri1[0],tri1[2]);


	b_inter=InterSegments(a,d_a,b,d_b,p,q);
	if (b_inter)
	{
		VcV(InterPoints[ number[0] ],q);
		number[0]++;
	}



}
PQP_REAL 
DisPointTriangle_Direction(PQP_REAL Point[3],PQP_REAL tri[3][3],PQP_REAL dirc[3],PQP_REAL norm[3], PQP_REAL D)
{
	PQP_REAL Interpoint[3], temp[3];

	bool b_inter=RayPlaneIntersection(Interpoint, Point, dirc,norm ,D);
	if (!b_inter)
	{
		return 1e+30;
	}
	
	bool b_inside = _insideTriangle(tri[0],tri[1], tri[2], Interpoint);

	if (b_inside)
	{
		VmV(temp,Point,Interpoint);
		return Vlength(temp);
	}
	return 1e+30;


   
}

PQP_REAL 
TriDist_Direction(PQP_REAL p[3], PQP_REAL q[3], 
				  PQP_REAL s[3][3],  PQP_REAL t[3][3],  PQP_REAL Direction[3], bool *bExisted)
{
	PQP_REAL e0_s[3],e2_s[3], norm_s[3];
	VmV(e0_s,s[1],s[0]);
	VmV(e2_s,s[2],s[0]);
	VcrossV(norm_s,e0_s,e2_s);
	Vnormalize(norm_s);
   
	PQP_REAL newt[3][3], Dir_t_s[3], D;

	VxS(Dir_t_s,Direction,1);
	D=-VdotV(norm_s,s[0]);
	bool bexit;
	for(int i=0; i<3; i++)
	{		
		bexit=RayPlaneIntersection(newt[i],t[i],Dir_t_s,norm_s,D);//projection triangle newt
		if (!bexit)
		{
			bExisted[0]=bexit;
			return 1e+30;
		}
			
	}
	PQP_REAL InterPoints[9][3];// Intersections between triangles
	int number;
	IntersectionTriangles_CoPlane(s, newt,InterPoints,&number);//intersection points between B and A'
    
	/*if (number==0)
	{
		bExisted[0]=bexit;
		return 1e+30;
	}*/

	PQP_REAL *Dist=new PQP_REAL[6+number];

	PQP_REAL norm_t[3],edge1[3],edge2[3];

	VmV(edge1, t[1], t[0]);
	VmV(edge2, t[2], t[0]);

	VcrossV(norm_t,edge1,edge2);
	Vnormalize(norm_t);
    VxS(Direction,Direction,-1);
	PQP_REAL D_t=-VdotV(norm_t,t[0]);

//////////////////////////////////////////////////////////////////////////////////////////////////////
/*	for (int i=0;i<3;i++)//A'
	{
		Dist[i] = DisPointTriangle_Direction(s[i],t,Direction,norm_t,D_t);//1e+30;

	}
	for (int i=0;i<3;i++)//B triangle
	{
		if (_insideTriangle(s[0],s[1],s[2],newt[i]))
		{
			Dist[i+3] =DisPointTriangle_Direction(newt[i],t,Direction,norm_t,D_t);// 1e+30;//
		}
		else
		{
			Dist[i+3] = 1e+30;
		   
		}
		

	}*/
//////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i=0;i<3;i++)//B triangle
	{

		if (_insideTriangle(newt[0],newt[1],newt[2],s[i]))
		{
			Dist[i] =DisPointTriangle_Direction(s[i],t,Direction,norm_t,D_t);// 1e+30;//
		}
		else
		{
			Dist[i] = 1e+30;

		}

		

	}
	for (int i=0;i<3;i++)//A'
	{
		if (_insideTriangle(s[0],s[1],s[2],newt[i]))
		{
			PQP_REAL temp[3];
			VmV(temp,newt[i],t[i]);
			Dist[i+3] = Vlength(temp);
		}
		else
		{
			Dist[i+3] = 1e+30;

		}


	}

	for (int i=0;i<number;i++)
	{
		Dist[i+6] = DisPointTriangle_Direction(InterPoints[i],t,Direction,norm_t,D_t);//1e+30;

	}
    PQP_REAL min_dist=1e+30;

	bExisted[0]=false;
	for (int i=0;i<6+number;i++)
	{
		if (min_dist>Dist[i])
		{
			min_dist=Dist[i];
			
		}
	}
	if (min_dist<1e+20)
	{
		bExisted[0]=true;
	}
	else
	{
		bExisted[0]=false;
	
	}
	delete Dist;
	return min_dist;

}