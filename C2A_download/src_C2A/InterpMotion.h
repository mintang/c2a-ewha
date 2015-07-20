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

  US Mail:             S. Gottschalk, E. Larsen
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

#ifndef INTERP_MOTION_H
#define INTERP_MOTION_H

#include "SWIFT_linalg.h"
#include "PQP_compile.h"
#include "bv.h"
#include "GMPcommon.h"
#include "PQP.h"


class CInterpMotion
{
public:
  CInterpMotion(GMP_INTERP_MODE itpMode, 
							 const PQP_REAL R0[3][3], const PQP_REAL T0[3],
							 const PQP_REAL R1[3][3], const PQP_REAL T1[3]);

  CInterpMotion(GMP_INTERP_MODE itpMode, 
							 const PQP_REAL q0[7], const PQP_REAL q1[7]);
  CInterpMotion();
  void SetSrcDst(const PQP_REAL q0[7], const PQP_REAL q1[7]);

  virtual ~CInterpMotion();
  virtual double computeTOC(PQP_REAL d, PQP_REAL r1);
  virtual void velocity(void);
  virtual void velocity_relate(PQP_REAL T_relate[3]);
  virtual string Print() {return "Motion Interpolation"; }

  // Interpolation: the results are represented as a quaternion
  // (q[0],... q[4]) Rotation, (q[5], q[6], q[7]) 
  virtual bool integrate(const double dt, PQP_REAL qua[7]) = 0;
  virtual bool integrate(const double dt, PQP_REAL R[3][3], PQP_REAL T[3]);

  
  virtual double computeTOC(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3]);//compute the TOC using directional motion bound, //use the BV's angular radius to compute TOC,
  virtual double compute_MotionBound(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3]);
  virtual double computeTOC_relate(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3]);
  virtual double computeTOC_MotionBound(PQP_REAL T[3],PQP_REAL R[3][3],  PQP_REAL d,BV *V,
						PQP_REAL S[3]);//compute the TOC using directional motion bound, //use the BV's angular radius to compute TOC,
  virtual double computeTOC_MotionBoundTri(PQP_REAL a,PQP_REAL T[3],PQP_REAL R[3][3], PQP_Model *o1, PQP_REAL d, BV *V, PQP_REAL p[3], PQP_REAL q[3]);
  virtual bool CAonRSS(PQP_REAL r1[3][3], PQP_REAL tt1[3],	  PQP_REAL r2[3][3],PQP_REAL tt2[3],  PQP_REAL R[3][3],
	  PQP_REAL T[3], BV *b1, BV *b2, 
	  PQP_REAL *mint,  PQP_REAL *distance);
  virtual bool  CAonNonAdjacentTriangles(PQP_REAL R[3][3], PQP_REAL T[3],   PQP_REAL triA[3][3],PQP_REAL triB[3][3],  PQP_REAL *mint, PQP_REAL *distance,int* Inter);

  SWIFT_Quaternion DeltaRt(PQP_REAL t);    // R_t
  SWIFT_Quaternion AbsoluteRt(PQP_REAL t); // R(t) = R_t * R(0)

  void LinearAngularVel(SWIFT_Triple &axis, PQP_REAL &angVel);
  

public:
  PQP_REAL ConservD(PQP_REAL d); // ConservD = d - security_distance_ratio * m_toc_delta 
  PQP_REAL m_toc_delta; // the d_delta for the CCD query



  GMP_INTERP_MODE m_itpMode;
  SWIFT_Transformation transform; // huangxin 12/01/2007
  SWIFT_Transformation transform_s; //huangxin 12/01/2007
  SWIFT_Transformation transform_t; //huangxin 12/01/2007
  SWIFT_Triple m_rel_tra; // The amount translation between q0 and q1 

  SWIFT_Transformation m_transform_t_inv; // the inverse of transform_t: avoid the comptuation of the inverse during each integration
  int m_time;
  int m_time1;


//protected:
  SWIFT_Triple cv; // linear translation  velocity
  SWIFT_Triple cv_relate;
  SWIFT_Triple m_axis; // the axis of the rotation, normalized
  PQP_REAL m_angVel; // the angular velocity around that axis: counter clockwise
					   // it might be negative, which means the clockwise rotation
  

};


class CInterpMotion_Linear : public CInterpMotion
{
public:
  virtual ~CInterpMotion_Linear();
  virtual void velocity(void);
  virtual void velocity_relate(PQP_REAL T_relate[3]);
  virtual bool integrate(const double dt, PQP_REAL qua[7]);

  CInterpMotion_Linear(const PQP_REAL R0[3][3], const PQP_REAL T0[3],
					   const PQP_REAL R1[3][3], const PQP_REAL T1[3]);

  CInterpMotion_Linear(const PQP_REAL q0[7], const PQP_REAL q1[7]);
  CInterpMotion_Linear();
  virtual string Print() {return "Linear Interpolation"; }

  virtual double computeTOC(PQP_REAL d, PQP_REAL r1);

 
  virtual double computeTOC(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3]);//compute the TOC using directional motion bound, //use the BV's angular radius to compute TOC, 
  virtual double compute_MotionBound(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3]);
  virtual double computeTOC_relate(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3]);
  virtual double computeTOC_MotionBound(PQP_REAL T[3],PQP_REAL R[3][3], PQP_REAL d, 
	                    BV *V, PQP_REAL S[3]);//using motion bound of RSS 
  virtual double computeTOC_MotionBoundTri( PQP_REAL a, PQP_REAL T[3],PQP_REAL R[3][3],PQP_Model *o1, PQP_REAL d, BV *V, PQP_REAL p[3], PQP_REAL q[3]);
  virtual bool CAonRSS(PQP_REAL r1[3][3], PQP_REAL tt1[3],	  PQP_REAL r2[3][3],PQP_REAL tt2[3], 	  PQP_REAL R[3][3],PQP_REAL T[3], 	  BV *b1, BV *b2, 
	          PQP_REAL *mint,  PQP_REAL *distance);
  virtual bool  CAonNonAdjacentTriangles(PQP_REAL R[3][3], PQP_REAL T[3],   PQP_REAL triA[3][3],PQP_REAL triB[3][3],  PQP_REAL *mint, PQP_REAL *distance , int* Inter);

};


#endif;