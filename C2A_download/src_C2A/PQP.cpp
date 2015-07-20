                                             
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
#include <iostream>
using namespace std;


#include <stdio.h>
#include <string.h>
#include "PQP.h"
#include "BVTQ.h"

#include "MatVec.h"
#include "GetTime.h"
#include "TriDist.h"
#include "InterpMotion.h"

#include "Build.h"



bool bPQPContact_Represenative = true;
double PQPContact_radius = 5;




enum BUILD_STATE
{ 
  PQP_BUILD_STATE_EMPTY,     // empty state, immediately after constructor
  PQP_BUILD_STATE_BEGUN,     // after BeginModel(), state for adding triangles
  PQP_BUILD_STATE_PROCESSED  // after tree has been built, ready to use
};

PQP_Model::PQP_Model()
{
  // no bounding volume tree yet

  b = 0;  
  num_bvs_alloced = 0;
  num_bvs = 0;

  // no tri list yet

  tris = 0;
  num_tris = 0;
  num_tris_alloced = 0;

  bMotionCoherence = false;
  motionTri = 0;
  motionBV = 0;

  level = 3;
//  numLocalTris = pow(2.0, level);
//  localTris = new int[numLocalTris];
  levelPID = 0;
  last_tri = 0;

  build_state = PQP_BUILD_STATE_EMPTY;
}

PQP_Model::~PQP_Model()
{
  if (b != NULL)
    delete [] b;
  if (tris != NULL)
    delete [] tris;
  if(trisConst != NULL)
	delete [] trisConst;

}

int
PQP_Model::BeginModel(int n)
{
  // reset to initial state if necessary

  if (build_state != PQP_BUILD_STATE_EMPTY) 
  {
    delete [] b;
    delete [] tris;
  
    num_tris = num_bvs = num_tris_alloced = num_bvs_alloced = 0;
  }

  // prepare model for addition of triangles

  if (n <= 0) n = 8;
  num_tris_alloced = n;
  tris = new PQP_Tri[n];
  if (!tris) 
  {
    fprintf(stderr, "PQP Error!  Out of memory for tri array on "
                    "BeginModel() call!\n");
    return PQP_ERR_MODEL_OUT_OF_MEMORY;  
  }

  // give a warning if called out of sequence

  if (build_state != PQP_BUILD_STATE_EMPTY)
  {
    fprintf(stderr,
            "PQP Warning! Called BeginModel() on a PQP_Model that \n"
            "was not empty. This model was cleared and previous\n"
            "triangle additions were lost.\n");
    build_state = PQP_BUILD_STATE_BEGUN;
    return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  build_state = PQP_BUILD_STATE_BEGUN;
  return PQP_OK;
}

void
PQP_Model::SetCenterOfMass(PQP_REAL c[3])
{
	VcV(com, c);
}
PQP_REAL
PQP_Model::Max_Cross_Product(PQP_REAL dir[3])
{
	PQP_REAL t[3];
	PQP_REAL cp[3];
	PQP_REAL v=0.0, v_;
	for (int i=0; i<num_tris; i++)
	{
		t[0] = tris[i].p1[0];		
		t[1] = tris[i].p1[1];		
		t[2] = tris[i].p1[2];
		VcrossV(cp, t, dir);
		v_=(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);
		
		if( v_ > v ) v=v_;
		
		t[0] = tris[i].p2[0];		
		t[1] = tris[i].p2[1];		
		t[2] = tris[i].p2[2];
		VcrossV(cp, t, dir);
		v_=(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);
		
		if( v_ > v ) v=v_;
		
		t[0] = tris[i].p3[0];		
		t[1] = tris[i].p3[1];		
		t[2] = tris[i].p3[2];
		VcrossV(cp, t, dir);
		v_=(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);
		
		if( v_ > v ) 
			v = v_;
	}
	
	return sqrt(v);
}

int
PQP_Model::AddTri(const PQP_REAL *p1, 
                  const PQP_REAL *p2, 
                  const PQP_REAL *p3, 
                  int id)
{
  if (build_state == PQP_BUILD_STATE_EMPTY)
  {
    BeginModel();
  }
  else if (build_state == PQP_BUILD_STATE_PROCESSED)
  {
    fprintf(stderr,"PQP Warning! Called AddTri() on PQP_Model \n"
                   "object that was already ended. AddTri() was\n"
                   "ignored.  Must do a BeginModel() to clear the\n"
                   "model for addition of new triangles\n");
    return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
  }
        
  // allocate for new triangles

  if (num_tris >= num_tris_alloced)
  {
    PQP_Tri *temp;
    temp = new PQP_Tri[num_tris_alloced*2];
    if (!temp)
    {
      fprintf(stderr, "PQP Error!  Out of memory for tri array on"
	              " AddTri() call!\n");
      return PQP_ERR_MODEL_OUT_OF_MEMORY;  
    }
    memcpy(temp, tris, sizeof(PQP_Tri)*num_tris);
    delete [] tris;
    tris = temp;
    num_tris_alloced = num_tris_alloced*2;
  }
  
  // initialize the new triangle

  tris[num_tris].p1[0] = p1[0];
  tris[num_tris].p1[1] = p1[1];
  tris[num_tris].p1[2] = p1[2];

  tris[num_tris].p2[0] = p2[0];
  tris[num_tris].p2[1] = p2[1];
  tris[num_tris].p2[2] = p2[2];

  tris[num_tris].p3[0] = p3[0];
  tris[num_tris].p3[1] = p3[1];
  tris[num_tris].p3[2] = p3[2];

  tris[num_tris].id = id;

  num_tris += 1;

  return PQP_OK;
}
int
PQP_Model::AddTri(const PQP_REAL *p1, 
                  const PQP_REAL *p2, 
                  const PQP_REAL *p3, 
                  int id, int i1,int i2,int i3)
{
	if (build_state == PQP_BUILD_STATE_EMPTY)
	{
		BeginModel();
	}
	else if (build_state == PQP_BUILD_STATE_PROCESSED)
	{
		fprintf(stderr,"PQP Warning! Called AddTri() on PQP_Model \n"
			"object that was already ended. AddTri() was\n"
			"ignored.  Must do a BeginModel() to clear the\n"
			"model for addition of new triangles\n");
		return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
	}
	
	// allocate for new triangles
	
	if (num_tris >= num_tris_alloced)
	{
		PQP_Tri *temp;
		temp = new PQP_Tri[num_tris_alloced*2];
		if (!temp)
		{
			fprintf(stderr, "PQP Error!  Out of memory for tri array on"
				" AddTri() call!\n");
			return PQP_ERR_MODEL_OUT_OF_MEMORY;  
		}
		memcpy(temp, tris, sizeof(PQP_Tri)*num_tris);
		delete [] tris;
		tris = temp;
		num_tris_alloced = num_tris_alloced*2;
	}
	
	// initialize the new triangle
	
	tris[num_tris].p1[0] = p1[0];
	tris[num_tris].p1[1] = p1[1];
	tris[num_tris].p1[2] = p1[2];
	
	tris[num_tris].p2[0] = p2[0];
	tris[num_tris].p2[1] = p2[1];
	tris[num_tris].p2[2] = p2[2];
	
	tris[num_tris].p3[0] = p3[0];
	tris[num_tris].p3[1] = p3[1];
	tris[num_tris].p3[2] = p3[2];
	
	tris[num_tris].id = id;
	tris[num_tris].Index[0] = i1;
	tris[num_tris].Index[1] = i2;
	tris[num_tris].Index[2] = i3;
	
	num_tris += 1;
	
	return PQP_OK;
}
void ComputeAbsTra(PQP_REAL R[3][3], PQP_REAL T[3], PQP_Model *o, int bv_ind)
{

	BV *pBV = &o->b[bv_ind];

	// R1*(R2*p + T2) + T1 = R1*R2*p + R1*T2 + T1
    MxM(pBV->R_abs,R, pBV->R);
#if PQP_BV_TYPE & OBB_TYPE
    MxVpV(pBV->To_abs, R, pBV->To,T);
#else
    MxVpV(pBV->Tc_abs, R, pBV->Tr,T);
#endif
		
	
	int l = pBV->Leaf();

	// recursive to compute the BV of the child
	if(!l)
	{
		int c1 = pBV->first_child;
		int c2 = c1 + 1;
#if PQP_BV_TYPE & OBB_TYPE
		ComputeAbsTra(pBV->R_abs, pBV->To_abs, 
					  o, c1);

		ComputeAbsTra(pBV->R_abs, pBV->To_abs, 
					  o, c2);
#else
		ComputeAbsTra(pBV->R_abs, pBV->Tc_abs, 
					  o, c1);

		ComputeAbsTra(pBV->R_abs, pBV->Tc_abs, 
					  o, c2);
#endif
	}
}

int
PQP_Model::EndModel()
{
	trisConst = new PQP_Tri[num_tris];
	for(int i = 0; i < num_tris; i++)
	{
		trisConst[i] = tris[i];
	}
	
	if (build_state == PQP_BUILD_STATE_PROCESSED)
  {
    fprintf(stderr,"PQP Warning! Called EndModel() on PQP_Model \n"
                   "object that was already ended. EndModel() was\n"
                   "ignored.  Must do a BeginModel() to clear the\n"
                   "model for addition of new triangles\n");
    return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // report error is no tris

  if (num_tris == 0)
  {
    fprintf(stderr,"PQP Error! EndModel() called on model with"
                   " no triangles\n");
    return PQP_ERR_BUILD_EMPTY_MODEL;
  }

  // shrink fit tris array 

  if (num_tris_alloced > num_tris)
  {
    PQP_Tri *new_tris = new PQP_Tri[num_tris];
    if (!new_tris) 
    {
      fprintf(stderr, "PQP Error!  Out of memory for tri array "
                      "in EndModel() call!\n");
      return PQP_ERR_MODEL_OUT_OF_MEMORY;  
    }
    memcpy(new_tris, tris, sizeof(PQP_Tri)*num_tris);
    delete [] tris;
    tris = new_tris;
    num_tris_alloced = num_tris;
  }

  //compute the center of mass and the radius of the model
  ComputeCenterOfMass();
  ComputeRadius();
  
  // create an array of BVs for the model

  b = new BV[2*num_tris - 1];
  if (!b)
  {
    fprintf(stderr,"PQP Error! out of memory for BV array "
                   "in EndModel()\n");
    return PQP_ERR_MODEL_OUT_OF_MEMORY;
  }
  num_bvs_alloced = 2*num_tris - 1;
  num_bvs = 0;

  // we should build the model now.

  build_model(this);
  build_state = PQP_BUILD_STATE_PROCESSED;

  last_tri = tris;

  PQP_REAL R[3][3];
  PQP_REAL T[3];

  R[0][0] = 1; R[0][1] = 0; R[0][2] = 0;
  R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
  R[2][0] = 0; R[2][1] = 0; R[2][2] = 1;

  T[0] = T[1] = T[2] = 0;

  // starting from the root
  ComputeAbsTra(R, T, this, 0);


  return PQP_OK;
}


int
PQP_Model::MemUsage(int msg)
{
  int mem_bv_list = sizeof(BV)*num_bvs;
  int mem_tri_list = sizeof(PQP_Tri)*num_tris;

  int total_mem = mem_bv_list + mem_tri_list + sizeof(PQP_Model);

  if (msg) 
  {
    fprintf(stderr,"Total for model %x: %d bytes\n", this, total_mem);
    fprintf(stderr,"BVs: %d alloced, take %d bytes each\n", 
            num_bvs, sizeof(BV));
    fprintf(stderr,"Tris: %d alloced, take %d bytes each\n", 
            num_tris, sizeof(PQP_Tri));
  }
  
  return total_mem;
}

//  COLLIDE STUFF
//
//--------------------------------------------------------------------------

PQP_CollideResult::PQP_CollideResult()
{
  pairs = 0;
  num_pairs = num_pairs_alloced = 0;
  num_bv_tests = 0;
  num_tri_tests = 0;
}

PQP_CollideResult::~PQP_CollideResult()
{
  delete [] pairs;
}

void
PQP_CollideResult::FreePairsList()
{
  num_pairs = num_pairs_alloced = 0;
  delete [] pairs;
  pairs = 0;
}

// may increase OR reduce mem usage
void
PQP_CollideResult::SizeTo(int n)
{
  CollisionPair *temp;

  if (n < num_pairs) 
  {
    fprintf(stderr, "PQP Error: Internal error in "
                    "'PQP_CollideResult::SizeTo(int n)'\n");
    fprintf(stderr, "       n = %d, but num_pairs = %d\n", n, num_pairs);
    return;
  }
  
  temp = new CollisionPair[n];
  memcpy(temp, pairs, num_pairs*sizeof(CollisionPair));
  delete [] pairs;
  pairs = temp;
  num_pairs_alloced = n;
  return;
}

void
PQP_CollideResult::Add(int a, int b)
{
  if (num_pairs >= num_pairs_alloced) 
  {
    // allocate more

    SizeTo(num_pairs_alloced*2+8);
  }

  // now proceed as usual

  pairs[num_pairs].id1 = a;
  pairs[num_pairs].id2 = b;
  num_pairs++;
}

// TRIANGLE OVERLAP TEST
       
inline
PQP_REAL
max_PQP(PQP_REAL a, PQP_REAL b, PQP_REAL c)
{
  PQP_REAL t = a;
  if (b > t) t = b;
  if (c > t) t = c;
  return t;
}

inline
PQP_REAL
min_PQP(PQP_REAL a, PQP_REAL b, PQP_REAL c)
{
  PQP_REAL t = a;
  if (b < t) t = b;
  if (c < t) t = c;
  return t;
}

int
project6(PQP_REAL *ax, 
         PQP_REAL *p1, PQP_REAL *p2, PQP_REAL *p3, 
         PQP_REAL *q1, PQP_REAL *q2, PQP_REAL *q3)
{
  PQP_REAL P1 = VdotV(ax, p1);
  PQP_REAL P2 = VdotV(ax, p2);
  PQP_REAL P3 = VdotV(ax, p3);
  PQP_REAL Q1 = VdotV(ax, q1);
  PQP_REAL Q2 = VdotV(ax, q2);
  PQP_REAL Q3 = VdotV(ax, q3);
  
  PQP_REAL mx1 = max_PQP(P1, P2, P3);
  PQP_REAL mn1 = min_PQP(P1, P2, P3);
  PQP_REAL mx2 = max_PQP(Q1, Q2, Q3);
  PQP_REAL mn2 = min_PQP(Q1, Q2, Q3);

  if (mn1 > mx2) return 0;
  if (mn2 > mx1) return 0;
  return 1;
}

// very robust triangle intersection test
// uses no divisions
// works on coplanar triangles
int 
TriContact(PQP_REAL *P1, PQP_REAL *P2, PQP_REAL *P3,
           PQP_REAL *Q1, PQP_REAL *Q2, PQP_REAL *Q3) 
{

  // One triangle is (p1,p2,p3).  Other is (q1,q2,q3).
  // Edges are (e1,e2,e3) and (f1,f2,f3).
  // Normals are n1 and m1
  // Outwards are (g1,g2,g3) and (h1,h2,h3).
  //  
  // We assume that the triangle vertices are in the same coordinate system.
  //
  // First thing we do is establish a new c.s. so that p1 is at (0,0,0).

  PQP_REAL p1[3], p2[3], p3[3];
  PQP_REAL q1[3], q2[3], q3[3];
  PQP_REAL e1[3], e2[3], e3[3];
  PQP_REAL f1[3], f2[3], f3[3];
  PQP_REAL g1[3], g2[3], g3[3];
  PQP_REAL h1[3], h2[3], h3[3];
  PQP_REAL n1[3], m1[3];

  PQP_REAL ef11[3], ef12[3], ef13[3];
  PQP_REAL ef21[3], ef22[3], ef23[3];
  PQP_REAL ef31[3], ef32[3], ef33[3];
  
  p1[0] = P1[0] - P1[0];  p1[1] = P1[1] - P1[1];  p1[2] = P1[2] - P1[2];
  p2[0] = P2[0] - P1[0];  p2[1] = P2[1] - P1[1];  p2[2] = P2[2] - P1[2];
  p3[0] = P3[0] - P1[0];  p3[1] = P3[1] - P1[1];  p3[2] = P3[2] - P1[2];
  
  q1[0] = Q1[0] - P1[0];  q1[1] = Q1[1] - P1[1];  q1[2] = Q1[2] - P1[2];
  q2[0] = Q2[0] - P1[0];  q2[1] = Q2[1] - P1[1];  q2[2] = Q2[2] - P1[2];
  q3[0] = Q3[0] - P1[0];  q3[1] = Q3[1] - P1[1];  q3[2] = Q3[2] - P1[2];
  
  e1[0] = p2[0] - p1[0];  e1[1] = p2[1] - p1[1];  e1[2] = p2[2] - p1[2];
  e2[0] = p3[0] - p2[0];  e2[1] = p3[1] - p2[1];  e2[2] = p3[2] - p2[2];
  e3[0] = p1[0] - p3[0];  e3[1] = p1[1] - p3[1];  e3[2] = p1[2] - p3[2];

  f1[0] = q2[0] - q1[0];  f1[1] = q2[1] - q1[1];  f1[2] = q2[2] - q1[2];
  f2[0] = q3[0] - q2[0];  f2[1] = q3[1] - q2[1];  f2[2] = q3[2] - q2[2];
  f3[0] = q1[0] - q3[0];  f3[1] = q1[1] - q3[1];  f3[2] = q1[2] - q3[2];
  
  VcrossV(n1, e1, e2);
  VcrossV(m1, f1, f2);

  VcrossV(g1, e1, n1);
  VcrossV(g2, e2, n1);
  VcrossV(g3, e3, n1);
  VcrossV(h1, f1, m1);
  VcrossV(h2, f2, m1);
  VcrossV(h3, f3, m1);

  VcrossV(ef11, e1, f1);
  VcrossV(ef12, e1, f2);
  VcrossV(ef13, e1, f3);
  VcrossV(ef21, e2, f1);
  VcrossV(ef22, e2, f2);
  VcrossV(ef23, e2, f3);
  VcrossV(ef31, e3, f1);
  VcrossV(ef32, e3, f2);
  VcrossV(ef33, e3, f3);
  
  // now begin the series of tests

  if (!project6(n1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(m1, p1, p2, p3, q1, q2, q3)) return 0;
  
  if (!project6(ef11, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef12, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef13, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef21, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef22, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef23, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef31, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef32, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef33, p1, p2, p3, q1, q2, q3)) return 0;

  if (!project6(g1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(g2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(g3, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h3, p1, p2, p3, q1, q2, q3)) return 0;

  return 1;
}

PQP_REAL
TriDistance_Dirsection(PQP_REAL R[3][3], PQP_REAL T[3], PQP_Tri *t1, PQP_Tri *t2,
			PQP_REAL p[3], PQP_REAL q[3], PQP_REAL direction[3])
{
	// transform tri 2 into same space as tri 1

	PQP_REAL tri1[3][3], tri2[3][3];

	VcV(tri1[0], t1->p1);
	VcV(tri1[1], t1->p2);
	VcV(tri1[2], t1->p3);
	MxVpV(tri2[0], R, t2->p1, T);
	MxVpV(tri2[1], R, t2->p2, T);
	MxVpV(tri2[2], R, t2->p3, T);



	return TriDist(p,q,tri1,tri2);
}

PQP_REAL
TriDistance(PQP_REAL R[3][3], PQP_REAL T[3], PQP_Tri *t1, PQP_Tri *t2,
            PQP_REAL p[3], PQP_REAL q[3])
{
  // transform tri 2 into same space as tri 1

  PQP_REAL tri1[3][3], tri2[3][3];

  VcV(tri1[0], t1->p1);
  VcV(tri1[1], t1->p2);
  VcV(tri1[2], t1->p3);
  MxVpV(tri2[0], R, t2->p1, T);
  MxVpV(tri2[1], R, t2->p2, T);
  MxVpV(tri2[2], R, t2->p3, T);
                                
  return TriDist(p,q,tri1,tri2);
}

PQP_REAL
TriDistance(PQP_REAL R[3][3], PQP_REAL T[3], PQP_Tri *t1, PQP_Tri *t2,
            PQP_REAL p[3], PQP_REAL q[3], ContactFeature &f1, ContactFeature &f2, bool &bCollided)
{
	// transform tri 2 into same space as tri 1
	
	PQP_REAL tri1[3][3], tri2[3][3];
	
	VcV(tri1[0], t1->p1);
	VcV(tri1[1], t1->p2);
	VcV(tri1[2], t1->p3);
	MxVpV(tri2[0], R, t2->p1, T);
	MxVpV(tri2[1], R, t2->p2, T);
	MxVpV(tri2[2], R, t2->p3, T);
	
	return TriDist(p,q,tri1,tri2, f1, f2, bCollided);
}

void
CollideRecurse(PQP_CollideResult *res,
               PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
               PQP_Model *o1, int b1, 
               PQP_Model *o2, int b2, int flag)
{
  // first thing, see if we're overlapping

  res->num_bv_tests++;

  if (!BV_Overlap(R, T, o1->child(b1), o2->child(b2))) return;

  // if we are, see if we test triangles next

  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();

  if (l1 && l2) 
  {
    res->num_tri_tests++;

#if 1
    // transform the points in b2 into space of b1, then compare

    PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];
    PQP_REAL q1[3], q2[3], q3[3];
    PQP_REAL *p1 = t1->p1;
    PQP_REAL *p2 = t1->p2;
    PQP_REAL *p3 = t1->p3;    
    MxVpV(q1, res->R, t2->p1, res->T);
    MxVpV(q2, res->R, t2->p2, res->T);
    MxVpV(q3, res->R, t2->p3, res->T);
    if (TriContact(p1, p2, p3, q1, q2, q3)) 
    {
      // add this to result

      res->Add(t1->id, t2->id);
	  o1->last_tri = t1;
	  o2->last_tri = t2;

    }
#else
    PQP_REAL p[3], q[3];

    PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];

    if (TriDistance(res->R,res->T,t1,t2,p,q) == 0.0)
    {
      // add this to result

      res->Add(t1->id, t2->id);
    }
#endif

    return;
  }

  // we dont, so decide whose children to visit next

  PQP_REAL sz1 = o1->child(b1)->GetSize();
  PQP_REAL sz2 = o2->child(b2)->GetSize();

  PQP_REAL Rc[3][3],Tc[3],Ttemp[3];
   if (l2 || (!l1 && (sz1 > sz2)))
	  {
		  int c1 = o1->child(b1)->first_child;
		  int c2 = c1 + 1;

		  MTxM(Rc,o1->child(c1)->R,R);
#if PQP_BV_TYPE & OBB_TYPE
		  VmV(Ttemp,T,o1->child(c1)->To);
#else
		  VmV(Ttemp,T,o1->child(c1)->Tr);
#endif
		  MTxV(Tc,o1->child(c1)->R,Ttemp);
		  CollideRecurse(res,Rc,Tc,o1,c1,o2,b2,flag);

		  if ((flag == PQP_FIRST_CONTACT) && (res->num_pairs > 0)) return;

		  MTxM(Rc,o1->child(c2)->R,R);
#if PQP_BV_TYPE & OBB_TYPE
		  VmV(Ttemp,T,o1->child(c2)->To);
#else
		  VmV(Ttemp,T,o1->child(c2)->Tr);
#endif
		  MTxV(Tc,o1->child(c2)->R,Ttemp);
		  CollideRecurse(res,Rc,Tc,o1,c2,o2,b2,flag);
	  }
	  else 
	  {
		  int c1 = o2->child(b2)->first_child;
		  int c2 = c1 + 1;

		  MxM(Rc,R,o2->child(c1)->R);
#if PQP_BV_TYPE & OBB_TYPE
		  MxVpV(Tc,R,o2->child(c1)->To,T);
#else
		  MxVpV(Tc,R,o2->child(c1)->Tr,T);
#endif
		  CollideRecurse(res,Rc,Tc,o1,b1,o2,c1,flag);

		  if ((flag == PQP_FIRST_CONTACT) && (res->num_pairs > 0)) return;

		  MxM(Rc,R,o2->child(c2)->R);
#if PQP_BV_TYPE & OBB_TYPE
		  MxVpV(Tc,R,o2->child(c2)->To,T);
#else
		  MxVpV(Tc,R,o2->child(c2)->Tr,T);
#endif
		  CollideRecurse(res,Rc,Tc,o1,b1,o2,c2,flag);
	  }
  
  return;
    
 
}

int 
PQP_Collide(PQP_CollideResult *res,
            PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
            PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
            int flag)
{
  double t1 = GetTime();

  // make sure that the models are built

  if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;
  if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;

  // clear the stats

  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // don't release the memory, but reset the num_pairs counter

  res->num_pairs = 0;

 

  
  // Okay, compute what transform [R,T] that takes us from cs1 to cs2.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part

  MTxM(res->R,R1,R2);
  PQP_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);  
  MTxV(res->T, R1, Ttemp);
  
  // compute the transform from o1->child(0) to o2->child(0)

  PQP_REAL Rtemp[3][3], R[3][3], T[3];

  MxM(Rtemp,res->R,o2->child(0)->R);
  MTxM(R,o1->child(0)->R,Rtemp);

#if PQP_BV_TYPE & OBB_TYPE
  MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
  VmV(Ttemp,Ttemp,o1->child(0)->To);
#else
  MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
  VmV(Ttemp,Ttemp,o1->child(0)->Tr);
#endif

  MTxV(T,o1->child(0)->R,Ttemp);

  // now start with both top level BVs  

  PQP_REAL *p1 = o1->last_tri->p1;
  PQP_REAL *p2 = o1->last_tri->p2;
  PQP_REAL *p3 = o1->last_tri->p3; 

  PQP_REAL q1[3], q2[3], q3[3];
  MxVpV(q1, res->R, o2->last_tri->p1, res->T);
  MxVpV(q2, res->R, o2->last_tri->p2, res->T);
  MxVpV(q3, res->R, o2->last_tri->p3, res->T);
  if ((flag == PQP_FIRST_CONTACT) &&TriContact(p1, p2, p3, q1, q2, q3)) 
  {
	  // add this to result

	  res->Add(o1->last_tri->id, o2->last_tri->id);

  }
  if ((flag == PQP_FIRST_CONTACT) && (res->num_pairs > 0)) return PQP_OK; 

  CollideRecurse(res,R,T,o1,0,o2,0,flag);
  
  double t2 = GetTime();
  res->query_time_secs = t2 - t1;
  
  return PQP_OK; 
}


// Use the absolute transformation
int 
PQP_Collide_Subtree(PQP_CollideResult *res,
            PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1, int bv1,
            PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2, int bv2,
            int flag)
{
  double t1 = GetTime();

  // make sure that the models are built

  if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;
  if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;

  // clear the stats

  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // don't release the memory, but reset the num_pairs counter

  res->num_pairs = 0;
  
  // Okay, compute what transform [R,T] that takes us from cs1 to cs2.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part

  MTxM(res->R,R1,R2);
  PQP_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);  
  MTxV(res->T, R1, Ttemp);
  
  // compute the transform from o1->child(0) to o2->child(0)

  PQP_REAL Rtemp[3][3], R[3][3], T[3];

  // MxM(Rtemp,res->R,o2->child(0)->R);
  MxM(Rtemp,res->R,o2->child(bv2)->R_abs);
  // MTxM(R,o1->child(0)->R,Rtemp);
  MTxM(R,o1->child(bv1)->R_abs,Rtemp);

#if PQP_BV_TYPE & OBB_TYPE
  // MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
  MxVpV(Ttemp,res->R,o2->child(bv2)->To_abs,res->T);
  // VmV(Ttemp,Ttemp,o1->child(0)->To);
  VmV(Ttemp,Ttemp,o1->child(bv1)->To_abs);
#else
  // MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
  MxVpV(Ttemp,res->R,o2->child(bv2)->Tr_abs,res->T);
  // VmV(Ttemp,Ttemp,o1->child(0)->Tr);
  VmV(Ttemp,Ttemp,o1->child(bv1)->Tr_abs);
#endif

  // MTxV(T,o1->child(0)->R,Ttemp);
  MTxV(T,o1->child(bv1)->R_abs,Ttemp);

  // now start with both top level BVs  

  // CollideRecurse(res,R,T,o1,0,o2,0,flag);
  CollideRecurse(res,R,T,o1,bv1,o2,bv2,flag);
  
  double t2 = GetTime();
  res->query_time_secs = t2 - t1;


 
  return PQP_OK; 
}


#if PQP_BV_TYPE & RSS_TYPE // distance/tolerance only available with RSS
                           // unless an OBB distance test is supplied in 
                           // BV.cpp

// DISTANCE STUFF
//
//--------------------------------------------------------------------------

void
DistanceRecurse(PQP_DistanceResult *res,
                PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
                PQP_Model *o1, int b1,
                PQP_Model *o2, int b2)
{
  PQP_REAL sz1 = o1->child(b1)->GetSize();
  PQP_REAL sz2 = o2->child(b2)->GetSize();
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();

  if (l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.

    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];

    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
  
    if (d < res->distance) 
    {
      res->distance = d;

	  res->t1 = t1->id;
	  res->t2 = t2->id;
	  
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed 
                               // into c.s. 2 later
      o1->last_tri = t1;
      o2->last_tri = t2;


    }

    return;
  }

  // First, perform distance tests on the children. Then traverse 
  // them recursively, but test the closer pair first, the further 
  // pair second.

  int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
  PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

  if (l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1

    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child+1;
    c2 = b2;
    
    MTxM(R1,o1->child(a1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(a1)->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,o1->child(a1)->R,Ttemp);

    MTxM(R2,o1->child(c1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(c1)->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,o1->child(c1)->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = o2->child(b2)->first_child;
    c1 = b1;
    c2 = o2->child(b2)->first_child+1;

    MxM(R1,R,o2->child(a2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,o2->child(a2)->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif

    MxM(R2,R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,o2->child(c2)->Tr,T);
#else
    MxVpV(T2,R,o2->child(c2)->To,T);
#endif
  }

  res->num_bv_tests += 2;
  PQP_REAL S[3];
  

  PQP_REAL d1 = BV_Distance(R1, T1, o1->child(a1), o2->child(a2),S);
  PQP_REAL d2 = BV_Distance(R2, T2, o1->child(c1), o2->child(c2),S);

  if (d2 < d1)
  {
    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      DistanceRecurse(res, R2, T2, o1, c1, o2, c2);      
    }

    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      DistanceRecurse(res, R1, T1, o1, a1, o2, a2);
    }
  }
  else 
  {
    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      DistanceRecurse(res, R1, T1, o1, a1, o2, a2);
    }

    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      DistanceRecurse(res, R2, T2, o1, c1, o2, c2);      
    }
  }
}




void
ContactRecurse(PQP_ContactResult *res,
                PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
                PQP_Model *o1, int b1,
                PQP_Model *o2, int b2)
{
	//return if the models are colliding
  if (res->modelcollided) 
  {
	  return;
  }
  PQP_REAL sz1 = o1->child(b1)->GetSize();
  PQP_REAL sz2 = o2->child(b2)->GetSize();
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();

  if (l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.

    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];

	ContactFeature f1, f2;	
    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,res->modelcollided);
  
	if (res->modelcollided) 
	{
		return;
	}
	
	if (d <= res->ctolerance) 
	{
		//add the contact feature if < ctolerance
		f1.tri = t1->id;
		f2.tri = t2->id;
		//huangxin 071126
		f1.bid = b1;
		f2.bid = b2;
		
		res->Add(f1, f2, d, p, q);
		res->closer_than_ctolerance = true;
		
		VcV(res->p1, p);         // p already in c.s. 1
		VcV(res->p2, q);         // q must be transformed 

		o1->last_tri = t1;
		o2->last_tri = t2;		
	}


    return;
  }

  // First, perform distance tests on the children. Then traverse 
  // them recursively, but test the closer pair first, the further 
  // pair second.

  int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
  PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

  if (l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1

    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child+1;
    c2 = b2;
    
    MTxM(R1,o1->child(a1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(a1)->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,o1->child(a1)->R,Ttemp);

    MTxM(R2,o1->child(c1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(c1)->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,o1->child(c1)->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = o2->child(b2)->first_child;
    c1 = b1;
    c2 = o2->child(b2)->first_child+1;

    MxM(R1,R,o2->child(a2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,o2->child(a2)->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif

    MxM(R2,R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,o2->child(c2)->Tr,T);
#else
    MxVpV(T2,R,o2->child(c2)->To,T);
#endif
  }

  res->num_bv_tests += 2;
  PQP_REAL S[3];
  
  PQP_REAL d1 = BV_Distance(R1, T1, o1->child(a1), o2->child(a2),S);
  PQP_REAL d2 = BV_Distance(R2, T2, o1->child(c1), o2->child(c2),S);

  if (d2 < d1)
  {
    if ((d2 < (res->ctolerance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R2, T2, o1, c1, o2, c2);      
    }
	
	if (res->modelcollided) 
	{
		return;
	}

    if ((d1 < (res->ctolerance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R1, T1, o1, a1, o2, a2);
    }
  }
  else 
  {
    if ((d1 < (res->ctolerance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R1, T1, o1, a1, o2, a2);
    }

	if (res->modelcollided) 
	{
		return;
	}

    if ((d2 < (res->ctolerance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R2, T2, o1, c1, o2, c2);      
    }
  }
}



void
DistanceQueueRecurse(PQP_DistanceResult *res, 
                     PQP_REAL R[3][3], PQP_REAL T[3],
                     PQP_Model *o1, int b1,
                     PQP_Model *o2, int b2)
{
  BVTQ bvtq(res->qsize);

  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R,R);
  VcV(min_test.T,T);

  while(1) 
  {  
    int l1 = o1->child(min_test.b1)->Leaf();
    int l2 = o2->child(min_test.b2)->Leaf();
    
    if (l1 && l2) 
    {  
      // both leaves.  Test the triangles beneath them.

      res->num_tri_tests++;

      PQP_REAL p[3], q[3];

      PQP_Tri *t1 = &o1->tris[-o1->child(min_test.b1)->first_child - 1];
      PQP_Tri *t2 = &o2->tris[-o2->child(min_test.b2)->first_child - 1];

      PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
  
      if (d < res->distance)
      {
        res->distance = d;

		res->t1 = t1->id;
		res->t2 = t2->id;
		
        VcV(res->p1, p);         // p already in c.s. 1
        VcV(res->p2, q);         // q must be transformed 
                                 // into c.s. 2 later
        o1->last_tri = t1;
        o2->last_tri = t2;
      }
    }		 
    else if (bvtq.GetNumTests() == bvtq.GetSize() - 1) 
    {  
      // queue can't get two more tests, recur
      
      DistanceQueueRecurse(res,min_test.R,min_test.T,
                           o1,min_test.b1,o2,min_test.b2);
    }
    else 
    {  
      // decide how to descend to children
      
      PQP_REAL sz1 = o1->child(min_test.b1)->GetSize();
      PQP_REAL sz2 = o2->child(min_test.b2)->GetSize();

      res->num_bv_tests += 2;
 
      BVT bvt1,bvt2;
      PQP_REAL Ttemp[3];

      if (l2 || (!l1 && (sz1 > sz2)))	
      {  
        // put new tests on queue consisting of min_test.b2 
        // with children of min_test.b1 
      
        int c1 = o1->child(min_test.b1)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R,o1->child(c1)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,o1->child(c1)->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c1)->To);
#endif
        MTxV(bvt1.T,o1->child(c1)->R,Ttemp);
		PQP_REAL S[3];
		
        bvt1.d = BV_Distance(bvt1.R,bvt1.T,
                            o1->child(bvt1.b1),o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = c2;
        bvt2.b2 = min_test.b2;
        MTxM(bvt2.R,o1->child(c2)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,o1->child(c2)->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c2)->To);
#endif
        MTxV(bvt2.T,o1->child(c2)->R,Ttemp);
		
		
        bvt2.d = BV_Distance(bvt2.R,bvt2.T,
                            o1->child(bvt2.b1),o2->child(bvt2.b2),S);
      }
      else 
      {
        // put new tests on queue consisting of min_test.b1 
        // with children of min_test.b2
      
        int c1 = o2->child(min_test.b2)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R,min_test.R,o2->child(c1)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->Tr,min_test.T);
#else
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->To,min_test.T);
#endif
		PQP_REAL S[3];
		
        bvt1.d = BV_Distance(bvt1.R,bvt1.T,
                            o1->child(bvt1.b1),o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R,min_test.R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->Tr,min_test.T);
#else
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->To,min_test.T);
#endif
		
		
        bvt2.d = BV_Distance(bvt2.R,bvt2.T,
                            o1->child(bvt2.b1),o2->child(bvt2.b2),S);
      }

      bvtq.AddTest(bvt1);	
      bvtq.AddTest(bvt2);
    }

    if (bvtq.Empty())
    {
      break;
    }
    else
    {
      min_test = bvtq.ExtractMinTest();

      if ((min_test.d + res->abs_err >= res->distance) && 
         ((min_test.d * (1 + res->rel_err)) >= res->distance)) 
      {
        break;
      }
    }
  }  
}	

void
ContactQueueRecurse(PQP_ContactResult *res, 
                     PQP_REAL R[3][3], PQP_REAL T[3],
                     PQP_Model *o1, int b1,
                     PQP_Model *o2, int b2)
{
	if (res->modelcollided) 
	{
		return;
	}
  BVTQ bvtq(res->qsize);

  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R,R);
  VcV(min_test.T,T);

  while(1) 
  {  
    int l1 = o1->child(min_test.b1)->Leaf();
    int l2 = o2->child(min_test.b2)->Leaf();
    
    if (l1 && l2) 
    {  
      // both leaves.  Test the triangles beneath them.

      res->num_tri_tests++;

      PQP_REAL p[3], q[3];

      PQP_Tri *t1 = &o1->tris[-o1->child(min_test.b1)->first_child - 1];
      PQP_Tri *t2 = &o2->tris[-o2->child(min_test.b2)->first_child - 1];
	  
	  ContactFeature f1, f2;
      PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,res->modelcollided);
  
	  if (d <= res->ctolerance) 
	  {
		  //add the contact feature if < ctolerance
		  f1.tri = t1->id;
		  f2.tri = t2->id;
		  res->Add(f1, f2, d, p, q);
		  res->closer_than_ctolerance = true;
		  
		  VcV(res->p1, p);         // p already in c.s. 1
		  VcV(res->p2, q);         // q must be transformed 
		  
		  o1->last_tri = t1;
		  o2->last_tri = t2;		  
	  }
	  
//	  if (d < res->distance) 
	  //	  {
	  //		  //compute the minimum distance
	  //		  res->distance = d;
	  //		  
	  //		  // into c.s. 2 later
	  //		  o1->last_tri = t1;
	  //		  o2->last_tri = t2;
	  //		  	  
	  //	  }
	  

    }		 
    else if (bvtq.GetNumTests() == bvtq.GetSize() - 1) 
    {  
      // queue can't get two more tests, recur
      
      ContactQueueRecurse(res,min_test.R,min_test.T,
                           o1,min_test.b1,o2,min_test.b2);
	  if (res->modelcollided) 
	  {
		  return;
	  }
    }
    else 
    {  
      // decide how to descend to children
      
      PQP_REAL sz1 = o1->child(min_test.b1)->GetSize();
      PQP_REAL sz2 = o2->child(min_test.b2)->GetSize();

      res->num_bv_tests += 2;
 
      BVT bvt1,bvt2;
      PQP_REAL Ttemp[3];

      if (l2 || (!l1 && (sz1 > sz2)))	
      {  
        // put new tests on queue consisting of min_test.b2 
        // with children of min_test.b1 
      
        int c1 = o1->child(min_test.b1)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R,o1->child(c1)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,o1->child(c1)->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c1)->To);
#endif
        MTxV(bvt1.T,o1->child(c1)->R,Ttemp);
		PQP_REAL S[3];
		
        bvt1.d = BV_Distance(bvt1.R,bvt1.T,
                            o1->child(bvt1.b1),o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = c2;
        bvt2.b2 = min_test.b2;
        MTxM(bvt2.R,o1->child(c2)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,o1->child(c2)->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c2)->To);
#endif
        MTxV(bvt2.T,o1->child(c2)->R,Ttemp);
		
		
        bvt2.d = BV_Distance(bvt2.R,bvt2.T,
                            o1->child(bvt2.b1),o2->child(bvt2.b2),S);
      }
      else 
      {
        // put new tests on queue consisting of min_test.b1 
        // with children of min_test.b2
      
        int c1 = o2->child(min_test.b2)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R,min_test.R,o2->child(c1)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->Tr,min_test.T);
#else
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->To,min_test.T);
#endif
		PQP_REAL S[3];
		
        bvt1.d = BV_Distance(bvt1.R,bvt1.T,
                            o1->child(bvt1.b1),o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R,min_test.R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->Tr,min_test.T);
#else
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->To,min_test.T);
#endif
		
		
        bvt2.d = BV_Distance(bvt2.R,bvt2.T,
                            o1->child(bvt2.b1),o2->child(bvt2.b2),S);
      }

      bvtq.AddTest(bvt1);	
      bvtq.AddTest(bvt2);
    }

    if (bvtq.Empty())
    {
      break;
    }
    else
    {
      min_test = bvtq.ExtractMinTest();

      if ((min_test.d + res->abs_err >= res->ctolerance) && 
         ((min_test.d * (1 + res->rel_err)) >= res->ctolerance)) 
      {
        break;
      }
    }
  }  
}

void PrintTri(FILE* output, PQP_Tri &ftri)
{
	fprintf(output,"\n p1:%f, %f, %f", ftri.p1[0], ftri.p1[1],  ftri.p1[2]);
	fprintf(output,"\n p2:%f, %f, %f", ftri.p2[0], ftri.p2[1],  ftri.p2[2]);
	fprintf(output,"\n p3:%f, %f, %f", ftri.p3[0], ftri.p3[1],  ftri.p3[2]);
	
}

//
void DeleteDuplicatePairs(PQP_ContactResult *res, 
							PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
							PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2)
{
	int iDistNum = 0;
	int iPairsNum = res->num_cpairs;
	int *pPairIndex = new int[iPairsNum];
	
//	FILE* output1 = NULL;
	//	// fclose(fp);
	//	output1 = fopen("contactpairsPQP.txt", "a+");
	//	
	//	if (NULL == output1)
	//	{
	//		//			  cout << "Cannot open the write file: " << strCount << endl;
	//		return;		
	//	}
	
	
	int i,j;
	bool bDistinct;
	int maxPair = 100;

	if(iPairsNum > maxPair)
	{
		// When the parameter for the contact query is very large, the number iPairsNum could be huge: for instance 342000.
		// At this time, just truncate it. 
		static int iWarning = 0;
		if(iWarning < 10)
		{
			iWarning++;
			printf("Warning %d: Wrong with the paramters of contact query. The number of paris is %d\n", iWarning, iPairsNum);
		}

		iPairsNum = maxPair;
	}

	for(i = iPairsNum-1; i >=0; i--)
	{
		//the start of first iterative
		bDistinct = true;
		for(j = i - 1; j >= 0; j--)
		{
			//the start of second iterative
			//go on if the distance and type of the two pairs are the same
			if (fabs(res->cpairs[i].distance - res->cpairs[j].distance) < 1e-10) 
			{
				if ((res->cpairs[i].f1.type == res->cpairs[j].f1.type) &&
					(res->cpairs[i].f2.type == res->cpairs[j].f2.type)) 
				{
					// compare the case of V/F, F/V, E/E
					//1 the case of V/F
					if (res->cpairs[i].f1.type == 0) 
					{
						// compare whether the vertex and face are the same
						int iVert, jVert, iVertTri, jVertTri, iTri, jTri;
						iVert = res->cpairs[i].f1.fid[0];
						jVert = res->cpairs[j].f1.fid[0];
						iVertTri = res->cpairs[i].f1.tri;
						jVertTri = res->cpairs[j].f1.tri;

						PQP_REAL tri1[3][3], tri2[3][3];						
						VcV(tri1[0], o1->trisConst[iVertTri].p1);
						VcV(tri1[1], o1->trisConst[iVertTri].p2);
						VcV(tri1[2], o1->trisConst[iVertTri].p3);

						VcV(tri2[0], o1->trisConst[jVertTri].p1);
						VcV(tri2[1], o1->trisConst[jVertTri].p2);
						VcV(tri2[2], o1->trisConst[jVertTri].p3);
						
						iTri = res->cpairs[i].f2.tri;
						jTri = res->cpairs[j].f2.tri;
						if ((iTri == jTri) && (tri1[iVert][0] == tri2[jVert][0]) &&
							(tri1[iVert][1] == tri2[jVert][1] && (tri1[iVert][2] == tri2[jVert][2]))) 
						{
							bDistinct = false;
							break;
						}
//						fprintf(output1,"all the vertex:\n");

//						fprintf(output1,"The triangel ID is %d, %d, %d, %d\n", iVertTri, o1->trisConst[iVertTri].id, jVertTri, o1->trisConst[jVertTri].id);
//						PrintTri(output1, o1->tris[iVertTri]);
//						PrintTri(output1, o1->tris[jVertTri]);
						


//						fprintf(output1, "i:%d, j:%d, distance: %f, itri:%d = jTri:%d; tri12[iVert][0]:%f = %f; tri12[iVert][1]:%f = %f; tri12[iVert][2]:%f = %f.\n", i, j, res->cpairs[i].distance, iTri, jTri, 
//							tri1[iVert][0], tri2[jVert][0], tri1[iVert][1], tri2[jVert][1], tri1[iVert][2] , tri2[jVert][2]);

//						fprintf(output1,"\n");
					}

					//2 the case of F/V
					if (res->cpairs[i].f1.type == 2) 
					{
					    // compare whether the face and vertex are the same
						int iVert, jVert, iVertTri, jVertTri, iTri, jTri;
						iVert = res->cpairs[i].f2.fid[0];
						jVert = res->cpairs[j].f2.fid[0];
						iVertTri = res->cpairs[i].f2.tri;
						jVertTri = res->cpairs[j].f2.tri;
						
						PQP_REAL tri1[3][3], tri2[3][3];						
						VcV(tri1[0], o2->trisConst[iVertTri].p1);
						VcV(tri1[1], o2->trisConst[iVertTri].p2);
						VcV(tri1[2], o2->trisConst[iVertTri].p3);
						
						VcV(tri2[0], o2->trisConst[jVertTri].p1);
						VcV(tri2[1], o2->trisConst[jVertTri].p2);
						VcV(tri2[2], o2->trisConst[jVertTri].p3);
						
						iTri = res->cpairs[i].f1.tri;
						jTri = res->cpairs[j].f1.tri;
						if ((iTri == jTri) && (tri1[iVert][0] == tri2[jVert][0]) &&
							(tri1[iVert][1] == tri2[jVert][1] && (tri1[iVert][2] == tri2[jVert][2]))) 
						{
							bDistinct = false;
							break;
						}
						/*
						fprintf(output1, "i:%d, j:%d, distance: %f, itri:%d = jTri:%d; tri12[iVert][0]:%f = %f; tri12[iVert][1]:%f = %f; tri12[iVert][2]:%f = %f.\n", i, j, res->cpairs[i].distance, iTri, jTri, 
							tri1[iVert][0], tri2[jVert][0], tri1[iVert][1], tri2[jVert][1], tri1[iVert][2] , tri2[jVert][2]);
						fprintf(output1,"\n");
						*/
						
					}

					//3 the case of E/E
					if (res->cpairs[i].f1.type == 1) 
					{
						// compare whether the edge and edge are the same
						int iEdge1, iEdge2, iTri1, iTri2, jEdge1, jEdge2, jTri1, jTri2;
						iEdge1 = res->cpairs[i].f1.fid[0];
						iEdge2 = res->cpairs[i].f2.fid[0];
						iTri1 = res->cpairs[i].f1.tri;
						iTri2 = res->cpairs[i].f2.tri;

						jEdge1 = res->cpairs[j].f1.fid[0];
						jEdge2 = res->cpairs[j].f2.fid[0];												
						jTri1 = res->cpairs[j].f1.tri;
						jTri2 = res->cpairs[j].f2.tri;
						
						PQP_REAL iVtri1[3][3], iVtri2[3][3];						
						PQP_REAL jVtri1[3][3], jVtri2[3][3];		
						
						VcV(iVtri1[0], o1->trisConst[iTri1].p1);
						VcV(iVtri1[1], o1->trisConst[iTri1].p2);
						VcV(iVtri1[2], o1->trisConst[iTri1].p3);
						
						VcV(iVtri2[0], o2->trisConst[iTri2].p1);
						VcV(iVtri2[1], o2->trisConst[iTri2].p2);
						VcV(iVtri2[2], o2->trisConst[iTri2].p3);

						VcV(jVtri1[0], o1->trisConst[jTri1].p1);
						VcV(jVtri1[1], o1->trisConst[jTri1].p2);
						VcV(jVtri1[2], o1->trisConst[jTri1].p3);
						
						VcV(jVtri2[0], o2->trisConst[jTri2].p1);
						VcV(jVtri2[1], o2->trisConst[jTri2].p2);
						VcV(jVtri2[2], o2->trisConst[jTri2].p3);

//						bDistinct = false;
						
						
//						fprintf(output1,"all the vertex of edge:\n");
						
						//						fprintf(output1,"The triangel ID is %d, %d, %d, %d\n", iVertTri, o1->trisConst[iVertTri].id, jVertTri, o1->trisConst[jVertTri].id);
						//						PrintTri(output1, o1->tris[iVertTri]);
						//						PrintTri(output1, o1->tris[jVertTri]);
						
						
						
//						fprintf(output1, "i:%d, j:%d, distance: %f, itri:%d = jTri:%d; tri12[iVert][0]:%f = %f; tri12[iVert][1]:%f = %f; tri12[iVert][2]:%f = %f.\n", i, j, res->cpairs[i].distance, iTri, jTri, 
//						           tri1[iVert][0], tri2[jVert][0], tri1[iVert][1], tri2[jVert][1], tri1[iVert][2] , tri2[jVert][2]);
						
//						fprintf(output1,"\n");

						
						//compare the edge on model 1
						if ((((iVtri1[iEdge1][0] == jVtri1[jEdge1][0]) && (iVtri1[iEdge1][1] == jVtri1[jEdge1][1]) && (iVtri1[iEdge1][2] == jVtri1[jEdge1][2])) || 
											((iVtri1[(iEdge1+1)%3][0] == jVtri1[(jEdge1+1)%3][0]) && (iVtri1[(iEdge1+1)%3][1] == jVtri1[(jEdge1+1)%3][1]) && (iVtri1[(iEdge1+1)%3][2] == jVtri1[(jEdge1+1)%3][2]))) ||
											(((iVtri1[iEdge1][0] == jVtri1[(jEdge1+1)%3][0]) && (iVtri1[iEdge1][1] == jVtri1[(jEdge1+1)%3][1]) && (iVtri1[iEdge1][2] == jVtri1[(jEdge1+1)%3][2])) || 
											((iVtri1[(iEdge1+1)%3][0] == jVtri1[jEdge1][0]) && (iVtri1[(iEdge1+1)%3][1] == jVtri1[jEdge1][1]) && (iVtri1[(iEdge1+1)%3][2] == jVtri1[jEdge1][2]))))
								
						{
							//compare the edge on model 2
							if ((((iVtri2[iEdge2][0] == jVtri2[jEdge2][0]) && (iVtri2[iEdge2][1] == jVtri2[jEdge2][1]) && (iVtri2[iEdge2][2] == jVtri2[jEdge2][2])) || 
								((iVtri2[(iEdge2+1)%3][0] == jVtri2[(jEdge2+1)%3][0]) && (iVtri2[(iEdge2+1)%3][1] == jVtri2[(jEdge2+1)%3][1]) && (iVtri2[(iEdge2+1)%3][2] == jVtri2[(jEdge2+1)%3][2]))) ||
								(((iVtri2[iEdge2][0] == jVtri2[(jEdge2+1)%3][0]) && (iVtri2[iEdge2][1] == jVtri2[(jEdge2+1)%3][1]) && (iVtri2[iEdge2][2] == jVtri2[(jEdge2+1)%3][2])) || 
								((iVtri2[(iEdge2+1)%3][0] == jVtri2[jEdge2][0]) && (iVtri2[(iEdge2+1)%3][1] == jVtri2[jEdge2][1]) && (iVtri2[(iEdge2+1)%3][2] == jVtri2[jEdge2][2])))) 
							{
								bDistinct = false;
								break;
							}
							
						}

					}
				}
			}
						
		}
		//the end of second iterative

		if (bDistinct) 
		{
			pPairIndex[iDistNum++] = i;
		}
	}
	//the end of first iterative
	if (iDistNum) 
	{
		ContactPair * newpairs = new ContactPair[iDistNum];
		for(int k = 0; k < iDistNum; k++)
		{
			newpairs[k].f1 = res->cpairs[pPairIndex[k]].f1;
			newpairs[k].f2 = res->cpairs[pPairIndex[k]].f2;
			newpairs[k].distance = res->cpairs[pPairIndex[k]].distance;
			VcV(newpairs[k].p1, res->cpairs[pPairIndex[k]].p1);
			VcV(newpairs[k].p2, res->cpairs[pPairIndex[k]].p2);
			
		}
		delete [] res->cpairs;
		res->cpairs = newpairs;
		res->num_cpairs = iDistNum;
		res->num_cpairs_alloced = iDistNum;
		
	}
	delete []pPairIndex;

	if(bPQPContact_Represenative)
	{
		//report the representative pairs by distance clustering
		int repNum = 0;
		int *pRepIndex = new int[iDistNum];
		pRepIndex[0] = 0;
		repNum++;
		PQP_REAL p1[3], p2[3], p3[3], p4[3];
		PQP_REAL square_epislon = res->ctolerance * res->ctolerance * PQPContact_radius * PQPContact_radius;

		//cluster the pairs
		for( i = 1; i < iDistNum; i++)
		{
			// It is ok to use the coornidates relative to R1, T1. 
			// MxVpV(p1, R1, res->cpairs[i].p1, T1);
			// MxVpV(p2, R1, res->cpairs[i].p2, T1);

			VcV(p1, res->cpairs[i].p1);
			VcV(p2, res->cpairs[i].p2);

			bool bSame = false;
			for(int j = 0; j < repNum; j++)
			{
				// It is ok to use the coornidates relative to R1, T1. 
				//MxVpV(p3, R1, res->cpairs[pRepIndex[j]].p1, T1);
				//MxVpV(p4, R1, res->cpairs[pRepIndex[j]].p2, T1);

				VcV(p3, res->cpairs[pRepIndex[j]].p1);
				VcV(p4, res->cpairs[pRepIndex[j]].p2);

				// if ((VdistV2(p1, p3) <= epislon ) || (VdistV2(p2, p4) <= epislon )) 
				if (  (VdistV2(p1, p3) <= square_epislon ) && (VdistV2(p2, p4) <= square_epislon ) )
				{
	//				pRepIndex[repNum++] = i;
	//				break;
					bSame = true;
					double d1 = res->cpairs[i].distance;
					double d2 = res->cpairs[pRepIndex[j]].distance;

					if (d1 < d2)
					{
						pRepIndex[j] = i;
					}

					break;
				}
			}

			if (!bSame) 
			{
				pRepIndex[repNum++] = i;			
			}

		}

		//copy the representative pairs
		ContactPair * reppairs = new ContactPair[repNum];	
		for(i = 0; i < repNum; i++)
		{
			reppairs[i].f1 = res->cpairs[pRepIndex[i]].f1;
			reppairs[i].f2 = res->cpairs[pRepIndex[i]].f2;
			reppairs[i].distance = res->cpairs[pRepIndex[i]].distance;
			VcV(reppairs[i].p1, res->cpairs[pRepIndex[i]].p1);
			VcV(reppairs[i].p2, res->cpairs[pRepIndex[i]].p2);
		}

		delete [] res->cpairs;
		res->cpairs = reppairs;
		res->num_cpairs = repNum;
		res->num_cpairs_alloced = repNum;
		delete []pRepIndex;
	}
}

int 
PQP_Contact(PQP_ContactResult *res,
             PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
             PQP_REAL rel_err, PQP_REAL abs_err, PQP_REAL tolerance,
             int qsize)
{
	
	double time1 = GetTime();
	
	// make sure that the models are built
	
	if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;
	if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;
	
	// Okay, compute what transform [R,T] that takes us from cs2 to cs1.
	// [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
	// First compute the rotation part, then translation part
	
	MTxM(res->R,R1,R2);
	PQP_REAL Ttemp[3];
	VmV(Ttemp, T2, T1);  
	MTxV(res->T, R1, Ttemp);

	res->closer_than_ctolerance = 0;
	if (tolerance < 0.0) tolerance = 0;
	res->ctolerance = tolerance;
	res->modelcollided = false;
	
	// establish initial upper bound using last triangles which 
	// provided the minimum distance
	//	res->t1 = o1->last_tri->id;
	//	res->t2 = o2->last_tri->id;
//	ContactFeature f1, f2;	
//	PQP_REAL p[3],q[3];
//	res->distance = TriDistance(res->R,res->T,o1->last_tri,o2->last_tri,p,q, f1, f2,res->modelcollided);
	//	
	//	if res->distance < tolerance) 
	//	{
	//
	//		f1.tri = o1->last_tri->id;
	//		f2.tri = o2->last_tri->id;
	//		res->Add(f1, f2, res->distance);		
	//		res->closer_than_ctolerance = true;
	//
	//		VcV(res->p1,p);
	//		VcV(res->p2,q);
	//	}
	
//	res->distance = tolerance;
	
	// initialize error bounds
	
	res->abs_err = abs_err;
	res->rel_err = rel_err;
	
	// clear the stats
	
	res->num_bv_tests = 0;
	res->num_tri_tests = 0;
	
	// compute the transform from o1->child(0) to o2->child(0)
	
	PQP_REAL Rtemp[3][3], R[3][3], T[3];
	
	MxM(Rtemp,res->R,o2->child(0)->R);
	MTxM(R,o1->child(0)->R,Rtemp);
	
#if PQP_BV_TYPE & RSS_TYPE
	MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->Tr);
#else
	MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->To);
#endif
	MTxV(T,o1->child(0)->R,Ttemp);
	
	// choose routine according to queue size
	
	if (qsize <= 2)
	{
		ContactRecurse(res,R,T,o1,0,o2,0);    
	}
	else 
	{ 
		res->qsize = qsize;
		
		ContactQueueRecurse(res,R,T,o1,0,o2,0);
	}
	
	// res->p2 is in cs 1 ; transform it to cs 2
	
	PQP_REAL u[3];
	VmV(u, res->p2, res->T);
	MTxV(res->p2, res->R, u);
	
	//delete the duplicate pairs
	// if (res->CloserThanTolerance() && !res->modelcollided) 
	if (res->CloserThanTolerance() && !res->modelcollided) 
	// Change by Liangjun. Jan 12, 2007
	{
		// Always delete the duplicated pairs
		DeleteDuplicatePairs(res, R1, T1, o1, R2, T2, o2);
	}
	else
	{
		res->FreePairsList();
	}

	double time2 = GetTime();
	res->query_time_secs = time2 - time1;  
	
	return PQP_OK;
}

int 
PQP_Distance(PQP_DistanceResult *res,
             PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
             PQP_REAL rel_err, PQP_REAL abs_err,
             int qsize)
{
  
  double time1 = GetTime();
  
  // make sure that the models are built

  if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;
  if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;

  // Okay, compute what transform [R,T] that takes us from cs2 to cs1.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part

  MTxM(res->R,R1,R2);
  PQP_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);  
  MTxV(res->T, R1, Ttemp);
  
  // establish initial upper bound using last triangles which 
  // provided the minimum distance

  PQP_REAL p[3],q[3];
  res->distance = TriDistance(res->R,res->T,o1->last_tri,o2->last_tri,p,q);
  res->t1 = o1->last_tri->id;
  res->t2 = o2->last_tri->id;
		
  VcV(res->p1,p);
  VcV(res->p2,q);

  // initialize error bounds

  res->abs_err = abs_err;
  res->rel_err = rel_err;
  
  // clear the stats

  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // compute the transform from o1->child(0) to o2->child(0)

  PQP_REAL Rtemp[3][3], R[3][3], T[3];

  MxM(Rtemp,res->R,o2->child(0)->R);
  MTxM(R,o1->child(0)->R,Rtemp);
  
#if PQP_BV_TYPE & RSS_TYPE
  MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
  VmV(Ttemp,Ttemp,o1->child(0)->Tr);
#else
  MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
  VmV(Ttemp,Ttemp,o1->child(0)->To);
#endif
  MTxV(T,o1->child(0)->R,Ttemp);

  // choose routine according to queue size

  if (qsize <= 2)
  {
    DistanceRecurse(res,R,T,o1,0,o2,0);    
  }
  else 
  { 
    res->qsize = qsize;

    DistanceQueueRecurse(res,R,T,o1,0,o2,0);
  }

  // res->p2 is in cs 1 ; transform it to cs 2

  PQP_REAL u[3];
  VmV(u, res->p2, res->T);
  MTxV(res->p2, res->R, u);

  double time2 = GetTime();
  res->query_time_secs = time2 - time1;  

  return PQP_OK;
}

// Tolerance Stuff
//
//---------------------------------------------------------------------------
void 
ToleranceRecurse(PQP_ToleranceResult *res, 
                 PQP_REAL R[3][3], PQP_REAL T[3],
                 PQP_Model *o1, int b1, PQP_Model *o2, int b2)
{
  PQP_REAL sz1 = o1->child(b1)->GetSize();
  PQP_REAL sz2 = o2->child(b2)->GetSize();
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();

  if (l1 && l2) 
  {
    // both leaves - find if tri pair within tolerance
    
    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];

    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
    
    if (d <= res->tolerance)  
    {  
      // triangle pair distance less than tolerance

      res->closer_than_tolerance = 1;
      res->distance = d;
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed 
                               // into c.s. 2 later
    }

    return;
  }

  int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
  PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

  if (l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1

    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child+1;
    c2 = b2;
    
    MTxM(R1,o1->child(a1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(a1)->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,o1->child(a1)->R,Ttemp);

    MTxM(R2,o1->child(c1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(c1)->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,o1->child(c1)->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = o2->child(b2)->first_child;
    c1 = b1;
    c2 = o2->child(b2)->first_child+1;

    MxM(R1,R,o2->child(a2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,o2->child(a2)->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif
    MxM(R2,R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,o2->child(c2)->Tr,T);
#else
    MxVpV(T2,R,o2->child(c2)->To,T);
#endif
  }

  res->num_bv_tests += 2;
  PQP_REAL S[3];
  
  PQP_REAL d1 = BV_Distance(R1, T1, o1->child(a1), o2->child(a2),S);
  PQP_REAL d2 = BV_Distance(R2, T2, o1->child(c1), o2->child(c2),S);

  if (d2 < d1) 
  {
    if (d2 <= res->tolerance) ToleranceRecurse(res, R2, T2, o1, c1, o2, c2);
    if (res->closer_than_tolerance) return;
    if (d1 <= res->tolerance) ToleranceRecurse(res, R1, T1, o1, a1, o2, a2);
  }
  else 
  {
    if (d1 <= res->tolerance) ToleranceRecurse(res, R1, T1, o1, a1, o2, a2);
    if (res->closer_than_tolerance) return;
    if (d2 <= res->tolerance) ToleranceRecurse(res, R2, T2, o1, c1, o2, c2);
  }
}

void
ToleranceQueueRecurse(PQP_ToleranceResult *res,
                      PQP_REAL R[3][3], PQP_REAL T[3],
                      PQP_Model *o1, int b1,
                      PQP_Model *o2, int b2)
{
  BVTQ bvtq(res->qsize);
  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R,R);
  VcV(min_test.T,T);

  while(1)
  {  
    int l1 = o1->child(min_test.b1)->Leaf();
    int l2 = o2->child(min_test.b2)->Leaf();
    
    if (l1 && l2) 
    {  
      // both leaves - find if tri pair within tolerance
    
      res->num_tri_tests++;

      PQP_REAL p[3], q[3];

      PQP_Tri *t1 = &o1->tris[-o1->child(min_test.b1)->first_child - 1];
      PQP_Tri *t2 = &o2->tris[-o2->child(min_test.b2)->first_child - 1];

      PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
    
      if (d <= res->tolerance)  
      {  
        // triangle pair distance less than tolerance

        res->closer_than_tolerance = 1;
        res->distance = d;
        VcV(res->p1, p);         // p already in c.s. 1
        VcV(res->p2, q);         // q must be transformed 
                                 // into c.s. 2 later
        return;
      }
    }
    else if (bvtq.GetNumTests() == bvtq.GetSize() - 1)
    {  
      // queue can't get two more tests, recur
      
      ToleranceQueueRecurse(res,min_test.R,min_test.T,
                            o1,min_test.b1,o2,min_test.b2);
      if (res->closer_than_tolerance == 1) return;
    }
    else 
    {  
      // decide how to descend to children
      
      PQP_REAL sz1 = o1->child(min_test.b1)->GetSize();
      PQP_REAL sz2 = o2->child(min_test.b2)->GetSize();

      res->num_bv_tests += 2;
      
      BVT bvt1,bvt2;
      PQP_REAL Ttemp[3];

      if (l2 || (!l1 && (sz1 > sz2)))	
      {
	      // add two new tests to queue, consisting of min_test.b2
        // with the children of min_test.b1

        int c1 = o1->child(min_test.b1)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R,o1->child(c1)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,o1->child(c1)->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c1)->To);
#endif
		PQP_REAL S[3];
		
        MTxV(bvt1.T,o1->child(c1)->R,Ttemp);
        bvt1.d = BV_Distance(bvt1.R,bvt1.T,
                            o1->child(bvt1.b1),o2->child(bvt1.b2),S);

	      // init bv test 2

	      bvt2.b1 = c2;
	      bvt2.b2 = min_test.b2;
	      MTxM(bvt2.R,o1->child(c2)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
	      VmV(Ttemp,min_test.T,o1->child(c2)->Tr);
#else
	      VmV(Ttemp,min_test.T,o1->child(c2)->To);
#endif
	      MTxV(bvt2.T,o1->child(c2)->R,Ttemp);
        bvt2.d = BV_Distance(bvt2.R,bvt2.T,
                            o1->child(bvt2.b1),o2->child(bvt2.b2),S);
      }
      else 
      {
        // add two new tests to queue, consisting of min_test.b1
        // with the children of min_test.b2

        int c1 = o2->child(min_test.b2)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R,min_test.R,o2->child(c1)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->Tr,min_test.T);
#else
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->To,min_test.T);
#endif
		PQP_REAL S[3];
		
        bvt1.d = BV_Distance(bvt1.R,bvt1.T,
                            o1->child(bvt1.b1),o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R,min_test.R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->Tr,min_test.T);
#else
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->To,min_test.T);
#endif
        bvt2.d = BV_Distance(bvt2.R,bvt2.T,
                            o1->child(bvt2.b1),o2->child(bvt2.b2),S);
      }

      // put children tests in queue

      if (bvt1.d <= res->tolerance) bvtq.AddTest(bvt1);
      if (bvt2.d <= res->tolerance) bvtq.AddTest(bvt2);
    }

    if (bvtq.Empty() || (bvtq.MinTest() > res->tolerance)) 
    {
      res->closer_than_tolerance = 0;
      return;
    }
    else 
    {
      min_test = bvtq.ExtractMinTest();
    }
  }  
}	

int
PQP_Tolerance(PQP_ToleranceResult *res,
              PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
              PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
              PQP_REAL tolerance,
              int qsize)
{
	double time1 = GetTime();
	
	// make sure that the models are built
	
	if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;
	if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;
	
	// Compute the transform [R,T] that takes us from cs2 to cs1.
	// [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
	
	MTxM(res->R,R1,R2);
	PQP_REAL Ttemp[3];
	VmV(Ttemp, T2, T1);
	MTxV(res->T, R1, Ttemp);
	
	// set tolerance, used to prune the search
	
	if (tolerance < 0.0) tolerance = 0.0;
	res->tolerance = tolerance;
	
	// clear the stats
	
	res->num_bv_tests = 0;
	res->num_tri_tests = 0;
	
	// initially assume not closer than tolerance
	
	res->closer_than_tolerance = 0;
	
	// compute the transform from o1->child(0) to o2->child(0)
	
	PQP_REAL Rtemp[3][3], R[3][3], T[3];
	PQP_REAL S[3];
	
	
	MxM(Rtemp,res->R,o2->child(0)->R);
	MTxM(R,o1->child(0)->R,Rtemp);
#if PQP_BV_TYPE & RSS_TYPE
	MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->Tr);
#else
	MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->To);
#endif
	MTxV(T,o1->child(0)->R,Ttemp);
	
	// find a distance lower bound for trivial reject
	
	PQP_REAL d = BV_Distance(R, T, o1->child(0), o2->child(0),S);
	
	if (d <= res->tolerance)
	{
		// more work needed - choose routine according to queue size
		
		if (qsize <= 2) 
		{
			ToleranceRecurse(res, R, T, o1, 0, o2, 0);
		}
		else 
		{
			res->qsize = qsize;
			ToleranceQueueRecurse(res, R, T, o1, 0, o2, 0);
		}
	}
	
	// res->p2 is in cs 1 ; transform it to cs 2
	
	PQP_REAL u[3];
	VmV(u, res->p2, res->T);
	MTxV(res->p2, res->R, u);
	
	double time2 = GetTime();
	res->query_time_secs = time2 - time1;
	
	return PQP_OK;
}

PQP_ContactResult::PQP_ContactResult()
{
	cpairs = 0;
	num_cpairs = num_cpairs_alloced = 0;
//	num_bv_tests = 0;
//	num_tri_tests = 0;
}

PQP_ContactResult::~PQP_ContactResult()
{
	delete [] cpairs;
}

void
PQP_ContactResult::FreePairsList()
{
	num_cpairs = num_cpairs_alloced = 0;
	delete [] cpairs;
	cpairs = 0;
}
// may increase OR reduce mem usage
void
PQP_ContactResult::SizeTo(int n)
{
	ContactPair *temp;
	
	if (n < num_cpairs) 
	{
		fprintf(stderr, "PQP Error: Internal error in "
			"'PQP_ContactResult::SizeTo(int n)'\n");
		fprintf(stderr, "       n = %d, but num_cpairs = %d\n", n, num_cpairs);
		return;
	}
	
	temp = new ContactPair[n];
	memcpy(temp, cpairs, num_cpairs*sizeof(ContactPair));
	delete [] cpairs;
	cpairs = temp;
	num_cpairs_alloced = n;
	return;
}

void PQP_ContactResult::Add(ContactFeature f1, ContactFeature f2, PQP_REAL distance, PQP_REAL p[], PQP_REAL q[])
{
	// Check whether a feature pair is duplicated or not
	bool bDuplicated = false;
	if(bPQPContact_Represenative)
	{
		// The tolerance $\kappa_r$
		PQP_REAL square_epislon = ctolerance * ctolerance * PQPContact_radius * PQPContact_radius;

		//report the representative pairs by distance clustering
		for( int i = 0; i < num_cpairs; i++)
		{
			// It is ok to use the coornidates relative to R1, T1. 
			// MxVpV(p1, R1, res->cpairs[i].p1, T1);
			// MxVpV(p2, R1, res->cpairs[i].p2, T1);

			// if ((VdistV2(p1, p3) <= epislon ) || (VdistV2(p2, p4) <= epislon )) 
			if (  (VdistV2(p, cpairs[i].p1) <= square_epislon ) 
			   && (VdistV2(q, cpairs[i].p2) <= square_epislon ) )
			{
				bDuplicated = true;

				// replace the item by the new item
				if (cpairs[i].distance > distance)
				{
					cpairs[i].f1 = f1;
					cpairs[i].f2 = f2;
					cpairs[i].distance = distance;
					VcV(cpairs[i].p1, p);
					VcV(cpairs[i].p2, q);
				}
			}
		}
	}

	// now proceed as usual
	if(!bDuplicated)
	{
		if (num_cpairs >= num_cpairs_alloced) 
		{
			// allocate more
			
			SizeTo(num_cpairs_alloced*2+8);
		}
		cpairs[num_cpairs].f1 = f1;
		cpairs[num_cpairs].f2 = f2;
		cpairs[num_cpairs].distance = distance;
		VcV(cpairs[num_cpairs].p1, p);
		VcV(cpairs[num_cpairs].p2, q);
	
		num_cpairs++;
	}
}

#endif

