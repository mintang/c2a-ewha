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
/*************************************************************************\

  Copyright 2009	Computer Graphics Lab, 
									Ewha Womans University,
									Seoul, Korea.
  All Rights Reserved.

  The authors may be contacted via:

  Mailing Address:     Xinyu Zhang
                       Min Tang
                       Department of Computer Science & Engineering
                       Ewha Womans University
                       Seoul, 120-750, Korea.

  EMail:               zhangxy@ewha.ac.kr
                       tangmin@ewha.ac.kr


\**************************************************************************/

//#include "SWIFT_linalg.h"
#include "PQP_Tri.h"
#include "BV.h"
#include <list>
#include <string>
#include <algorithm>


using namespace std;

class PQP_Model
{
public:

  int build_state;

  PQP_Tri *tris;  
  PQP_Tri *trisConst; //the index of the tri will not change, add by huangxin 11/10/2007
  int num_tris;
  int num_tris_alloced;
  PQP_REAL Vtemp[3];

  BV *b;
  int num_bvs;
  int num_bvs_alloced;
  int maxdeep;
  PQP_REAL P[3];
  PQP_REAL Q[3];

  //use a single tri coherence
 bool bMotionCoherence;
 int motionTri;
 int motionBV;

 //use a local tree coherence, set local BV level typicall 3
 int level; //set local BV level typicall 3
 int levelPID;//the parent ID at level tree from the leaf
// int numLocalTris; // pow(2, level)
// int *localTris;//all the tris in the local BV tree

  PQP_Tri *last_tri;       // closest tri on this model in last distance test
  
  BV *child(int n) { return &b[n]; }

  PQP_Model();
  ~PQP_Model();

  int BeginModel(int num_tris = 8); // preallocate for num_tris triangles;
                                    // the parameter is optional, since
                                    // arrays are reallocated as needed
  int AddTri(const PQP_REAL *p1, const PQP_REAL *p2, const PQP_REAL *p3, 
             int i, int i1, int i2, int i3);
  int AddTri(const PQP_REAL *p1, const PQP_REAL *p2, const PQP_REAL *p3, 
	  int i);
  int EndModel();
  int MemUsage(int msg);  // returns model mem usage.  
                          // prints message to stderr if msg == TRUE
  
  PQP_REAL com[3]; // the center of Mass, huangxin
  PQP_REAL radius; // the radius of the model using com, huangxin

  PQP_REAL cwLength;
  PQP_REAL v1v2DifLength;

  
  PQP_REAL Max_Cross_Product(PQP_REAL dir[3]);
  void SetCenterOfMass(PQP_REAL c[3]);
  void ComputeCenterOfMass();
  void ComputeRadius();

};


// Add by Liangjun: to compute the aboslute transform
void ComputeAbsTra(PQP_REAL R[3][3], PQP_REAL T[3], PQP_Model *o, int bv_ind);


struct CollisionPair
{
  int id1;
  int id2;
};

struct PQP_CollideResult  
{
  // stats

  int num_bv_tests;
  int num_tri_tests;
  int num_frontnodes;
  double query_time_secs;

  // xform from model 1 to model 2

  PQP_REAL R[3][3];
  PQP_REAL T[3];

  int num_pairs_alloced;
  int num_pairs;
  CollisionPair *pairs;

  void SizeTo(int n);    
  void Add(int i1, int i2); 

  PQP_CollideResult();
  ~PQP_CollideResult();

  // statistics

  int NumBVTests() { return num_bv_tests; }
  int NumTriTests() { return num_tri_tests; }
  double QueryTimeSecs() { return query_time_secs; }

  // free the list of contact pairs; ordinarily this list is reused
  // for each query, and only deleted in the destructor.

  void FreePairsList(); 

  // query results

  int Colliding() { return (num_pairs > 0); }
  int NumPairs() { return num_pairs; }
  int Id1(int k) { return pairs[k].id1; }
  int Id2(int k) { return pairs[k].id2; }
};

#if PQP_BV_TYPE & RSS_TYPE // distance/tolerance are only available with RSS

struct PQP_DistanceResult 
{
  // stats

  int num_bv_tests;
  int num_tri_tests;
  double query_time_secs;

  // xform from model 1 to model 2

  PQP_REAL R[3][3];
  PQP_REAL T[3];

  PQP_REAL rel_err; 
  PQP_REAL abs_err; 

  PQP_REAL distance;
  PQP_REAL p1[3]; 
  PQP_REAL p2[3];
  int qsize;
  
  int t1;//the closest pair of t1
  int t2;//the closest pair of t2
  
  // statistics

  int NumBVTests() { return num_bv_tests; }
  int NumTriTests() { return num_tri_tests; }
  double QueryTimeSecs() { return query_time_secs; }

  // The following distance and points established the minimum distance
  // for the models, within the relative and absolute error bounds 
  // specified.
  // Points are defined: PQP_REAL p1[3], p2[3];

  PQP_REAL Distance() { return distance; }
  const PQP_REAL *P1() { return p1; }
  const PQP_REAL *P2() { return p2; }
};

struct ContactFeature
{
	int type;//0: vertex; 1: edge; 2: face
	int tri; //index of the triangle
	int fid[3]; //index of the vertex or edge according to the type
	int bid; //the BV of the triangle 		//huangxin 071126	
};
struct ContactPair
{
	ContactFeature f1;
	ContactFeature f2;
	double distance; // the diatnce between the two features
	PQP_REAL p1[3]; //the nearest point of f1 to f2
	PQP_REAL p2[3]; //the nearest point of f2 to f1
};
struct PQP_ContactResult 
{
	// stats
	
	int num_bv_tests;
	int num_tri_tests;
	double query_time_secs;
	
	// xform from model 1 to model 2
	
	PQP_REAL R[3][3];
	PQP_REAL T[3];
	
	PQP_REAL rel_err; 
	PQP_REAL abs_err; 
	
	
	
//	PQP_REAL distance;
	PQP_REAL p1[3]; 
	PQP_REAL p2[3];
	int qsize;
	
	//contact features
	ContactPair *cpairs;
	int num_cpairs_alloced;
	int num_cpairs;
	bool closer_than_ctolerance;
	bool modelcollided; //whether the two models are colliding
	PQP_REAL ctolerance;      
	
	void FreePairsList(); 
	void SizeTo(int n);
	void Add(ContactFeature f1, ContactFeature f2, PQP_REAL distance, PQP_REAL p[], PQP_REAL q[]); 
	
	// statistics
	
	int NumBVTests() { return num_bv_tests; }
	int NumTriTests() { return num_tri_tests; }
	double QueryTimeSecs() { return query_time_secs; }
	
	// The following distance and points established the minimum distance
	// for the models, within the relative and absolute error bounds 
	// specified.
	// Points are defined: PQP_REAL p1[3], p2[3];
	PQP_ContactResult();
	~PQP_ContactResult();
	bool CloserThanTolerance() { return closer_than_ctolerance; }
	
	
//	PQP_REAL Distance() { return distance; }
	const PQP_REAL *P1() { return p1; }
	const PQP_REAL *P2() { return p2; }
};


struct PQP_ToleranceResult 
{
	// stats
	
	int num_bv_tests;
	int num_tri_tests;
	double query_time_secs;
	
	// xform from model 1 to model 2
	
	PQP_REAL R[3][3];
	PQP_REAL T[3];
	
	int    closer_than_tolerance;   
	PQP_REAL tolerance;      
	
	PQP_REAL distance;
	PQP_REAL p1[3]; 
	PQP_REAL p2[3]; 
	int qsize;
	
	// statistics
	
	int NumBVTests() { return num_bv_tests; }
	int NumTriTests() { return num_tri_tests; }
	double QueryTimeSecs() { return query_time_secs; }
	
	// If the models are closer than ( <= ) tolerance, these points 
	// and distance were what established this.  Otherwise, 
	// distance and point values are not meaningful.
	
	PQP_REAL Distance() { return distance; }
	const PQP_REAL *P1() { return p1; }
	const PQP_REAL *P2() { return p2; }
	
	// boolean says whether models are closer than tolerance distance
	
	int CloserThanTolerance() { return closer_than_tolerance; }
};

struct ContactF
{
public:
	void IcI(int a[3],int b[3]);
	void VcV_L(PQP_REAL a[3],PQP_REAL b[3]);
	ContactF(int FeatureType_A, int FeatureType_B, int FeatureID_A[3], int FeatureID_B[3],
		int TriangleID_A,	int TriangleID_B,	PQP_REAL P_A[3],	PQP_REAL P_B[3],
		PQP_REAL Distance);
public:
	int FeatureType_A;
	int FeatureType_B;
	int FeatureID_A[3];
	int FeatureID_B[3];
	int TriangleID_A;
	int TriangleID_B;
	PQP_REAL P_A[3];
	PQP_REAL P_B[3];
	PQP_REAL Distance;
	friend bool operator > (const ContactF&,const ContactF&);
	friend bool operator < (const ContactF&,const ContactF&);
	
	
};







typedef list<ContactF> ContactF_LIST;

typedef list<ContactF>::iterator ContactF_IT;

typedef list<ContactF>::reverse_iterator REVERSE_ContactF_IT;

struct PQP_TOCResult
{
 // stats
 
 int num_bv_tests;
 int num_tri_tests;
 int num_frontnodes;
 double query_time_secs;
 
 ContactF_LIST cont_l;
 // xform from model 1 to model 2
 
 PQP_REAL R[3][3];
 PQP_REAL T[3];
 
 PQP_REAL rel_err;
 PQP_REAL abs_err;
 PQP_REAL UpboundTOC;
 PQP_REAL ini_dis; 
 PQP_REAL distance;
 PQP_REAL Treshold_Contact;
 PQP_REAL  ori_distance;
 PQP_REAL Tridistance;
 PQP_REAL p1[3];
 PQP_REAL p2[3];
 int qsize;
 PQP_REAL mint; //the minimum time of BV
 
 PQP_REAL BV1Radius;//the radius of BV1
 PQP_REAL BV2Radius;//the radius of BV2
 
 //toc and collision flag
 PQP_REAL toc;//toc >= 1.0 if collision free, else toc < 1.0
 bool collisionfree;//collision free flag
 int numCA; //the number of CA performed
 PQP_REAL R1_toc[3][3], T1_toc[3];
 PQP_REAL R2_toc[3][3], T2_toc[3];

 int num_contact;
 int CurrentLength;

public:
 
 int NumBVTests() { return num_bv_tests; }
 int NumTriTests() { return num_tri_tests; }
 double QueryTimeSecs() { return query_time_secs; }
 
 // The following distance and points established the minimum distance
 // for the models, within the relative and absolute error bounds
 // specified.
 // Points are defined: PQP_REAL p1[3], p2[3];
 
 PQP_REAL Distance() { return distance; }
 const PQP_REAL *P1() { return p1; }
 const PQP_REAL *P2() { return p2; }
};
#endif




















































































