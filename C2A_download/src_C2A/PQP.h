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

#ifndef PQP_H
#define PQP_H

#include "PQP_Compile.h"   
#include "PQP_Internal.h"                             
                      
//----------------------------------------------------------------------------
//
//  PQP API Return Values
//
//----------------------------------------------------------------------------

const int PQP_OK = 0; 
  // Used by all API routines upon successful completion except
  // constructors and destructors

const int PQP_ERR_MODEL_OUT_OF_MEMORY = -1; 
  // Returned when an API function cannot obtain enough memory to
  // store or process a PQP_Model object.

const int PQP_ERR_OUT_OF_MEMORY = -2;
  // Returned when a PQP query cannot allocate enough storage to
  // compute or hold query information.  In this case, the returned
  // data should not be trusted.

const int PQP_ERR_UNPROCESSED_MODEL = -3;
  // Returned when an unprocessed model is passed to a function which
  // expects only processed models, such as PQP_Collide() or
  // PQP_Distance().

const int PQP_ERR_BUILD_OUT_OF_SEQUENCE = -4;
  // Returned when: 
  //       1. AddTri() is called before BeginModel().  
  //       2. BeginModel() is called immediately after AddTri().  
  // This error code is something like a warning: the invoked
  // operation takes place anyway, and PQP does what makes "most
  // sense", but the returned error code may tip off the client that
  // something out of the ordinary is happenning.

const int PQP_ERR_BUILD_EMPTY_MODEL = -5; 
  // Returned when EndModel() is called on a model to which no
  // triangles have been added.  This is similar in spirit to the
  // OUT_OF_SEQUENCE return code, except that the requested operation
  // has FAILED -- the model remains "unprocessed", and the client may
  // NOT use it in queries.





const int PQP_ALL_CONTACTS = 1;  // find all pairwise intersecting triangles
const int PQP_FIRST_CONTACT = 2; // report first intersecting tri pair found

int 
PQP_Collide(PQP_CollideResult *result,
            PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
            PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
            int flag = PQP_ALL_CONTACTS);


int 
PQP_Collide_Subtree(PQP_CollideResult *result,
            PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1, int bv1, 
            PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2, int bv2,
            int flag = PQP_ALL_CONTACTS);


#if PQP_BV_TYPE & RSS_TYPE  // this is true by default,
                            // and explained in PQP_Compile.h

//----------------------------------------------------------------------------
//
//  PQP_DistanceResult
//
//  This saves and reports results from a distance query.  
//
//----------------------------------------------------------------------------
//
//  struct PQP_DistanceResult - declaration contained in PQP_Internal.h
//  {
//    // statistics
//  
//    int NumBVTests();
//    int NumTriTests();
//    PQP_REAL QueryTimeSecs();
//  
//    // The following distance and points established the minimum distance
//    // for the models, within the relative and absolute error bounds 
//    // specified.
//
//    PQP_REAL Distance();
//    const PQP_REAL *P1();  // pointers to three PQP_REALs
//    const PQP_REAL *P2();  
//  };

//----------------------------------------------------------------------------
//
//  PQP_Distance() - computes the distance between two PQP_Models
//
//
//  Declare a PQP_DistanceResult struct and pass its pointer to collect
//  distance information.
//
//  "rel_err" is the relative error margin from actual distance.
//  "abs_err" is the absolute error margin from actual distance.  The
//  smaller of the two will be satisfied, so set one large to nullify
//  its effect.
//
//  "qsize" is an optional parameter controlling the size of a priority
//  queue used to direct the search for closest points.  A larger queue
//  can help the algorithm discover the minimum with fewer steps, but
//  will increase the cost of each step. It is not beneficial to increase
//  qsize if the application has frame-to-frame coherence, i.e., the
//  pair of models take small steps between each call, since another
//  speedup trick already accelerates this situation with no overhead.
//
//  However, a queue size of 100 to 200 has been seen to save time in a
//  planning application with "non-coherent" placements of models.
//
//----------------------------------------------------------------------------

int 
PQP_Distance(PQP_DistanceResult *result, 
             PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
             PQP_REAL rel_err, PQP_REAL abs_err,
             int qsize = 2);


//----------------------------------------------------------------------------
//
//  PQP_ToleranceResult
//
//  This saves and reports results from a tolerance query.  
//
//----------------------------------------------------------------------------
//
//  struct PQP_ToleranceResult - declaration contained in PQP_Internal.h
//  {
//    // statistics
//  
//    int NumBVTests(); 
//    int NumTriTests();
//    PQP_REAL QueryTimeSecs();
//  
//    // If the models are closer than ( <= ) tolerance, these points 
//    // and distance were what established this.  Otherwise, 
//    // distance and point values are not meaningful.
//  
//    PQP_REAL Distance();
//    const PQP_REAL *P1();
//    const PQP_REAL *P2();
//  
//    // boolean says whether models are closer than tolerance distance
//  
//    int CloserThanTolerance();
//  };

//----------------------------------------------------------------------------
//
// PQP_Tolerance() - checks if distance between PQP_Models is <= tolerance
//
//
// Declare a PQP_ToleranceResult and pass its pointer to collect
// tolerance information.
//
// The algorithm returns whether the true distance is <= or >
// "tolerance".  This routine does not simply compute true distance
// and compare to the tolerance - models can often be shown closer or
// farther than the tolerance more trivially.  In most cases this
// query should run faster than a distance query would on the same
// models and configurations.
// 
// "qsize" again controls the size of a priority queue used for
// searching.  Not setting qsize is the current recommendation, since
// increasing it has only slowed down our applications.
//
//----------------------------------------------------------------------------

int
PQP_Tolerance(PQP_ToleranceResult *res, 
              PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
              PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
              PQP_REAL tolerance,
              int qsize = 2);


//----------------------------------------------------------------------------
//
//  PQP_ContactResult
//struct PQP_ContactResult 
//{
//	// stats
//	
//	int num_bv_tests;
//	int num_tri_tests;
//	double query_time_secs;
//	
//	// xform from model 1 to model 2
//	
//	PQP_REAL R[3][3];
//	PQP_REAL T[3];
//	
//	PQP_REAL rel_err; 
//	PQP_REAL abs_err; 
//	
//	//	PQP_REAL distance;
//	PQP_REAL p1[3]; 
//	PQP_REAL p2[3];
//	int qsize;
//	
//	//contact features
//	ContactPair *cpairs;
//	int num_cpairs_alloced;
//	int num_cpairs;
//	bool closer_than_ctolerance;
//	bool modelcollided; //whether the two models are colliding
//	PQP_REAL ctolerance;      
//	
//	void FreePairsList(); 
//	void SizeTo(int n);
//	void Add(ContactFeature f1, ContactFeature f2, PQP_REAL distance); 
//	
//	// statistics
//	
//	int NumBVTests() { return num_bv_tests; }
//	int NumTriTests() { return num_tri_tests; }
//	double QueryTimeSecs() { return query_time_secs; }
//	
//	// The following distance and points established the minimum distance
//	// for the models, within the relative and absolute error bounds 
//	// specified.
//	// Points are defined: PQP_REAL p1[3], p2[3];
//	PQP_ContactResult();
//	~PQP_ContactResult();
//	bool CloserThanTolerance() { return closer_than_ctolerance; }
//	
//	
//	//	PQP_REAL Distance() { return distance; }
//	const PQP_REAL *P1() { return p1; }
//	const PQP_REAL *P2() { return p2; }
//}
//  This saves and reports all the contact results from a contact query.  
//
//----------------------------------------------------------------------------


// Extension of PQP_1_3
// To report the closest features within some tolerance for polygonal soup models 
// Dec 2007 to Jan 2008
int 
PQP_Contact(PQP_ContactResult *result, 
             PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
             PQP_REAL rel_err, PQP_REAL abs_err, PQP_REAL tolerance,
             int qsize = 2);

// when the flag is enabled, PQP_Contact return a set of represenative pairs
extern bool bPQPContact_Represenative;


PQP_REAL
TriDistance(PQP_REAL R[3][3], PQP_REAL T[3], PQP_Tri *t1, PQP_Tri *t2,
            PQP_REAL p[3], PQP_REAL q[3]);

PQP_REAL
TriDistance(PQP_REAL R[3][3], PQP_REAL T[3], PQP_Tri *t1, PQP_Tri *t2,
            PQP_REAL p[3], PQP_REAL q[3], ContactFeature &f1, ContactFeature &f2, bool &bCollided);


#endif 
#endif