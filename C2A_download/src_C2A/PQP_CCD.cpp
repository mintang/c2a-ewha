#include "PQP.h"
#include <stdio.h>
#include <string.h>
#include "PQP_CCD.h"

#include "BVTQ.h"
#include "Build.h"
#include "MatVec.h"
#include "GetTime.h"
#include "TriDist.h"
#include "ListTri.h"
#include "InterpMotion.h"

#define   Max_Value 1e+30

/*the definition just for translation motion mode*/
//#define TranslationCCD    
#define  C2A_METHOD true




PQP_REAL tolerance_Cd=1;

// CCD algorithm based on conservative advancement for polygonal soup models 

// toc
double TOCStepRecurse_Dis(PQP_REAL tolerance_d, CInterpMotion *objmotion1, CInterpMotion *objmotion2,PQP_TOCResult *res, double R[3][3], double T[3], 
        PQP_Model *o1, int b1,
        PQP_Model *o2, int b2,
        double delta)
{
	
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();
  

  PQP_REAL r1[3][3], r2[3][3], tt1[3], tt2[3], result1[3], result2[3], temp1[3], temp2[3],cwc[3];

  (objmotion1->transform).Rotation().Get_Value(r1);
  (objmotion2->transform).Rotation().Get_Value(r2);
 
  (objmotion1->transform).Translation().Get_Value(tt1);
  (objmotion2->transform).Translation().Get_Value(tt2);
  
  objmotion1->m_axis.Get_Value(cwc);

  if (l1 && l2)
  { 
    double p[3], q[3];

    PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];

    PQP_REAL dTri = TriDistance(res->R,res->T,t1,t2,p,q);//TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,bCollided) ;
	
	if (dTri <= res->distance)
	{
		res->distance=dTri;
		double S[3];
		MxV(temp1, r1, p);
		VpV(result1, temp1, tt1);
		MxV(temp2, r1, q);
		VpV(result2, temp2, tt1);
		VmV(S,result2,result1);
		VcV(res->p1,p);
		VcV(res->p2,q); 

		//PQP_REAL mint = objmotion1->computeTOC(dTri, o1->child(b1)->angularRadius, S); //  

		PQP_REAL motionbound1 = objmotion1->compute_MotionBound(dTri, o1->child(b1)->angularRadius, S); // 
		PQP_REAL motionbound2 = objmotion2->compute_MotionBound(dTri, o2->child(b2)->angularRadius, S); // 

		PQP_REAL mint = (dTri-0.01*tolerance_d)/(motionbound1+motionbound2);//  

		if (mint <= res->mint) 
			res->mint = mint;

		res->distance=dTri;	
		o1->last_tri = t1;
		o2->last_tri = t2;

	
	}
	  res->num_tri_tests++;
	  res->num_frontnodes++;
	  res->ori_distance=dTri;

      return dTri;  
  }

  // First, perform distance tests on the children. Then traverse
  // them recursively, but test the closer pair first, the further
  // pair second.

  int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
  double R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

  double sz1 = o1->child(b1)->GetSize();
  double sz2 = o2->child(b2)->GetSize();

  if (l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1

    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child+1;
    c2 = b2;
   
    MTxM(R1,o1->child(a1)->R,R);
    VmV(Ttemp,T,o1->child(a1)->Tr);
    MTxV(T1,o1->child(a1)->R,Ttemp);

    MTxM(R2,o1->child(c1)->R,R);
    VmV(Ttemp,T,o1->child(c1)->Tr);
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
    MxVpV(T1,R,o2->child(a2)->Tr,T);

    MxM(R2,R,o2->child(c2)->R);
    MxVpV(T2,R,o2->child(c2)->Tr,T);
  }

  PQP_REAL a1_clo[3], a2_clo[3];
  PQP_REAL c1_clo[3], c2_clo[3];
  bool bValid_A_clo;
  bool bValid_C_clo;
  PQP_REAL S[3];// direction of distance between two minimum distance points

  double d2 = BV_Distance(R2, T2, o1->child(c1), o2->child(c2), c1_clo, c2_clo, bValid_C_clo,S);

  
   if(!bValid_C_clo)
   {
  	  //compute the normal  
  	 MxV(temp1, r1, o1->child(c1)->com);
  	 VpV(result1, temp1, tt1);
  	 MxV(temp2, r2, o2->child(c2)->com);
  	 VpV(result2, temp2, tt2);  
  	 VmV(S,result2,result1);
  }
  else
  {
  	 MxV(temp1, o1->child(c1)->R_loc, S);
  	 MxV(S, r1, temp1);
  } 
   
  PQP_REAL minta,mintb;
  PQP_REAL motionbound1,motionbound2;
  
 // minta = objmotion1->computeTOC(d2, o1->child(c1)->angularRadius,S);//

  motionbound1 = objmotion1->computeTOC_MotionBound(T2,R2,d2,o1->child(c1),S);
  motionbound2 = objmotion2->computeTOC_MotionBound(T2,R2,d2,o2->child(c2),S);
  minta = fabs(d2 - 0.01* tolerance_d)/(motionbound1+motionbound2);

  /////////////////////////////////////////////////////////////////////////////////////////////

  double d1 = BV_Distance(R1, T1, o1->child(a1), o2->child(a2), a1_clo, a2_clo, bValid_A_clo,S);
  
  if(!bValid_A_clo)
  {
  	 		  
  	 MxV(temp1, r1, o1->child(a1)->com);
  	 VpV(result1, temp1, tt1);
  	 MxV(temp2, r2, o2->child(a2)->com);
  	 VpV(result2, temp2, tt2);  
  	 VmV(S,result2,result1);
  }
  else
  {
  	 MxV(temp1, o1->child(a1)->R_loc, S);
  	 MxV(S, r1, temp1);
   
  }	 

 // mintb = objmotion1->computeTOC(d1, o1->child(a1)->angularRadius,S); // 
  motionbound1 = objmotion1->computeTOC_MotionBound(T1,R1,d1,o1->child(a1),S);
  motionbound2 = objmotion2->computeTOC_MotionBound(T1,R1,d1,o2->child(a2),S);
  mintb = fabs(d1 - 0.01* tolerance_d)/(motionbound1+motionbound2);

  res->num_bv_tests += 2;
 
  if (d2 < d1)
  {	
  
	  if (minta<res->UpboundTOC && ((d2 < (res->distance - res->abs_err)) || 
		  (d2*(1 + res->rel_err) < res->distance))) 
	  {  
		  TOCStepRecurse_Dis(tolerance_d,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);  
	
	  }
	  else
	  {
		  res->num_frontnodes++;

		  if (minta < res->mint)	
		  {
			  res->mint=minta;
		  }		 
	  }
	 
	  if (mintb<res->UpboundTOC &&((d1 < (res->distance - res->abs_err)) || 
	 	 (d1*(1 + res->rel_err) < res->distance))) 
	  {  
		TOCStepRecurse_Dis(tolerance_d,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
	  }
	  else
	  {
		  res->num_frontnodes++;

		  if (mintb < res->mint)
		  {     
	 		  res->mint=mintb;
		  }

	  }
	
  }
  else 
  {
	  if (mintb<res->UpboundTOC &&((d1 < (res->distance - res->abs_err)) || 
		  (d1*(1 + res->rel_err) < res->distance))) 
	  { 
    			  
		  TOCStepRecurse_Dis(tolerance_d,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
	  }
	  else
	  {	
		  res->num_frontnodes++;
	
		  if (mintb < res->mint)
		  {     
			  res->mint=mintb;
		  }
	  }
	  
	
	  if (minta<res->UpboundTOC &&((d2 < (res->distance - res->abs_err)) || 
		  (d2*(1 + res->rel_err) < res->distance))) 
	  {	 
	  
		 
		  TOCStepRecurse_Dis(tolerance_d,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);
	  }
	  else
	  {
		  res->num_frontnodes++;
		 
		  if (minta < res->mint)
		  {     
			  res->mint=minta;
		  }
		  
	  }

  
  }
  return PQP_OK;
  
}


//*********************************************************************************************************************//
//  CA_Translation                                                                                                     //
//*********************************************************************************************************************//



double TOCStepRecurse_Dis_Translation(PQP_REAL tolerance_t, CInterpMotion *objmotion1, CInterpMotion *objmotion2,PQP_TOCResult *res, double R[3][3], double T[3], 
						  PQP_Model *o1, int b1,
						  PQP_Model *o2, int b2,
						  double delta)
{

	int l1 = o1->child(b1)->Leaf();
	int l2 = o2->child(b2)->Leaf();


	PQP_REAL r1[3][3], r2[3][3], tt1[3], tt2[3],cwc[3];

	(objmotion1->transform).Rotation().Get_Value(r1);
	(objmotion2->transform).Rotation().Get_Value(r2);

	(objmotion1->transform).Translation().Get_Value(tt1);
	(objmotion2->transform).Translation().Get_Value(tt2);

	objmotion1->m_axis.Get_Value(cwc);

	if (l1 && l2)
	{ 
		PQP_REAL tri1[3][3], tri2[3][3];

		PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
		PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];
		
		VcV(tri1[0], t1->p1);
		VcV(tri1[1], t1->p2);
		VcV(tri1[2], t1->p3);
		MxVpV(tri2[0], res->R, t2->p1, res->T);
		MxVpV(tri2[1], res->R, t2->p2, res->T);
		MxVpV(tri2[2], res->R, t2->p3, res->T);
		PQP_REAL dist=1e+30;
		
		int Inter=1;
        bool CAtri_flag;


		CAtri_flag=objmotion1->CAonNonAdjacentTriangles(r1, tt1, 
			tri1,tri2,&res->mint, &dist,&Inter);


		res->numCA+=Inter;
		
		
		if (dist<res->distance)
		{	
		
			res->distance=dist;
			o1->last_tri = t1;
			o2->last_tri = t2;

		}
		res->num_tri_tests++;		

		return res->distance;  
	}


	int a1,a2,c1,c2;  
	double R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

	double sz1 = o1->child(b1)->GetSize();
	double sz2 = o2->child(b2)->GetSize();

	if (l2 || (!l1 && (sz1 > sz2)))
	{
		// visit the children of b1

		a1 = o1->child(b1)->first_child;
		a2 = b2;
		c1 = o1->child(b1)->first_child+1;
		c2 = b2;

		MTxM(R1,o1->child(a1)->R,R);
		VmV(Ttemp,T,o1->child(a1)->Tr);
		MTxV(T1,o1->child(a1)->R,Ttemp);

		MTxM(R2,o1->child(c1)->R,R);
		VmV(Ttemp,T,o1->child(c1)->Tr);
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
		MxVpV(T1,R,o2->child(a2)->Tr,T);

		MxM(R2,R,o2->child(c2)->R);
		MxVpV(T2,R,o2->child(c2)->Tr,T);
	}



	
	double d1, d2;
	PQP_REAL minta,mintc;


    d1=d2=Max_Value;
	minta=mintc=res->mint;
	bool CARSS_flag_a=objmotion1->CAonRSS(r1, tt1, r2, tt2, R1, T1, o1->child(a1), o2->child(a2),  &minta, &d1);

	bool CARSS_flag_c=objmotion1->CAonRSS(r1, tt1, r2, tt2, R2, T2, o1->child(c1), o2->child(c2),  &mintc, &d2);

	res->num_bv_tests+=2;

	if (d2 < d1)
	{	

		if (mintc<res->mint && ((d2 < (res->distance - res->abs_err)) || 
			(d2*(1 + res->rel_err) < res->distance))) 
		{  
			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);  

		}
		



		if (minta<res->mint &&((d1 < (res->distance - res->abs_err)) || 
			(d1*(1 + res->rel_err) < res->distance))) 
		{  
			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
		}
		

	}
	else 
	{
		if (minta<res->mint &&((d1 < (res->distance - res->abs_err)) || 
			(d1*(1 + res->rel_err) < res->distance))) 
		{ 

			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
		}
		

		if (mintc<res->mint &&((d2 < (res->distance - res->abs_err)) || 
			(d2*(1 + res->rel_err) < res->distance))) 
		{	 


			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);
		}		


	}
	return PQP_OK;

}



void
TOCStepRecurse_Dis_contact(PQP_TOCResult *res,
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


    PQP_REAL p[3], q[3];

    PQP_Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    PQP_Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];
	ContactFeature f1,f2;
	bool bCollided=false;
    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,bCollided) ;
  
    if (d <= res->distance) 
    {
	
			int FeatureType_A = f1.type+1;//1-vertex;2-Edge;3-Face;
			int FeatureType_B = f2.type+1;//1-vertex;2-Edge;3-Face;
			int FeatureID_A[3];
			int FeatureID_B[3];
			
			if (f1.type==0)
			{
				FeatureID_A[0]=t1->Index[f1.fid[0]];
				
			}
			if(f1.type==1)
			{
				FeatureID_A[0]=t1->Index[f1.fid[0]];
				FeatureID_A[1]=t1->Index[(f1.fid[0]+1)%3];		
				
			}
			if(f1.type==2)
			{
				FeatureID_A[0]=t1->Index[0];
				FeatureID_A[1]=t1->Index[1];	
				FeatureID_A[2]=t1->Index[2];		
				
			}
			
			
			if(f2.type==0)
			{
				FeatureID_B[0]=t2->Index[f2.fid[0]];
				
			}
			if(f2.type==1)
			{
				FeatureID_B[0]=t2->Index[f2.fid[0]];
				FeatureID_B[1]=t2->Index[(f2.fid[0]+1)%3];		
				
			}
			if(f2.type==2)
			{
				FeatureID_B[0]=t2->Index[0];
				FeatureID_B[1]=t2->Index[1];	
				FeatureID_B[2]=t2->Index[2];		
				
			}
	
			
			int TriangleID_A = -o1->child(b1)->first_child - 1;
			int TriangleID_B = -o2->child(b2)->first_child - 1;
			
			PQP_REAL P_A[3];
			VcV(P_A,p);
			
			PQP_REAL P_B[3];
			PQP_REAL temp[3];
			VmV(temp,q,res->T);
			MTxV(P_B,res->R,temp);
			
			
			
			
			PQP_REAL Distance=d;
			
			res->cont_l.push_front(ContactF(FeatureType_A, FeatureType_B, FeatureID_A, FeatureID_B,
				TriangleID_A, TriangleID_B,  P_A,  P_B, Distance));

			if(res->Tridistance>d)
			{
				res->Tridistance=d;

			}
			
			res->num_contact++;

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

    VmV(Ttemp,T,o1->child(a1)->Tr);


    MTxV(T1,o1->child(a1)->R,Ttemp);

    MTxM(R2,o1->child(c1)->R,R);

    VmV(Ttemp,T,o1->child(c1)->Tr);


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

    MxVpV(T1,R,o2->child(a2)->Tr,T);



    MxM(R2,R,o2->child(c2)->R);

    MxVpV(T2,R,o2->child(c2)->Tr,T);

    

  }

  PQP_REAL S[3];
  

  PQP_REAL d1 = BV_Distance(R1, T1, o1->child(a1), o2->child(a2),S);
  PQP_REAL d2 = BV_Distance(R2, T2, o1->child(c1), o2->child(c2),S);

  if (d2 < d1)
  {
    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R2, T2, o1, c1, o2, c2);      
    }

    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R1, T1, o1, a1, o2, a2);
    }
  }
  else 
  {
    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R1, T1, o1, a1, o2, a2);
    }

    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R2, T2, o1, c1, o2, c2);      
    }
  }
}

int 
PQP_TOCStep_Contact(CInterpMotion *objmotion1, CInterpMotion *objmotion2, PQP_TOCResult *res,
			PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
			PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
			double threshold)
{


	MTxM(res->R,R1,R2);
	double Ttemp[3];
	VmV(Ttemp, T2, T1);  
	MTxV(res->T, R1, Ttemp);

	
	double Rtemp[3][3], R[3][3], T[3];
	
	MxM(Rtemp,res->R,o2->child(0)->R);
	MTxM(R,o1->child(0)->R,Rtemp);
	
	MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->Tr);
	MTxV(T,o1->child(0)->R,Ttemp);

	PQP_REAL dTri;
	dTri=res->distance;

	res->distance = threshold;
	res->Tridistance = res->distance;
		
    res->rel_err=0;
	res->abs_err=0;
	
	   
	res->mint = 1;

    TOCStepRecurse_Dis_contact(res, R,T,o1,0,o2,0);

	res->distance = res->Tridistance;//res->Tridistance;
	

	return PQP_OK;


	
}

// every CA iterative
int 
PQP_TOCStep(CInterpMotion *objmotion1, CInterpMotion *objmotion2, PQP_TOCResult *res,
			PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
			PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
			PQP_REAL tolerance_d,PQP_REAL tolerance_t,
			int qsize)
{
    double Ttemp[3];
   
	MTxM(res->R,R1,R2);
	
	VmV(Ttemp, T2, T1);  
	MTxV(res->T, R1, Ttemp);

	
	double Rtemp[3][3], R[3][3], T[3];
	
	MxM(Rtemp,res->R,o2->child(0)->R);
	MTxM(R,o1->child(0)->R,Rtemp);
	
	MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->Tr);
	MTxV(T,o1->child(0)->R,Ttemp);


	PQP_REAL p[3],q[3];


	PQP_Tri *t1, *t2;
	t1=o1->last_tri;
	t2=o2->last_tri;

#ifdef TranslationCCD
	PQP_REAL r1[3][3],tt1[3];
	(objmotion1->transform).Rotation().Get_Value(r1);	

	(objmotion1->transform).Translation().Get_Value(tt1);


    PQP_REAL dTri,mint;

	PQP_REAL tri1[3][3], tri2[3][3];
	int Inter;
	

	VcV(tri1[0], t1->p1);
	VcV(tri1[1], t1->p2);
	VcV(tri1[2], t1->p3);
	MxVpV(tri2[0], res->R, t2->p1, res->T);
	MxVpV(tri2[1], res->R, t2->p2, res->T);
	MxVpV(tri2[2], res->R, t2->p3, res->T);
	mint=1.0;
	dTri=Max_Value;
	bool CAtri_flag=objmotion1->CAonNonAdjacentTriangles(r1, tt1, 
		tri1,tri2,&mint, &dTri,&Inter);
	if (CAtri_flag)
	{
		res->mint=mint;
		res->distance=dTri;
	}
	else
	{
		res->mint=1.0;
		res->distance=Max_Value;
	   
	}


#else


	PQP_REAL dTri = res->distance = TriDistance(res->R,res->T,t1,t2,p,q) ;

	if (res->numCA==0)
	{ 	
		res->mint=1;
	}
	if (!C2A_METHOD)
	{
		res->abs_err=0;
		res->rel_err=0;

	}
	else
	{
		//traverse the BVH Tree to compute the toc
		//controlled

		if (res->mint<=0.001||res->distance<=0.1||res->numCA>4)
		{ 
			res->abs_err=0;
			res->rel_err=0;

		}
		else
		{

			res->abs_err=Max_Value;

			if(res->numCA<=3)
			{   
				res->rel_err=2;

			}
			else
			{
				res->rel_err=0.5;
			}   	

		}


	}

	
	   
	res->mint = 1; 
#endif
    
#ifdef TranslationCCD 

	res->abs_err=0;
	res->rel_err=0;

	res->numCA = 0;

    TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R,T,o1,0,o2,0,qsize);
    if (res->num_tri_tests>0)
    {
		res->numCA= int(res->numCA/res->num_tri_tests);
    }
	else
	{
	    res->numCA=0;
	}
	
	t1=o1->last_tri;
	t2=o2->last_tri;
	objmotion1->integrate(res->mint, R1, T1);
	objmotion2->integrate(res->mint, R2, T2); 

	MTxM(res->R,R1,R2);

	VmV(Ttemp, T2, T1);  
	MTxV(res->T, R1, Ttemp);

	res->distance=TriDistance(res->R,res->T,t1,t2,p,q) ;

#else
	TOCStepRecurse_Dis(tolerance_d,objmotion1, objmotion2, res, R,T,o1,0,o2,0,qsize);
#endif
	
	return PQP_OK;
	
}

//get the contact features
PQP_REAL PQP_QueryContact(CInterpMotion *objmotion1, CInterpMotion *objmotion2,
								PQP_TOCResult *res,  PQP_Model *o1, PQP_Model *o2,
								double threshold)
{

	PQP_REAL toc = 0;

	res->num_contact=0;
	res->UpboundTOC=1;
 
	PQP_REAL R1[3][3];
	PQP_REAL T1[3];
	PQP_REAL R2[3][3];
	PQP_REAL T2[3];

	objmotion1->transform.Rotation().Get_Value(R1);
 	objmotion1->transform.Translation().Get_Value(T1);
	objmotion2->transform.Rotation().Get_Value(R2);
	objmotion2->transform.Translation().Get_Value(T2);	
	
	PQP_TOCStep_Contact(objmotion1, objmotion2, res, R1,T1,o1,
			R2,T2,o2,threshold);

	return PQP_OK;
		
		
}

//calculate the toc 
PQP_REAL PQP_QueryTimeOfContact(CInterpMotion *objmotion1, CInterpMotion *objmotion2,
								PQP_TOCResult *res,  PQP_Model *o1, PQP_Model *o2,
								PQP_REAL tolerance_d, PQP_REAL tolerance_t, int qsize)
{
 
	PQP_REAL mint;
	PQP_REAL toc = 0;

	res->num_bv_tests = 0;
	res->num_tri_tests = 0;
	res->num_frontnodes =0;
	res->num_contact=0;
	res->UpboundTOC=1;
    res->numCA=0;
	PQP_REAL R1[3][3];
	PQP_REAL T1[3];
	PQP_REAL R2[3][3];
	PQP_REAL T2[3];

	objmotion1->transform.Rotation().Get_Value(R1);
 	objmotion1->transform.Translation().Get_Value(T1);
	objmotion2->transform.Rotation().Get_Value(R2);
	objmotion2->transform.Translation().Get_Value(T2);

    PQP_REAL lamda = 0.0, lastLamda=0, dlamda=0.0;
	PQP_REAL dist=0;
	int nItrs = 0;

		
	PQP_TOCStep(objmotion1, objmotion2, res, R1,T1,o1,
		R2,T2,o2,tolerance_d,tolerance_t,qsize);

	nItrs=0;

	dist = res->distance; 
	mint = res->mint;
#ifdef TranslationCCD
	res->toc=mint;
	if (res->toc>=1.0)
	{
		res->collisionfree = true; 
	}
	else
	{
		res->collisionfree = false;
		
	}
	toc=mint;

	objmotion1->integrate(toc, res->R1_toc, res->T1_toc);
	objmotion2->integrate(toc, res->R2_toc, res->T2_toc);


	return toc;

#endif
		
	res->numCA = 1;
	if(dist  <= tolerance_d  || mint <= tolerance_t)
	{

		
	
		if(toc>=1)
		{
			res->toc=1.0; 
		    res->collisionfree = true; 
		}
		else
		{
			res->collisionfree = false; 

			res->toc = 0;  

			objmotion1->integrate(res->toc, res->R1_toc, res->T1_toc);
			objmotion2->integrate(res->toc, res->R2_toc, res->T2_toc);
		   
		}
	
	}
		
	while (dist > tolerance_d) //not close enough
	{
		nItrs++;

		if(nItrs > 150)
		{
			break;
		}
			// 1. If dist > path_max, no collision during t= [0, 1] 
			// 2. If dist < path_max, look for TOC
		if( mint >= 1.0)
		{
			res->collisionfree = true;
			res->toc = 0;
			toc = 0;

			return 0;
				
		}
		else
		{
			dlamda = mint;
			if (dlamda < tolerance_t)
			{												
				break;
			}	
			lamda += dlamda;   				
			if (lamda >= 1.0)
			{
				res->collisionfree = true;
				res->toc = 0;
				toc = 0;
				return 0;

			}
			lastLamda = lamda;				
			res->numCA++;   
			objmotion1->integrate(lamda, R1, T1);
			objmotion2->integrate(lamda, R2, T2); 
			res->UpboundTOC=1.0-lamda;				

			PQP_TOCStep(objmotion1, objmotion2, res, R1,T1,o1,
					R2,T2,o2,tolerance_d,tolerance_t,qsize);

			dist=res->distance;
			mint = res->mint; 

		}
	}

	if(dist == 0 && mint<0)
	{
		toc = lastLamda- dlamda ;
	}
	else
	{
		toc = lastLamda;
	}

	res->collisionfree = false;
			
	if(toc>=1-tolerance_t)
		{toc=0;}
	res->toc = toc;
	//compute the result R, T at TOC 
	objmotion1->integrate(toc, res->R1_toc, res->T1_toc);//make res->R_toc equals R1
	objmotion2->integrate(toc, res->R2_toc, res->T2_toc);	
		
	return toc;
		
		
}

//calculate the toc 
PQP_REAL PQP_QueryTimeOfContact_TranslateCA(CInterpMotion *objmotion1, CInterpMotion *objmotion2,
								PQP_TOCResult *res,  PQP_Model *o1, PQP_Model *o2,
								PQP_REAL tolerance_d, PQP_REAL tolerance_t, int qsize)
{


	int tmpi;
	int tmpj;

	PQP_REAL mint;
	PQP_REAL toc = 0;

	res->num_bv_tests = 0;
	res->num_tri_tests = 0;
	res->num_contact=0;
	res->UpboundTOC=1;
	res->numCA=0;
	PQP_REAL R1[3][3];
	PQP_REAL T1[3];
	PQP_REAL R2[3][3];
	PQP_REAL T2[3];

	objmotion1->transform.Rotation().Get_Value(R1);
	objmotion1->transform.Translation().Get_Value(T1);
	objmotion2->transform.Rotation().Get_Value(R2);
	objmotion2->transform.Translation().Get_Value(T2);
	PQP_REAL lamda = 0.0, lastLamda=0, dlamda=0.0;
	PQP_REAL dist=0;
	int nItrs = 0;
	bool T_CA=false;



	PQP_TOCStep(objmotion1, objmotion2, res, R1,T1,o1,
		R2,T2,o2,tolerance_d,qsize);

	nItrs=0;

	dist = res->distance; 
	mint = res->mint;

	res->numCA = 1;


	if(dist  <= tolerance_d  || mint <= tolerance_t)
	{
		res->collisionfree = false; 			
		res->toc = tolerance_t;
		for(tmpi=0; tmpi<3; tmpi++)
		{	 
			for( tmpj=0; tmpj<3; tmpj++)
			{
				res->R1_toc[tmpi][tmpj] = R1[tmpi][tmpj];
				res->R2_toc[tmpi][tmpj] = R2[tmpi][tmpj];
			}
			res->T1_toc[tmpi] = T1[tmpi];
			res->T2_toc[tmpi] = T2[tmpi];
		}
		if(toc>=1)
			toc=0;	

	}


	while (dist > tolerance_d) //not close enough
	{
		nItrs++;

		if(nItrs > 150)
		{
			break;
		}
		// 1. If dist > path_max, no collision during t= [0, 1] 
		// 2. If dist < path_max, look for TOC
		if( mint >= 1.0)
		{
			res->collisionfree = true;
			res->toc = 0;
			toc = 0;

			return 0;

		}
		else
		{
			dlamda = mint;
			if (dlamda < tolerance_t)
			{												
				break;
			}	
			lamda += dlamda;   				
			if (lamda >= 1.0)
			{
				res->collisionfree = true;
				res->toc = 0;
				toc = 0;
				return 0;

			}
			lastLamda = lamda;				
			res->numCA++;   
			objmotion1->integrate(lamda, R1, T1);
			objmotion2->integrate(lamda, R2, T2); 
			res->UpboundTOC=1.0-lamda;				

			PQP_TOCStep(objmotion1, objmotion2, res, R1,T1,o1,
				R2,T2,o2,tolerance_d,tolerance_t,qsize);

			dist=res->distance;
			mint = res->mint; 

		}
	}

	if(dist == 0 && mint<0)
	{
		toc = lastLamda- dlamda ;
	}
	else
	{
		toc = lastLamda;
	}

	res->collisionfree = false;

	if(toc>=1-tolerance_t)
	{toc=0;}
	res->toc = toc;
	//compute the result R, T at TOC 
	objmotion1->integrate(toc, res->R1_toc, res->T1_toc);//make res->R_toc equals R1


	return toc;


}





void PQP_Model::ComputeCenterOfMass()
{
    int i;
    PQP_REAL area_x2;
    PQP_REAL total_area;
	
    com[0] = 0.0;
    com[1] = 0.0;
    com[2] = 0.0;

    total_area = 0.0;
    for( i = 0; i < num_tris; i++ ) 
	{
		PQP_REAL E0[3], E1[3], S[3];
		VmV(E0, tris[i].p1, tris[i].p2);
		VmV(E1, tris[i].p1, tris[i].p3);
		VcrossV(S, E0, E1);
		area_x2=sqrt(S[0]*S[0] + S[1]*S[1] + S[2]*S[2]);
        total_area += area_x2;
		
        com[0] += (tris[i].p1[0] + tris[i].p2[0] + tris[i].p3[0]) * area_x2;
        com[1] += (tris[i].p1[1] + tris[i].p2[1] + tris[i].p3[1]) * area_x2;
		com[2] += (tris[i].p1[2] + tris[i].p2[2] + tris[i].p3[2]) * area_x2;
    }
	com[0] = com[0] / (3.0 * total_area);
	com[1] = com[1] / (3.0 * total_area);
	com[2] = com[2] / (3.0 * total_area);
	//    com /= 3.0 * total_area;
}

void PQP_Model::ComputeRadius()
{
    int i;
	radius=0.0;
	PQP_REAL d;
    for( i = 0; i < num_tris; i++ )
	{
		PQP_REAL p[3], pr[3];

		p[0] = tris[i].p1[0];
		p[1] = tris[i].p1[1];
		p[2] = tris[i].p1[2];		
		VmV(pr, p, com);
		d = Vlength(pr);		
		if( d > radius ) radius = d;

		p[0] = tris[i].p2[0];
		p[1] = tris[i].p2[1];
		p[2] = tris[i].p2[2];		
		VmV(pr, p, com);
		d = Vlength(pr);	
		if( d > radius ) radius = d;
		
		p[0] = tris[i].p3[0];
		p[1] = tris[i].p3[1];
		p[2] = tris[i].p3[2];		
		VmV(pr, p, com);
		d = Vlength(pr);	
		if( d > radius ) radius = d;
		
    }
    radius;
}












































































































































































































































































































