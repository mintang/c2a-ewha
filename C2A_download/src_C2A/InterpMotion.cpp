#include "interpmotion.h"
#include "MatVec.h"
#include "gmpParam.h"
#include <iostream>
#include <string>
using namespace std;
#include "gl/glut.h"
#include "TriDist.h"


void tmp_PQP_2_SWIFT_vec( const PQP_REAL T[3], SWIFT_Real T_new[3])
{
	int i;
	for(i = 0; i < 3; i++)
	{
		T_new[i] = T[i];
	}
}
void tmp_PQP_2_SWIFT_mat(const PQP_REAL R[3][3], SWIFT_Real R_new[12])
{
	int i,j, k;
	k = 0;
	for( i = 0; i < 3; i++)
	{
		for( j = 0; j < 3; j++)
		{
			R_new[k] = 	R[i][j];	
			k++;
		}
	}
}

CInterpMotion::CInterpMotion(GMP_INTERP_MODE itpMode, 
							 const PQP_REAL R0[3][3], const PQP_REAL T0[3],
							 const PQP_REAL R1[3][3], const PQP_REAL T1[3])
{

	SWIFT_Real A_R1_r[12], A_T1_r[3];
	SWIFT_Real A_R2_r[12], A_T2_r[3];

	tmp_PQP_2_SWIFT_mat(R0, A_R1_r);
	tmp_PQP_2_SWIFT_mat(R1, A_R2_r);

	tmp_PQP_2_SWIFT_vec(T0, A_T1_r);
	tmp_PQP_2_SWIFT_vec(T1, A_T2_r);

	m_itpMode = itpMode;
	transform.Set_Value(A_R1_r, A_T1_r);
	transform_s.Set_Value(A_R1_r, A_T1_r);
	transform_t.Set_Value(A_R2_r, A_T2_r);
	m_transform_t_inv.Invert(transform_t);

}

CInterpMotion::CInterpMotion()
{
}
double CInterpMotion::computeTOC_relate(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3])
{
	PQP_REAL time=0.0;

	return time;
}


CInterpMotion::CInterpMotion(GMP_INTERP_MODE itpMode, 
							 const PQP_REAL q0[7], const PQP_REAL q1[7])
{
	m_itpMode = itpMode;
	SetSrcDst(q0, q1);
	velocity();
}

double CInterpMotion::computeTOC(PQP_REAL d, PQP_REAL r1)
{
 PQP_REAL time=0;
 
 return time;

}

void CInterpMotion::velocity(void)
{
}
void CInterpMotion::velocity_relate(PQP_REAL T_relate[3])
{


}
CInterpMotion::~CInterpMotion()
{
}

void CInterpMotion::SetSrcDst(const PQP_REAL q0[7], const PQP_REAL q1[7])
{
	SWIFT_Quaternion s_q0(q0[1], q0[2], q0[3], q0[0]);
	SWIFT_Triple t0(q0[4], q0[5], q0[6]);

	transform.Set_Rotation(s_q0);
	transform.Set_Translation(t0);

	transform_s = transform;

	SWIFT_Quaternion s_q1(q1[1], q1[2], q1[3], q1[0]);
	SWIFT_Triple t1(q1[4], q1[5], q1[6]);
	transform_t.Set_Rotation(s_q1);
	transform_t.Set_Translation(t1);
	m_transform_t_inv.Invert(transform_t);
}


void CInterpMotion::LinearAngularVel(SWIFT_Triple &axis, SWIFT_Real& angVel)
{
#ifdef USE_QUATERNION_DIFF
	SWIFT_Quaternion orn0	= m1.Quaternion();
	SWIFT_Quaternion orn1a	= m2.Quaternion();
	SWIFT_Quaternion orn1	= orn0.Farthest(orn1a);
	SWIFT_Quaternion dorn	= orn1 % orn0.Inverse();
#else
	SWIFT_Matrix33 dmat = transform_t.Rotation() * ( transform_s.Rotation().Inverse());
	SWIFT_Quaternion dorn= dmat.Quaternion();
#endif


	angVel = dorn.Angle();


	axis.Set_Value(dorn[0],dorn[1],dorn[2]);


	SWIFT_Real len = axis.Length_Sq();
	if (len < 0.00000001f)
	{


		axis = SWIFT_Triple(1.f,0.f,0.f);
	}
	else
		axis /= sqrt(len);

}


SWIFT_Quaternion CInterpMotion::DeltaRt(SWIFT_Real t)
{
	SWIFT_Real   ang  = 0.5f*m_angVel *t;
	SWIFT_Triple axis = m_axis *( sin(ang));

	return SWIFT_Quaternion(axis.X(),axis.Y(),axis.Z(),cos(ang));
}

SWIFT_Quaternion CInterpMotion::AbsoluteRt(SWIFT_Real t)
{
	SWIFT_Quaternion d_rt = DeltaRt(t);
	SWIFT_Quaternion orn0 = transform_s.Quaternion();

	return d_rt % orn0;
}



SWIFT_Real CInterpMotion::ConservD(SWIFT_Real d)
{
	SWIFT_Real retD;

	retD = d - GMP_CCD_SECURITY_DISTANCE_RATIO * m_toc_delta;

	return retD;
}


bool CInterpMotion::integrate(const double dt_input, PQP_REAL R[3][3], PQP_REAL T[3])
{
	PQP_REAL qua[7];
	integrate(dt_input, qua);

	transform.Rotation().Get_Value(R);
	transform.Translation().Get_Value(T);

	return true;
}

double CInterpMotion::computeTOC(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3])
{
	PQP_REAL time=0;
	return time;
}
double CInterpMotion::compute_MotionBound(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3])
{
	PQP_REAL bound=0.0;
	

	return bound;
}



double CInterpMotion::computeTOC_MotionBound(PQP_REAL T[3],PQP_REAL R[3][3],PQP_REAL d, BV *V,
								  PQP_REAL S[3])
{
	PQP_REAL time=0;
	return time;
}
double CInterpMotion::computeTOC_MotionBoundTri(PQP_REAL a,PQP_REAL T[3],PQP_REAL R[3][3], PQP_Model *o1,PQP_REAL d, BV *V,
											 PQP_REAL p[3], PQP_REAL q[3])
{
	PQP_REAL time=0;
	return time;
}


bool  CInterpMotion::CAonRSS(PQP_REAL r1[3][3], PQP_REAL tt1[3],
									PQP_REAL r2[3][3],PQP_REAL tt2[3], 
									PQP_REAL R[3][3],PQP_REAL T[3], 
									BV *b1, BV *b2, 
									PQP_REAL *mint,  PQP_REAL *distance)
{
	bool flag=true;
	return flag;

}

bool  CInterpMotion::CAonNonAdjacentTriangles(PQP_REAL R[3][3], PQP_REAL T[3], 
											  PQP_REAL triA[3][3],PQP_REAL triB[3][3],
											  PQP_REAL *mint, PQP_REAL *distance, int* Inter)
{
	bool flag=true;
	return flag;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Linear Motion

//compute TOC without motion bound
double CInterpMotion_Linear::computeTOC(PQP_REAL d, PQP_REAL r1)
{
	PQP_REAL time;

	PQP_REAL w_max =  r1 * fabs(m_angVel); //m_angVel;  the angular velocity around that axis: counter clockwise
					                       // it might be negative, which means the clockwise rotation

	SWIFT_Triple v3 = cv;//cv linear translation  velocity
	PQP_REAL v_max = v3.Length();
	PQP_REAL path_max= w_max + v_max;

	time = ConservD(d) / path_max;// ConservD = d - security_distance_ratio * m_toc_delta 

	return time;
}

CInterpMotion_Linear::~CInterpMotion_Linear()
{
 
}

CInterpMotion_Linear::CInterpMotion_Linear()
{
}

void CInterpMotion_Linear::velocity_relate(PQP_REAL T_relate[3])
{
	cv_relate.Set_Value(T_relate[0],T_relate[1],T_relate[2]);
  
}
void CInterpMotion_Linear::velocity()
{

	cv = (transform_t.Translation() - transform_s.Translation());
	LinearAngularVel(m_axis, m_angVel);
}


CInterpMotion_Linear::CInterpMotion_Linear(const PQP_REAL R0[3][3], const PQP_REAL T0[3],
										   const PQP_REAL R1[3][3], const PQP_REAL T1[3]) : CInterpMotion(GMP_IM_LINEAR, R0, T0, R1, T1)
{
	velocity();
}

CInterpMotion_Linear::CInterpMotion_Linear(const PQP_REAL q0[7], const PQP_REAL q1[7]) : CInterpMotion(GMP_IM_LINEAR, q0, q1)
{
	velocity();
}


bool CInterpMotion_Linear::integrate(const double dt_input, PQP_REAL qua[7])
{
	double dt = dt_input;

	if(dt > 1)
		dt = 1;

	transform.Set_Translation(transform_s.Translation() + cv* dt);

	SWIFT_Quaternion predictedOrn = AbsoluteRt(dt);
	transform.Set_Rotation(predictedOrn);

	// Rotation
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Translation
	qua[4] = transform.Translation().X();
	qua[5] = transform.Translation().Y();
	qua[6] = transform.Translation().Z();

	return true;
}

bool  CInterpMotion_Linear::CAonRSS(PQP_REAL r1[3][3], PQP_REAL tt1[3],
									PQP_REAL r2[3][3],PQP_REAL tt2[3], 
									PQP_REAL R[3][3],PQP_REAL T[3], 
									BV *b1, BV *b2, 
									PQP_REAL *mint,  PQP_REAL *distance)
{
	PQP_REAL c1_clo[3], c2_clo[3], temp1[3], result1[3], result2[3], temp2[3];

	bool bValid_C_clo;
	PQP_REAL S[3];

	double d = BV_Distance(R, T, b1, b2, c1_clo, c2_clo, bValid_C_clo,S);


	if(!bValid_C_clo)
	{
		//compute the normal  
		MxV(temp1, r1, b1->com);
		VpV(result1, temp1, tt1);
		MxV(temp2, r2, b2->com);
		VpV(result2, temp2, tt2);  
		VmV(S,result2,result1);
	}
	else
	{
		MxV(temp1, b1->R_loc, S);
		MxV(S, r1, temp1);
	} 
	PQP_REAL cvc[3], Vel[3];
	cv_relate.Get_Value(cvc);

	MTxV(temp1, r1, cvc);
	MTxV(Vel, b1->R_loc, temp1);

	PQP_REAL tocf = computeTOC_relate(d, b1->angularRadius,S);

	PQP_REAL  total_toc=tocf;
	int nIters=1;

	while ((d>=m_toc_delta)&&(total_toc<=mint[0])&&(nIters<50))
	{
		VpVxS(temp2,T,Vel,-total_toc);			

		d=BV_Distance(R, temp2, b1, b2, c1_clo, c2_clo, bValid_C_clo,S);
		if (d==0)
		{
			break;
		}

		if(!bValid_C_clo)
		{
			//compute the normal  
	    	MxV(temp1, r1, b1->com);
			VpV(temp2, temp1, tt1);
			VpVxS(result1, temp2, cvc,total_toc);

			MxV(temp2, r2, b2->com);
			VpV(result2, temp2, tt2);  

			VmV(S,result2,result1);



		}
		else
		{
			MxV(temp1, b1->R_loc, S);
			MxV(S, r1, temp1);
		} 

		tocf = computeTOC_relate(d, b1->angularRadius,S);

		nIters++;

		if(tocf<m_toc_delta) 
			break;

		total_toc+=tocf;

	}

	if (total_toc<1.0)
	{
		mint[0]=total_toc;
		distance[0]=(total_toc-m_toc_delta)*Vlength(cvc);
		if (distance[0]<0)
		{
           distance[0]=0;

		}
		return true;

	}

	return false;


}

bool  CInterpMotion_Linear::CAonNonAdjacentTriangles(PQP_REAL R[3][3], PQP_REAL T[3], 
													 PQP_REAL triA[3][3],PQP_REAL triB[3][3],
													 PQP_REAL *mint, PQP_REAL *distance, int* Inter)
{
	PQP_REAL p[3],q[3],triA_t[3][3];

	PQP_REAL Vel[3];
	PQP_REAL cvc[3];

	cv_relate.Get_Value(cvc);
	MTxV(Vel, R, cvc);

	PQP_REAL d=TriDist(p,q,triA,triB);

	PQP_REAL n[3],total_toc;
	VmV(n,q,p);
	Vnormalize(n);		
    int nIters=1;
	PQP_REAL u= (VdotV(Vel, n));

	if (u<=0)
	{
		u=1e-30;
	}
		
	if (d==0 )//&& u>SMALL_NUM
	{
		Inter[0]=nIters;
		mint[0]=0.0;
		distance[0]=0.0;
		return true;
	}

		PQP_REAL dt=ConservD(d)/u;//
		

		if (dt>=mint[0])
		{
			Inter[0]=nIters;
			
			return false;

		}
		total_toc=dt;
	

		while ((d>m_toc_delta)&&(total_toc<=mint[0])&&(nIters<50))
		{
			for(int i=0;i<3;i++)
			{
				VpVxS(triA_t[i],triA[i],Vel,total_toc);
				
			}				

			d=TriDist(p,q,triA_t,triB);
			if (d==0.0)
			{
				break;
			}

		
			VmV(n,q,p);

			Vnormalize(n);
			u=(VdotV(Vel, n));
			if (u<=0)
			{
				u=1e-30;
			}
			PQP_REAL tofc=ConservD(d)/u;//
			nIters++;
            
			if(tofc<m_toc_delta) 
				break;
            total_toc+=tofc;
			

		}
		Inter[0]=nIters;

		if (total_toc<=mint[0]&&total_toc>=0.0)//
		{
			mint[0]=total_toc;
			distance[0]=(total_toc)*Vlength(cvc)+m_toc_delta;
			return true;

		}

		return false;

}

double CInterpMotion_Linear::computeTOC_relate(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3])
{
	PQP_REAL time;
	PQP_REAL cwc[3], cross[3];
	PQP_REAL v_max,w_max;
	//set the direction

	Vnormalize(normal);


	if (m_angVel==0)
	{
		w_max=0;
		PQP_REAL cvc[3];
		cv_relate.Get_Value(cvc);
		v_max= fabs(VdotV(cvc, normal));
		if (v_max<0)
		{
			v_max=0;
		}
	}
	else
	{
		m_axis.Get_Value(cwc);
		cwc[0] *= m_angVel;
		cwc[1] *= m_angVel;
		cwc[2] *= m_angVel;
		VcrossV(cross, cwc, normal);

		//compute the angular motion
		w_max =  r1 * Vlength(cross);	

		//compute the translation motion
		PQP_REAL cvc[3];
		cv.Get_Value(cvc);
		v_max= fabs(VdotV(cvc, normal));
		if (v_max<0)
		{
			v_max=0;
		}

		//angular motion + translation motion


	}
	PQP_REAL path_max= w_max + v_max;


	if(path_max<=0)
	{
		path_max=1e-30;
	}

	//return path_max;
	time = ConservD(d) / path_max;
	if (time<=0)
	{
		time = 0;
	}

	return time;
}

double CInterpMotion_Linear::compute_MotionBound(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3])
{
	
	PQP_REAL cwc[3], cross[3];
	PQP_REAL v_max,w_max;
	//set the direction

	Vnormalize(normal);


	if (m_angVel==0)
	{
		w_max=0;
		PQP_REAL cvc[3];
		cv.Get_Value(cvc);
		v_max= fabs(VdotV(cvc, normal));
		if (v_max<0)
		{
			v_max=0;
		}
	}
	else
	{
		m_axis.Get_Value(cwc);
		cwc[0] *= m_angVel;
		cwc[1] *= m_angVel;
		cwc[2] *= m_angVel;
		VcrossV(cross, cwc, normal);

		//compute the angular motion
		w_max =  r1 * Vlength(cross);	

		//compute the translation motion
		PQP_REAL cvc[3];
		cv.Get_Value(cvc);
		v_max= fabs(VdotV(cvc, normal));
		if (v_max<0)
		{
			v_max=0;
		}

		//angular motion + translation motion


	}
	PQP_REAL path_max= w_max + v_max;


	if(path_max<=0)
	{
		path_max=1e-30;
	}

	return path_max;
	
}

double CInterpMotion_Linear::computeTOC(PQP_REAL d, PQP_REAL r1,PQP_REAL normal[3])
{
	PQP_REAL time;
	PQP_REAL cwc[3], cross[3];
    PQP_REAL v_max,w_max;
	//set the direction

	Vnormalize(normal);


	if (m_angVel==0)
	{
		w_max=0;
		PQP_REAL cvc[3];
		cv.Get_Value(cvc);
		v_max= fabs(VdotV(cvc, normal));
		if (v_max<0)
		{
			v_max=0;
		}
	}
	else
	{
		m_axis.Get_Value(cwc);
		cwc[0] *= m_angVel;
		cwc[1] *= m_angVel;
		cwc[2] *= m_angVel;
		VcrossV(cross, cwc, normal);

		//compute the angular motion
		w_max =  r1 * Vlength(cross);	

		//compute the translation motion
		PQP_REAL cvc[3];
		cv.Get_Value(cvc);
		v_max= fabs(VdotV(cvc, normal));
		if (v_max<0)
		{
			v_max=0;
		}

		//angular motion + translation motion


	}
	PQP_REAL path_max= w_max + v_max;


	if(path_max<=0)
	{
		path_max=1e-30;
	}

  
 	time = ConservD(d) / path_max;
	if (time<=0)
	{
		time = 0;
	}
	
	return time;
}


double CInterpMotion_Linear::computeTOC_MotionBound(PQP_REAL T[3],PQP_REAL R[3][3], PQP_REAL d, BV *V, PQP_REAL S[3])
{
	
	PQP_REAL N[3], cwc[3], cross[3] , w_max, C1[3],C2[3],C3[3],C4[3],a1,a2,a3,a4,a;

	//set the direction
	VcV(N,S);
	Vnormalize(N);
	m_axis.Get_Value(cwc);
	VcrossV(cross, cwc, N);
	VcrossV(C1,V->Corner[0],T);
	VcrossV(C2,V->Corner[1],T);
	VcrossV(C3,V->Corner[2],T);
	VcrossV(C4,V->Corner[3],T);	
		
	a1 = Vlength(C1);
	a2 = Vlength(C2);
	a3 = Vlength(C3);
	a4 = Vlength(C4);
		
	a = a1;
	if (a2 > a) a = a2;
	if (a3 > a) a = a3;
	if (a4 > a) a = a4;   
		
	

	PQP_REAL temp=V->r+a;

	if (temp>V->angularRadius)
	{
	  temp=V->angularRadius;

	}
	
	
    w_max = (temp) * Vlength(cross) * m_angVel;

   
	PQP_REAL cvc[3];
	cv.Get_Value(cvc);
	PQP_REAL v_max = fabs(VdotV(cvc, N));
	
	//angular motion + translation motion
	PQP_REAL path_max = v_max + w_max;

	return path_max;



}
double CInterpMotion_Linear::computeTOC_MotionBoundTri(PQP_REAL a, PQP_REAL T[3],PQP_REAL R[3][3], PQP_Model *o1, PQP_REAL d, BV *V, PQP_REAL p[3], PQP_REAL q[3])
{
	PQP_REAL time;
	PQP_REAL N[3], cwc[3], cross[3], w_max;

	//set the direction
	VmV(N, q, p);
	Vnormalize(N);
	m_axis.Get_Value(cwc);

	VcrossV(cross, cwc, N);


    w_max = a* Vlength(cross)*m_angVel;

    	
	PQP_REAL cvc[3];
	cv.Get_Value(cvc);
	PQP_REAL v_max = fabs(VdotV(cvc, N));
	

	PQP_REAL path_max= v_max+w_max;

	time = ConservD(d) / path_max;
	
	return time;

}





















