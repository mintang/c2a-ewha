#include "C2A/interpmotion.h"
#include "MatVec.h"

#include <iostream>
#include <string>

#include "gl/glut.h"
#include "TriDist.h"


void tmp_PQP_2_SWIFT_vec( const PQP_REAL T[3], Real T_new[3])
{
	int i;
	for(i = 0; i < 3; i++)
	{
		T_new[i] = T[i];
	}
}
void tmp_PQP_2_SWIFT_mat(const PQP_REAL R[3][3], Real R_new[12])
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

inline 
void Mat_V(PQP_REAL i[3],PQP_REAL j[3],PQP_REAL k[3],const PQP_REAL R[3][3])
{
	for (int l=0;l<3;l++)		
	{
		i[l] = R[0][l];
		j[l] = R[1][l];
		k[l] = R[2][l];
	}
	

}
CScrewMotion::CScrewMotion(
				   const PQP_REAL R0[3][3], 
				   const PQP_REAL T0[3],						 
				   const PQP_REAL R1[3][3], 
				   const PQP_REAL T1[3])
{

	//change rotation matrix to quaternion
	PQP_REAL q0[7],q1[7];

	Matrix3x3 R0_temp,R1_temp;
	R0_temp.Set_Value(R0);
	R1_temp.Set_Value(R1);


	Quaternion Q0,Q1;

	Q0 = R0_temp.Quaternion_();
	Q1 = R1_temp.Quaternion_();

	//R0_temp.Set_Value(Q0);
	//R1_temp.Set_Value(Q1);


	Q0.Get_Value(q0);
	Q1.Get_Value(q1);

	for (int i=0;i<3;i++)
	{
		q0[i+4] = T0[i];
		q1[i+4] = T1[i];
	}

	CScrewMotion screw(q0, q1); 
	b = screw.b;
	d = screw.d;
	VcV(p,screw.p);
	VcV(S,screw.S);


}
CScrewMotion::CScrewMotion()
{
	b=0;
	d=0;
	pure_translation = false;
  
}
CScrewMotion::CScrewMotion(const PQP_REAL q0[7], 
						   const PQP_REAL q1[7])
{
	Quaternion s_q0(-q0[0], -q0[1], -q0[2], q0[3]);

	Quaternion s_q1(q1[0], q1[1], q1[2], q1[3]);

	Quaternion s_temp=s_q0 % s_q1;


	double s    = 1<s_temp.W()?1:s_temp.W();
	double sign = s < 0 ? -1:1;
	double a    = (fabs(s - 1) <= 1e-40 ||
		fabs(s + 1) <= 1e-40) ? (2 * sign) : (sign * acos(2 * s * s - 1) / sqrt(1 - s * s));

	double tagent[3];
	tagent[0] = a * s_temp.X(); 
	tagent[1] = a * s_temp.Y(); 
	tagent[2] = a * s_temp.Z(); 

	VcV(S,tagent);
	b=Vlength(S);

	Vnormalize(S);
	PQP_REAL Ttemp[3],temp[3];


	VmV(Ttemp,q1+4,q0+4);

	if (b<1e-30)
	{
		pure_translation = true;
		b=0;
		d=Vlength(Ttemp);
		VcV(S,Ttemp);
		VcV(p,q0+4);		
	}
	else
	{
		pure_translation = false;

		VpV(temp,q0+4,q1+4);

		PQP_REAL STtemp[3],t[3];
		VcrossV(STtemp,S,Ttemp);
		VxS(t,STtemp,1/tan(b/2));

		VpV(STtemp,t,temp);
		VxS(p,STtemp,double(1/2.0));

		d=VdotV(Ttemp,S);
	}

}

CInterpMotion::CInterpMotion(GMP_INTERP_MODE itpMode, 
							 const PQP_REAL R0[3][3], const PQP_REAL T0[3],
							 const PQP_REAL R1[3][3], const PQP_REAL T1[3])
{

	Real A_R1_r[12], A_T1_r[3];
	Real A_R2_r[12], A_T2_r[3];

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
void CInterpMotion::velocity_screw(const PQP_REAL R0[3][3],const PQP_REAL T0[3],const PQP_REAL R1[3][3],const PQP_REAL T1[3])
{
}

void CInterpMotion::velocity_screw(const PQP_REAL q0[7],const PQP_REAL q1[7])
{
}
CInterpMotion::~CInterpMotion()
{
}

void CInterpMotion::SetSrcDst(const PQP_REAL q0[7], const PQP_REAL q1[7])
{
	Quaternion s_q0(q0[1], q0[2], q0[3], q0[0]);
	Coord3D t0(q0[4], q0[5], q0[6]);

	quaternion_s.Set_Value(q0[1], q0[2], q0[3], q0[0]);

	transform.Set_Rotation(s_q0);
	transform.Set_Translation(t0);

	transform_s = transform;

	Quaternion s_q1(q1[1], q1[2], q1[3], q1[0]);

    quaternion_t.Set_Value(q1[1], q1[2], q1[3], q1[0]);


	Coord3D t1(q1[4], q1[5], q1[6]);
	transform_t.Set_Rotation(s_q1);
	transform_t.Set_Translation(t1);
	m_transform_t_inv.Invert(transform_t);
}


void CInterpMotion::LinearAngularVelocity(Coord3D &axis, PQP_REAL& angVel)
{
#ifdef USE_QUATERNION_DIFF
	Quaternion s_q0(-quaternion_s.X(), -quaternion_s.Y(), -quaternion_s.Z(), quaternion_s.W());

	Quaternion s_q1(quaternion_t.X(), quaternion_t.Y(), quaternion_t.Z(), quaternion_t.W());

	
#else
	quaternion_s = transform_s.Rotation().Quaternion_(); 
	quaternion_t = transform_t.Rotation().Quaternion_(); 
	Quaternion s_q0, s_q1;
	s_q0.Set_Value(-quaternion_s.X(), -quaternion_s.Y(), -quaternion_s.Z(), quaternion_s.W());
	s_q1 =  quaternion_t;

#endif

	Quaternion s_temp=s_q0 % s_q1;


	double s    = 1<s_temp.W()?1:s_temp.W();
	double sign = s < 0 ? -1:1;
	double a    = (fabs(s - 1) <= 1e-40 ||
		fabs(s + 1) <= 1e-40) ? (2 * sign) : (sign * acos(2 * s * s - 1) / sqrt(1 - s * s));

	double tagent[3];
	tagent[0] = a * s_temp.X(); 
	tagent[1] = a * s_temp.Y(); 
	tagent[2] = a * s_temp.Z(); 

	axis.Set_Value(tagent);
	angVel=Vlength(tagent);

	PQP_REAL len = axis.Length_Sq();
	if (len < 0.00000001f)
	{
		axis.Set_Value(1.f,0.f,0.f);
	}
	else
		axis /= sqrt(len);


}


Quaternion CInterpMotion::DeltaRt(Real t)
{
	Real   ang  = 0.5f*m_angVel *t;
	Coord3D axis = m_axis *( sin(ang));

	return Quaternion(axis.X(),axis.Y(),axis.Z(),cos(ang));
}

Quaternion CInterpMotion::AbsoluteRt(Real t)
{
	Quaternion d_rt = DeltaRt(t);
	Quaternion orn0 = transform_s.Quaternion_();

	return orn0% d_rt ;
}

#define GMP_CCD_SECURITY_DISTANCE_RATIO 0.1

Real CInterpMotion::ConservD(Real d)
{
	Real retD;

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



double CInterpMotion::computeTOC(PQP_REAL d, PQP_REAL r1, PQP_REAL S[3])
{
	PQP_REAL time=0;
	return time;
}

double CInterpMotion::computeTOC_MotionBound(PQP_REAL T[3],PQP_REAL d, C2A_BV *V,
								  PQP_REAL N[3])
{
	PQP_REAL path_max=0;
	return path_max;
}


bool  CInterpMotion::CAonRSS(PQP_REAL r1[3][3], PQP_REAL tt1[3],
									PQP_REAL r2[3][3],PQP_REAL tt2[3], 
									PQP_REAL R[3][3],PQP_REAL T[3], 
									C2A_BV *b1, C2A_BV *b2, 
									PQP_REAL *mint,  PQP_REAL *distance)
{
	bool flag=true;
	return flag;

}

bool  CInterpMotion::CAonNonAdjacentTriangles(PQP_REAL R[3][3], PQP_REAL T[3], 
											  PQP_REAL triA[3][3],PQP_REAL triB[3][3],
											  PQP_REAL *mint, PQP_REAL *distance)
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
    Coord3D v3;
#ifdef SCREW_MOTION
	v3=cv_screw;
#else
	v3=cv;
#endif

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

void CInterpMotion_Linear::velocity_screw(const PQP_REAL q0[7],const PQP_REAL q1[7])
{


	CScrewMotion screw(q0,q1);
	m_axis.Set_Value(screw.S);
	m_angVel=screw.b;


	point_screw.Set_Value(screw.p);

	VcV(m_screw.p,screw.p);
	VcV(m_screw.S,screw.S);
	m_screw.d = screw.d;
	m_screw.b = screw.b;
	m_screw.pure_translation = screw.pure_translation;

	PQP_REAL temp1[3],temp2[3],temp3[3];
	if (screw.pure_translation)
	{
		VxS(temp2,screw.S,screw.d);
		cv.Set_Value(temp2);
		cv_screw.Set_Value(temp2);
	}
	else
	{
		VcrossV(temp1,screw.S,screw.p);

		VxS(temp3,temp1,screw.b);
		VxS(temp2,screw.S,screw.d);
		VmV(temp1,temp2,temp3);
		cv.Set_Value(temp2);



		Coord3D temp4,cv_temp;
		cv_temp.Set_Value(temp1);
		temp4=transform_s.Translation();
		temp4.Get_Value(temp2);

		VcrossV(temp3,m_screw.S,temp2);

		VxS(temp2,temp3,m_screw.b);

		temp4.Set_Value(temp2);
		cv_screw = cv_temp + temp4;

	}




}

void CInterpMotion_Linear::velocity_screw(const PQP_REAL R0[3][3],const PQP_REAL T0[3],const PQP_REAL R1[3][3],const PQP_REAL T1[3])
{


	CScrewMotion screw(R0,T0,R1,T1);
    m_axis.Set_Value(screw.S);
	m_angVel=screw.b;

	PQP_REAL temp1[3],temp2[3],temp3[3];
	if (screw.pure_translation)
	{
		VxS(temp2,screw.S,screw.d);
		cv.Set_Value(temp2);
	}
	else
	{
		VcrossV(temp1,screw.S,screw.p);

		VxS(temp3,temp1,screw.b);
		VxS(temp2,screw.S,screw.d);
		VmV(temp1,temp2,temp3);
		cv.Set_Value(temp2);



 		Coord3D temp4,cv_temp;
		cv_temp.Set_Value(temp1);
		temp4=transform_s.Translation();
		temp4.Get_Value(temp2);

		VcrossV(temp3,m_screw.S,temp2);

		VxS(temp2,temp3,m_screw.b);

		temp4.Set_Value(temp2);
		cv_screw = cv_temp + temp4;
	
	}

	point_screw.Set_Value(screw.p);

	VcV(m_screw.p,screw.p);
	VcV(m_screw.S,screw.S);
	m_screw.d = screw.d;
	m_screw.b = screw.b;
	m_screw.pure_translation = screw.pure_translation;


}

void CInterpMotion_Linear::velocity()
{

	cv = (transform_t.Translation() - transform_s.Translation());
	LinearAngularVelocity(m_axis, m_angVel);
}


CInterpMotion_Linear::CInterpMotion_Linear(const PQP_REAL R0[3][3], const PQP_REAL T0[3],
										   const PQP_REAL R1[3][3], const PQP_REAL T1[3]) : CInterpMotion(GMP_IM_LINEAR, R0, T0, R1, T1)
{
#ifdef SCREW_MOTION
	velocity_screw(R0,T0,R1,T1);
#else
	velocity();
#endif
}

CInterpMotion_Linear::CInterpMotion_Linear(const PQP_REAL q0[7], const PQP_REAL q1[7]) : CInterpMotion(GMP_IM_LINEAR, q0, q1)
{
#ifdef SCREW_MOTION

	velocity_screw(q0,q1);

#else
	velocity();
#endif
}


bool CInterpMotion_Linear::integrate(const double dt_input, PQP_REAL qua[7])
{
	double dt = dt_input;

	if(dt > 1)
		dt = 1;

#ifdef SCREW_MOTION
	

	Matrix3x3 temp;
	temp.Set_Value_exp(m_screw.S,m_screw.b*dt);
	Coord3D temp1;

	temp1=transform_s.Translation()-point_screw;

	PQP_REAL temp_val[3][3],temp1_val[3],temp2[3],temp3[3];
	temp.Get_Value(temp_val);
	temp1.Get_Value(temp1_val);

	MxV(temp2,temp_val,temp1_val);


	VpVxS(temp3,temp2,m_screw.S,dt*m_screw.d);

	temp1.Set_Value(temp3);

	transform.Set_Translation(temp1  + point_screw );




#else
	transform.Set_Translation(transform_s.Translation() + cv* dt);

#endif
	
	Quaternion predictedOrn = AbsoluteRt(dt);
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
									C2A_BV *b1, C2A_BV *b2, 
									PQP_REAL *mint,  PQP_REAL *distance)
{
	PQP_REAL c1_clo[3], c2_clo[3], temp1[3], result1[3], result2[3], temp2[3];

	bool bValid_C_clo;
	PQP_REAL S[3];

	double d = C2A_BV_Distance(R, T, b1, b2, c1_clo, c2_clo, bValid_C_clo,S);

	MxV(temp1, b1->R_loc, S);
	MxV(S, r1, temp1);
	 
	PQP_REAL cvc[3], Vel[3];
	cv.Get_Value(cvc);

	MTxV(temp1, r1, cvc);
	MTxV(Vel, b1->R_loc, temp1);

	PQP_REAL tocf = computeTOC(d, b1->angularRadius,S);

	PQP_REAL  total_toc=tocf;
	int nIters=1;

	while ((d>=m_toc_delta)&&(total_toc<=mint[0])&&(nIters<50))
	{
		VpVxS(temp2,T,Vel,-total_toc);			

		d=C2A_BV_Distance(R, temp2, b1, b2, c1_clo, c2_clo, bValid_C_clo,S);
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

		tocf = computeTOC(d, b1->angularRadius,S);

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
													 PQP_REAL *mint, PQP_REAL *distance)
{
	PQP_REAL p[3],q[3],triA_t[3][3];

	PQP_REAL Vel[3];
	PQP_REAL cvc[3];
	cv.Get_Value(cvc);
	MTxV(Vel, R,cvc);

	PQP_REAL d=TriDist(p,q,triA,triB);

	PQP_REAL n[3],total_toc;
	VmV(n,q,p);
	Vnormalize(n);		

	PQP_REAL u= (VdotV(Vel, n));

	if (u<=0)
	{
		u=1e-30;
	}
		
	if (d==0)//&& u>SMALL_NUM
	{
		mint[0]=0.0;
		distance[0]=0.0;
		return true;
	}

		PQP_REAL dt=ConservD(d)/u;//
		int nIters=1;

		if (dt>=mint[0])
		{
			
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

		if (total_toc<=mint[0]&&total_toc>=0.0)//
		{
			mint[0]=total_toc;
			distance[0]=(total_toc)*Vlength(cvc)+m_toc_delta;
			return true;

		}

		return false;

}


double CInterpMotion_Linear::computeTOC(PQP_REAL d, PQP_REAL r1, PQP_REAL S[3])
{
	
	PQP_REAL cwc[3], cross[3];
    PQP_REAL v_max,w_max;

	Vnormalize(S);


	if (m_angVel==0)
	{
		w_max=0;
		PQP_REAL cvc[3];

#ifdef SCREW_MOTION

		PQP_REAL cc[3],temp1[3],temp2[3],temp;
		VxS(cc,m_screw.S,m_screw.d);
		Coord3D temp4=transform_s.Translation();
		temp4.Get_Value(temp1);
		VmV(temp1,temp1,m_screw.p);
		VcrossV(temp2,temp1,m_screw.S);
		v_max = (VdotV(cc, S));
		temp=r1+Vlength(temp2);
		w_max = temp * Vlength(cross);

#else
		cv.Get_Value(cvc);
#endif	

		v_max= VdotV(cvc, S);
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
		VcrossV(cross, cwc, S);

		//compute the angular motion
		w_max =  r1 * Vlength(cross);	

		//compute the translation motion
		PQP_REAL cvc[3];

#ifdef SCREW_MOTION
		PQP_REAL cc[3],temp1[3],temp2[3],temp;
		VxS(cc,m_screw.S,m_screw.d);
		Coord3D temp4=transform_s.Translation();
		temp4.Get_Value(temp1);
		VmV(temp1,temp1,m_screw.p); 
		VcrossV(temp2,temp1,m_screw.S);
		v_max = (VdotV(cc, S));
		temp=r1+Vlength(temp2);
		w_max = temp * Vlength(cross);
#else
		cv.Get_Value(cvc);
		v_max= VdotV(cvc, S);
#endif	
		
		if (v_max<0)
		{
			v_max=0;
		}

		//angular motion + translation motion


	}
	PQP_REAL path_max= w_max + v_max;

	if(path_max==0)
	{
		path_max=1e-30;
	}

	return path_max;
}


double CInterpMotion_Linear::computeTOC_MotionBound(PQP_REAL T[3], PQP_REAL d, C2A_BV *V, PQP_REAL N[3])
{	
	PQP_REAL  cwc[3], cross[3], w_max;
//	PQP_REAL  C1[3],C2[3],C3[3],C4[3],a1,a2,a3,a4,a;

	//set the direction
	Vnormalize(N);
	m_axis.Get_Value(cwc);
	VcrossV(cross, cwc, N);

	//VcrossV(C1,V->Corner[0],cwc);
	//VcrossV(C2,V->Corner[1],cwc);
	//VcrossV(C3,V->Corner[2],cwc);
	//VcrossV(C4,V->Corner[3],cwc);	
		
	//a1=Vlength(C1);
	//a2=Vlength(C2);
	//a3=Vlength(C3);
	//a4=Vlength(C4);
		
	//a= a1;
	//if (a2 > a) a = a2;
	//if (a3 > a) a = a3;
	//if (a4 > a) a = a4;   
	

     PQP_REAL temp;
	//temp=V->r+a;
//	if (temp>V->angularRadius)
	{
		temp=V->angularRadius;		

	}

    w_max = (temp)* Vlength(cross)*m_angVel;
  
	PQP_REAL cvc[3], v_max;

#ifdef SCREW_MOTION

	PQP_REAL cc[3],temp1[3],temp2[3],temp5;
	VxS(cc,m_screw.S,m_screw.d);
	Coord3D temp4=transform_s.Translation();
	temp4.Get_Value(temp1);
	VmV(temp1,temp1,m_screw.p);
	VcrossV(temp2,temp1,m_screw.S);
	v_max = VdotV(cc, N);
	temp5 = V->angularRadius + Vlength(temp2);
	w_max = temp5 * Vlength(cross);

#else
	cv.Get_Value(cvc);
	v_max = VdotV(cvc, N);
#endif	

	
	if (v_max<0)
	{
		v_max=0;
	}
	
	//angular motion + translation motion
	PQP_REAL path_max= v_max+w_max;
	if (path_max<=0)
	{
		path_max=1e-30;
	}
	
	return path_max;

}




















