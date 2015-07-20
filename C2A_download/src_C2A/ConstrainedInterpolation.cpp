#include "interpmotion.h"
#include "ConstrainedInterpolation.h"
#include "MatVec.h"
#include "gmpParam.h"
#include <iostream>
#include <string>
using namespace std;
#include "gl/glut.h"

// #define EPS_PLANE_PARRAL (sin(10*pi/180))^2 ~~~ 0.0302
// #define EPS_PLANE_PARRAL 0.0302
SWIFT_Real EPS_PLANE_PARRAL = 0.0302;


////////////////////////////////////////////////////////////////////////////////////////////
// Linear angular velocity
// The shortest distance between the contact feature is linearly interpolated. 
SWIFT_Real eps_normalize = 0.001;
// SWIFT_Real eps_normalize = 0.00000001f;
int gDegenerate = 0;

bool  MyNormalize(SWIFT_Triple &n)
{
	SWIFT_Real len = n.Length_Sq();
	if (len < eps_normalize)
	{
		// GMP_WARNING("The normal is NULL; Degenerate Situations");
		// GMP_ASSERT("The normal is NULL; Degenerate Situations");

		if(len < 0.00000001f)
			n.Set_Value(0, 0, 1); // set an arbitrary axis
		else
			n.Normalize();

		gDegenerate++;
		return false;
	}
	else
	{
		n.Normalize();
		return true;
	}
}

/*!
  Returns the determinant of the 3x3 submatrix specified by the row and
  column indices.
 */

SWIFT_Real MyDeterminant(const SWIFT_Matrix33 & mat,
						 int r1, int r2, int r3,
					     int c1, int c2, int c3) 
{
/* #if COIN_EXTRA_DEBUG
  // Check indices.
  if (r1<0 || r1>3 || r2<0 || r2>3 || r3<0 || r3>3 ||
      c1<0 || c1>3 || c2<0 || c2>3 || c3<0 || c3>3) {
    SoDebugError::post("SbMatrix::det3",
                       "At least one idx out of bounds [0, 3]. ");
  }
  if (r1==r2 || r1==r3 || r2==r3 ||
      c1==c2 || c1==c3 || c2==c3)
    SoDebugError::post("SbMatrix::det3", "Indices should be distinct.");
#endif // COIN_EXTRA_DEBUG
*/

  // More or less directly from "Advanced Engineering Mathematics"
  // (E. Kreyszig), 6th edition.
  SWIFT_Real a11 = mat.Value()[3*r1+c1];
  SWIFT_Real a12 = mat.Value()[3*r1+c2];
  SWIFT_Real a13 = mat.Value()[3*r1+c3];
  SWIFT_Real a21 = mat.Value()[3*r2+c1];
  SWIFT_Real a22 = mat.Value()[3*r2+c2];
  SWIFT_Real a23 = mat.Value()[3*r2+c3];
  SWIFT_Real a31 = mat.Value()[3*r3+c1];
  SWIFT_Real a32 = mat.Value()[3*r3+c2];
  SWIFT_Real a33 = mat.Value()[3*r3+c3];

  SWIFT_Real M11 = a22 * a33 - a32 * a23;
  SWIFT_Real M21 = -(a12 * a33 - a32 * a13);
  SWIFT_Real M31 = a12 * a23 - a22 * a13;

  return (a11 * M11 + a21 * M21 + a31 * M31);
}


// e_rob = {t0, t1 ), e_obs = {s0, s1}
// The normal is from e_obs to e_rob
SWIFT_Real LineLineSignDist(const SWIFT_Triple & t0, const SWIFT_Triple & t1, 
						const SWIFT_Triple & s0, const SWIFT_Triple & s1)
{
	SWIFT_Triple nor = (s1 - s0) % (t1- t0);
	MyNormalize(nor);
	// nor.Normalize();
	return nor * (t0-s0);
}


void
LinePoints(SWIFT_Real VEC[3], 
		  SWIFT_Real X[3], SWIFT_Real Y[3],             // closest points
          const SWIFT_Real P[3], const SWIFT_Real A[3], // seg 1 origin, vector
          const SWIFT_Real Q[3], const SWIFT_Real B[3]) // seg 2 origin, vector
{
  SWIFT_Real T[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;
  SWIFT_Real TMP[3];

  VmV(T,Q,P);
  A_dot_A = VdotV(A,A);
  B_dot_B = VdotV(B,B);
  A_dot_B = VdotV(A,B);
  A_dot_T = VdotV(A,T);
  B_dot_T = VdotV(B,T);

  // t parameterizes ray P,A 
  // u parameterizes ray Q,B 

  SWIFT_Real t,u;

  // compute t for the closest point on ray P,A to
  // ray Q,B

  SWIFT_Real denom = A_dot_A*B_dot_B - A_dot_B*A_dot_B;

  t = (A_dot_T*B_dot_B - B_dot_T*A_dot_B) / denom;

  // clamp result so t is on the segment P,A

  // if ((t < 0) || isnan(t)) t = 0; else if (t > 1) t = 1;

  // find u for point on ray Q,B closest to point at t

  u = (t*A_dot_B - B_dot_T) / B_dot_B;

  // if u is on segment Q,B, t and u correspond to 
  // closest points, otherwise, clamp u, recompute and
  // clamp t 
/*
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
*/
    VpVxS(Y, Q, B, u);

 /*   if ((t <= 0) || isnan(t)) {
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
  */  VpVxS(X, P, A, t);
      VcrossV(VEC, A, B);
      if (VdotV(VEC, T) < 0) {
        VxS(VEC, VEC, -1);
      }
   /* }
  } */
}



void LineLineNearPts(const SWIFT_Triple & A_worl_1, const SWIFT_Triple & A_worl_2, 
				 const SWIFT_Triple & B_1,      const SWIFT_Triple & B_2,  
				 SWIFT_Triple & near_p1,        SWIFT_Triple & near_p2)
{
	SWIFT_Real VEC[3];
	SWIFT_Real X[3], Y[3]; // closest points
	SWIFT_Real P[3], A[3]; // seg 1 origin, vector
	SWIFT_Real Q[3], B[3]; // seg 2 origin, vector

	P[0] = A_worl_1.X(); P[1] = A_worl_1.Y(); P[2] = A_worl_1.Z(); 

	SWIFT_Triple A_dir = A_worl_2 - A_worl_1;
	MyNormalize(A_dir);
	// A_dir.Normalize();
	A[0] = A_dir.X(); A[1] = A_dir.Y(); A[2] = A_dir.Z();

	Q[0] = B_1.X(); Q[1] = B_1.Y(); Q[2] = B_1.Z(); 

	SWIFT_Triple B_dir = B_2 - B_1;
	MyNormalize(B_dir);
	// B_dir.Normalize();
	B[0] = B_dir.X(); B[1] = B_dir.Y(); B[2] = B_dir.Z();

	LinePoints(VEC, X, Y, P, A, Q, B);

	near_p1.Set_Value(X[0], X[1], X[2]);
	near_p2.Set_Value(Y[0], Y[1], Y[2]);
}


CItpMot_OneContact::CItpMot_OneContact(const PQP_REAL q0[7], const PQP_REAL q1[7] ) : CInterpMotion(GMP_IM_LINEAR_ONE_CONTACT, q0, q1)
{
	m_contact_type = GMP_ITP_ONE_CONTACT_V_F;

	m_vf_pt.Set_Value(0, 0, 0);
	m_rob_m_0_0.Set_Value(0, 0, 0);
	m_rob_m_1_1.Set_Value(0, 0, 0);
	m_obs_p0.Set_Value(0, 0, 0);
	m_obs_p1.Set_Value(0, 0, 0);

	velocity();
}

CItpMot_OneContact::CItpMot_OneContact()
{
}


CItpMot_OneContact::~CItpMot_OneContact()
{
	delete m_pMultConstraint;
}

void CItpMot_OneContact::velocity(void)
{
	LinearAngularVel(m_axis, m_angVel);
}

bool CItpMot_OneContact::integrate(const double dt_input, PQP_REAL qua[7])
{
	double t = dt_input;

	if(t > 1)
		t = 1;

	if(m_contact_type == GMP_ITP_ONE_CONTACT_V_F)
	{
		vf_integrate(t, qua);
	}
	else if(m_contact_type == GMP_ITP_ONE_CONTACT_F_V)
	{
		fv_integrate_SPM08(t, qua);
	}
	else if(m_contact_type == GMP_ITP_ONE_CONTACT_E_E)
	{
		ee_integrate(t, qua);
	}

	return true;
}

string CItpMot_OneContact::Print()
{
	char oStr[256] = " ";

	if(m_contact_type == GMP_ITP_ONE_CONTACT_V_F)
		strcat(oStr, "V-F");
	else if(m_contact_type == GMP_ITP_ONE_CONTACT_F_V)
		strcat(oStr, "F-V");
	else if(m_contact_type == GMP_ITP_ONE_CONTACT_E_E)
		strcat(oStr, "E-E");
	else
		strcat(oStr, "Unknown Wrong ");

	return oStr;
}

double fv_sign_0 = 1;
double fv_sign_1 = 1;

// should be changed RSS formulae though the equation here is same
void CItpMot_OneContact::fv_integrate_SPM08(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);

	// Rotation
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Transformation
	SWIFT_Real dist_0, dist_1, dist_t;
	SWIFT_Triple t0;

	// Note: m_0_1 should be computed by: R_1*m_0_1 + T_1 = m_1_1 
	// m_1_1 is known
	SWIFT_Triple m_0_t, m_0_1; 
	SWIFT_Triple tra;

	dist_0 = m_fv_obs_pt.Dist(m_fv_rob_m_0_0); 
	dist_1 = m_fv_obs_pt.Dist(m_fv_rob_m_1_1); 

	// There is still a bug to handle the situation when the sign is negative. 
	dist_0 *= fv_sign_0;

	if(m_fv_sign_0_flip)
		fv_sign_0 = -1;

	if(m_fv_sign_1_flip)
		fv_sign_1 = -1;

	dist_1 *= fv_sign_1;

	m_0_1 = transform_s * (m_transform_t_inv * m_fv_rob_m_1_1);

	dist_t = (1-t)*dist_0 + t*dist_1;
	m_0_t  = (1-t)*m_fv_rob_m_0_0 + t* m_0_1;

	tra = transform_s.Translation() - dist_t * m_fv_rob_nor_0 - m_0_t;
	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!
	tra = d_rt.Rotate(tra);
	tra = tra + m_fv_obs_pt;

	qua[4] = tra.X();
	qua[5] = tra.Y();
	qua[6] = tra.Z();

	transform.Set_Translation(tra);
}

void CItpMot_OneContact::SetupEE(SWIFT_Triple A0, SWIFT_Triple A1,
							     SWIFT_Triple B0, SWIFT_Triple B1)
{
	m_pt_A_0 = A0;
	m_pt_A_1 = A1;
	m_pt_B_0 = B0;
	m_pt_B_1 = B1;

	// Can be precomptued. Be careful DeltaRt is correct
	m_rel_tra = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());

	SWIFT_Triple A_0_new, A_1_new;

	A_0_new = DeltaRt(1).Rotate(m_pt_A_0)+m_rel_tra;
	A_1_new = DeltaRt(1).Rotate(m_pt_A_1)+m_rel_tra;

	m_d_0 = LineLineSignDist(m_pt_A_0, m_pt_A_1, m_pt_B_0, m_pt_B_1); 
	m_d_1 = LineLineSignDist(A_0_new, A_1_new, 
					     m_pt_B_0, m_pt_B_1); 

	LineLineNearPts(m_pt_A_0, m_pt_A_1, m_pt_B_0, m_pt_B_1,  
			m_near_A_0_0,
			m_near_B_0);

	LineLineNearPts(A_0_new, A_1_new, m_pt_B_0, m_pt_B_1,  
			m_near_A_1_1,
			m_near_B_1);

	// Transform back from configuration q1 to q0
	m_near_A_0_1 = DeltaRt(-1).Rotate(m_near_A_1_1 - m_rel_tra);
}

void CItpMot_OneContact::ee_integrate(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Can be precomptued. Be careful DeltaRt is correct
	m_rel_tra = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());

	// The goal is to compute the delta_tra
	// The abs_tra = Delta_Rot * q0.T + delta_tra
	SWIFT_Triple delta_tra;
	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!
	SWIFT_Real d_t;
	SWIFT_Triple near_A_0_t, near_B_t, N_t;

	d_t = (1-t)*m_d_0 + t*m_d_1;
	N_t = (m_pt_B_1 - m_pt_B_0) % d_rt.Rotate(m_pt_A_1-m_pt_A_0);
	MyNormalize(N_t);

	near_B_t = (1-t)*m_near_B_0 + t*m_near_B_1;
	near_A_0_t = (1-t)*m_near_A_0_0 + t*m_near_A_0_1;

	delta_tra = d_t * N_t + near_B_t - d_rt.Rotate(near_A_0_t);

	SWIFT_Triple abs_tra;
	abs_tra = DeltaRt(t).Rotate(transform_s.Translation()) + delta_tra;

	qua[4] = abs_tra.X();
	qua[5] = abs_tra.Y();
	qua[6] = abs_tra.Z();
}


double ee_sign_0 = 1;
double ee_sign_1 = 1;
// bool bOld_Normal_Method = true; // when use the "old normal" method with ee_sign_1 = -1 could get the same result as the new normal method !!! unbelievable 
bool bOld_Normal_Method = false;



// It is buggy; to remove
void CItpMot_OneContact::ee_integrate_SPM08_temp(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);

	// Rotation
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Transformation
	SWIFT_Real dist_0, dist_1, dist_t;
	SWIFT_Triple t0, N0;

	// Note: m_0_1 should be computed by: R_1*m_0_1 + T_1 = m_1_1 
	// m_1_1 is known
	SWIFT_Triple m_0_t, m_0_1; 
	SWIFT_Triple tra;

	dist_0 = m_rob_m_0_0.Dist(m_obs_p0); 
	dist_1 = m_rob_m_1_1.Dist(m_obs_p1); 

	// There is still a bug to handle the situation when the sign is negative. 
	dist_0 *= ee_sign_0;
	dist_1 *= ee_sign_1;

	m_0_1 = transform_s * (m_transform_t_inv * m_rob_m_1_1);
	//N0    = m_0_1 -  m_rob_m_0_0;
	//N0    = N0 % (m_obs_p1 -  m_obs_p0);
	N0    = m_rob_m_0_0 - m_obs_p0;
	N0.Normalize();

	dist_t = (1-t)*dist_0 + t*dist_1;
	m_0_t  = (1-t)*m_rob_m_0_0 + t* m_0_1;

	if(bOld_Normal_Method)
	{
		tra = transform_s.Translation() + dist_t * N0 - m_0_t;
		SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!
		tra = d_rt.Rotate(tra);
	}
	else // new
	{
/*		tra = transform_s.Translation() + dist_t * N0 - m_0_t;
		SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!
		tra = d_rt.Rotate(tra);
*/
		tra = transform_s.Translation() - m_0_t;
		SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!
		tra = d_rt.Rotate(tra);
		// The normal is computed based on cross product
		SWIFT_Triple NT = DeltaRt(t).Rotate(m_rob_m_0_0-m_0_1) % (m_obs_p1 - m_obs_p0);
		NT.Normalize();
		tra = tra + dist_t * NT;
	}







	tra = tra + (1-t)*m_obs_p0+t*m_obs_p1;

	qua[4] = tra.X();
	qua[5] = tra.Y();
	qua[6] = tra.Z();

	transform.Set_Translation(tra);


  // For one contact V/F
  //SWIFT_Triple m_basePt
  //SWIFT_Triple m_tra_dir; // constant translation velocity of the basePt

}

void CItpMot_OneContact::vf_integrate(const double dt, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(dt);
	transform.Set_Rotation(predictedOrn);

	// Rotation
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Transformation
	SWIFT_Triple tra;
	SWIFT_Triple p, p0, p1;

	// Global coordinate
	p = m_vf_pt; // Local coordinate of the point p
	p0 = transform_s * p;
	p1 = transform_t * p;

	tra = (1-dt) * p0 + dt * p1;
	tra = tra -  transform.Rotation() * p; // !!! p not p0

	qua[4] = tra.X();
	qua[5] = tra.Y();
	qua[6] = tra.Z();

	transform.Set_Translation(tra);


  // For one contact V/F
  //SWIFT_Triple m_basePt
  //SWIFT_Triple m_tra_dir; // constant translation velocity of the basePt
}


/////////////////////////////////////////////////////////////
CItpMot_TwoContacts::CItpMot_TwoContacts(const PQP_REAL q0[7], const PQP_REAL q1[7] ) : CInterpMotion(GMP_IM_LINEAR_TWO_CONTACTS, q0, q1)
{
	m_contact_type = GMP_ITP_TWO_CONTACTS_VF_VF;

	//m_vf_pt.Set_Value(0, 0, 0);
	//m_rob_m_0_0.Set_Value(0, 0, 0);
	//m_rob_m_1_1.Set_Value(0, 0, 0);
	//m_obs_p0.Set_Value(0, 0, 0);
	//m_obs_p1.Set_Value(0, 0, 0);

	velocity();
}

CItpMot_TwoContacts::CItpMot_TwoContacts()
{
}


CItpMot_TwoContacts::~CItpMot_TwoContacts()
{
}

void CItpMot_TwoContacts::velocity(void)
{
	LinearAngularVel(m_axis, m_angVel);
}

bool CItpMot_TwoContacts::integrate(const double dt_input, PQP_REAL qua[7])
{
	double t = dt_input;

	if(t > 1)
		t = 1;

	if(m_contact_type == GMP_ITP_TWO_CONTACTS_VF_VF)
	{
		vf_vf_integrate_new(t, qua);
	}
	else if(m_contact_type == GMP_ITP_TWO_CONTACTS_FV_FV)
	{
		fv_fv_integrate(t, qua);
	}
	else if(m_contact_type == GMP_ITP_TWO_CONTACTS_EE_EE)
	{
		ee_ee_integrate(t, qua);
	}
	//else if(m_contact_type == GMP_ITP_ONE_CONTACT_E_E)
	//{
	//	ee_integrate(t, qua);
	//}

	return true;
}

/*
double fv_sign_0 = 1;
double fv_sign_1 = 1;

void CItpMot_OneContact::fv_integrate(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);

	// Rotation
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Transformation
	SWIFT_Real dist_0, dist_1, dist_t;
	SWIFT_Triple t0;

	// Note: m_0_1 should be computed by: R_1*m_0_1 + T_1 = m_1_1 
	// m_1_1 is known
	SWIFT_Triple m_0_t, m_0_1; 
	SWIFT_Triple tra;

	dist_0 = m_fv_obs_pt.Dist(m_fv_rob_m_0_0); 
	dist_1 = m_fv_obs_pt.Dist(m_fv_rob_m_1_1); 

	// There is still a bug to handle the situation when the sign is negative. 
	dist_0 *= fv_sign_0;

	if(m_fv_sign_flip)
		fv_sign_1 = -1;

	dist_1 *= fv_sign_1;

	m_0_1 = transform_s * (m_transform_t_inv * m_fv_rob_m_1_1);

	dist_t = (1-t)*dist_0 + t*dist_1;
	m_0_t  = (1-t)*m_fv_rob_m_0_0 + t* m_0_1;

	tra = transform_s.Translation() - dist_t * m_fv_rob_nor_0 - m_0_t;
	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!
	tra = d_rt.Rotate(tra);
	tra = tra + m_fv_obs_pt;

	qua[4] = tra.X();
	qua[5] = tra.Y();
	qua[6] = tra.Z();

	transform.Set_Translation(tra);
}


double ee_sign_0 = 1;
double ee_sign_1 = 1;

void CItpMot_OneContact::ee_integrate(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);

	// Rotation
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Transformation
	SWIFT_Real dist_0, dist_1, dist_t;
	SWIFT_Triple t0, N0;

	// Note: m_0_1 should be computed by: R_1*m_0_1 + T_1 = m_1_1 
	// m_1_1 is known
	SWIFT_Triple m_0_t, m_0_1; 
	SWIFT_Triple tra;

	dist_0 = m_rob_m_0_0.Dist(m_obs_p0); 
	dist_1 = m_rob_m_1_1.Dist(m_obs_p1); 

	// There is still a bug to handle the situation when the sign is negative. 
	dist_0 *= ee_sign_0;
	dist_1 *= ee_sign_1;

	m_0_1 = transform_s * (m_transform_t_inv * m_rob_m_1_1);
	//N0    = m_0_1 -  m_rob_m_0_0;
	//N0    = N0 % (m_obs_p1 -  m_obs_p0);
	N0    = m_rob_m_0_0 - m_obs_p0;
	N0.Normalize();

	dist_t = (1-t)*dist_0 + t*dist_1;
	m_0_t  = (1-t)*m_rob_m_0_0 + t* m_0_1;

	tra = transform_s.Translation() + dist_t * N0 - m_0_t;
	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!
	tra = d_rt.Rotate(tra);
	tra = tra + (1-t)*m_obs_p0+t*m_obs_p1;

	qua[4] = tra.X();
	qua[5] = tra.Y();
	qua[6] = tra.Z();

	transform.Set_Translation(tra);


  // For one contact V/F
  //SWIFT_Triple m_basePt
  //SWIFT_Triple m_tra_dir; // constant translation velocity of the basePt

}
*/

void CItpMot_TwoContacts::vf_vf_integrate(const double dt, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(dt);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Rotation
	// Constant: Only need to be computed once
	SWIFT_Triple z_dir;
	z_dir = m_vf_vf_f1_nor % m_vf_vf_f2_nor;
	z_dir.Normalize(); // Normalize: to check whether its lenght is less than a threshold: tolerance!!!
	SWIFT_Triple T0, T1, Tt_local, Tt_global;
	SWIFT_Real c;

	T0.Set_Value(0, 0, 0);
	T1 = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());

	c = T1 * z_dir;

	// The distance of the closest feature to the face at q0 and q1
	SWIFT_Real dist_pt1_s, dist_pt1_t;
	SWIFT_Real dist_pt2_s, dist_pt2_t;
	SWIFT_Triple pt_1_t, pt_2_t;

	pt_1_t = DeltaRt(1).Rotate(m_vf_vf_pt_1) + T1;
	pt_2_t = DeltaRt(1).Rotate(m_vf_vf_pt_2) + T1;

	dist_pt1_s = m_vf_vf_f1_nor * m_vf_vf_pt_1 + m_vf_vf_f1_d;
	dist_pt1_t = m_vf_vf_f1_nor * pt_1_t + m_vf_vf_f1_d;

	dist_pt2_s = m_vf_vf_f2_nor * m_vf_vf_pt_2 + m_vf_vf_f2_d;
	dist_pt2_t = m_vf_vf_f2_nor * pt_2_t + m_vf_vf_f2_d;

	// Varing w.r.t. t
	SWIFT_Quaternion d_rt = DeltaRt(dt); // This is the relative amount of rotation !!!
	SWIFT_Real s1, s2; // varing w.r.t t

	s1 = m_vf_vf_f1_nor * d_rt.Rotate(m_vf_vf_pt_1) + m_vf_vf_f1_d - ((1-dt)*dist_pt1_s + dt*dist_pt1_t); 
	s2 = m_vf_vf_f2_nor * d_rt.Rotate(m_vf_vf_pt_2) + m_vf_vf_f2_d - ((1-dt)*dist_pt2_s + dt*dist_pt2_t); 

	SWIFT_Real theta; // constant
	theta = acos(m_vf_vf_f1_nor * m_vf_vf_f2_nor); // Might a bug here



	SWIFT_Matrix33 frame;
	SWIFT_Triple x_axis, y_axis, z_axis;
	x_axis = m_vf_vf_f1_nor;
	z_axis = z_dir;
	y_axis = z_axis % x_axis;
	y_axis.Normalize();

	// The frame is due to the concept of coordinate transformation. Set column
	// Set cols, since you want the pt with local coordinate (1, 0, 0) be x_axis with global coordinate
	frame.Set_Value_Cols(x_axis, y_axis, z_axis);

	Tt_local.Set_X(-s1);
	Tt_local.Set_Y(-s2/sin(theta) - s1 / tan(theta));
	Tt_local.Set_Z(dt * c);

	Tt_global = frame * Tt_local;

	SWIFT_Triple tra;
	SWIFT_Triple tmp_old_tra = d_rt.Rotate(transform_s.Translation());

	tra = d_rt.Rotate(transform_s.Translation()) + Tt_global;
	tra = tmp_old_tra + Tt_global;

	qua[4] = tra.X();
	qua[5] = tra.Y();
	qua[6] = tra.Z();

	transform.Set_Translation(tra);

	vf_vf_integrate_new(dt, qua);


/*
	// Transformation
	SWIFT_Triple tra;
	SWIFT_Triple p, p0, p1;

	// Global coordinate
	p = m_vf_pt; // Local coordinate of the point p
	p0 = transform_s * p;
	p1 = transform_t * p;

	tra = (1-dt) * p0 + dt * p1;
	tra = tra -  transform.Rotation() * p; // !!! p not p0

	qua[4] = tra.X();
	qua[5] = tra.Y();
	qua[6] = tra.Z();

	transform.Set_Translation(tra);
*/

  // For one contact V/F
  //SWIFT_Triple m_basePt
  //SWIFT_Triple m_tra_dir; // constant translation velocity of the basePt
}

void CItpMot_TwoContacts::vf_vf_integrate_new(const double t, PQP_REAL qua[7])
{

  m_pt_0 = m_vf_vf_pt_1;   // The position of the vertex (at q_0 for robot) (for obstacle, ready for using) 
  m_f0_nor = m_vf_vf_f1_nor; // The normal of the face (at q_0 for robot) (for obstacle, ready for using) 
  m_f0_d = m_vf_vf_f1_d;

  m_pt_1 = m_vf_vf_pt_2;
  m_f1_nor = m_vf_vf_f2_nor;
  m_f1_d = m_vf_vf_f2_d;





	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Can be precomptued. Be careful DeltaRt is correct
	m_rel_tra = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());

	// The goal is to compute the delta_tra
	// The abs_tra = Delta_Rot * q0.T + delta_tra
	SWIFT_Triple delta_tra;
	SWIFT_Triple base;
	SWIFT_Matrix33 frame;
	SWIFT_Real s0_t, s1_t, lambda;
	SWIFT_Triple N0_t, N1_t, dirIntersect; // The normal of intersection of the planes 0 and 1

	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!

	N0_t = m_f0_nor;
	N1_t = m_f1_nor;

	// To set the frame
	// TODO: for VF/VF or FV/FV case, the inverse can be precomputed
	dirIntersect = N0_t % N1_t;
	dirIntersect.Normalize();
	// a bug frame.Set_Value_Cols(N0_t, N1_t, dirIntersect); 
	frame.Set_Value_Rows(N0_t, N1_t, dirIntersect); // The frame here is due to the linear system. So set rows

	SWIFT_Real dist_pt0_s, dist_pt0_t; 
	SWIFT_Real dist_pt1_s, dist_pt1_t; 
	// The distance of the closest feature to the face at q0 and q1
	SWIFT_Triple pt_0_t, pt_1_t;
	pt_0_t = DeltaRt(1).Rotate(m_pt_0) + m_rel_tra; // should be the new point at t=1!!!
	pt_1_t = DeltaRt(1).Rotate(m_pt_1) + m_rel_tra;

	dist_pt0_s = m_f0_nor * m_pt_0 + m_f0_d;
	dist_pt0_t = m_f0_nor * pt_0_t + m_f0_d;

	dist_pt1_s = m_f1_nor * m_pt_1 + m_f1_d;
	dist_pt1_t = m_f1_nor * pt_1_t + m_f1_d;

	s0_t = m_f0_nor * d_rt.Rotate(m_pt_0) + m_f0_d - ((1-t)*dist_pt0_s + t*dist_pt0_t);
	s1_t = m_f1_nor * d_rt.Rotate(m_pt_1) + m_f1_d - ((1-t)*dist_pt1_s + t*dist_pt1_t);

	lambda = m_rel_tra * dirIntersect;

	base = frame.Inverse() * SWIFT_Triple(-1*s0_t, -1*s1_t, 0);
	delta_tra = base + dirIntersect * t * lambda;

	SWIFT_Triple abs_tra;
	SWIFT_Triple tmp_tra = d_rt.Rotate(transform_s.Translation());

	// abs_tra = DeltaRt(t).Rotate(transform_s.Translation()) + delta_tra;
	// abs_tra = d_rt.Rotate(transform_s.Translation()) + delta_tra;
	abs_tra = tmp_tra + delta_tra;

	qua[4] = abs_tra.X();
	qua[5] = abs_tra.Y();
	qua[6] = abs_tra.Z();

	transform.Set_Translation(abs_tra);
}

SWIFT_Real CItpMot_TwoContacts::tmp_lambda()
{
	PQP_REAL qua[7];
	SWIFT_Real t = 1;
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Can be precomptued. Be careful DeltaRt is correct
	m_rel_tra = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());

	// The goal is to compute the delta_tra
	// The abs_tra = Delta_Rot * q0.T + delta_tra
	SWIFT_Triple delta_tra;
	SWIFT_Triple base;
	SWIFT_Matrix33 frame;
	SWIFT_Real s0_t, s1_t, lambda;
	SWIFT_Triple N0_t, N1_t, dirIntersect; // The normal of the fi

	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!

	N0_t = -1*d_rt.Rotate(m_f0_nor);
	N1_t = -1*d_rt.Rotate(m_f1_nor);

	SWIFT_Real dist_pt0_s, dist_pt0_t; 
	SWIFT_Real dist_pt1_s, dist_pt1_t; 

	// The distance of the point to the transformed plane
	dist_pt0_s = m_f0_nor * m_pt_0                              + m_f0_d;
	dist_pt0_t = DeltaRt(1).Rotate(m_f0_nor)*(m_pt_0-m_rel_tra) + m_f0_d; // should be the new normal at t=1!!!

	dist_pt1_s = m_f1_nor * m_pt_1                              + m_f1_d;
	dist_pt1_t = DeltaRt(1).Rotate(m_f1_nor)*(m_pt_1-m_rel_tra) + m_f1_d; // should be the new normal at t=1!!!

	s0_t = m_f0_d - N0_t * m_pt_0 - ((1-t)*dist_pt0_s + t*dist_pt0_t);
	s1_t = m_f1_d - N1_t * m_pt_1 - ((1-t)*dist_pt1_s + t*dist_pt1_t);

	/*
	{
		// A big question here
		// It is incorrect
		// lambda = m_axis * m_rel_tra;
		SWIFT_Triple tmp_n0 = DeltaRt(1).Rotate(m_f0_nor);
		SWIFT_Triple tmp_n1 = DeltaRt(1).Rotate(m_f1_nor);
		SWIFT_Triple tmp_n2 = tmp_n0 % tmp_n1;
		tmp_n2.Normalize();

		lambda = tmp_n2 * m_rel_tra;
	}
	*/

	// To set the frame
	// TODO: for VF/VF or FV/FV case, the inverse can be precomputed
	dirIntersect = N0_t % N1_t;
	dirIntersect.Normalize();
	frame.Set_Value_Rows(N0_t, N1_t, dirIntersect);

	base = frame.Inverse() * SWIFT_Triple(-1*s0_t, -1*s1_t, 0);

	lambda = (m_rel_tra - base) * dirIntersect;

	return lambda;
	/*
	delta_tra = base + dirIntersect * t * lambda;

	SWIFT_Triple abs_tra;
	abs_tra = DeltaRt(t).Rotate(transform_s.Translation()) + delta_tra;

	qua[4] = abs_tra.X();
	qua[5] = abs_tra.Y();
	qua[6] = abs_tra.Z();

	transform.Set_Translation(abs_tra);
	*/
}


void CItpMot_TwoContacts::fv_fv_integrate(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Can be precomptued. Be careful DeltaRt is correct
	m_rel_tra = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());

	// The goal is to compute the delta_tra
	// The abs_tra = Delta_Rot * q0.T + delta_tra
	SWIFT_Triple delta_tra;
	SWIFT_Triple base;
	SWIFT_Matrix33 frame;
	SWIFT_Real s0_t, s1_t, lambda;
	SWIFT_Triple N0_t, N1_t, dirIntersect; // The normal of the fi

	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!

	// A trick here:
	// \N \cdot ( \Delta\R(t)^{-1} ( \p - \delta(t)) + d 
	// \Delta\R(t)\N \cdot (\p - \delta(t)) + d // by making use of \R(t) is orthogonal
	N0_t = -1*d_rt.Rotate(m_f0_nor);
	N1_t = -1*d_rt.Rotate(m_f1_nor);

	SWIFT_Real dist_pt0_s, dist_pt0_t; 
	SWIFT_Real dist_pt1_s, dist_pt1_t; 

	// The distance of the point to the transformed plane
	dist_pt0_s = m_f0_nor * m_pt_0                              + m_f0_d;
	dist_pt0_t = DeltaRt(1).Rotate(m_f0_nor)*(m_pt_0-m_rel_tra) + m_f0_d; // should be the new normal at t=1!!!

	dist_pt1_s = m_f1_nor * m_pt_1                              + m_f1_d;
	dist_pt1_t = DeltaRt(1).Rotate(m_f1_nor)*(m_pt_1-m_rel_tra) + m_f1_d; // should be the new normal at t=1!!!

	s0_t = m_f0_d - N0_t * m_pt_0 - ((1-t)*dist_pt0_s + t*dist_pt0_t);
	s1_t = m_f1_d - N1_t * m_pt_1 - ((1-t)*dist_pt1_s + t*dist_pt1_t);

	{
		// A big question here
		// It is incorrect
		// lambda = m_axis * m_rel_tra;
		SWIFT_Triple tmp_n0 = DeltaRt(1).Rotate(m_f0_nor);
		SWIFT_Triple tmp_n1 = DeltaRt(1).Rotate(m_f1_nor);
		SWIFT_Triple tmp_n2 = tmp_n0 % tmp_n1;
		tmp_n2.Normalize();

		lambda = tmp_n2 * m_rel_tra;

		// lambda = tmp_lambda();









	}

	// To set the frame
	// TODO: for VF/VF or FV/FV case, the inverse can be precomputed
	dirIntersect = N0_t % N1_t;
	dirIntersect.Normalize();
	// frame.Set_Value_Cols(N0_t, N1_t, dirIntersect);
	frame.Set_Value_Rows(N0_t, N1_t, dirIntersect); // The frame here is due to the linear system. So set row

	base = frame.Inverse() * SWIFT_Triple(-1*s0_t, -1*s1_t, 0);
	delta_tra = base + dirIntersect * t * lambda;

	SWIFT_Triple abs_tra;
	abs_tra = DeltaRt(t).Rotate(transform_s.Translation()) + delta_tra;

	qua[4] = abs_tra.X();
	qua[5] = abs_tra.Y();
	qua[6] = abs_tra.Z();

	transform.Set_Translation(abs_tra);
}


void CItpMot_TwoContacts::ee_ee_integrate(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();

	// Can be precomptued. Be careful DeltaRt is correct
	m_rel_tra = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());

	// The goal is to compute the delta_tra
	// The abs_tra = Delta_Rot * q0.T + delta_tra
	SWIFT_Triple delta_tra;
	SWIFT_Triple base;
	SWIFT_Matrix33 frame;
	SWIFT_Real s0_t, s1_t;
	SWIFT_Triple N0_t, N1_t, dirIntersect; // The normal of the fi

	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!

	SWIFT_Triple base_n0, base_n1;

	// The normal is from e_obs to e_rob
	base_n0 = (m_pt_3-m_pt_2) % (m_pt_1-m_pt_0);
	base_n0.Normalize();
	base_n1 = (m_pt_7-m_pt_6) % (m_pt_5-m_pt_4);
	base_n1.Normalize();

	SWIFT_Real lambda;
	{
		// A big question here
		// It is incorrect
		// lambda = m_axis * m_rel_tra;

		// Rotate the normal is incorrect!!!; // R(t)V1 % V2 != R(t) (V1%v2) !!!
		// SWIFT_Triple tmp_n0 = DeltaRt(1).Rotate(base_n0);
		// SWIFT_Triple tmp_n1 = DeltaRt(1).Rotate(base_n1);
		SWIFT_Triple tmp_n0 = (m_pt_3-m_pt_2) % DeltaRt(1).Rotate(m_pt_1 - m_pt_0);
		SWIFT_Triple tmp_n1 = (m_pt_7-m_pt_6) % DeltaRt(1).Rotate(m_pt_5 - m_pt_4);

		SWIFT_Triple tmp_n2 = tmp_n0 % tmp_n1;
		tmp_n2.Normalize();

		lambda = tmp_n2 * m_rel_tra;
	}

	// R(t)V1 % V2 != R(t) (V1%v2) !!!
	N0_t = d_rt.Rotate(base_n0);
	N0_t = (m_pt_3-m_pt_2) % d_rt.Rotate(m_pt_1 - m_pt_0);
	N0_t.Normalize(); // To normalize it !!!
	N1_t = d_rt.Rotate(base_n1);
	N1_t = (m_pt_7-m_pt_6) % d_rt.Rotate(m_pt_5 - m_pt_4);
	N1_t.Normalize();

	SWIFT_Real dist_ee0_s, dist_ee0_t; 
	SWIFT_Real dist_ee1_s, dist_ee1_t; 

	// The distance of the point to the transformed plane
/*	dist_pt0_s = m_f0_nor * m_pt_0                              + m_f0_d;
	dist_pt0_t = DeltaRt(1).Rotate(m_f0_nor)*(m_pt_0-m_rel_tra) + m_f0_d; // should be the new normal at t=1!!!

	dist_pt1_s = m_f1_nor * m_pt_1                              + m_f1_d;
	dist_pt1_t = DeltaRt(1).Rotate(m_f1_nor)*(m_pt_1-m_rel_tra) + m_f1_d; // should be the new normal at t=1!!!
*/
	// Todo: to compute the distance between the edge
	// Can make use the normal
	dist_ee0_s = LineLineSignDist(m_pt_0, m_pt_1, m_pt_2, m_pt_3);
	// Transform the vertex first
	// It is not the most efficient way
	dist_ee0_t = LineLineSignDist(DeltaRt(1).Rotate(m_pt_0)+m_rel_tra, 
						      DeltaRt(1).Rotate(m_pt_1)+m_rel_tra, 
							  m_pt_2, m_pt_3);

	dist_ee1_s = LineLineSignDist(m_pt_4, m_pt_5, m_pt_6, m_pt_7);
	dist_ee1_t = LineLineSignDist(DeltaRt(1).Rotate(m_pt_4)+m_rel_tra, 
						      DeltaRt(1).Rotate(m_pt_5)+m_rel_tra, 
							  m_pt_6, m_pt_7);

	// Rob
	s0_t = N0_t * (d_rt.Rotate(m_pt_0) - m_pt_2) - ((1-t)*dist_ee0_s + t*dist_ee0_t);
	s1_t = N1_t * (d_rt.Rotate(m_pt_4) - m_pt_6) - ((1-t)*dist_ee1_s + t*dist_ee1_t);

	// To set the frame
	// TODO: for VF/VF or FV/FV case, the inverse can be precomputed
	dirIntersect = N0_t % N1_t;
	dirIntersect.Normalize();
	// frame.Set_Value_Cols(N0_t, N1_t, dirIntersect);
	frame.Set_Value_Rows(N0_t, N1_t, dirIntersect); // The frame here is due to the linear system. So set row

	base = frame.Inverse() * SWIFT_Triple(-1*s0_t, -1*s1_t, 0);
	delta_tra = base + dirIntersect * t * lambda;

	SWIFT_Triple abs_tra;
	abs_tra = DeltaRt(t).Rotate(transform_s.Translation()) + delta_tra;

	qua[4] = abs_tra.X();
	qua[5] = abs_tra.Y();
	qua[6] = abs_tra.Z();

	transform.Set_Translation(abs_tra);
}

///////////////////////
CItpMot_MultContacts::CItpMot_MultContacts(const PQP_REAL q0[7], const PQP_REAL q1[7] ) : CInterpMotion(GMP_IM_ALLMODE, q0, q1)
{
	velocity();
	m_nConstraints = 0;
	m_DistConsts = NULL;
}


void CItpMot_MultContacts::velocity(void)
{
	LinearAngularVel(m_axis, m_angVel);

	// Can be precomptued. Be careful DeltaRt is correct
	m_rel_tra = transform_t.Translation() - DeltaRt(1).Rotate(transform_s.Translation());
	m_rel_rot = DeltaRt(1);
}

bool CItpMot_MultContacts::integrate(const double t, PQP_REAL qua[7])
{
	bool bRet;
	if(m_nConstraints == 2)
		bRet = TwoConsts_integrate(t, qua);
	else if(m_nConstraints == 3)
		bRet = ThreeConsts_integrate(t, qua);

	return bRet;
}


// Does not work very well, maybe you should change the change in the direction of $x$ and $y$ is the direction of Z
//bool bTestNewDistancFun = true;
//SWIFT_Real d_ncr_full, d_ncr = 1;


// To seperate the code
bool CItpMot_MultContacts::Degenerate2D(const SWIFT_Triple& N_0, 
									    const SWIFT_Triple& N_1)
{
	return true;
}


bool CItpMot_MultContacts::TwoConsts_integrate(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();


	// The goal is to compute the delta_tra
	// The abs_tra = Delta_Rot * q0.T + delta_tra
	SWIFT_Triple delta_tra;
	SWIFT_Triple base;
	SWIFT_Matrix33 frame;
	SWIFT_Real s0_t, s1_t, lambda;
	SWIFT_Triple N0_t, N1_t, dirIntersect; // The normal of intersection of the planes 0 and 1

	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!

	m_DistConsts[0]->Compute_N_t(t, d_rt, N0_t);
	m_DistConsts[1]->Compute_N_t(t, d_rt, N1_t);
	// N1_t = m_f1_nor;
	// N0_t = m_f0_nor;

	// To set the frame
	// TODO: for VF/VF or FV/FV case, the inverse can be precomputed
	dirIntersect = N0_t % N1_t;
	//if(bTestNewDistancFun)
	//	d_ncr = dirIntersect.Length();

	if(!MyNormalize(dirIntersect))
		return false;


	// dirIntersect.Normalize();
	// a bug frame.Set_Value_Cols(N0_t, N1_t, dirIntersect); 
	frame.Set_Value_Rows(N0_t, N1_t, dirIntersect); // The frame here is due to the linear system. So set rows

	/*
	SWIFT_Real dist_pt0_s, dist_pt0_t; 
	SWIFT_Real dist_pt1_s, dist_pt1_t; 
	// The distance of the closest feature to the face at q0 and q1
	SWIFT_Triple pt_0_t, pt_1_t;
	pt_0_t = DeltaRt(1).Rotate(m_pt_0) + m_rel_tra; // should be the new point at t=1!!!
	pt_1_t = DeltaRt(1).Rotate(m_pt_1) + m_rel_tra;

	dist_pt0_s = m_f0_nor * m_pt_0 + m_f0_d;
	dist_pt0_t = m_f0_nor * pt_0_t + m_f0_d;

	dist_pt1_s = m_f1_nor * m_pt_1 + m_f1_d;
	dist_pt1_t = m_f1_nor * pt_1_t + m_f1_d;

	s0_t = m_f0_nor * d_rt.Rotate(m_pt_0) + m_f0_d - ((1-t)*dist_pt0_s + t*dist_pt0_t);
	s1_t = m_f1_nor * d_rt.Rotate(m_pt_1) + m_f1_d - ((1-t)*dist_pt1_s + t*dist_pt1_t);
	*/

	s0_t = m_DistConsts[0]->Compute_s_t(t, d_rt);
	s1_t = m_DistConsts[1]->Compute_s_t(t, d_rt);


	SWIFT_Triple N_cr_full, N_0_full, N_1_full;

	m_DistConsts[0]->Compute_N_t(1, m_rel_rot, N_0_full);
	m_DistConsts[1]->Compute_N_t(1, m_rel_rot, N_1_full);
	N_cr_full = N_0_full % N_1_full;

	//if(bTestNewDistancFun)
	//	d_ncr_full = N_cr_full.Length();


	// To check whether the constraints are parallel
	if(N_cr_full.Length_Sq() < EPS_PLANE_PARRAL) // a parallel constraint
		return false;
	if(!MyNormalize(N_cr_full))
		return false;


	lambda = m_rel_tra * N_cr_full;

	if(fabs(MyDeterminant(frame, 0, 1, 2, 0, 1, 2)) < 0.00001f)
	{
		GMP_ASSERT("Degnerate sitatuions in 2 Distance Constraints");
		return false;
	}

	base = frame.Inverse() * SWIFT_Triple(-1*s0_t, -1*s1_t, 0);
	delta_tra = base + dirIntersect * t * lambda ;


	// delta_tra = base + dirIntersect * t * lambda / d_ncr;
	// printf("t=%.3f, d_ncr=%.3f, old t=%.3f, new t=%.3f\n", t, dirIntersect.Length(),  t * lambda / d_ncr);
	//printf("%.3f %.3f %.3f %.3f\n", t, d_ncr,  t * lambda, t * lambda * d_ncr);



	//if(bTestNewDistancFun)
	//	delta_tra = base + dirIntersect * t * lambda * d_ncr * d_ncr / (d_ncr_full *d_ncr_full);
	//else
	//	delta_tra = base + dirIntersect * t * lambda;

	SWIFT_Triple abs_tra;
	SWIFT_Triple tmp_tra = d_rt.Rotate(transform_s.Translation());

	// abs_tra = DeltaRt(t).Rotate(transform_s.Translation()) + delta_tra;
	// abs_tra = d_rt.Rotate(transform_s.Translation()) + delta_tra;
	abs_tra = tmp_tra + delta_tra;

	qua[4] = abs_tra.X();
	qua[5] = abs_tra.Y();
	qua[6] = abs_tra.Z();

	transform.Set_Translation(abs_tra);

	return true;
}

bool CItpMot_MultContacts::Degenerate3D(const SWIFT_Triple& N_0, 
				  const SWIFT_Triple& N_1,
				  const SWIFT_Triple& N_2, 
				  const SWIFT_Matrix33 &conSystem)
{
	if(m_nConstraints != 3)
		GMP_ASSERT("Wrong here for Degenerate3D!!!");

	// Rank is 3, return true;
	if(fabs(MyDeterminant(conSystem, 0, 1, 2, 0, 1, 2)) > 0.01f)
		return false;

	SWIFT_Triple N[3];

	N[0] = N_0;
	N[1] = N_1;
	N[2] = N_2;

	for(int i=0; i<3; i++)
	{
		if(MyNormalize(N[i] % N[(i+1)%3]))
		{
			m_thrDeg_Valid_0 = i; 
			m_thrDeg_Valid_1 = (i+1)%3;

			m_thrDeg_nValid = 2;

			return true; // a degenerate case 
		}
	}

	m_thrDeg_Valid_0 = 0;
	m_thrDeg_nValid = 1;

	return true; // a degenerate case 
}

bool CItpMot_MultContacts::ThreeConsts_integrate(const double t, PQP_REAL qua[7])
{
	SWIFT_Quaternion predictedOrn = AbsoluteRt(t);
	transform.Set_Rotation(predictedOrn);
	qua[0] = predictedOrn.W();
	qua[1] = predictedOrn.X();
	qua[2] = predictedOrn.Y();
	qua[3] = predictedOrn.Z();


	// The goal is to compute the delta_tra
	// The abs_tra = Delta_Rot * q0.T + delta_tra
	SWIFT_Triple delta_tra;
	SWIFT_Matrix33 system;
	SWIFT_Triple N0_t, N1_t, N2_t;
	SWIFT_Real s0_t, s1_t, s2_t;

	SWIFT_Quaternion d_rt = DeltaRt(t); // This is the relative amount of rotation !!!

	m_DistConsts[0]->Compute_N_t(t, d_rt, N0_t);
	m_DistConsts[1]->Compute_N_t(t, d_rt, N1_t);
	m_DistConsts[2]->Compute_N_t(t, d_rt, N2_t);

	s0_t = m_DistConsts[0]->Compute_s_t(t, d_rt);
	s1_t = m_DistConsts[1]->Compute_s_t(t, d_rt);
	s2_t = m_DistConsts[2]->Compute_s_t(t, d_rt);

	system.Set_Value_Rows(N0_t, N1_t, N2_t); // The frame here is due to the linear system. So set rows

	// Degenerate situation 
	if(Degenerate3D(N0_t, N1_t, N2_t, system))
		return false;
	else
		delta_tra = system.Inverse() * SWIFT_Triple(-s0_t, -s1_t, -s2_t);

	SWIFT_Triple abs_tra;
	SWIFT_Triple tmp_tra = d_rt.Rotate(transform_s.Translation());

	// abs_tra = DeltaRt(t).Rotate(transform_s.Translation()) + delta_tra;
	// abs_tra = d_rt.Rotate(transform_s.Translation()) + delta_tra;
	abs_tra = tmp_tra + delta_tra;

	qua[4] = abs_tra.X();
	qua[5] = abs_tra.Y();
	qua[6] = abs_tra.Z();

	transform.Set_Translation(abs_tra);

	return true;
}

string CItpMot_MultContacts::Print()
{
	char oStr[256];

	sprintf(oStr, "%d DistConst: ", m_nConstraints);

	for(int i=0; i<m_nConstraints; i++)
	{
		strcat(oStr, ", ");
		if(m_DistConsts[i]->m_distConType == GMP_DC_V_F)
			strcat(oStr, "V-F");
		else if(m_DistConsts[i]->m_distConType == GMP_DC_F_V)
			strcat(oStr, "F-V");
		else if(m_DistConsts[i]->m_distConType == GMP_DC_E_E)
			strcat(oStr, "E-E");
		else
			strcat(oStr, "Unknown");
	}
	return oStr;
}

void CItpMot_MultContacts::DrawFeatures(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale)
{
	for(int i=0; i<m_nConstraints; i++)
		m_DistConsts[i]->Draw(robR, robT, obsR, obsT, scale);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance Constraint

SWIFT_Real CVF_DC::Compute_s_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt)
{
	SWIFT_Real s_t;
	SWIFT_Triple pt_t;

	// TODO: 
	// Can be passed here
	// m_rel_tra
	// d_rt
	// DeltaRt(1)
	CItpMot_MultContacts *pItp = (CItpMot_MultContacts *)m_pOwner_Itp;
	pt_t = pItp->m_rel_rot.Rotate(m_A_pt) + pItp->m_rel_tra; // should be the new point at t=1!!!

	SWIFT_Real dist_pt_s = m_B_f_nor * m_A_pt + m_B_f_d;
	SWIFT_Real dist_pt_t = m_B_f_nor * pt_t + m_B_f_d;

	s_t = m_B_f_nor * d_rt.Rotate(m_A_pt) + m_B_f_d - ((1-t)*dist_pt_s + t*dist_pt_t);

	return s_t;
}

void CVF_DC::Compute_N_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt, 
						 SWIFT_Triple& N_t)
{
	N_t = m_B_f_nor;
}


SWIFT_Real CFV_DC::Compute_s_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt)
{
	SWIFT_Real s_t;
	CItpMot_MultContacts *pItp = (CItpMot_MultContacts *)m_pOwner_Itp;

	// A trick here:
	// \N \cdot ( \Delta\R(t)^{-1} ( \p - \delta(t)) + d 
	// \Delta\R(t)\N \cdot (\p - \delta(t)) + d // by making use of \R(t) is orthogonal
	SWIFT_Real dist_pt_s, dist_pt_t; 

	// The distance of the point to the transformed plane
	dist_pt_s = m_A_f_nor * m_B_pt								           + m_A_f_d;
	dist_pt_t = pItp->m_rel_rot.Rotate(m_A_f_nor)*(m_B_pt-pItp->m_rel_tra) + m_A_f_d; // should be the new normal at t=1!!!

	// Inefficient !!!
	SWIFT_Triple N_t;
	Compute_N_t(t, d_rt, N_t);
	s_t = m_A_f_d - N_t * m_B_pt - ((1-t)*dist_pt_s + t*dist_pt_t);

	return s_t;
}

void CFV_DC::Compute_N_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt, 
						 SWIFT_Triple& N_t)
{
	N_t = -1*d_rt.Rotate(m_A_f_nor);
}

SWIFT_Real CEE_DC::Compute_s_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt)
{
	SWIFT_Real s_t;
	SWIFT_Real dist_ee_s, dist_ee_t; 

	CItpMot_MultContacts *pItp = (CItpMot_MultContacts *)m_pOwner_Itp;

	// Todo: precompute 
	dist_ee_s = LineLineSignDist(m_A_pt_0, m_A_pt_1, m_B_pt_0, m_B_pt_1);
	// Transform the vertex first 
	// It is not the most efficient way
	dist_ee_t = LineLineSignDist(pItp->m_rel_rot.Rotate(m_A_pt_0)+pItp->m_rel_tra, 
							     pItp->m_rel_rot.Rotate(m_A_pt_1)+pItp->m_rel_tra, 
							     m_B_pt_0, m_B_pt_1);

	// Inefficient !!!
	SWIFT_Triple N_t;
	Compute_N_t(t, d_rt, N_t);
	s_t = N_t * (d_rt.Rotate(m_A_pt_0) - m_B_pt_0) - ((1-t)*dist_ee_s + t*dist_ee_t);
	
	return s_t;
}

void CEE_DC::Compute_N_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt, 
						 SWIFT_Triple& N_t)
{
	// R(t)V1 % V2 != R(t) (V1%v2) !!!
	N_t = (m_B_pt_1-m_B_pt_0) % d_rt.Rotate(m_A_pt_1-m_A_pt_0);
	MyNormalize(N_t);
	// N_t.Normalize(); // To normalize it !!!
}

void CEE_DC::Draw(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale)
{
	glPushAttrib(GL_ENABLE_BIT);
	glPushAttrib(GL_LINE_BIT);

	glLineWidth(scale*60);
	glDisable(GL_LIGHTING);

	SWIFT_Triple pos0, pos1;
	SWIFT_Matrix33 R;
	SWIFT_Triple T;

	pos0 = robR * m_A_pt_0 + robT;
	pos1 = robR * m_A_pt_1 + robT;
	glColor3f(1, 0, 1);
	glBegin(GL_LINES);
		glVertex3dv(pos0.Value());
		glVertex3dv(pos1.Value());
	glEnd();

	pos0 = obsR * m_B_pt_0 + obsT;
	pos1 = obsR * m_B_pt_1 + obsT;
	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
		glVertex3dv(pos0.Value());
		glVertex3dv(pos1.Value());
	glEnd();

	glPopAttrib();
	glPopAttrib();
}


void DrawVertex_AsSphere(const SWIFT_Triple &lcs, double scale)
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslated(lcs.X(), lcs.Y(), lcs.Z());
	glutSolidSphere(0.5*scale, 60, 60);
	glPopMatrix();
}


void CVF_DC::Draw(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale)
{
	// Robot vertex
	SWIFT_Triple rob_pt = robR * m_A_pt + robT;
	glColor3f(1, 0, 1);
	DrawVertex_AsSphere(rob_pt, scale);

	// Obstacle Face
	SWIFT_Triple proj;

	proj = rob_pt - (m_B_f_nor * rob_pt + m_B_f_d) * m_B_f_nor;
	glColor3f(0, 0, 1);
	DrawVertex_AsSphere(proj, scale);

	glBegin(GL_LINES);
		glVertex3dv(proj.Value());
		glVertex3dv(rob_pt.Value());
	glEnd();
}

// DrawPrim is better for visualization
// pRob->DrawPrim(*q_iter, SWIFT_VERTEX, pOneContact->m_vf_pt_id, gUsingPQP_ForPDG);
void CFV_DC::Draw(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale)
{
	// Obs vertex
	glColor3f(0, 0, 1);
	SWIFT_Triple obs_pt = obsR * m_B_pt + obsT;
	DrawVertex_AsSphere(obs_pt, scale);

	// Robot Face
	SWIFT_Triple proj;

	glColor3f(1, 0, 1);
	proj = obs_pt - (m_A_f_nor * obs_pt + m_A_f_d) * m_A_f_nor;
	DrawVertex_AsSphere(proj, scale);
	
	glBegin(GL_LINES);
		glVertex3dv(proj.Value());
		glVertex3dv(obs_pt.Value());
	glEnd();
}

SWIFT_Real DistanceInterpolant(SWIFT_Real t)
{
	GMP_ASSERT("F-V for single constraint is wrong");
	return t;
}
