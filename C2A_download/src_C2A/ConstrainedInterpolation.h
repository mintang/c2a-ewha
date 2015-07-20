#ifndef CONSTRAINED_MOTION_INTERPOLATION
#define CONSTRAINED_MOTION_INTERPOLATION

// Linear angular velocity
// The shortest distance between the contact feature is linearly interpolated. 

class CItpMot_MultContacts;

enum GMP_ITP_ONE_CONTACT
{
	GMP_ITP_ONE_CONTACT_V_F = 0,   // 
	GMP_ITP_ONE_CONTACT_F_V = 1,   // 
	GMP_ITP_ONE_CONTACT_E_E = 2,   // 
};

void LineLineNearPts(const SWIFT_Triple & A_worl_1, const SWIFT_Triple & A_worl_2, 
				 const SWIFT_Triple & B_1,      const SWIFT_Triple & B_2,  
				 SWIFT_Triple & near_p1,        SWIFT_Triple & near_p2);


class CItpMot_OneContact : public CInterpMotion
{
public:
  CItpMot_OneContact(const PQP_REAL q0[7], const PQP_REAL q1[7]);
  CItpMot_OneContact();


  GMP_ITP_ONE_CONTACT m_contact_type; 

  virtual ~CItpMot_OneContact();
  virtual void velocity(void);
  virtual bool integrate(const double dt, PQP_REAL qua[7]);
  virtual string Print();

  void SetupEE(SWIFT_Triple A0, SWIFT_Triple A1,
			   SWIFT_Triple B0, SWIFT_Triple B1);

  // For ee
  // Refer to the paper
  SWIFT_Triple m_pt_A_0, m_pt_A_1; // Robot's edge: the position of two end poitns at q_0
  SWIFT_Triple m_pt_B_0, m_pt_B_1; // Obstacle's edge: the position of two end poitns

  // Pre-computed according to the input
  SWIFT_Triple m_near_A_0_0, m_near_A_1_1; // Can be pre-computed
  SWIFT_Triple m_near_A_0_1; 
  SWIFT_Triple m_near_B_0, m_near_B_1; // Can be pre-computed
  SWIFT_Real m_d_0, m_d_1;

  // there are issues with the old implementaiton of EE
  //   1. m_obs_p0 may be same as m_obs_p1, which results in the bad normal computaiton 
  //   2. sign distances
  // same notation                 // for E/E situation
  SWIFT_Triple m_rob_m_0_0, m_rob_m_1_1; // vertices on the line of the robot
  SWIFT_Triple m_obs_p0, m_obs_p1; // vertices on the line of the obstacle


  SWIFT_Triple m_vf_pt;     // for V/F case 
  int m_vf_pt_id;

  // for F/V case 
  SWIFT_Triple m_fv_obs_pt; // the point of the obstacle
  SWIFT_Triple m_fv_rob_nor_0; // N0
  SWIFT_Triple m_fv_rob_m_0_0, m_fv_rob_m_1_1; // the closest points on the plane containing the face: 
  bool m_fv_sign_0_flip; // It is tricky. 
  bool m_fv_sign_1_flip; // For signed distance


  // the code should be more merged with multiple constraints
  CItpMot_MultContacts *m_pMultConstraint;



private:
  void vf_integrate(const double dt, PQP_REAL qua[7]);
  void ee_integrate(const double dt, PQP_REAL qua[7]);
  void fv_integrate_SPM08(const double dt, PQP_REAL qua[7]);

  // It is buggy; to remove
  void ee_integrate_SPM08_temp(const double dt, PQP_REAL qua[7]);
};


enum GMP_ITP_TWO_CONTACTS
{
	GMP_ITP_TWO_CONTACTS_VF_VF = 0,   // 
	GMP_ITP_TWO_CONTACTS_FV_FV = 1,   //
	GMP_ITP_TWO_CONTACTS_EE_EE = 2,   
	// GMP_ITP_ONE_CONTACT_VF_V = 1,   // 
	// GMP_ITP_ONE_CONTACT_E_E = 2,   // 
};

class CItpMot_ThreeContacts : public CInterpMotion
{
public:
  //CItpMot_ThreeContacts(const PQP_REAL q0[7], const PQP_REAL q1[7]);
  //CItpMot_ThreeContacts();
};

class CItpMot_TwoContacts : public CInterpMotion
{
public:
  CItpMot_TwoContacts(const PQP_REAL q0[7], const PQP_REAL q1[7]);
  CItpMot_TwoContacts();

  GMP_ITP_TWO_CONTACTS m_contact_type; 

  virtual ~CItpMot_TwoContacts();
  virtual void velocity(void);
  virtual bool integrate(const double dt, PQP_REAL qua[7]);

  SWIFT_Triple m_vf_vf_pt_1,   m_vf_vf_pt_2;     // for VF, vf, case 
  SWIFT_Triple m_vf_vf_f1_nor, m_vf_vf_f2_nor;   
  SWIFT_Real   m_vf_vf_f1_d,   m_vf_vf_f2_d;   


  SWIFT_Triple m_pt_0;   // The position of the vertex (at q_0 for robot) (for obstacle, ready for using) 
  SWIFT_Triple m_f0_nor; // The normal of the face (at q_0 for robot) (for obstacle, ready for using) 
  SWIFT_Real   m_f0_d;

  SWIFT_Triple m_pt_1;
  SWIFT_Triple m_f1_nor;
  SWIFT_Real   m_f1_d;


  SWIFT_Triple m_pt_2, m_pt_3;


  SWIFT_Triple m_pt_4, m_pt_5;
  SWIFT_Triple m_pt_6, m_pt_7;



private:
  void vf_vf_integrate(const double dt, PQP_REAL qua[7]);
  void vf_vf_integrate_new(const double dt, PQP_REAL qua[7]);
  void fv_fv_integrate(const double t, PQP_REAL qua[7]);
  void ee_ee_integrate(const double t, PQP_REAL qua[7]);


  SWIFT_Real tmp_lambda();


  // void ee_integrate(const double dt, PQP_REAL qua[7]);
  // void fv_integrate(const double dt, PQP_REAL qua[7]);
};

//////////////////////////////////////
// Distance Constraint and CItpMot_MultContacts

enum GMP_DIST_CONSTRAINT
{
	GMP_DC_V_F = 0,   // 
	GMP_DC_F_V= 1,   // 
	GMP_DC_E_E = 2,   // 
};

class CDistanceConstraint
{
public:
	GMP_DIST_CONSTRAINT m_distConType;
	CInterpMotion *m_pOwner_Itp;

	virtual SWIFT_Real Compute_s_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt) = 0;
	virtual void Compute_N_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt,
					 SWIFT_Triple& N_t) = 0;
	virtual void Draw(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale) = 0;
};

class CVF_DC : public CDistanceConstraint
{
public:
	CVF_DC()
	{
		m_distConType = GMP_DC_V_F;
	}

	SWIFT_Triple m_A_pt; // the point of the robot at configuration q0
	SWIFT_Triple m_B_f_nor; // The face of the obstacle
	SWIFT_Real m_B_f_d;

	SWIFT_Real Compute_s_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt);
	void Compute_N_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt,
					 SWIFT_Triple& N_t);
	virtual void Draw(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale);
};

class CFV_DC : public CDistanceConstraint
{
public:
	CFV_DC()
	{
		m_distConType = GMP_DC_F_V;
	}
	SWIFT_Triple m_A_f_nor; // The face of the robot at configuration q0
	SWIFT_Real   m_A_f_d;
	SWIFT_Triple m_B_pt; // the point of the obstacle

	SWIFT_Real Compute_s_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt);
	void Compute_N_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt,
					 SWIFT_Triple& N_t);
	virtual void Draw(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale);
};

class CEE_DC : public CDistanceConstraint
{
public:
	CEE_DC()
	{
		m_distConType = GMP_DC_E_E;
	}

    // rob edge = m_A_pt_0, m_A_pt_1
    // obs edge = m_B_pt_0, m_A_pt_1
	SWIFT_Triple m_A_pt_0, m_A_pt_1;
	SWIFT_Triple m_B_pt_0, m_B_pt_1;

	SWIFT_Real Compute_s_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt);
	void Compute_N_t(const SWIFT_Real t, const SWIFT_Quaternion &d_rt,
					 SWIFT_Triple& N_t);

	virtual void Draw(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale);
};

// at t=0, return 0
// at t=1, return 1
SWIFT_Real DistanceInterpolant(SWIFT_Real t);



// Can be 2 or 3
class CItpMot_MultContacts : public CInterpMotion
{
public:
  CItpMot_MultContacts(const PQP_REAL q0[7], const PQP_REAL q1[7]);
  CItpMot_MultContacts()
  {
		m_nConstraints = 0;
		m_DistConsts = NULL;
  }

  virtual ~CItpMot_MultContacts()
  {
	  for(int i=0; i<m_nConstraints; i++)
		  delete m_DistConsts[i];
	  delete m_DistConsts;
  }

  virtual void velocity(void);
  virtual bool integrate(const double dt, PQP_REAL qua[7]);

  virtual void DrawFeatures(const SWIFT_Matrix33 &robR, const SWIFT_Triple &robT, 
				  const SWIFT_Matrix33 &obsR, const SWIFT_Triple &obsT, 
				  double scale);




  virtual string Print();

  int m_nConstraints;
  CDistanceConstraint **m_DistConsts;

public:
  // the relative translation between q0 and q1
  // the relative rotation between q0 and q1
  // Setup by function velocity()
  SWIFT_Triple m_rel_tra;
  SWIFT_Quaternion m_rel_rot;

  // For handle the degenerate situation of three constraints
	int m_thrDeg_nValid;
	int m_thrDeg_Valid_0;
	int m_thrDeg_Valid_1;



private:
  bool TwoConsts_integrate(const double t, PQP_REAL qua[7]);
  bool ThreeConsts_integrate(const double t, PQP_REAL qua[7]);

  bool Degenerate2D(const SWIFT_Triple& N_0, 
				    const SWIFT_Triple& N_1);

  bool Degenerate3D(const SWIFT_Triple& N_0, 
				    const SWIFT_Triple& N_1,
				    const SWIFT_Triple& N_2, 
				    const SWIFT_Matrix33 &conSystem);



};


#endif