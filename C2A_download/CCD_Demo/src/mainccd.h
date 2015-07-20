#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include "model.h"
#include "PQP.h"
#include "PQP_CCD.h"
#include "InterpMotion.h"
#include "MatVec.h"
#include <iostream>
#include <fstream>
#include "Transform.h"
#include "lut.h"

#define GetContactNumber false


void readXFormFrames( const char* filename1, const char* filename2 );
void outputPerformance(const char* filename );
void outputXformTOCs(const char* filename );

struct TransRT
{
	PQP_REAL R[3][3];	
	PQP_REAL T[3];	
};
TransRT xformframes[2][1000];
TransRT xformTOCs[1000];
float tframes[1000];
float tocframes[1000];
float distframes[1000];
int iterframes[1000];
float contframes[1000];
int FrontNodes[1000];

static const char* timingFileName = "../results/timingPQPCCD_Dir.txt";
PQP_CollideResult *cres;
//PQP_DistanceResult dres;
PQP_TOCResult dres;

LUT object2_lut, object1_lut;




void outputPerformance(const char* filename )
{
	  

	ofstream fout(filename);


	fout<<endl<<"Frame No.| Timing |  Iterative |   TOC   |   Distance   | Contact  | FrontNode";
	fout<<endl;

	for(int i=0; i<nframes-1; i++ )
	{
        if(tocframes[i]==0.0)
			fout<<i<<"   \t "<<tframes[i]*1000<<"   \t "<<iterframes[i]<<"   \t "<<tocframes[i]<<"    \t "<<distframes[i]<<"      \t "<<contframes[i]<<"      \t "<<FrontNodes[i]<<endl;
		else
			fout<<i<<"   \t "<<tframes[i]*1000<<"   \t"<<iterframes[i]<<"   \t "<<tocframes[i]<<"   \t "<<distframes[i]<<"      \t "<<contframes[i]<<"      \t "<<FrontNodes[i]<<endl;
	}
	cout<<"Stop!"<<endl;

}

void outputXformTOCs(const char* filename )
{
	ofstream fout(filename);
	TransRT *xform;

	for(int i=0; i<nframes-1; i++ )
	{
		xform=&xformTOCs[i];

		fout<<i<<" "<<xform->R[0][0]<<" "
					<<xform->R[0][1]<<" "
					<<xform->R[0][2]<<" "
					<<xform->R[1][0]<<" "
					<<xform->R[1][1]<<" "
					<<xform->R[1][2]<<" "
					<<xform->R[2][0]<<" "
					<<xform->R[2][1]<<" "
					<<xform->R[2][2]<<" "
					<<xform->T[0]<<" "
					<<xform->T[1]<<" "
					<<xform->T[2]<<endl;
	}
	cout<<"Finish outputting xforms at TOC!"<<endl;
}

inline void Transform2PQP(Transform *t, PQP_REAL R[3][3], PQP_REAL T[3])
{
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			R[i][j] = t->basis[i][j];
		}
		T[i] = t->origin[i];
	}	
}

inline void PQP2Transform(PQP_REAL R[3][3], PQP_REAL T[3], Transform &t)
{
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			t.basis[i][j]= R[i][j];
		}
		t.origin[i]= T[i];
	}	
}

 inline  void output_list_New(ContactF_LIST ConL)
	{
		ContactF_IT ContactFit;
		for (ContactFit = ConL.begin(); ContactFit != ConL.end(); ContactFit++)
		{
			cout << (*ContactFit).FeatureType_A << '\t' << (*ContactFit).FeatureType_B
				<<'\t' << (*ContactFit).Distance<<'\n' 
				<< (*ContactFit).TriangleID_A<<'\t' << (*ContactFit).FeatureID_A[0]<<'\t'<< (*ContactFit).FeatureID_A[1]<<'\t' << (*ContactFit).FeatureID_A[2]<<'\n'
				<< (*ContactFit).TriangleID_B<<'\t' << (*ContactFit).FeatureID_B[0]<<'\t'<< (*ContactFit).FeatureID_B[1]<<'\t' << (*ContactFit).FeatureID_B[2]<< endl;
		}
	}


//CA toc
RESULT CA(Transform *trans00,Transform *trans01,  PQP_Model * obj1_tested,
		  Transform *trans10, Transform *trans11,  PQP_Model * obj2_tested,
		  Transform & trans0, Transform & trans1, PQP_REAL &toc, int& nItrs, int & NTr)
{
	// linear velocity and angular velocity
	Vector3D	cv0, cw0, cv1, cw1;

	trans00->velocity(trans01, 1.0, cv0, cw0);
	trans10->velocity(trans11, 1.0, cv1, cw1);

	PQP_REAL Vel_relate[3];
	Vector3D temp1,temp2,temp;
	temp1 = trans01->origin-trans00->origin;
	temp2 = trans11->origin-trans10->origin;
	temp  = temp1 - temp2;

	Vel_relate[0]=temp.x;
	Vel_relate[1]=temp.y;
	Vel_relate[2]=temp.z;


	nItrs=0;
	Transform interpolatedXform[2];
	
	PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3];

	PQP_REAL dR1[12], dR2[12];	
	
	Transform2PQP(trans00, R1, T1);
	Transform2PQP(trans10, R2, T2);

	int k = 0;
	int i,j;

	for( i = 0; i < 3; i++)
	{
		for( j = 0; j < 3; j++)
		{
			dR1[k] = R1[i][j];			
			dR2[k] = R2[i][j];	
			k++;
		}
	}

	

	PQP_REAL R1e[3][3], T1e[3], R2e[3][3], T2e[3];	

	// get the rotation matrix and translational vector about the final state
	Transform2PQP(trans01, R1e, T1e);
	Transform2PQP(trans11, R2e, T2e);
	
	k = 0;
	for( i = 0; i < 3; i++)
	{
		for( j = 0; j < 3; j++)
		{
			dR1[k] = R1e[i][j];			
			dR2[k] = R2e[i][j];	
			k++;
		}
	}
  

   //perform toc query

	CInterpMotion_Linear motion1(R1, T1, R1e, T1e);
	CInterpMotion_Linear motion2(R2, T2, R2e, T2e);

	motion1.velocity_relate(Vel_relate);


	PQP_REAL t_delta = 0.0001;
	PQP_REAL d_delta = 0.0001;

	motion1.m_toc_delta = d_delta;
	motion2.m_toc_delta = d_delta;
	motion1.m_time=t_delta;
    motion2.m_time=t_delta;

    //calculate the toc 
	PQP_QueryTimeOfContact( &motion1, &motion2, &dres,obj1_tested, 
						 obj2_tested, d_delta, t_delta, 0);
	
	

    dres.cont_l.clear();

	if(!dres.collisionfree)
	{
		PQP_REAL qua[7];
		Quaternion q;
		Coord3D org;		

		motion1.integrate(dres.toc, qua);
		q.w = qua[0]; q.x = qua[1]; q.y = qua[2]; q.z = qua[3]; 
		org.x = qua[4]; org.y = qua[5]; org.z = qua[6];
		trans0.setRotation(q);
		trans0.setOrigin(org);

		motion2.integrate(dres.toc, qua);
		q.w = qua[0]; q.x = qua[1]; q.y = qua[2]; q.z = qua[3]; 
		org.x = qua[4]; org.y = qua[5]; org.z = qua[6];
		trans1.setRotation(q);
		trans1.setOrigin(org);

		  //used to get the contact feature
	
		if (GetContactNumber)
		{	
			double threshold = 2*dres.distance+0.001;		
        	PQP_QueryContact(&motion1,&motion2,&dres,obj1_tested,obj2_tested,threshold);  

		}		 	
	
		
	}
 //	dres.cont_l.sort();
 // output_list_New(dres.cont_l);

	toc = dres.toc;
    NTr = dres.num_contact;
	nItrs = dres.numCA;	   

	
	return TOCFound;
	
}




































































































































































































































































































































































































































































































