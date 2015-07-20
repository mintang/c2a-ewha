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
#ifndef PQP_CCD_H
#define PQP_CCD_H

class CInterpMotion ;


//CInterpMotion &objmotion1, CInterpMotion &objmotion2,
PQP_REAL PQP_QueryTimeOfContact(CInterpMotion *objmotion1, CInterpMotion *objmotion2,
								PQP_TOCResult *res, PQP_Model *o1, PQP_Model *o2,
								PQP_REAL tolerance_d, PQP_REAL tolerance_t, int qsize=2);
PQP_REAL PQP_QueryContact(CInterpMotion *objmotion1, CInterpMotion *objmotion2,
								PQP_TOCResult *res, PQP_Model *o1, PQP_Model *o2,
								double threshold);

int
PQP_TOCStep(CInterpMotion *objmotion1, CInterpMotion *objmotion2, PQP_TOCResult *res, 
			PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
			PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
			PQP_REAL tolerance,PQP_REAL tolerance_t,
			int qsize = 2);

 

#endif