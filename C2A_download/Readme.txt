C2A is a continuous collision detection (CCD) library written in C++ for rigid polygon soup models.
The code is written using Visual C++ 8.0 (or Intel c++ compiler 11.0 for better performance)in MS Windows and is not tested on other platforms yet.
To run the code, VC++ 8.0 is necessary.
//...........................................................................................................................//
Part files of PQP and Blender are also needed. 

    Please download:
    PQP at http://www.cs.unc.edu/%7Egeom/SSV/index.html
    

1.Copy the src folder in PQP root path into C2A_download root path.

2.copy all the files inside src_C2A to the src folder.
//...........................................................................................................................//

In C2A.sln, two projects exist:

PQP_CCD: continous collision detection based on PQP soft package, and form PQP_CCD.lib
The main functions are:

PQP_REAL PQP_QueryTimeOfContact(CInterpMotion *objmotion1, CInterpMotion *objmotion2,
				PQP_TOCResult *res,  PQP_Model *o1, PQP_Model *o2,
				PQP_REAL tolerance_d, PQP_REAL tolerance_t, int qsize) 
//which function is used to execute conservative advancment algorithm

PQP_TOCStep(CInterpMotion *objmotion1, CInterpMotion *objmotion2, PQP_TOCResult *res,
			PQP_REAL R1[3][3], PQP_REAL T1[3], PQP_Model *o1,
			PQP_REAL R2[3][3], PQP_REAL T2[3], PQP_Model *o2,
			PQP_REAL tolerance,PQP_REAL tolerance_t,
			int qsize)
//which function is used to set parameter and call BVTT function for each CA iterative

double TOCStepRecurse_Dis(PQP_REAL tolerance_t, CInterpMotion *objmotion1, CInterpMotion *objmotion2,PQP_TOCResult *res, double R[3][3], double T[3], 
        PQP_Model *o1, int b1, PQP_Model *o2, int b2, double delta)
//Controlled CA function by BVTT



CCD_Demo: the Demo fucntion for showing how to use PQP_CCD.lib.
The result of C2A is saved in result folder under root path.
Macro GetContactNumber is defined in mainccd.h file to switch for getting contact features or not.

//...........................................................................................................................//

The project also referenced the follow open sources:
FAST: 
     http://graphics.ewha.ac.kr/FAST/
SWIFT++: 
     http://www.cs.unc.edu/~geom/SWIFT++/
Bullet Physics lib:  
     http://www.bulletphysics.org/

//...........................................................................................................................//

If you have any questions or find any bugs, please feel free to connect with Min Tang by EMail: tangmin@ewha.ac.kr.

