#ifndef GMP_COMMON_H
#define GMP_COMMON_H

enum GMP_INTERP_MODE
{
	GMP_IM_EULER   = 0,   // linear interpolation of the euler angles.
	GMP_IM_LINEAR  = 1,   // constant rotation velocity and constant translation velocity of the origin.
	GMP_IM_SCREW   = 2,   // screw motion
	GMP_IM_SLERP   = 3,   // SLERP motion
	GMP_IM_LINEAR_ONE_CONTACT = 4, // Constained motion (v constaint) 
	GMP_IM_LINEAR_TWO_CONTACTS = 5,
	GMP_IM_LINEAR_THR_CONTACTS = 6,
	GMP_IM_ALLMODE = 7,
};


extern char *GMP_INTERP_MODE_Desc[GMP_IM_ALLMODE]; 

void GMP_ASSERT(const char *str);
void GMP_WARNING(const char *str);
bool IsIMCMI(GMP_INTERP_MODE imMode); // check whether the imMode is CMI

#define GMP_IM_METHOD_DEAFULT GMP_IM_LINEAR

#endif
