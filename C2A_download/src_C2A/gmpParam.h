#ifndef GMP_PARAM_H
#define GMP_PARAM_H

// Setting parameter in GMP
// Dec 5, 2007

// CCD should return a contact with distance less than eps and larger than RATIO*eps
#define GMP_CCD_SECURITY_DISTANCE_RATIO 0.1

#define GMP_CCD_DELTA_D 0.0001
#define GMP_CCD_DELTA_T 0.0001

// #define GMP_CCD_DELTA_D 0.1
// #define GMP_CCD_DELTA_T 0.1
// #define CCD_AT_MOST_EPSILON 0.1


// From GPDCompCD.h
#define CCD_AT_MOST_EPSILON 0.001
#define CCD_AT_LEAST_EPSILON 0.00005
#define  CD_TOLERANCE  0 // CCD_AT_LEAST_EPSILON /2 

extern int SWIFT_CONTACT_TOLEARNCE_BUG;



#define CONTACT_TOL_3D 5e-3

#define DISCREETE_CCD 100




#endif
