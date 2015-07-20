/******************************************************************
** Copyright (c) 2003 Computer Graphics & Imaging Laboratory
** Filename£º		Transform.h
** Creator£º		Zhang Xinyu
** Date£º			2003-10-15
** Description£º	Declare Transform class
**
** Version£º		1.0.1
**-----------------------------------------------------------------------------

**Revision history£º
**
** Reviser£º			date£º				Description£º
**-----------------------------------------------------------------------------
******************************************************************/

#ifndef TRANSFORM_H
#define TRANSFORM_H


#include "def.h"
#include "Matrix3x3.h"

typedef enum { 
	IDENTITY	= 0x0,
	TRANSLATION = 0x01,
	ROTATION    = 0x02,
	RIGID       = TRANSLATION | ROTATION,  
	SCALING     = 0x04,
	LINEAR      = ROTATION | SCALING,
	AFFINE      = TRANSLATION | LINEAR
} XFormTYPE;

/** Defines a 4 x 4 transform matrix.  */
class Transform{

public:
	/** Default constructor. Initializes the coordinate to identity. */
	Transform();
	
	/** Constructor. Initializes the coordinate. */
	Transform(const Matrix3x3& b, const Vector3D& c = Vector3D(0.0, 0.0, 0.0), XFormTYPE type = AFFINE);

	/** Constructor. Initializes the coordinate. */
	Transform(const Quaternion& q, const Vector3D& c = Vector3D(0.0, 0.0, 0.0), XFormTYPE type = RIGID);

	/** Default deconstructor.*/
	~Transform();


	/** Sets the current Xform according to the given one 
		\param t The Xform to copy
		\return A reference to this object
	*/
	Transform& operator =(Transform& t);

	/** Returns the multiply of current Xform and \a Xform
		\param t The Xform to multiply
		\return A new Xform object with the result
	*/
	Transform operator *(Transform& t);

	/** Returns sum of the given vector and the multiply of current Xform and the \a vector 
		\param x The vector to multiply
		\return A new vector object with the result
	*/
	Vector3D operator *(Vector3D& x);

	/** Returns sum of the given vector and the multiply of current Xform and the \a vector 
		\param x The vector to multiply
		\return A new vector object with the result
	*/
	Vector3D operator()(Vector3D& x);

	/** Returns the multiply result (1/X1)X2
		\param t The transform to multiply
		\return A new transform object with the result
	*/
	Transform multInverseLeft(Transform& t);

	/** Return the inverse  of the current transform.
	*/
	Transform inverse();

	/** Set the origin of the current transform.
		\param o	The given origin	
	*/
	void setOrigin(Coord3D o);

	/** Set the origin of the current transform.
		\param x	The given origin x	
		\param y	The given origin y	
		\param z	The given origin z	
	*/
	void setOrigin(double x, double y, double z);

	/** Return the origin of the current transform.
	*/
	Coord3D getOrigin(void);

	/** Return the origin of the current transform.
	*/
	void getOrigin(double* o);

	/** Set the rotation matrix of the current transform.
		\param m	The given rotation matrix	
	*/
	void setRotation(Quaternion q);

	/** Return the rotation matrix of the current transform.
	*/
	Quaternion getRotation(void);

	/** Set the rotation matrix of the current transform.
		\param m	The given rotation matrix	
	*/
	void setBasis(Matrix3x3 m);

	/** Return the rotation matrix of the current transform.
	*/
	Matrix3x3 getBasis(void);

	/** Return the rotation matrix of the current transform.
	*/
	void getBasis(double *b);

	/** Set the current transform identity.
	*/
	void setIdentity(void);

	/** Check whether the current transform is identity.
	*/
	bool isIdentity(void);

	/** Set the elements with OpenGL matrix
		\param m The given OpenGL matrix	
	*/
	void setFromOpenGLMatrix(const double *m);

	/** Return the elements in the order of OpenGL matrix
		\param m The return OpenGL matrix	
	*/
	void getOpenGLMatrix(double *m);

	/** Scale the Xform with the given vector
		\param s The given vector	
	*/
	void scale(const Vector3D& s) const;

	/** Integrate a transform using the given linear velocity and angular velocity.
		\param v	The given linear velocity	
		\param w	The given angular velocity
		\param dt	The given time interval
	*/
	Transform integrate(const Vector3D* v, const Vector3D* w, const double dt);

	/** Compute the velocities from the current and the other given transform.
		\param xform	The given transform	
		\param dt		The given time interval
		\param v		The linear velocity returned	
		\param w		The angular velocity returned
	*/
	void velocity(Transform* xform, double dt, Vector3D& v, Vector3D& w);


public:
	Matrix3x3 basis;	/* The rotation components */
	Vector3D  origin;	/* The translation components */
	XFormTYPE type;
} ;




#endif
