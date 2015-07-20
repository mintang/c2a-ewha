/******************************************************************
** Copyright (c) 2003 Computer Graphics & Imaging Laboratory
** Filename£º		Quaternion.h
** Creator£º		Zhang Xinyu
** Date£º			2005-02-23
** Description£º	Declare Quaternion class
**
** Version£º		1.0.1
**-----------------------------------------------------------------------------

**Revision history£º
**
** Reviser£º			date£º				Description£º
**-----------------------------------------------------------------------------
******************************************************************/

#ifndef QUATERNION_H
#define QUATERNION_H

#include <math.h>
#include "Coord.h"

#ifdef _DEBUG
#include <iostream>
using namespace std;
#endif

/** \nosubgrouping */

/** Defines a Quaternion.  */
class Quaternion  
{

public:
	/** The x coordinate */
    double x;
	/** The y coordinate */
	double y; 
	/** The z coordinate */
	double z;
	/** The w coordinate */
	double w;

public:
/************************************************/
/** \name Contruction and Deconstruction. A quaternion is defined by q=w+xi+yj+zk. */
//@{

	/** Default constructor. Initializes the quaternion to (0, 0, 0, 0). */
	Quaternion();

	/** Constructor. Initializes the quaternion with another quaternion. 
		\param q The given Quaternion	
	*/
	Quaternion(const Quaternion& q);

	/** Constructor. Initializes the quaternion with another quaternion. 
		\param x The x coordinate	
		\param y The y coordinate	
		\param z The z coordinate	
		\param w The w coordinate	
	*/
	Quaternion(const double x, const double y, const double z, const double w);

	/** Constructor. Initializes the quaternion with a rotation axis and angle. 
		\param axis		The rotation axis	
		\param angle	The rotation angle	
	*/
	Quaternion(Vector3D axis, double angle);

	Quaternion(const double& yaw, const double& pitch, const double& roll);

	/** Default deconstructor. */
	~Quaternion();

//@}

/************************************************/
/** \name Set/Get coordinates.*/
//@{

	/** Set the coordinate
		\param x The x coordinate	
		\param y The y coordinate
		\param z The z coordinate
		\param w The w coordinate
	*/
	void setCoords(double x = 0.0, double y = 0.0, double z = 0.0, double w=0.0);

	/** Set the coordinate
		\param q The given Quaternion	
	*/
	void setCoords(const Quaternion& q);

	/** Set the coordinate with rotation
		\param axis		The rotation axis	
		\param angle	The rotation angle	
	*/
	void setRotation(const Vector3D& axis, const double& angle);

	/** Set the coordinate with three Euler angles
		\param yaw		The first angle	
		\param pitch	The second angle		
		\param roll		The third angle	
	*/
	void setEuler(const double& yaw, const double& pitch, const double& roll);

	/** Copy \a q to current quaternion. 
		\param q The quaternion to copy
		\return A reference to this object
	*/
	Quaternion& operator =(const Quaternion &q);

	/** Adds \a q to current quaternion. 
		\param q The quaternion to add
		\return A reference to this object
	*/
	Quaternion& operator +=(const Quaternion &q);

	/** Subtract \a q from current quaternion. 
		\param coord The quaternion to subtract
		\return A reference to this object
	*/
	Quaternion& operator -=(const Quaternion &q);

	/** Adds \a q to current quaternion
		\param q The quaternion to add
		\return A new Coord object with the result
	*/
	Quaternion operator +(const Quaternion &q) const;

	/** Subtract \a q from current quaternion
		\param q The quaternion to subtract
		\return A new Coord object with the result
	*/
	Quaternion operator -(const Quaternion &q) const;

	/** Returns the current quaternion multiplied by -1
		\return A new quaternion object with the result
	*/
	Quaternion operator -() const;

	/** Compares 2 quaternions and returns true if they are the same
		\param q The quaternion to compare to
		\retval true The 2 quaternions are the same
		\retval false The quaternions are different
	*/
	bool operator ==(const Quaternion &q) const;

	/** Multiplies quaternion by \a d
		\param d The parameter to multiply by
		\return A reference to this object
	*/
	Quaternion& operator *=(double d);

	/** Multiplies quaternion by \a q
		\param q The parameter to multiply by
		\return A reference to this object
	*/
	Quaternion& operator*=(const Quaternion& q);

	/** Divides quaternion by \a d
		\param d The parameter to divide by
		\return A reference to this object
	*/
	Quaternion& operator /=(double d);

	/** Multiplies quaternion by \a d
		\param d The parameter to multiply by
		\return A new Coord object with the result
	*/
	Quaternion operator *(double d) const;

	/** Multiplies quaternion by \a q
		\param q The parameter to multiply by
		\return A new Coord object with the result
	*/
	Quaternion operator *(Quaternion q) const;

	/** Multiplies quaternion by \a q
		\param q The parameter to multiply by
		\return A new Coord object with the result
	*/
	friend Quaternion operator*(const Vector3D& v, const Quaternion& q);

	/** Divides quaternion by \a d
		\param d The parameter to divide by
		\return A new Coord object with the result
	*/
	Quaternion operator /(double d) const;

	/** Returns the x/y/z/w coordinate (depends on \a i). 
		\param i Specifies the x/y/z/w coordinate (0 for x, 1 for y, 2 for z, 3 for w)
		\return The current x/y/z/w coordinate
	*/
	double& operator [](int i) ;

	/** Returns the x/y/z/w coordinate (depends on \a i). 
		\param i Specifies the x/y/z/w coordinate (0 for x, 1 for y, 2 for z, 3 for w)
		\return The current x/y/z/w coordinate
	*/
	double operator [](int i) const;


//@}

/************************************************/
/** \name Operation .*/
//@{

	/** Returns the square of the norm of the quaternion.
		It is also the dot product between the quaternion and itself 
		\return The norm of the quaternion
	*/
	double sqabs() const;

	/** Returns the norm of the quaternion.
		It is also the dot product between the quaternion and itself 
		\return The square of the norm of the quaternion
	*/
	double abs() const;

	/** Returns a normalized vector (same direction as original vector, but with length=1) 
		\return Normalized vector
	*/
	Quaternion unit() const;

	/** Normalizes the vector to be of length=1 */
    void normalize();

	/** Dot multiply a quaternion
 		\param q The quaternion to multiply
	*/
    double dot(const Quaternion& q) const;

	/** Return the inverse quaternion
	*/
    Quaternion inverse() const;

	/** Returns the angle between the current quaternion and the given quaternion.
 		\param q The given quaternion
	*/
	double angle(const Quaternion& q) const;

	/** Returns the angle between the current quaternion and the given quaternion.
 		\param q The given quaternion
	*/
	double angle() const;
	
	/** Get the farthest quaternion
 		\param q The given quaternion to compare
	*/
	Quaternion farthest( const Quaternion& q) const;

	/** Get the slerp of the quaternion
 		\param q The given quaternion to compare
 		\param t The given parameter
	*/
	Quaternion slerp(const Quaternion& q, const double& t) const;
//@}


};


#endif
