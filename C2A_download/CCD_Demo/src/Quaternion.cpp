/******************************************************************
** Copyright (c) 2003 Computer Graphics & Imaging Laboratory
** Filename£º		Quaternion.cpp
** Creator£º		Zhang Xinyu
** Date£º			2005-02-23
** Description£º	Define Quaternion class
**
** Version£º		1.0.1
**-----------------------------------------------------------------------------

**Revision history£º
**
** Reviser£º			date£º				Description£º
**-----------------------------------------------------------------------------
******************************************************************/

 
#include "Quaternion.h"
#include "assert.h"


/** Default constructor. Initializes the quaternion to (0, 0, 0, 0). */
Quaternion::Quaternion()
{
	x=0; y=0; z=0; w=0; 
}

/** Constructor. Initializes the quaternion with another quaternion. 
	\param q The given Quaternion	
*/
Quaternion::Quaternion(const Quaternion& q)
{
	x=q.x; y=q.y; z=q.z; w=q.w; 
}

/** Constructor. Initializes the quaternion with another quaternion. 
	\param x The x coordinate	
	\param y The y coordinate	
	\param z The z coordinate	
	\param w The w coordinate	
*/
Quaternion::Quaternion(const double x, const double y, const double z, const double w)
{
	this->x=x; this->y=y; this->z=z; this->w=w; 
}

/** Constructor. Initializes the quaternion with a rotation axis and angle. 
	\param axis		The rotation axis	
	\param angle	The rotation angle	
*/
Quaternion::Quaternion(Vector3D axis, double angle)
{
	setRotation(axis, angle);
}

/** Constructor. Initializes the quaternion with three rotation angles. 
	\param yaw		The first angle	
	\param pitch	The second angle		
	\param roll		The third angle	
*/
Quaternion::Quaternion(const double& yaw, const double& pitch, const double& roll)
{ 
	setEuler(yaw, pitch, roll); 
}

/** Default deconstructor. */
Quaternion::~Quaternion()
{

}

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
void Quaternion::setCoords(double x , double y , double z , double w)
{
	this->x=x; this->y=y; this->z=z; this->w=w; 
}

/** Set the coordinate
	\param q The given Quaternion	
*/
void Quaternion::setCoords(const Quaternion& q)
{
	x=q.x; y=q.y; z=q.z; w=q.w; 
}

/** Set the coordinate with rotation
	\param axis		The rotation axis	
	\param angle	The rotation angle	
*/
void Quaternion::setRotation(const Vector3D& axis, const double& angle)
{
	double d = axis.sqrabs();
	assert(d != double(0.0));
	double s = sin(angle * double(0.5)) / d;
	setCoords(axis[0] * s, axis[1] * s, axis[2] * s, cos(angle * double(0.5)));
}

/** Set the coordinate with three Euler angles
	\param yaw		The first angle	
	\param pitch	The second angle		
	\param roll		The third angle	
*/
void Quaternion::setEuler(const double& yaw, const double& pitch, const double& roll)
{
	double halfYaw	 =	double(yaw) * 0.5;  
	double halfPitch =	double(pitch) * 0.5;  
	double halfRoll  =	double(roll) * 0.5;  
	double cosYaw	 =	cos(halfYaw);
	double sinYaw	 =	sin(halfYaw);
	double cosPitch	 =	cos(halfPitch);
	double sinPitch  =	sin(halfPitch);
	double cosRoll   =	cos(halfRoll);
	double sinRoll   =	sin(halfRoll);

	setCoords(cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw,
		cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw,
		sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw,
		cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw);
}

/** Copy \a q to current quaternion. 
	\param q The quaternion to copy
	\return A reference to this object
*/
Quaternion& Quaternion::operator =(const Quaternion &q)
{
	x=q.x;	y=q.y;	z=q.z; w=q.w;

	return (*this);
}

/** Adds \a q to current quaternion. 
	\param q The quaternion to add
	\return A reference to this object
*/
Quaternion& Quaternion::operator +=(const Quaternion &q)
{
	x+=q.x;	y+=q.y;	z+=q.z; w+=q.w;
	return (*this);
}

/** Subtract \a q from current quaternion. 
	\param coord The quaternion to subtract
	\return A reference to this object
*/
Quaternion& Quaternion::operator -=(const Quaternion &q)
{
	x-=q.x;	y-=q.y;	z-=q.z; w-=q.w;
	return (*this);
}

/** Multiplies quaternion by \a d
	\param d The parameter to multiply by
	\return A reference to this object
*/
Quaternion& Quaternion::operator *=(double d)
{
	x*=d;	y*=d;	z*=d; w*=d;
	return (*this);
}

/** Multiplies quaternion by \a q
	\param q The parameter to multiply by
	\return A reference to this object
*/
Quaternion& Quaternion::operator*=(const Quaternion& q)
{
	setCoords(	w * q.x + x * q.w + y * q.z - z * q.y,
				w * q.y + y * q.w + z * q.x - x * q.z,
				w * q.z + z * q.w + x * q.y - y * q.x,
				w * q.w - x * q.x - y * q.y - z * q.z);
	return (*this);
}

/** Divides quaternion by \a d
	\param d The parameter to divide by
	\return A reference to this object
*/
Quaternion& Quaternion::operator /=(double d)
{
	x/=d;	y/=d;	z/=d; w/=d;
	return (*this);
}

/** Adds \a q to current quaternion
	\param q The quaternion to add
	\return A new Coord object with the result
*/
Quaternion Quaternion::operator +(const Quaternion &q) const
{
	return Quaternion(x+q.x, y+q.y, z+q.z, w+q.w);
}

/** Subtract \a q from current quaternion
	\param q The quaternion to subtract
	\return A new Coord object with the result
*/
Quaternion Quaternion::operator -(const Quaternion &q) const
{
	return Quaternion(x-q.x, y-q.y, z-q.z, w-q.w);
}

/** Multiplies quaternion by \a d
	\param d The parameter to multiply by
	\return A new Coord object with the result
*/
Quaternion Quaternion::operator *(double d) const
{
	return Quaternion(x*d, y*d, z*d, w*d);
}

/** Multiplies quaternion by \a q
	\param q The parameter to multiply by
	\return A new Coord object with the result
*/
Quaternion Quaternion::operator *(Quaternion q) const
{
	return Quaternion( 
		w * q[0] + x * q[3] + y * q[2] - z * q[1],
		w * q[1] + y * q[3] + z * q[0] - x * q[2],
		w * q[2] + z * q[3] + x * q[1] - y * q[0],
		w * q[3] - x * q[0] - y * q[1] - z * q[2]);
}


/** Multiplies quaternion by \a q
	\param q The parameter to multiply by
	\return A new Coord object with the result
*/
Quaternion operator*(const Vector3D& v, const Quaternion& q)
{
	return Quaternion(  q[3] * v.x	+ q[2] * v.y  - q[1] * v.z,
						q[3] * v.y	+ q[0] * v.z  - q[2] * v.x,
						q[3] * v.z	+ q[1] * v.x  - q[0] * v.y,
					   -q[0]* v.x	- q[1] * v.y  - q[2] * v.z); 
}

/** Divides quaternion by \a d
	\param d The parameter to divide by
	\return A new Coord object with the result
*/
Quaternion Quaternion::operator /(double d) const
{
	return Quaternion(x/d, y/d, z/d, w/d);
}
/** Returns the current quaternion multiplied by -1
	\return A new quaternion object with the result
*/
Quaternion Quaternion::operator -() const
{
	return Quaternion(-x, -y, -z, -w);
}

/** Compares 2 quaternions and returns true if they are the same
	\param q The quaternion to compare to
	\retval true The 2 quaternions are the same
	\retval false The quaternions are different
*/
bool Quaternion::operator ==(const Quaternion &q) const
{
	return (x==q.x&&y==q.y&&z==q.z&&w==q.w)? true:false;
}

/** Returns the x/y/z/w coordinate (depends on \a i). 
	\param i Specifies the x/y/z/w coordinate (0 for x, 1 for y, 2 for z, 3 for w)
	\return The current x/y/z/w coordinate
*/
double& Quaternion::operator [](int i)
{
	switch(i)
	{
	case 0:return x;
	case 1:return y;
	case 2:return z;
	case 3:return w;
	default:
		return z;
	}
}

/** Returns the x/y/z/w coordinate (depends on \a i). 
	\param i Specifies the x/y/z/w coordinate (0 for x, 1 for y, 2 for z, 3 for w)
	\return The current x/y/z/w coordinate
*/
double Quaternion::operator [](int i) const
{
	switch(i)
	{
	case 0:return x;
	case 1:return y;
	case 2:return z;
	case 3:return w;
	default:
		return z;
	}
}


//@}

/************************************************/
/** \name Operation .*/
//@{

/** Returns the square of the norm of the quaternion.
	It is also the dot product between the quaternion and itself 
	\return The norm of the quaternion
*/
double Quaternion::sqabs() const
{
	return x*x+y*y+z*z+w*w;
}

/** Returns the square of the norm of the quaternion.
	It is also the dot product between the quaternion and itself 
	\return The square of the norm of the quaternion
*/
double Quaternion::abs() const
{
	return sqrt(x*x+y*y+z*z+w*w);
}

/** Returns a normalized vector (same direction as original vector, but with length=1) 
	\return Normalized vector
*/
Quaternion Quaternion::unit() const
{
	double absval=sqrt(x*x+y*y+z*z+w*w);
	if(absval>0.0)
		return Quaternion(x/absval,y/absval,z/absval,w/absval);
	else
		return Quaternion(0.0, 0.0, 0.0, 0.0 );
}

/** Normalizes the vector to be of length=1 */
void Quaternion::normalize()
{
	double absval=sqrt(x*x+y*y+z*z);
	if(absval> 0.0)
	{
		x/=absval;
		y/=absval;
		z/=absval;
		w/=absval;
	}
}

/** Dot multiply a quaternion
 	\param q The quaternion to multiply
*/
double Quaternion::dot(const Quaternion& q) const
{
	return x*q.x+y*q.y+z*q.z+w*q.w;
}

/** Return the inverse quaternion
*/
Quaternion Quaternion::inverse() const
{
	return Quaternion(x, y, z, -w);
}

/** Returns the angle between the current quaternion and the given quaternion.
 	\param q The given quaternion
*/
double Quaternion::angle(const Quaternion& q) const 
{
	double s = sqrt(sqabs() * q.sqabs());
	assert(s != double(0.0));
	return acos(dot(q) / s);
}

/** Returns the angle .
*/
double Quaternion::angle() const 
{
	return 2.0 * acos(w);
}

/** Get the farthest quaternion
 	\param q The given quaternion to compare
*/
Quaternion Quaternion::farthest( const Quaternion& q) const 
{
	Quaternion diff,sum;
	diff = *this - q;
	sum = *this + q;
	if( diff.dot(diff) > sum.dot(sum) )
		return q;
	return (-q);
}

/** Get the slerp of the quaternion
 	\param q The given quaternion to compare
 	\param t The given parameter
*/
Quaternion Quaternion::slerp(const Quaternion& q, const double& t) const
{
	double theta = angle(q);
	if (theta != double(0.0))
	{
		double d = double(1.0) / sin(theta);
		double s0 = sin((double(1.0) - t) * theta);
		double s1 = sin(t * theta);   
		return Quaternion((x * s0 + q.x * s1) * d,
			(y * s0 + q.y * s1) * d,
			(z * s0 + q.z * s1) * d,
			(w * s0 + q.w * s1) * d);
	}
	else
	{
		return *this;
	}
}

//@}



