/******************************************************************
** Copyright (c) 2003 Computer Graphics & Imaging Laboratory
** Filename£º		Coord.cpp
** Creator£º		Zhang Xinyu
** Date£º			2003-09-16
** Description£º	Define Coord3D and Coord2D classes
**
** Version£º		1.0.1
**-----------------------------------------------------------------------------

**Revision history£º
**
** Reviser£º			date£º				Description£º
**-----------------------------------------------------------------------------
******************************************************************/
//Refer to meshmaker 5.0
 
#include "Coord.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif


/** Default constructor. Initializes the coordinate to (0, 0, 0). */
Coord3D::Coord3D()
{
	x=0.0;y=0.0;z=0.0;
}

/** Constructor. Initializes the coordinate to (x, y, z). 
	\param x The x coordinate	
	\param y The y coordinate
	\param z The z coordinate
*/
Coord3D::Coord3D(double x, double y , double z )
{
	this->x=x;this->y=y;this->z=z;
}

/** Set the coordinate.
	\param x The x coordinate	
	\param y The y coordinate
	\param z The z coordinate
*/
void Coord3D::setCoords(double x , double y , double z )
{
	this->x=x; this->y=y; this->z=z;
}

/** Set the coordinate
	\param coord The coordinate	
*/
void Coord3D::setCoords(const Coord3D& coord)
{
	x=coord.x; y=coord.y; z=coord.z;
}

/** Set the coordinate to (0.0, 0.0, 0.0)
*/
void Coord3D::zero()
{
	x=0.0;y=0.0;z=0.0;
}

/** Copy \a coord to current coordinate. 
	\param coord The coordinate to copy
	\return A reference to this object
*/
Coord3D& Coord3D::operator =(const Coord3D &coord)
{
	x=coord.x; y=coord.y; z=coord.z;
	return (*this);
}

/** Adds \a coord to current coordinate. 
	\param coord The coordinate to add
	\return A reference to this object
*/
Coord3D& Coord3D::operator +=(const Coord3D & coord)
{
	x+=coord.x; y+=coord.y; z+=coord.z;
	return *this;
}

/** Subtract \a coord from current coordinate. 
	\param coord The coordinate to subtract
	\return A reference to this object
*/
Coord3D& Coord3D::operator -=(const Coord3D & coord)
{
	x-=coord.x; y-=coord.y; z-=coord.z;
	return *this;
}

/** Adds \a coord to current coordinate
	\param coord The coordinate to add
	\return A new Coord3D object with the result
*/
Coord3D Coord3D::operator +(const Coord3D & coord) const
{
	Coord3D rtcoord(x+coord.x, y+coord.y, z+coord.z);
	return rtcoord;
}

/** Subtract \a coord from current coordinate
	\param coord The coordinate to subtract
	\return A new Coord3D object with the result
*/
Coord3D Coord3D::operator -(const Coord3D & coord) const
{
	Coord3D rtcoord(x-coord.x, y-coord.y, z-coord.z);
	return rtcoord;
}

/** Returns the current coordinate multiplied by -1
	\return A new Coord3D object with the result
*/
Coord3D Coord3D::operator -() const
{
	Coord3D rtcoord(-x, -y, -z);
	return rtcoord;
}

/** Compares 2 coordinates and returns true if they are the same
	\param coord The coordinate to compare to
	\retval true The 2 coordinates are the same
	\retval false The coordinates are different
*/
bool Coord3D::operator ==(const Coord3D & coord) const
{
	return (x==coord.x&&y==coord.y&&z==coord.z)? true:false;
}

/** Multiplies coordinate by \a d
	\param d The parameter to multiply by
	\return A reference to this object
*/
Coord3D& Coord3D::operator *=(double d)
{
	x*=d; y*=d; z*=d;
	return *this;
}

/** Divides coordinate by \a d
	\param d The parameter to divide by
	\return A reference to this object
*/
Coord3D& Coord3D::operator /=(double d)
{
	x/=d; y/=d; z/=d;
	return *this;

}

/** Multiplies coordinate by \a d
	\param d The parameter to multiply by
	\return A new Coord3D object with the result
*/
Coord3D Coord3D::operator *(double d) const
{
	Coord3D rtcoord(x*d, y*d, z*d);
	return rtcoord;
}

/** Multiplies 3 x 3 matrix by \a m
	\param m The matrix parameter to multiply by
	\return A new Coord3D object with the result
*/
//Coord3D Coord3D::operator *(const Matrix3x3& m) const
//{
//	return Coord3D(m.tdot(0, (*this)), m.tdot(1, (*this)), m.tdot(2, (*this)));
//}

/** Divides coordinate by \a d
	\param d The parameter to divide by
	\return A new Coord3D object with the result
*/
Coord3D Coord3D::operator /(double d) const
{
	Coord3D rtcoord(x/d, y/d, z/d);
	return rtcoord;
}

/** Returns the x/y/z coordinate (depends on \a i). 
	\param i Specifies the x/y/z coordinate (0 for x, 1 for y, 2 for z)
	\return The current x/y/z coordinate
*/
double& Coord3D::operator [](int i)
{
	switch(i)
	{
	case 0:return x;
	case 1:return y;
	case 2:return z;
	default: return z;
	}
}

/** Returns the x/y/z coordinate (depends on \a i). 
	\param i Specifies the x/y/z coordinate (0 for x, 1 for y, 2 for z)
	\return The current x/y/z coordinate
*/
double Coord3D::operator [](int i) const
{
	switch(i)
	{
	case 0:return x;
	case 1:return y;
	case 2:return z;
	default:
		return z;
	}
}

/** Returns the index of the longest coordinate.
	\return The index of the longest coordinate.
*/
int Coord3D::longestCoord() const
{
	/* compute and index to the largest component of D */
	
	int index=0;
	double max=fabs(x);
	double b=fabs(y);
	double c=fabs(z);
	if(b>max) max=b,index=1;
	if(c>max) max=c,index=2;

	return index;
}

/** Returns the dot product between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The dot product result
*/
double Coord3D::dot(const Coord3D& coord)
{
	return this->x*coord.x + this->y*coord.y+ this->z*coord.z;
}

/** Returns the cross product between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The cross product result
*/
Coord3D Coord3D::cross(const Coord3D& coord)
{
	Coord3D rtcoord(	this->y*coord.z-coord.y*this->z, 
					this->z*coord.x-coord.z*this->x, 
					this->x*coord.y-coord.x*this->y);
	return rtcoord;
}

/** Returns the angle between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The angle result.
*/
double Coord3D::vectorAngle(const Coord3D& coord)
{
	double result;
	double acosx=dot(coord)/(this->sqrabs()*coord.sqrabs());
	
	if( fabs(fabs(acosx)-1.0)  <1.0E-6)
	{
		if(acosx<0.0  )
			acosx=-1.0;
		else
			acosx=1.0;

	}

	result=acos(acosx);
	
	return result;
}

/** Returns the norm of the vector.
	It is also the dot product between the vector and itself 
	\return The norm of the vector
*/
double Coord3D::sqrabs() const
{
	return sqrt(x*x+y*y+z*z);
}

/** Returns the square of the norm of the vector.
	It is also the dot product between the vector and itself 
	\return The square of the norm of the vector
*/
double Coord3D::abs() const
{
	return fabs(x*x+y*y+z*z);
}

/** Returns a normalized vector (same direction as original vector, but with length=1) 
	\return Normalized vector
*/
Coord3D Coord3D::unit() const
{
	double absval=sqrt(x*x+y*y+z*z);
	if(absval>0.0)
		return Coord3D(x/absval,y/absval,z/absval);
	else
		return Coord3D(0.0, 0.0, 0.0 );
}

/** Normalizes the vector to be of length=1 */
void Coord3D::normalize()
{
	double absval=sqrt(x*x+y*y+z*z);
	if(absval> 0.0)
	{
		x/=absval;
		y/=absval;
		z/=absval;
	}
}

/** Returns a vector perpendicular to the current coord. 
	\return a vector perpendicular to the current coord.
*/
Coord3D Coord3D::perpendicular()
{
	Coord3D result;
	if (fabs(x) < fabs(y)) {
	  if (fabs(x) < fabs(z)) {
		 result.x = 1.0f-x*x;
		 result.y = -x*y;
		 result.z = -x*z;
		 return result;
	  }
	}
	else {
	  if (fabs(y) < fabs(z)) {
		 result.x  =-y*x;
		 result.y = 1.0f-y*y;
		 result.z = -y*z;
		 return result;
	  }
	}
	result.x = -z*x;
	result.y = -z*y;
	result.z = 1.0f-z*z;
	result.normalize();
	return result;
}

/** Returns the dot product between the 2 coordinates
	\param a First coordinate
	\param b Second coordinate
	\return The dot product result
*/
double dot(const Coord3D& a, const Coord3D& b)
{
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

/** Returns the cross product between the 2 coordinates
	\param a First coordinate
	\param b Second coordinate
	\return The cross product result
*/
Coord3D cross(const Coord3D& a, const Coord3D& b)
{
	Coord3D rtcoord(a.y*b.z-b.y*a.z, a.z*b.x-b.z*a.x, a.x*b.y-b.x*a.y);
	return rtcoord;
}

/** Returns the angle  between between the 2 coordinates.
	\param a First coordinate
	\param b Second coordinate
	\return The angle result
*/
double vectorAngle(const Coord3D& a, const Coord3D& b)
{
	float alen = (float)sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
	float blen = (float)sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
	Coord3D c(a.y*b.z - a.z*b.y, a.z*b.x -a.x*b.z, a.x*b.y - a.y*b.x);
	float clen = (float)sqrt( c.x*c.x + c.y*c.y + c.z*c.z);
	float d = clen / (alen * blen);
	if((a.x*b.x + a.y*b.y + a.z*b.z)>0.0)
		return asin(d);
    else 
        return PI - asin(d);
}




/** Default constructor. Initializes the coordinate to (0, 0). */
Coord2D::Coord2D()
{
	x=0.0; y=0.0;
}

	/** Constructor. Initializes the coordinate to (x, y). 
		\param x The x coordinate	
		\param y The y coordinate
	*/
Coord2D::Coord2D(double x, double y )
{
	this->x=x;this->y=y;
}

	/** Set the coordinate
		\param x The x coordinate	
		\param y The y coordinate
	*/
void Coord2D::setCoords(double x , double y )
{
	this->x=x; this->y=y;
}

	/** Set the coordinate
		\param coord The coordinate	
	*/
void Coord2D::setCoords(const Coord2D& coord)
{
	this->x=coord.x; this->y=coord.y;
}

/** Copy \a coord to current coordinate. 
	\param coord The coordinate to copy
	\return A reference to this object
*/
Coord2D& Coord2D::operator =(const Coord2D &coord)
{
	x=coord.x; y=coord.y;
	return (*this);
}

	/** Adds \a coord to current coordinate. 
		\param coord The coordinate to add
		\return A reference to this object
	*/
Coord2D& Coord2D::operator +=(const Coord2D & coord2d)
{
	x+=coord2d.x;y+=coord2d.y;
	return *this;
}

	/** Subtract \a coord from current coordinate. 
		\param coord The coordinate to subtract
		\return A reference to this object
	*/
Coord2D& Coord2D::operator -=(const Coord2D & coord2d)
{
	x-=coord2d.x;y-=coord2d.y;
	return *this;
}

	/** Adds \a coord to current coordinate
		\param coord The coordinate to add
		\return A new Coord2D object with the result
	*/
Coord2D Coord2D::operator +(const Coord2D & coord2d) const
{
	Coord2D rtcoord2d(x+coord2d.x, y+coord2d.y);
	return rtcoord2d;

}

	/** Subtract \a coord from current coordinate
		\param coord The coordinate to subtract
		\return A new Coord2D object with the result
	*/
Coord2D Coord2D::operator -(const Coord2D & coord2d) const
{
	Coord2D rtcoord2d(x-coord2d.x, y-coord2d.y);
	return rtcoord2d;
}

	/** Returns the current coordinate multiplied by -1
		\return A new Coord2D object with the result
	*/
Coord2D Coord2D::operator -() const
{
	Coord2D rtcoord2d(-x,-y);
	return rtcoord2d;
}

	/** Compares 2 coordinates and returns true if they are the same
		\param coord The coordinate to compare to
		\retval true The 2 coordinates are the same
		\retval false The coordinates are different
	*/
bool Coord2D::operator ==(const Coord2D & coord2d) const
{
	return x==coord2d.x&&y==coord2d.y;
}

	/** Multiplies coordinate by \a d
		\param d The parameter to multiply by
		\return A reference to this object
	*/
Coord2D& Coord2D::operator *=(double d)
{
	x*=d;y*=d;
	return *this;
}

	/** Divides coordinate by \a d
		\param d The parameter to divide by
		\return A reference to this object
	*/
Coord2D& Coord2D::operator /=(double d)
{
	x/=d;y/=d;
	return *this;
}

	/** Multiplies coordinate by \a d, returns a new Coord2D object
		\param d The parameter to multiply by
		\return A new Coord2D object with the result
	*/
Coord2D Coord2D::operator *(double d) const
{
	Coord2D rtcoord2d(x*d, y*d);
	return rtcoord2d;
}

	/** Divides coordinate by \a d, returns a new Coord2D object
		\param d The parameter to divide by
		\return A new Coord2D object with the result
	*/
Coord2D Coord2D::operator /(double d) const
{
	Coord2D rtcoord2d(x/d, y/d);
	return rtcoord2d;
}

	/** Returns the x/y coordinate (depends on \a i). 
		\param i Specifies the x/y coordinate (0 for x, 1 for y)
		\return The current x/y coordinate
	*/
double& Coord2D::operator [](int i) 
{
	switch(i)
	{
	case 0:return x;
	case 1:return y;
	default:return y;
	}
}

	/** Returns the x/y/z coordinate (depends on \a i). 
		\param i Specifies the x/y/z coordinate (0 for x, 1 for y)
		\return The current x/y coordinate
	*/
double Coord2D::operator [](int i) const
{
	switch(i)
	{
	case 0:return x;
	case 1:return y;
	default:return 0;
	}
}

/** Returns the index of the longest coordinate.
	\return The index of the longest coordinate.
*/
int Coord2D::longestCoord() const
{
	return (fabs(x)>fabs(y))?0:1;
}
/** Returns the norm of the vector.
	It is also the dot product between the vector and itself 
	\return The norm of the vector
*/
double Coord2D::sqrabs() const
{
	return sqrt(x*x+y*y);
}

/** Returns the square of the length of the vector.
	It is also the square root of the dot product between the vector and itself 
	\return The length of the vector
*/
double Coord2D::abs() const
{
	return fabs(x*x+y*y);
}

/** Returns a normalized vector (same direction as original vector, but with length=1) 
	\return Normalized vector
*/
Coord2D Coord2D::unit() const
{
	double absval=sqrt(x*x+y*y);
	Coord2D rtcoord(x/absval,y/absval);
	return rtcoord;
}

/** Normalizes the vector to be of length=1 */
void Coord2D::normalize()
{
	double absval=sqrt(x*x+y*y);
	x/=absval;y/=absval;
}


/** Returns the dot product between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The dot product result
*/
double Coord2D::dot(const Coord2D& coord)
{
	return this->x*coord.x + this->y*coord.y;
}

/** Returns the cross product between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The cross product result
*/
Coord3D Coord2D::cross(const Coord2D& coord)
{
	Coord3D rtcoord(	0, 	0, this->x*coord.y-coord.x*this->y);
	return rtcoord;
}

/** Returns the angle between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The angleresult.
*/
double Coord2D::vectorAngle2D(const Coord2D& coord)
{
	float alen = (float)sqrt(x*x + y*y);
	float blen = (float)sqrt(coord.x*coord.x + coord.y*coord.y);
	float d = (float)(fabs(x*coord.y - y*coord.x) / (alen * blen));
	if((x*coord.x + y*coord.y )>0.0)
		return asin(d);
    else 
        return PI - asin(d);
}


/** Returns the dot product between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The dot product result
*/
double dot(const Coord2D& a, const Coord2D& b)
{
	return a.x*b.x+a.y*b.y;
}

/** Returns the cross product between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The cross product result
*/
Coord3D cross(const Coord2D& a, const Coord2D& b)
{
	Coord3D rtcoord(	0, 	0, a.x*b.y-b.x*a.y);
	return rtcoord;
}

/** Returns the angle between the current coordinates and the given coordinates.
	\param coord The given coordinates.
	\return The angleresult.
*/
double vectorAngle2D(const Coord2D& a, const Coord2D& b)
{
	float alen = (float)sqrt(a.x*a.x + a.y*a.y);
	float blen = (float)sqrt(b.x*b.x + b.y*b.y);
	float d = (float)(fabs(a.x*b.y - a.y*b.x) / (alen * blen));
	if((a.x*b.x + a.y*b.y )>0.0)
		return asin(d);
    else 
        return PI - asin(d);
}

