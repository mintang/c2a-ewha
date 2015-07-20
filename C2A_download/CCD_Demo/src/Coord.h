/******************************************************************
** Copyright (c) 2003 Computer Graphics & Imaging Laboratory
** Filename£º		Coord.h
** Creator£º		Zhang Xinyu
** Date£º			2003-09-16
** Description£º	Declare Coord3D and Coord2D classes
**
** Version£º		1.0.1
**-----------------------------------------------------------------------------

**Revision history£º
**
** Reviser£º			date£º				Description£º
**-----------------------------------------------------------------------------
******************************************************************/
//Refer to meshmaker 5.0

#ifndef COORD_H
#define COORD_H

#include <math.h>
#include "def.h"

inline double Degree2Radian(double degree)
{
	return (degree / 180.0 * PI);
}

inline double Radian2Degree(double radian)
{
	return (radian / PI * 180.0);
}


/** Defines a 3D coordinate. Every component (x, y, and z) is a double. 
	The coordinate can also be used as a vector, starting at (0,0,0) and ending at (x,y,z). */
class Coord3D  
{
public:
	/** The x coordinate */
    double x;
	/** The y coordinate */
	double y; 
	/** The z coordinate */
	double z;

public:
	/** Default constructor. Initializes the coordinate to (0, 0, 0). */
	Coord3D();
	/** Constructor. Initializes the coordinate to (x, y, z). 
		\param x The x coordinate	
		\param y The y coordinate
		\param z The z coordinate
	*/
	Coord3D(double x, double y = 0.0, double z = 0.0);
	/** Set the coordinate
		\param x The x coordinate	
		\param y The y coordinate
		\param z The z coordinate
	*/
	void setCoords(double x = 0.0, double y = 0.0, double z = 0.0);

	/** Set the coordinate
		\param coord The coordinate	
	*/
	void setCoords(const Coord3D& coord);

	/** Set the coordinate to (0.0, 0.0, 0.0)
	*/
	void zero();

	/** Copy \a coord to current coordinate. 
		\param coord The coordinate to copy
		\return A reference to this object
	*/
	Coord3D& operator =(const Coord3D &coord);
	/** Adds \a coord to current coordinate. 
		\param coord The coordinate to add
		\return A reference to this object
	*/
	Coord3D& operator +=(const Coord3D &coord);
	/** Subtract \a coord from current coordinate. 
		\param coord The coordinate to subtract
		\return A reference to this object
	*/
	Coord3D& operator -=(const Coord3D &coord);
	/** Adds \a coord to current coordinate
		\param coord The coordinate to add
		\return A new Coord3D object with the result
	*/
	Coord3D operator +(const Coord3D &coord) const;
	/** Subtract \a coord from current coordinate
		\param coord The coordinate to subtract
		\return A new Coord3D object with the result
	*/
	Coord3D operator -(const Coord3D &coord) const;
	/** Returns the current coordinate multiplied by -1
		\return A new Coord3D object with the result
	*/
	Coord3D operator -() const;

	/** Compares 2 coordinates and returns true if they are the same
		\param coord The coordinate to compare to
		\retval true The 2 coordinates are the same
		\retval false The coordinates are different
	*/
	bool operator ==(const Coord3D &coord) const;

	/** Multiplies coordinate by \a d
		\param d The parameter to multiply by
		\return A reference to this object
	*/
	Coord3D& operator *=(double d);
	/** Divides coordinate by \a d
		\param d The parameter to divide by
		\return A reference to this object
	*/
	Coord3D& operator /=(double d);
	/** Multiplies coordinate by \a d
		\param d The parameter to multiply by
		\return A new Coord3D object with the result
	*/
	Coord3D operator *(double d) const;
	/** Multiplies 3 x 3 matrix by \a m
		\param m The matrix parameter to multiply by
		\return A new Coord3D object with the result
	*/
	//Coord3D operator *(const Matrix3x3& m) const;
	/** Divides coordinate by \a d
		\param d The parameter to divide by
		\return A new Coord3D object with the result
	*/
	Coord3D operator /(double d) const;

	/** Returns the x/y/z coordinate (depends on \a i). 
		\param i Specifies the x/y/z coordinate (0 for x, 1 for y, 2 for z)
		\return The current x/y/z coordinate
	*/
	double& operator [](int i) ;
	/** Returns the x/y/z coordinate (depends on \a i). 
		\param i Specifies the x/y/z coordinate (0 for x, 1 for y, 2 for z)
		\return The current x/y/z coordinate
	*/
	double operator [](int i) const;

	/** Returns the index of the longest coordinate.
		\return The index of the longest coordinate.
	*/
	int longestCoord() const;

	/** Returns the dot product between the current coordinates and the given coordinates.
		\param coord The given coordinates.
		\return The dot product result
	*/
	double dot(const Coord3D& coord);

	/** Returns the cross product between the current coordinates and the given coordinates.
		\param coord The given coordinates.
		\return The cross product result
	*/
	Coord3D cross(const Coord3D& coord);

	/** Returns the angle between the current coordinates and the given coordinates.
		\param coord The given coordinates.
		\return The angle result.
	*/
	double vectorAngle(const Coord3D& coord);

	/** Returns the norm of the vector.
		It is also the dot product between the vector and itself 
		\return The norm of the vector
	*/
	double sqrabs() const;

	/** Returns the square of the norm of the vector.
		It is also the dot product between the vector and itself 
		\return The square of the norm of the vector
	*/
	double abs() const;

	/** Returns a normalized vector (same direction as original vector, but with length=1) 
		\return Normalized vector
	*/
	Coord3D unit() const;

	/** Normalizes the vector to be of length=1 */
    void normalize();


	/** Returns a vector perpendicular to the current coord. 
		\return a vector perpendicular to the current coord.
	*/
	Coord3D perpendicular();

	/** Returns the dot product between the 2 coordinates
		\param a First coordinate
		\param b Second coordinate
		\return The dot product result
	*/
	friend double dot(const Coord3D& a, const Coord3D& b);

	/** Returns the cross product between the 2 coordinates
		\param a First coordinate
		\param b Second coordinate
		\return The cross product result
	*/
	friend Coord3D cross(const Coord3D& a, const Coord3D& b);

	/** Returns the angle  between between the 2 coordinates.
		\param a First coordinate
		\param b Second coordinate
		\return The angle result
	*/
	friend double vectorAngle(const Coord3D& a, const Coord3D& b);
};

/** Defines a 2D coordinate. Every component (x and y) is a double. 
	The coordinate can also be used as a vector, starting at (0,0) and ending at (x,y). */

class Coord2D 
{ 
public:
	/** The x coordinate */
    double x;
	/** The y coordinate */
	double y; 

public:
	/** Default constructor. Initializes the coordinate to (0, 0). */
	Coord2D();
	/** Constructor. Initializes the coordinate to (x, y). 
		\param x The x coordinate	
		\param y The y coordinate
	*/
	Coord2D(double x, double y = 0.0);
	/** Set the coordinate
		\param x The x coordinate	
		\param y The y coordinate
	*/
	void setCoords(double x = 0.0, double y = 0.0);

	/** Set the coordinate
		\param coord The coordinate	
	*/
	void setCoords(const Coord2D& coord);

	/** Copy \a coord to current coordinate. 
		\param coord The coordinate to copy
		\return A reference to this object
	*/
	Coord2D& operator =(const Coord2D &coord);
	/** Adds \a coord to current coordinate. 
		\param coord The coordinate to add
		\return A reference to this object
	*/
	Coord2D& operator +=(const Coord2D &coord);
	/** Subtract \a coord from current coordinate. 
		\param coord The coordinate to subtract
		\return A reference to this object
	*/
	Coord2D& operator -=(const Coord2D &coord);
	/** Adds \a coord to current coordinate
		\param coord The coordinate to add
		\return A new Coord2D object with the result
	*/
	Coord2D operator +(const Coord2D &coord) const;
	/** Subtract \a coord from current coordinate
		\param coord The coordinate to subtract
		\return A new Coord2D object with the result
	*/
	Coord2D operator -(const Coord2D &coord) const;
	/** Returns the current coordinate multiplied by -1
		\return A new Coord2D object with the result
	*/
	Coord2D operator -() const;
	/** Compares 2 coordinates and returns true if they are the same
		\param coord The coordinate to compare to
		\retval true The 2 coordinates are the same
		\retval false The coordinates are different
	*/
	bool operator ==(const Coord2D &coord) const;
	/** Multiplies coordinate by \a d
		\param d The parameter to multiply by
		\return A reference to this object
	*/
	Coord2D& operator *=(double d);
	/** Divides coordinate by \a d
		\param d The parameter to divide by
		\return A reference to this object
	*/
	Coord2D& operator /=(double d);
	/** Multiplies coordinate by \a d, returns a new Coord2D object
		\param d The parameter to multiply by
		\return A new Coord2D object with the result
	*/
	Coord2D operator *(double d) const;
	/** Divides coordinate by \a d, returns a new Coord2D object
		\param d The parameter to divide by
		\return A new Coord2D object with the result
	*/
	Coord2D operator /(double d) const;

	/** Returns the x/y coordinate (depends on \a i). 
		\param i Specifies the x/y coordinate (0 for x, 1 for y)
		\return The current x/y coordinate
	*/
	double& operator [](int i) ;
	/** Returns the x/y/z coordinate (depends on \a i). 
		\param i Specifies the x/y/z coordinate (0 for x, 1 for y)
		\return The current x/y coordinate
	*/
	double operator [](int i) const;

	/** Returns the index of the longest coordinate.
		\return The index of the longest coordinate.
	*/
	int longestCoord() const;

	/** Returns the norm of the vector.
		It is also the dot product between the vector and itself 
		\return The norm of the vector
	*/
	double sqrabs() const;
	/** Returns the square of the length of the vector.
		It is also the square root of the dot product between the vector and itself 
		\return The length of the vector
	*/
	double abs() const;
	/** Returns a normalized vector (same direction as original vector, but with length=1) 
		\return Normalized vector
	*/
	Coord2D unit() const;

	/** Normalizes the vector to be of length=1 */
    void normalize();

	
	/** Returns the dot product between the current coordinates and the given coordinates.
		\param coord The given coordinates.
		\return The dot product result
	*/
	double dot(const Coord2D& coord);

	/** Returns the cross product between the current coordinates and the given coordinates.
		\param coord The given coordinates.
		\return The cross product result
	*/
	Coord3D cross(const Coord2D& coord);
	
	/** Returns the angle between the current coordinates and the given coordinates.
		\param coord The given coordinates.
		\return The angleresult.
	*/
	double vectorAngle2D(const Coord2D& coord);

	/** Returns the dot product between the 2 coordinates
		\param a First coordinate
		\param b Second coordinate
		\return The dot product result
	*/
	friend double dot(const Coord2D& a, const Coord2D& b);

	/** Returns the cross product between the 2 coordinates
		\param a First coordinate
		\param b Second coordinate
		\return The cross product result
	*/
	friend Coord3D cross(const Coord2D& a, const Coord2D& b);

	/** Returns the angle  between between the 2 coordinates.
		\param a First coordinate
		\param b Second coordinate
		\return The angle result
	*/
	friend double vectorAngle2D(const Coord2D& a, const Coord2D& b);
  
};


typedef Coord3D Position;
typedef Coord3D Vector3D;


#endif
