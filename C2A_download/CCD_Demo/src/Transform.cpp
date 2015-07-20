/******************************************************************
** Copyright (c) 2003 Computer Graphics & Imaging Laboratory
** Filename£º		Transform.cpp
** Creator£º		Zhang Xinyu
** Date£º			2003-10-15
** Description£º	Define Transform class
**
** Version£º		1.0.1
**-----------------------------------------------------------------------------

**Revision history£º
**
** Reviser£º			date£º				Description£º
**-----------------------------------------------------------------------------
******************************************************************/



#include "Transform.h"
#include "math.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define ANGULARMOTIONTRESHOLD 0.78539815 //PI/4

const double XFORM_EPSILON=1e-15;

Transform::Transform()
{
	basis.setIdentity();
	origin.zero();
	type=IDENTITY;
}

Transform::Transform(const Matrix3x3& b, const Vector3D& c, XFormTYPE type_)
	: basis(b),	origin(c),	type(type_)
{

}
Transform::Transform(const Quaternion& q, const Vector3D& c , XFormTYPE type_)
	: basis(q),	origin(c),	type(type_)
{

}

Transform::~Transform()
{
	
}

Transform& Transform::operator =(Transform& t)
{
	setOrigin(t.getOrigin());
	setBasis(t.getBasis());

	return *this;
}

Transform Transform::operator *(Transform& t)
{
	Transform Xform;
	Xform.setOrigin( (*this) * t.getOrigin() );
	Xform.setBasis ( getBasis()  * t.getBasis() );
	Xform.type= (XFormTYPE)(type | t.type);
	return Xform;
}

Vector3D Transform::operator()(Vector3D& x)
{
	return Vector3D(basis[0].dot(x) + origin[0], 
					basis[1].dot(x) + origin[1], 
					basis[2].dot(x) + origin[2]);
}

Vector3D Transform::operator *(Vector3D& x)
{
	return Vector3D(basis[0].dot(x) + origin[0], 
					basis[1].dot(x) + origin[1], 
					basis[2].dot(x) + origin[2]);
}

Transform Transform::multInverseLeft(Transform& t) 
{
	Transform Xform;

	Vector3D v = t.origin - origin;
	if (type & SCALING) {
		Matrix3x3 inv = basis.inverse();
		Xform.basis = inv * t.basis;
		Xform.origin = inv * v;
	}
	else {
		Xform.basis = multTransposeLeft(basis, t.basis);
		Xform.origin = v * basis;
	}
	Xform.type = (XFormTYPE) (type | t.type);

	return Xform;
}

/** Set the origin of the current transform.
	\param o	The given origin	
*/
void Transform::setOrigin(Coord3D o)
{
	origin=o;
}

/** Set the origin of the current transform.
	\param x	The given origin x	
	\param y	The given origin y	
	\param z	The given origin z	
*/
void Transform::setOrigin(double x, double y, double z)
{
	origin.setCoords(x, y, z);
}

/** Return the origin of the current transform.
*/
Coord3D Transform::getOrigin(void)
{
	return origin;
}

/** Return the origin of the current transform.
*/
void Transform::getOrigin(double* o)
{
	o[0]=origin[0]; o[1]=origin[1]; o[2]=origin[2];
}

/** Set the rotation matrix of the current transform.
	\param m	The given rotation matrix	
*/
void Transform::setRotation(Quaternion q)
{
	basis.setRotation(q);
}

/** Return the rotation matrix of the current transform.
*/
Quaternion Transform::getRotation(void)
{
	return basis.getRotation();
}

/** Set the rotation matrix of the current transform.
	\param m	The given rotation matrix	
*/
void Transform::setBasis(Matrix3x3 m)
{
	basis=m;
}

/** Return the rotation matrix of the current transform.
*/
Matrix3x3 Transform::getBasis(void)
{
	return basis;
}

/** Return the rotation matrix of the current transform.
*/
void Transform::getBasis(double *b)
{
	b[0]=basis[0][0]; b[1]=basis[0][1]; b[2]=basis[0][2];
	b[3]=basis[1][0]; b[4]=basis[1][1]; b[5]=basis[1][2];
	b[6]=basis[2][0]; b[7]=basis[2][1]; b[8]=basis[2][2];
}

/** Set the current transform identity.
*/
void Transform::setIdentity(void)
{
	basis.setIdentity();
	origin.zero();
}

/** Set the current transform identity.
*/
bool Transform::isIdentity(void)
{
	return type==IDENTITY;
}

void Transform::setFromOpenGLMatrix(const double *m)
{
	basis.setFromOpenGLSubMatrix(m);
	origin[0] = m[12];
	origin[1] = m[13];
	origin[2] = m[14];
}

void Transform::getOpenGLMatrix(double *m) 
{
	basis.getOpenGLSubMatrix(m);
	m[12] = origin[0];
	m[13] = origin[1];
	m[14] = origin[2];
	m[15] = 1.0;
}

/** Scale the Xform with the given vector
	\param s The given vector	
*/
void Transform::scale(const Vector3D& s) const
{
	basis.scaled(s); 
}

Transform Transform::inverse()
{
	if (type)
	{
		Matrix3x3 inv = (type & SCALING) ? basis.inverse(): basis.transpose();

		return Transform(inv, inv * (-origin), type);
	}

	return *this;
}


Transform Transform::integrate(const Vector3D* v, const Vector3D* w, const double dt)
{
	Transform Xform;
	Xform.setOrigin(getOrigin() + (*v)* dt);

#ifdef QUATERNION_DERIVATIVE
	Quaternion orn = getRotation();
	orn += ((*v) * orn) * (dt * 0.5f);
	orn.normalize();
#else
	//exponential map
	Vector3D axis;
	double	angle = w->sqrabs();
	if( (angle*dt) > ANGULARMOTIONTRESHOLD)
	{
		angle = (ANGULARMOTIONTRESHOLD)/dt;
	}

	//limit the angular motion
	if ( angle < 0.001f )
	{
		// use Taylor's expansions of sync function
		axis   = (*w)*( 0.5f*dt-(dt*dt*dt)*(0.020833333333)*angle*angle );
	}
	else
	{
		// sync(fAngle) = sin(c*fAngle)/t
		axis   = (*w)*( sin(0.5f*angle*dt)/angle );
	}
	Quaternion dorn (axis[0],axis[1],axis[2],cos( angle*dt*0.5f ));
	Quaternion orn0 = getRotation();

	Quaternion predictedOrn = dorn * orn0;
#endif
	Xform.setRotation(predictedOrn);

	return Xform;
}

// get velocity by the now state and the final state of model
// Xform is the final state 
// dt time interval
// v the translational velocity
// w the rotation velocity

void Transform::velocity(Transform* xform, double dt, Vector3D& v, Vector3D& w)
{
	//linear velocity
	v = (xform->getOrigin() - getOrigin()) / dt;

#ifdef USE_QUATERNION_DIFF
	Quaternion orn0 = getRotation();
	Quaternion orn1a = xform->getRotation();
	Quaternion orn1 = orn0.farthest(orn1a);
	Quaternion dorn = orn1 * orn0.inverse();
#else
	Matrix3x3 dmat = xform->getBasis()* ( getBasis().inverse());;// 
	Quaternion dorn= dmat.getRotation();
#endif

	Vector3D axis;
	double  angle = dorn.angle();
	axis = Vector3D(dorn[0],dorn[1],dorn[2]);

	//check for axis length
	double len = axis.abs();
	if (len < 0.001f)
		axis = Vector3D(1.f,0.f,0.f);
	else
		axis /= sqrtf(len);

	//angular velocity	
	w = axis * angle / dt;

}
