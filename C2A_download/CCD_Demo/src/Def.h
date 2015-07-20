/******************************************************************
** Copyright (c) 2005 Computer Graphics Laboratory
**					  Ewha Womans University, Seoul, Korea
** Filename:		CCDChecker.h
** Creator:			Zhang Xinyu
** Date:			2005-11-23
** Description:		Declare CCDChecker class
**
** Version:			1.0.1
**-----------------------------------------------------------------------------

**Revision history:
**
** Reviser:			date:				Description:
**-----------------------------------------------------------------------------
******************************************************************/
//Refer to DGPShop

#ifndef DEFINITION_H
#define DEFINITION_H

//#pragma warning(disable:4786)
//#pragma warning(disable:4251)
//#pragma warning(disable:4996)

#include <windows.h>
#include <vector>
#include <GL/glut.h>
#include <fstream>
#include <iostream>
//#include "./LinearMath/btTransform.h"
//#include "./LinearMath/btTransformUtil.h"


#define DEBUG_GRAPHICS
#define COMPUTE_EXE_TIME
#define SAMPLE_RATE 44100

typedef enum 
{ 
	OK,						//Succeed
	TOCFound,				//A time of collision is found
	CollisionFound,			//A collision is found
	CollisionFree,			//A collision is found
	CollisionNotFound,		//No collision is found
	IntersectionFound,		//An intersection is found at the current position
	PenetrationFound,		//A penetration is found at the current position
	MaxInterationExceed,	//A maximum interation has been exceeded
	LargeDifferenceFound,	//A large difference between in-between motion and dynamic motion
	AddingFileErrorFound	//A penetration is found at the current position
} RESULT;

typedef enum 
{
	NCAForNonConvex, 
	TRCAForNonConvex,
	NCAForConvex,
	BisectionForNonConvex
} CCDALGORITHM;
/*

class TimeOfCollision {
public:
	float t;
	float d;
	btVector3 n;
	btTransform TransformA;
	btTransform TransformB;
};


class NodeOfVelocity {
public:
	btVector3 vA;
	btVector3 vB;
	btVector3 wA;
	btVector3 wB;
	float radiusA;
	float radiusB;
	float vpmax;
	int	  level;
	NodeOfVelocity* leftleftChild;
	NodeOfVelocity* leftrightChild;
	NodeOfVelocity* rightleftChild;
	NodeOfVelocity* rightrightChild;
	NodeOfVelocity* parent;
};*/

//static float EPSILON_DISTANCE=0.001f;
#define EPSILON_TIME 0.001
#define EPSILON_DELTA_TIME 0.02
#define MAX_ITERATIONS 50
#define EPSILON_DISTANCE 0.001
#define EPSILON_DELTA_DISTANCE 0.00002 //=EPSILON_DISTANCE/MAX_ITERATIONS

static GLfloat mat_ambient_diffuse_0_i[]	={0.5f, 0.0f, 0.0f, 1.0f};
static GLfloat mat_ambient_diffuse_0_m[]	={1.0f, 0.5f, 0.1f, 0.5f};
static GLfloat mat_ambient_diffuse_0_toc[]	={1.0f, 0.0f, 0.0f, 1.0f};
static GLfloat mat_ambient_diffuse_0_front[]={0.2f, 0.2f, 0.0f, 1.0f};
static GLfloat mat_ambient_diffuse_1_i[]	={0.0f, 0.5f, 0.0f, 1.0f};
static GLfloat mat_ambient_diffuse_1_m[]	={0.0f, 1.0f, 0.1f, 0.5f};
static GLfloat mat_ambient_diffuse_1_toc[]	={0.0f, 1.0f, 0.0f, 1.0f};
static GLfloat mat_ambient_diffuse_1_front[]={0.0f, 0.2f, 0.2f, 1.0f};

#define ROTATE_SENSITIVITY		0.2
#define ZOOM_SENSITIVITY		0.005
#define MIN_ZOOM_FACTOR			0.0
#define MAX_ZOOM_FACTOR			100.0

#define CURSOR_SIZE_PIXELS		20
#define MAX_STRING_SIZE			255

typedef enum {SOLID=0, WIREFRAME} RENDERMODE;
typedef enum {NOTHING, TRANSLATE, SCALE, ROTATE, PAN, ZOOM, SPIN, PICKING, PAINT, NONE} TRANSFORM_MODE;

typedef enum {R=0, G, B} COLOR;

#ifndef PI
#define PI	3.1415926
#endif

#ifndef PI_M_2
#define PI_M_2	6.2831852
#endif

#ifndef PI_S_2
#define PI_S_2	1.5707963
#endif


#endif

