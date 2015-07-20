#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include "model.h"
#include "PQP.h"
#include "MatVec.h"
#include "Stopwatch.h"
#include <iostream>
#include <fstream>
#include "Transform.h"



#pragma warning(disable: 4996)

int nframes;
int iframe=0;
int nTests=20;

#include "mainccd.h"


PQP_Model  *object1_tested,*object2_tested;
Model      *object1_drawn, *object2_drawn;



StopwatchWin32 timer;
float timing;

int animating = 1;
int mode;
double beginx, beginy;
double dis = 10.0, azim = 0.0, elev = 0.0;
double ddis = 0.0, dazim = 0.0, delev = 0.0;

int step = 0;
int number_of_steps;
int query_type = 1;
double tolerance = .05;



static const char* modelFile1 = "../models/torusknot_3K.tri";
static const char* modelFile2 = "../models/torusknot_3K.tri";


static const char* aniPath1 = "../models/torusknot1.ani";
static const char* aniPath2 = "../models/torusknot2.ani";


void
init_viewer_window()
{
  GLfloat Ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };  
  GLfloat Diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };  
  GLfloat Specular[] = { 0.1f, 0.1f, 0.1f, 1.0f };   
  GLfloat SpecularExp[] = { 50 };              
  GLfloat Emission[] = { 0.1f, 0.1f, 0.1f, 1.0f };   

  glMaterialfv(GL_FRONT, GL_AMBIENT, Ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, Diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, Specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, SpecularExp);
  glMaterialfv(GL_FRONT, GL_EMISSION, Emission);

  glMaterialfv(GL_BACK, GL_AMBIENT, Ambient);
  glMaterialfv(GL_BACK, GL_DIFFUSE, Diffuse);
  glMaterialfv(GL_BACK, GL_SPECULAR, Specular);
  glMaterialfv(GL_BACK, GL_SHININESS, SpecularExp);
  glMaterialfv(GL_BACK, GL_EMISSION, Emission);

  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

  glEnable(GL_COLOR_MATERIAL);
  glutInitWindowSize(700, 700); 
  GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
   
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glShadeModel(GL_FLAT);
  glClearColor(1, 1, 1, 0);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-350,350,-350,350,-350,350);   
 //   glOrtho(-35,35,-35,35,-35,35); 
  glMatrixMode(GL_MODELVIEW);

}

void
cb_mouse(int _b, int _s, int _x, int _y)
{
  if (_s == GLUT_UP)
  {
    dis += ddis;
    if (dis < .1) dis = .1;
    azim += dazim;
    elev += delev;
    ddis = 0.0;
    dazim = 0.0;
    delev = 0.0;
    return;
  }

  if (_b == GLUT_RIGHT_BUTTON)
  {
    mode = 0;
    beginy = _y;
    return;
  }
  else
  {
    mode = 1;
    beginx = _x;
    beginy = _y;
  }

}

void
cb_motion(int _x, int _y)
{
  if (mode == 0)
  {
    ddis = dis * (double)(_y - beginy)/0.5;
  }
  else
  {
    dazim = (_x - beginx)/5;
    delev = (_y - beginy)/5;      
  }  
  glutPostRedisplay();
}

void cb_keyboard(unsigned char key, int x, int y) 
{
  switch(key) 
  {
  case 'q': 
    delete object1_drawn; 
    delete object2_drawn; 
    delete object1_tested;
    delete object2_tested;
	
    exit(0);
  case '0': query_type = 0; break;
  case '1': query_type = 1; break;
  case 'p': outputPerformance(timingFileName); break;
  case 't': 
  default: {animating = 1 - animating;
	  printf("NowFrame: %d\n",iframe);}
  }

  glutPostRedisplay();
}




void
BeginDraw()
{
	
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glLoadIdentity(); 
  glTranslatef(0.0, 0.0, -(dis+ddis));
  glRotated(elev+delev, 1.0, 0.0, 0.0);
  glRotated(azim+dazim, 0.0, 1.0, 0.0);
  glRotated(90.0,-1.0,0.0,0.0);
}

void
EndDraw()
{
  glFlush();
  glutSwapBuffers();
}

inline void glVertex3v(float V[3]) { glVertex3fv(V); }
inline void glVertex3v(double V[3]) { glVertex3dv(V); }

PQP_REAL toc;
int nItr;
int NTr;

void
cb_display()//Display
{  

  BeginDraw();
  double oglm[16];
  if(iframe>=nframes)	{
	  iframe = 0;
  }
  int step1, step2,step3,step4;
  if(iframe<101){
	  step1=0;		step2=1;
  }
  else if (iframe<202)
  {
	  step1=101;		step2=102;
  }else{
	  step1=202;		step2=203;
  }


  step3=iframe;		step4=iframe+1;
  if (iframe==101) {
	  step3=iframe;		step4=iframe;
  }
  if (iframe==202) {
	  step3=iframe;		step4=iframe;
  }
  if (iframe==302)
  {
	  step3=iframe;		step4=iframe;
  }

  TransRT tr0;

  Transform trans0;
  Transform trans1;
  Transform trans00;
  Transform trans01;
  Transform trans10;
  Transform trans11;
  for(int i=0; i<3; i++)
  {
	  for(int j=0; j<3; j++)
	  {
		  trans00.basis[i][j]= xformframes[0][step1].R[i][j];
		  trans01.basis[i][j]= xformframes[0][step2].R[i][j];
		  trans10.basis[i][j]= xformframes[1][step3].R[i][j];
		  trans11.basis[i][j]= xformframes[1][step4].R[i][j];
	  }
	  trans00.origin[i]= xformframes[0][step1].T[i];
	  trans01.origin[i]= xformframes[0][step2].T[i];
	  trans10.origin[i]= xformframes[1][step3].T[i];
	  trans11.origin[i]= xformframes[1][step4].T[i];
  }

  switch (query_type)
  {
  case 0:

    // draw model 1

	//red-initial A
    glColor3f(1.0,0.0,0);                 // setup color and transform
    MVtoOGL(oglm,xformframes[0][step1].R,xformframes[0][step1].T);
    glPushMatrix();
    glMultMatrixd(oglm);
    object1_drawn->Draw();             // do gl rendering
    glPopMatrix();                    // restore transform

	//blue final A
    glColor3f(0,0,1.0);                 // setup color and transform
    MVtoOGL(oglm,xformframes[0][step2].R,xformframes[0][step2].T);
    glPushMatrix();
    glMultMatrixd(oglm);
    object1_drawn->Draw();             // do gl rendering
    glPopMatrix();                    // restore transform

    //draw model 2
	//yellow rotating B 
    glColor3f(1,1,0.0);                 // setup color and transform
    MVtoOGL(oglm,xformframes[1][iframe].R,xformframes[1][iframe].T);
    glPushMatrix();
    glMultMatrixd(oglm);
    object2_drawn->Draw();
    glPopMatrix();

    break;

  case 1:
    // perform continuous collision detection
  
	PQP_CollideResult res_collide;
	
	int flag=2;

	
	timer.Reset(); 
	timer.Start();
    for(int iTest=0;iTest<nTests;iTest++) 
	{		  
		
		CA( &trans00, &trans01, object1_tested, 
			&trans10, &trans11, object2_tested, 
			trans0, trans1, toc, nItr ,NTr);  

	}


	timer.Stop();
	tframes[iframe]=timer.GetTime()/nTests;


	iterframes[iframe] = nItr;
	tocframes[iframe]  = toc;
	distframes[iframe] = dres.Distance();
	contframes[iframe] = NTr;
	FrontNodes[iframe] = dres.num_frontnodes;

    glColor3f(1, 0,0);                 // setup color and transform
    MVtoOGL(oglm,xformframes[0][step1].R,xformframes[0][step1].T);
    glPushMatrix();
    glMultMatrixd(oglm);
    object1_drawn->Draw();             // do gl rendering
    glPopMatrix();                    // restore transform

	//blue final A	
 
 

	glColor3f(0, 0,1);                 // setup color and transform
	MVtoOGL(oglm,xformframes[0][step2].R,xformframes[0][step2].T);
	glPushMatrix();
	glMultMatrixd(oglm);
	object1_drawn->Draw();             // do gl rendering
	glPopMatrix();                    // restore transform

  	//green- toc A
	glColor3f(0, 1,0);	
	if (dres.collisionfree)// 
	{
 
     MVtoOGL(oglm,xformframes[0][step2].R,xformframes[0][step2].T);
     glPushMatrix();
     glMultMatrixd(oglm);
     object1_drawn->Draw();            
     glPopMatrix();     

	 Transform2PQP(&trans11, tr0.R, tr0.T);
	 glColor3f(1, 1, 0);                 // setup color and transform
	 MVtoOGL(oglm,tr0.R,tr0.T);
	 glPushMatrix();
	 glMultMatrixd(oglm);
	 object2_drawn->Draw();
	 glPopMatrix();
		
	}
	else
	{
		// Fetch R and T of the contact configuration
		Transform2PQP(&trans0, tr0.R, tr0.T);
	
		MVtoOGL(oglm,tr0.R,tr0.T);
		glPushMatrix();
		glMultMatrixd(oglm);
		object1_drawn->Draw();
		glPopMatrix();

		Transform2PQP(&trans1, tr0.R, tr0.T);
		glColor3f(1, 1, 0);                 // setup color and transform
		MVtoOGL(oglm,tr0.R,tr0.T);
		glPushMatrix();
		glMultMatrixd(oglm);
		object2_drawn->Draw();
		glPopMatrix();
    
	}

	//yellow rotating B 
    // draw model 2



    break;
  }

  EndDraw();

  iframe++;
  if(animating)
	glutPostRedisplay();
}

void initOpengl()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);

	// initialize OpenGL lighting
	GLfloat matAmb[4] =    {1.0, 1.0, 1.0, 1.0};
	GLfloat matDiff[4] =   {1.0, 0.1, 0.2, 1.0};
	GLfloat matSpec[4] =   {1.0, 1.0, 1.0, 1.0};

	GLfloat lightpos1[] = {1.0, -1.0, -1.0, 0};
	GLfloat lightpos2[] = {-1.0, 1.0, 1.0, 0};

	GLfloat lightAmb[4] =  {0.0, 0.0, 0.0, 1.0};
	GLfloat lightDiff[4] = {1.0, 1.0, 1.0, 1.0};
	GLfloat lightSpec[4] = {1.0, 1.0, 1.0, 1.0};

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmb);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiff);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matSpec);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 60.0);

	glLightfv(GL_LIGHT0, GL_POSITION, lightpos1);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiff);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);

	glLightfv(GL_LIGHT1, GL_POSITION, lightpos2);
	glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmb);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiff);
	glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpec);

//	glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL_EXT, GL_SEPARATE_SPECULAR_COLOR_EXT);
	GLfloat black[] =  {0.0, 0.0, 0.0, 1.0};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, black);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
}

int main(int argc, char **argv)
{
  // init glut
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA | GLUT_STENCIL);

  glutInitWindowSize(720, 480);
  glutInitWindowPosition(550,400);
  glutCreateWindow("CCD using PQP");

  // set OpenGL graphics state -- material props, perspective, etc.

  init_viewer_window();

  // set the callbacks

  glutDisplayFunc(cb_display);
  glutMouseFunc(cb_mouse);
  glutMotionFunc(cb_motion);  
  glutKeyboardFunc(cb_keyboard);

  // create models

  FILE *fp;
  int ntris, nvertics,i;
  double a,b,c;
  char str[10];
  PQP_REAL p1[3],p2[3],p3[3];
  PQP_REAL *p;
  int i1,i2,i3;
  
  //reading tri files

  // model 1 
  printf("reading the file: %s\n",modelFile1);

  object1_drawn = new Model(modelFile1);
  object1_tested = new PQP_Model();
  fp = fopen(modelFile1,"r");
  fscanf(fp,"%s",str);
  fscanf(fp,"%d",&nvertics);
  fscanf(fp,"%d",&ntris);

  printf("building the BVH for the model: %s\n",modelFile1);
  object1_tested->BeginModel();
  p=new PQP_REAL[3*nvertics];
  for (i = 0; i < nvertics; i++)
  {
	  fscanf(fp,"%lf %lf %lf",&a,&b,&c);
	  p[3*i+0] = (PQP_REAL)a;
	  p[3*i+1] = (PQP_REAL)b;
	  p[3*i+2] = (PQP_REAL)c;
	
  }
  
  for (i = 0; i < ntris; i++)
  {
    fscanf(fp,"%d %d %d",&i1,&i2,&i3);
    p1[0] = p[i1*3];
    p1[1] = p[i1*3+1];
    p1[2] = p[i1*3+2];

	p2[0] = p[i2*3];
    p2[1] = p[i2*3+1];
    p2[2] = p[i2*3+2];
	
	p3[0] = p[i3*3];
    p3[1] = p[i3*3+1];
    p3[2] = p[i3*3+2];
  
    object1_tested->AddTri(p1,p2,p3,i,i1,i2,i3);
  }
  object1_tested->EndModel(); 


  object1_lut.Create(object1_tested);

  fclose(fp);
  delete []p;

  // model 2

  printf("reading the file: %s\n",modelFile2);
  object2_drawn = new Model(modelFile2);
  object2_tested = new PQP_Model();

  
  fp = fopen(modelFile2,"r");
  fscanf(fp,"%s",str);
  fscanf(fp,"%d",&nvertics);
  fscanf(fp,"%d",&ntris);
  
  printf("building the BVH for the model: %s\n",modelFile1);
  object2_tested->BeginModel();
  p=new PQP_REAL[3*nvertics];
  for (i = 0; i < nvertics; i++)
  {
	  fscanf(fp,"%lf %lf %lf",&a,&b,&c);
	  p[3*i+0] = (PQP_REAL)a;
	  p[3*i+1] = (PQP_REAL)b;
	  p[3*i+2] = (PQP_REAL)c;	  
	  
  }  

  for (i = 0; i < ntris; i++)
  {
	  fscanf(fp,"%d %d %d",&i1,&i2,&i3);
	  p1[0] = p[i1*3];
	  p1[1] = p[i1*3+1];
	  p1[2] = p[i1*3+2];
	  
	  p2[0] = p[i2*3];
	  p2[1] = p[i2*3+1];
	  p2[2] = p[i2*3+2];
	  
	  p3[0] = p[i3*3];
	  p3[1] = p[i3*3+1];
	  p3[2] = p[i3*3+2];
	  
	  object2_tested->AddTri(p1,p2,p3,i,i1,i2,i3);
  }
  
 object2_tested->EndModel();  

 object2_lut.Create(object2_tested);

 fclose(fp);
 delete []p;
  
 readXFormFrames(aniPath1, aniPath2);


  // print instructions
  printf("PQP Demo - Falling:\n"
         "Press:\n"
         "0 - initial configurations, just animation\n"
         "    the initial A is shown in red\n"
         "    the initial B is shown in blue\n"
         "    the initial/final A is shown in yellow.\n"
         "1 - continus collision query\n"
         "    the toc A is shown in green.\n"
         "p - profiling into %s\n"
         "any other key to toggle animation on/off\n", timingFileName);

  // Enter the main loop.

  glutMainLoop();
}



void readXFormFrames( const char* filename1, const char* filename2 )
{
	ifstream fin1(filename1);
	ifstream fin2(filename2);
	float nf;
	float i;
	TransRT *xform;
	char c[5];

	fin1>>nf>>c;
	nframes=nf;
	while(!fin1.eof())
	{
		fin1>>i>>c;
		xform=&xformframes[0][(int)i];
		fin1>>xform->R[0][0]>>xform->R[1][0]>>xform->R[2][0]
			>>xform->R[0][1]>>xform->R[1][1]>>xform->R[2][1]
			>>xform->R[0][2]>>xform->R[1][2]>>xform->R[2][2];
		fin1>>xform->T[0]>>xform->T[1]>>xform->T[2];
	}
	cout<<"Read animation file: <"<<filename1<<">"<<endl;

	fin2>>nf>>c;
	while(!fin2.eof())
	{
		fin2>>i>>c;
		xform=&xformframes[1][(int)i];
		fin2>>xform->R[0][0]>>xform->R[1][0]>>xform->R[2][0]
			>>xform->R[0][1]>>xform->R[1][1]>>xform->R[2][1]
			>>xform->R[0][2]>>xform->R[1][2]>>xform->R[2][2];
		fin2>>xform->T[0]>>xform->T[1]>>xform->T[2];
	}
	cout<<"Read animation file: <"<<filename2<<">"<<endl;
}




















































































