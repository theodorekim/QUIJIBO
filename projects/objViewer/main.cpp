/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/

//------------------------------------------------------------------------------
// GL interface elements are from:
//------------------------------------------------------------------------------
// GLVU : Copyright 1997 - 2002 
//        The University of North Carolina at Chapel Hill
//------------------------------------------------------------------------------
// Permission to use, copy, modify, distribute and sell this software and its 
// documentation for any purpose is hereby granted without fee, provided that 
// the above copyright notice appear in all copies and that both that copyright 
// notice and this permission notice appear in supporting documentation. 
// Binaries may be compiled with this software without any royalties or 
// restrictions. 
//
// The University of North Carolina at Chapel Hill makes no representations 
// about the suitability of this software for any purpose. It is provided 
// "as is" without express or implied warranty.

#include <iostream>
#include <cstring>
#include "OBJ.h"
//#include "BOX.h"
#include "FIELD_3D.h"

// In Windows, pop up a preview window
#include <glvu.h>
#include <snapshot.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// GLOBALS
//////////////////////////////////////////////////////////////////////////////
//int res = 64;
int res = 32;
//int res = 16;
OBJ objFile;
bool animate = false;
int DIV = 1048576;
char *divisor = "M";
int WIDTH = 4;
int startingFree;
float zSlice = 0.5f;

enum VIEWING_MODE { CAMERA_ADJUST = 0, TRANSLATING = 1, ROTATING = 2};
enum WHICH_AXIS { AXIS_X = 0, AXIS_Y = 1, AXIS_Z = 2};

//bool spinning = true;
int mode = CAMERA_ADJUST;
int rotationAxis = AXIS_X;
int translationAxis = AXIS_X;

GLVU glvu;
VEC3 center;
VEC3 lengths;
//BOX box;
FIELD_3D box;

GLuint objDisplayList = 0;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void buildObjDisplayList()
{
  objDisplayList = glGenLists(1);
  glNewList(objDisplayList, GL_COMPILE);
    objFile.draw();
  glEndList();
}

//////////////////////////////////////////////////////////////////////////////
// print which mode we're currently in
//////////////////////////////////////////////////////////////////////////////
void printMode()
{
  char buffer[256];
  if (mode == CAMERA_ADJUST)
    sprintf(buffer, "Camera viewing mode");
  if (mode == TRANSLATING)
  {
    char axis;
    if (translationAxis == AXIS_X)
      axis = 'X';
    if (translationAxis == AXIS_Y)
      axis = 'Y';
    if (translationAxis == AXIS_Z)
      axis = 'Z';
    sprintf(buffer, "Translating axis %c", axis);
  }
  if (mode == ROTATING)
  {
    char axis;
    if (rotationAxis == AXIS_X)
      axis = 'X';
    if (rotationAxis == AXIS_Y)
      axis = 'Y';
    if (rotationAxis == AXIS_Z)
      axis = 'Z';
    sprintf(buffer, "Rotating axis %c", axis);
  }
  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  //glvu.Print(5,700 - 15,buffer);
  glvu.Print(5,5,buffer);
}

//////////////////////////////////////////////////////////////////////////////
// draw the coordinate axes
//////////////////////////////////////////////////////////////////////////////
void drawAxes()
{
  glLineWidth(3.0f);
  glBegin(GL_LINES);
    // x axis is red
    glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    // y axis is green 
    glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 10.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    
    // z axis is blue
    glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
  glEnd();
  glLineWidth(1.0f);
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED DRAWING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userDisplayFunc()
{
  // get camera settings
  Camera* camera = glvu.GetCurrentCam();
  glvuVec3f Eye, Lookat, Up;
  camera->GetLookAtParams(&Eye, &Lookat, &Up);

  // repack camera settings into arrays
  float eye[] = {Eye.x, Eye.y, Eye.z};
  float look[3];
  look[0] = Lookat.x - eye[0];
  look[1] = Lookat.y - eye[1];
  look[2] = Lookat.z - eye[2];
  float magnitude = 1.0f / sqrt(look[0] * look[0] + look[1] * look[1] + look[2] * look[2]);
  look[0] *= magnitude;
  look[1] *= magnitude;
  look[2] *= magnitude;

  VEC3F boxCenter = box.center();
  glvu.BeginFrame();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
    // draw coordinate axes
    glPushMatrix();
      glTranslatef(boxCenter[0], boxCenter[1], boxCenter[2]);
      VEC3F axis;
      Real angle;
      box.rotation().axisAngle(axis, angle);
      glRotatef(angle, axis[0], axis[1], axis[2]);
      drawAxes();
    glPopMatrix();

    glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
    if (objDisplayList == 0)
      buildObjDisplayList();
    else
      glCallList(objDisplayList);
    //objFile.draw();

    // draw a box outline
    glColor4f(1.0f, 1.0f, 1.0f, 10.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
      box.drawBoundingBox();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);

    // draw the center
    glColor4f(1.0f, 0.0f, 0.0f, 10.0f);
    glPointSize(10);
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);
      glVertex3f(boxCenter[0], boxCenter[1], boxCenter[2]);
    glEnd();
    glEnable(GL_LIGHTING);

    // print which picking mode we're currently in
    printMode();

  glvu.EndFrame();
}

void printBox()
{
  VEC3F center  = box.center();
  VEC3F lengths = box.lengths();
  QUATERNION rotation = box.rotation();
  cout << "field center  = (" << center[0] << "," << center[1] << "," << center[2] << ")" << endl;
  cout << "field lengths = (" << lengths[0] << "," << lengths[1] << "," << lengths[2] << ")" << endl;
  cout << "rotation = (" << rotation[0] << "," << rotation[1] << "," << rotation[2] << "," << rotation[3] << ")" << endl;
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED KEYBOARD HANDLING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userKeyboardFunc(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 's':
      {
        box.scaleLengths(0.95);
        break;
      }
    case 'S':
      {
        box.scaleLengths(1.05);
        break;
      }
    case 'r':
      mode = ROTATING;
      break;
    case 't':
      mode = TRANSLATING;
      break;
    case ' ':
      mode = CAMERA_ADJUST;
      printBox();
      break;
    case 'v':
      {
        Camera* camera = glvu.GetCurrentCam();
        glvuVec3f eye;
        glvuVec3f ref;
        glvuVec3f up;
        camera->GetLookAtParams(&eye, &ref, &up);
        cout << " Eye(" << eye[0] << ", " << eye[1] << ", " << eye[2]  << "),";
        cout << " LookAtCntr(" << ref[0] << ", " << ref[1] << ", " << ref[2]    << "),";
        cout << " Up(" << up[0] << ", " << up[1] << ", " << up[2]  << ");" << endl;
        cout << "<lookat target = \"" << ref[0] << "," << ref[1] << ","<< ref[2] << "\" " 
             << "origin = \"" << eye[0] << "," << eye[1] << "," << eye[2] << "\" "
             << "up = \"" << up[0] << ", " << up[1] << "," << up[2] << "\"/>" << endl;
        break;
      }
    case 'q':
    case 'Q':
      exit(0);
      break;
  }

  // if we're spinning the camera, all done
  if (mode == CAMERA_ADJUST)
  {
    glvu.Keyboard(Key,x,y);
    return;
  }

  if (mode == TRANSLATING)
  {
    switch(Key)
    {
      case 'x':
        {
          translationAxis = AXIS_X;
          break;
        }
      case 'y':
        {
          translationAxis = AXIS_Y;
          break;
        }
      case 'z':
        {
          translationAxis = AXIS_Z;
          break;
        }
    };
  }

  if (mode == ROTATING)
  {
    Real delta = 0.01;
    switch(Key)
    {
      case 'x':
        {
          rotationAxis = AXIS_X;
          break;
        }
      case 'y':
        {
          rotationAxis = AXIS_Y;
          break;
        }
      case 'z':
        {
          rotationAxis = AXIS_Z;
          break;
        }
    };
  }


  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED MOUSE MOTION ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userMotionFunc(int x, int y)
{
  if (mode == CAMERA_ADJUST)
  {
    glvu.SetInertiaEnabled(true);
    glvu.Motion(x,y);
    return;
  }
  if (mode == TRANSLATING)
  {
    glvu.SetInertiaEnabled(false);

    int deltaX = glvu.GetMouseDeltaX(x);
    Real fraction = (Real)deltaX / 700;
    if (translationAxis == AXIS_X)
      box.translateX(fraction);
    if (translationAxis == AXIS_Y)
      box.translateY(fraction);
    if (translationAxis == AXIS_Z)
      box.translateZ(fraction);
    return;
  }

  if (mode == ROTATING)
  {
    glvu.SetInertiaEnabled(false);

    int deltaX = glvu.GetMouseDeltaX(x);
    Real fraction = (Real)deltaX / 700;
    if (rotationAxis == AXIS_X)
      box.rotateX(fraction);
    if (rotationAxis == AXIS_Y)
      box.rotateY(fraction);
    if (rotationAxis == AXIS_Z)
      box.rotateZ(fraction);
    return;
  }
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED IDLE ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userIdleFunc()
{
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  glvu.Init("GLVU Window",
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            100,100,640,480);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {0.7, 0.6, 0.6, 1.0};
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glClearColor(0,0,0,0);

  glutMotionFunc(userMotionFunc);
  glutDisplayFunc(userDisplayFunc);
  glutKeyboardFunc(userKeyboardFunc);
  glutIdleFunc(userIdleFunc);

  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        Eye(0.5,0.5,3), LookAtCntr(0.5,0.5,0.5), Up(0,1,0);
        //Eye(-2.60453, -3.57762, 2.87471), LookAtCntr(-2.04369, -2.81653, 2.54882), Up(0.147047, 0.295801, 0.943862);
        //LookAtCntr(0.681017,0.344339, 0.841651), Eye(1.06275, -0.147389, 1.62427), Up(-0.180728, -0.870101, -0.458544);
        //Eye(1.02965, 0.247607, 0.31065), LookAtCntr(0.204498, -0.0502375, 0.79068), Up(-0.382829, -0.330058, -0.86285);
        //Eye(1.0611, 0.144932, 0.338933), LookAtCntr(0.181573, 0.136956, 0.814712), Up(-0.445088, -0.339855, -0.828491);
        //Eye(0.975366, 0.0865425, 0.430146), LookAtCntr(0.0329559, 0.23162, 0.731513), Up(-0.0605996, -0.960187, 0.272722);

  // reset center
  VEC3 objCenter = objFile.center();
  glvuVec3f glvuCenter(objCenter[0], objCenter[1], objCenter[2]);

  LookAtCntr = glvuCenter;
  Eye = glvuCenter + glvuVec3f(0,0,2.5);

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye,LookAtCntr,Up, Yfov,Aspect, Near,Far);

  glvu.SetWorldCenter(glvuCenter);

  // recenter the box to the center of the mesh
  cout << " OBJ center: " << objCenter << endl;
  box.setCenter(objCenter);

  glvu.SetMoveSpeed(0.1);

  /*
//box.center() = VEC3F(-1.00856,1.02805,2.63454);
box.setCenter(VEC3F(-1.00856,1.02805,2.63454));
box.setLengths(VEC3F(0.302194,0.302194,0.302194));
box.rotation() = QUATERNION(0.893501,0.449061,0,0);
*/

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
// generate test from implicit functions
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    cout << " Usage: " << argv[0] << " <OBJ filename> " << endl;
    return 0;
  }

  //objFile.Load("./meshes/bunny_small_zoom_0.res.50.iterations.3_translationz.obj");
  //objFile.Load("./meshes/bunny_zoom_0.res.1100.iterations.3.obj");
  objFile.Load(argv[1]);
  objFile.ComputeVertexNormals();

  center = VEC3(0,0,0);
  //center = VEC3(0,0,1.41);
  lengths = VEC3(3.06186, 3.06186, 3.06186);
  //box = BOX(center, lengths);  
  box = FIELD_3D(10,10,10,center, lengths);  

  glutInit(&argc, argv);
  glvuWindow();
  return 0;
}
