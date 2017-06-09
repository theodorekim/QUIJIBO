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
#include "BOX.h"

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

bool spinning = true;

float fov = 45;
float aspect = 1.0;
float near = 0.001;
float far = 10.0;

GLVU glvu;
VEC3 center;
VEC3 lengths;
BOX box;

glvuVec3f glvuMin, glvuMax;

VEC3 objCenter;

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

  glvu.BeginFrame();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
    // draw coordinate axes
    glPushMatrix();
      glTranslatef(objCenter[0], objCenter[1], objCenter[2]);
      drawAxes();
    glPopMatrix();

    glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
    glColor4f(0.5f, 0.5f, 0.5f, 1.0f);

      glTranslatef(0.5, 0.5, 0.5);
      glScalef(0.25, 0.25, 0.25);
      glRotatef(-89.0, 1, 0, 0);
      glRotatef(-175.0, 0, 1, 0);
      glTranslatef(-0.1, -0.8, 0.3);
      glTranslatef(0,0,-1.41);
    objFile.draw();

    // draw a box outline
    glColor4f(1.0f, 1.0f, 1.0f, 10.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
      box.draw();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);

    // draw the center
    glColor4f(1.0f, 0.0f, 0.0f, 10.0f);
    glPointSize(10);
    VEC3F boxCenter = box.center();
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);
      glVertex3f(boxCenter[0], boxCenter[1], boxCenter[2]);
    glEnd();
    glEnable(GL_LIGHTING);

  glvu.EndFrame();
}

//////////////////////////////////////////////////////////////////////////////
// USER-PROVIDED KEYBOARD HANDLING ROUTINE
//////////////////////////////////////////////////////////////////////////////
void userKeyboardFunc(unsigned char Key, int x, int y)
{
  Camera* camera = glvu.GetCurrentCam();
  glvuVec3f eye;
  glvuVec3f ref;
  glvuVec3f up;
  camera->GetLookAtParams(&eye, &ref, &up);

  switch(Key)
  {
    case ' ':
      spinning = !spinning;
      break;
    case 'v':
      {
        cout << " Eye(" << eye[0] << ", " << eye[1] << ", " << eye[2]  << "),";
        cout << " LookAtCntr(" << ref[0] << ", " << ref[1] << ", " << ref[2]    << "),";
        cout << " Up(" << up[0] << ", " << up[1] << ", " << up[2]  << ");" << endl;
        cout << "<lookat target = \"" << ref[0] << "," << ref[1] << ","<< ref[2] << "\" " 
             << "origin = \"" << eye[0] << "," << eye[1] << "," << eye[2] << "\" "
             << "up = \"" << up[0] << ", " << up[1] << "," << up[2] << "\"/>" << endl;
        cout << "<float name=\"fov\" value=\"" << fov << "\"/>" << endl;
        break;
      }
    case 'f':
      fov -= 1.0;
      glvu.SetAllCams(glvuMin, glvuMax, eye, ref, up, fov, aspect, near, far);
      break;
    case 'F':
      fov += 1.0;
      glvu.SetAllCams(glvuMin, glvuMax, eye, ref, up, fov, aspect, near, far);
      break;
    case 'q':
    case 'Q':
      exit(0);
      break;
  }

  // if we're spinning the camera, all done
  if (spinning)
  {
    glvu.Keyboard(Key,x,y);
    return;
  }

  static int steps = 0;
  
  switch(Key)
  {
    case 's':
      {
        box.scale(0.95);
        break;
      }
    case 'S':
      {
        box.scale(1.05);
        break;
      }
    case 'x':
      {
        box.translateX(-0.025);
        break;
      }
    case 'X':
      {
        box.translateX(0.025);
        break;
      }

    case 'y':
      {
        box.translateY(-0.025);
        break;
      }
    case 'Y':
      {
        box.translateY(0.025);
        break;
      }

    case 'z':
      {
        box.translateZ(-0.025);
        break;
      }
    case 'Z':
      {
        box.translateZ(0.025);
        break;
      }
  };

  VEC3F center  = box.center();
  VEC3F lengths = box.lengths();
  cout << "center  = VEC3F(" << center[0] << "," << center[1] << "," << center[2] << ");" << endl;
  cout << "lengths = VEC3F(" << lengths[0] << "," << lengths[1] << "," << lengths[2] << ");" << endl;

  glutPostRedisplay();
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

  glutDisplayFunc(userDisplayFunc);
  glutKeyboardFunc(userKeyboardFunc);
  glutIdleFunc(userIdleFunc);

  // reset center
  //VEC3 objCenter = objFile.center();
  objCenter = objFile.center();

  objCenter += VEC3(0,0,-1.41);
  objCenter += VEC3(-0.1, -0.8, 0.3);
  objCenter = MATRIX3::rotation(VEC3(0,1,0), -175.0 / 360.0 * 2.0 * M_PI) * objCenter;
  objCenter = MATRIX3::rotation(VEC3(1,0,0), -89.0 / 360.0 * 2.0 * M_PI) * objCenter;
  objCenter *= 0.25;
  objCenter += VEC3(0.5, 0.5, 0.5);

  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        //Eye(0.5,0.5,3), LookAtCntr(0.5,0.5,0.5), Up(0,1,0);
        //Eye(-2.60453, -3.57762, 2.87471), LookAtCntr(-2.04369, -2.81653, 2.54882), Up(0.147047, 0.295801, 0.943862);
        //LookAtCntr(0.681017,0.344339, 0.841651), Eye(1.06275, -0.147389, 1.62427), Up(-0.180728, -0.870101, -0.458544);
        //Eye(1.02965, 0.247607, 0.31065), LookAtCntr(0.204498, -0.0502375, 0.79068), Up(-0.382829, -0.330058, -0.86285);
        //Eye(1.0611, 0.144932, 0.338933), LookAtCntr(0.181573, 0.136956, 0.814712), Up(-0.445088, -0.339855, -0.828491);
        //Eye(0.975366, 0.0865425, 0.430146), LookAtCntr(0.0329559, 0.23162, 0.731513), Up(-0.0605996, -0.960187, 0.272722);
        //Eye(0.848018, 0.105181, 0.412587), LookAtCntr(-0.0891528, -0.118941, 0.679945), Up(-0.0267212, 0.810215, 0.585525);
        //Eye(0.81415, 0.0598335, 0.437854), LookAtCntr(-0.183826, 0.112008, 0.401492), Up(0.0253124, 0.850417, 0.525507);
        //Eye(0.819703, 0.0661519, 0.39502), LookAtCntr(-0.102349, 0.00123858, 0.776606), Up(0.140805, 0.862045, 0.486882);
        //Eye(0.830368, 0.0599492, 0.400172), LookAtCntr(-0.120662, 0.0558909, 0.709244), Up(0.153079, 0.862497, 0.482356);
        Eye(0.760062, 0.0628541, 0.445488), LookAtCntr(-0.186215, 0.0501026, 0.768596), Up(0.153079, 0.862496, 0.482356);

  //glvuVec3f glvuCenter(objCenter[0], objCenter[1], objCenter[2]);
  //LookAtCntr = glvuCenter;
  //Eye = glvuCenter + glvuVec3f(0,0,2.5);

  //fov = 45;
  //fov = 37;
  fov = 45;
  //fov = 20;
  //fov = 45;
  aspect = 640.0 / 480.0;
  near = 0.001f;
  far = 10.0f;
  glvuMin = ModelMin;
  glvuMax = ModelMax;
  glvu.SetAllCams(ModelMin, ModelMax, Eye,LookAtCntr,Up, fov, aspect, near,far);

  //glvu.SetMoveSpeed(0.01);
  glvu.SetMoveSpeed(0.0001);
  
  glvuVec3f center(objCenter[0], objCenter[1], objCenter[2]);
  glvu.SetWorldCenter(center);
  glvu.SetInertiaEnabled(false);

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
  box = BOX(center, lengths);  

  glutInit(&argc, argv);
  glvuWindow();
  return 0;
}
