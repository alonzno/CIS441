/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#if USING_ROOM_100_MAC
#include "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk/System/Library/Frameworks/OpenGL.framework/Headers/gl.h"
#endif

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

//I add
#include <tgmath.h>
#include <set>


class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double         fieldValue[3]; // always between 0 and 1
      double         normals[3][3];
};

//
// Function: GetTriangles
//
// Purpose: reads triangles from the geometry file: proj1e_geometry.vtk.
//
std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        tris[idx].fieldValue[0] = (color_ptr[ptIds[0]]-1)/5.;
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        tris[idx].fieldValue[1] = (color_ptr[ptIds[1]]-1)/5.;
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
        tris[idx].fieldValue[2] = (color_ptr[ptIds[2]]-1)/5.;
    }

    return tris;
}

//
// Function: GetColorMap
//
// Purpose: returns a 256x3 array of colors
//
unsigned char *
GetColorMap(void)
{
    unsigned char controlPts[8][3] =
    {
        {  71,  71, 219 },
        {   0,   0,  91 },
        {   0, 255, 255 },
        {   0, 127,   0 },
        { 255, 255,   0 },
        { 255,  96,   0 },
        { 107,   0,   0 },
        { 224,  76,  76 },
    };
    int textureSize = 256;
    unsigned char *ptr = new unsigned char[textureSize*3];
    int nControlPts = 8;
    double amountPerPair = ((double)textureSize-1.0)/(nControlPts-1.0);
    for (int i = 0 ; i < textureSize ; i++)
    {
        int lowerControlPt = (int)(i/amountPerPair);
        int upperControlPt = lowerControlPt+1;
        if (upperControlPt >= nControlPts)
            upperControlPt = lowerControlPt; // happens for i == textureSize-1

        double proportion = (i/amountPerPair)-lowerControlPt;
        for (int j = 0 ; j < 3 ; j++)
            ptr[3*i+j] = controlPts[lowerControlPt][j]
                       + proportion*(controlPts[upperControlPt][j]-
                                     controlPts[lowerControlPt][j]);
    }

    return ptr;
}


class vtk441Mapper : public vtkOpenGLPolyDataMapper
{
  protected:
   GLuint displayList;
   bool   initialized;
   std::vector<Triangle> tris;
   unsigned char * c_map;


  public:
   vtk441Mapper()
   {
     initialized = false;
     tris = GetTriangles();
     c_map = GetColorMap();
   }
    
   void
   RemoveVTKOpenGLStateSideEffects()
   {
     float Info[4] = { 0, 0, 0, 1 };
     glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Info);
     float ambient[4] = { 1,1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
     float diffuse[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
     float specular[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
   }


   void SetupLight(void)
   {
       glEnable(GL_LIGHTING);
       glEnable(GL_LIGHT0);
       GLfloat diffuse0[4] = { 0.8, 0.8, 0.8, 1 };
       GLfloat ambient0[4] = { 0.2, 0.2, 0.2, 1 };
       GLfloat specular0[4] = { 0.0, 0.0, 0.0, 1 };
       GLfloat pos0[4] = { 1, 2, 3, 0 };
       glLightfv(GL_LIGHT0, GL_POSITION, pos0);
       glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
       glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
       glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
       glDisable(GL_LIGHT1);
       glDisable(GL_LIGHT2);
       glDisable(GL_LIGHT3);
       glDisable(GL_LIGHT5);
       glDisable(GL_LIGHT6);
       glDisable(GL_LIGHT7);
   }
};

/* DrawBox
* 
* Draws a 20x20x20 white box centered on (0,0,0)
*
*/
void DrawBox(void) {
    unsigned char paths[16] = {0, 1, 3, 2,
                               6, 2, 0, 4,
                               5, 4, 6, 7,
                               3, 7, 5, 1};

    glDisable(GL_TEXTURE_1D);

    for (int i = 0; i < 4; i++) {
        glBegin(GL_LINE_STRIP);
        glColor3ub(255, 255, 255);
        for (int j = 0; j < 4; j++) {
            float a[3] = {10.0f, 10.0f, 10.0f};
            if (paths[i*4 + j] & 4) a[0] = -10.0f;
            if (paths[i*4 + j] & 2) a[1] = -10.0f;
            if (paths[i*4 + j] & 1) a[2] = -10.0f;
            glVertex3f(a[0], a[1], a[2]);
        }
        glEnd();
    }
}

class vtk441MapperPart1 : public vtk441Mapper
{
 public:
   static vtk441MapperPart1 *New();
   
   
   /*
    *
    *  REVISE THIS
    *
    */

   virtual void RenderPiece(vtkRenderer *ren, vtkActor *act)
   {
      
      RemoveVTKOpenGLStateSideEffects();
      
      glEnable(GL_COLOR_MATERIAL);
      glDisable(GL_TEXTURE_1D);
      SetupLight();

      // Write all the triangles
      glBegin(GL_TRIANGLES);
      {
        unsigned char *buff;

        for (Triangle t: tris) {
            for (int i = 0; i < 3; i++) {
                buff = &c_map[3 * (int) roundf(255 * t.fieldValue[i])];
                glColor3ub(buff[0], buff[1], buff[2]);
                glNormal3f(t.normals[i][0], t.normals[i][1], t.normals[i][2]);
                glVertex3f(t.X[i], t.Y[i], t.Z[i]);
            }
        }
      } 
      glEnd();

      // Drawing the white box
      DrawBox();
   }
};

vtkStandardNewMacro(vtk441MapperPart1);

class vtk441MapperPart2 : public vtk441Mapper
{
 public:
   static vtk441MapperPart2 *New();
   
   GLuint displayList;
   bool   initialized;

   vtk441MapperPart2()
   {
     initialized = false;
   }
   
   /*
    *
    *  REVISE THIS
    *
    */
   
   virtual void RenderPiece(vtkRenderer *ren, vtkActor *act)
   {
       RemoveVTKOpenGLStateSideEffects();
       SetupLight();
      
       glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, c_map);
       glEnable(GL_COLOR_MATERIAL);
       //glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
       glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
       glEnable(GL_TEXTURE_1D);

      // Write all the triangles
      glBegin(GL_TRIANGLES);
      {
        unsigned char *buff;

        for (Triangle t: tris) {
            for (int i = 0; i < 3; i++) {
                glTexCoord1f(t.fieldValue[i]);
                glNormal3f(t.normals[i][0], t.normals[i][1], t.normals[i][2]);
                glVertex3f(t.X[i], t.Y[i], t.Z[i]);
            }
        }
      } 
      glEnd();
      
      // Drawing the white box
      DrawBox();
   }
};

vtkStandardNewMacro(vtk441MapperPart2);


int main()
{
  // Dummy input so VTK pipeline mojo is happy.
  //
  vtkSmartPointer<vtkSphereSource> sphere =
    vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetThetaResolution(100);
  sphere->SetPhiResolution(50);

  // The mapper is responsible for pushing the geometry into the graphics
  // library. It may also do color mapping, if scalars or other attributes
  // are defined. 
  //
  vtkSmartPointer<vtk441MapperPart1> win1Mapper =
    vtkSmartPointer<vtk441MapperPart1>::New();
  win1Mapper->SetInputConnection(sphere->GetOutputPort());

  vtkSmartPointer<vtkActor> win1Actor =
    vtkSmartPointer<vtkActor>::New();
  win1Actor->SetMapper(win1Mapper);

  vtkSmartPointer<vtkRenderer> ren1 =
    vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtk441MapperPart2> win2Mapper =
    vtkSmartPointer<vtk441MapperPart2>::New();
  win2Mapper->SetInputConnection(sphere->GetOutputPort());

  vtkSmartPointer<vtkActor> win2Actor =
    vtkSmartPointer<vtkActor>::New();
  win2Actor->SetMapper(win2Mapper);

  vtkSmartPointer<vtkRenderer> ren2 =
    vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);
  ren1->SetViewport(0, 0, 0.5, 1);
  renWin->AddRenderer(ren2);
  ren2->SetViewport(0.5, 0, 1.0, 1);

  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);

  // Add the actors to the renderer, set the background and size.
  //
  bool doWindow1 = true;
  if (doWindow1)
     ren1->AddActor(win1Actor);
  ren1->SetBackground(0.0, 0.0, 0.0);
  bool doWindow2 = true;
  if (doWindow2)
      ren2->AddActor(win2Actor);
  ren2->SetBackground(0.0, 0.0, 0.0);
  renWin->SetSize(1200, 600);

  // Set up the lighting.
  //
  vtkRenderer *rens[2] = { ren1, ren2 };
  for (int i = 0 ; i < 2 ; i++)
  {
     rens[i]->GetActiveCamera()->SetFocalPoint(0,0,0);
     rens[i]->GetActiveCamera()->SetPosition(0,0,70);
     rens[i]->GetActiveCamera()->SetViewUp(0,1,0);
     rens[i]->GetActiveCamera()->SetClippingRange(20, 120);
     rens[i]->GetActiveCamera()->SetDistance(70);
  }
  
  // This starts the event loop and invokes an initial render.
  //
  ((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
  iren->Initialize();
  iren->Start();

  return EXIT_SUCCESS;
}




