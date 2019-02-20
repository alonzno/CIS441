#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

#include <math.h>

#define X_SCALE 1786
#define Y_SCALE 1344

//#define Y_SCALE 1786
//#define X_SCALE 1344

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

using std::cerr;
using std::endl;

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

inline bool eq(double a, double b) {
    //Change this
    return a == b;
    //return fabs(a - b) < 0.1;

}


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      unsigned char color[3];

    void copy(const Triangle &other) {
        memcpy(X, other.X, 3 * sizeof(double));
        memcpy(Y, other.Y, 3 * sizeof(double));
        memcpy(color, other.color, 3 * sizeof(unsigned char));
    }

  // would some methods for transforming the triangle in place be helpful?
    const void print() {
        cerr << "Triangle " << this << ":" << endl;
        for (int i = 0; i < 3; i++) {
            cerr << X[i] << "\t" << Y[i] << endl;
        }
    }
    
    //Place the vertices in CW order, v1 is top left, v2 is bottom/top right
    //v3 is top(left),
    void sortVertices() {
        double temp;
        if (Y[0] < Y[1]) {
            temp = X[0];
            X[0] = X[1];
            X[1] = temp;
            
            temp = Y[0];
            Y[0] = Y[1];
            Y[1] = temp;
        }
        
        if (Y[1] < Y[2]) {
            temp = X[1];
            X[1] = X[2];
            X[2] = temp;
            
            temp = Y[1];
            Y[1] = Y[2];
            Y[2] = temp;
        }
        
        if (Y[0] < Y[1]) {
            temp = X[0];
            X[0] = X[1];
            X[1] = temp;
            
            temp = Y[0];
            Y[0] = Y[1];
            Y[1] = temp;
        }
        
        //Case of downward tri
        if (eq(Y[0], Y[1])) {
            if (X[0] > X[1]) {
               temp = X[0];
               X[0] = X[1];
               X[1] = temp;
               
               temp = Y[0];
               Y[0] = Y[1];
               Y[1] = temp;
            }
        }
       
        //Case of upward tri
        if (eq(Y[1], Y[2])) {
            if (X[1] < X[2]) {
               temp = X[1];
               X[1] = X[2];
               X[2] = temp;
               
               temp = Y[1];
               Y[1] = Y[2];
               Y[2] = temp;
            }
        }
    }
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
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
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }
    cerr << "Done reading" << endl;

    return tris;
}

// Code that I've added

void writePoint(int x, int y, const unsigned char *color, unsigned char *buffer) {
    // Perform bounds checking
    if (x < 0 || x >= X_SCALE || y < 0 || y >= Y_SCALE) return;
     
    int loc = 3 * (y*X_SCALE + x);
    buffer[loc + 0] = color[0];
    buffer[loc + 1] = color[1];
    buffer[loc + 2] = color[2];

}

void rasterUpTriangle(const Triangle &t, unsigned char *buffer){
    if (t.Y[0] < 0.0) {
        return;
    }
    double row_min = ceil_441(t.Y[2]);
    double row_max = floor_441(t.Y[0]);
    
    if (eq(t.X[0], t.X[1])) {
        double l_m = (t.Y[2] - t.Y[0]) / (t.X[2] - t.X[0]);
        double l_b = t.Y[0] - l_m * t.X[0];
        //only move left
        double scan_XL;
        int left_end;
        int right_end = floor_441(t.X[1]);

        for (int i = row_min; i <= row_max; i++) {
            scan_XL = (i - l_b) / l_m;
            left_end  = (int) ceil_441(scan_XL); 

            //Fill Buffer for a row
            for (int j = left_end; j <= right_end; j++) {
                writePoint(j, i, t.color, buffer);
            }
        }
        return;
    }
    if (eq(t.X[0], t.X[2])) {
        double r_m = (t.Y[1] - t.Y[0]) / (t.X[1] - t.X[0]);
        double r_b = t.Y[1] - r_m * t.X[1];
        //only move right
        double scan_XR;
        int right_end;
        int left_end = ceil_441(t.X[0]);
        //cerr << row_min << "\t" << row_max << endl;

        for (int i = row_min; i <= row_max; i++) {
            scan_XR = (i - r_b) / r_m;
            right_end = (int) floor_441(scan_XR);

            for (int j = left_end; j <= right_end; j++) {
                if (i < 0) break; // Don't draw below screen
                writePoint(j, i, t.color, buffer);
            }
        }
        return;
    }


    double l_m = (t.Y[2] - t.Y[0]) / (t.X[2] - t.X[0]);
    double r_m = (t.Y[1] - t.Y[0]) / (t.X[1] - t.X[0]);
    
    double l_b = t.Y[0] - l_m * t.X[0];
    double r_b = t.Y[1] - r_m * t.X[1];

    double scan_XL, scan_XR;
    int right_end, left_end;

    for (int i = row_min; i <= row_max; i++) {
        
        // Bound Checking left right
        scan_XR = (i - r_b) / r_m;
        scan_XL = (i - l_b) / l_m;

        right_end = (int) floor_441(scan_XR);
        left_end  = (int) ceil_441(scan_XL); 

        //Fill Buffer for a row
        for (int j = left_end; j <= right_end; j++) {
            if (i < 0) break; // Don't draw below screen
            writePoint(j, i, t.color, buffer);
        }
    }
}

void rasterDownTriangle(const Triangle &t, unsigned char *buffer) {
    double row_min = ceil_441(t.Y[2]);
    double row_max = floor_441(t.Y[0]);
    
    if (eq(t.X[1], t.X[2])) {
        double l_m = (t.Y[2] - t.Y[0]) / (t.X[2] - t.X[0]);
        double l_b = t.Y[0] - l_m * t.X[0];
        //only move left
        double scan_XL;
        int left_end;
        int right_end = floor_441(t.X[1]);

        for (int i = row_min; i <= row_max; i++) {
            scan_XL = (i - l_b) / l_m;
            left_end  = (int) ceil_441(scan_XL); 

            //Fill Buffer for a row
            for (int j = left_end; j <= right_end; j++) {
                writePoint(j, i, t.color, buffer);
            }
        }
        return;
    }
    if (eq(t.X[0], t.X[2])) {
        double r_m = (t.Y[1] - t.Y[2]) / (t.X[1] - t.X[2]);
        double r_b = t.Y[1] - r_m * t.X[1];
        //only move right
        double scan_XR;
        int right_end;
        int left_end = ceil_441(t.X[0]);
        //cerr << row_min << "\t" << row_max << endl;

        for (int i = row_min; i <= row_max; i++) {
            scan_XR = (i - r_b) / r_m;
            right_end = (int) floor_441(scan_XR);

            for (int j = left_end; j <= right_end; j++) {
                if (i < 0) break; // Don't draw below screen
                writePoint(j, i, t.color, buffer);
            }
        }
        return;
    }


    double l_m = (t.Y[2] - t.Y[0]) / (t.X[2] - t.X[0]);
    double r_m = (t.Y[1] - t.Y[2]) / (t.X[1] - t.X[2]);
    
    double l_b = t.Y[0] - l_m * t.X[0];
    double r_b = t.Y[1] - r_m * t.X[1];

    double scan_XL, scan_XR;
    int right_end, left_end;

    for (int i = row_min; i <= row_max; i++) {
        
        // Bound Checking left right
        scan_XR = (i - r_b) / r_m;
        scan_XL = (i - l_b) / l_m;

        right_end = (int) floor_441(scan_XR);
        left_end  = (int) ceil_441(scan_XL); 

        //Fill Buffer for a row
        for (int j = left_end; j <= right_end; j++) {
            if (i < 0) break; // Don't draw below screen
            writePoint(j, i, t.color, buffer);
        }
    }
}

void rasterTriangle(Triangle &t, unsigned char *buffer) {
    // Here we will first check if the triangle is degenerate
    // e.g. is the triangle just a line or perhaps just a point.
    if ((eq(t.X[0], t.X[1]) && eq(t.Y[0], t.Y[1])) ||
        (eq(t.X[1], t.X[2]) && eq(t.Y[1], t.Y[2])) ||
        (eq(t.X[2], t.X[0]) && eq(t.Y[2], t.Y[0]))) {
        
        return;
    }

    // Make sure that the triangle is on the screen at all
    
    //Make sure vertices are in a predicatable order
    t.sortVertices();
    
    //Here we will determine if the triangle is an uptriangle, 
    //a downtriangle, or neither, and treat it accordingly.
    if (eq(t.Y[1], t.Y[2])) rasterUpTriangle(t, buffer);
    else if (eq(t.Y[0], t.Y[1])) rasterDownTriangle(t, buffer);
    else {
        //cerr << "splitting" << endl;
        Triangle top, bottom;
        top.copy(t);
        bottom.copy(t);

        double new_x = t.X[0]
                    +((t.Y[1]-t.Y[0])/(t.Y[2]-t.Y[0]))
                    * (t.X[2]-t.X[0]);
        //double m = (t.Y[0] - t.Y[2]) / (t.X[0] - t.X[2]);
        //double b = (t.Y[0] - t.X[0] * m);
        //double old = (t.Y[0] - b) / m;
        //cerr << old_x << "\t" << new_x << endl;
        
        top.X[2] = new_x;
        top.Y[2] = top.Y[1];

        bottom.X[0] = new_x;
        bottom.Y[0] = bottom.Y[1];
       
        //cerr << "ANOTHER ONE" << endl;
        //top.print();
        //bottom.print();

        rasterTriangle(top, buffer);
        rasterTriangle(bottom, buffer);
    }
}


//END Code that I've added


int main()
{
   vtkImageData *image = NewImage(X_SCALE, Y_SCALE);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = X_SCALE*Y_SCALE;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = X_SCALE;
   screen.height = Y_SCALE;
   


   // Rasterize all triangles 
   for (Triangle t: triangles) {
        //cerr << "In" << endl;
        rasterTriangle(t, buffer);
        //cerr << "out" << endl;
    }

   WriteImage(image, "allTriangles");
}
