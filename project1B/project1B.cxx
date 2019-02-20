#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

#include <math.h>

#define X_SCALE 1000
#define Y_SCALE 1000

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

  // would some methods for transforming the triangle in place be helpful?
    void print() {
        cerr << "Triangle " << this << ":" << endl;
        for (int i = 0; i < 3; i++) {
            cerr << X[i] << "\t" << Y[i] << endl;
        }
    }
    
    //Place the vertices in CW order, v1 is top left, v2 is bottom(left),
    // v3 is bottom/top right
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
        if (Y[0] == Y[1]) {
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
        if (Y[1] == Y[2]) {
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
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       if (i == 50)
           rv[i].X[firstPt] = -10;
       rv[i].Y[firstPt] = posJ+10*(idxJ+1);
       rv[i].X[(firstPt+1)%3] = posI+99;
       rv[i].Y[(firstPt+1)%3] = posJ+10*(idxJ+1);
       rv[i].X[(firstPt+2)%3] = posI+i;
       rv[i].Y[(firstPt+2)%3] = posJ;
       if (i == 5)
          rv[i].Y[(firstPt+2)%3] = -50;
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;
}

// Code that I've added
void printDebug(std::vector<Triangle> tris) {
    int c = 0;
    for (Triangle t: tris) {
        std::cout << "Triangle " << c << ": " << std::endl;
        c++;

        for (int i = 0 ; i < 3; i++) {
            std::cout << t.X[i] << "\t" << t.Y[i] << std::endl;
        }
        std::cout << "R: " << (int) t.color[0] << "\t" << std::endl;
        std::cout << "G: " << (int) t.color[1] << "\t" << std::endl;
        std::cout << "B: " << (int) t.color[2] << "\t" << std::endl;
        std::cout << std::endl;
    }
}

void rasterTriangle(Triangle t, unsigned char *buffer) {
    // Here we will first check if the triangle is degenerate
    // e.g. is the triangle just a line or perhaps just a point.
    
    //Here we will determine if the triangle is an uptriangle, 
    //a downtriangle, or neither, and treat it accordingly.
}

void rasterUpTriangle(Triangle t, unsigned char *buffer){
    double row_min = ceil_441(t.Y[2]);
    double row_max = floor_441(t.Y[0]);
    
    row_max = fmin(Y_SCALE - 1, row_max);   // Don't Draw Above Screen

    //Avoids divide by zero
    double left_slope = (t.X[0] - t.X[1]) / (t.Y[0] - t.Y[1]);
    double right_slope = (t.X[2] - t.X[0]) / (t.Y[2] - t.Y[0]);
    
    double scan_XL = t.X[0];
    double scan_XR = t.X[0];

    int loc;
    int right_end, left_end;

    for (int i = row_max; i >= row_min; i--) {
        
        // Bound Checking left right
        right_end = (int) fmin(X_SCALE-1, floor_441(scan_XR));
        left_end  = (int) fmax(0.0, ceil_441(scan_XL)); 

        //Fill Buffer for a row
        for (int j = left_end; j <= right_end; j++) {
            if (i < 0) break; // Don't draw below screen

            loc = 3*i*X_SCALE + 3*j;
            buffer[loc + 0] = t.color[0];
            buffer[loc + 1] = t.color[1];
            buffer[loc + 2] = t.color[2];
        }
        scan_XL -= left_slope;
        scan_XR += right_slope;
    }
}

void rasterDownTriangle(Triangle t, unsigned char *buffer) {
    double row_min = ceil_441(t.Y[2]);
    double row_max = floor_441(t.Y[0]);
    
    row_max = fmin(Y_SCALE - 1, row_max);   // Don't Draw Above Screen

    //Avoids divide by zero
    double left_slope = (t.X[2] - t.X[0]) / (t.Y[2] - t.Y[0]);
    double right_slope = (t.X[1] - t.X[2]) / (t.Y[1] - t.Y[2]);
    
    double scan_XL = t.X[2];
    double scan_XR = t.X[2];

    int loc;
    int right_end, left_end;

    for (int i = row_min; i <= row_max; i++) {
        
        // Bound Checking left right
        right_end = (int) fmin(X_SCALE-1, floor_441(scan_XR));
        left_end  = (int) fmax(0.0, ceil_441(scan_XL)); 

        //Fill Buffer for a row
        for (int j = left_end; j <= right_end; j++) {
            if (i < 0) break; // Don't draw below screen

            loc = 3*i*X_SCALE + 3*j;
            buffer[loc + 0] = t.color[0];
            buffer[loc + 1] = t.color[1];
            buffer[loc + 2] = t.color[2];
        }
        scan_XL += left_slope;
        scan_XR += right_slope;
    }
}

//END Code that I've added


int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   // Rasterize all triangles 
   for (Triangle t: triangles) {
        t.sortVertices();
        rasterDownTriangle(t, buffer);
    }

   WriteImage(image, "allTriangles");
}
