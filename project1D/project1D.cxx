#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

#include <math.h>

#define X_SCALE 1000
#define Y_SCALE 1000

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
      double    X[3];
      double    Y[3];
      double    Z[3];
      double    colors[3][3];

    void copy(const Triangle &other) {
        memcpy(X, other.X, 3 * sizeof(double));
        memcpy(Y, other.Y, 3 * sizeof(double));
        memcpy(Z, other.Z, 3 * sizeof(double));
        memcpy(colors[0], other.colors[0], 3 * sizeof(double));
        memcpy(colors[1], other.colors[1], 3 * sizeof(double));
        memcpy(colors[2], other.colors[2], 3 * sizeof(double));
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
    void swap(int i, int j) {
        double temp, temp_c[3];
        temp = X[j];
        X[j] = X[i];
        X[i] = temp;
        
        temp = Y[j];
        Y[j] = Y[i];
        Y[i] = temp;
        
        temp = Z[j];
        Z[j] = Z[i];
        Z[i] = temp;
        
        memcpy(temp_c, colors[j], 3 * sizeof(double)); 
        memcpy(colors[j], colors[i], 3 * sizeof(double)); 
        memcpy(colors[i], temp_c, 3 * sizeof(double)); 
        
    }
    void sortVertices() {
        double temp;
        if (Y[0] < Y[1]) {
            swap(0, 1);
        }
        
        if (Y[1] < Y[2]) {
            swap(1, 2);
        }
        
        if (Y[0] < Y[1]) {
            swap(0, 1);
        }
        
        //Case of downward tri
        if (eq(Y[0], Y[1])) {
            if (X[0] > X[1]) {
                swap(0, 1);
            }
            /*  0 1
             *   2
             */
        }
       
        //Case of upward tri
        if (eq(Y[1], Y[2])) {
            if (X[1] < X[2]) {
                swap(1, 2);
            }
            /*   0 
             *  2 1
             */
        }
    }
};


double lerp_d(double a, double b, double c, double d, double x) {
    double t = (x - a) / (b - a);
    return c + t * (d - c);
}

class Screen
{
  public:
      unsigned char   *buffer;
      double          *z_buffer;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
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
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
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
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}
/* std::vector<Triangle>
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
*/

// Code that I've added

void writePoint(int x, int y, const Triangle &t, Screen &screen, double x_l, double x_r, bool up) {
    // Perform bounds checking
    if (x < 0 || x >= X_SCALE || y < 0 || y >= Y_SCALE) return;
     
    int loc = 3 * (y*X_SCALE + x);
    int loc_z = (y*X_SCALE + x);

    double t_l, t_r, r, g, b, z;

    /*  0 1
     *   2
     */
    /*   0 
     *  2 1
     */

    if (up) {
        t_l = lerp_d(t.Y[0], t.Y[2], t.colors[0][0], t.colors[2][0], y); 
        t_r = lerp_d(t.Y[0], t.Y[1], t.colors[0][0], t.colors[1][0], y); 
        r   = lerp_d(x_l, x_r, t_l, t_r, x);
        
        t_l = lerp_d(t.Y[0], t.Y[2], t.colors[0][1], t.colors[2][1], y); 
        t_r = lerp_d(t.Y[0], t.Y[1], t.colors[0][1], t.colors[1][1], y); 
        g   = lerp_d(x_l, x_r, t_l, t_r, x);
        
        t_l = lerp_d(t.Y[0], t.Y[2], t.colors[0][2], t.colors[2][2], y); 
        t_r = lerp_d(t.Y[0], t.Y[1], t.colors[0][2], t.colors[1][2], y); 
        b   = lerp_d(x_l, x_r, t_l, t_r, x);
        
        t_l = lerp_d(t.Y[0], t.Y[2], t.Z[0], t.Z[2], y); 
        t_r = lerp_d(t.Y[0], t.Y[1], t.Z[0], t.Z[1], y); 
        z   = lerp_d(x_l, x_r, t_l, t_r, x);
    }
    else {
        t_l = lerp_d(t.Y[0], t.Y[2], t.colors[0][0], t.colors[2][0], y); 
        t_r = lerp_d(t.Y[1], t.Y[2], t.colors[1][0], t.colors[2][0], y); 
        r   = lerp_d(x_l, x_r, t_l, t_r, x);
        
        t_l = lerp_d(t.Y[0], t.Y[2], t.colors[0][1], t.colors[2][1], y); 
        t_r = lerp_d(t.Y[1], t.Y[2], t.colors[1][1], t.colors[2][1], y); 
        g   = lerp_d(x_l, x_r, t_l, t_r, x);
        
        t_l = lerp_d(t.Y[0], t.Y[2], t.colors[0][2], t.colors[2][2], y); 
        t_r = lerp_d(t.Y[1], t.Y[2], t.colors[1][2], t.colors[2][2], y); 
        b   = lerp_d(x_l, x_r, t_l, t_r, x);
        
        t_l = lerp_d(t.Y[0], t.Y[2], t.Z[0], t.Z[2], y); 
        t_r = lerp_d(t.Y[1], t.Y[2], t.Z[1], t.Z[2], y); 
        z   = lerp_d(x_l, x_r, t_l, t_r, x);
    }
    
    if (screen.z_buffer[loc_z] < z) {
        screen.z_buffer[loc_z] = z;
        screen.buffer[loc + 0] = ceil_441(r * 255);
        screen.buffer[loc + 1] = ceil_441(g * 255);
        screen.buffer[loc + 2] = ceil_441(b * 255);
    }
}

void rasterUpTriangle(const Triangle &t, Screen &screen){
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
                writePoint(j, i, t, screen, scan_XL, t.X[1], true);
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
                writePoint(j, i, t, screen, t.X[0], scan_XR, true);
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
            writePoint(j, i, t, screen, scan_XL, scan_XR, true);
        }
    }
}

void rasterDownTriangle(const Triangle &t, Screen &screen) {
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
                writePoint(j, i, t, screen, scan_XL, t.X[1], false);
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
                writePoint(j, i, t, screen, t.X[0], scan_XR, false);
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
            writePoint(j, i, t, screen, scan_XL, scan_XR, false);
        }
    }
}

void rasterTriangle(Triangle &t, Screen &screen) {
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
    if (eq(t.Y[1], t.Y[2])) rasterUpTriangle(t, screen);
    else if (eq(t.Y[0], t.Y[1])) rasterDownTriangle(t, screen);
    else {
        Triangle top, bottom;
        top.copy(t);
        bottom.copy(t);

        double new_x = t.X[0]
                    +((t.Y[1]-t.Y[0])/(t.Y[2]-t.Y[0]))
                    * (t.X[2]-t.X[0]);
        double new_z = lerp_d(t.Y[0], t.Y[2], t.Z[0], t.Z[2], t.Y[1]);
        
        top.X[2] = new_x;
        top.Y[2] = top.Y[1];
        top.Z[2] = new_z;

        top.colors[2][0] = lerp_d(t.Y[0], t.Y[2], t.colors[0][0], t.colors[2][0], t.Y[1]);
        top.colors[2][1] = lerp_d(t.Y[0], t.Y[2], t.colors[0][1], t.colors[2][1], t.Y[1]);
        top.colors[2][2] = lerp_d(t.Y[0], t.Y[2], t.colors[0][2], t.colors[2][2], t.Y[1]);

        bottom.X[0] = new_x;
        bottom.Y[0] = bottom.Y[1];
        bottom.Z[0] = new_z;
       
        bottom.colors[0][0] = lerp_d(t.Y[0], t.Y[2], t.colors[0][0], t.colors[2][0], t.Y[1]);
        bottom.colors[0][1] = lerp_d(t.Y[0], t.Y[2], t.colors[0][1], t.colors[2][1], t.Y[1]);
        bottom.colors[0][2] = lerp_d(t.Y[0], t.Y[2], t.colors[0][2], t.colors[2][2], t.Y[1]);
        //cerr << "ANOTHER ONE" << endl;
        //top.print();
        //bottom.print();

        rasterTriangle(top, screen);
        rasterTriangle(bottom, screen);
    }
}


//END Code that I've added


int main()
{
   vtkImageData *image = NewImage(X_SCALE, Y_SCALE);
   int npixels = X_SCALE*Y_SCALE;

   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   double *z_buffer = new double[npixels];

   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;
   for (int i = 0 ; i < npixels ; i++)
       z_buffer[i] = -1.0;


   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.z_buffer = z_buffer;
   screen.width = X_SCALE;
   screen.height = Y_SCALE;
   


   // Rasterize all triangles 
   for (Triangle t: triangles) {
        //cerr << "In" << endl;
        rasterTriangle(t, screen);
        //cerr << "out" << endl;
    }

   WriteImage(image, "allTriangles");
}
