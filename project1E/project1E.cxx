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
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

using std::cerr;
using std::endl;

// Helper Functions

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

inline bool eq(double a, double b) {
    return a == b;  // To be changed if floating point accuracy issue
}

double lerp_d(const double a, const double b, const double c, const double d, const double x) {
    double t = (x - a) / (b - a);
    return c + t * (d - c);
}

double lerp_d_t(const double a, const double b, const double t) {
    return a + t * (b - a);
}

double *cross_3(const double *a, const double *b) {
    double *c = (double *) malloc(3 * sizeof(double));
    
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    
    return c;
    
}

double dot_3(const double *a, const double *b) {
    double total = 0.0;
    
    for (int i =0; i < 3; i++) {
        total += a[i] * b[i];
    }
    
    return total;
}

double * sub_3(const double *a, const double *b) {
    double *c = (double *) malloc(3 * sizeof(double));
    
    for (int i = 0; i < 3; i++) {
        c[i] = a[i] - b[i];
    }
    
    return c;
}

void normalize_3(double *v) {
    double length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    for (int i = 0; i < 3; i++) {
        v[i] /= length;
    }
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

/*
 *  Triangle Class
 */

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

    const void print() {
        cerr << "Triangle " << this << ":" << endl;
        for (int i = 0; i < 3; i++) {
            cerr << X[i] << "\t" << Y[i] << endl;
        }
    }
    
    //Place the vertices in CW order, v0 is top left, v1 is bottom/top right
    //v2 is bottom(left),
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

/*
 *  Screen Class
 */

class Screen
{
  public:
      unsigned char   *buffer;
      double          *z_buffer;
      int width, height;

};

/*
 *  Matrix Class
 */

class Matrix
{
public:
    double          A[4][4]= {{0.0, 0.0, 0.0, 0.0},
                              {0.0, 0.0, 0.0, 0.0},
                              {0.0, 0.0, 0.0, 0.0},
                              {0.0, 0.0, 0.0, 0.0}};  // A[i][j] means row i, column j
    
    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }
    
    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
    + ptIn[1]*A[1][0]
    + ptIn[2]*A[2][0]
    + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
    + ptIn[1]*A[1][1]
    + ptIn[2]*A[2][1]
    + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
    + ptIn[1]*A[1][2]
    + ptIn[2]*A[2][2]
    + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
    + ptIn[1]*A[1][3]
    + ptIn[2]*A[2][3]
    + ptIn[3]*A[3][3];
}

/*
 *  Camera Class
 */

class Camera
{
public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    
    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(void);
};

Matrix Camera::ViewTransform() {
    Matrix m;
    m.A[0][0] = 1.0 / tan(angle / 2.0);
    m.A[1][1] = 1.0 / tan(angle / 2.0);
    m.A[2][2] = (far + near) / (far - near);
    m.A[2][3] = -1.0;
    m.A[3][2] = (2.0 * far * near) / (far - near);
    return m;
}

Matrix Camera::CameraTransform() {
    Matrix m;
    double *w = sub_3(position, focus);
    normalize_3(w);
    double *u = cross_3(up, w);
    normalize_3(u);
    double *v = cross_3(w, u);
    normalize_3(v);
    
    double t[3] = {-position[0], -position[1], -position[2]};
    double udt = dot_3(u, t);
    double vdt = dot_3(v, t);
    double wdt = dot_3(w, t);
    
    m.A[0][0] = u[0];
    m.A[1][0] = u[1];
    m.A[2][0] = u[2];
    m.A[3][0] = udt;
    
    m.A[0][1] = v[0];
    m.A[1][1] = v[1];
    m.A[2][1] = v[2];
    m.A[3][1] = vdt;
    
    m.A[0][2] = w[0];
    m.A[1][2] = w[1];
    m.A[2][2] = w[2];
    m.A[3][2] = wdt;
    
    m.A[3][3] = 1.0;
    
    return m;
}

Matrix Camera::DeviceTransform() {
    Matrix m;
    m.A[0][0] = X_SCALE / 2.0;
    m.A[3][0] = X_SCALE / 2.0;
    
    m.A[1][1] = Y_SCALE / 2.0;
    m.A[3][1] = Y_SCALE / 2.0;
    
    m.A[2][2] = 1.0;
    m.A[3][3] = 1.0;
    
    return m;
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

//  Reader method
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
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
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
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

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
    
    // TODO replace the lerp_d with lerp_d_t
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

        rasterTriangle(top, screen);
        rasterTriangle(bottom, screen);
    }
}

void produceImage(int width, int height, std::vector<Triangle> triangles, int frame, int nframes) {
    vtkImageData *image = NewImage(width, height);
    int npixels = width * height;
    
    unsigned char *buffer =
    (unsigned char *) image->GetScalarPointer(0,0,0);
    double *z_buffer = new double[npixels];
    
    for (int i = 0 ; i < npixels*3 ; i++)
        buffer[i] = 0;
    for (int i = 0 ; i < npixels ; i++)
        z_buffer[i] = -INFINITY;
    
    
    // Apply Camera
    Camera c = GetCamera(frame, nframes);
    
    //DEBUG 
    /*
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 0;
    c.position[1] = 40;
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    */
    Matrix m = Matrix::ComposeMatrices(c.CameraTransform(), c.ViewTransform());
    m = Matrix::ComposeMatrices(m, c.DeviceTransform());
    std::vector<Triangle> new_triangles;
    
    /* //DEBUG
    c.DeviceTransform().Print(cerr);
    c.ViewTransform().Print(cerr);
    c.CameraTransform().Print(cerr);
    m.Print(cerr);
    */
   
    int i = 0; //DEBUG
    
    for (Triangle t: triangles) {
        Triangle new_triangle(t);
        for (int j = 0; j < 3; j++) {
            const double old_pt[4] = {t.X[j], t.Y[j], t.Z[j], 1};
            //const double old_pt[4] = {0, 36.4645, 36.4645, 1}; //DEBUG
            double       new_pt[4];
            m.TransformPoint(old_pt, new_pt);
            
            //DEBUG
            /*
            cerr << "Tranformed " << old_pt[0] << ", " 
                << old_pt[1] << ", " 
                << old_pt[2] << " to "
                << new_pt[0]/new_pt[3] << ", "
                << new_pt[1]/new_pt[3] << ", "
                << new_pt[2]/new_pt[3] << endl;
            if ( ++i > 1) {
                return;
            }
            */

            new_triangle.X[j] = new_pt[0] / new_pt[3];
            new_triangle.Y[j] = new_pt[1] / new_pt[3];
            new_triangle.Z[j] = new_pt[2] / new_pt[3];
        }
        new_triangles.push_back(new_triangle);
    }
    
    
    Screen screen;
    screen.buffer = buffer;
    screen.z_buffer = z_buffer;
    screen.width = width;
    screen.height = height;
    
    // Rasterize all triangles
    for (Triangle t: new_triangles) {
        rasterTriangle(t, screen);
    }
    char buff[100];
    snprintf(buff, sizeof(buff), "frame%.3d", frame);
    WriteImage(image, buff);
}


int main()
{
    std::vector<Triangle> triangles = GetTriangles();
    for (int i = 0; i < 1000; i += 250) {
        produceImage(X_SCALE, Y_SCALE, triangles, i, 1000);
    }
    //produceImage(X_SCALE, Y_SCALE, triangles, 0, 1000);
    
    /*
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
        rasterTriangle(t, screen);
    }

   WriteImage(image, "allTriangles");
     */
    
    
}
