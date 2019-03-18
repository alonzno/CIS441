/**
 * CIS 541 Winter 2019
 * Project G
 * He He
 */

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <cmath>

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

double mag_3(double *v) {
    double length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    return length;
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

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

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

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    // use your own project 1 code for those three functions
    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(int n, int m);
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

Matrix Camera::DeviceTransform(int n, int m) {
    Matrix M;
    M.A[0][0] = n / 2.0;
    M.A[3][0] = n / 2.0;
    
    M.A[1][1] = m / 2.0;
    M.A[3][1] = m / 2.0;
    
    M.A[2][2] = 1.0;
    M.A[3][3] = 1.0;
    
    return M;
}

// note that this is different from project 1 starter code
Camera
GetCamera(int frame, int nframes)
{
    double t = (double)frame/(double)nframes;
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 80 * sin(2*M_PI*t);
    c.position[1] = 30;
    c.position[2] = 80 * cos(2*M_PI*t);
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

class Screen
{
    public:
        unsigned char   *buffer;
        double          *zbuffer;
        int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

void colorPixel(Screen * screen, int r, int c, unsigned char color[3], double depth){
    // implement this function, or use your own project 1E code.
    if (r < 0 || c < 0 || r >= screen -> height || c >= screen -> width) return;
    if (depth > screen -> zbuffer[r * screen -> width + c]) {
        //cerr << "Writing " << r << " " << c << endl;
        //cerr << 1+color[0] << " " << 1+color[1] << " "  << 1+color[2] << endl;
        screen -> buffer[3 * (r * screen -> width + c) + 0] = color[0];
        screen -> buffer[3 * (r * screen -> width + c) + 1] = color[1];
        screen -> buffer[3 * (r * screen -> width + c) + 2] = color[2];
    }
}

class Line{
public:
    double X[2];
    double Y[2];
    double Z[2];
    unsigned char   color[3];

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }
    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Draw(Screen * screen){
        int x0, y0, x1, y1;
        if(Y[0] < Y[1]){
            x0 = X[0];
            y0 = Y[0];
            x1 = X[1];
            y1 = Y[1];
        }else{
            x0 = X[1];
            y0 = Y[1];
            x1 = X[0];
            y1 = Y[0];
        }

        int e;
        int dx = x1 - x0;
        int dy = y1 - y0;
        int x = x0;
        int y = y0;
        double dz = Z[1] - Z[0];

        int stepx = 1;
        int stepy = 1;
        if (dx<0){
            dx = - dx;
            stepx = -1;
        }
        if (dy<0){
            dy = - dy;
            stepy = -1;
        }

        int i;
        // first octant
        if(dy<dx){
            e = 2 * dy - dx;
            for(i=0; i<dx; i++){
                colorPixel(screen, y, x, color, i*dz/dx+Z[0]);
                while(e > 0){
                    y += stepy;
                    e = e - 2 * dx;
                }
                x += stepx;
                e = e + 2 * dy;
            }
        // second octant
        }else{
            e = 2 * dx - dy;
            for(i=0; i<dy; i++){
                colorPixel(screen, y, x, color, i*dz/dy+Z[0]);
                while(e > 0){
                    x += stepx;
                    e = e - 2 * dy;
                }
                y += stepy;
                e = e + 2 * dx;
            }
        }
    }

    void Print(){
        std::cout << "v1: " << X[0] << ", " << Y[0] << ", " << Z[0] << '\n';
        std::cout << "v2: " << X[1] << ", " << Y[1] << ", " << Z[1] << '\n';
    }
};

void Bezier_Divide(double pX[4], double pY[4], double pZ[4], double qX[4], double qY[4], double qZ[4], double rX[4], double rY[4], double rZ[4]){
    // your code goes here
    qX[0] = pX[0];
    qY[0] = pY[0];
    qZ[0] = pZ[0];
    
    qX[1] = 0.5 * (pX[0] + pX[1]);
    qY[1] = 0.5 * (pY[0] + pY[1]);
    qZ[1] = 0.5 * (pZ[0] + pZ[1]);

    qX[2] = 0.5 * qX[1] + 0.25 * (pX[1] + pX[2]);
    qY[2] = 0.5 * qY[1] + 0.25 * (pY[1] + pY[2]);
    qZ[2] = 0.5 * qZ[1] + 0.25 * (pZ[1] + pZ[2]);

    rX[3] = pX[3];
    rY[3] = pY[3];
    rZ[3] = pZ[3];

    rX[2] = 0.5 * (pX[2] + pX[3]);
    rY[2] = 0.5 * (pY[2] + pY[3]);
    rZ[2] = 0.5 * (pZ[2] + pZ[3]);

    rX[1] = 0.5 * rX[2] + 0.25 * (pX[1] + pX[2]);
    rY[1] = 0.5 * rY[2] + 0.25 * (pY[1] + pY[2]);
    rZ[1] = 0.5 * rZ[2] + 0.25 * (pZ[1] + pZ[2]);

    qX[3] = 0.5 * (qX[2] + rX[1]);
    qY[3] = 0.5 * (qY[2] + rY[1]);
    qZ[3] = 0.5 * (qZ[2] + rZ[1]);

    rX[0] = qX[3];
    rY[0] = qY[3];
    rZ[0] = qZ[3];

    return;
}

class BezierCurve{
public:
    double X[4];
    double Y[4];
    double Z[4];
    unsigned char color[3];

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }
    
    bool eval() {
        
        return sqrt((X[3]-X[0]) * (X[3]-X[0]) +
                    (Y[3]-Y[0]) * (Y[3]-Y[0]) +
                    (Z[3]-Z[0]) * (Z[3]-Z[0]) ) < 100;
        
        double x1[3] = {X[0], Y[0], Z[0]};
        double x2[3] = {X[3], Y[3], Z[3]};
        
        double t0[3] = {X[1], Y[1], Z[1]};
        double t1[3] = {X[2], Y[2], Z[2]};
        
        double d = mag_3(cross_3(sub_3(t0,x1), sub_3(t0, x2))) / mag_3(sub_3(x2,x1));
               d += mag_3(cross_3(sub_3(t1,x1), sub_3(t1, x2))) / mag_3(sub_3(x2,x1));
        
        double EPSILON = 10;
        
        return d < EPSILON;
    }

    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Print(){
        printf("p0: %lf, %lf, %lf\n", X[0], Y[0], Z[0]);
        printf("p1: %lf, %lf, %lf\n", X[1], Y[1], Z[1]);
        printf("p2: %lf, %lf, %lf\n", X[2], Y[2], Z[2]);
        printf("p3: %lf, %lf, %lf\n", X[3], Y[3], Z[3]);
    }

    void RotateY(BezierCurve * dst, double angleDegree){
        int i;
        for(i=0; i<4; i++){
            dst->X[i] = X[i] * cos(2*M_PI*angleDegree/360)
                        + Z[i] * sin(2*M_PI*angleDegree/360);
            dst->Y[i] = Y[i];
            dst->Z[i] = Z[i] * cos(2*M_PI*angleDegree/360)
                        - X[i] * sin(2*M_PI*angleDegree/360);
        }
    }

    void Draw(Screen * screen){
        //cerr << "Writing" << endl;
        if (eval()) {
            for (int i = 0; i < 3; i++) {
                Line l;
                l.SetPointByIndex(0, X[i], Y[i], Z[i]);
                l.SetPointByIndex(1, X[i+1], Y[i+1], Z[i+1]);
                l.SetColor(color[0], color[1], color[2]);
                l.Draw(screen);
            }
        }
        else {
            BezierCurve c0;
            BezierCurve c1;
            
            Bezier_Divide(X, Y, Z,
                          c0.X, c0.Y, c0.Z,
                          c1.X, c1.Y, c1.Z);
            c0.SetColor(color[0], color[1], color[2]);
            c1.SetColor(color[0], color[1], color[2]);
            
            c0.Draw(screen);
            c1.Draw(screen);
        }
    }
};

class BezierSurface{
public:
    double X[16];
    double Y[16];
    double Z[16];
    unsigned char color[3];

    void SetPoints(double x[16], double y[16], double z[16]){
        int i;
        for(i=0; i<16; i++){
            X[i] = x[i];
            Y[i] = y[i];
            Z[i] = z[i];
        }
    }

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }

    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Print(){
        printf("v0: %lf, %lf, %lf\n", X[0], Y[0], Z[0]);
        printf("v1: %lf, %lf, %lf\n", X[1], Y[1], Z[1]);
        printf("v2: %lf, %lf, %lf\n", X[2], Y[2], Z[2]);

    }

    void Draw(Screen * screen, int divisions){
        // your code goes here
        if (divisions == 0) {
            //cerr << "GOing 0" << endl;
            for (int i = 0; i < 4; i++) {
                BezierCurve b;
                b.SetColor(color[0], color[1], color[2]);
                for (int j = 0; j < 4; j++) {
                    //cerr << "GOo1" << endl;
                    //cerr << X[i*4+j] << " " <<  Y[i*4+j] << " " << Z[i*4+j] << endl;
                    b.SetPointByIndex(j, X[i*4+j], Y[i*4+j], Z[i*4+j]);
                }
                //cerr << "exit" << endl;
                b.Draw(screen);
                //cerr << "exita" << endl;
            }
            for (int i = 0; i < 4; i++) {
                BezierCurve b;
                b.SetColor(color[0], color[1], color[2]);
                for (int j = 0; j < 4; j++) {
                    //cerr << "GOo2" << endl;
                    b.SetPointByIndex(j, X[j*4+i], Y[j*4+i], Z[j*4+i]);
                }
                b.Draw(screen);
            }
        }
        else {
            //cerr << "GOing 1" << endl;
            BezierSurface b[4];
            for (int i = 0; i < 4; i++) {
                b[i].SetColor(color[0], color[1], color[2]);
            }
            double nX[4][8], nY[4][8], nZ[4][8];
            for (int i = 0; i < 4; i++) {
                Bezier_Divide(&(X[i*4]), &(Y[i*4]), &(Z[i*4]),
                              &(nX[i][0]),  &(nY[i][0]), &(nZ[i][0]),
                              &(nX[i][0])+8, &(nY[i][0])+8, &(nZ[i][0])+8);
            }
            double oX[8][4], oY[8][4], oZ[8][4];
            for(int i=0; i < 4; ++i)
                for(int j=0; j < 8; ++j) {
                    oX[j][i] = nX[i][j];
                    oY[j][i] = nY[i][j];
                    oZ[j][i] = nZ[i][j];
            }
            
            for (int i = 0; i < 4; i++) {
                Bezier_Divide(&(oX[i][0]), &(oY[i][0]), &(oZ[i][0]),
                              &(b[0].X[4*i]), &(b[0].Y[4*i]), &(b[0].Z[4*i]),
                              &(b[1].X[4*i]), &(b[1].Y[4*i]), &(b[1].Z[4*i]));
            }
            for (int i = 4; i < 8; i++) {
                Bezier_Divide(&(oX[i][0]), &(oY[i][0]), &(oZ[i][0]),
                              &(b[2].X[4*(i-4)]), &(b[2].Y[4*(i-4)]), &(b[2].Z[4*(i-4)]),
                              &(b[3].X[4*(i-4)]), &(b[3].Y[4*(i-4)]), &(b[3].Z[4*(i-4)]));
            }
            for (int i = 0; i < 4; i++) {
                b[i].Draw(screen, divisions - 1);
            }
            
        }

        return;
    }
};

// returns an array of BezierCurves of size 6
BezierCurve * getLeaves(){
    BezierCurve *c = (BezierCurve*)malloc(sizeof(BezierCurve)*6);
    int i;
    for(i=0; i<6; i++)
        c[i].SetColor(30, 215, 97);

    c[0].SetPointByIndex(0, 0, -10, 0);
    c[0].SetPointByIndex(1, 14, -7.2, 0);
    c[0].SetPointByIndex(2, 9.8, -3, 2.8);
    c[0].SetPointByIndex(3, 14, 4, 14);

    c[1].SetPointByIndex(0, 0, -10, 0);
    c[1].SetPointByIndex(1, 0, -7.2, 14);
    c[1].SetPointByIndex(2, 2.8, -3, 9.8);
    c[1].SetPointByIndex(3, 14, 4, 14);

    c[0].RotateY(&(c[2]), 120);
    c[1].RotateY(&(c[3]), 120);
    c[0].RotateY(&(c[4]), 240);
    c[1].RotateY(&(c[5]), 240);

    return &(c[0]);
}

// returns an array of BezierSurfaces of size 2
BezierSurface * getSurfaces(){
    BezierSurface *s = (BezierSurface*)malloc(sizeof(BezierSurface)*2);

    s[0].SetPointByIndex(0, 0.0, 0.0, 0.0); // first row
	s[0].SetPointByIndex(1, 0.0, 3, 5);
	s[0].SetPointByIndex(2, 0.0, 7.5, 10);
	s[0].SetPointByIndex(3, 0.0, 1.5, 15.0);
	s[0].SetPointByIndex(4, 5, 12, 0.0); // second row
	s[0].SetPointByIndex(5, 5, -1.5, 5);
	s[0].SetPointByIndex(6, 5, 0.0, 10);
	s[0].SetPointByIndex(7, 5, 1.5, 15.0);
	s[0].SetPointByIndex(8, 10, 4.5, 0.0); // third row
	s[0].SetPointByIndex(9, 10, 12, 5);
	s[0].SetPointByIndex(10, 10, 13.5, 10);
	s[0].SetPointByIndex(11, 10, 7.5, 15.0);
	s[0].SetPointByIndex(12, 15.0, 6, 0.0); // fourth row
	s[0].SetPointByIndex(13, 15.0, 3, 5);
	s[0].SetPointByIndex(14, 15.0, 7.5, 10);
	s[0].SetPointByIndex(15, 15.0, 15.0, 15.0);
    s[0].SetColor(51, 133, 229);

    s[1].SetPointByIndex(0, 0.0, -3, 0.0); // first row
	s[1].SetPointByIndex(1, 0.0, -3, 5);
	s[1].SetPointByIndex(2, 0.0, -3, 10);
	s[1].SetPointByIndex(3, 0.0, -3, 15);
	s[1].SetPointByIndex(4, 5, -3, 0.0); // second row
	s[1].SetPointByIndex(5, 5, -3, 5);
	s[1].SetPointByIndex(6, 5, -3, 10);
	s[1].SetPointByIndex(7, 5, -3, 15);
	s[1].SetPointByIndex(8, 10, -3, 0.0); // third row
	s[1].SetPointByIndex(9, 10, -3, 5);
	s[1].SetPointByIndex(10, 10, -3, 10);
	s[1].SetPointByIndex(11, 10, -3, 15);
	s[1].SetPointByIndex(12, 15, -3, 0.0); // fourth row
	s[1].SetPointByIndex(13, 15, -3, 5);
	s[1].SetPointByIndex(14, 15, -3, 10);
	s[1].SetPointByIndex(15, 15, -3, 15);
    s[1].SetColor(31, 179, 83);

    return &(s[0]);
}

int main()
{
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer =
        (unsigned char *) image->GetScalarPointer(0,0,0);

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;
    screen.zbuffer = (double*)malloc(sizeof(double) * screen.width * screen.height);

    // uncooment the following two lines to get the curves and surfaces
    BezierCurve *c;
    BezierSurface *s = getSurfaces();

    // Camera rotates around y-axis
    // It comes back to the original position after 100 iterations
    int i;
    for(i=0; i<100; i++){
        int j;
        int npixels = screen.width * screen.height;
        for (int j = 0 ; j < npixels*3 ; j++)
            screen.buffer[j] = 0;
        for(int j = 0; j < npixels; j++)
            screen.zbuffer[j] = -1;
        
        c = getLeaves();
        s = getSurfaces();
        
        Camera camera = GetCamera(i*10, 1000);
        Matrix ct = camera.CameraTransform();
        Matrix vt = camera.ViewTransform();
        Matrix intermediate = Matrix::ComposeMatrices(ct, vt);
        Matrix dt = camera.DeviceTransform(screen.width, screen.height);
        Matrix total = Matrix::ComposeMatrices(intermediate, dt);
        
        bool leaves = true;
        if (!leaves) {
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 16; j++) {
                    //cerr <<"here" << endl;
                    double in[4] = {s[i].X[j], s[i].Y[j], s[i].Z[j], 1.0};
                    double out[4];
                    total.TransformPoint(in, out);
                    //cerr <<"there" << endl;
                    s[i].SetPointByIndex(j, out[0]/out[3], out[1]/out[3], out[2]/out[3]);
                    //cerr << out[0]/out[3] << " " << out[1]/out[3]  << " " << out[2]/out[3] << " "  << out[3]/out[3] << endl;

                }
            }
            for (int i = 0 ; i < 2; i++) {
                s[i].Draw(&screen, 2);
                cerr << "Writtne"<<endl;
            }
        }
        else {
            // draw your curves and surfaces here
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 4; j++) {
                    double in[4] = {c[i].X[j], c[i].Y[j], c[i].Z[j], 1.0};
                    double out[4];

                    total.TransformPoint(in, out);

                    //cerr << out[0]/out[3] << " " << out[1]/out[3]  << " " << out[2]/out[3] << " "  << out[3]/out[3] << endl;

                    c[i].SetPointByIndex(j, out[0]/out[3], out[1]/out[3], out[2]/out[3]);
                }
            }
            for (int i = 0 ; i < 6; i++) {
                c[i].Draw(&screen);
            }
        }

        char name[20];
        // make sure you have a directory named images
        // so that writing images won't cause you errors
        sprintf(name, "./images/frame%02d", i);
        WriteImage(image, name);
        //cerr << "Written" << endl;

    }
}
