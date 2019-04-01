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
#include <math.h>


#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <fstream>
#include <stdlib.h>     /* atof */

using std::cerr;
using std::endl;

vtkImageData* NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

void normalize(double* v, double* result, int size){
    double m = 0;
    for (int i=0; i<size; i++){
        m = m + v[i] * v[i];
    }
    m = sqrt(m);
    for (int j=0; j<size; j++){
        result[j] = v[j] / m;
    }
}

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
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

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
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


    double          U[3];
    double          V[3];
    double          W[3];
    double          T[3];

    // use your own project 1 code for those three functions
    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(int n, int m);

    void            config(){

        // u = up x (O-focus) => up x w
        // v = (O-focus) x u  => w x u
        // w = (O-focus)
        // O

        for (int i=0; i<3; i++){
            W[i] = position[i] - focus[i];
            T[i] = - position[i];
        }

        normalize(W, W, 3);

        U[0] = up[1] * W[2] - up[2] * W[1];
        U[1] = up[2] * W[0] - up[0] * W[2];
        U[2] = up[0] * W[1] - up[1] * W[0];
        normalize(U, U, 3);

        V[0] = W[1] * U[2] - W[2] * U[1];
        V[1] = W[2] * U[0] - W[0] * U[2];
        V[2] = W[0] * U[1] - W[1] * U[0];
        normalize(V, V, 3);

    }
};

Matrix Camera::ViewTransform(){
    Matrix m;
    for (int i = 0 ; i < 4 ; i++){
        for (int j = 0 ; j < 4 ; j++)
        {
            m.A[i][j] = 0;
        }
    }
    m.A[0][0] = 1 / tan(angle / 2);
    m.A[1][1] = 1 / tan(angle / 2);
    m.A[2][2] = (far + near) / (far - near);
    m.A[2][3] = -1;
    m.A[3][2] = 2 * far * near / (far - near);

    return m;
}

Matrix Camera::CameraTransform(void) {
    Matrix m;

    double ut = 0;
    double vt = 0;
    double wt = 0;
    for (int i=0; i<3; i++){
        ut = ut + U[i] * T[i];
        vt = vt + V[i] * T[i];
        wt = wt + W[i] * T[i];
    }

    m.A[0][0] = U[0];
    m.A[1][0] = U[1];
    m.A[2][0] = U[2];
    m.A[3][0] = ut;

    m.A[0][1] = V[0];
    m.A[1][1] = V[1];
    m.A[2][1] = V[2];
    m.A[3][1] = vt;

    m.A[0][2] = W[0];
    m.A[1][2] = W[1];
    m.A[2][2] = W[2];
    m.A[3][2] = wt;

    m.A[0][3] = 0;
    m.A[1][3] = 0;
    m.A[2][3] = 0;
    m.A[3][3] = 1;

    return m;
}

Matrix Camera::DeviceTransform(int n, int m) {
    Matrix mat;

    for (int i = 0 ; i < 4 ; i++){
        for (int j = 0 ; j < 4 ; j++)
        {
            mat.A[i][j] = 0;
        }
    }
    mat.A[0][0] = n/2;
    mat.A[1][1] = m/2;
    mat.A[2][2] = 1;
    mat.A[3][0] = n/2;
    mat.A[3][1] = m/2;
    mat.A[3][3] = 1;

    return mat;
}

// note that this is different from project 1 starter code
Camera GetCamera(int frame, int nframes)
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

    if (r > screen->width || r <0){
        return ;
    }
    if (c > screen->height || c < 0){
        return ;
    }
    
    double zval = screen->zbuffer[r*screen->height + c];

    if (depth > zval) {
        int pixelval = (r * screen->height + c ) *3;
        screen->buffer[pixelval] = color[0];
        screen->buffer[pixelval+1] = color[1];
        screen->buffer[pixelval+2] = color[2];
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


void casteljau(double p[4], double q[4], double r[4]){

    q[0] = p[0];
    r[3] = p[3];
    q[1] = 0.5 * (p[0] + p[1]);
    r[2] = 0.5 * (p[2] + p[3]);
    q[2] = 0.5 * q[1] + 0.25 * (p[1] + p[2]);
    r[1] = 0.5 * r[2] + 0.25 * (p[1] + p[2]);
    q[3] = 0.5 * (q[2] + r[1]);
    r[0] = q[3];

}

void Bezier_Divide(double pX[4], double pY[4], double pZ[4], double qX[4], double qY[4], double qZ[4], double rX[4], double rY[4], double rZ[4]){
    
    casteljau(pX, qX, rX);
    casteljau(pY, qY, rY);
    casteljau(pZ, qZ, rZ);
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2){
    double result = 0;
    result = (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2);
    return sqrt(result);
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

        double qX[4],qY[4],qZ[4],rX[4],rY[4],rZ[4];
        Bezier_Divide(X, Y, Z, qX, qY, qZ, rX, rY, rZ);

        double size1 = distance(qX[0], qY[0], qZ[0], qX[3], qY[3], qZ[3]);
        double size2 = distance(rX[0], rY[0], rZ[0], rX[3], rY[3], rZ[3]);
        if (size1 <2 || size2 < 2){
            for (int i=0; i<3; i++){
                Line l1, l2;

                l1.SetPointByIndex(0, qX[i], qY[i], qZ[i]);
                l1.SetPointByIndex(1, qX[i+1], qY[i+1], qZ[i+1]);
                l1.SetColor(color[0], color[1], color[2]);
                
                l2.SetPointByIndex(0, rX[i], rY[i], rZ[i]);
                l2.SetPointByIndex(1, rX[i+1], rY[i+1], rZ[i+1]);
                l2.SetColor(color[0], color[1], color[2]);

                l1.Draw(screen);
                l2.Draw(screen);
            }
        }
        else{
            BezierCurve c1, c2;
            for (int j=0; j<4; j++){
                c1.SetPointByIndex(j, qX[j], qY[j], qZ[j]);
                c2.SetPointByIndex(j, rX[j], rY[j], rZ[j]);
            }
            c1.SetColor(color[0], color[1], color[2]);
            c2.SetColor(color[0], color[1], color[2]);

            c1.Draw(screen);
            c2.Draw(screen);
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
        printf("v3: %lf, %lf, %lf\n", X[3], Y[3], Z[3]);
        printf("v4: %lf, %lf, %lf\n", X[4], Y[4], Z[4]);
        printf("v5: %lf, %lf, %lf\n", X[5], Y[5], Z[5]);
        printf("v6: %lf, %lf, %lf\n", X[6], Y[6], Z[6]);
        printf("v7: %lf, %lf, %lf\n", X[7], Y[7], Z[7]);
        printf("v8: %lf, %lf, %lf\n", X[8], Y[8], Z[8]);
        printf("v9: %lf, %lf, %lf\n", X[9], Y[9], Z[9]);
        printf("v10: %lf, %lf, %lf\n", X[10], Y[10], Z[10]);
        printf("v11: %lf, %lf, %lf\n", X[11], Y[11], Z[11]);
        printf("v12: %lf, %lf, %lf\n", X[12], Y[12], Z[12]);
        printf("v13: %lf, %lf, %lf\n", X[13], Y[13], Z[13]);
        printf("v14: %lf, %lf, %lf\n", X[14], Y[14], Z[14]);
        printf("v15: %lf, %lf, %lf\n", X[15], Y[15], Z[15]);

    }

    void Draw(Screen * screen, int divisions){

        if (divisions ==0){
            // split surface into curves and draw.
            for (int i=0; i<4;i++){

                double qX[4],qY[4],qZ[4],rX[4],rY[4],rZ[4];
                double rowX[4], rowY[4], rowZ[4];
                double colX[4], colY[4], colZ[4];

                for (int j=0; j<4; j++){
                    rowX[j] = X[4*i +j];
                    rowY[j] = Y[4*i +j];
                    rowZ[j] = Z[4*i +j];
                    colX[j] = X[4*j +i];
                    colY[j] = Y[4*j +i];
                    colZ[j] = Z[4*j +i];
                }
                BezierCurve c1, c2, c3, c4;
                Bezier_Divide(rowX, rowY, rowZ, qX, qY, qZ, rX, rY, rZ);
                for (int k=0; k<4; k++){
                    c1.SetPointByIndex(k, qX[k], qY[k], qZ[k]);
                    c2.SetPointByIndex(k, rX[k], rY[k], rZ[k]);
                }
                
                Bezier_Divide(colX, colY, colZ, qX, qY, qZ, rX, rY, rZ);
                for (int k=0; k<4; k++){
                    c3.SetPointByIndex(k, qX[k], qY[k], qZ[k]);
                    c4.SetPointByIndex(k, rX[k], rY[k], rZ[k]);
                }

                c1.SetColor(color[0], color[1], color[2]);
                c2.SetColor(color[0], color[1], color[2]);
                c3.SetColor(color[0], color[1], color[2]);
                c4.SetColor(color[0], color[1], color[2]);

                c1.Draw(screen);
                c2.Draw(screen);
                c3.Draw(screen);
                c4.Draw(screen);

            }
        }
        else{
            // split surface and recursion
            divisions--;
            BezierSurface s1, s2;
            BezierSurface subS1, subS2, subS3, subS4;

            subS1.SetColor(color[0], color[1], color[2]);
            subS2.SetColor(color[0], color[1], color[2]);
            subS3.SetColor(color[0], color[1], color[2]);
            subS4.SetColor(color[0], color[1], color[2]);
            

            double qX[4],qY[4],qZ[4],rX[4],rY[4],rZ[4];
            double qsX[16], qsY[16], qsZ[16], rsX[16], rsY[16], rsZ[16];


            // split into 2 sub surfaces
            for (int i=0; i<4;i++){
                double rowX[4], rowY[4], rowZ[4];
                for (int j=0; j<4; j++){
                    rowX[j] = X[4*i +j];
                    rowY[j] = Y[4*i +j];
                    rowZ[j] = Z[4*i +j];
                }
                Bezier_Divide(rowX, rowY, rowZ, qX, qY, qZ, rX, rY, rZ);

                for (int k=0; k<4; k++){
                    qsX[4*i +k] = qX[k];
                    qsY[4*i +k] = qY[k];
                    qsZ[4*i +k] = qZ[k];
                    rsX[4*i +k] = rX[k];
                    rsY[4*i +k] = rY[k];
                    rsZ[4*i +k] = rZ[k];
                }
            }
            s1.SetPoints(qsX, qsY, qsZ);
            s2.SetPoints(rsX, rsY, rsZ);


            // Columns
            for (int i=0; i<4;i++){
                double colX[4], colY[4], colZ[4];
                for (int j=0; j<4; j++){
                    colX[j] = s1.X[4*j +i];
                    colY[j] = s1.Y[4*j +i];
                    colZ[j] = s1.Z[4*j +i];
                }
                Bezier_Divide(colX, colY, colZ, qX, qY, qZ, rX, rY, rZ);
                for (int k=0; k<4; k++){
                    qsX[4*i +k] = qX[k];
                    qsY[4*i +k] = qY[k];
                    qsZ[4*i +k] = qZ[k];
                    rsX[4*i +k] = rX[k];
                    rsY[4*i +k] = rY[k];
                    rsZ[4*i +k] = rZ[k];
                }
            }
            subS1.SetPoints(qsX, qsY, qsZ);
            subS2.SetPoints(rsX, rsY, rsZ);


            for (int i=0; i<4;i++){
                double colX[4], colY[4], colZ[4];
                for (int j=0; j<4; j++){
                    colX[j] = s2.X[4*j +i];
                    colY[j] = s2.Y[4*j +i];
                    colZ[j] = s2.Z[4*j +i];
                }
                Bezier_Divide(colX, colY, colZ, qX, qY, qZ, rX, rY, rZ);
                for (int k=0; k<4; k++){
                    qsX[4*i +k] = qX[k];
                    qsY[4*i +k] = qY[k];
                    qsZ[4*i +k] = qZ[k];
                    rsX[4*i +k] = rX[k];
                    rsY[4*i +k] = rY[k];
                    rsZ[4*i +k] = rZ[k];
                }
            }
            subS3.SetPoints(qsX, qsY, qsZ);
            subS4.SetPoints(rsX, rsY, rsZ);

            subS1.Draw(screen, divisions);
            subS2.Draw(screen, divisions);
            subS3.Draw(screen, divisions);
            subS4.Draw(screen, divisions);

        }

    }
};

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

        Camera camera = GetCamera(i*10, 1000);
        camera.config();
        Matrix ct = camera.CameraTransform();
        Matrix vt = camera.ViewTransform();
        Matrix intermediate = Matrix::ComposeMatrices(ct, vt);
        Matrix dt = camera.DeviceTransform(screen.width, screen.height);
        Matrix total = Matrix::ComposeMatrices(intermediate, dt);

        BezierSurface * s = getSurfaces();
        
        for (int sidx=0; sidx<2; sidx++){
            // cout<<i<<"\n";
            // s[sidx].Print();
            for (int n = 0; n<16; n++){
                double v[4];
                double transformed[4];
                v[0] = s[sidx].X[n];
                v[1] = s[sidx].Y[n];
                v[2] = s[sidx].Z[n];
                v[3] = 1;
                normalize(v, v, 4);

                total.TransformPoint(v, transformed);
                s[sidx].SetPointByIndex(n, transformed[0] / transformed[3], transformed[1] / transformed[3], transformed[2] / transformed[3]);
            }
            // cout<<"\n";
            // s[sidx].Print();
            // cout<<"\n\n";
            s[sidx].Draw(&screen, 3);
        }

        char name[20];
        // make sure you have a directory named images
        // so that writing images won't cause you errors
        sprintf(name, "./images/frame%02d", i);
        // sprintf(name, "./images/frame%02d", i);
        WriteImage(image, name);

    }
}
