#include <cuComplex.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#define h_x 0.01
#define h_y 0.005
#define N 100 //размер расчетной плоскости по X
#define M 400 //размер расчетной плоскости по Y
#define tmax 5 // нужный момент времени
double c(double x, double y);
double v(double t, double x);
void save(const char *prefix, int step, double *a);
int ind (int i, int j) {
return j*(N+1) + i;
}
int main(int argc, char **argv){//int main( int argc, char * argv [] ){
int i;
double t=h_x/2;
double *a = new double[(N+1) * (M+1)];
// U = (double**)malloc((N) * sizeof(double*));
/* // for (i=0;i<N;i++){
U[i]=(double*)malloc((M) * sizeof(double));
}
*/
std::cout << "Bang!" << std::endl;
for (int step = 0; step < 1000; step++){
for (int i = 0; i < N+1; i++)
for (int j = 0; j < M+1; j++){
double x = h_x*i;
double y = h_y*j;
a[ind(i,j)] = x + y + step / 1000.;
}
if (step == (step / 10) * 10 )
save("output", step, a);
}
}
//---------------------------------------------------------------------------//
template<typename T>
static void put(std::fstream &f, const T value) {
union {
char buf[sizeof(T)];
T val;
} helper;
helper.val = value;
std::reverse(helper.buf, helper.buf + sizeof(T));
f.write(helper.buf, sizeof(T));
}
void save(const char *prefix, int step, double *a){
char buffer[50];
sprintf(buffer, "%s.%05d.vtk", prefix, step);
std::fstream f(buffer, std::ios::out);
if (!f){
std::cerr << "Unable to open file " << buffer << std::endl;
return;
}
f << "# vtk DataFile Version 3.0" << std::endl;
f << "U data" << std::endl;
f << "BINARY" << std::endl;
f << "DATASET STRUCTURED_POINTS" << std::endl;
f << "DIMENSIONS " << N+1 << " " << M+1 << " 1" << std::endl;
f << "SPACING " << h_x << " " << h_y << " 1" << std::endl;
f << "ORIGIN 0 0 0" << std::endl;
f << "POINT_DATA " << (N+1) * (M+1) << std::endl;
f << "SCALARS u double" << std::endl;
f << "LOOKUP_TABLE default" << std::endl;
for (int j = 0 ; j < M+1; j++){
for (int i = 0; i < N+1; i++)
put(f, a[ind(i,j)]);
}
}
double v(double t, double x) {
if (5*t<2*M_PI){
return sin(5*t)*exp(-30*(x-0.5)*(x-0.5));
}
else return 0.0;
}
double c(double x, double y) {
if (1.0<y<=1.2){
return 0.8;
}
else if( (0.5<y<=0.8)&(0.2<x<=0.5) ){
return 1.0;
}
else return 0.5;
} 
