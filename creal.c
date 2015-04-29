#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#define N 200//размер расчетной прямой
#define w 5 //частота
#define tmax 10 // нужный момент времени

void EMatr(int n, double** );
double gf(int m);
double cf(double x);
void FMatr(int n,double h, double* );	
void DMatr(int n, double** );
void KMatr(int n, double** );
void MMatr(int n, double** );

int main(){
    double h=1/N;//шаг по прямой
    double t=h/2 ;//шаг по времени
    double Nt=tmax/t ;// шагов по времени
    int k;
    k=(int)pow(2,(int)ceil(log(N-1)/log(2)));
    int i,j;
    printf("%d", k);
    double **K, **D, **M, *F;
    int n;
    //scanf("%d",&n);
    complex** A;
    n=k;
    K = (double**)malloc((n) * sizeof(double*));
    M = (double**)malloc((n) * sizeof(double*));
    F = (double*)malloc((n)  * sizeof(double));
    D = (double**)malloc((n) * sizeof(double*));
    A = (complex**)malloc((n)* sizeof(complex*));
    if ((!A)&(!D)&(!M)&(!F)&(!D)) {
 	printf("error");
 	exit(1);
 }
    for (i=0;i<n;i++){
        K[i]=(double*)malloc((n) * sizeof(double));
        M[i]=(double*)malloc((n) * sizeof(double));
        D[i]=(double*)malloc((n) * sizeof(double));
        A[i]=(complex*)malloc((n) * sizeof(complex));
    }
    //EMatr(n,E);
    FMatr(n,h,F );	
    DMatr(n, D);
    KMatr(n, K);
    MMatr(n, M);

    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            A[i][j]=w*w*M[i][j]+K[i][j] +w*D[i][j]*I;
        }
        //printf("\n");
    }
    
return 0; 
}

void DMatr(int n, double** M){
    M[n-1][n-1]=1/cf(1);
    }
void FMatr(int n, double h, double* M){
    M[0]=w*w*h*gf(0.5)/6+1/h;
    }
double gf(int m){
    double h=1/N;//шаг по прямой
    return 1/cf(m*h)/cf(m*h);
    }
double cf(double x){
    return (0.1 +3.6*(x-0.5)*(x-0.5));
    }
void EMatr(int n,double** M){
    int i;
    for (i=0;i<n;i++){
        M[i][i]=1;
        }
    }
void KMatr(int n, double** M){
    int i;
    for (i=0;i<(n-1);i++){
        M[i][i]=2;
        M[i+1][i]=-1;
        M[i][i+1]=-1;
}
    M[n-1][n-1]=1;
    }
void MMatr(int n, double **M ){      
    int i;
    for (i=0;i<(n-1);i++){
        M[i][i]=2*(gf(i+0.5)+gf(i+1.5));
        M[i+1][i]=gf(i+1.5);
        M[i][i+1]=gf(i+1.5);
        }

    M[n-1][n-1]=2*(gf(n+0.5));
    }    