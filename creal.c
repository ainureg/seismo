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
double FMatr(int n,double h);	
double DMatr(int n);
double KMatr(int n);
double MMatr(int n);

int main(){
    double h=1/N;//шаг по прямой
    double t=h/2 ;//шаг по времени
    double Nt=tmax/t ;// шагов по времени
    int k;
    int p,ptemp;
    ptemp=(int)ceil(log(N-1)/log(2));
    p=(int)ptemp;
    k=(int)pow(2,(int)ceil(log(N-1)/log(2)));
    int i,j;
    printf("%d", k);
    double **E;
    int n;
    scanf("%d",&n);     			//FMatr(n,h);     //n=3;
    
    E = (double**)malloc((n) * sizeof(double*));
    if (!E) {
 	printf("error");
 	exit(1);
 }

    for (i=0;i<n;i++){
        E[i]=(double*)malloc((n) * sizeof(double));
    }
    EMatr(n,E);
//    printf("%f",k);
 
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            printf("%lf  ", E[i][j]);
        }
        printf("\n");
    }
    
return 0; 
}

double DMatr(int n){
    double** M = (double**)malloc((n*n) * sizeof(double));
    int i,j;
    
    for (i=0;i<n;i++){
        M[i]=(double*)malloc((n*n) * sizeof(double));
}
    M[n][n]=1/cf(1);
    return **M;
    
    }
double FMatr(int n, double h){
    

    double** M = (double**)malloc((n*n) * sizeof(double));
    int i,j;
    
    for (i=0;i<n;i++){
        M[i]=(double*)malloc((n*n) * sizeof(double));
}
    
        M[1][1]=w*w*h*gf(0.5)/6+1/h;


    
    return **M;
    
    }
double gf(int m){
    double h=1/N;//шаг по прямой
    return 1/cf(m*h)/cf(m*h);
}
double cf(double x){
    return (0.1 +3.6*(x-0.5)*(x-0.5));
}
void EMatr(int n,double** M){
    int i,j;
    for (i=0;i<n;i++){
        M[i][i]=1;
        }
    }
double KMatr(int n){
    

    double** M = (double**)malloc((n*n) * sizeof(double));
    int i,j;
    
    for (i=0;i<n;i++){
        M[i]=(double*)malloc((n*n) * sizeof(double));
}
    
    
    for (i=0;i<n;i++){
        M[i][i]=2;
        M[i-1][i]=-1;
        M[i][i-1]=-1;
}
    M[1][1]=1;
    M[n][n-1]=-1;
    M[n-1][n]=-1;
    M[n][n]=1;

    return **M;
    }
double MMatr(int n){      

    double** M = (double**)malloc((n*n) * sizeof(double));
    int i,j;
    
    for (i=0;i<n;i++){
        M[i]=(double*)malloc((n*n) * sizeof(double));
}
    
    
    for (i=1;i<(n-1);i++){
        M[i][i]=2*(gf(i-0.5)+gf(i+0.5));
        M[i-1][i]=gf(i-0.5);
        M[i][i-1]=gf(i-0.5);
}

M[1][1]=2*(gf(1-0.5)+gf(1+0.5));
M[n][n-1]=gf(n-0.5);
M[n-1][n]=gf(n-0.5);
M[n][n]=2*(gf(n-0.5));
    return **M;
    }    