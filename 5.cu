#include <cuComplex.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#define N 400 //размер расчетной прямой
#define tmax 10 // нужный момент времени

float mu=100;
float const w = 5;
float h=(float) 2/N;//шаг по прямой;
void KMatr(int n, cuFloatComplex** M);
cuFloatComplex prm (float x);
cuFloatComplex sigm(float x);
float gf(float m);
float cf(float x);
void MMatr(int n, cuFloatComplex** M);
void FMatr(int n, cuFloatComplex* M);
int solve(int n, cuFloatComplex *top, cuFloatComplex *mid, cuFloatComplex *bot, cuFloatComplex *b, cuFloatComplex *x);

int main( int argc, char *  argv []  ){
    
    
    int k,i,n,j;
    k=(int)pow(2,(int)ceil(log(N-1)/log(2)));
    n=k+1;
    cuFloatComplex **K, **M, *F, **A;
    K = (cuFloatComplex**)malloc((n) * sizeof(cuFloatComplex*));
    M = (cuFloatComplex**)malloc((n) * sizeof(cuFloatComplex*));
    F = (cuFloatComplex*) malloc((n) * sizeof(cuFloatComplex)) ;
    A = (cuFloatComplex**)malloc((n) * sizeof(cuFloatComplex*));
    
    for (i=0;i<n;i++){
        K[i]=(cuFloatComplex*)malloc((n) * sizeof(cuFloatComplex));
        M[i]=(cuFloatComplex*)malloc((n) * sizeof(cuFloatComplex));
        A[i]=(cuFloatComplex*)malloc((n) * sizeof(cuFloatComplex));
    }
    
     KMatr(N, K);
     MMatr(N, M);
     FMatr(N, F);
     for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            K[i][j]=cuCdivf(K[i][j],make_cuFloatComplex(h,0)) ;
            M[i][j]=cuCmulf(M[i][j],make_cuFloatComplex(h/6,0));
            } 
        }
     for (i=0;i<N;i++){
       for (j=0;j<N;j++){
           A[i][j]=cuCaddf( cuCmulf( make_cuFloatComplex( -w*w,0), M[i][j]), K[i][j]);
        }        
     }
     cuFloatComplex *a,*b,*c,*u;
    a = (cuFloatComplex*)malloc((n)  * sizeof(cuFloatComplex));
    b = (cuFloatComplex*)malloc((n)  * sizeof(cuFloatComplex));
    c = (cuFloatComplex*)malloc((n)  * sizeof(cuFloatComplex));
    u = (cuFloatComplex*)malloc((n)  * sizeof(cuFloatComplex));
    
     for (i=1;i<(N-1);i++){
        a[i+1]=A[i+1][i];
        b[i]=A[i][i];
        c[i-1]=A[i][i-1];
    }
    
    a[0]=make_cuFloatComplex(0,0);
    a[1]=A[1][0];
    b[0]=A[0][0];
    b[N-1]=A[N-1][N-1];
    c[N-2]=A[N-1][N-2];
    c[N-1]=make_cuFloatComplex(0,0);
     
    for (j=N;j<(n);j++) {
            a[j]=make_cuFloatComplex(0,0);
            b[j]=make_cuFloatComplex(1,0);
            c[j]=make_cuFloatComplex(0,0);
            F[j]=make_cuFloatComplex(0,0);
      }
      
      
    int l;
    l=solve (n,a,b,c,F,u);
}     

int solve(int n, cuFloatComplex *top, cuFloatComplex *mid, cuFloatComplex *bot, cuFloatComplex *b, cuFloatComplex *x){
    return 0;
}

void MMatr(int n, cuFloatComplex** M){
   int i;
   cuFloatComplex t1,t2;
   for (i=0;i<(n-1);i++) {
        t1=cuCmulf(make_cuFloatComplex( gf(i+1.5),0), prm(i+1.5));
        t2=cuCmulf( make_cuFloatComplex(gf(i+0.5),0) , prm(i+0.5));
        
       
        M[i][i]=cuCmulf(make_cuFloatComplex(2,0), cuCaddf( t1 ,t2 ));
        M[i+1][i]=t1;
        M[i][i+1]=t1;
        }
    M[n-1][n-1]=cuCmulf(make_cuFloatComplex( 2*gf(n-0.5),0), prm(n-0.5)) ;
    }
cuFloatComplex prm ( float x){
    return cuCdivf(cuCaddf( make_cuFloatComplex(w,0), cuCmulf(make_cuFloatComplex(0,-1), sigm(x*h) )) , make_cuFloatComplex(w,0));
    }        
cuFloatComplex sigm(float x){
if (x <= 1)
    return  make_cuFloatComplex (0,0);
else
    return make_cuFloatComplex( mu*(x-1)*(x-1),0);
    }    
float gf(float m){
    return ( ( 1/cf(m*h) ) / cf(m*h) );
    }       
float cf(float x){
    if (x<=1)
    return (0.1 +3.6*(x-0.5)*(x-0.5));
    else
    return 1;
     }    
void KMatr(int n, cuFloatComplex ** M){
    int i;
    cuFloatComplex t1,t2;

    for (i=0;i<(n-1);i++){
        t1=cuCdivf(make_cuFloatComplex(1,0), prm(i+0.5));
        t2=cuCdivf(make_cuFloatComplex(1,0), prm(i+1.5));
        
        
        M[i][i]= cuCaddf(t1,t2);
        M[i+1][i]=cuCdivf(make_cuFloatComplex(-1,0), prm(i+1.5));
        M[i][i+1]=cuCdivf(make_cuFloatComplex(-1,0), prm(i+1.5));
}
    M[n-1][n-1]=cuCdivf(make_cuFloatComplex(1,0), prm(n-0.5));
    }     
void FMatr(int n, cuFloatComplex * M){
    M[0]=cuCaddf( cuCmulf(make_cuFloatComplex(w*w*h*gf(0.5)/6,0) , prm(0.5) ) , cuCmulf(make_cuFloatComplex(1/h,0),prm(0.5)));
    }    