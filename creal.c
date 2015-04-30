#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#define N 200//размер расчетной прямой
#define w 5 //частота
#define tmax 10 // нужный момент времени


void EMatr(int n, double** );
double gf(double m);
double cf(double x);
void FMatr(int n,double h, double* );	
void DMatr(int n, double** );
void KMatr(int n, double** );
void MMatr(int n, double** );

int main(){
    double h;
    double N1=200;
    h=1/N1;//шаг по прямой
    //printf("%lf\n",h);
    double t=h/2 ;//шаг по времени
    double Nt=tmax/t ;// шагов по времени
    int k;
    k=(int)pow(2,(int)ceil(log(N-1)/log(2)));
    int i,j;
    double **K, **D, **M, *F;
    int n,n1;
    complex *a,*b,*c,*u;
    //scanf("%d",&n);
    complex** A;
    n1=200;
    n=(k+1);
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
    FMatr(n1,h,F);	
    DMatr(n1, D);
    KMatr(n1, K);
    MMatr(n1, M);

    
    for (i=0;i<n1;i++){
        for (j=0;j<n1;j++){
            K[i][j]=K[i][j]/h;
            M[i][j]=M[i][j]*h/6;
            //D[i][j]=D[i][j];
        }        
    }
    for (i=0;i<n1;i++){
       for (j=0;j<n1;j++){
           A[i][j]=-w*w*M[i][j]+K[i][j] +w*D[i][j]*I;
       }        
    }
//     for (i=0;i<4;i++){
//         for (j=0;j<4;j++){
//             printf("A[%d][%d]=%lf    ",i,j,creal(A[i][j]));
//             }        
//         printf("\n");
//    }
    
    
    
    a = (complex*)malloc((n)  * sizeof(complex));
    b = (complex*)malloc((n)  * sizeof(complex));
    c = (complex*)malloc((n)  * sizeof(complex));
    u = (complex*)malloc((n)  * sizeof(complex));
    
   
    for (i=1;i<(n1-1);i++){
        a[i+1]=A[i+1][i];
        b[i]=A[i][i];
        c[i-1]=A[i][i-1];
    }
    a[0]=0;
    a[1]=A[1][0];
    a[n1-1]= A[n1-1][n1-2];
    b[0]=A[0][0];
    b[n1-1]=A[n1-1][n1-1];
    c[n1-2]=A[n1-1][n1-2];
    c[n1-1]=0;
     
    for (j=n1;j<(n);j++) {
            a[j]=0;
            b[j]=1;
            c[j]=0;
            F[j]=0;
      }
//     for (i=195;i<201;i++){
//          printf("F[%d]=%lf   b[%d]=%lf  c[%d]=%lf   a[%d] = %lf \n",i,creal(F[i]),i,creal(b[i]),i,creal(c[i]), i,creal(a[i]));
//         }
      
    int s=1;
    complex Ctemp,Atemp;
    int q;
    q=0;
	
    while (s<=(k)/2){
	//the first and the last eq-s
	Ctemp=-c[0]/b[s];
	b[0]= b[0]+Ctemp*a[s];
	F[0]=F[0]+Ctemp*F[s];
	c[0]=Ctemp*c[s];
	
	Atemp=-a[k]/b[k-s];
	b[k]=b[k]+Atemp*a[k-s];
	F[k]=F[k]+Atemp*F[k-s];
	a[k]=Atemp*a[k-s];
	
	j=2*s;

	while (j<=(k-1)){
	    //others eq-s
	    Atemp=-a[j]/b[j-s];
	    Ctemp=-c[j]/b[j+s];
	    //#print(c(Atemp, Ctemp))
	    
	    b[j]= b[j]+Atemp*c[j-s]+Ctemp*a[j+s];
	    F[j]= F[j]+Atemp*F[j-s]+Ctemp*F[j+s];
	    a[j]= a[j-s]*Atemp;
	    c[j]= c[j+s]*Ctemp;
	    j=j+2*s;
	    
	    }
	s=s*2;
            }
            //now we can calculate 
    u[0]=(F[0]*b[k]-c[0]*F[k])/(b[0]*b[k]-a[k]*c[0]);  
    u[k]=(F[k]*b[0]-a[k]*F[0])/(b[0]*b[k]-a[k]*c[0]);
//     printf ("UUU=%lf+%lf*I\n",creal( u[0]),cimag(u[0]) ) ;
//     printf("%d\n",q);
    s=(k)/2;
    while (s>=1){
        j=s;
        while (j<(k)){
            q++;
            u[j]=(-F[j]-a[j]*u[j-s]-c[j]*u[j+s])/b[j];
            j=j+2*s;
        }
    s=s/2;
    }
    for (i=195;i<205;i++){
        printf("u[%d] = %lf+%lf \n",i,creal(u[i]),cimag(u[i]));
    }
    printf("%d\n",q);
         
return 0 ;
}

void DMatr(int n, double** M){
    M[n-1][n-1]=1/cf(1);
    }
void FMatr(int n, double h, double* M){
    M[0]=w*w*h*gf(0.5)/6+1/h;
    }
double gf(double m){
    //double h=1/N;//шаг по прямой
    return ( ( 1/cf(m/N) ) / cf(m/N) );
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