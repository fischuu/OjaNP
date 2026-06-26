#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

/* NOTE (fixes_j / CRAN compliance):
 * - Replaced Calloc/Free macros with R_Calloc/R_Free (R-allocated memory)
 * - Added explicit C23-compliant function prototypes (removed outdated extern decls)
 *   Previously used `extern double det(); void dp(), nextp();` style,
 *   which is implicit int in old C - flagged by modern compilers.
 */
/* Function Prototypes inline with C23 */
double det(double a[], int n);
double sgn(double x);
double sp(double x[], int k, int p[], double obsvec[], int a[], double *d);
double sq(double x[], int k, int q[], double obsvec[], double *d);
void dp(double x[], int n, int k, int p[], int a[], double cf[]);
void eq(double x[], int n, int k, int q[], double cf[]);
void nextp(int n, int k, int p[], int *last);

/******************************************************************/
/******************************************************************/
void ojacrnk(double *x, int *nrow, int *ncol, double *rn)     
{ 
  int i,j,m,n,k,nk,cnt=0,lc,h,*p,*a;
  double sp2,da,*dp2,*obsvec;
  
  n=*nrow;
  k=*ncol;
  nk=(*nrow)*(*ncol);
  //kk=(*ncol)*(*ncol);

  p = (int*) R_Calloc(k,int); 
  if (p == NULL) error("memory allocation failed");
  a = (int*) R_Calloc(k,int);
  if (a == NULL) error("memory allocation failed");
  dp2 = (double*) R_Calloc(k,double);
  if (dp2 == NULL) error("memory allocation failed");
  obsvec = (double*) R_Calloc(k,double);
  if (obsvec == NULL) error("memory allocation failed");

  /* Initialization */
  for(i=0;i<nk;i++)rn[i]=0;
  for(h=0; h<k; h++)a[h]=1;

  for(i=0; i<n; i++){
    cnt=0;
    lc=1;
    for(m=0; m<k; m++)
    { 
      p[m]=m+1;
      obsvec[m]=x[i*k+m];
    }
    while(lc!=0)
    {
      cnt = cnt+1;
      sp2=sp(x,k,p,obsvec,a,&da);
      dp(x,n,k,p,a,dp2);     
      for(j=0; j<k; j++)rn[i*k+j]=rn[i*k+j]+sp2*dp2[j];
      nextp(n,k,p,&lc);
    }
  }

  R_Free(p);
  R_Free(a);
  R_Free(dp2);
  R_Free(obsvec);

  for(i=0;i<nk;i++)rn[i]/=cnt;
}   

/******************************************************************/
/******************************************************************/
void ojasrnk(double *x, int *nrow, int *ncol, double *rn)     
{  
  int i,j,m,n,k,nk,*p,cnt=0,lc,*a,h,na,hh;   
  double sp2,*dp2,*obsvec,da;
 
  n=*nrow;
  k=*ncol;
  nk=(*nrow)*(*ncol);
  //kk=(*ncol)*(*ncol);

  p = (int*) R_Calloc(k,int);
  if (p == NULL) error("memory allocation failed");
  a = (int*) R_Calloc(k,int);
  if (a == NULL) error("memory allocation failed");
  dp2 = (double*) R_Calloc(k,double);
  if (dp2 == NULL) error("memory allocation failed");
  obsvec = (double*) R_Calloc(k,double);
  if (obsvec == NULL) error("memory allocation failed");

  /* Initialization */
  for(i=0;i<nk;i++)rn[i]=0;

  for(i=0; i<n; i++){
    cnt=0;
    lc=1;
    for(m=0; m<k; m++)
    {
      p[m]=m+1;
      obsvec[m]=x[i*k+m];
    }
    na = (int) pow(2.,k);

    while(lc!=0)
    {
      for(h=0; h<k; h++)a[h] = -1;
      for(hh=0; hh<na; hh++)
      {
        cnt = cnt+1;
        sp2=sp(x,k,p,obsvec,a,&da);
        dp(x,n,k,p,a,dp2);   
        for(j=0; j<k; j++)rn[i*k+j]=rn[i*k+j]+sp2*dp2[j];
        if(hh<(na-1)){
          h=k-1;
          while(a[h]!=-1){
            a[h]=-1;
            h=h-1;
          }
          a[h]=1;
        }
      }
      nextp(n,k,p,&lc);
    }
  }

  R_Free(p);
  R_Free(a);
  R_Free(dp2);
  R_Free(obsvec);

  for(i=0;i<nk;i++)rn[i]/=cnt;
}   



/******************************************************************/
/******************************************************************/
void ojasrnk2(double *x, int *nrow, int *ncol, double *rn)     
{  
  int i,j,m,n,k,nk,*p,cnt=0,lc,*a,h;   
  double sp2,*dp2,*obsvec,*x2,da;
 
  n=*nrow;
  k=*ncol;
  nk=(*nrow)*(*ncol);
  //kk=(*ncol)*(*ncol);

  x2      = (double*)R_Calloc((2*n)*k,double);
  if (x2 == NULL) error("memory allocation failed");
  p       = (int*)R_Calloc(k,int);
  if (p == NULL) error("memory allocation failed");
  dp2     = (double*)R_Calloc(k,double);
  if (dp2 == NULL) error("memory allocation failed");
  obsvec  = (double*)R_Calloc(k,double);
  if (obsvec == NULL) error("memory allocation failed");
  a       = (int*)R_Calloc(k,int);
  if (a == NULL) error("memory allocation failed");

  for(i=0; i<(2*n); i++){
    for(j=0; j<k; j++){
      if(i<n)
        x2[i*k+j]=x[i*k+j];
      else
        x2[i*k+j]= -x[(i-n)*k+j];
    }
  }

  /* Initialization */
  for(i=0;i<nk;i++)rn[i]=0;
  for(h=0; h<k; h++)a[h]=1;

  for(i=0; i<n; i++){
    cnt=0;
    lc=1;
    for(m=0; m<k; m++)
    { 
      p[m]=m+1;
      obsvec[m]=x[i*k+m];
    }
    while(lc!=0)
    {
      cnt = cnt+1;
      sp2=sp(x2,k,p,obsvec,a,&da);
      dp(x2,2*n,k,p,a,dp2);     
      for(j=0; j<k; j++)rn[i*k+j]=rn[i*k+j]+sp2*dp2[j];
      nextp(2*n,k,p,&lc);
    }
  }

  R_Free(x2);
  R_Free(p);
  R_Free(dp2);
  R_Free(obsvec);
  R_Free(a);

  for(i=0;i<nk;i++)rn[i]/=cnt;
}   



/******************************************************************/
/******************************************************************/
void ojasn(double *a, int *nrow, int *ncol, double *sn)     
{  
  int i,j,m,n,k,nk,*q,cnt=0,lc;   
  double sq2,*eq2,*obsvec,da;
   
  n=*nrow;
  k=*ncol;
  nk=(*nrow)*(*ncol);
  //kk=(*ncol)*(*ncol);

  q       = (int*)R_Calloc(k,int);
  if (q == NULL) error("memory allocation failed");
  eq2     = (double*)R_Calloc(k,double);
  if (eq2 == NULL) error("memory allocation failed");
  obsvec  = (double*)R_Calloc(k,double);
  if (obsvec == NULL) error("memory allocation failed");

  /* Initialization */
  for(i=0;i<nk;i++)
    sn[i]=0.0;


  for(i=0; i<n; i++){
    cnt=0;
    lc=1;
    for(m=0; m<k; m++)
    { 
      q[m]=m+1;
      obsvec[m]=a[i*k+m];
    }   
    while(lc!=0)
    {
      cnt += 1;
      sq2 = sq(a,k,q,obsvec,&da);
      eq(a,n,k,q,eq2);     
      for(j=0; j<k; j++)sn[i*k+j]=sn[i*k+j]+sq2*eq2[j];
      nextp(n,k-1,q,&lc);
    }
  }

  R_Free(q);
  R_Free(eq2);
  R_Free(obsvec);

  for(i=0;i<nk;i++)sn[i]/=cnt; 
}   


/******************************************************************/
/******************************************************************/
double sgn(double x)
{
  float ans;

  if(x<0)ans=-1.;
  else if(x>0)ans=1.;
  else ans=0.;
  return(ans);
}

/******************************************************************/
/******************************************************************/

double sp(double x[], int k, int p[], double obsvec[], int a[], double *d)
{  
  int i, j;
  double *x2,maxabsj=0.0,absj,g;

  x2 = (double*)R_Calloc(k*k,double);
  if (x2 == NULL) error("memory allocation failed");

  for(i=0; i<k; i++){
    for(j=0; j<k; j++){
      x2[i*k+j]=a[i]* x[(p[i]-1)*k+j]-obsvec[j];
      if((i==0)&(j==0))maxabsj=fabs(x2[0]);
      absj=fabs(x2[i*k+j]);
      if(absj>maxabsj)maxabsj=absj;
    }
  }
  
  g=pow(maxabsj,k);
  *d=det(x2,k);
  if((fabs(*d)/g)<1.0e-10) *d=0.0;

  R_Free(x2);

  return(sgn(*d));	   /* sign of the determinant */
}



/******************************************************************/
/******************************************************************/
double sq(double x[], int k, int q[], double obsvec[], double *d)
{  
  int i, j;
  double *x2,maxabsj=0.0,absj,g;

  x2 = (double*)R_Calloc(k*k,double);
  if (x2 == NULL) error("memory allocation failed");

  for(j=0; j<k; j++){
    x2[j]=obsvec[j];
    if(j==0)maxabsj=fabs(x2[0]);
    absj=fabs(x2[j]);
    if(absj>maxabsj)maxabsj=absj;
  }
  
  for(i=1; i<k; i++){
    for(j=0; j<k; j++){
      x2[i*k+j]=x[(q[i-1]-1)*k+j];
      absj=fabs(x2[i*k+j]);
      if(absj>maxabsj)maxabsj=absj;
    }
  }
 
  g=pow(maxabsj,k);
  *d=det(x2,k);
  if((fabs(*d)/g)<1.0e-10) *d=0.0;

  R_Free(x2); 

  return(sgn(*d));	   /* sign of the determinant */
}



/******************************************************************/
/******************************************************************/
void dp(double x[], int n, int k, int p[], int a[], double cf[])
{
  //int i, j, h, u, v;     2009.05.08
  int i, j, u, v;
  double *x2;

  x2 = (double*)R_Calloc(k*k,double);
  if (x2 == NULL) error("memory allocation failed");

  for(u=0;u<k;u++) /* u:th cofactor */
  {
    for(i=0;i<k;i++) /* i:th row */
    {
      v=1;
      x2[i*k]=1.0; /* i:th row and 0:th column */
      for(j=0;j<k;j++) 
      {
        if(j!=u) /* don't take the u:th component */
        {
          x2[i*k+v]=a[i]*x[(p[i]-1)*k+j]; 
          v += 1;
        }
      }
    }
    cf[u]=pow(-1.,u+1)*det(x2,k);
  }
 
  R_Free(x2);
  
}


/******************************************************************/
/******************************************************************/
void eq(double x[], int n, int k, int q[], double cf[])
{
 // int i, j, h, u, v;   2009.05.08
  int i, j, u, v;
  double *x2;

  x2 = (double*)R_Calloc(k*k,double);
  if (x2 == NULL) error("memory allocation failed");

  for(u=0;u<k;u++) /* u:th cofactor */
  {
    for(i=0;i<(k-1);i++) /* i:th row */
    {
      v=0;
      for(j=0;j<k;j++) 
      {
        if(j!=u) /* don't take the u:th component */
        {
          x2[i*(k-1)+v] = x[(q[i]-1)*k+j]; 
          v += 1;
        }
      }
    }
    cf[u]=pow(-1.,u)*det(x2,k-1); 
  }

  R_Free(x2);

}



/******************************************************************/
/******************************************************************/
void nextp(int n, int k, int p[], int *last)
{  
  int i, lc;
       
  lc = k;
  if(p[k-1] < n){
    p[k-1] = p[k-1] + 1;
  }
  else{
    while(1){
      lc = lc - 1;
      if(lc==0)
        break;
      if(p[lc-1] < (n-k+lc))
        break;
    }
    if(lc > 0){
      p[lc-1] = p[lc-1] + 1;
      for(i=(lc+1);i<=k;i++)
	p[i-1] = p[i-2]+1;
    }
  }
  *last=lc;
}
