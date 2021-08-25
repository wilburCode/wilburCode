#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <Hyper.h>
using namespace std;
namespace iret {

Hyper::Hyper(void) {
}

Hyper::Hyper(long n) {
   long i;
   nobj=n;

   log_num=new double[n+1];
   for(i=1;i<n+1;i++)log_num[i]=log10((double)i);

   log_fac=new double[n+1];
   *log_fac=*(log_fac+1)=0;
   for(long i=2;i<n+1;i++)*(log_fac+i)=*(log_fac+i-1)+log_num[i];
}

Hyper::Hyper(long n,double epx) {
   long i;
   nobj=n;

   log_num=new double[n+1];
   for(i=1;i<n+1;i++)log_num[i]=log10((double)i);

   log_fac=new double[n+1];
   *log_fac=*(log_fac+1)=0;
   for(i=2;i<n+1;i++)*(log_fac+i)=*(log_fac+i-1)+log_num[i];

   if(epx>0){eps=-log(epx);epz=epx;}
   else {cout << "Error, eps>0 is necessary!" << endl;exit(0);}
}


Hyper::~Hyper() {
  delete [] log_fac;
  if(log_num!=NULL)delete [] log_num;
}

double Hyper::nlog_pval(long n_st,long n_s, long n_t, long N){
   double x,y;
   long i,j,m;
  
   if(N>nobj){cout << "Error in size nobj = " << nobj << " and N = " << N << endl;exit(0);}
   m=(n_s<n_t)?n_s:n_t;
   x=log_prob(n_st,n_s,n_t,N);
   for(i=n_st+1;i<=m;i++){
      y=log_prob(i,n_s,n_t,N);
      x=addl(x,y);
   }
   return(-x);
}

double Hyper::nlog_pval_appx(long n_st,long n_s, long n_t, long N){
   double x,y,*mt;
   long i,j,m;
  
   if(N>nobj){cout << "Error in size nobj = " << nobj << " and N = " << N << endl;exit(0);}
   m=(n_s<n_t)?n_s:n_t;
   x=log_prob(n_st,n_s,n_t,N);
   if(n_st>=((double)n_s)*((double)n_t)/((double)N)){
      mt=log_num+m-n_st;
      for(i=n_st+1;i<=m;i++){
         y=log_prob(i,n_s,n_t,N);
         if(y+eps+*(mt--)>x){
            x=addl(x,y);
         }
         else break;
      }
      return(-x);
   }
   else return(this->nlog_pval(n_st,n_s,n_t,N));
}

double Hyper::nlog_pval_appx2(long n_st,long n_s, long n_t, long N){
   double x,y,*mt,nx,ny,nr,nn,xa,xb,xc;
   long i,j,m;

   if(N>nobj){cout << "Error in size nobj = " << nobj << " and N = " << N << endl;exit(0);}
   m=(n_s<n_t)?n_s:n_t;
   x=log_prob(n_st,n_s,n_t,N);
   if(n_st>=((double)n_s)*((double)n_t)/((double)N)){
      nx=n_s-n_st+1;ny=n_t-n_st+1;nr=n_st;nn=N-n_s-n_t+n_st;
      xa=xb=1.0;
      for(i=n_st+1;i<=m;i++){
         nx--;ny--,nr++;nn++;
         xc=nx*ny/(nr*nn);
         xb*=xc;
         xa+=xb;
         if((xb<epz)&&(xc<=0.5))break;
      }
      return(-x-log10(xa));
   }
   else return(this->nlog_pval(n_st,n_s,n_t,N));
}

double Hyper::log_prob(long n_st,long n_s, long n_t, long N){

return log_fac[n_s] + log_fac[n_t] + log_fac[N - n_s] + log_fac[N - n_t]
-log_fac[N] -log_fac[n_st]- log_fac[n_s - n_st] - log_fac[n_t -n_st]
- log_fac[N -n_s -n_t + n_st];
}

double Hyper::addl(double x,double y)
{
double xt,u;
u=x-y;
if(u<=0){
   xt=pow(10.0,u);
   return(y+log10(1.0+xt));
         }
else {
   xt=pow(10.0,-u);
   return(x+log10(1.0+xt));
      }
}

double Hyper::log_binomCoeff(long m,long n){
   return(log_fac[n]-log_fac[m]-log_fac[n-m]);
}

double Hyper::logOdds(long n_st,long n_s,long n_t,long N){
   long k;
   double xx,yy,zz,uu,vv;

   if(n_t>n_s){
     k=n_t;
     n_t=n_s;
     n_s=k;
   }
   xx=(double)n_st;
   yy=(double)n_s;
   zz=(double)n_t;
   uu=(double)N;

   if(xx>yy*zz/uu){
     vv=-yy*log_num[n_s+2]+xx*log_num[n_st+1]+(yy-xx)*log_num[n_s-n_st+1];
     vv-=log_binomCoeff(n_t-n_st,N-n_s);
     vv+=log_binomCoeff(n_t,N);
   }
   else {
     vv=yy*log_num[n_s+2]-xx*log_num[n_st+1]-(yy-xx)*log_num[n_s-n_st+1];
     vv+=log_binomCoeff(n_t-n_st,N-n_s);
     vv-=log_binomCoeff(n_t,N);
   }
   return(vv);
}

double Hyper::HlogOdds(long n_st,long n_s,long n_t,long N){
   long k;
   double xx,yy,zz,uu,vv;

   if(n_t>n_s){
     k=n_t;
     n_t=n_s;
     n_s=k;
   }
   xx=(double)n_st;
   yy=(double)n_s;
   zz=(double)n_t;
   uu=(double)N;

   if(xx>yy*zz/uu){
     k=(long)ceil(yy*zz/xx);
     vv=log_binomCoeff(n_t-n_st,k-n_s);
     vv-=log_binomCoeff(n_t,k);
     vv-=log_binomCoeff(n_t-n_st,N-n_s);
     vv+=log_binomCoeff(n_t,N);
   }
   else {
     k=(long)ceil((uu-yy)*zz/(zz-xx));
     vv=-log_binomCoeff(n_st,k-N+n_s);
     vv+=log_binomCoeff(n_t,k);
     vv+=log_binomCoeff(n_st,n_s);
     vv-=log_binomCoeff(n_t,N);
   }
   return(vv);
}

//SHyper

SHyper::SHyper(void) {
}

SHyper::SHyper(double epx) {
   if(epx>0)eps=-log(epx);
   else {cout << "Error, eps>0 is necessary!" << endl;exit(0);}
}


SHyper::~SHyper() {
}

double SHyper::nlog_pval(long n_st,long n_s, long n_t, long N){
   double x,y;
   long i,j,m;
  
   m=1+((n_s<n_t)?n_s:n_t);
   set_log_prob(n_st,n_s,n_t,N);
   x=xi;
   for(i=n_st+1;i<m;i++){
      log_prob(i);
      x=addl(x,xi);
   }
   return(-x);
}

double SHyper::nlog_pval_appx(long n_st,long n_s, long n_t, long N){
   double x,y;
   long i,j,m;
  
   m=1+((n_s<n_t)?n_s:n_t);
   set_log_prob(n_st,n_s,n_t,N);
   x=xi;
   if(n_st>=((double)n_s)*((double)n_t)/((double)N)){
      for(i=n_st+1;i<m;i++){
         log_prob(i);
         if(xi+eps+log10((double)(m-i))>x){
            x=addl(x,xi);
         }
         else break;
      }
      return(-x);
   }
   else return(this->nlog_pval(n_st,n_s,n_t,N));
}

double SHyper::log_prob(long k){
   if(k==i1)return(xi);
   else if(k==i1+1){
      i2--;
      i3--;
      i1++;
      xi+=log10((double)i2)+log10((double)i3)-\
          log10((double)i1)-log10((double)i4);
      i4++;
      return(xi);
   }
   else return(0);
}

void SHyper::set_log_prob(long n_st,long n_s, long n_t, long N){
   i1=n_st;
   i2=n_s+1-i1;
   i3=n_t+1-i1;
   i4=N-n_s-n_t+1+i1;

   xi=log_fac(i2,n_s)+log_fac(i3,n_t)+\
      log_fac(i4,N-n_s)-log_fac(2,i1)-log_fac(N-n_t+1,N);
}

double SHyper::log_fac(long l,long u){
   long i;
   double sum=0;
   for(i=l;i<=u;i++)sum+=log10((double)i);
   return(sum);
}

double SHyper::addl(double x,double y)
{
double xt,u;
u=x-y;
if(u<=0){
   xt=pow(10.0,u);
   return(y+log10(1.0+xt));
         }
else {
   xt=pow(10.0,-u);
   return(x+log10(1.0+xt));
      }
}

}
