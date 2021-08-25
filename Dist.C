#include <iostream>
#include <fstream>
#include <cstdlib>
#include <fcntl.h>
#include <sys/mman.h>
#include <cmath>
#include <cstring>
#include <cassert>
#include "Dist.h"

using namespace std;
namespace iret {


double lphi(double x){
   long i;
   double y,z,u,s,xi;

   if(fabs(x)<10.0){
      y=sqrt(2.0);
      if(x>=0)return(log((1.0+erf(x/y))/2.0));
      else return(log(erfc(-x/y)/2.0));
   }
   else if(x>0){
      y=sqrt(2.0*3.14159)*x;
      y=exp(-x*x/2.0)/y;
      z=1.0;
      u=-1.0/(x*x);
      s=1.0;
      for(i=1;i<6;i++){
        xi=(double)i;
        z=z*xi*u;
        s+=z;
      }
      return(log(1.0-y*s));
   }
   else {
      y=-x*x/2.0-log(2.0*3.14159)/2.0-log(fabs(x));
      z=1.0;
      u=-1.0/(x*x);
      s=1.0;
      for(i=1;i<6;i++){
        xi=(double)i;
        z=z*xi*u;
        s+=z;
      }
      return(y+log(s));
   }
}


LBin::LBin(long n){
   log_fac=new double[n+1];
   *log_fac=*(log_fac+1)=0;
   for(long i=2;i<n+1;i++)*(log_fac+i)=*(log_fac+i-1)+log((double)i);
}

LBin::~LBin() {
  if(log_fac) delete [] log_fac;
}

double LBin::combination(long n, long r){
     return(log_fac[n] - (log_fac[r] + log_fac[n-r]));
}

double LBin::log_binomial(long n, long r, double p){
     double x; 
     x = combination(n, r) + 1.0 * r*log(p);
     x= x + 1.0*(n-r)*log(1.0 -p);
     return (x);
}


double LBin::log_binomial_pval(long n, long r, double p){
     double x,y;
     long i;
     x= log_binomial(n, 0, p);
     for(i=1;i<=r;i++){ 
      y = log_binomial(n, i, p);
      x = addl(x,y);
     }
     return(x); 

}

double LBin::addl(double x,double y){
   double xt,u;
   u=x-y;
   if(u<=0){
      xt=exp(u);
      return(y+log(1.0+xt));
   }
   else {
      xt=exp(-1.0*u);
      return(x+log(1.0+xt));
   }
}

//Binomial Confidence Limit

Binomial::Binomial(double errs, double errp,long maxN){
  long i,j,k,m;
  double xx;

  eps=errs;
  epp=errp;
  m=maxN/1000+1;
  xfac=new double[m];
  xfac[0]=0;
  if(m==1){return;}
  xx=0;
  k=0;
  for(i=1;i<m;i++){
     for(j=1;j<=1000;j++){
        xx+=log((double)(k+j));
     }
     xfac[i]=xx;
     k+=1000;
  }
}

Binomial::~Binomial(void) {
}

double Binomial::fact(long m){
   long i,j,k,n;
   double xx,yy;

   i=m/1000;
   j=m%1000;
   n=m-j;
   xx=0;
   for(k=1;k<=j;k++){
      xx+=log((double)(n+k));
   }
   return(xfac[i]+xx);
}

double Binomial::combnr(long n,long r){
   return(fact(n)-fact(r)-fact(n-r));
}
 
double Binomial::binomial_pval_appx(long n, long r, double p){
     long i;
     double xx, yy,f;


     if(r) f=(1.0-p)/p;
     else {
       yy= exp(n*log(1.0 -p));
       return(yy);

     }

     yy=0.0;
     xx = combnr(n,r) + r*log(p);
     xx= xx + (n-r)*log(1.0 -p);
     xx= exp(xx);

     yy+=xx;
     i=r;
     while((i*xx)>=epp && i>0){
        xx=((i*f)/(double)(n-i+1))*xx;
        yy+=xx;
        i--;
     }
     return(yy);

}

double Binomial::upper_limit(long n, long r, double lim){
   long i;
   double p, qn, xx;
   double x,y,z;
   qn=lim;
   p=(double)r/(double)n;

   xx=this->binomial_pval_appx(n, r, p);
   if(xx<=qn){cout << "check the lim" << endl;return(xx);}
   if(n==r) return(p);
   x=p;
   y=1.0;
   while(y-x>eps){
      p=(y+x)/2;
      xx=this->binomial_pval_appx(n, r, p);
      if(qn<xx)x=p;
      else y=p;
   }
   return(x);
}

double Binomial::lower_limit(long n, long r, double lim){
   long i;
   double p,q,x;
   x=upper_limit(n,n-r,lim);
   return(1.0-x);
}

Wmw::Wmw(void){
   kn=sn=0;
   ka=sa=0;
   sr2=sqrt(2.0);
}

Wmw::~Wmw(void){
}

void Wmw::setup(long s,long k,long n){
   long i,j,flag=1;

   //Clean up from previous work
   if(ka){
      if((ka<k)||(sa<s)){
         for(i=1;i<ka+1;i++){
            delete [] sx[i];
            delete [] tx[i];
         }
         delete [] sx;
         delete [] tx;
      }
      else flag=0;
   }

   sn=s;
   kn=k;
   zn=n;

   ux=1.0;
   ct=sn-(kn*(kn-1))/2;
   ct=(zn<ct)?zn:ct;
   if(zn>ct){
      i=zn;
      while(i>ct){
         ux*=((double)(i-kn))/((double)i);
         i--;
      }
   }
   if(ct>=kn){
      if(flag){
         sx=new double*[kn+1];
         tx=new double*[kn+1];
         for(i=1;i<kn+1;i++){
            sx[i]=new double[sn+1];
            tx[i]=new double[sn+1];
         }
         ka=kn;
         sa=sn;
      }
      for(i=1;i<kn+1;i++){
         for(j=1;j<sn+1;j++){
            *(sx[i]+j)=*(tx[i]+j)=0.0;
         }
      }
      *(sx[kn]+sn)=ux;
   }
   wx=0.0;
}

void Wmw::calc(long n){
   long i,j,m,v;
   double xx,yy,zz;

   //k is 1 case
   m=((n<sn)?n:sn)+1;
   for(j=1;j<m;j++)wx+=*(sx[1]+j)*((double)j)/((double)n);
   for(j=m;j<sn+1;j++)wx+=*(sx[1]+j);
   //k is >1 case
   for(i=2;i<kn+1;i++){
      if(i>n)continue;
      m=(i*(i+1))/2;
      if(i==n){
         for(j=m;j<sn+1;j++)wx+=*(sx[i]+j);
      }
      else {
         xx=((double)(n-i))/((double)n);
         yy=((double)i)/((double)n);
         v=n+(i*(i-1))/2;
         for(j=m;j<v;j++){
            zz=*(sx[i]+j);
            if(zz){
               *(tx[i]+j)+=zz*xx;
            }
         }
         for(j=v;j<sn+1;j++){
            zz=*(sx[i]+j);
            if(zz){
               *(tx[i]+j)+=zz*xx;
               *(tx[i-1]+j-n)+=zz*yy;
            }
         }
      }
   }

   //Reset arrays sx and tx
   for(i=1;i<kn+1;i++){
      if(i>n)continue;
      for(j=1;j<sn+1;j++){
         *(sx[i]+j)=*(tx[i]+j);
         *(tx[i]+j)=0.0;
      }
   }
}

double Wmw::wmwcomp(long s,long k,long n){
   long i;

   setup(s,k,n);
   if(ct<kn)return(wx);
   i=ct;
   while(i>0){
      calc(i);
      i--;
   }
   return(wx);
}

double Wmw::wmwappx(long s,long k,long n){
   double xs=(double)s, xk=(double)k, xn=(double)n;
   double xx=(xs-xk*(xn+1)/2.0)/sqrt(xk*(xn-xk)*(xn+1)/12.0);
   if(xx<=0)return(0.5*erfc(-xx/sr2));
   else return(0.5+0.5*erf(xx/sr2));
}

double Wmw::wmwpval(long s,long k,long n){
   if((k>=10)&&(n-k>=10))return(wmwappx(s,k,n));
   else return(wmwcomp(s,k,n));
}

}

