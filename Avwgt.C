#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "Avwgt.h"
#include "Exmat.h"
using namespace std;
namespace iret {

Avwgt::Avwgt(long n,long gn) {
   nobj=n;
   gobj=gn;
   long z=(gn<n-gn)?n-gn:gn;
   long zz=(long)ceiling(z/6.0);
   ntp=(z<zz)?zz:z;
   z=(long)ceiling(100.0*log10((double)z));
   ntp=(ntp<z)?z:ntp;
}

Avwgt::~Avwgt() {
}

Exmat *Avwgt::onedr(long alp,long bet){
   long i,j,k;
   double xx,yy,zz;
   Exmat *pS=new Exmat(1000,ntp,20);
   Exmat A(1000,ntp,20);
   Exmat B(1000,ntp,20);
   Exmat C(1000,ntp,20);
   A.load(alp);
   A.mul(alp+1);
   B.load(1);
   for(i=alp+2;i<alp+bet;i++){
      A.mul(i);
      B.mul(i-alp);
      if(!(i%100)){
         A.div(B);
         B.load(1);
      }
   }
   if((i-1)%100){
      A.div(B);
   }
   A.mul(-1);
   B.copy(A);
   C.load(alp*alp);
   B.div(C);
   pS->copy(B);
   for(i=1;i<bet;i++){
      A.mul(i-bet);
      C.load(i);
      A.div(C);
      B.copy(A);
      C.load((alp+i)*(alp+i));
      B.div(C);
      pS->add(B);
   }
   pS->debug("final");
   return(pS);
}

double Avwgt::ave_wgt(long np,long nq){
   long nnp=gobj-np;
   long nnq=(nobj-gobj)-nq;
   Exmat *pS1=onedr(np,nnp);
   Exmat *pS2=onedr(nnq,nq);
   Exmat *pT1=onedr(nnp,np);
   Exmat *pT2=onedr(nq,nnq);
   pS1->add(*pS2);
   pS1->sbt(*pT1);
   pS1->sbt(*pT2);
   pS1->debug("sum");
   double xx=pS1->val();
   delete pS1;
   delete pS2;
   delete pT1;
   delete pT2;
   return(xx);
}

double Avwgt::mxl_wgt(long np,long nq){
   long nnp=gobj-np;
   long nnq=(nobj-gobj)-nq;
   double psm=((double)np)/((double)(np+nnp));
   double psk=((double)nq)/((double)(nq+nnq));
   return(log(psm*(1.0-psk)/(psk*(1.0-psm))));
}

}
