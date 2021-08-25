#ifndef SPMRK_H
#define SPMRK_H
#include <fstream>
#include <iostream>
#include <DataObj.h>
#include <Heap.h>
#include <runn.h>

using namespace std;
namespace iret {

template<class Y,class Z>
class  SpMrk {
   public:
      SpMrk(long n); //Array size
     ~SpMrk(void);
       //Assumes CreateLenCos called before
   Ordr<Y,Z> *Skim(Y n); //Skims off the top n items 
       //but not more than sck
   Ordr<Y,Z> *SkimC(Z Th); //Skims off items scoring >=Th 

   inline void SetVal(Y ix,Z sx){
      xx[ix]=sx;
      if(!mrk[ix]){
         stk[sck++]=ix;
         mrk[ix]=true;
      }
   }
   inline void AddVal(Y ix,Z sx){
      xx[ix]+=sx;
      if(!mrk[ix]){
         stk[sck++]=ix;
         mrk[ix]=true;
      }
   }
   inline void Zero(void){
      long i;
      for(i=0;i<sck;i++){
         stk[i]=0;
         mrk[i]=false;
      }
      sck=0;
   }

   //Data
   long ns; //Array sizes
   Z *xx; //Data array
   Y *stk; //stack of visited indices
   long sck; //stack counter
   bool *mrk; //array to mark visited indices
};

template<class Y,class Z>
SpMrk<Y,Z>::SpMrk(long n){
   long i;
   ns=n;
   xx=new Z[ns];
   stk=new Y[ns];
   mrk = new bool[ns];
   for(i=0;i<ns;i++){
      xx[i]=0;
      mrk[i]=false;
   }
}

template<class Y,class Z>
SpMrk<Y,Z>::~SpMrk(void){
   delete [] xx;
   delete [] stk;
   delete [] mrk;
}

template<class Y,class Z> 
Ordr<Y,Z> *SpMrk<Y,Z>::Skim(Y n){
   Y i,j,k;
   Z *ra=new Z[sck];
   Y *rb=new Y[sck];
   Y *ord=new Y[sck];
   Heap2<Z,Y> H2(n,ra,rb);
   for(i=0;i<sck;i++){
      j=stk[i];
      H2.add_pointF(xx[j],j);
      ord[i]=i;
   }
   k=H2.orderF();
   hRort(k,rb,ord);
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>;
   pOrd->ix=k;
   pOrd->idx=rb;
   pOrd->order=ord;
   pOrd->score=ra;
   Y *ax=pOrd->invert();
   pOrd->order=ax;
   delete [] ord;
   return(pOrd);
}

template<class Y,class Z> 
Ordr<Y,Z> *SpMrk<Y,Z>::SkimC(Z Th){
   Y i,j,k,m=0;
   Z *ra=new Z[sck];
   Y *rb=new Y[sck];
   Y *ord=new Y[sck];
   Heap2<Z,Y> H2(sck,ra,rb);
   for(i=0;i<sck;i++){
      j=stk[i];
      if(xx[j]>=Th){
         H2.add_pointF(xx[j],j);
         ord[m]=m++;
      }
   }
   k=H2.orderF();
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>;
   if(!k){
      pOrd->ix=k;
      pOrd->idx=rb;
      pOrd->order=ord;
      pOrd->score=ra;
      return(pOrd);
   }
   hRort(k,rb,ord);
   pOrd->ix=k;
   pOrd->idx=rb;
   pOrd->order=ord;
   pOrd->score=ra;
   Y *ax=pOrd->invert();
   pOrd->order=ax;
   delete [] ord;
   return(pOrd);
}

}

#endif
