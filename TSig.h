#ifndef TSIG_H
#define TSIG_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <DataObj.h>
#include <Hyper.h>
#include <XPost.h>

using namespace std;
namespace iret {

//template class TSig 
template<class Y,class Z>
class TSig : public XPost<Y,Z> {
public:
   TSig(const char *nspost);//name of XPost set
   TSig(const char *nspost,const char *pnam);//name of XPost set
      //pnam is either used in "path_pnam" as place to find path
      //or if begins with ':' is followed by path itself.
   ~TSig(void);

      //Counting functions.
   void init_cnt(void); //Sets pdoc uniform.
      //Sets memory for tx and sx.
      //must call gopen_db_map() of XPost first
   void countSX1(Indx<Y> *cnd); //adds all counts in cnd to s1.
   void countSX2(Indx<Y> *cnd); //adds all counts in cnd to s2.
   Z termValue(Y ui); //Divides postings and performs counts
      //zeroes ui's counts to remove it from the computation
      //then process counts to get term value and rezeroes all arrays
      //returns term value. Uses crossvalidation. 
   Z setValue(Indx<Y> *gnd); //Divides set and performs counts
      //then process counts to get term value and rezeroes all arrays
      //returns term value. Does crossvalication
   Z newValue(Y ui); //Uses postings and performs counts
      //processes terms counts to get term value
      //returns term value. Computes an expectation based on a
      //prior of 1 sig. term/250 terms expected.
   Z nsetValue(Indx<Y> *gnd); //preceding
      //except has no ui. Just uses the gnd. Makes
      //same prior=1/250 assumption.
   //Information measures
   double infoPair(Y ui,Y uj); //computes context information between the term pair
      //Same computation as done in the morphological paper
   double *infoRank(Y ui); //Produces the info (without taking log) for all terms
      //versus the term given by ui and returns a pointer at the score array of size nwrd 
   double infoPair2(Y ui,Y uj); //computes context information between the term pair
      //Computation as done in the morphological paper modified by HyperProb on 1st term.
   double *infoRank2(Y ui); //Produces the info (without taking log) as in inforPair2
      //versus the term given by ui and returns a pointer at the score array of size nwrd 
   //Mean relatedness
   double reldPair(Y ui,Y uj); //computes mean relatedness based on HyperG probs
   double *reldRank(Y ui); //Produces mean relatedness for all terms verus
      //a term given by ui and returns a pointer at the score array of size nwrd
   //Context overlap
   double contxPair(Y ui,Y uj,double &xi,double &xj,double &xc); //Returns the 2*xc/(xi+xj) (Dice coefficient)
      //xi is sum of squares for the first term HGProbs, xj for the second and xc is sum of product
      //of first term times second term
   void contxComp(Y ui,Y uj,double &xi,double &xj,double &xc); 
      //xi is newValue of ui, xj is newValue of uj and xc is sum of min 
      //of first term compared with second term summands
   void contxComp2(Y ui,Y uj,double &xi,double &xj,double &xc); 
      //xi is weight*prob from ui, xj is weight*prob from uj and xc is sum of min 
      //of first term compared with second term summands

   void Set_Term_wt(Z cut); //Set the mrk array 1 if wt>=cut
   void Set_Term_alp(Z cut); //Set the mrk array 1 if alpha>=cut
   void Set_Term_chi(Z cut); //Set the mrk array 1 if chi sq>=cut
   void Set_Term_muti(Z cut); //Set the mrk array 1 if muti>=cut
      //Used this on Rebase and found 0.02897 worked well (involved the
      //approx. 32k best terms).

   void Set_Term(Indx<Y> *pTrm);//Set the mrk array to process terms 
   			// listed in pTrm
   void Set_Term_freq(Y nm); //Set the mrk array
     //If term_freq<nm then 0, else 1. Based on count array tx.
   void Set_Term_freq2(Y nm); //Set the mrk array
     //If term_freq<nm then 0, else 1. Based on freq array
   void Set_Term_between(Y nm,Y bg); //Set the mrk array
     //If term_freq<nm or term_freq>bg then 0, else 1. Based on count array tx.
   void Set_Term_between2(Y nm,Y bg); //Set the mrk array
     //If term_freq<nm or term_freq>bg then 0, else 1. Based on freq array
   //Special mrk array functions
   void Set_All(int n); //Creates and sets mrk array to be n everywhere
   void Set_Char(char c,int n); //Sets mrk n if character c in string
   void Set_NChar(char c,int n); //Sets mrk n if character c not in string
   void Set_String(char *str,int n); //Sets mrk n if str substring of string
   void Set_NString(char *str,int n); //Sets mrk n if str not substring of string

   //Relate to the terms processing
   Z eps; //Size of 1/ndoc.
   Z nnx;
   Z nsx;
   Z *tx;
   Z *sx;
   Z *s1; //holds counts for the first set
   Y *m1; //Marks term numbers already seen
   Y *q1; //Que for term numbers seen in the count
   Y cq1; //Number of terms in q1
   Z *s2;
   Y *m2;
   Y *q2;
   Y cq2; //Number of terms in q2
   Y *mrk; //Marks which terms to include in process
   Y tmk;  //Number of words marked.
   Hyper *pHp; //for Hyper G calculations
};

template<class Y,class Z>
TSig<Y,Z>::TSig(const char *namspost) : XPost<Y,Z>(namspost){
   mrk=NULL;
   tx=NULL;
   s1=NULL;
   m1=NULL;
   q1=NULL;
   s2=NULL;
   m2=NULL;
   q2=NULL;
}

template<class Y,class Z>
TSig<Y,Z>::TSig(const char *namspost,const char *pnam) : XPost<Y,Z>(namspost,pnam){
   mrk=NULL;
   tx=NULL;
   s1=NULL;
   m1=NULL;
   q1=NULL;
   s2=NULL;
   m2=NULL;
   q2=NULL;
}

template<class Y,class Z>
TSig<Y,Z>::~TSig(){
   if(mrk!=NULL)delete [] mrk;
   if(tx!=NULL)delete [] tx;
   if(s1!=NULL)delete [] s1;
   if(m1!=NULL)delete [] m1;
   if(q1!=NULL)delete [] q1;
   if(s2!=NULL)delete [] s2;
   if(m2!=NULL)delete [] m2;
   if(q2!=NULL)delete [] q2;
}

template<class Y,class Z>
void TSig<Y,Z>::init_cnt(void){
   Y i;
   if(tx==NULL)tx=new Z[this->nwrd];
   if(s1==NULL)s1=new Z[this->nwrd];
   if(m1==NULL)m1=new Y[this->nwrd];
   if(q1==NULL)q1=new Y[this->nwrd];
   if(s2==NULL)s2=new Z[this->nwrd];
   if(m2==NULL)m2=new Y[this->nwrd];
   if(q2==NULL)q2=new Y[this->nwrd];
   for(Y i=0;i<this->nwrd;i++){
      *(tx+i)=this->freq[i];
      s1[i]=s2[i]=0;
      m1[i]=m2[i]=q1[i]=q2[i]=0;
   }
   pHp=new Hyper(this->ndoc);
}

template<class Y,class Z>
void TSig<Y,Z>::countSX1(Indx<Y> *cnd){
   Y i,j,k;

   cq1=0;
   for(i=0;i<cnd->ix;i++){
      this->readp_db(cnd->idx[i]);
      for(k=0;k<this->nw;k++){
         j=*(this->nwd+k);
         (*(s1+j))++;
         if(!m1[j]){
            m1[j]=1;
            q1[cq1++]=j;
         }
      }
      this->mark(i,1000,"docs");
   }
}

template<class Y,class Z>
void TSig<Y,Z>::countSX2(Indx<Y> *cnd){
   Y i,j,k;

   cq2=0;
   for(i=0;i<cnd->ix;i++){
      this->readp_db(cnd->idx[i]);
      for(k=0;k<this->nw;k++){
         j=*(this->nwd+k);
         (*(s2+j))++;
         if(!m2[j]){
            m2[j]=1;
            q2[cq2++]=j;
         }
      }
      this->mark(i,1000,"docs");
   }
}

template<class Y,class Z>
Z TSig<Y,Z>::termValue(Y ui){
   Y i,j,k,flag;
   Z xx,yy,zz,u1,u2,v1,v2,sum=0;
   Z nsx1,nsx2,nnx1,nnx2,n_t,n_st,max,min;
   Z pt,qt,wt;

   Indx<Y> *gnd=this->readp(ui);
   Indx<Y> *gnd1=gnd->Subsample(gnd->ix/2,987456);
   if(!gnd1){cout << "Error in gnd1!" << endl;exit(0);}
   Indx<Y> *gnd2=gnd->cbool_Butnot(gnd1);
   if(!gnd2){cout << "Error in gnd2!" << endl;exit(0);}
   countSX1(gnd1);
   countSX2(gnd2);
   s1[ui]=s2[ui]=0;
   u1=u2=0;
   
   nsx1=(Z)gnd1->ix;
   nsx2=(Z)gnd2->ix;
   nnx1=this->ndoc-nsx2;
   nnx2=this->ndoc-nsx1;
   Z frc1=nsx1/nnx1;
   Z frc2=nsx2/nnx2;
   Z diff=(Z)(gnd->ix-this->ndoc);
   for(i=0;i<cq1;i++){
      j=q1[i];
      n_t=tx[j]-s2[j];
      if(n_t<nnx1){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s1[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            u2+=wt*s2[j];
         }
      }
      n_t=tx[j]-s1[j];
      if(s2[j]&&(n_t<nnx2)){
         min=(n_t<nsx2)?n_t:nsx2;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s2[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc2)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc2)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx2;
         qt=(n_t-n_st)/(nnx2-nsx2);
         if(pt>qt){
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            u1+=wt*s1[j];
         }
      }
      if(s2[j]){
         s2[j]=m2[j]=0;
      }
      s1[j]=m1[j]=0;
   }
   for(i=0;i<cq2;i++){
      j=q2[i];
      if(!m2[j])continue;
      else s2[j]=m2[j]=0;
   }
   delete gnd1;
   delete gnd2;
   sum=0.5*(u1/nsx1+u2/nsx2);      
   return(sum);   
}

template<class Y,class Z>
Z TSig<Y,Z>::newValue(Y ui){
   Y i,j,k,flag,k1,k2;
   Z xx,yy,zz,u1,u2,v1,v2,sum=0;
   Z nsx1,nsx2,nnx1,nnx2,n_t,n_st,max,min;
   Z pt,qt,wt,vv;

   Indx<Y> *gnd=this->readp(ui);
   countSX1(gnd);
   s1[ui]=0;
   u1=0;
   k1=gnd->ix-1;
   k2=this->ndoc-1;
   
   nsx1=(Z)gnd->ix;
   nnx1=this->ndoc;
   Z frc1=nsx1/nnx1;
   Z diff=(Z)(nsx1-nnx1);
   for(i=0;i<cq1;i++){
      j=q1[i];
      n_t=tx[j];
      if((n_t<nnx1)&&(s1[j]>1.0)){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s1[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            vv=pHp->HlogOdds(rnd(s1[j])-1,k1,rnd(tx[j])-1,k2)-5.4;
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            u1+=wt*s1[j]/(1.0+pow(10.0,-vv));
         }
      }
      s1[j]=m1[j]=0;
   }
   return(u1/gnd->ix);   
}

template<class Y,class Z>
Z TSig<Y,Z>::setValue(Indx<Y> *gnd){
   Y i,j,k,flag;
   Z xx,yy,zz,u1,u2,v1,v2,sum=0;
   Z nsx1,nsx2,nnx1,nnx2,n_t,n_st,max,min;
   Z pt,qt,wt;

   Indx<Y> *gnd1=gnd->Subsample(gnd->ix/2,987456);
   if(!gnd1){cout << "Error in gnd1!" << endl;exit(0);}
   Indx<Y> *gnd2=gnd->cbool_Butnot(gnd1);
   if(!gnd2){cout << "Error in gnd2!" << endl;exit(0);}
   countSX1(gnd1);
   countSX2(gnd2);
   u1=u2=0;
   
   nsx1=(Z)gnd1->ix;
   nsx2=(Z)gnd2->ix;
   nnx1=this->ndoc-nsx2;
   nnx2=this->ndoc-nsx1;
   Z frc1=nsx1/nnx1;
   Z frc2=nsx2/nnx2;
   Z diff=(Z)(gnd->ix-this->ndoc);
   for(i=0;i<cq1;i++){
      j=q1[i];
      n_t=tx[j]-s2[j];
      if(n_t<nnx1){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s1[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            u2+=wt*s2[j];
         }
      }
      n_t=tx[j]-s1[j];
      if(s2[j]&&(n_t<nnx2)){
         min=(n_t<nsx2)?n_t:nsx2;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s2[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc2)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc2)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx2;
         qt=(n_t-n_st)/(nnx2-nsx2);
         if(pt>qt){
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            u1+=wt*s1[j];
         }
      }
      if(s2[j]){
         s2[j]=m2[j]=0;
      }
      s1[j]=m1[j]=0;
   }
   for(i=0;i<cq2;i++){
      j=q2[i];
      if(!m2[j])continue;
      else s2[j]=m2[j]=0;
   }
   delete gnd1;
   delete gnd2;
   sum=0.5*(u1/nsx1+u2/nsx2);      
   return(sum);   
}

template<class Y,class Z>
Z TSig<Y,Z>::nsetValue(Indx<Y> *gnd){
   Y i,j,k,flag,k1,k2;
   Z xx,yy,zz,u1,u2,v1,v2,sum=0;
   Z nsx1,nsx2,nnx1,nnx2,n_t,n_st,max,min;
   Z pt,qt,wt,vv;

   countSX1(gnd);
   u1=0;
   k1=gnd->ix-1;
   k2=this->ndoc-1;

   nsx1=(Z)gnd->ix;
   nnx1=this->ndoc;
   Z frc1=nsx1/nnx1;
   Z diff=(Z)(nsx1-nnx1);
   for(i=0;i<cq1;i++){
      j=q1[i];
      n_t=tx[j];
      if((n_t<nnx1)&&(s1[j]>1.0)){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s1[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            vv=pHp->HlogOdds(rnd(s1[j])-1,k1,rnd(tx[j])-1,k2)-5.4;
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            u1+=wt*s1[j]/(1.0+pow(10.0,-vv));
         }
      }
      s1[j]=m1[j]=0;
   }
   return(u1/gnd->ix);
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_wt(Z cut){
   Z n_t,n_st,min;
   Y nstw=0, flag;
   Z xx,wtt,frc,pt,qt,rtx;
   if(mrk==NULL)mrk=new Y[this->nwrd];

   frc=nsx/nnx;

   for(Y i=0;i<this->nwrd;i++){
      n_t=*(tx+i);
      if(n_t){
         min=(n_t<nsx)?n_t:nsx;
         n_st=*(sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-eps >n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st==0){
            if(eps<n_t*frc)n_st=eps;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         pt =n_st/nsx;
         qt =(n_t-n_st)/(nnx-nsx);
         //calculate wt for this term
         wtt=log(pt*(1.0-qt))-log(qt*(1.0-pt));

         if(fabs(wtt)>=cut){
            mrk[i]=1;
            nstw++;
         }
         else mrk[i]=0;
      }
      else {
         mrk[i]=0;
      }
   }
   if(this->pflag)cout << "Number of weight marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_alp(Z cut){
   Z n_t,n_st,min;
   Y nstw=0, flag;
   Z xx,frc,pt,qt,rtx;
   if(mrk==NULL)mrk=new Y[this->nwrd];

   frc=nsx/nnx;

   for(Y i=0;i<this->nwrd;i++){
      n_t=*(tx+i);
      if(n_t){
         min=(n_t<nsx)?n_t:nsx;
         n_st=*(sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-eps >n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st==0){
            if(eps<n_t*frc)n_st=eps;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         pt =n_st/nsx;
         qt =(n_t-n_st)/(nnx-nsx);
         rtx =n_t/nnx;

         xx=-n_st*log(rtx/pt) - (n_t-n_st)*log(rtx/qt);
         xx-=(nsx-n_st)*log((1.0-rtx)/(1.0-pt));
         xx-=(nnx-nsx-n_t+n_st)*log((1.0-rtx)/(1.0-qt));
         if(xx>=cut){
            mrk[i]=1;
            nstw++;
         }
         else mrk[i]=0;

      }
      else {
         mrk[i]=0;
      }
   }
   if(this->pflag)cout << "Number of alpha marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_chi(Z cut){
   Z n_t,n_st,min;
   Y nstw=0, flag;
   Z xx,yy,max,zz,uu,frc,thr,diff;
   Z sum;
   if(mrk==NULL)mrk=new Y[this->nwrd];

   frc=nsx/nnx;
   diff=nsx-nnx;

   for(Y i=0;i<this->nwrd;i++){
      n_t=*(tx+i);
      if(n_t&&(n_t<nnx)){
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-eps>n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+eps<n_t*frc)n_st=max+eps;
            else flag=0; 
         }
      }
      else flag=0;

      if(flag){
         xx=frc*n_t;
         thr=eps*cut*xx*(1.0-xx/((Z)nnx));
         xx-=(Z)n_st;
         if(xx*xx>=thr){
            mrk[i]=1;
            nstw++;
         }
         else mrk[i]=0;
      }
      else {
         mrk[i]=0;
      }
   }
   if(this->pflag)cout << "Number of chi sq marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_muti(Z cut){
   Z n_t,n_st,min;
   Y nstw=0, flag;
   Z xx,yy,max,zz,uu,frc,thr,diff;
   Z sum;
   if(mrk==NULL)mrk=new Y[this->nwrd];

   frc=nsx/nnx;
   diff=nsx-nnx;

   for(Y i=0;i<this->nwrd;i++){
      n_t=*(tx+i);
      if(n_t&&(n_t<nnx)){
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-eps>n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+eps<n_t*frc)n_st=max+eps;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         sum=n_st*(log(nnx*n_st)-log(n_t*nsx));
         sum+=(nsx-n_st)*(log(nnx*(nsx-n_st)-log((nnx-n_t)*nsx)));
         sum+=(n_t-n_st)*(log(nnx*(n_t-n_st))-log(n_t*(nnx-nsx)));
         sum+=(nnx-nsx-n_t+n_st)*(log(nnx*(nnx-nsx-n_t+n_st))-log((nnx-nsx)*(nnx-n_t)));
         if(sum>=cut){
            mrk[i]=1;
            nstw++;
         }
         else mrk[i]=0;
      }
      else {
         mrk[i]=0;
      }
   }
   if(this->pflag)cout << "Number of mutual info marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term(Indx<Y> *pTrm){
   Y i;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++) mrk[i]=0;
   
   for(i=0; i<pTrm->ix; i++) mrk[pTrm->idx[i]]=1;
   if(this->pflag)cout << "Number of marked terms= " << pTrm->ix << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_freq(Y nm){
   Y i,nstw=0;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if(tx[i]<nm)mrk[i]=0;
      else {
         mrk[i]=1;
         nstw++;
      }
   }
   if(this->pflag)cout << "Number of freq marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_freq2(Y nm){
   Y i,nstw=0;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if(this->freq[i]<nm)mrk[i]=0;
      else {
         mrk[i]=1;
         nstw++;
      }
   }
   if(this->pflag)cout << "Number of freq marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_between(Y nm,Y bg){
   Y i,nstw=0;
   Z xx=eps*nm,yy=eps*bg;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if((tx[i]<xx)||(tx[i]>yy))mrk[i]=0;
      else {
         mrk[i]=1;
         nstw++;
      }
   }
   if(this->pflag)cout << "Number of freq marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Term_between2(Y nm,Y bg){
   Y i,nstw=0;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if((this->freq[i]<nm)||(this->freq[i]>bg))mrk[i]=0;
      else {
         mrk[i]=1;
         nstw++;
      }
   }
   if(this->pflag)cout << "Number of freq marked terms= " << nstw << endl;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_All(int n){
   Y i;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++)mrk[i]=n;
}

template<class Y,class Z>
void TSig<Y,Z>::Set_Char(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strchr(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void TSig<Y,Z>::Set_String(char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strstr(this->show(i),str))mrk[i]=n;
   }
}

template<class Y,class Z>
void TSig<Y,Z>::Set_NChar(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strchr(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void TSig<Y,Z>::Set_NString(char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strstr(this->show(i),str))mrk[i]=n;
   }
}

template<class Y,class Z>
double TSig<Y,Z>::infoPair(Y ui,Y uj){
   long i,j,k,u;
   countSX1(this->readp(ui));
   countSX2(this->readp(uj));
   i=0;
   for(j=0;j<this->nwrd;j++){
      i+=this->freq[j];
   }
   double ttx=(double)i;
   double pr1,pr2,prb=0;
   i=j=0;
   for(k=0;k<cq1;k++){  
      u=q1[k];
      i+=rnd(s1[u]);
      prb+=(s1[u]/this->freq[u])*s2[u];
      s1[u]=0;
      m1[u]=0;
   }
   pr1=(double)i;  
   for(k=0;k<cq2;k++){  
      u=q2[k];
      j+=rnd(s2[u]);
      s2[u]=0;
      m2[u]=0;
   }
   pr2=(double)j;  
   return(log(prb*ttx/(pr1*pr2)));
}

template<class Y,class Z>
double* TSig<Y,Z>::infoRank(Y uj){
   long i,j,k,fq,mq,iq,*fqq,*irb,sb,ui;
   double zz,xx,yy,sy,su,dx,ttx,*pv,*exc,eps=0.01,alp;
   double *pr1,*prx,*pry;

   pr1=new double[this->nwrd];
   prx=new double[this->nwrd];
   pry=new double[this->nwrd];

   sb=0;
   for(ui=0;ui<this->nwrd;ui++){
      prx[ui]=pry[ui]=0;
      sb+=this->freq[ui];
   }
   ttx=(double)sb;
   countSX1(this->readp(uj));
   for(ui=0;ui<this->nwrd;ui++){
      pr1[ui]=s1[ui]/this->freq[ui];
      if(m1[ui]){s1[ui]=0;m1[ui]=0;}
   }
   cout << "Setup complete" << endl;
   for(ui=0;ui<this->ndoc;ui++){
      this->readp_db(ui);
      for(i=0;i<this->nw;i++){
         k=this->nwd[i];
         prx[k]+=this->nw;
         pry[k]+=pr1[k];
         for(j=i+1;j<this->nw;j++){
            pry[k]+=pr1[this->nwd[j]];
            pry[this->nwd[j]]+=pr1[k];
         }
      }
      mark(1,ui,1000,"docs");
   }
   dx=prx[uj];
   for(ui=0;ui<this->nwrd;ui++){
      pry[ui]=pry[ui]*ttx/(dx*prx[ui]);
   }
   delete [] pr1;
   delete [] prx;
   return(pry);
}

template<class Y,class Z>
double TSig<Y,Z>::infoPair2(Y ui,Y uj){
   long i,j,k,m,u;
   double xx,yy,zz;

   m=this->freq[ui];
   countSX1(this->readp(ui));
   countSX2(this->readp(uj));
   double *phg=new double[this->nwrd];
   double ttx=0;
   for(j=0;j<this->nwrd;j++){
      yy=pHp->HlogOdds(rnd(s1[j]),m,this->freq[j],this->ndoc);
      phg[j]=1.0/(1.0+pow(10.0,-yy+5.4));
      ttx+=phg[j]*this->freq[j];
   }
   double pr1,pr2,prb=0;
   j=0;
   pr1=0;
   for(k=0;k<cq1;k++){
      u=q1[k];
      pr1+=yy=s1[u]*phg[u];
      prb+=(yy/this->freq[u])*s2[u];
      s1[u]=0;
      m1[u]=0;
   }
   pr2=0;
   for(k=0;k<cq2;k++){
      u=q2[k];
      pr2+=s2[u]*phg[u];
      s2[u]=0;
      m2[u]=0;
   }
   return(log(prb*ttx/(pr1*pr2)));
}

template<class Y,class Z>
double* TSig<Y,Z>::infoRank2(Y uj){
   long i,j,k,fq,mq,iq,*fqq,*irb,sb,ui,m;
   double zz,xx,yy,sy,su,dx,ttx,*pv,*exc,eps=0.01,alp;
   double *pr1,*prx,*pry;

   m=this->freq[uj];
   pr1=new double[this->nwrd];
   prx=new double[this->nwrd];
   pry=new double[this->nwrd];
   double *phg=new double[this->nwrd];

   countSX1(this->readp(uj));
   ttx=0;
   for(ui=0;ui<this->nwrd;ui++){
      prx[ui]=pry[ui]=0;
      yy=pHp->HlogOdds(rnd(s1[ui]),m,this->freq[ui],this->ndoc);
      phg[ui]=1.0/(1.0+pow(10.0,-yy+5.4));
      ttx+=phg[ui]*this->freq[ui];
      pr1[ui]=phg[ui]*s1[ui]/this->freq[ui];
      if(m1[ui]){s1[ui]=0;m1[ui]=0;}
   }
   cout << "Setup complete" << endl;
   for(ui=0;ui<this->ndoc;ui++){
      this->readp_db(ui);
      for(i=0;i<this->nw;i++){
         k=this->nwd[i];
         prx[k]+=phg[k];
         pry[k]+=pr1[k];
         for(j=i+1;j<this->nw;j++){
            pry[k]+=pr1[this->nwd[j]];
            pry[this->nwd[j]]+=pr1[k];
            prx[k]+=phg[this->nwd[j]];
            prx[this->nwd[j]]+=phg[k];
         }
      }
      mark(1,ui,1000,"docs");
   }
   dx=prx[uj];
   for(ui=0;ui<this->nwrd;ui++){
      pry[ui]=pry[ui]*ttx/(dx*prx[ui]);
   }
   delete [] pr1;
   delete [] prx;
   return(pry);
}

template<class Y,class Z>
double TSig<Y,Z>::reldPair(Y ui,Y uj){
   long i,j,k,u,m,n;
   Indx<Y> *pInd=this->readp(ui);
   m=pInd->ix;
   countSX1(pInd);
   pInd=this->readp(uj);
   n=pInd->ix;
   countSX2(pInd);
   double pr2,prb=0,xx,yy;
   j=0;
   for(k=0;k<cq2;k++){  
      u=q2[k];
      j=rnd(s1[u]);
      xx=pHp->HlogOdds(j,m,this->freq[u],this->ndoc);
      yy=1.0/(1.0+pow(10.0,-xx+5.4));
      prb+=yy*s2[u];
      s2[u]=0;
      m2[u]=0;
   }
   pr2=(double)(n+m);  
   for(k=0;k<cq1;k++){  
      u=q1[k];
      s1[u]=0;
      m1[u]=0;
   }
   return(prb/pr2);
}

template<class Y,class Z>
double* TSig<Y,Z>::reldRank(Y uj){
   long i,j,k,u,ui,m,n;
   double zz,xx,yy,sy,su,dx,ttx,*pv,*exc,eps=0.01,alp;
   double *pr1,*prx,*pry;

   pr1=new double[this->nwrd];
   pry=new double[this->nwrd];

   for(ui=0;ui<this->nwrd;ui++){
      pr1[ui]=pry[ui]=0;
   }
   Indx<Y> *pInd=this->readp(uj);
   m=pInd->ix;
   countSX1(pInd);
   for(k=0;k<cq1;k++){  
      u=q1[k];
      j=rnd(s1[u]);
      xx=pHp->HlogOdds(j,m,this->freq[u],this->ndoc);
      pr1[u]=1.0/(1.0+pow(10.0,-xx+5.4));
      s1[u]=0;m1[u]=0;
   }
   cout << "Setup complete" << endl;
   for(ui=0;ui<this->ndoc;ui++){
      this->readp_db(ui);
      for(i=0;i<this->nw;i++){
         k=this->nwd[i];
         pry[k]+=pr1[k];
         for(j=i+1;j<this->nw;j++){
            pry[k]+=pr1[this->nwd[j]];
            pry[this->nwd[j]]+=pr1[k];
         }
      }
      mark(1,ui,1000,"docs");
   }
   for(ui=0;ui<this->nwrd;ui++){
      pry[ui]=pry[ui]/(this->freq[ui]+m);
   }
   delete [] pr1;
   return(pry);
}

template<class Y,class Z>
double TSig<Y,Z>::contxPair(Y ui,Y uj,double &xi,double &xj,double &xc){
   long i,j,k,m,n,u;
   double xx,yy,zz;
   xi=xj=xc=0;

   m=this->freq[ui];
   n=this->freq[uj];
   countSX1(this->readp(ui));
   countSX2(this->readp(uj));
   for(j=0;j<this->nwrd;j++){
      yy=pHp->HlogOdds(rnd(s1[j]),m,this->freq[j],this->ndoc);
      xx=1.0/(1.0+pow(10.0,-yy+5.4));
      yy=pHp->HlogOdds(rnd(s2[j]),n,this->freq[j],this->ndoc);
      zz=1.0/(1.0+pow(10.0,-yy+5.4));
      xi+=xx*xx;
      xj+=zz*zz;
      xc+=xx*zz;
      if(m1[j]){m1[j]=0;s1[j]=0;}
      if(m2[j]){m2[j]=0;s2[j]=0;}
   }
   return(2.0*xc/(xi+xj));
}

template<class Y,class Z>
void TSig<Y,Z>::contxComp(Y ui,Y uj,double &xi,double &xj,double &xc){
   Y i,j,k,flag,k1,k2;
   Z xx,yy,zz,u1,u2,v1,v2,sum=0;
   Z nsx1,nsx2,nnx1,nnx2,n_t,n_st,max,min;
   Z pt,qt,wt,vv;

   Indx<Y> *gnd=this->readp(ui);
   countSX1(gnd);
   xi=xj=xc=0;
   k1=gnd->ix-1;
   k2=this->ndoc-1;

   nsx1=(Z)gnd->ix;
   nnx1=this->ndoc;
   Z frc1=nsx1/nnx1;
   Z diff=(Z)(nsx1-nnx1);
   for(i=0;i<cq1;i++){
      j=q1[i];
      n_t=tx[j];
      if((n_t<nnx1)&&(s1[j]>1.0)){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s1[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            vv=pHp->HlogOdds(rnd(s1[j])-1,k1,rnd(tx[j])-1,k2)-2.4;
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            xi+=s1[j]=wt*s1[j]/(gnd->ix*(1.0+pow(10.0,-vv)));
         }
         else s1[j]=0;
      }
      else s1[j]=0;
   }
   gnd=this->readp(uj);
   countSX2(gnd);
   k1=gnd->ix-1;
   k2=this->ndoc-1;

   nsx1=(Z)gnd->ix;
   nnx1=this->ndoc;
   frc1=nsx1/nnx1;
   diff=(Z)(nsx1-nnx1);
   for(i=0;i<cq2;i++){
      j=q2[i];
      n_t=tx[j];
      if((n_t<nnx1)&&(s2[j]>1.0)){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s2[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            vv=pHp->HlogOdds(rnd(s2[j])-1,k1,rnd(tx[j])-1,k2)-2.4;
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            xj+=u1=wt*s2[j]/(gnd->ix*(1.0+pow(10.0,-vv)));
            xc+=(u1<s1[j])?u1:s1[j];
         }
      }
      s2[j]=m2[j]=0;
   }
   for(i=0;i<cq1;i++){
      j=q1[i];
      s1[j]=m1[j]=0;
   }
}


template<class Y,class Z>
void TSig<Y,Z>::contxComp2(Y ui,Y uj,double &xi,double &xj,double &xc){
   Y i,j,k,flag,k1,k2;
   Z xx,yy,zz,u1,u2,v1,v2,sum=0;
   Z nsx1,nsx2,nnx1,nnx2,n_t,n_st,max,min;
   Z pt,qt,wt,vv;

   Indx<Y> *gnd=this->readp(ui);
   countSX1(gnd);
   xi=xj=xc=0;
   k1=gnd->ix-1;
   k2=this->ndoc-1;

   nsx1=(Z)gnd->ix;
   nnx1=this->ndoc;
   Z frc1=nsx1/nnx1;
   Z diff=(Z)(nsx1-nnx1);
   for(i=0;i<cq1;i++){
      j=q1[i];
      n_t=tx[j];
      if((n_t<nnx1)&&(s1[j]>1.0)){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s1[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            vv=pHp->HlogOdds(rnd(s1[j])-1,k1,rnd(tx[j])-1,k2)-2.4;
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            xi+=s1[j]=wt/(1.0+pow(10.0,-vv));
         }
         else s1[j]=0;
      }
      else s1[j]=0;
   }
   gnd=this->readp(uj);
   countSX2(gnd);
   k1=gnd->ix-1;
   k2=this->ndoc-1;

   nsx1=(Z)gnd->ix;
   nnx1=this->ndoc;
   frc1=nsx1/nnx1;
   diff=(Z)(nsx1-nnx1);
   for(i=0;i<cq2;i++){
      j=q2[i];
      n_t=tx[j];
      if((n_t<nnx1)&&(s2[j]>1.0)){
         min=(n_t<nsx1)?n_t:nsx1;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=s2[j];
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc1)n_st-=1.0;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc1)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         pt=n_st/nsx1;
         qt=(n_t-n_st)/(nnx1-nsx1);
         if(pt>qt){
            vv=pHp->HlogOdds(rnd(s2[j])-1,k1,rnd(tx[j])-1,k2)-2.4;
            wt=logf(pt*(1.0-qt))-logf(qt*(1.0-pt));
            xj+=u1=wt/(1.0+pow(10.0,-vv));
            xc+=(u1<s1[j])?u1:s1[j];
         }
      }
      s2[j]=m2[j]=0;
   }
   for(i=0;i<cq1;i++){
      j=q1[i];
      s1[j]=m1[j]=0;
   }
}

}

#endif
