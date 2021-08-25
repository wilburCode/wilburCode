#ifndef HUBR_H
#define HUBR_H

#include <iostream>
#include <fstream>
#include <CMark.h>

using namespace std;
namespace iret {

//template class HubeX
template<class Y,class Z>
class HubeX : public CMark<Y,Z> {
public:
   HubeX(void); //for copies
   HubeX(const char *nspost);//name of XPost set
   HubeX(const char *nspost,const char *pnam);//name of XPost set
      //pnam is used in "path_pnam" where to find path, or if ':'
      //is first char then remainder is the path
   ~HubeX(void);
   void SetMem(HubeX &Hx); //sets data and memory to Hx data and memory

   //Learning functions
   void Set_Docs(Indx<Y> *gdd,Indx<Y> *bdd); //Set cls array
     //gdd good docs, bdd bad docs
   void Set_Lambda_Norm(Z lam); //1.0E-7 is good lam for Rebase
   void Set_Lambda_Norm_lw(Z lam); //1.0E-5 is good lam for Rebase
   void Learn(void); //Optimizer
     //im is improvement threshold and iter the iteration limit
   void Learn_S(void); //Optimizer, skips 9/10 of time for already correct
     //is pointer at the partials array
     //is pointer at the partials array, xu multiplies weight of + class
   void Learn_lw(void); 
     //Uses local weights
   void Learn_S_lw(void); 
     //Uses local weights
   Z Func(void); //Computes and returns function value
   Z Func_lw(void); //Computes and returns function value
   Z Norm(void); //Computes the two norm of weight vector
     //without the th value included.
   Z *ScoreAll(void);
     //must call gopen_map() of XPost first 
   Z *ScoreSet(Indx<Y> *ind);
      //must call gopen_db_map() of XPost first
   Z *ScoreAll_lw(void);
     //must call gopen_map() of XPost first 
     //Uses local weights
   Z *ScoreSet_lw(Indx<Y> *ind);
      //must call gopen_db_map() of XPost first
     //Uses local weights
   void Save(int n); //Saves the set of weights wt. Marks with n
   void Load(int n); //Loads the set of weights wt marked with n
   void Release(int n); //Unmaps the weight vector mapped by load


   //Derivative function
   void Set_deriv(void); //Sets pdoc to derivative contribution for
      //that doc as it depends on rt.
   void Sift_terms(Y fl,Z th); //Looks for term pairs with
      //frequency of one of the terms at least fl and with a partial
      //derivative at least th in absolute value and prints into a file
      //with extension "dr".

   //Relate to the documents
   Y tdoc; //Number of training documents 
   Y *cls; //Marks the class with 1 or -1

   //Relate to the terms processing
   Z lambda; //lambda
   Z *wt; //weights of the process
   Z *wd; //derivative factor 
   Z *lm; //independent variables of function
   Z th; //Threshhold
   Y *sr;  //Maps the tdoc training set
   Z xnm; //norm of weight vector
   Z xf; //sum of violations
   Z tf; //True function value
 
   Z *sco;
};

template<class Y,class Z>
HubeX<Y,Z>::HubeX(void) : CMark<Y,Z>("null"){
   wt=NULL;
   wd=NULL;
   sr=NULL;
   sco=NULL;
   cls=NULL;
}

template<class Y,class Z>
HubeX<Y,Z>::HubeX(const char *namspost) : CMark<Y,Z>(namspost){
   wt=NULL;
   wd=NULL;
   sr=NULL;
   sco=NULL;
   cls=NULL;
}

template<class Y,class Z>
HubeX<Y,Z>::HubeX(const char *namspost,const char *pnam) : CMark<Y,Z>(namspost,pnam){
   wt=NULL;
   wd=NULL;
   sr=NULL;
   sco=NULL;
   cls=NULL;
}

template<class Y,class Z>
HubeX<Y,Z>::~HubeX(){
   //if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   if(sco!=NULL)delete [] sco;
   if(cls!=NULL)delete [] cls;
}

template<class Y,class Z>
void HubeX<Y,Z>::SetMem(HubeX<Y,Z> &Hx){
   XPost<Y,Z>::SetMem((XPost<Y,Z>&)Hx);
   wt=Hx.wt;
   th=Hx.th;
}

template<class Y,class Z>
void HubeX<Y,Z>::Set_Docs(Indx<Y> *gdd,Indx<Y> *bdd){
   Y i;

   if(cls!=NULL)delete [] cls;
   cls=new Y[this->ndoc];
   for(i=0;i<this->ndoc;i++)cls[i]=0;
   for(i=0;i<gdd->ix;i++){
      cls[gdd->idx[i]]=1;
   }
   for(i=0;i<bdd->ix;i++){
      cls[bdd->idx[i]]=-1;
   }
   tdoc=gdd->ix+bdd->ix;
}

template<class Y,class Z>
void HubeX<Y,Z>::Set_Lambda_Norm(Z lam){
   Y i,j,k; 
   Z sum=0,xx,yy,zz,cnx;
         
   this->gopen_db_map();
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         this->readp_db(i);
         xx=0;
         for(j=0;j<this->nw;j++){
            if(this->mrk[this->nwd[j]]){
               xx++;
            }
         }
         sum+=sqrt(xx);
      }
   }
   sum/=(Z)tdoc;
   lambda=lam*sum*sum;
   if(this->pflag)cout << "lambda= " << lambda << endl;
   lambda*=(Z)tdoc;
}

template<class Y,class Z>
void HubeX<Y,Z>::Set_Lambda_Norm_lw(Z lam){
   Y i,j,k;
   Z sum=0,xx,yy,zz,cnx;

   this->gopen_db_map();
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         this->readp_db(i);
         this->readz_db(i);
         xx=0;
         for(j=0;j<this->nw;j++){
            if(this->mrk[this->nwd[j]]){
               xx+=this->lwt[j]*this->lwt[j];
            }
         }
         sum+=sqrt(xx);
      }
   }
   sum/=(Z)tdoc;
   lambda=lam*sum*sum;
   if(this->pflag)cout << "lambda= " << lambda << endl;
   lambda*=(Z)tdoc;
}

template<class Y,class Z>
void HubeX<Y,Z>::Learn(void){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,mode,ct,gt,cu;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta;
   Z dth,f,of,xrch,xvt,ux,uy;
   Z fac,minf,rch,rfac;
   Z ss,dl,tho,*wto,llm;
   Z efu,efd,ominf,xfac;

   Store<Z> St(this->nwrd+1);
   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
   wd=new Z[this->nwrd];
   sr=new Y[tdoc];

   th=0;
   for(i=0;i<this->nwrd;i++){
      wt[i]=0;
   }

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         sr[n++]=i;
      }
   }
   llm=0;
   of=ominf=minf=1.0E100;
   fac=2.0/3.0;
   xfac=2.0/3.0;
   rfac=3.0/2.0;
   mode=0;
   efu=0;

   gt=ct=cu=0;
   flag=1;
 
 while(flag){
   gt++;ct++;
   //Deal with regularizer
   f=th*th;
   for(i=0;i<this->nwrd;i++){
      wd[i]=lambda*(xx=wt[i]);
      f+=xx*xx;
   }
   f*=0.5*lambda;
   dth=lambda*th;

   //Deal with main loss function
   for(i=0;i<tdoc;i++){
      j=sr[i];
      yy=cls[j];
      this->readp_db(j);
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         if(xx>-1.0){
            zz=1.0-xx;
            uu=-2.0*yy*zz;
            f+=zz*zz;
         }
         else {
            f+=-4.0*xx;
            uu=-4.0*yy;
         }
         dth+=uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wd[iu]+=uu;
            }
         }
      }
   }
   //Track the minimum f seen
   if(f>of)cu++;
   if(f<minf){minf=f;wt[this->nwrd]=th;St.CopyIn(wt);}
   if(cu==2){
      switch(mode){
         case 0: cu=ct=0;
                 if(minf==ominf){
                    mode=1;
                    fac*=xfac;
                 }
                 else ominf=minf;
                 break;
         case 1: efd=(ominf-minf)/ct;
                 ominf=minf;
                 mode=2;
                 fac*=rfac;
                 cu=ct=0;
                 break;
         case 2: efu=(ominf-minf)/ct;
                 ominf=minf;
                 cu=ct=0;
                 if(efd>=efu){
                    fac*=xfac;
                 }
                 mode=0;
      }
   }
   of=f;     

   del=f-minf+fac*minf;

   if(fac<0.05)break;
   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   eta=del/grad;
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];

   if(this->pflag)cout << endl << gt << "  " << cu << " mode " << mode << endl;
   if(this->pflag)cout << "    f " << f << " minf " << minf << " fac " << fac << endl;
   if(this->pflag)cout << " grad " << grad << " eta " << eta << endl; 
 }


   St.CopyOut(wt);
   th=wt[this->nwrd];
   delete [] wd;
   wd=NULL;
   delete [] sr;
   sr=NULL;
   if(this->pflag)cout << "Minf " << minf << endl;
}

template<class Y,class Z>
void HubeX<Y,Z>::Learn_S(void){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,flag2,gt,ct,cu,mode;
   int *wk,wkk;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta;
   Z dth,f,of,xrch,*wn,ux,uy;
   Z fac,minf,rfac,llm;
   Z efu,efd,ominf,xfac;


   Store<Z> St(this->nwrd+1);
   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
   wd=new Z[this->nwrd];
   sr=new Y[tdoc];
   wk=new int[tdoc];

   th=0;
   for(i=0;i<this->nwrd;i++)wt[i]=0;

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         wk[n]=0;
         sr[n++]=i;
      }
   } 
   of=ominf=minf=1.0E100;
   fac=2.0/3.0;
   xfac=2.0/3.0;
   rfac=3.0/2.0;
   mode=0;
   efu=0;

   gt=ct=cu=0;
   flag=1;

 while(flag){
   gt++;ct++;
   //Deal with regularizer
   f=th*th;
   for(i=0;i<this->nwrd;i++){
      wd[i]=lambda*(xx=wt[i]);
      f+=xx*xx;
   }
   f*=0.5*lambda;
   dth=lambda*th;

   //Deal with main loss function
   for(i=0;i<tdoc;i++){
      if((wkk=wk[i])>=0){
      j=sr[i];
      yy=cls[j];
      this->readp_db(j);
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         wk[i]=10;
         if(xx>-1.0){
            zz=1.0-xx;
            uu=-2.0*yy*zz;
            f+=zz*zz;
         }
         else {
            f+=-4.0*xx;
            uu=-4.0*yy;
         }
         dth+=uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wd[iu]+=uu;
            }
         }
      }
      else {
         if(wkk>0)wk[i]--;
         else wk[i]=-5-zrand(10);
      }
      }
      else wk[i]++;
   }
   //Track the minimum f seen
   if(f>of)cu++;
   if(f<minf){minf=f;wt[this->nwrd]=th;St.CopyIn(wt);}
   if(cu==2){
      switch(mode){
         case 0: cu=ct=0;
                 if(minf==ominf){
                    mode=1;
                    fac*=xfac;
                 }
                 else ominf=minf;
                 break;
         case 1: efd=(ominf-minf)/ct;
                 ominf=minf;
                 mode=2;
                 fac*=rfac;
                 cu=ct=0;
                 break;
         case 2: efu=(ominf-minf)/ct;
                 ominf=minf;
                 cu=ct=0;
                 if(efd>=efu){
                    fac*=xfac;
                 }
                 mode=0;
      }
   }
   of=f;

   del=f-minf+fac*minf;

   if(fac<0.05)break;
   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   eta=del/grad;
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];

   cout << endl << gt << "  " << cu << " mode " << mode << endl;
   cout << "    f " << f << " minf " << minf << " fac " << fac << endl;
   cout << " grad " << grad << " eta " << eta << endl;
 }

   St.CopyOut(wt);
   th=wt[this->nwrd];
   delete [] wd;
   wd=NULL;
   delete [] sr;
   sr=NULL;
   cout << "Minf " << minf << endl;
}

template<class Y,class Z>
void HubeX<Y,Z>::Learn_lw(void){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,flag2,ct,gt,cu,mode;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta;
   Z dth,f,of,xrch,*wn,ux,uy;
   Z fac,minf,rfac,llm;
   Z efu,efd,ominf,xfac;

   Store<Z> St(this->nwrd+1);
   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
   wd=new Z[this->nwrd];
   sr=new Y[tdoc];

   th=0;
   for(i=0;i<this->nwrd;i++)wt[i]=0;

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         sr[n++]=i;
      }
   }

   llm=0;
   of=ominf=minf=1.0E100;
   fac=2.0/3.0;
   xfac=2.0/3.0;
   rfac=3.0/2.0;
   mode=0;
   efu=0;

   gt=ct=cu=0;
   flag=1;
 
 while(flag){
   gt++;ct++;
   //Deal with regularizer
   f=th*th;
   for(i=0;i<this->nwrd;i++){
      wd[i]=lambda*(xx=wt[i]);
      f+=xx*xx;
   }
   f*=0.5*lambda;
   dth=lambda*th;

   //Deal with main loss function
   for(i=0;i<tdoc;i++){
      j=sr[i];
      yy=cls[j];
      this->readp_db(j);
      this->readz_db(j);
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu]*this->lwt[iz];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         if(xx>-1.0){
            zz=1.0-xx;
            uu=-2.0*yy*zz;
            f+=zz*zz;
         }
         else {
            f+=-4.0*xx;
            uu=-4.0*yy;
         }
         dth+=uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wd[iu]+=uu*this->lwt[iz];
            }
         }
      }
   }
   //Track the minimum f seen
   if(f>of)cu++;
   if(f<minf){minf=f;wt[this->nwrd]=th;St.CopyIn(wt);}
   if(cu==2){
      switch(mode){
         case 0: cu=ct=0;
                 if(minf==ominf){
                    mode=1;
                    fac*=xfac;
                 }
                 else ominf=minf;
                 break;
         case 1: efd=(ominf-minf)/ct;
                 ominf=minf;
                 mode=2;
                 fac*=rfac;
                 cu=ct=0;
                 break;
         case 2: efu=(ominf-minf)/ct;
                 ominf=minf;
                 cu=ct=0;
                 if(efd>=efu){
                    fac*=xfac;
                 }
                 mode=0;
      }
   }
   of=f;

   del=f-minf+fac*minf;

   if(fac<0.05)break;
   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   eta=del/grad;
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];

   if(this->pflag)cout << endl << gt << "  " << cu << " mode " << mode << endl;
   if(this->pflag)cout << "    f " << f << " minf " << minf << " fac " << fac << endl;
   if(this->pflag)cout << " grad " << grad << " eta " << eta << endl;
 }

   St.CopyOut(wt);
   th=wt[this->nwrd];
   delete [] wd;
   wd=NULL;
   delete [] sr;
   sr=NULL;
   if(this->pflag)cout << "Minf " << minf << endl;
}

template<class Y,class Z>
void HubeX<Y,Z>::Learn_S_lw(void){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,flag2,ct,gt,cu,mode;
   int *wk,wkk;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta;
   Z dth,f,of,xrch,*wn,ux,uy;
   Z fac,minf,rfac,llm;
   Z efu,efd,ominf,xfac;

   Store<Z> St(this->nwrd+1);
   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
   wd=new Z[this->nwrd];
   sr=new Y[tdoc];
   wk=new int[tdoc];

   th=0;
   for(i=0;i<this->nwrd;i++)wt[i]=0;

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         wk[n]=0;
         sr[n++]=i;
      }
   }

   llm=0;
   of=ominf=minf=1.0E100;
   fac=2.0/3.0;
   xfac=2.0/3.0;
   rfac=3.0/2.0;
   mode=0;
   efu=0;

   gt=ct=cu=0;
   flag=1;
 
 while(flag){
   gt++;ct++;
   //Deal with regularizer
   f=th*th;
   for(i=0;i<this->nwrd;i++){
      wd[i]=lambda*(xx=wt[i]);
      f+=xx*xx;
   }
   f*=0.5*lambda;
   dth=lambda*th;

   //Deal with main loss function
   for(i=0;i<tdoc;i++){
      if((wkk=wk[i])>=0){
      j=sr[i];
      yy=cls[j];
      this->readp_db(j);
      this->readz_db(j);
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu]*this->lwt[iz];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         wk[i]=10;
         if(xx>-1.0){
            zz=1.0-xx;
            uu=-2.0*yy*zz;
            f+=zz*zz;
         }
         else {
            f+=-4.0*xx;
            uu=-4.0*yy;
         }
         dth+=uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wd[iu]+=uu*this->lwt[iz];
            }
         }
         if(wkk<5)wk[i]=5;
      }
      else {
         if(wkk>0)wk[i]--;
         else wk[i]=-5-zrand(10);
      }
      }
      else wk[i]++;
   }
   //Track the minimum f seen
   if(f>of)cu++;
   if(f<minf){minf=f;wt[this->nwrd]=th;St.CopyIn(wt);}
   if(cu==2){
      switch(mode){
         case 0: cu=ct=0;
                 if(minf==ominf){
                    mode=1;
                    fac*=xfac;
                 }
                 else ominf=minf;
                 break;
         case 1: efd=(ominf-minf)/ct;
                 ominf=minf;
                 mode=2;
                 fac*=rfac;
                 cu=ct=0;
                 break;
         case 2: efu=(ominf-minf)/ct;
                 ominf=minf;
                 cu=ct=0;
                 if(efd>=efu){
                    fac*=xfac;
                 }
                 mode=0;
      }
   }
   of=f;

   del=f-minf+fac*minf;

   if(fac<0.05)break;
   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   eta=del/grad;
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];

   cout << endl << gt << "  " << cu << " mode " << mode << endl;
   cout << "    f " << f << " minf " << minf << " fac " << fac << endl;
   cout << " grad " << grad << " eta " << eta << endl;
 }

   St.CopyOut(wt);
   th=wt[this->nwrd];
   delete [] wd;
   wd=NULL;
   delete [] sr;
   sr=NULL;
   cout << "Minf " << minf << endl;
}

template<class Y,class Z>
Z *HubeX<Y,Z>::ScoreAll(void){
   Y j,n;
   Z sum;
   
   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];
      
   for(j=0;j<this->ndoc;j++){
      this->readp_db(j);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]];
      sco[j]=sum;
      this->mark(j+1,100,"docs scored");
   }
   return(sco);
}  

template<class Y,class Z>
Z *HubeX<Y,Z>::ScoreSet(Indx<Y> *ind){
   Y i,n;
   Z sum;

   if(sco!=NULL)delete [] sco;
   sco=new Z[ind->ix];

   for(i=0;i<ind->ix;i++){
      this->readp_db(ind->idx[i]);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]];
      sco[i]=sum;
      this->mark(i+1,100,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z *HubeX<Y,Z>::ScoreAll_lw(void){
   Y j,n;
   Z sum;

   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];

   for(j=0;j<this->ndoc;j++){
      this->readp_db(j);
      this->readz_db(j);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]]*this->lwt[n];
      sco[j]=sum;
      this->mark(j+1,100,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z *HubeX<Y,Z>::ScoreSet_lw(Indx<Y> *ind){
   Y i,j,n;
   Z sum;

   if(sco!=NULL)delete [] sco;
   sco=new Z[ind->ix];

   for(i=0;i<ind->ix;i++){
      j=ind->idx[i];
      this->readp_db(j);
      this->readz_db(j);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]]*this->lwt[n];
      sco[i]=sum;
      this->mark(i+1,100,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z HubeX<Y,Z>::Norm(void){
   Y i;
   Z xx,sw=th*th;
   for(i=0;i<this->nwrd;i++){
      if(this->mrk[i]){
         xx=wt[i];
         sw+=xx*xx;
      }
   }
   return(sqrt(sw));
}

template<class Y,class Z>
Z HubeX<Y,Z>::Func(void){
   Y i,iz,ui;
   Z sum,xx,yy,zz;
   xnm=th*th;
   xf=0;
   for(i=0;i<this->nwrd;i++){
      if(this->mrk[i]){
         xx=wt[i];
         xnm+=xx*xx;
      }
   }
   this->gopen_db_map();
   for(i=0;i<this->ndoc;i++){
      if(yy=cls[i]){
         sum=th;
         this->readp_db(i);
         for(iz=0;iz<this->nw;iz++){
            ui=this->nwd[iz];
            if(this->mrk[ui])sum+=wt[ui];
         }
         xx=yy*sum;
         if(xx<1.0){
            if(xx>-1.0)xf+=(1.0-xx)*(1.0-xx);
            else xf+=-4.0*xx;
         }
      }
   }
   tf=xf+0.5*lambda*xnm;
   return(tf);
}

template<class Y,class Z>
Z HubeX<Y,Z>::Func_lw(void){
   Y i,iz,ui;
   Z fs=0,sum,xx,yy,zz,sw=th*th;
   xnm=th*th;
   xf=0;
   for(i=0;i<this->nwrd;i++){
      if(this->mrk[i]){
         xx=wt[i];
         xnm+=xx*xx;
      }
   }  
   this->gopen_db_map();
   for(i=0;i<this->ndoc;i++){
      if(yy=cls[i]){
         sum=th;
         this->readp_db(i);
         this->readz_db(i);
         for(iz=0;iz<this->nw;iz++){
            ui=this->nwd[iz];
            if(this->mrk[ui])sum+=wt[ui]*this->lwt[iz];
         }
         xx=yy*sum;
         if(xx<1.0){
            if(xx>-1.0)xf+=(1.0-xx)*(1.0-xx);
            else xf+=-4.0*xx;
         }
      } 
   }
   tf=xf+0.5*lambda*xnm;
   return(tf);
}

template<class Y,class Z>
void HubeX<Y,Z>::Set_deriv(void){
   Y i,j,iz,iu,n,*sr;
   Z xx,yy,zz,uu,sum;

   sr=new Y[tdoc];
   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         sr[n++]=i;
      }
   }

   for(j=0;j<this->ndoc;j++){
      this->pdoc[j]=0;
   }

   for(i=0;i<tdoc;i++){
      j=sr[i];
      yy=cls[j];
      this->readp_db(j);
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         if(xx>-1.0){
            zz=1.0-xx;
            uu=-2.0*yy*zz;
         }
         else {
            uu=-4.0*yy;
         }
         this->pdoc[j]=uu;
      }
   }
   ofstream *pfout=this->get_Ostr("deriv",ios::out);
   pfout->write((char*)this->pdoc,this->ndoc*sizeof(Z));
   this->dst_Ostr(pfout);
}

template<class Y,class Z>
void HubeX<Y,Z>::Sift_terms(Y fl,Z th){
   Y i,j;
   Indx<Y> *pUnd;
   ofstream *pfout = this->get_Ostr("dr",ios::out);
   for(i=0;i<this->nwrd;i++){
      if(this->mrk[i]&&(this->freq[i]>=fl)){
         pUnd=this->readp(i);
         this->zerot();
         countTX(pUnd);
         for(j=0;j<i;j++){
            if((this->freq[j]<fl)&&(fabs(this->tx[j])>=th)){
               *pfout << this->tx[j] << " " << j << " " << i << endl;
            }
         }
         for(j=i+1;j<this->nwrd;j++){
            if(fabs(this->tx[j])>=th){
               *pfout << this->tx[j] << " " << i << " " << j << endl;
            }
         }
      }
      this->mark(i,10000,"terms");
   }
   this->dst_Ostr(pfout);
}

template<class Y,class Z>
void HubeX<Y,Z>::Save(int n){
   this->put_Nnum(n,"zh",this->ndoc,this->nwrd);
   this->bin_Writ(n,"weight",this->nwrd*sizeof(Z),(char*)wt);
   this->bin_Writ(n,"thresh",sizeof(Z),(char*)&th);
}

template<class Y,class Z>
void HubeX<Y,Z>::Load(int n){
   long i,j,k;
   this->get_Nnum(n,"zh",i,k);
   this->ndoc=i;
   this->nwrd=k;
   wt=(Z*)this->get_Mmap(n,"weight");
   ifstream *pfin=this->get_Istr(n,"thresh",ios::in);
   pfin->read((char*)&th,sizeof(Z));
   this->dst_Istr(pfin);
}

template<class Y,class Z>
void HubeX<Y,Z>::Release(int n){
    this->dst_Mmap(n,"weight",(char*&)wt);
    wt=NULL;
}

}
#endif
