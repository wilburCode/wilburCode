#ifndef SVMG_H
#define SVMG_H

#include <iostream>
#include <fstream>
#include <DataObj.h>
#include <XPost.h>
#include <Hyper.h>
#include <Elev.h>
#include <CMark.h>

using namespace std;
namespace iret {

//template class SvmGX
template<class Y,class Z>
class SvmGX : public CMark<Y,Z> {
public:
   SvmGX(const char *nspost);//name of XPost set
   SvmGX(const char *nspost,const char *pnam);//name of XPost set
      //pnam is either used in "path_pnam" as place to find path
      //or if begins with ':' is followed by path itself.
   ~SvmGX(void);

   //Learning functions
   void Set_Docs(Indx<Y> *gdd,Indx<Y> *bdd); //Set cls array
     //gdd good docs, bdd bad docs

   void LearnSGD(Y iter,Y ui,CValdX<Y> *Cval); //Optimization based on lam=0 and early stopping
     //iter is number of learning cycles, ui is branch to learn and test
   void LearnSGDTest(Y iter,Y ui,CValdX<Y> *Cval); //Optimization based on lam=0 and early stopping
     //iter is number of learning cycles, ui is branch to learn and test
   void LearnSGD_lw(Y iter,Y ui,CValdX<Y> *Cval); //Optimization based on lam=0 and early stopping
     //iter is number of learning cycles, ui is branch to learn and test
   void LearnSGDO(Y iter,Z neta,Z peta,CValdX<Y> *Cval); //Optimization based on lam=0 and early stopping
     //Assumes ui=0. Experimental
   void LearnSGD(Y iter); //Simple SGD like update. Optimize with lam=0 and early stopping
     //iter is number of cycles. Find iter=8 is good for large MEDLINE collections, iter=9 for smaller
     //training sets.
   void LearnSGD_lw(Y iter); //Simple SGD like update. Optimize with lam=0 and early stopping
     //iter is number of cycles. 
   void DocAnal(Y nf); //nf is file number for id purposes. Analyzes the training docs and records
     //sizes of descisive sets for all in a file.

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

   //Debug functions
   double Debug_Func(void); //Computes and returns function value
      //Puts out the partial sums

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
SvmGX<Y,Z>::SvmGX(const char *namspost) : CMark<Y,Z>(namspost){
   wt=NULL;
   wd=NULL;
   sr=NULL;
   sco=NULL;
   cls=NULL;
}

template<class Y,class Z>
SvmGX<Y,Z>::SvmGX(const char *namspost,const char *pnam) : CMark<Y,Z>(namspost,pnam){
   wt=NULL;
   wd=NULL;
   sr=NULL;
   sco=NULL;
   cls=NULL;
}

template<class Y,class Z>
SvmGX<Y,Z>::~SvmGX(){
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   if(sco!=NULL)delete [] sco;
   if(cls!=NULL)delete [] cls;
}

template<class Y,class Z>
void SvmGX<Y,Z>::Set_Docs(Indx<Y> *gdd,Indx<Y> *bdd){
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
void SvmGX<Y,Z>::LearnSGD(Y iter,Y ui,CValdX<Y> *Cval){
   Y i,j,k,n,sc,iz,iu,ix;
   Y gt,ik;
   Z xx,yy,cv,eta=0.002;
   Z sum,zz,uu,rp,rpo;
   Z zm,vax,zmo;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd];
   sr=new Y[tdoc];

   th=0;
   for(i=0;i<this->nwrd;i++)wt[i]=0;

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         sr[n++]=i;
      }
   }
   rpo=tdoc;
   zmo=(Z)tdoc;
for(ik=0;ik<iter;ik++){
   vax=0;
   //Randomize order
   for(i=tdoc-1;i>0;i--){
      k=(Y)zrand((long)i+1);
      j=sr[i];
      sr[i]=sr[k];
      sr[k]=j;
   }
   //Deal with main loss function
   zm=0;
   rp=0;
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
         uu=-yy;
         zm+=1.0-xx;
         th-=eta*uu;
         n=1;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wt[iu]-=eta*uu;
               n++;
            }
         }
         rp+=eta*n;
      }
     this->mark(i,1000000,"docs");
   }
   cout << "Pseudo loss " << zm <<  " violators " << rp << " ratio " << (zmo-zm)/rpo << endl;
   rpo=rp;
   zmo=zm;
   vax=xx*xx;
   this->pflag=0;
   ScoreAll();
   DocAnal(ik);
   //training performance
   Indx<Y> *test=new Indx<Y>(Cval->pWTR[ui]);
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(test->ix,test,this->sco);
   EvalX<Y,Z> Elv(this->ndoc);
   Elv.load(Cval->pGTR[ui]);
   Z av=Elv.ave_prec(pOrd);
   delete pOrd;
   delete test;
   //testing performance
   test=new Indx<Y>(Cval->pWTS[ui]);
   pOrd=new Ordr<Y,Z>(test->ix,test,this->sco);
   EvalX<Y,Z> Ets(this->ndoc);
   Ets.load(Cval->pGTS[ui]);
   Z bv=Ets.ave_prec(pOrd);
   delete pOrd;
   delete test;
   cout << ik << "\t"<< bv << "\t" << av <<endl;
 }
   delete [] sr;
   sr=NULL;

}

template<class Y,class Z>
void SvmGX<Y,Z>::LearnSGDTest(Y iter,Y ui,CValdX<Y> *Cval){
   Y i,j,k,n,sc,iz,iu,rp,ix,rpo;
   Y gt,ik,xt,rn,rs;
   Z xx,yy,cv,eta=0.002;
   Z sum,zz,uu,lim=0.0;
   Z zm,vax,*vxx,xm,xv;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd];
   sr=new Y[tdoc];
   vxx=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)vxx[i]=1.0;

   th=0;
   for(i=0;i<this->nwrd;i++)wt[i]=0;

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         sr[n++]=i;
      }
   }
   rpo=tdoc;

xt=0;
rs=0;
for(ik=0;ik<iter;ik++){
   //Randomize order
   for(i=tdoc-1;i>0;i--){
      k=(Y)zrand((long)i+1);
      j=sr[i];
      sr[i]=sr[k];
      sr[k]=j;
   }
   //Deal with main loss function
   xm=xv=0;
   zm=0;
   rp=0;
   rn=0;
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
         uu=-yy;
         if(yy==1)rp++;
         else rn++;
         zm+=1.0-xx;
         th-=eta*uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wt[iu]-=eta*uu;
            }
         }
      }
      this->mark(i,1000000,"docs");
   }
   rs+=rn-rp;
   //xt+=rp;
   //cout << "Psuedo Loss " << zm << " violators " << rp << " total time " << xt << endl;
   cout << "Diff " << rs << endl;
   this->pflag=0;
   ScoreAll();
   //training performance
   Indx<Y> *test=new Indx<Y>(Cval->pWTR[ui]);
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(test->ix,test,this->sco);
   EvalX<Y,Z> Elv(this->ndoc);
   Elv.load(Cval->pGTR[ui]);
   Z av=Elv.ave_prec(pOrd);
   delete pOrd;
   delete test;
   //testing performance
   test=new Indx<Y>(Cval->pWTS[ui]);
   pOrd=new Ordr<Y,Z>(test->ix,test,this->sco);
   EvalX<Y,Z> Ets(this->ndoc);
   Ets.load(Cval->pGTS[ui]);
   Z bv=Ets.ave_prec(pOrd);
   delete pOrd;
   delete test;
   cout << ik << "\t"<< bv << "\t" << av <<endl;
 }
   delete [] sr;
   sr=NULL;
   delete [] vxx;

}

template<class Y,class Z>
void SvmGX<Y,Z>::LearnSGD_lw(Y iter,Y ui,CValdX<Y> *Cval){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y gt,ik;
   Z xx,yy,cv,eta=0.008;
   Z sum,zz,uu;
   Z zm,vax;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd];
   sr=new Y[tdoc];

   th=0;
   for(i=0;i<this->nwrd;i++)wt[i]=0;

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         sr[n++]=i;
      }
   }

for(ik=0;ik<iter;ik++){
   vax=0;
   //Randomize order
   for(i=tdoc-1;i>0;i--){
      k=(Y)zrand((long)i+1);
      j=sr[i];
      sr[i]=sr[k];
      sr[k]=j;
   }
   //Deal with main loss function
   zm=0;
   rp=0;
   for(i=0;i<tdoc;i++){
      j=sr[i];
      yy=cls[j];
      readp_db(j);
      readz_db(j);
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu]*this->lwt[iz];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         uu=-yy;
         zm+=1.0-xx;
         th-=eta*uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wt[iu]-=eta*uu*this->lwt[iz];
            }
         }
         rp++;
      }
     this->mark(i,1000000,"docs");
   }
   xx=Norm();
   cout << "Pseudo loss " << zm <<  " violators " << rp << " norm2 " << xx*xx << " ** " << zm+lambda*xx*xx << endl;
   vax=xx*xx;
   this->pflag=0;
   ScoreAll_lw();
   //training performance
   Indx<Y> *test=new Indx<Y>(Cval->pWTR[ui]);
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(test->ix,test,this->sco);
   EvalX<Y,Z> Elv(this->ndoc);
   Elv.load(Cval->pGTR[ui]);
   Z av=Elv.ave_prec(pOrd);
   delete pOrd;
   delete test;
   //testing performance
   test=new Indx<Y>(Cval->pWTS[ui]);
   pOrd=new Ordr<Y,Z>(test->ix,test,this->sco);
   EvalX<Y,Z> Ets(this->ndoc);
   Ets.load(Cval->pGTS[ui]);
   Z bv=Ets.ave_prec(pOrd);
   delete pOrd;
   delete test;
   cout << ik << "\t"<< bv << "\t" << av <<endl;
 }
   delete [] sr;
   sr=NULL;

}

template<class Y,class Z>
void SvmGX<Y,Z>::LearnSGDO(Y iter,Z neta,Z peta,CValdX<Y> *Cval){
   Y i,j,k,n,m,sc,iz,iu,rn,rp;
   Y gt,*sn,*sp,mx,ix,jx;
   Z xx,yy,cv;
   Z sum,zz,uu;
   Z zn,zp;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   wt=new Z[this->nwrd];
   sn=new Y[tdoc];
   sp=new Y[tdoc];

   th=0;
   for(i=0;i<this->nwrd;i++)wt[i]=0;

   n=m=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]<0)sn[n++]=i;
      else if(cls[i]>0)sp[m++]=i;
   }
   mx=(m<n)?n:m;

ix=jx=0;
for(gt=1;gt<=iter;gt++){
   zn=zp=0;
   rn=rp=0;
   for(i=0;i<mx;i++){
      //deal with neg class
      j=sn[ix++];
      yy=-1;
      readp_db(j);
      //Deal with main loss function
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         uu=-yy;
         zn+=1.0-xx;
         th-=neta*uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wt[iu]-=neta*uu;
            }
         }
         rn++;
      }
      if(!(ix%n)){
         xshuffle(n,sn);
         ix=0;
     }

      //deal with pos class
      j=sp[jx++];
      yy=1;
      readp_db(j);
      //Deal with main loss function
      sum=th;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wt[iu];
         }
      }
      //Loss zz calculation from p=sum
      xx=yy*sum;
      if(xx<1.0){
         uu=-yy;
         zp+=1.0-xx;
         th-=peta*uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wt[iu]-=peta*uu;
            }
         }
         rp++;
      }
      if(!(jx%m)){
         xshuffle(m,sp);
         jx=0;
     }
     this->mark(i,1000000,"docs");
   }
   cout << "Pseudo - loss " << zn <<  " - violators " << rn << endl;
   cout << "Pseudo + loss " << zp <<  " + violators " << rp << endl;
   this->pflag=0;
   ScoreAll();
   Indx<Y> *test=new Indx<Y>(Cval->pWTS[0]);
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(test->ix,test,this->sco);
   EvalX<Y,Z> Elv(this->ndoc);
   Elv.load(Cval->pGTS[0]);
   Z av=Elv.ave_prec(pOrd);
   delete pOrd;
   delete test;
   cout << "## "<<av <<endl;
 }
   delete [] sn;
   delete [] sp;

}

template<class Y,class Z>
void SvmGX<Y,Z>::LearnSGD(Y iter){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,mode,ct,gt,cu;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta=0.002;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
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

   for(gt=0;gt<iter;gt++){
      xshuffle(tdoc,sr);
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
            th+=yy*eta;
            for(iz=0;iz<this->nw;iz++){
               if(this->mrk[iu=this->nwd[iz]]){
                  wt[iu]+=yy*eta;
               }
            }
         }
      }
      if(this->pflag)cout << "iter " << gt+1 << endl;
    }
    delete [] sr;
    sr=NULL;
}

template<class Y,class Z>
void SvmGX<Y,Z>::LearnSGD_lw(Y iter){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,mode,ct,gt,cu;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta=0.008;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
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

   for(gt=0;gt<iter;gt++){
      xshuffle(tdoc,sr);
      for(i=0;i<tdoc;i++){
         j=sr[i];
         yy=cls[j];
         readp_db(j);
         readz_db(j);
         sum=th;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               sum+=wt[iu]*this->lwt[iz];
            }
         }
         //Loss zz calculation from p=sum
         xx=yy*sum;
         if(xx<1.0){
            th+=yy*eta;
            for(iz=0;iz<this->nw;iz++){
               if(this->mrk[iu=this->nwd[iz]]){
                  wt[iu]+=yy*eta*this->lwt[iz];
               }
            }
         }
      }
      if(this->pflag)cout << "iter " << gt+1 << endl;
    }
    delete [] sr;
    sr=NULL;
}


template<class Y,class Z>
void SvmGX<Y,Z>::DocAnal(Y nf){
  Y i,j,k;
  Z xx,yy,zz;

  this->gopen_db_map();
  ofstream *pfout=this->get_Ostr(nf,"dsc");

  for(i=0;i<this->ndoc;i++){
     if(cls[i]>0){
        this->readp_db(i);
        xx=yy=0;
        for(j=0;j<this->nw;j++){
           if(wt[this->nwd[j]]<0)xx-=wt[this->nwd[j]];
        }
        j=0;
        while((j<this->nw)&&(yy<=xx)){
           if(wt[this->nwd[j]]>0)yy+=wt[this->nwd[j]];
           j++;
        }
        *pfout << i << "\t" << 1 << "\t" << j << endl;
     }
     else if(cls[i]<0){
        this->readp_db(i);
        xx=yy=0;
        for(j=0;j<this->nw;j++){
           if(wt[this->nwd[j]]>0)xx-=wt[this->nwd[j]];
        }
        j=0;
        while((j<this->nw)&&(yy>=xx)){
           if(wt[this->nwd[j]]<0)yy+=wt[this->nwd[j]];
           j++;
        }
        *pfout << i << "\t" << -1 << "\t" << j << endl;
     }
     mark(this->pflag,i,1000,"A docs");
  }
  this->dst_Ostr(pfout);
}

template<class Y,class Z>
Z *SvmGX<Y,Z>::ScoreAll(void){
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
Z *SvmGX<Y,Z>::ScoreSet(Indx<Y> *ind){
   Y i,n;
   Z sum;

   if(sco!=NULL)delete [] sco;
   sco=new Z[ind->ix];

   for(i=0;i<ind->ix;i++){
      readp_db(ind->idx[i]);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]];
      sco[i]=sum;
      this->mark(i+1,100,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z *SvmGX<Y,Z>::ScoreAll_lw(void){
   Y j,n;
   Z sum;

   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];

   for(j=0;j<this->ndoc;j++){
      readp_db(j);
      readz_db(j);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]]*this->lwt[n];
      sco[j]=sum;
      this->mark(j+1,100,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z *SvmGX<Y,Z>::ScoreSet_lw(Indx<Y> *ind){
   Y i,j,n;
   Z sum;

   if(sco!=NULL)delete [] sco;
   sco=new Z[ind->ix];

   for(i=0;i<ind->ix;i++){
      j=ind->idx[i];
      readp_db(j);
      readz_db(j);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]]*this->lwt[n];
      sco[i]=sum;
      this->mark(i+1,100,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z SvmGX<Y,Z>::Norm(void){
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
Z SvmGX<Y,Z>::Func(void){
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
         readp_db(i);
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
   tf=xf/lambda+0.5*xnm;
   return(tf);
}

template<class Y,class Z>
Z SvmGX<Y,Z>::Func_lw(void){
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
         readp_db(i);
         readz_db(i);
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
   tf=xf/lambda+0.5*xnm;
   return(tf);
}

template<class Y,class Z>
void SvmGX<Y,Z>::Set_deriv(void){
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
      readp_db(j);
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
void SvmGX<Y,Z>::Save(int n){
   this->put_Nnum(n,"zh",this->ndoc,this->nwrd);
   this->bin_Writ(n,"weight",this->nwrd*sizeof(Z),(char*)wt);
   this->bin_Writ(n,"thresh",sizeof(Z),(char*)&th);
}

template<class Y,class Z>
void SvmGX<Y,Z>::Load(int n){
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
void SvmGX<Y,Z>::Release(int n){
    this->dst_Mmap(n,"weight",(char*&)wt);
}

}
#endif
