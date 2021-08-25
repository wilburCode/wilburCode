#ifndef SVMR_H
#define SVMR_H

#include <iostream>
#include <fstream>
#include <DataObj.h>
#include <XPost.h>
#include <Elev.h>
#include <CMark.h>

using namespace std;
namespace iret {

//template class SvmX
template<class Y,class Z>
class SvmX : public CMark<Y,Z> {
public:
   SvmX(const char *nspost);//name of XPost set
   SvmX(const char *nspost,const char *pnam);//name of XPost set
      //pnam is either used in "path_pnam" as place to find path
      //or if begins with ':' is followed by path itself.
   ~SvmX(void);

   //Learning functions
   void Set_Docs(Indx<Y> *gdd,Indx<Y> *bdd); //Set cls array
     //gdd good docs, bdd bad docs
   void Set_lambdaTest(Y rr); //Number of rounds rr (Zhang paper)
   void Set_Cnx(Z cx); //Sets parameter cnx
   void Set_Cnx_Norm(Z cx); //1.0E7 is good cx for Rebase
   void Set_Cnx_Norm_lw(Z cx); //10,000 is good cx for Rebase
   void Learn(void); //Optimizer
   void LearnProbeL(Y cy,Y lv); //Optimizer, cy the number of cycles per level
     //lv is number of levels. Does the line search
   void LearnSProbeL(Y cy,Y lv); //Optimizer, cy the number of cycles per level
     //lv is number of levels. Does the line search
     //Skip version of algorithm
   void LearnPlus(Z xu); //Optimizer, xu multiplies weight of + class
   void Learn_S(void); //Optimizer, skips 9/10 of time for already correct
     //is pointer at the partials array
   void Learn_S2(Z eps); //Optimizer, skips 9/10 of time for already correct
     //is pointer at the partials array
     //eps is the depth parameter
   void Learn_SPlus(Z xu); //Optimizer, skips 9/10 of time for already correct
     //is pointer at the partials array, xu multiplies weight of + class
   void Learn_lw(void); 
     //Uses local weights
   void LearnPlus_lw(Z xu); //Optimizer, xu multiplies weight of + class
   void Learn_S_lw(void); 
     //Uses local weights
   void LearnPos_lw(Z xu); 
     //Uses local weights and the weights it produces are all positive
     //Optimizer, xu multiplies weight of + class
   void Learn_SGrad(Y cyc); //Optimizer, skips 9/10 of time for already correct
     //is pointer at the partials array
     //attempt to use the gradient to steer the process
     //cyc is number of special increases between xfac decrements
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
SvmX<Y,Z>::SvmX(const char *namspost) : CMark<Y,Z>(namspost){
   wt=NULL;
   wd=NULL;
   sr=NULL;
   sco=NULL;
   cls=NULL;
}

template<class Y,class Z>
SvmX<Y,Z>::SvmX(const char *namspost,const char *pnam) : CMark<Y,Z>(namspost,pnam){
   wt=NULL;
   wd=NULL;
   sr=NULL;
   sco=NULL;
   cls=NULL;
}

template<class Y,class Z>
SvmX<Y,Z>::~SvmX(){
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   if(sco!=NULL)delete [] sco;
   if(cls!=NULL)delete [] cls;
}

template<class Y,class Z>
void SvmX<Y,Z>::Set_Docs(Indx<Y> *gdd,Indx<Y> *bdd){
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
void SvmX<Y,Z>::Set_lambdaTest(Y rr){
   Y i,j,k; 
   Z sum=0,xx,yy,zz,cnx;

   this->gopen_db_map();
   yy=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         yy++;
         this->readp_db(i);
         xx=0;
         for(j=0;j<this->nw;j++){
            if(this->mrk[this->nwd[j]]){
               xx++;
            }
         }
         sum+=xx;
      }
   }
   sum/=yy;
   cout << "Ssquared " << sum << " number T " << rnd(yy) << endl;
   lambda=2.0*tdoc/(yy*rr*0.002); //Note tdoc and yy cancel
   cin >> k;
}

template<class Y,class Z>
void SvmX<Y,Z>::Set_Cnx(Z cx){
   Y i;
   lambda=tdoc/cx;
}        

template<class Y,class Z>
void SvmX<Y,Z>::Set_Cnx_Norm(Z cx){
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
   lambda=tdoc*sum*sum/cx;
   cnx=cx/(sum*sum);
   if(this->pflag)cout << "cnx= " << cnx << endl;
}

template<class Y,class Z>
void SvmX<Y,Z>::Set_Cnx_Norm_lw(Z cx){
   Y i,j,k;
   Z sum=0,xx,yy,zz,cnx;

   this->gopen_db_map();
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         this->readp_db(i);
         readz_db(i);
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
   lambda=tdoc*sum*sum/cx;
   cnx=cx/(sum*sum);
   if(this->pflag)cout << "cnx= " << cnx << endl;
}

template<class Y,class Z>
void SvmX<Y,Z>::Learn(void){
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
         uu=-yy;
         f+=1.0-xx;
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
void SvmX<Y,Z>::LearnProbeL(Y cy,Y lv){
   Y i,j,k,m,n,iz,iu,rp,ix,fm;
   Y flag,mode,ct,gt,xcy,up,xlv=0;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta;
   Z dth,f,of,xrch,xvt,ux,uy;
   Z fac,minf,rch,rfac,glam,sglam;
   Z ss,dl,tho,*wto,llm;
   Z efu,efd,ominf,xfac,proj;
   Z gzro,gone,geas,gsum,edg,xgeas;
   Z *sc,*sd,*ea,*sm,*cd;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
   wd=new Z[this->nwrd];
   sr=new Y[tdoc];
   sc=new Z[tdoc];
   sd=new Z[tdoc];
   ea=new Z[tdoc];
   Store<Z> St(this->nwrd+1);
   if(tdoc<2000)fm=2;
   else fm=tdoc/1000;
   cd=new Z[fm];
   ofstream fout("bounce.log",ios::out|ios::app);
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
   dw=0.5*fm;
   proj=0;
   of=ominf=minf=1.0E100;
   fac=0.5;
   up=1;

   gt=ct=xcy=0;
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
      ea[i]=0;
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
      sc[i]=xx;
      if(xx<1.0){
         uu=-yy;
         f+=1.0-xx;
         dth+=uu;
         for(iz=0;iz<this->nw;iz++){
            if(this->mrk[iu=this->nwd[iz]]){
               wd[iu]+=uu;
            }
         }
      }
   }
   //Track the minimum f seen
   if(f>of){xcy++;up=1;}
   if(f<minf){minf=f;fout << "          CurMin " << minf << endl;wt[this->nwrd]=th;St.CopyIn(wt);}
   if(xcy>=cy){
      fout << "switch" << endl;
      xcy=0;
      xlv++;
      fac*=0.5;
      if(xlv>=lv){flag=0;up=1;}
   }
   of=f;
   cout << "    f " << f << " minf " << minf << endl << endl;

if(up){
   fout << "proj " << proj << endl;
   //regularizer change rate
   gzro=-th*dth;
   gone=dth*dth;
   ct=0;
   for(i=0;i<this->nwrd;i++){
      gzro-=wt[i]*wd[i];
      gone+=wd[i]*wd[i];
   }
   //wd scores
   geas=0;
   gsum=0;
   for(i=0;i<tdoc;i++){
      j=sr[i];
      yy=cls[j];
      this->readp_db(j);
      sum=dth;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wd[iu];
         }
      }
      sd[i]=xx=yy*sum;
      if(fabs(xx)>1.0E-10){
         wx=sc[i];
         zz=(wx-1.0)/xx;
         if(zz>0){
            ea[i]=zz;
            geas+=zz;
            ct++;
            if(wx<1.0){gsum+=xx;sd[i]=-sd[i];}
         }
         else if(zz==0){
            if(xx>0){
               ea[i]=1.0E-20;
            }
         }
         else if(xx>0)gsum+=xx;
      }
   }
   geas/=dw*ct;
   xgeas=(fm-1.0)*geas;
   for(i=0;i<fm;i++)cd[i]=0;
   for(i=0;i<tdoc;i++){
      xx=ea[i];
      if(xx>0){
         if(xx>xgeas)k=fm-1;
         else k=(Y)floor(xx/geas);
         cd[k]+=sd[i];
      }
   }
   gsum+=gzro;
   i=0;
cout << "gsum " << gsum << endl;
   while((gsum<0)&&(i<fm)){
      gsum+=gone*geas+cd[i++];
   }
   if((i==fm)&&(gsum<0)){
      eta=geas*fm;
cout << "type 1 eta " << eta << endl;
   }
   else {
      gsum-=gone*geas+cd[--i];
      j=0;
      edg=i*geas;
      for(k=0;k<tdoc;k++){
         xx=ea[k];
         if(xx>0){
            if(xx>xgeas)m=fm-1;
            else m=(Y)floor(xx/geas);
            if(m==i){
               ea[j]=xx;
               sd[j]=sd[k];
               j++;
            }
         }
      }
      hSort(j,ea,sd);
      k=0;
      xx=gsum;
      while((k<j)&&(xx<0)){
         gsum+=sd[k];
         xx=gsum+(ea[k++]-edg)*gone;
      }
      if((k==j)&&(xx<0)){
         eta=edg-gsum/gone;
cout << "type 2 j " << j << " eta " << eta << endl;
      }
      else {
         xx-=sd[--k];
         if(xx<0)eta=ea[k];
         else {gsum-=sd[k];eta=edg-gsum/gone;}
cout << "type 3 k j " << k << " " << j << " eta " << eta << endl;
      }
   }
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];
   up=0;
}
else {
   del=f-minf+fac*minf;
   proj=minf-fac*minf;

   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   eta=del/grad;
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];
   cout << " grad " << grad << " eta " << eta << endl;
}
   cout << gt << "  " << xcy << endl;
 }

   delete [] wd;
   wd=NULL;
   delete [] sr;
   sr=NULL;
   delete [] sd;
   delete [] sc;
   delete [] ea;
   delete [] cd;
   cout << "Minf " << minf << endl;
   fout << "          FinalMin " << minf << endl << endl << endl;
   fout.close();
   fout.clear();
   St.CopyOut(wt);
   th=wt[this->nwrd];
}

template<class Y,class Z>
void SvmX<Y,Z>::LearnSProbeL(Y cy,Y lv){
   Y i,j,k,m,n,iz,iu,rp,ix,fm;
   Y flag,mode,ct,gt,xcy,up,xlv=0;
   int *wk,wkk;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta;
   Z dth,f,of,xrch,xvt,ux,uy;
   Z fac,minf,rch,rfac,glam,sglam;
   Z ss,dl,tho,*wto,llm;
   Z efu,efd,ominf,xfac,proj;
   Z gzro,gone,geas,gsum,edg,xgeas;
   Z *sc,*sd,*ea,*sm,*cd;

   this->gopen_db_map();
   if(wt!=NULL)delete [] wt;
   if(wd!=NULL)delete [] wd;
   if(sr!=NULL)delete [] sr;
   wt=new Z[this->nwrd+1];
   wd=new Z[this->nwrd];
   sr=new Y[tdoc];
   sc=new Z[tdoc];
   sd=new Z[tdoc];
   ea=new Z[tdoc];
   wk=new int[tdoc];
   Store<Z> St(this->nwrd+1);
   if(tdoc<2000)fm=2;
   else fm=tdoc/1000;
   cd=new Z[fm];
   ofstream fout("bounce.log",ios::out|ios::app);
   th=0;
   for(i=0;i<this->nwrd;i++){
      wt[i]=0;
   }

   n=0;
   for(i=0;i<this->ndoc;i++){
      if(cls[i]){
         wk[n]=0;
         sr[n++]=i;
      }
   }
   dw=0.5*fm;
   proj=0;
   of=ominf=minf=1.0E100;
   fac=0.5;
   up=1;

   gt=ct=xcy=0;
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
      ea[i]=0;
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
      sc[i]=xx;
      if(xx<1.0){
         wk[i]=10;
         uu=-yy;
         f+=1.0-xx;
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
   if(f>of){xcy++;up=1;}
   if(f<minf){minf=f;fout << "          CurMin " << minf << endl;wt[this->nwrd]=th;St.CopyIn(wt);}
   if(xcy>=cy){
      fout << "switch" << endl;
      xcy=0;
      xlv++;
      fac*=0.5;
      if(xlv>=lv){flag=0;up=1;}
   }
   of=f;
   cout << "    f " << f << " minf " << minf << endl << endl;

if(up){
   for(i=0;i<tdoc;i++){
      if(wk[i]<=0){
      ea[i]=0;
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
      sc[i]=xx;
      if(xx<1.0)wk[i]=10;
      }
   }
   fout << "proj " << proj << endl;
   //regularizer change rate
   gzro=-th*dth;
   gone=dth*dth;
   ct=0;
   for(i=0;i<this->nwrd;i++){
      gzro-=wt[i]*wd[i];
      gone+=wd[i]*wd[i];
   }
   //wd scores
   geas=0;
   gsum=0;
   for(i=0;i<tdoc;i++){
      j=sr[i];
      yy=cls[j];
      this->readp_db(j);
      sum=dth;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=wd[iu];
         }
      }
      sd[i]=xx=yy*sum;
      if(fabs(xx)>1.0E-10){
         wx=sc[i];
         zz=(wx-1.0)/xx;
         if(zz>0){
            ea[i]=zz;
            geas+=zz;
            ct++;
            if(wx<1.0){gsum+=xx;sd[i]=-sd[i];}
         }
         else if(zz==0){
            if(xx>0){
               ea[i]=1.0E-20;
            }
         }
         else if(xx>0)gsum+=xx;
      }
   }
   geas/=dw*ct;
   xgeas=(fm-1.0)*geas;
   for(i=0;i<fm;i++)cd[i]=0;
   for(i=0;i<tdoc;i++){
      xx=ea[i];
      if(xx>0){
         if(xx>xgeas)k=fm-1;
         else k=(Y)floor(xx/geas);
         cd[k]+=sd[i];
      }
   }
   gsum+=gzro;
   i=0;
cout << "gsum " << gsum << endl;
   while((gsum<0)&&(i<fm)){
      gsum+=gone*geas+cd[i++];
   }
   if((i==fm)&&(gsum<0)){
      eta=geas*fm;
cout << "type 1 eta " << eta << endl;
   }
   else {
      gsum-=gone*geas+cd[--i];
      j=0;
      edg=i*geas;
      for(k=0;k<tdoc;k++){
         xx=ea[k];
         if(xx>0){
            if(xx>xgeas)m=fm-1;
            else m=(Y)floor(xx/geas);
            if(m==i){
               ea[j]=xx;
               sd[j]=sd[k];
               j++;
            }
         }
      }
      hSort(j,ea,sd);
      k=0;
      xx=gsum;
      while((k<j)&&(xx<0)){
         gsum+=sd[k];
         xx=gsum+(ea[k++]-edg)*gone;
      }
      if((k==j)&&(xx<0)){
         eta=edg-gsum/gone;
cout << "type 2 j " << j << " eta " << eta << endl;
      }
      else {
         xx-=sd[--k];
         if(xx<0)eta=ea[k];
         else {gsum-=sd[k];eta=edg-gsum/gone;}
cout << "type 3 k j " << k << " " << j << " eta " << eta << endl;
      }
   }
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];
   up=0;
}
else {
   del=f-minf+fac*minf;
   proj=minf-fac*minf;

   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   eta=del/grad;
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];
   cout << " grad " << grad << " eta " << eta << endl;
}
   cout << gt << "  " << xcy << endl;
 }

   delete [] wd;
   wd=NULL;
   delete [] sr;
   sr=NULL;
   delete [] sd;
   delete [] sc;
   delete [] ea;
   delete [] cd;
   cout << "Minf " << minf << endl;
   fout << "          FinalMin " << minf << endl << endl << endl;
   fout.close();
   fout.clear();
   St.CopyOut(wt);
   th=wt[this->nwrd];
}

template<class Y,class Z>
void SvmX<Y,Z>::LearnPlus(Z xu){
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
         if(yy>0){
            uu=-yy*xu;
            f+=(1.0-xx)*xu;
         }
         else {
            uu=-yy;
            f+=1.0-xx;
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
void SvmX<Y,Z>::Learn_S(void){
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
         uu=-yy;
         f+=1.0-xx;
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

   cout << endl << gt << "  " << cu << endl;
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
void SvmX<Y,Z>::Learn_S2(Z eps){
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
         uu=-yy;
         f+=1.0-xx;
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

   if(fac<eps)break;
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
void SvmX<Y,Z>::Learn_SPlus(Z xu){
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
         if(yy>0){
            uu=-yy*xu;
            f+=(1.0-xx)*xu;
         }
         else {
            uu=-yy;
            f+=1.0-xx;
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
                 llm-=xfac*(minf-llm);
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
   cout << "Minf " << minf << endl;
}

template<class Y,class Z>
void SvmX<Y,Z>::Learn_lw(void){
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
         f+=1.0-xx;
         uu=-yy;
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
void SvmX<Y,Z>::LearnPos_lw(Z xu){
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
   for(i=0;i<this->nwrd;i++)wt[i]=-5.0;

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
      xx=exp((double)wt[i]);
      wd[i]=lambda*xx*xx;
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
      sum=0;
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            sum+=exp((double)wt[iu])*this->lwt[iz];
         }
      }
      //Loss zz calculation from p=sum
      if(yy<0){
         if(sum>1.0){
            f+=(sum-1.0)*(sum-1.0);
            for(iz=0;iz<this->nw;iz++){
               if(this->mrk[iu=this->nwd[iz]]){
                  wd[iu]+=2.0*(sum-1.0)*exp((double)wt[iu])*this->lwt[iz];
               }
            }
         }
      }
      else {
         if(sum<3.0){
            f+=(3.0-sum)*(3.0-sum)*xu;
            for(iz=0;iz<this->nw;iz++){
               if(this->mrk[iu=this->nwd[iz]]){
                  wd[iu]+=2.0*(sum-3.0)*exp((double)wt[iu])*this->lwt[iz]*xu;
               }
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
   if(gt>1000)break;
   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   eta=del/grad;
   if((ux=th-eta*dth)<2.0)th=ux;
   for(i=0;i<this->nwrd;i++){
      if((ux=wt[i]-eta*wd[i])<2.0)wt[i]=ux;
   }

   if(this->pflag)cout << endl << gt << "  " << cu << " mode " << mode << endl;
   if(this->pflag)cout << "    f " << f << " minf " << minf << " fac " << fac << endl;
   if(this->pflag)cout << " grad " << grad << " eta " << eta << endl;
 }

   delete [] wd;
   wd=NULL;
   delete [] sr;
   sr=NULL;
   St.CopyOut(wt);
   for(i=0;i<this->nwrd;i++)wt[i]=exp(wt[i]);
   if(this->pflag)cout << "Minf " << minf << endl;
}

template<class Y,class Z>
void SvmX<Y,Z>::LearnPlus_lw(Z xu){
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
         if(yy>0){
            f+=(1.0-xx)*xu;
            uu=-yy*xu;
         }
         else {
            f+=1.0-xx;
            uu=-yy;
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
void SvmX<Y,Z>::Learn_S_lw(void){
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
         uu=-yy;
         f+=1.0-xx;
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
void SvmX<Y,Z>::Learn_SGrad(Y cyc){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,flag2,gt,ct,cu,mode;
   int *wk,wkk;
   Z xx,yy,cv,dv,dw,wx,blk;
   Z sum,zz,uu,del,grad,eta;
   Z dth,f,of,xrch,*wdo,ux,uy;
   Z fac,minf,rfac,llm,dtho;
   Z efu,efd,ominf,xfac,sgrado,sgrad;

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
   of=minf=1.0E100;
   xfac=0.5;

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
         uu=-yy;
         f+=1.0-xx;
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

   if(f<minf){
      minf=f;
      wt[this->nwrd]=th;
      St.CopyIn(wt);
   }
   else {
      if(f>of)cu++;
      of=f;
   }
   if(cu==cyc){xfac*=2.0/3.0;cu=0;}
   if(xfac<1E-2)break;
   grad=dth*dth;
   for(i=0;i<this->nwrd;i++)grad+=wd[i]*wd[i];
   del=f-minf+xfac*minf;
   eta=del/grad;
   th-=eta*dth;
   for(i=0;i<this->nwrd;i++)wt[i]-=eta*wd[i];
   cout << gt << "  " << cu << endl;
   cout << "    f " << f << " minf " << minf << endl;
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
Z *SvmX<Y,Z>::ScoreAll(void){
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
Z *SvmX<Y,Z>::ScoreSet(Indx<Y> *ind){
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
Z *SvmX<Y,Z>::ScoreAll_lw(void){
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
Z *SvmX<Y,Z>::ScoreSet_lw(Indx<Y> *ind){
   Y i,j,n;
   Z sum;

   if(sco!=NULL)delete [] sco;
   sco=new Z[ind->ix];

   for(i=0;i<ind->ix;i++){
      j=ind->idx[i];
      this->readp_db(j);
      readz_db(j);
      sum=th;
      for(n=0;n<this->nw;n++)sum+=wt[this->nwd[n]]*this->lwt[n];
      sco[i]=sum;
      this->mark(i+1,100,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z SvmX<Y,Z>::Norm(void){
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
Z SvmX<Y,Z>::Func(void){
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
   tf=xf/lambda+0.5*xnm;
   return(tf);
}

template<class Y,class Z>
Z SvmX<Y,Z>::Func_lw(void){
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
   tf=xf/lambda+0.5*xnm;
   return(tf);
}

template<class Y,class Z>
void SvmX<Y,Z>::Set_deriv(void){
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
void SvmX<Y,Z>::Save(int n){
   this->put_Nnum(n,"zh",this->ndoc,this->nwrd);
   this->bin_Writ(n,"weight",this->nwrd*sizeof(Z),(char*)wt);
   this->bin_Writ(n,"thresh",sizeof(Z),(char*)&th);
}

template<class Y,class Z>
void SvmX<Y,Z>::Load(int n){
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
void SvmX<Y,Z>::Release(int n){
    this->dst_Mmap(n,"weight",(char*&)wt);
}

}
#endif
