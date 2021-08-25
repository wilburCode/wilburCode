#ifndef BAYEX_H
#define BAYEX_H

#include <iostream>
#include <fstream>
#include <Elev.h>
#include <Dist.h>
#include <CMark.h>

namespace iret {

template<class Y,class Z>
class BayeX : public CMark<Y,Z> {
   public:
      BayeX(const char *nmspost); //nmspost is SPost object name
      BayeX(const char *nmxpost,const char *pnam); //nmxpost is XPost object name
         //pnam points at name used for find path of from path_*pnam file. If begins
         //with ":" then ":" is ignored and remainder is the path
      ~BayeX();
      void gopen_BayeX(void); //Opens for operating.
         //initializes all pdoc to 1.0

      //Weighting functions.
      void weightSall(void); //Weights assigned to all terms by MBM 
      void weightSmrk(void); //Weights assigned to all terms by MBM 
      void weightDmrk(Y docNum); //Weights assigned to all terms of docNum by MBM 
      void weightRmrk(void); //Weights assigned to all terms by MBM. Inverts s and t
      void weightPQmrk(void); //Weights assigned to all terms by MBM. Uses the p-q strategy
      void weightSidf(void); //Weights all terms by inverse document frequency
           //weighting. Assumes that tx has full counts in it.
      void weightSidfMrk(Y docNum); //idf weights assigned to all terms of docNum  
           //that are marked
      void weightSmut(void); //Weights all terms by reciprocal of document frequency
           //weighting. Assumes that tx has full counts in it. Based on mutual info.
      void weightBall(void); //Weights assigned to all terms by BLodds 
      inline double BLodds(Z nst,Z nt,Z ns,Z num){ //BLodds(Bernoulli log odds ratio)
         double axx,pt,rt,rs;

         axx=0;
         rs=(double)ns/((double)num);
         rt=(double)nt/((double)num);
      
         if(nst){
            pt=(double)nst/((double)num);
            axx+=nst*log(pt/(rs*rt));
         }
         if(nt-nst){
            pt=(double)(nt-nst)/((double)num);
            axx+=(nt-nst)*log(pt/((1.0-rs)*rt));
         }
         if(ns-nst){
            pt=(double)(ns-nst)/((double)num);
            axx+=(ns-nst)*log(pt/((1.0-rt)*rs));
         }
         if(num-ns-nt+nst){
            pt=(double)(num-ns-nt+nst)/((double)num);
            axx+=(num-ns-nt+nst)*log(pt/((1.0-rs)*(1.0-rt)));
         }
         return(axx);
      }          

      //Scoring function - scores only terms with 1 in mrk array
      Z *ScoreAll(void); 
      Z *ScoreSet(Indx<Y> *ind);
         //Both require gopen_db_map to be called first
      Z *ScoreAll(Indx<Y> *pTer); //Scores the list of terms, but only if 1 mrk array
         //Requires gopen_map to be called first
      Z ScorePair(Indx<Y> *pTer,Indx<Y> *pSer); //Scores overlapping terms with weg if mrk=1
      void ScorePairM(Z *scx,Indx<Y> *pTer,Indx<Y> *pSer); //Scores overlapping terms with weg 
         //if mrk>0. scx dim 6, B-1, G-2, M-3, R-4, T-5. requires mrk array specially set
      void DScorePairM(Z *scx,Indx<Y> *pTer,Indx<Y> *pSer); //Scores overlapping terms with weg 
         //if mrk>0. scx dim 6, B-1, G-2, M-3, R-4, T-5. requires mrk array specially set. For debugging
      void SetMrkSp(long *tcx); //Sets mark array for ScorePairM function and collects total counts
         //in tcx array. Must have called gopen_lexos and gopen_map.
      Ordr<Y,Z> *Skim(Y n);
      //Allows to save and load weights after learning
      void Save(int n); //Saves the set of weights weg. Marks with n
      void Load(int n); //Loads the set of weights weg marked with n
      void Release(int n); //Unmaps the weight vector mapped by load

      //Global data
      Z *weg; //For weights, nwrd size.
      Z *axx; //For alpha values, nwrd size.
      Z *csg; //For individual cs numbers.
      Z cs; //Additve constant for correct scoring of docs.
      Z *sco;
};

template<class Y,class Z>
BayeX<Y,Z>::BayeX(const char *namspost) : CMark<Y,Z>(namspost){
   weg=NULL;
   axx=NULL;
   csg=NULL;
   sco=NULL;
}

template<class Y,class Z>
BayeX<Y,Z>::BayeX(const char *namspost,const char *pnam) : CMark<Y,Z>(namspost,pnam){
   weg=NULL;
   axx=NULL;
   csg=NULL;
   sco=NULL;
}

template<class Y,class Z>
BayeX<Y,Z>::~BayeX(){
   if(weg!=NULL)delete [] weg;
   if(axx!=NULL)delete [] axx;
   if(csg!=NULL)delete [] csg;
   if(sco!=NULL)delete [] sco;
}

template<class Y,class Z>
void BayeX<Y,Z>::gopen_BayeX(void){
   Y i;

   this->gopen_db_map();
   this->init_cnt();
   for(i=0;i<this->ndoc;i++)this->pdoc[i]=1.0;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightSall(void){
   Y nstw=0,flag,i;
   Z min,max,diff,n_t,n_st;
   Z xx,frc,pt,qt,rt;

   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   if(axx!=NULL)delete [] axx;
   axx=new Z[this->nwrd];
   cs=0;
   frc=((Z)this->nsx)/((Z)this->nnx);
   diff=this->nsx-this->nnx;

   for(i=0;i<this->nwrd;i++){
      n_t=*(this->tx+i);
      if(n_t&&(n_t<this->nnx)){
         min=(n_t<this->nsx)?n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(this->sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st -=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         nstw++;
         pt =(Z)n_st/(Z)this->nsx;
         qt =(Z)(n_t-n_st)/(Z)(this->nnx-this->nsx);
         rt =(Z)n_t/(Z)this->nnx;
         //calculate wt for this term
         *(weg+i)=log(pt*(1.0-qt))-log(qt*(1.0-pt));

         xx=-(Z)n_st*log(rt/pt) - (Z)(n_t-n_st)*log(rt/qt);
         xx-=(Z)(this->nsx-n_st)*log((1.0-rt)/(1.0-pt));
         xx-=(Z)(this->nnx-this->nsx-n_t+n_st)*log((1.0-rt)/(1.0-qt));
         *(axx+i)=xx;

         xx=log((Z)(this->nsx-n_st))+log((Z)(this->nnx-this->nsx));
         xx-=(log((Z)this->nsx)+log((Z)(this->nnx-this->nsx-n_t+n_st)));
         *(csg+i)=xx;
         cs+=xx;
      }
      else {
         *(weg+i)=0;
         *(csg+i)=0;
         *(axx+i)=0;
      }
      this->mark(i,10000,"terms weighted");
   }
if(this->pflag)cout << "cs= " << cs << endl;
if(this->pflag)cout << "Number of weighted terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightBall(void){
   Y nstw=0,flag,i;
   Z min,max,diff,n_t,n_st;
   Z xx,frc,pt,qt,rt;

   if(axx!=NULL)delete [] axx;
   axx=new Z[this->nwrd];
   cs=0;
   frc=((Z)this->nsx)/((Z)this->nnx);
   diff=this->nsx-this->nnx;

   for(i=0;i<this->nwrd;i++){
      n_t=*(this->tx+i);
      if(n_t&&(n_t<this->nnx)){
         min=(n_t<this->nsx)?n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(this->sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st -=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         nstw++;
         *(axx+i)=BLodds(n_st,n_t,this->nsx,this->nnx);
      }
      else {
         *(axx+i)=0;
      }
      this->mark(i,10000,"terms weighted");
   }
if(this->pflag)cout << "Number of weighted terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightSmrk(void){
   Y nstw=0,flag,i;
   Z min,max,diff,n_t,n_st;
   Z xx,frc,pt,qt,rt;

   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   if(axx!=NULL)delete [] axx;
   axx=new Z[this->nwrd];
   cs=0;
   frc=((Z)this->nsx)/((Z)this->nnx);
   diff=this->nsx-this->nnx;

   for(i=0;i<this->nwrd;i++){
      if((!this->mrk[i])||(this->sx[i]<2)){weg[i]=csg[i]=axx[i]=0;continue;}
      n_t=*(this->tx+i);
      if(n_t&&(n_t<this->nnx)){
         min=(n_t<this->nsx)?n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(this->sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st -=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         nstw++;
         pt =(Z)n_st/(Z)this->nsx;
         qt =(Z)(n_t-n_st)/(Z)(this->nnx-this->nsx);
         rt =(Z)n_t/(Z)this->nnx;
         //calculate wt for this term
         *(weg+i)=log(pt*(1.0-qt))-log(qt*(1.0-pt));

         xx=-(Z)n_st*log(rt/pt) - (Z)(n_t-n_st)*log(rt/qt);
         xx-=(Z)(this->nsx-n_st)*log((1.0-rt)/(1.0-pt));
         xx-=(Z)(this->nnx-this->nsx-n_t+n_st)*log((1.0-rt)/(1.0-qt));
         *(axx+i)=xx;

         xx=log((Z)(this->nsx-n_st))+log((Z)(this->nnx-this->nsx));
         xx-=(log((Z)this->nsx)+log((Z)(this->nnx-this->nsx-n_t+n_st)));
         *(csg+i)=xx;
         cs+=xx;
      }
      else {
         *(weg+i)=0;
         *(csg+i)=0;
         *(axx+i)=0;
      }
      this->mark(i,10000,"terms weighted");
   }
if(this->pflag)cout << "cs= " << cs << endl;
if(this->pflag)cout << "Number of weighted terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightDmrk(Y docNum){
   Y nstw=0,flag,i,j;
   Z min,max,diff,n_t,n_st;
   Z xx,frc,pt,qt,rt;

   this->readp_db(docNum);
   frc=((Z)this->nsx)/((Z)this->nnx);
   diff=this->nsx-this->nnx;

   for(j=0;j<this->nw;j++){
      i=this->nwd[j];
      if((!this->mrk[i])||(this->sx[i]<2)){weg[i]=csg[i]=axx[i]=0;continue;}
      n_t=*(this->tx+i);
      if(n_t&&(n_t<this->nnx)){
         min=(n_t<this->nsx)?n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(this->sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st -=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         nstw++;
         pt =(Z)n_st/(Z)this->nsx;
         qt =(Z)(n_t-n_st)/(Z)(this->nnx-this->nsx);
         rt =(Z)n_t/(Z)this->nnx;
         //calculate wt for this term
         *(weg+i)=log(pt*(1.0-qt))-log(qt*(1.0-pt));

         xx=-(Z)n_st*log(rt/pt) - (Z)(n_t-n_st)*log(rt/qt);
         xx-=(Z)(this->nsx-n_st)*log((1.0-rt)/(1.0-pt));
         xx-=(Z)(this->nnx-this->nsx-n_t+n_st)*log((1.0-rt)/(1.0-qt));
         *(axx+i)=xx;

         cs-=*(csg+i);
         xx=log((Z)(this->nsx-n_st))+log((Z)(this->nnx-this->nsx));
         xx-=(log((Z)this->nsx)+log((Z)(this->nnx-this->nsx-n_t+n_st)));
         *(csg+i)=xx;
         cs+=xx;
      }
      else {
         *(weg+i)=0;
         *(csg+i)=0;
         *(axx+i)=0;
      }
      this->mark(i,10000,"terms weighted");
   }
if(this->pflag)cout << "cs= " << cs << endl;
if(this->pflag)cout << "Number of weighted terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightRmrk(void){
   Y nstw=0,flag,i;
   Z min,max,diff,n_t,n_st;
   Z xx,frc,pt,qt,rt;
   Binomial Bn(1.0E-8,1.0E-8,2000000);

   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   if(axx!=NULL)delete [] axx;
   axx=new Z[this->nwrd];
   cs=0;
   frc=((Z)this->nsx)/((Z)this->nnx);
   diff=this->nsx-this->nnx;

   for(i=0;i<this->nwrd;i++){
      if((!this->mrk[i])||(this->sx[i]<2)){weg[i]=csg[i]=axx[i]=0;continue;}
      n_t=*(this->tx+i);
      if(*(this->sx+i)==n_t){weg[i]=csg[i]=axx[i]=0;continue;}
      if(n_t&&(n_t<this->nnx)){
         min=(n_t<this->nsx)?n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(this->sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st -=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         nstw++;
         pt =(Z)n_st/(Z)this->nsx;
         qt =(Z)(n_t-n_st)/(Z)(this->nnx-this->nsx);
         rt =(Z)n_t/(Z)this->nnx;
         //calculate wt for this term
         *(weg+i)=log(pt*(1.0-qt))-log(qt*(1.0-pt));
         
         //xx=log((Z)n_st)+log((Z)this->nnx)-log((Z)n_t)-log((Z)(this->nsx));
         if(n_t<2000000){
            *(axx+i)=Bn.lower_limit(n_t,n_st,0.01);
         }
         else *(axx+i)=0;

         xx=log((Z)(this->nsx-n_st))+log((Z)(this->nnx-this->nsx));
         xx-=(log((Z)this->nsx)+log((Z)(this->nnx-this->nsx-n_t+n_st)));
         *(csg+i)=xx;
         cs+=xx;
      }
      else {
         *(weg+i)=0;
         *(csg+i)=0;
         *(axx+i)=0;
      }
      this->mark(i,10000,"terms weighted");
   }
if(this->pflag)cout << "cs= " << cs << endl;
if(this->pflag)cout << "Number of weighted terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightPQmrk(void){
   Y nstw=0,flag,i;
   Z min,max,diff,n_t,n_st;
   Z xx,frc,pt,qt,rt;

   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   if(axx!=NULL)delete [] axx;
   axx=new Z[this->nwrd];
   cs=0;
   frc=((Z)this->nsx)/((Z)this->nnx);
   diff=this->nsx-this->nnx;

   for(i=0;i<this->nwrd;i++){
      if((!this->mrk[i])||(this->sx[i]<2)){weg[i]=csg[i]=axx[i]=0;continue;}
      n_t=*(this->tx+i);
      if(n_t&&(n_t<this->nnx)){
         min=(n_t<this->nsx)?n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(this->sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st -=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         nstw++;
         pt =(Z)n_st/(Z)n_t;
         //pt =(Z)n_st/(Z)this->nsx;
         //qt =(Z)(n_t-n_st)/(Z)(this->nnx-this->nsx);
         qt =(Z)(this->nsx-n_st)/(Z)(this->nnx-n_t);
         //rt =(Z)n_t/(Z)this->nnx;
         //calculate wt for this term
         *(weg+i)=pt-qt;

        /*
         xx=-(Z)n_st*log(rt/pt) - (Z)(n_t-n_st)*log(rt/qt);
         xx-=(Z)(this->nsx-n_st)*log((1.0-rt)/(1.0-pt));
         xx-=(Z)(this->nnx-this->nsx-n_t+n_st)*log((1.0-rt)/(1.0-qt));
         *(axx+i)=xx;

         xx=log((Z)(this->nsx-n_st))+log((Z)(this->nnx-this->nsx));
         xx-=(log((Z)this->nsx)+log((Z)(this->nnx-this->nsx-n_t+n_st)));
         */
         xx=0;
         *(csg+i)=xx;
         cs+=xx;
      }
      else {
         *(weg+i)=0;
         *(csg+i)=0;
         *(axx+i)=0;
      }
      this->mark(i,10000,"terms weighted");
   }
if(this->pflag)cout << "cs= " << cs << endl;
if(this->pflag)cout << "Number of weighted terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightSidf(void){
   Y i,j;
   Z xx=(Z)this->nnx,yy;
   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if(yy=*(this->tx+i)){
         if(yy<20)yy=20.0;
         *(weg+i)=log(xx/yy);
      }
      else {
         *(weg+i)=0;
         this->mrk[i]=0;
      }
      csg[i]=0;
   }
   cs=0;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightSidfMrk(Y docNum){
   Y i,j;
   Z xx=(Z)this->ndoc,yy;
   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   for(i=0;i<this->nwrd;i++)weg[i]=csg[i]=0;
   this->readp_db(docNum);
   for(i=0;i<this->nw;i++){
      j=this->nwd[i];
      if(this->mrk[j]&&(yy=*(this->freq+j))){
         if(yy<20)yy=20.0;
         *(weg+i)=log(xx/yy);
      }
      else {
         *(weg+i)=0;
      }
      csg[j]=0;
   }
   cs=0;
}

template<class Y,class Z>
void BayeX<Y,Z>::weightSmut(void){
   Y i,j;
   Z xx=(Z)this->nnx,yy;
   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if(yy=*(this->tx+i)){
         *(weg+i)=1.0/yy;
      }
      else {
         *(weg+i)=0;
         this->mrk[i]=0;
      }
      csg[i]=0;
   }
   cs=0;
}

template<class Y,class Z>
Z *BayeX<Y,Z>::ScoreAll(void){
   Y i,j,n;
   Z xx;

   cs=0;
   for(i=0;i<this->nwrd;i++){
      if(this->mrk[i])cs+=csg[i];
   }
   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];
   for(n=0;n<this->ndoc;n++){
      xx=cs;
      this->readp_db(n);
      for(i=0;i<this->nw;i++){
         j=this->nwd[i];
         if(this->mrk[j])xx+=weg[j];
      }
      sco[n]=xx;
      this->mark(n+1,1000,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z *BayeX<Y,Z>::ScoreAll(Indx<Y> *pTer){
   Y i,j,k,n;
   Z xx;
   Indx<Y> *pInd;

   cs=0;
   for(i=0;i<pTer->ix;i++){
      j=pTer->idx[i];
      if(this->mrk[j])cs+=csg[j];
   }
   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=cs;
   for(i=0;i<pTer->ix;i++){
      n=pTer->idx[i];
      if(this->mrk[n]){
         pInd=this->readp(n);
         for(j=0;j<pInd->ix;j++){
            k=pInd->idx[j];
            sco[k]+=weg[n];
         }
      }
      this->mark(i,1,"terms scored");
   }
   return(sco);
}

template<class Y,class Z>
Z *BayeX<Y,Z>::ScoreSet(Indx<Y> *ind){
   Y i,j,n,u;
   Z xx;

   cs=0;
   for(i=0;i<this->nwrd;i++){
      if(this->mrk[i])cs+=csg[i];
   }
   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];
   for(u=0;u<ind->ix;u++){
      n=ind->idx[u];
      xx=cs;
      this->readp_db(n);
      for(i=0;i<this->nw;i++){
         j=this->nwd[i];
         if(this->mrk[j])xx+=weg[j];
      }
      sco[n]=xx;
      this->mark(u+1,1000,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z BayeX<Y,Z>::ScorePair(Indx<Y> *pTer,Indx<Y> *pSer){
   Y i=0,j=0,k,m,n;
   Z xx=0;

   n=pSer->idx[j];
   while(i<pTer->ix){
      m=pTer->idx[i];
      if(m<n)i++;
      else if(n<m){
         while(j+1<pSer->ix){
            j++;
            if((n=pSer->idx[j])>=m)break;
         }
         if(n==m){
            if(this->mrk[n])xx+=weg[n];
  //if(this->mrk[n]){cout << "scored: " << this->show(n) << endl;}
            i++;j++;
            if(j<pSer->ix)n=pSer->idx[j];
         }
         else if(m<n)i++;
         else i=pTer->ix;
      }
      else {
         if(this->mrk[n])xx+=weg[n];
  //if(this->mrk[n]){cout << "scored: " << this->show(n) << endl;}
         i++;j++;
         if(j<pSer->ix)n=pSer->idx[j];
      }
   }
   return(xx);
}

template<class Y,class Z>
void BayeX<Y,Z>::ScorePairM(Z *scx,Indx<Y> *pTer,Indx<Y> *pSer){
   Y i=0,j=0,k,m,n;
   scx[1]=scx[2]=scx[3]=scx[4]=scx[5]=0;

   n=pSer->idx[j];
   while(i<pTer->ix){
      m=pTer->idx[i];
      if(m<n)i++;
      else if(n<m){
         while(j+1<pSer->ix){
            j++;
            if((n=pSer->idx[j])>=m)break;
         }
         if(n==m){
            if(this->mrk[n])scx[this->mrk[n]]+=weg[n];
            i++;j++;
            if(j<pSer->ix)n=pSer->idx[j];
         }
         else if(m<n)i++;
         else i=pTer->ix;
      }
      else {
         if(this->mrk[n])scx[this->mrk[n]]+=weg[n];
         i++;j++;
         if(j<pSer->ix)n=pSer->idx[j];
      }
   }
}

template<class Y,class Z>
void BayeX<Y,Z>::DScorePairM(Z *scx,Indx<Y> *pTer,Indx<Y> *pSer){
   Y i=0,j=0,k,m,n;
   scx[1]=scx[2]=scx[3]=scx[4]=scx[5]=0;

   n=pSer->idx[j];
   while(i<pTer->ix){
      m=pTer->idx[i];
      if(m<n)i++;
      else if(n<m){
         while(j+1<pSer->ix){
            j++;
            if((n=pSer->idx[j])>=m)break;
         }
         if(n==m){
            if(this->mrk[n])scx[this->mrk[n]]+=weg[n];
            if(this->mrk[n]){cout << "scored: " << this->show(n) << endl;}
            i++;j++;
            if(j<pSer->ix)n=pSer->idx[j];
         }
         else if(m<n)i++;
         else i=pTer->ix;
      }
      else {
         if(this->mrk[n])scx[this->mrk[n]]+=weg[n];
         if(this->mrk[n]){cout << "scored: " << this->show(n) << endl;}
         i++;j++;
         if(j<pSer->ix)n=pSer->idx[j];
      }
   }
}

template<class Y,class Z>
void BayeX<Y,Z>::SetMrkSp(long *tcx){
   Y i,j,k;
   char *pch;
   tcx[1]=tcx[2]=tcx[3]=tcx[4]=tcx[5]=0;
   if(this->mrk==NULL)this->mrk=new Y[this->nwrd];

   for(i=0;i<this->nwrd;i++){
      k=this->freq[i];
      pch=this->show(i);
      j=strlen(pch);
      if((j<4)||(pch[j-2]!='!')){
         this->mrk[i]=1;
         tcx[1]+=k;
      }
      else {
         switch(pch[j-1]){
            case 'j': this->mrk[i]=1;
                      tcx[1]+=k;
                      break;
            case 'f': this->mrk[i]=1;
                      tcx[1]+=k;
                      break;
            case 'g': this->mrk[i]=2;
                      tcx[2]+=k;
                      break;
            case 'm': this->mrk[i]=3;
                      tcx[3]+=k;
                      break;
            case 'r': this->mrk[i]=4;
                      tcx[4]+=k;
                      break;
            case 't': this->mrk[i]=5;
                      tcx[5]+=k;
                      break;
            default : this->mrk[i]=0;
         }
      }
      this->mark(i,100000,"terms");
   }
}
        
template<class Y,class Z>
void BayeX<Y,Z>::Save(int n){
   this->put_Nnum(n,"zb",this->ndoc,this->nwrd);
   this->bin_Writ(n,"weightb",this->nwrd*sizeof(Z),(char*)weg);
   this->bin_Writ(n,"threshb",sizeof(Z),(char*)&cs);
}

template<class Y,class Z>
void BayeX<Y,Z>::Load(int n){
   long i,j,k;
   this->get_Nnum(n,"zb",i,k);
   this->ndoc=i;
   this->nwrd=k;
   weg=(Z*)this->get_Mmap(n,"weightb");
   ifstream *pfin=this->get_Istr(n,"threshb",ios::in);
   pfin->read((char*)&cs,sizeof(Z));
   this->dst_Istr(pfin);
}

template<class Y,class Z>
void BayeX<Y,Z>::Release(int n){
    this->dst_Mmap(n,"weightb",(char*&)weg);
}

template<class Y,class Z> 
Ordr<Y,Z> *BayeX<Y,Z>::Skim(Y n){
   if(sco==NULL)return(NULL);
   Z sx;
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(n,this->ndoc,sco);

   Y i=0;
   pOrd->ind(i,sx);
   if(sx>0)return(pOrd);
   else return(NULL);
}

//Homogeneous Bayes class where Z is the float type throughout

template<class Y,class Z>
class BayeQ : public XPost<Y,Z> {
   public:
      BayeQ(const char *nmspost); //nmspost is SPost object name
         //Creates space for tx, sx, and ax.
      ~BayeQ();
      void gopen_BayeQ(void); //Opens for operating.
         //initializes all pdoc to 1.0
      void gclose_BayeQ(void); //Closes for operating.

      //Counting functions.
      void zerot(void); //Zeroes tx array.
      void freq_tx(void); //Fills tx array from freq data.
      void zeros(void); //Zeroes sx array.
      void countDoc(Y i); //Adds the counts from doc i to
         //tx using read function.
      void counsDoc(Y i); //Adds the counts from doc i to
         //sx using read function.
      void counbDoc(Y i); //Adds the counts from doc i to
         //sx & tx using read function.
      //Set reads
      void countTX(Indx<Y> *cnd); //adds all counts in cnd to tx.
      void countSX(Indx<Y> *cnd); //adds all counts in cnd to sx.
      void countBX(Indx<Y> *cnd); //adds all counts in cnd to sx & tx.

      //Weighting functions.
      void weightSall(void); //Weights assigned to all terms by MBM 
      void weightSidf(void); //Weights all terms by inverse document frequency
           //weighting. Assumes that tx has full counts in it.
      void Set_axx_chi(void); //Computes the chi square value and puts
           //in axx array.

      //Marking functions
      void Set_mrk(Indx<Y> *pTrm);//Set the mrk array to process terms 
         // listed in pTrm
      void Set_mrk_freq(Y nm); //Set the mrk arry for terms that
         //satisfy freq > nm
      void Set_mrk_weg(Z cut); //Set the mrk arry for terms that
         //satisfy |weg| > cut.
      void Set_mrk_axx(Z cut); //Set the mrk arry for terms that
         //satisfy axx > cut.
      void Cut_mrk_weg(Z cut); //Trim the mrk arry for terms that
         //satisfy |weg| <= cut.
      //Special mrk array functions
      void Set_All(int n); //Creates and sets mrk array to be n everywhere
      void Set_Char(char c,int n); //Sets mrk n if character c in string
      void Set_NChar(char c,int n); //Sets mrk n if character c not in string
      void Set_String(char *str,int n); //Sets mrk n if str substring of string
      void Set_NString(char *str,int n); //Sets mrk n if str not substring of string

      //Scoring function - scores only terms with 1 in mrk array
      Z *ScoreAll(void); 
      Z *ScoreSet(Indx<Y> *ind);

      //Global data
      Z *tx; //For total counts based on current set. nwrd size.
      Z *sx; //For total counts based on subject set. nwrd size.
      Z nnx; //Total documents for purposes of weighting.
      Z nsx; //Number of documents in the subject area for weighting.
      Z *pdoc; //Weight of documents
      Z *weg; //For weights, nwrd size.
      Z *axx; //For alpha values, nwrd size.
      Z *csg; //For individual cs numbers.
      Z cs; //Additve constant for correct scoring of docs.
      Z *sco;
      Y *mrk; //Marks which terms to include for scoring
};

template<class Y,class Z>
BayeQ<Y,Z>::BayeQ(const char *namspost) : XPost<Y,Z>(namspost){
   weg=NULL;
   axx=NULL;
   csg=NULL;
   mrk=NULL;
   sco=NULL;
   tx=NULL;
   sx=NULL;
   pdoc=NULL;
}

template<class Y,class Z>
BayeQ<Y,Z>::~BayeQ(){
   if(weg!=NULL)delete [] weg;
   if(axx!=NULL)delete [] axx;
   if(csg!=NULL)delete [] csg;
   if(mrk!=NULL)delete [] mrk;
   if(sco!=NULL)delete [] sco;
   if(tx!=NULL)delete [] tx;
   if(sx!=NULL)delete [] sx;
   if(pdoc!=NULL)delete [] pdoc;
}

template<class Y,class Z>
void BayeQ<Y,Z>::gopen_BayeQ(void){
   Y i;

   this->gopen_db_map();

   if(tx!=NULL)delete [] tx;
   tx=new Z[this->nwrd];

   if(sx!=NULL)delete [] sx;
   sx=new Z[this->nwrd];

   if(mrk!=NULL)delete [] mrk;
   mrk=new Y[this->nwrd];

   if(pdoc!=NULL)delete [] pdoc;
   pdoc=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)pdoc[i]=1.0;
}

template<class Y,class Z>
void BayeQ<Y,Z>::gclose_BayeQ(void){
   this->gclose_db_map();
   if(tx!=NULL){delete [] tx;tx=NULL;}
   if(sx!=NULL){delete [] sx;sx=NULL;}
   if(mrk!=NULL){delete [] mrk;mrk=NULL;}
   if(pdoc!=NULL){delete [] pdoc;pdoc=NULL;}
}

template<class Y,class Z>
void BayeQ<Y,Z>::zerot(void){
   Y i;
   nnx=0.0;
   for(i=0;i<this->nwrd;i++)*(tx+i)=0;
}

template<class Y,class Z>
void BayeQ<Y,Z>::freq_tx(void){
   Y i;
   nnx=(Z)this->ndoc;
   this->gopen_map();
   for(i=0;i<this->nwrd;i++)*(tx+i)=this->freq[i];
   this->gclose_map();
}

template<class Y,class Z>
void BayeQ<Y,Z>::zeros(void){
   Y i;
   nsx=0.0;
   for(Y i=0;i<this->nwrd;i++)*(sx+i)=0;
}

template<class Y,class Z>
void BayeQ<Y,Z>::countDoc(Y i){
   Y j,k;
   Z pt;
   pt=pdoc[i];

   readp_db(i);
   for(k=0;k<this->nw;k++){
      j=*(this->nwd+k);
      (*(tx+j))+=pt;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::counsDoc(Y i){
   Y j,k;
   Z pt;
   pt=pdoc[i];

   readp_db(i);
   for(k=0;k<this->nw;k++){
      j=*(this->nwd+k);
      (*(sx+j))+=pt;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::counbDoc(Y i){
   Y j,k;
   Z pt;
   pt=pdoc[i];

   readp_db(i);
   for(k=0;k<this->nw;k++){
      j=*(this->nwd+k);
      (*(sx+j))+=pt;
      (*(tx+j))+=pt;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::countTX(Indx<Y> *cnd){
   Y i,j;

   for(i=0;i<cnd->ix;i++){
      j=cnd->idx[i];
      this->countDoc(j);
      nnx+=pdoc[j];
      this->mark(i,1000,"Docs T counted");
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::countSX(Indx<Y> *cnd){
   long i,j;

   for(i=0;i<cnd->ix;i++){
      j=cnd->idx[i];
      this->counsDoc(j);
      nsx+=pdoc[j];
      this->mark(i,1000,"Docs S counted");
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::countBX(Indx<Y> *cnd){
   Y i,j;

   for(i=0;i<cnd->ix;i++){
      j=cnd->idx[i];
      this->counbDoc(j);
      nnx+=pdoc[j];
      nsx+=pdoc[j];
      this->mark(i,1000,"Docs B counted");
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::weightSall(void){
   Y nstw=0,flag,i;
   Z min,max,diff,n_t,n_st;
   Z xx,frc,pt,qt,rt;

   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   if(csg!=NULL)delete [] csg;
   csg=new Z[this->nwrd];
   if(axx!=NULL)delete [] axx;
   axx=new Z[this->nwrd];
   cs=0;
   frc=((Z)nsx)/((Z)nnx);
   diff=nsx-nnx;

   for(i=0;i<this->nwrd;i++){
      if((!mrk[i])||(sx[i]<2)){weg[i]=csg[i]=axx[i]=0;continue;}
      //if(!mrk[i]){weg[i]=0;continue;}
      n_t=*(tx+i);
      if(n_t&&(n_t<nnx)){
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st -=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         nstw++;
         pt =(Z)n_st/(Z)nsx;
         qt =(Z)(n_t-n_st)/(Z)(nnx-nsx);
         rt =(Z)n_t/(Z)nnx;
         //calculate wt for this term
         *(weg+i)=log(pt*(1.0-qt))-log(qt*(1.0-pt));

         xx=-(Z)n_st*log(rt/pt) - (Z)(n_t-n_st)*log(rt/qt);
         xx-=(Z)(nsx-n_st)*log((1.0-rt)/(1.0-pt));
         xx-=(Z)(nnx-nsx-n_t+n_st)*log((1.0-rt)/(1.0-qt));
         *(axx+i)=xx;

         xx=log((Z)(nsx-n_st))+log((Z)(nnx-nsx));
         xx-=(log((Z)nsx)+log((Z)(nnx-nsx-n_t+n_st)));
         *(csg+i)=xx;
         //if(fabs(weg[i])>6.0)cs+=xx;
         cs+=xx;
      }
      else {
         *(weg+i)=0;
         *(csg+i)=0;
         *(axx+i)=0;
         //mrk[i]=0;
      }
      this->mark(i,10000,"terms weighted");
   }
if(this->pflag)cout << "cs= " << cs << endl;
if(this->pflag)cout << "Number of weighted terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeQ<Y,Z>::weightSidf(void){
   Y i,j;
   Z xx=(Z)nnx,yy;
   if(weg!=NULL)delete [] weg;
   weg=new Z[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if(!mrk[i])continue;
      if(yy=*(tx+i)){
         if(yy<20)yy=20.0;
         *(weg+i)=log(xx/yy);
      }
      else {
         *(weg+i)=0;
         mrk[i]=0;
      }
   }
   cs=0;
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_axx_chi(void){
   Y nstw=0,flag,i;
   Z min,max,diff,n_t,n_st;
   Z xx,yy,uu,frc,pt,qt,rt;

   if(axx!=NULL)delete [] axx;
   axx=new Z[this->nwrd];
   cs=0;
   frc=((Z)nsx)/((Z)nnx);
   diff=nsx-nnx;

   for(i=0;i<this->nwrd;i++){
      n_t=*(tx+i);
      if(n_t&&(n_t<nnx)){
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(sx+i);
         flag=1;
         if(n_st==min){
            if(n_st-1>n_t*frc)n_st-=1;
            else flag=0;
         }
         else if(n_st<=max){
            if(max+1<n_t*frc)n_st=max+1;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         xx=n_st;
         yy=n_t*nsx/nnx;
         uu=(xx-yy)*(xx-yy)/yy;
         xx=n_t-n_st;
         yy=n_t*(nnx-nsx)/nnx;
         uu+=(xx-yy)*(xx-yy)/yy;
         xx=nsx-n_st;
         yy=(nnx-n_t)*nsx/nnx;
         uu+=(xx-yy)*(xx-yy)/yy;
         xx=nnx-n_t-nsx+n_st;
         yy=(nnx-n_t)*(nnx-nsx)/nnx;
         uu+=(xx-yy)*(xx-yy)/yy;
         axx[i]=uu;
      }
      else axx[i]=0;
   }
   if(this->pflag)cout << "chi terms " << nstw << endl;
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_mrk(Indx<Y> *pTrm){
   Y i;
   if(mrk!=NULL)delete [] mrk;
   mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++) mrk[i]=0;

   for(i=0; i<pTrm->ix; i++) mrk[pTrm->idx[i]]=1;
   if(this->pflag)cout << "Number of this->marked terms= " << pTrm->ix << endl;
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_mrk_freq(Y nm){
   Y i,nstw=0;
   if(mrk!=NULL)delete [] mrk;
   mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if(rnd(tx[i])<nm)mrk[i]=0;
      else {
         mrk[i]=1;
         nstw++;
      }
   }
   if(this->pflag)cout << "Number of freq this->marked terms= " << nstw << endl;
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_mrk_weg(Z cut){
   Y i;
   if(mrk!=NULL)delete [] mrk;
   mrk=new Y[this->nwrd];

   for(i=0;i<this->nwrd;i++){
      if(fabs(weg[i])<cut)mrk[i]=0;
      else mrk[i]=1;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::Cut_mrk_weg(Z cut){
   Y i;

   for(i=0;i<this->nwrd;i++){
      if(fabs(weg[i])<cut)mrk[i]=0;
      else mrk[i]=1;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_mrk_axx(Z cut){
   Y i;
   if(mrk!=NULL)delete [] mrk;
   mrk=new Y[this->nwrd];

   for(i=0;i<this->nwrd;i++){
      if(axx[i]<cut)mrk[i]=0;
      else mrk[i]=1;
   }
}

template<class Y,class Z>
Z *BayeQ<Y,Z>::ScoreAll(void){
   Y i,j,n;
   Z xx;

   cs=0;
   for(i=0;i<this->nwrd;i++){
      if(mrk[i])cs+=csg[i];
   }
   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];
   for(n=0;n<this->ndoc;n++){
      xx=cs;
      readp_db(n);
      for(i=0;i<this->nw;i++){
         j=this->nwd[i];
         if(mrk[j])xx+=weg[j];
      }
      sco[n]=xx;
      this->mark(n+1,1000,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
Z *BayeQ<Y,Z>::ScoreSet(Indx<Y> *ind){
   Y i,j,n,u;
   Z xx;

   cs=0;
   for(i=0;i<this->nwrd;i++){
      if(mrk[i])cs+=csg[i];
   }
   if(sco!=NULL)delete [] sco;
   sco=new Z[this->ndoc];
   for(u=0;u<ind->ix;u++){
      n=ind->idx[u];
      xx=cs;
      readp_db(n);
      for(i=0;i<this->nw;i++){
         j=this->nwd[i];
         if(mrk[j])xx+=weg[j];
      }
      sco[n]=xx;
      this->mark(u+1,1000,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_All(int n){
   Y i;
   if(mrk!=NULL)delete [] mrk;
   mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++)mrk[i]=n;
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_Char(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strchar(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_String(char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strstr(this->show(i),str))mrk[i]=n;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_NChar(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strchar(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void BayeQ<Y,Z>::Set_NString(char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strstr(this->show(i),str))mrk[i]=n;
   }
}

}
#endif
