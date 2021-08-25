#ifndef DTREEX_H
#define DTREEX_H

#include <iostream>
#include <fstream>
#include <DataObj.h>
#include <Isgrid.h>
#include <XPost.h>
#include <Dist.h>

using namespace std;
namespace iret {
template<class Y,class Z>
class TreeBaseXPost: public XPost<Y,Z>{
public:
   TreeBaseXPost(const char *nampspost);
   ~TreeBaseXPost(void);
   void init(Y n);//after gopen_Postg() 
   Y get_max_alpha(Z cut);
   Y get_max_delta(Z cut);
   Y get_max_gain(Z cut);
   Y get_min_partition(Z cut);
   void countST(Indx<Y> *ind1,Indx<Y> *ind2);
   void countDoc(Y i);
   void counbDoc(Y i);  
   void zerot(void);
   void zeros(void);
   Z ptxx; //1/n
  
   Z *tx;
   Z *sx;
   Z nnx;
   Z nsx;
   Z *dt;
   Y ncdoc;
   Z lim; //Confidential Limit 

};
//

template <class Y, class Z>
TreeBaseXPost<Y,Z>::TreeBaseXPost(const char *namspost) : XPost<Y,Z>(namspost){
  tx=NULL;
  sx=NULL;
}

template <class Y, class Z>
TreeBaseXPost<Y,Z>::~TreeBaseXPost(){
   if(tx) delete [] tx;
   if(sx) delete [] sx;
}

template <class Y, class Z>
void TreeBaseXPost<Y,Z>::init(Y n){
   ptxx=(Z)n;
   tx=new Z[this->nwrd];
   sx=new Z[this->nwrd];
   lim=0.25;
}

template <class Y, class Z>
Y TreeBaseXPost<Y,Z>::get_max_alpha(Z cut){
   Z n_t,n_st,min, max, diff;
   Y flag;
   Y nstw=0,sid=-1;
   Z xx,frc,pt,qt,rt, xxi,xxo,xxstd, eps;
   Z alx=-1000000.0,wstd;
   diff=nsx-nnx;   
   frc=(Z)nsx/(Z)nnx;
   eps=(Z)nnx/(Z)ncdoc;
  
   for(Y i=0;i<this->nwrd;i++){
      if(!(n_t=*(tx+i)))continue;
      
      if(n_t&&(n_t<nnx)){
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(sx+i);
         flag=1;
         if(n_st+eps>min){
            if(n_st-eps>n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st<max+eps){
            if(max+eps<n_t*frc)n_st=max+eps;
            else flag=0;
         }
      }
     
      else flag=0;

      if(flag){
         
         pt =(Z)n_st/(Z)nsx;
         qt =(Z)(n_t-n_st)/(Z)(nnx-nsx);
         rt =(Z)n_t/(Z)nnx;
         
         xx=-n_st*log(rt/pt) - (n_t-n_st)*log(rt/qt);
         xx-=(nsx-n_st)*log((1.0-rt)/(1.0-pt));
         xx-=(nnx-nsx-n_t+n_st)*log((1.0-rt)/(1.0-qt));
         if(alx<xx){
            wstd=log(pt*(1.0-qt))-log(qt*(1.0-pt));
            if(fabs(wstd)>=cut){
                alx=xx;
                sid=i;
            }
            
         }
         
         
         nstw++;
        
      }
      tx[i]=sx[i]=0;
   }
   
   return(sid);
}

template <class Y, class Z>
Y TreeBaseXPost<Y,Z>::get_max_delta(Z cut){
   Z n_t,n_st,min, max, diff;
   Y nstw=0,sid=-1,flag;
   Z xx,rt_l,rt_r,rt;
   Z alx=-1000000.0,wstd;
   Z frc, eps;
   
   diff=(Z)nsx-nnx;   
   frc=(Z)nsx/(Z)nnx;
   eps=(Z)nnx/(Z)ncdoc;

   for(Y i=0;i<this->nwrd;i++){
      if(!(n_t=*(tx+i))) {
         tx[i]=sx[i]=0;
         continue;
      }
      if(n_t&&(n_t<nnx)){
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(sx+i);
         flag=1;
         if(n_st+eps>min){
            if(n_st-eps>n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st<max+eps){
            if(max+eps<n_t*frc)n_st=max+eps;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){

         rt_l = (1.0/ptxx)*(n_st-(n_st*n_st)/n_t);
         rt_r = (1.0/ptxx)*((nsx-n_st) -(Z)((nsx-n_st)*(nsx-n_st))/(Z)(nnx-n_t));
         rt = (1.0/ptxx)*(nsx -(Z)(nsx*nsx)/(Z)nnx);
         
         xx=rt-rt_l-rt_r;
         if(alx<xx){
            wstd=rt;
            if(wstd>=cut){
                alx=xx;
                sid=i;
            }
            
         }
     }
     tx[i]=sx[i]=0;
   }
   return(sid);
}

template <class Y, class Z>
Y TreeBaseXPost<Y,Z>::get_max_gain(Z cut){
   Z n_t,n_st,min, max;
   Y nstw=0,sid=-1,flag;
   Z xx,px_1, px_2, py_1, py_2, info, info_x, p1, p2;
   Z alx=-1000000.0,wstd;
   Z diff, frc,eps;
   Y i;
   if(nsx==0 || nsx==nnx) {
      for(i=0;i<this->nwrd;i++) tx[i]=sx[i]=0.0;
      return(sid);
   }

   diff=(Z)nsx-(Z)nnx;   
   frc=(Z)nsx/(Z)nnx;
   eps=(Z)nnx/(Z)ncdoc;

   p1=(Z)nsx/(Z)nnx;
   p2=(Z)(nnx-nsx)/(Z)nnx;
   info=-p1*(log(p1)/log(2.0))-p2*(log(p2)/log(2));

   for(i=0;i<this->nwrd;i++){
      if(!(n_t=*(tx+i))) {
         tx[i]=sx[i]=0;
         continue;
      }
      if(n_t&&(n_t<nnx)){
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         n_st=*(sx+i);
         flag=1;
         if(n_st+eps>min){
            if(n_st-eps>n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st<max+eps){
            if(max+eps<n_t*frc)n_st=max+eps;
            else flag=0;
         }
      }
      else flag=0;
      
      if(flag){
         px_1= (Z)n_st/(Z)n_t;
         px_2= (Z)(n_t-n_st)/(Z)n_t;
         py_1= (Z)(nsx-n_st)/(Z)(nnx-n_t);
         py_2= (Z)(nnx-nsx-n_t+n_st)/(Z)(nnx-n_t);
         info_x=((Z)n_t/(Z)nnx)*(-px_1*(log(px_1)/log(2))-px_2*(log(px_2)/log(2)));
         info_x+=((Z)(nnx-n_t)/(Z)nnx)*(-py_1*(log(py_1)/log(2))-py_2*(log(py_2)/log(2)));
         
         xx=info-info_x;
         if(alx<xx){
            if(info>=cut){
                alx=xx;
                sid=i;
            }
            
         }
     }
     tx[i]=sx[i]=0;
   }
   return(sid);
}

template <class Y, class Z>
Y TreeBaseXPost<Y,Z>::get_min_partition(Z cut){
   Z rt_l, rt_r, rt, frc;
   Z max,min, diff, eps, tmax=0.0;
   Y flag;
   diff=(Z)nsx-(Z)nnx;
   frc=(Z)(nsx)/(Z)(nnx);
   eps=(Z)nnx/(Z)ncdoc;
   Z n_t,n_st;
   Z w00,w01,w10,w11, w_gd,w_bd,w0,w1,w;
   Y sid=-1;
   Y i;
   if(nsx==0 || nsx==nnx) {
      for(i=0;i<this->nwrd;i++) tx[i]=sx[i]=0.0;
      return(sid);
   }
   w_gd=nsx;
   w_bd=nnx-nsx;
   w0=2.0*sqrt(w_gd*w_bd);
   if(w0<=cut){
      for(i=0;i<this->nwrd;i++) tx[i]=sx[i]=0.0;
      return(sid);
   }
   for(i=0;i<this->nwrd;i++){
      n_t=*(tx+i);
      if(n_t&&(n_t<nnx)){
         n_st=*(sx+i);
         min=(n_t<nsx)?n_t:nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         flag=1;
         if(n_st+eps>min){
            if(n_st-eps>n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st<max+eps){
            if(max+eps<n_t*frc)n_st=max+eps;
            else flag=0;
         }
      }
      else flag=0;
      if(flag){
         w00=(nnx-nsx-n_t+n_st);
         w01=(n_t-n_st);
         w10=(nsx-n_st);
         w11=(n_st);
         w1=2.0*(sqrt(w11*w01)+sqrt(w10*w00));
         w=w0-w1;
         if(tmax<w) {
            tmax=w;
            sid=i;
         }  
      }
      tx[i]=sx[i]=0.0;

   }
   return(sid);


}

template <class Y, class Z>
void TreeBaseXPost<Y,Z>::countST(Indx<Y> *ind1,Indx<Y> *ind2){
   Y i;
   nsx=nnx=0.0;
   ncdoc=0;
   if(ind1){ 
      for(i=0;i<ind1->ix;i++){
         this->counbDoc(ind1->idx[i]);
         nsx+=dt[ind1->idx[i]];
      }
      nnx=nsx;
      ncdoc=ind1->ix;
   }
   if(ind2){
      for(i=0;i<ind2->ix;i++){
         this->countDoc(ind2->idx[i]);
         nnx+=dt[ind2->idx[i]];
      }
      ncdoc+=ind2->ix;  
   }
  
}

template <class Y, class Z>
void TreeBaseXPost<Y,Z>::countDoc(Y i){
   Y j;
   Z pt;
   pt=dt[i];
   this->readp_db(i);
   for(i=0;i<this->nw;i++){
      j=*(this->nwd+i);
      (*(tx+j))+=pt;
   }
}


template <class Y, class Z>
void TreeBaseXPost<Y,Z>::counbDoc(Y i){
   Y j;
   Z pt;
   pt=dt[i];
   this->readp_db(i);
   for(i=0;i<this->nw;i++){
      j=*(this->nwd+i);
      (*(sx+j))+=pt;
      (*(tx+j))+=pt;
   }
}

template <class Y, class Z>
void TreeBaseXPost<Y,Z>::zerot(void){
   for(Y i=0;i<this->nwrd;i++)*(tx+i)=0.0;
}

template <class Y, class Z>
void TreeBaseXPost<Y,Z>::zeros(void){
   for(Y i=0;i<this->nwrd;i++)*(sx+i)=0.0;
}

template <class KNodx,class Y,class Z>
class Param_XPost {
public:
   Param_XPost(KNodx *pKn);
   ~Param_XPost();
   void Set_alpha(void);
   void Set_delta(void);
   void Set_gain(Z lim);
   void Set_boost();
   Z n_st;
   Z n_t;
   Z nnx;
   Z nsx;
   Y ndoc;
   Z win;
   Z wout;
   Z alpha;
   Z wstd;
   Y ncdoc;  
   Z ptxx;   

};

template<class Y,class Z>
class KNode_XPost {
public:
   KNode_XPost(Indx<Y> *ind1, Indx<Y> *ind2);
   ~KNode_XPost();
   
   virtual Y Split(Z cut, Y level,TreeBaseXPost<Y,Z> *pBb)=0;
   virtual void Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx)=0;
   virtual Z Prune(void)=0;
   virtual void debug(void)=0;
   Indx<Y> *snd; //PoYer at good set 
   Indx<Y> *nnd; //PoYer at bad set
   Y split; //term number for splitting
   Z alpha; //represents improvement from splitting this node
   Z win; //score to left
   Z wout; //score to right
   Z wstd; //goodness of this node
   Z rdiff; //variable for pruning
   Z lfcnt; //number of leaves below this node
   Z n_st;
   Z n_t;
   Z nsx;
   Z nnx;
   Z ptxx; //sum of dt array for everything
   Y flio; //0 if left child, 1 if right child  
   Y dep; //depth in tree, root at 0   
};
template <class Y, class Z>
class KNode_XPost_CT : public KNode_XPost<Y,Z>  {
public:
   KNode_XPost_CT(Indx<Y> *ind1, Indx<Y> *ind2);
   ~KNode_XPost_CT();
   Y Split(Z cut,Y level, TreeBaseXPost<Y,Z> *pBb);
   void Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx);
   Z Prune(void);
   void debug(void);
   KNode_XPost_CT *up;
   KNode_XPost_CT *bin;
   KNode_XPost_CT *bout;
};

template <class Y, class Z>
class KNode_XPost_Cart : public KNode_XPost<Y,Z>  {
public:
   KNode_XPost_Cart(Indx<Y> *ind1, Indx<Y> *ind2);
   ~KNode_XPost_Cart();
   
    Y Split(Z cut,Y level, TreeBaseXPost<Y,Z> *pBb);
   void Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx);
   Z Prune(void);
   void debug(void);
   KNode_XPost_Cart *up;
   KNode_XPost_Cart *bin;
   KNode_XPost_Cart *bout;
};

template <class Y, class Z>
class KNode_XPost_C45 : public KNode_XPost<Y,Z>  {
public:
   KNode_XPost_C45(Indx<Y> *ind1, Indx<Y> *ind2);
   ~KNode_XPost_C45();
   Y Split(Z cut,Y level, TreeBaseXPost<Y,Z> *pBb);
   void Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx);
   Z Prune(void);
   void debug(void);
   KNode_XPost_C45 *up;
   KNode_XPost_C45 *bin;
   KNode_XPost_C45 *bout;
};

template <class Y, class Z>
class KNode_XPost_Boost : public KNode_XPost<Y,Z> {
public:
   KNode_XPost_Boost(Indx<Y> *ind1, Indx<Y> *ind2);
   ~KNode_XPost_Boost();
   Y Split(Z cut,Y level, TreeBaseXPost<Y,Z> *pBb);
   void Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx);
   Z Prune(void);
   void debug(void);
   KNode_XPost_Boost *up;
   KNode_XPost_Boost *bin;
   KNode_XPost_Boost *bout;
};




//
//Node
template <class Y, class Z>
KNode_XPost<Y,Z>::KNode_XPost(Indx<Y> *ind1, Indx<Y> *ind2){
   this->snd=ind1;
   nnd=ind2;
   
}

template <class Y, class Z>
KNode_XPost<Y,Z>::~KNode_XPost(){
   if(this->snd != NULL) delete this->snd;
   if(nnd != NULL) delete nnd;
}


template <class Y, class Z>
void KNode_XPost_CT<Y,Z>::debug(void){
   cout <<" address of the this->snd: " << this->snd <<endl;
   cout <<" address of the nnd: " << this->nnd <<endl;
   cout <<" address of the up: " <<  this->up <<endl;
   cout <<" address of the bin: " << this->bin <<endl;
   cout <<" address of the bout: " << this->bout <<endl;
   cout <<" term number : " << this->split <<endl;
   cout <<" alpha value : " << this->alpha <<endl;
   cout <<" win value: " << this->win <<endl;
   cout <<" wout value : " << this->wout <<endl;
   cout <<" wstd value : " << this->wstd <<endl;
}

template <class Y, class Z>
KNode_XPost_CT<Y,Z>::KNode_XPost_CT(Indx<Y> *ind1, Indx<Y> *ind2) : KNode_XPost<Y,Z>(ind1, ind2){

}

template <class Y, class Z>
KNode_XPost_CT<Y,Z>::~KNode_XPost_CT() {

}

template <class Y, class Z>
Y KNode_XPost_CT<Y,Z>::Split(Z cut, Y level,TreeBaseXPost<Y,Z> *pBb){
   
   Y flag, min, max, diff,i;
   Z pt, qt, rt, frc;
   if(this->dep==level) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }
   pBb->countST(this->snd,this->nnd);
   this->split=pBb->get_max_alpha(cut);
   
   if(this->split<0) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }
   Indx<Y> *sbin, *nbin, *sbout, *nbout,*pInd;
   if(this->snd){
      pInd=pBb->readp(this->split);
      sbin=this->snd->cbool_And(pInd);
      sbout=this->snd->cbool_Butnot(pInd);
   }
   else {
      sbin=NULL;
      sbout=NULL;
   }
   if(this->nnd){
      pInd=pBb->readp(this->split);
      nbin=this->nnd->cbool_And(pInd);
      nbout=this->nnd->cbool_Butnot(pInd);
   }
   else {
      nbin=NULL;
      nbout=NULL;
   }
   
   bin=new KNode_XPost_CT(sbin,nbin);
   bin->up=this;
   bin->flio=0;
   bin->dep=this->dep+1;
   
      

   bout=new KNode_XPost_CT(sbout,nbout);
   bout->up=this;
   bout->flio=1;
   bout->dep=this->dep+1;

   
   
   this->n_st=this->n_t=0.0;
   if(bin->snd){
      for(i=0;i<bin->snd->ix;i++) 
         (this->n_st)+=pBb->dt[bin->snd->idx[i]];
   }
   this->n_t=this->n_st;
   if(bin->nnd){
      for(i=0;i<bin->nnd->ix;i++) (this->n_t)+=pBb->dt[bin->nnd->idx[i]];
   }

   this->nsx=this->nnx=0.0;
   if(this->snd){
      for(i=0;i<this->snd->ix;i++) this->nsx+=pBb->dt[this->snd->idx[i]];
   }
   
   this->nnx=this->nsx;
   if(this->nnd){
      for(i=0;i<this->nnd->ix;i++) this->nnx+=pBb->dt[this->nnd->idx[i]];
   }     
   
   
   Param_XPost<KNode_XPost_CT<Y,Z>,Y,Z> *Pm=new Param_XPost<KNode_XPost_CT<Y,Z>, Y,Z>(this);
   Pm->ptxx=pBb->ptxx;
   Pm->Set_alpha();
   this->win=Pm->win;
   this->wout=Pm->wout;
   this->wstd=Pm->wstd;
   this->alpha=Pm->alpha;
   delete Pm;
   if(this->snd) { 
     delete this->snd;
     this->snd=NULL;
   }
  
   if(this->nnd) {
      delete this->nnd;
      this->nnd=NULL;
   }
   return(1);         
}

template <class Y, class Z>
Z KNode_XPost_CT<Y,Z>::Prune(void){
   Z xx=0;
   this->rdiff=this->alpha;
   this->lfcnt=2.0;
   if(bin){
      this->rdiff+=bin->alpha;
      this->lfcnt+=(bin->lfcnt-1.0);
   }
   if(bout){
      this->rdiff+=bout->alpha;
      this->lfcnt+=(bout->lfcnt-1.0);
   }  
   return(this->rdiff/(this->lfcnt-1.0));

}

template <class Y, class Z>
void KNode_XPost_CT<Y,Z>::Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx){
   Y i,j,k;
   Indx<Y> *ind, *ond, *pInd;
   if(this->wstd && this->snd){
      pInd=pBx->readp(this->split);
      ind=this->snd->cbool_And(pInd);
      ond=this->snd->cbool_Butnot(pInd);
      
      delete this->snd;
      this->snd=NULL;
      
      if(ind){ 
         for(i=0;i<ind->ix;i++) sxx[ind->idx[i]]+=this->win;
         
      }
      if(ond){
         for(i=0;i<ond->ix;i++) sxx[ond->idx[i]]+=this->wout;
      }
      if(bin) bin->snd=ind;
      else delete ind;
      
      if(bout) bout->snd=ond;
      else delete ond;      
      
   }
} 





//----Begin of the Cart----------------------------

template <class Y, class Z>
void KNode_XPost_Cart<Y,Z>::debug(void){
   cout <<" address of the this->snd: " <<  this->snd <<endl;
   cout <<" address of the this->nnd: " <<  this->nnd <<endl;
   cout <<" address of the up: " <<  this->up <<endl;
   cout <<" address of the bin: " <<this->bin <<endl;
   cout <<" address of the bout: " << this->bout <<endl;
   cout <<" term number : " << this->split <<endl;
   cout <<" alpha value : " << this->alpha <<endl;
   cout <<" win value: " << this->win <<endl;
   cout <<" wout value : " << this->wout <<endl;
   cout <<" wstd value : " << this->wstd <<endl;
}

template <class Y, class Z>
KNode_XPost_Cart<Y,Z>::KNode_XPost_Cart(Indx<Y> *ind1, Indx<Y> *ind2) : KNode_XPost<Y,Z>(ind1, ind2){

}

template <class Y, class Z>
KNode_XPost_Cart<Y,Z>::~KNode_XPost_Cart() {

}

template <class Y, class Z>
Y KNode_XPost_Cart<Y,Z>::Split(Z cut, Y level,TreeBaseXPost<Y,Z> *pBb){
   Y i;
   if(this->dep==level) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }

   pBb->countST(this->snd,this->nnd);
   this->split=pBb->get_max_delta(cut);
   
   if(this->split<0) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }
   Indx<Y> *sbin, *nbin, *sbout, *nbout,*pInd;
   if(this->snd){
      pInd=pBb->readp(this->split);
      sbin=this->snd->cbool_And(pInd);
      sbout=this->snd->cbool_Butnot(pInd);
   }
   else {
      sbin=NULL;
      sbout=NULL;
   }
   if(this->nnd){
      pInd=pBb->readp(this->split);
      nbin=this->nnd->cbool_And(pInd);
      nbout=this->nnd->cbool_Butnot(pInd);
   }
   else {
      nbin=NULL;
      nbout=NULL;
   }

   this->bin=new KNode_XPost_Cart(sbin,nbin);
   this->bin->up=this;
   this->bin->flio=0;
   bin->dep=this->dep+1;

   this->bout=new KNode_XPost_Cart(sbout,nbout);
   this->bout->up=this;
   this->bout->flio=1;
   this->bout->dep=this->dep+1;

   this->n_st=this->n_t=0.0;
   if(bin->snd){
      for(i=0;i<bin->snd->ix;i++) 
         (this->n_st)+=pBb->dt[bin->snd->idx[i]];
   }
   this->n_t=this->n_st;
   if(bin->nnd){
      for(i=0;i<bin->nnd->ix;i++) this->n_t+=pBb->dt[bin->nnd->idx[i]];
   }

   this->nsx=this->nnx=0.0;
   if(this->snd){
      for(i=0;i<this->snd->ix;i++) this->nsx+=pBb->dt[this->snd->idx[i]];
   }
   
   this->nnx=this->nsx;
   if(this->nnd){
      for(i=0;i<this->nnd->ix;i++) this->nnx+=pBb->dt[this->nnd->idx[i]];
   }     



   Param_XPost<KNode_XPost_Cart<Y,Z>,Y,Z> *Pm=new Param_XPost<KNode_XPost_Cart<Y,Z>,Y,Z>(this);
   Pm->ptxx=pBb->ptxx;
   Pm->Set_delta();
   
   this->win=Pm->win;
   this->wout=Pm->wout;
   this->wstd=Pm->wstd;
   this->alpha=Pm->alpha;
   delete Pm;

   if(this->snd) { 
     delete this->snd;
     this->snd=NULL;
   }
  
   if(this->nnd) {
      delete this->nnd;
      this->nnd=NULL;
   }
   return(1);   
}

template <class Y, class Z>
void KNode_XPost_Cart<Y,Z>::Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx){
   Y i,j,k;
   Indx<Y> *ind, *ond,*pInd;
   if(this->wstd && this->snd){
      pInd=pBx->readp(this->split);      
      ind=this->snd->cbool_And(pInd);
      ond=this->snd->cbool_Butnot(pInd);
      delete this->snd;
      this->snd=NULL;
      
      if(bin) bin->snd=ind;
      else if(ind){
         for(i=0;i<ind->ix;i++) sxx[ind->idx[i]]+=this->win;
         delete ind;
         
      }
      if(bout) bout->snd=ond;
      else if(ond){
         for(i=0;i<ond->ix;i++) sxx[ond->idx[i]]+=this->wout;
         delete ond;
      }
   }
} 


template <class Y, class Z>
Z KNode_XPost_Cart<Y,Z>::Prune(void){
   this->rdiff=this->alpha;
   this->lfcnt=2.0;
   if(bin){
      this->rdiff+=bin->alpha;
      this->lfcnt+=bin->lfcnt-1.0;
   }
   if(bout){
      this->rdiff+=bout->alpha;
      this->lfcnt+=bout->lfcnt-1.0;
   }  
   return(this->rdiff/(this->lfcnt-1.0));

}







//-------------------------End of Cart--------------------------




//Begin of C45//
template <class Y, class Z>
void KNode_XPost_C45<Y,Z>::debug(void){
   cout <<" address of the this->snd: " <<  this->snd <<endl;
   cout <<" address of the this->nnd: " <<  this->nnd <<endl;
   cout <<" address of the up: " <<  up <<endl;
   cout <<" address of the bin: " <<  bin <<endl;
   cout <<" address of the bout: " << bout <<endl;
   cout <<" term number : " << this->split <<endl;
   cout <<" alpha value : " << this->alpha <<endl;
   cout <<" win value: " << this->win <<endl;
   cout <<" wout value : " << this->wout <<endl;
   cout <<" wstd value : " << this->wstd <<endl;
}


template <class Y, class Z>
KNode_XPost_C45<Y,Z>::KNode_XPost_C45(Indx<Y> *ind1, Indx<Y> *ind2) : KNode_XPost<Y,Z>(ind1, ind2){

}

template <class Y, class Z>
KNode_XPost_C45<Y,Z>::~KNode_XPost_C45() {

}

template <class Y, class Z>
Y KNode_XPost_C45<Y,Z>::Split(Z cut, Y level,TreeBaseXPost<Y,Z> *pBb){
   Y i;
   if(this->dep==level) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }

   Z lim=pBb->lim;
   pBb->countST(this->snd,this->nnd);
   this->split=pBb->get_max_gain(cut);
  
   if(this->split<0) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }
      
   Indx<Y> *sbin, *nbin, *sbout, *nbout,*pInd;
   if(this->snd){
      pInd=pBb->readp(this->split);
      sbin=this->snd->cbool_And(pInd);
      sbout=this->snd->cbool_Butnot(pInd);
   }
   else {
      sbin=NULL;
      sbout=NULL;
   }
   if(this->nnd){
      pInd=pBb->readp(this->split);
      nbin=this->nnd->cbool_And(pInd);
      nbout=this->nnd->cbool_Butnot(pInd);
   }
   else {
      nbin=NULL;
      nbout=NULL;
   }

   bin=new KNode_XPost_C45<Y,Z>(sbin,nbin);
   bin->up=this;
   bin->flio=0;
   bin->dep=this->dep+1;

   bout=new KNode_XPost_C45<Y,Z>(sbout,nbout);
   bout->up=this;
   bout->flio=1;
   bout->dep=this->dep+1;
   
   this->n_st=this->n_t=0.0;
   if(bin->snd){
      for(i=0;i<bin->snd->ix;i++) 
         this->n_st+=pBb->dt[bin->snd->idx[i]];
   }
   this->n_t=this->n_st;
   if(bin->nnd){
      for(i=0;i<bin->nnd->ix;i++) this->n_t+=pBb->dt[bin->nnd->idx[i]];
   }

   this->nsx=this->nnx=0.0;
   if(this->snd){
      for(i=0;i<this->snd->ix;i++) this->nsx+=pBb->dt[this->snd->idx[i]];
   }
   
   this->nnx=this->nsx;
   if(this->nnd){
      for(i=0;i<this->nnd->ix;i++) this->nnx+=pBb->dt[this->nnd->idx[i]];
   }     

 
   Param_XPost<KNode_XPost_C45<Y,Z>,Y,Z> *Pm=new Param_XPost<KNode_XPost_C45<Y,Z>,Y,Z>(this);
   Pm->ptxx=pBb->ptxx;
   Pm->Set_gain(lim);
   this->win=Pm->win;
   this->wout=Pm->wout;
   this->wstd=Pm->wstd;
   this->alpha=Pm->alpha;
   delete Pm;
   
   if(this->snd) { 
     delete this->snd;
     this->snd=NULL;
   }
  
   if(this->nnd) {
      delete this->nnd;
      this->nnd=NULL;
   }
   return(1);   
}

template <class Y, class Z>
Z KNode_XPost_C45<Y,Z>::Prune(void){
   this->rdiff=this->alpha;
   if(bin) this->rdiff+=bin->alpha;
   if(bout) this->rdiff+=bout->alpha;
   return(this->rdiff/this->wstd);

}


template <class Y, class Z>
void KNode_XPost_C45<Y,Z>::Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx){
   Y i,j,k;
   Indx<Y> *ind, *ond, *pInd;
   if(this->wstd && this->snd){
      pInd=pBx->readp(this->split);
      ind=this->snd->cbool_And(pInd);
      ond=this->snd->cbool_Butnot(pInd);
      delete this->snd;
      this->snd=NULL;
      
      if(bin) bin->snd=ind;
      else if(ind){
         for(i=0;i<ind->ix;i++) sxx[ind->idx[i]]+=this->win;
         delete ind;
      }
      if(bout) bout->snd=ond;
      else if(ond){
         for(i=0;i<ond->ix;i++) sxx[ond->idx[i]]+=this->wout;
         delete ond;
      }
   }
} 



//End of C45




//Begin of Boost


template <class Y, class Z>
void KNode_XPost_Boost<Y,Z>::debug(void){
   cout <<" address of the this->snd: " <<  this->snd <<endl;
   cout <<" address of the this->nnd: " <<  this->nnd <<endl;
   cout <<" address of the up: " <<  up <<endl;
   cout <<" address of the bin: " << bin <<endl;
   cout <<" address of the bout: " << bout <<endl;
   cout <<" term number : " << this->split <<endl;
   cout <<" alpha value : " << this->alpha <<endl;
   cout <<" win value: " << this->win <<endl;
   cout <<" wout value : " << this->wout <<endl;
   cout <<" wstd value : " << this->wstd <<endl;
}

template <class Y, class Z>
KNode_XPost_Boost<Y,Z>::KNode_XPost_Boost(Indx<Y> *ind1, Indx<Y> *ind2) : KNode_XPost<Y,Z>(ind1, ind2){

}

template <class Y, class Z>
KNode_XPost_Boost<Y,Z>::~KNode_XPost_Boost() {

}


template <class Y, class Z>
Y KNode_XPost_Boost<Y,Z>::Split(Z cut,Y level,TreeBaseXPost<Y,Z> *pBb){
   
   Y flag, min, max, diff,i,nn;
   Z pt, qt, rt, frc;
   if(this->dep==level) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }
   pBb->countST(this->snd,this->nnd);
   this->split=pBb->get_min_partition(cut);
   if(this->split<0) {
      if(this->snd) { 
        delete this->snd;
        this->snd=NULL;
      }
      if(this->nnd) {
         delete this->nnd;
         this->nnd=NULL;
      }
      return(0);
   }
   Indx<Y> *sbin, *nbin, *sbout, *nbout,*pInd;
   if(this->snd){
      pInd=pBb->readp(this->split);
      sbin=this->snd->cbool_And(pInd);
      sbout=this->snd->cbool_Butnot(pInd);
      
   }
   else {
      sbin=NULL;
      sbout=NULL;
   }
   if(this->nnd){
      pInd=pBb->readp(this->split);
      nbin=this->nnd->cbool_And(pInd);
      nbout=this->nnd->cbool_Butnot(pInd);
      
   }
   else {
      nbin=NULL;
      nbout=NULL;
   }
   
   bin=new KNode_XPost_Boost(sbin,nbin);
   bin->up=this;
   bin->flio=0;
   bin->dep=this->dep+1;

   bout=new KNode_XPost_Boost(sbout,nbout);
   bout->up=this;
   bout->flio=1;
   bout->dep=this->dep+1;
   
   this->n_st=this->n_t=0.0;
   if(bin->snd){
      for(i=0;i<bin->snd->ix;i++) 
         this->n_st+=pBb->dt[bin->snd->idx[i]];
   }
   this->n_t=this->n_st;
   if(bin->nnd){
      for(i=0;i<bin->nnd->ix;i++) this->n_t+=pBb->dt[bin->nnd->idx[i]];
   }
   this->nsx=this->nnx=0.0;
   if(this->snd){
      for(i=0;i<this->snd->ix;i++) this->nsx+=pBb->dt[this->snd->idx[i]];
   }
   
   this->nnx=this->nsx;
   if(this->nnd){
      for(i=0;i<this->nnd->ix;i++) this->nnx+=pBb->dt[this->nnd->idx[i]];
   }     
   Param_XPost<KNode_XPost_Boost<Y,Z>,Y,Z> *Pm=new Param_XPost<KNode_XPost_Boost<Y,Z>,Y,Z>(this);
   Pm->ptxx=pBb->ptxx;
   
   
   Pm->Set_boost();
   this->alpha=Pm->alpha;
   this->win=Pm->win;
   this->wout=Pm->wout;
   this->wstd=Pm->wstd;
   delete Pm; 
   if(this->snd) { 
     delete this->snd;
     this->snd=NULL;
   }
  
   if(this->nnd) {
      delete this->nnd;
      this->nnd=NULL;
   }
   return(1);         
}

template <class Y, class Z>
Z KNode_XPost_Boost<Y,Z>::Prune(void){
   Z xx=0;
   this->rdiff=this->alpha;
   this->lfcnt=2.0; 
   if(bin){
      this->rdiff+=bin->alpha;
      this->lfcnt+=(bin->lfcnt-1.0);
   }
   if(bout){
      this->rdiff+=bout->alpha;
      this->lfcnt+=(bout->lfcnt-1.0);
   }  
   return(this->rdiff/(this->lfcnt-1.0));

}

template <class Y, class Z>
void KNode_XPost_Boost<Y,Z>::Score(Z *sxx, TreeBaseXPost<Y,Z> *pBx){
   Y i,j,k;
   Indx<Y> *ind, *ond,*pInd;
   if(this->wstd && this->snd){
      pInd=pBx->readp(this->split);
      ind=this->snd->cbool_And(pInd);
      ond=this->snd->cbool_Butnot(pInd);
      delete this->snd;
      this->snd=NULL;
      if(bin) bin->snd=ind;
      else if(ind) {
        for(i=0;i<ind->ix;i++) sxx[ind->idx[i]]+=this->win;
        delete ind;
      }
      if(bout) bout->snd=ond;
      else if(ond){
        for(i=0;i<ond->ix;i++) sxx[ond->idx[i]]+=this->wout;
        delete ond;
      }
      
   }

} 

//End of Node
//

template <class KNodx,class Y,class Z>
class Dtree_XPost : public FBase {
public:
   Dtree_XPost(Y nd, TreeBaseXPost<Y,Z> *pBx);
   Dtree_XPost(Y nd, TreeBaseXPost<Y,Z> *pBx, Y flag);
   ~Dtree_XPost();
   Y Build_Tree(Z cut, Y level,Indx<Y> *sub, Indx<Y> *nub);
   void Dest_Tree(void);
   void Score_Tree(Indx<Y> *test);
   Z Prune_Tree(void); //return (alpha_min+delta)
   Y Prune_Tree(Z level);
   void gopen_write(const char *nam);
   void write_Tree(void); //write tree to a file
   void write_Tree_Binary(void);
   void gclose_write(void);
   void Convert_AllScore_Boost(Isgrid *pIsg, Z eps);
   void Convert_LeafScore_Boost(Isgrid *pIsg, Z eps);
   Y nnode;
   Y ndoc;
   TreeBaseXPost<Y,Z> *pBb;
   KNodx *root;
   Z *sco;
   Y pflag;
   Y ntree; //Counts number of trees written out.
   ofstream *pfout; //PoYs to file being written.
};





template <class KNodx, class Y, class Z> 
Dtree_XPost<KNodx,Y,Z>::Dtree_XPost(Y nd,TreeBaseXPost<Y,Z> *pBx) : FBase("ktree","null"){
   ndoc=nd;
   pBb=pBx;
   root=NULL;
   sco=NULL;
   pflag=1;
}

template <class KNodx, class Y, class Z> 
Dtree_XPost<KNodx,Y,Z>::Dtree_XPost(Y nd,TreeBaseXPost<Y,Z> *pBx, Y flag) : FBase("ktree","null"){
   ndoc=nd;
   pBb=pBx;
   root=NULL;
   sco=NULL;
   pflag=flag;
}

template <class KNodx, class Y, class Z>
Dtree_XPost<KNodx,Y,Z>::~Dtree_XPost(){
   Y flag;
   if(sco!=NULL) delete [] sco;
   KNodx *kkn, *kko;
   if(!root) return;
   kkn=root;
   
   if(root->bin){
      kkn=root->bin;
   }else if(root->bout){
      kkn=root->bout;
   }
   else{
      delete root;
      return;
   }      
   while( kkn != root ){
      flag=1;
      while(flag){
         while(kkn->bin){
            kkn=kkn->bin;
         }
         if(kkn->bout) kkn=kkn->bout;
         else flag=0;
      }
      
      kko=kkn;  
      flag=kkn->flio;
      while(flag){
         kkn=kkn->up;
         kkn->bout=NULL;
         delete kko;
         kko=kkn;
         flag=kkn->flio; 
      }
      
      if(kkn !=root){
         kkn=kkn->up;
         kkn->bin=NULL;
         delete kko;
         if(kkn->bout) kkn=kkn->bout;
      }
   }
   delete root;
}

template <class KNodx, class Y, class Z>
Y Dtree_XPost<KNodx,Y,Z>::Build_Tree(Z cut, Y level, Indx<Y> *subb, Indx<Y> *nubb){
   Indx<Y> *sub=new Indx<Y>(subb);
   Indx<Y> *nub=new Indx<Y>(nubb);
   Y flag;
   nnode=1;
   KNodx *kko, *kkn;
   pBb->zerot();
   pBb->zeros(); 
   root = new KNodx(sub, nub);
   kkn=root; 
   kkn->flio=0;
   kkn->dep=0;
   flag=0;
   if(kkn->Split(cut,level,pBb)==0) return 0;
   kkn=kkn->bin;
   kko=NULL;
   nnode+=2;
   while( kkn != root ){
      while(kkn->Split(cut,level,pBb)){
         kkn=kkn->bin;
         nnode+=2;
         
      }
      kko=kkn;
      flag=kkn->flio;
      while(flag){
         kkn=kkn->up;
         if(kko){
            kkn->bout=NULL;
            delete kko;
            nnode--;
            kko=NULL; 
         }  
         flag=kkn->flio;
      }
      if(kkn !=root){
         kkn=kkn->up;
         if(kko){
            kkn->bin=NULL;
            delete kko;
            nnode--;
            kko=NULL;
         }  
         kkn=kkn->bout;
      }
   }
   return 1;
}

template <class KNodx, class Y, class Z>
void Dtree_XPost<KNodx,Y,Z>::Dest_Tree(void){
   Y flag;
   KNodx *kkn, *kko;
   if(root->split<0) {
      delete root; 
      return;
   }
   if(!root) return;
   kkn=root;

   if(root->bin){
      kkn=root->bin;
   }else if(root->bout){
      kkn=root->bout;
   }
   else{
      delete root;
      return;
   }
   while( kkn != root ){
      flag=1;
      while(flag){
         while(kkn->bin){
            kkn=kkn->bin;
         }
         if(kkn->bout) kkn=kkn->bout;
         else flag=0;
      }

      kko=kkn;
      flag=kkn->flio;
      while(flag){
         kkn=kkn->up;
         kkn->bout=NULL;
         delete kko;
         kko=kkn;
         flag=kkn->flio;
      }

      if(kkn !=root){
         kkn=kkn->up;
         kkn->bin=NULL;
         delete kko;
         if(kkn->bout) kkn=kkn->bout;
      }
   }
   delete root;
   root=NULL;
}

template <class KNodx, class Y, class Z>
void Dtree_XPost<KNodx,Y,Z>::Score_Tree(Indx<Y> *temp){
   Indx<Y> *test=new Indx<Y>(temp);
   Y flag,i, hflag;
   if(sco!=NULL) {
      delete [] sco;
      sco=NULL;
   }
   sco=new Z[ndoc];
   for(i=0;i<ndoc;i++) sco[i]=0.0;
   KNodx *kkn, *kko;
   kkn=NULL;
   if((root->split) < 0) return;
   if(!root) return;
   root->snd=test;
   root->Score(sco, pBb);
   if(root->bin){
      kkn=root->bin;
      hflag=0;
      kkn->Score(sco, pBb);   
   }else if(root->bout){
      kkn=root->bout;
      hflag=0;
      kkn->Score(sco, pBb);
   }
       
   while( kkn != root && kkn ){
      flag=1;
      while(flag){
         
         while(kkn->bin && !hflag){
            
            kkn=kkn->bin;
            hflag=0;
            kkn->Score(sco, pBb);
         }
         if(kkn->bout) {
            kkn=kkn->bout;
            hflag=0;
            kkn->Score(sco, pBb);
         }
         else flag=0;
      }
      
       
      flag=kkn->flio;
      while(flag){
    
         kkn=kkn->up;
         hflag=2;
         flag=kkn->flio; 
      }
      
      if(kkn !=root){
         kkn=kkn->up;
         hflag=1;
         if(kkn->bout){
            kkn=kkn->bout;
            hflag=0;
            kkn->Score(sco, pBb);
         } 

     }
   }
   
}

template <class KNodx, class Y, class Z>
Z Dtree_XPost<KNodx,Y,Z>::Prune_Tree(void){
   
   Y flag,i, hflag;
   Y fl_rem;
   Z alpha_min=1000000.0,xx, delta=1.0E-20;
   KNodx *kkn, *kko, *kkm=NULL;
   
   if(!root) {
      cout <<"root node is null"<<endl;
      return(0.0);
   }   
   if(root->bin){
      kkn=root->bin;
      hflag=0;
   }else if(root->bout){
      kkn=root->bout;
      hflag=0;
   }else{
      kkn=root;
      xx=kkn->Prune();
      if(xx < alpha_min) {
         alpha_min=xx;
         kkm=kkn;
         
      }
   }
   
   
   while( kkn != root ){
      flag=1;
      
      while(flag){
         
         while(kkn->bin && !hflag){
            kkn=kkn->bin;
            hflag=0;
         }
         if(kkn->bout) {
            kkn=kkn->bout;
            hflag=0;
         }
         else if(!hflag){
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }
            flag=0;
         }
         else{
            //cout <<"Is it possible?"<<" "<<kkn->bin<<" "<<hflag<<endl;
            if(!kkn->bin) {cout <<"should not be NULL"<<endl; exit(0);}
            if(hflag!=1) {cout <<"The current node should come from the left child"<<endl; exit(0);}
            flag=0;
         }
                           
         
      }
            
       
      flag=kkn->flio;
      while(flag){
         
         kkn=kkn->up;
         flag=kkn->flio;
         hflag=2; 
         if(!kkn->bin){
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }

         }
         
         else{
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }

         }
         
      }
       
   
      if(kkn !=root){
         
         kkn=kkn->up;
         hflag=1;
         if(kkn->bout){
            kkn=kkn->bout;
            hflag=0;
         } else{
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }

         }
         
      }
      
   }
   
   if(!kkm) { cout <<"Unexpected NULL kkm"<<endl; exit(0);}
   fl_rem=kkm->flio;
   kkm->flio=0;
   kkn=kkm;
   
         

   if(kkm->bin){
     
      kkn=kkm->bin;
   }else if(kkm->bout){
     
      kkn=kkm->bout;
   }
   else{
      
      if(kkm==root) {
         cout <<"Attempt to prune root"<<endl;
         return(0.0);
      }
      kkn=kkm->up;
      if(fl_rem) kkn->bout=NULL;
      else kkn->bin=NULL;
      delete kkm;
      return(alpha_min+delta);
      
   } 
       
   while( kkn != kkm ){
      
      flag=1;
      while(flag){
         while(kkn->bin){
            kkn=kkn->bin;
         }
         if(kkn->bout) kkn=kkn->bout;
         else flag=0;
         
      }
      
      kko=kkn;  
      flag=kkn->flio;
      while(flag){
         
         kkn=kkn->up;
         kkn->bout=NULL;
         delete kko;
         kko=kkn;
         flag=kkn->flio; 
      }
      
      if(kkn !=kkm){
         kkn=kkn->up;
         kkn->bin=NULL;
         delete kko;
         if(kkn->bout) kkn=kkn->bout;
      }
      
   }
   
   if(kkm==root) {
      cout <<"Attempt to back up from the root"<<endl;
      return(0.0);
   }
   
   kkn=kkm->up;
   if(fl_rem) kkn->bout=NULL;
   else kkn->bin=NULL;
   delete kkm;

     
   cout <<"value of minimum alpha= "<<alpha_min<<endl;
   return(alpha_min+delta);



}


  

template <class KNodx, class Y, class Z>
Y Dtree_XPost<KNodx,Y,Z>::Prune_Tree(Z level){
  
   Y flag,i, hflag;
   Y fl_rem;
   Z alpha_min=1000000.0,xx, xxx;
   
   
   KNodx *kkn, *kko, *kkm=NULL;
   
   if(!root) {
      cout <<"NULL root"<<endl;
      return(0);
   }   
   if(root->bin){
      kkn=root->bin;
      hflag=0;
   }else if(root->bout){
      kkn=root->bout;
      hflag=0;
   }else{
      kkn=root;
      xx=kkn->Prune();
      if(xx < alpha_min) {
         alpha_min=xx;
         kkm=kkn;
         
      }
   }
   
   
   while( kkn != root ){
      flag=1;
      while(flag){
         while(kkn->bin && !hflag){
            kkn=kkn->bin;
            hflag=0;
         }
         if(kkn->bout) {
            kkn=kkn->bout;
            hflag=0;
         }
         else if(!hflag){
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }
            flag=0;
         }
         else{
            if(!kkn->bin) {cout <<"should not be NULL"<<endl; exit(0);}
            if(hflag!=1) {cout <<"The current node should come from the left child"<<endl; exit(0);}
            flag=0;
         }
         
         
      }
            
       
      flag=kkn->flio;
      while(flag){
         kkn=kkn->up;
         flag=kkn->flio;
         hflag=2; 
         if(!kkn->bin){
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }

         }
         else{
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }

         }
      }
        
   
      if(kkn !=root){
         kkn=kkn->up;
         hflag=1;
         if(kkn->bout){
            kkn=kkn->bout;
            hflag=0;
         } else{
            xx=kkn->Prune();
            if(xx < alpha_min) {
               alpha_min=xx;
               kkm=kkn;
               
            }

         }
         
      }
   }
  
   if(alpha_min > level) {
      cout <<"Minimum Alpha = " << alpha_min<<endl;
      cout <<"Level to prune = " << level<<endl;
      return(0);
   }
   if(!kkm) { cout <<"Unexpected NULL kkm"<<endl; exit(0);}
   fl_rem=kkm->flio;
   kkm->flio=0;
   kkn=kkm;
   
      

   if(kkm->bin){
      kkn=kkm->bin;
   }else if(kkm->bout){
      kkn=kkm->bout;
   }
   else{
      if(kkm==root) {
         cout <<"Attempt to prune root"<<endl;
         return(0);
      }
      kkn=kkm->up;
      if(fl_rem) kkn->bout=NULL;
      else kkn->bin=NULL;
      delete kkm;
      return(1);
   }      
   while( kkn != kkm ){
      flag=1;
      while(flag){
         while(kkn->bin){
            kkn=kkn->bin;
         }
         if(kkn->bout) kkn=kkn->bout;
         else flag=0;
      }
      
      kko=kkn;  
      flag=kkn->flio;
      while(flag){
         kkn=kkn->up;
         kkn->bout=NULL;
         delete kko;
         kko=kkn;
         flag=kkn->flio; 
      }
      
      if(kkn !=kkm){
         kkn=kkn->up;
         kkn->bin=NULL;
         delete kko;
         if(kkn->bout) kkn=kkn->bout;
      }
   }
   kkn=kkm->up;
   if(fl_rem) kkn->bout=NULL;
   else kkn->bin=NULL;
   delete kkm;

     
   cout <<"value of minimum alpha= "<<alpha_min<<endl;
   return(1);



}

template <class KNodx, class Y, class Z>
void Dtree_XPost<KNodx,Y,Z>::gopen_write(const char *nam){
   set_name(nam);
   pfout=get_Ostr("t",ios::out);
   ntree=0;
}

template <class KNodx, class Y, class Z>
void Dtree_XPost<KNodx,Y,Z>::gclose_write(void){
   dst_Ostr(pfout);
   pfout=get_Ostr("n",ios::out);
   *pfout << ntree << endl;
   dst_Ostr(pfout);
}

template <class KNodx,class Y, class Z>
void Dtree_XPost<KNodx,Y,Z>::write_Tree(void){
   Y i;
   Y *pr=new Y[nnode];
   Y *lc=new Y[nnode];
   Y *rc=new Y[nnode];
   ntree++;
   for(i=0;i<nnode;i++) {
      lc[i]=-1;
      rc[i]=-1;
   }
   *pfout <<nnode<<endl;
   
   Y flag, hflag, cnod=0, curr=0;
  

   KNodx *kkn, *kko;
   kkn=NULL;
   
   if(!root) { cout <<"Tree is NULL"<<endl; return;}
   
   *pfout<<curr<<" "<<root->split<<" "<<root->win<<"  "<<root->wout<<endl;
   pr[curr]=0;
  
   if(root->bin){
      lc[curr]=++cnod;
      kkn=root->bin;
      pr[cnod]=curr;
      curr=cnod;
      *pfout<<curr <<" "<<kkn->split<<" "<<kkn->win<<"  "<<kkn->wout<<endl;
      hflag=0;
   }else if(root->bout){
      rc[curr]=++cnod;
      kkn=root->bout;
      pr[cnod]=curr;
      curr=cnod;
      *pfout<<curr<<" "<<kkn->split<<" "<<kkn->win<<"  "<<kkn->wout<<endl;
      hflag=0;
   }
       
   while( kkn != root && kkn ){
      flag=1;
      
      while(flag){
         
         while(kkn->bin && !hflag){
            lc[curr]=++cnod;
            kkn=kkn->bin;
            pr[cnod]=curr;
            curr=cnod;
            *pfout<<curr<<" "<<kkn->split<<" "<<kkn->win<<"  "<<kkn->wout<<endl;
            hflag=0;
         }
         if(kkn->bout) {
            rc[curr]=++cnod;
            kkn=kkn->bout;
            pr[cnod]=curr;
            curr=cnod;
            *pfout<<curr<<" "<<kkn->split<<" "<<kkn->win<<"  "<<kkn->wout<<endl;
            hflag=0;
         }
         else flag=0;
      }
      
       
      flag=kkn->flio;
      while(flag){
         curr=pr[curr];
         kkn=kkn->up;
         hflag=2;
         flag=kkn->flio; 
      }
      
      if(kkn !=root){
         curr=pr[curr];
         kkn=kkn->up;
         hflag=1;
         if(kkn->bout){
            rc[curr]=++cnod;
            kkn=kkn->bout;
            pr[cnod]=curr;
            curr=cnod;
            *pfout<<curr<<" "<<kkn->plit<<" "<<kkn->win<<"  "<<kkn->wout<<endl;
            hflag=0;
         } 

     }
   }
   for(i=0;i<nnode;i++) *pfout<<pr[i]<<" "<<lc[i]<<" "<<rc[i]<<endl;
  

}

template <class KNodx, class Y,class Z>
void Dtree_XPost<KNodx,Y,Z>::write_Tree_Binary(void){
   Y i;
   Y *pr=new Y[nnode];
   Y *lc=new Y[nnode];
   Y *rc=new Y[nnode];
   ntree++;
   for(i=0;i<nnode;i++) {
      lc[i]=-1;
      rc[i]=-1;
   }
   pfout->write((char*)&nnode,sizeof(Y));
   
   
   Y flag, hflag, cnod=0, curr=0;
  

   KNodx *kkn, *kko;
   kkn=NULL;
   
   if(!root) { cout <<"Tree is NULL"<<endl; return;}
   
   pfout->write((char*)&(root->split),sizeof(Y));
   pfout->write((char*)&(root->win),sizeof(Z));
   pfout->write((char*)&(root->wout),sizeof(Z));
   
   pr[curr]=0;
  
   if(root->bin){
      lc[curr]=++cnod;
      kkn=root->bin;
      pr[cnod]=curr;
      curr=cnod;
      pfout->write((char*)&(kkn->split),sizeof(Y));
      pfout->write((char*)&(kkn->win),sizeof(Z));
      pfout->write((char*)&(kkn->wout),sizeof(Z));
      
      hflag=0;
   }else if(root->bout){
      rc[curr]=++cnod;
      kkn=root->bout;
      pr[cnod]=curr;
      curr=cnod;
      pfout->write((char*)&(kkn->split),sizeof(Y));
      pfout->write((char*)&(kkn->win),sizeof(Z));
      pfout->write((char*)&(kkn->wout),sizeof(Z));
      
      hflag=0;
   }
       
   while( kkn != root && kkn ){
      flag=1;
      
      while(flag){
         
         while(kkn->bin && !hflag){
            lc[curr]=++cnod;
            kkn=kkn->bin;
            pr[cnod]=curr;
            curr=cnod;
            pfout->write((char*)&(kkn->split),sizeof(Y));
            pfout->write((char*)&(kkn->win),sizeof(Z));
            pfout->write((char*)&(kkn->wout),sizeof(Z));
            
            hflag=0;
         }
         if(kkn->bout) {
            rc[curr]=++cnod;
            kkn=kkn->bout;
            pr[cnod]=curr;
            curr=cnod;
            pfout->write((char*)&(kkn->split),sizeof(Y));
            pfout->write((char*)&(kkn->win),sizeof(Z));
            pfout->write((char*)&(kkn->wout),sizeof(Z));
            
            hflag=0;
         }
         else flag=0;
      }
      
       
      flag=kkn->flio;
      while(flag){
         curr=pr[curr];
         kkn=kkn->up;
         hflag=2;
         flag=kkn->flio; 
      }
      
      if(kkn !=root){
         curr=pr[curr];
         kkn=kkn->up;
         hflag=1;
         if(kkn->bout){
            rc[curr]=++cnod;
            kkn=kkn->bout;
            pr[cnod]=curr;
            curr=cnod;
            pfout->write((char*)&(kkn->split),sizeof(Y));
            pfout->write((char*)&(kkn->win),sizeof(Z));
            pfout->write((char*)&(kkn->wout),sizeof(Z));
            
            hflag=0;
         } 

     }
   }
   for(i=0;i<nnode;i++) { 
      pfout->write((char*)(lc+i),sizeof(Y));
      pfout->write((char*)(rc+i),sizeof(Y));
      
   }

}

template <class KNodx, class Y, class Z>
void Dtree_XPost<KNodx,Y,Z>::Convert_AllScore_Boost(Isgrid *pIsg, Z eps){
   Z deps=1.0-eps;
   Y flag,i, hflag;
   Z score=0.0, xx, px;
   KNodx *kkn, *kko;
   kkn=NULL;
   if(!root) return;
   
   if(root->bin){
      kkn=root->bin;
      score+=root->win;
      hflag=0;
   }else if(root->bout){
      xx=score+root->win;
      px=pIsg->val_1df(xx);
      if(px<eps)px=eps;
      if(px>deps)px=deps;
      root->win=0.5*log(px/(1.0-px));
      kkn=root->bout;
      score+=root->wout;
      hflag=0;
   }else{
      xx=score+root->win;
      px=pIsg->val_1df(xx);
      if(px<eps)px=eps;
      if(px>deps)px=deps;
      root->win=0.5*log(px/(1.0-px));      

      xx=score+root->wout;
      px=pIsg->val_1df(xx);
      if(px<eps)px=eps;
      if(px>deps)px=deps;
      root->wout=0.5*log(px/(1.0-px));
   }
       
   while( kkn != root && kkn ){
      flag=1;
      
      while(flag){
         
         while(kkn->bin && !hflag){
            score+=kkn->win;
            kkn=kkn->bin;
            hflag=0;
         }
         if(!hflag){
            xx=score+kkn->win;
            px=pIsg->val_1df(xx);
            if(px<eps)px=eps;
            if(px>deps)px=deps;
            kkn->win=0.5*log(px/(1.0-px));
         }
         if(kkn->bout) {
            score+=kkn->wout;
            kkn=kkn->bout;
            hflag=0;
         }
         else {
            xx=score+kkn->wout;
            px=pIsg->val_1df(xx);
            if(px<eps)px=eps;
            if(px>deps)px=deps;
            kkn->wout=0.5*log(px/(1.0-px));
            flag=0;
         }
      }
      
       
      flag=kkn->flio;
      while(flag){
         kkn=kkn->up;
         score-=kkn->wout;
         hflag=2;
         flag=kkn->flio; 
      }
      
      if(kkn !=root){
         kkn=kkn->up;
         score-=kkn->win;
         hflag=1;
         if(kkn->bout){
            score+=kkn->wout;
            kkn=kkn->bout;
            hflag=0;
         } 
         
      }
   }
}

template <class KNodx, class Y, class Z>
void Dtree_XPost<KNodx,Y,Z>::Convert_LeafScore_Boost(Isgrid *pIsg, Z eps){
   Z deps=1.0-eps;
   Y flag,i, hflag;
   Z px;
   KNodx *kkn, *kko;
   kkn=NULL;
   if(!root) return;
   
   if(root->bin){
      kkn=root->bin;
      hflag=0;
   }else if(root->bout){
      px=pIsg->val_1df(root->win);
      if(px<eps)px=eps;
      if(px>deps)px=deps;
      root->win=0.5*log(px/(1.0-px));
      kkn=root->bout;
      hflag=0;
   }else{
      px=pIsg->val_1df(root->win);
      if(px<eps)px=eps;
      if(px>deps)px=deps;
      root->win=0.5*log(px/(1.0-px));      

      px=pIsg->val_1df(root->wout);
      if(px<eps)px=eps;
      if(px>deps)px=deps;
      root->wout=0.5*log(px/(1.0-px));
   }
       
   while( kkn != root && kkn ){
      flag=1;
      while(flag){
         
         while(kkn->bin && !hflag){
            kkn=kkn->bin;
            hflag=0;
         }
         if(!hflag){
            px=pIsg->val_1df(kkn->win);
            if(px<eps)px=eps;
            if(px>deps)px=deps;
            kkn->win=0.5*log(px/(1.0-px));
         }
         if(kkn->bout) {
            kkn=kkn->bout;
            hflag=0;
         }
         else {
            px=pIsg->val_1df(kkn->wout);
            if(px<eps)px=eps;
            if(px>deps)px=deps;
            kkn->wout=0.5*log(px/(1.0-px));
            flag=0;
         }
      }
      
       
      flag=kkn->flio;
      while(flag){
         kkn=kkn->up;
         hflag=2;
         flag=kkn->flio; 
      }
      
      if(kkn !=root){
         kkn=kkn->up;
         hflag=1;
         if(kkn->bout){
            kkn=kkn->bout;
            hflag=0;
         } 
         
      }
   }
}




template <class KNodx,class Y, class Z>
Param_XPost<KNodx,Y,Z>::Param_XPost(KNodx *pKn){
    n_st=pKn->n_st;
    n_t=pKn->n_t;
    this->nsx=pKn->nsx;
    this->nnx=pKn->nnx;
    ncdoc=pKn->snd->ix+pKn->nnd->ix;
}
template <class KNodx,class Y, class Z>
Param_XPost<KNodx,Y,Z>::~Param_XPost(){
}

template <class KNodx,class Y, class Z>
void Param_XPost<KNodx,Y,Z>::Set_alpha(void){
      Z pt, qt, rt, frc, min, max, diff, eps;
      Y flag;
      
      diff=(Z)(this->nsx-this->nnx);
      frc=(Z)(this->nsx)/(Z)(this->nnx);
      eps=(Z)this->nnx/(Z)ncdoc;

      if(n_t&&(this->n_t<this->nnx)){
         min=(this->n_t<this->nsx)?this->n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         flag=1;
         if(n_st+eps>min){
            if(this->n_st-eps>this->n_t*frc)this->n_st -=eps;
            else flag=0;
         }
         else if(this->n_st<max+eps){
            if(max+eps<n_t*frc)this->n_st=max+eps;
            else flag=0;
         }
      }
      if(flag){
         
         pt =(Z)this->n_st/(Z)this->nsx;
         qt =(Z)(this->n_t-this->n_st)/(Z)(this->nnx-this->nsx);
         rt =(Z)this->n_t/(Z)this->nnx;
         
         this->win=log(pt)-log(qt);
         this->wout=log(1.0-pt)-log(1.0-qt);
         this->alpha=-this->n_st*log(rt/pt) - (this->n_t-this->n_st)*log(rt/qt);
         this->alpha-=(this->nsx-this->n_st)*log((1.0-rt)/(1.0-pt));
         this->alpha-=(this->nnx-this->nsx-this->n_t+this->n_st)*log((1.0-rt)/(1.0-qt));
         this->wstd=this->win-this->wout;

      }

}

template <class KNodx,class Y, class Z>
void Param_XPost<KNodx,Y,Z>::Set_delta(void){
      Y flag;
      Z rt_l, rt_r, rt;
      Z diff, frc,max,min, eps;
      diff=(Z)(this->nsx-this->nnx);
      frc=(Z)(this->nsx)/(Z)(this->nnx);
      eps=(Z)this->nnx/(Z)ncdoc;
     
      if(this->n_t&&(this->n_t<this->nnx)){
         min=(this->n_t<this->nsx)?this->n_t:this->nsx;
         max=this->n_t+diff;
         max=(0<max)?max:0;
         flag=1;
         if(this->n_st+eps>min){
            if(this->n_st-eps>this->n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(this->n_st<max+eps){
            if(max+eps<this->n_t*frc)this->n_st=max+eps;
            else flag=0;
         }
      } else flag=0;

      if(flag){
         rt_l=(1.0/ptxx)*(this->n_st-(Z)(this->n_st*this->n_st)/(Z)this->n_t);
         rt_r=(1.0/ptxx)*((this->nsx-this->n_st) -(Z)((this->nsx-this->n_st)*(this->nsx-this->n_st))/(Z)(this->nnx-this->n_t));
         rt =(1.0/ptxx)*(this->nsx -(Z)(this->nsx*this->nsx)/(Z)this->nnx);
         
         this->win=(Z)this->n_st/(Z)this->n_t;
         this->wout=(Z)(this->nsx-this->n_st)/(Z)(this->nnx-this->n_t);
         this->alpha=rt-rt_l-rt_r;
         this->wstd=rt;
      }
     
}

template <class KNodx,class Y, class Z>
void Param_XPost<KNodx,Y,Z>::Set_gain(Z lim){
      Y R, flag, nsx_rnd, n_st_rnd,n_t_rnd, nnx_rnd, min_rnd;
      Z p, eps=0.01, rt, rt_l, rt_r;
      Binomial b(1.0E-8,1.0E-8,1000000000);
      Z diff, frc, eps1, max, min;
      diff=(Z)(this->nsx-this->nnx);
      frc=(Z)(this->nsx)/(Z)(this->nnx);
      eps1=(Z)this->nnx/(Z)ncdoc;
     
      if(this->n_t&&(this->n_t<this->nnx)){
         min=(this->n_t<this->nsx)?this->n_t:this->nsx;
         max=this->n_t+diff;
         max=(0<max)?max:0;
         flag=1;
         if(this->n_st+eps>min){
            if(this->n_st-eps1>this->n_t*frc)this->n_st -=eps1;
            else flag=0;
         }
         else if(this->n_st<max+eps){
            if(max+eps1<n_t*frc)this->n_st=max+eps1;
            else flag=0;
         }
      } else flag=0;

      if(flag){
         this->win=(Z)this->n_st/(Z)this->n_t;
         this->wout=(Z)(this->nsx-this->n_st)/(Z)(this->nnx-this->n_t);
        
         //cout <<this->nnx<<" "<<ncdoc <<endl;
         nsx_rnd=rnd(((1.0*ncdoc)/this->nnx)*this->nsx);
         n_st_rnd=rnd(((1.0*ncdoc)/this->nnx)*this->n_st);
         n_t_rnd=rnd(((1.0*ncdoc)/this->nnx)*this->n_t);
         nnx_rnd=ncdoc;
         
         
         min_rnd=(nsx_rnd<(nnx_rnd-nsx_rnd))?nsx_rnd:(nnx_rnd-nsx_rnd);
         p=b.upper_limit(nnx_rnd,min_rnd,lim);
         rt=nnx_rnd*p;
         //cout <<"# of errors="<<min_rnd<<"|"<<rt<<endl;

         min_rnd=(n_st_rnd<(n_t_rnd-n_st_rnd))?n_st_rnd:(n_t_rnd - n_st_rnd);
         p=b.upper_limit(n_t_rnd, min_rnd, lim);
         rt_l=n_t_rnd*p;
         //cout <<"# of left child errors="<<min_rnd<<"|"<<rt_l<<endl;

         min_rnd=((nsx_rnd-n_st_rnd)<(nnx_rnd-nsx_rnd-n_t_rnd+n_st_rnd))?(nsx_rnd-n_st_rnd):(nnx_rnd-nsx_rnd-n_t_rnd+n_st_rnd);
         p=b.upper_limit(nnx_rnd-n_t_rnd, min_rnd, lim);
         rt_r=(nnx_rnd-n_t_rnd)*p;
         //cout <<"# of right child errors="<<min_rnd<<"|"<<rt_r<<endl;
         
         alpha=rt-rt_l-rt_r;
         //cout <<"# of error difference="<<alpha<<endl;
         wstd=this->nnx;

      }
     
}

template <class KNodx,class Y, class Z>
void Param_XPost<KNodx,Y,Z>::Set_boost(void){
      Z w00,w01,w10,w11;
      Z rt_l, rt_r, rt, frc;

      Z pt, qt,max,min, diff, eps;
      Y flag;
      
      diff=this->nsx-this->nnx;
      frc=(Z)(this->nsx)/(Z)(this->nnx);
      eps=(Z)this->nnx/(Z)ncdoc;
      
      if(n_t&&(n_t<this->nnx)){
         min=(n_t<this->nsx)?n_t:this->nsx;
         max=n_t+diff;
         max=(0<max)?max:0;
         flag=1;
         if(n_st+eps>min){
            if(n_st-eps>n_t*frc)n_st -=eps;
            else flag=0;
         }
         else if(n_st<max+eps){
            if(max+eps<n_t*frc)n_st=max+eps;
            else flag=0;
         }
      }
      else flag=0;

      if(flag){
         
         w00=this->nnx-this->nsx-n_t+n_st;
         w01=n_t-n_st;
         w10=this->nsx-n_st;
         w11=n_st;
         
         win=(1.0/2.0)*log(w11/w01);
         wout=(1.0/2.0)*log(w10/w00);
         
          
         rt_l=(1.0/ptxx)*(n_st-(n_st*n_st)/n_t);
         rt_r=(1.0/ptxx)*((this->nsx-n_st) -(Z)((this->nsx-n_st)*(this->nsx-n_st))/(Z)(this->nnx-n_t));
         rt =(1.0/ptxx)*(this->nsx -(Z)(this->nsx*this->nsx)/(Z)this->nnx);
         
         
         alpha=rt-rt_l-rt_r;
         wstd=rt;

      }

}
    
}
#endif
