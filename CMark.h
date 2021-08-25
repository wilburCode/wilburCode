#ifndef CMRK_H
#define CMRK_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <DataObj.h>
#include <XPost.h>

using namespace std;
namespace iret {

//template class CMark 
template<class Y,class Z>
class CMark : public XPost<Y,Z> {
public:
   CMark(const char *nspost);//name of XPost set
   CMark(const char *nspost,const char *pnam);//name of XPost set
      //pnam is either used in "path_pnam" as place to find path
      //or if begins with ':' is followed by path itself.
   ~CMark(void);

      //Counting functions.
   void init_cnt(void); //Sets pdoc uniform.
      //Sets memory for tx and sx.
      //must call gopen_db_map() of XPost first
   void zerot(void); //Zeroes tx array.
   void freq_tx(void); //Puts freq values in tx array
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
   void Set_String(const char *str,int n); //Sets mrk n if str substring of string
   void Set_NString(const char *str,int n); //Sets mrk n if str not substring of string

   void Create_Inverse(void); //replaces mrk[i]!=0 by its count among !=0 in order
      //also sets up the array *inv as an inverse mapping so inv[mrk[i]]=i.
   //Relate to the terms processing
   Z eps; //Size of 1/ndoc.
   Z nnx;
   Z nsx;
   Z *tx;
   Z *sx;
   Z *pdoc;
   Y *mrk; //Marks which terms to include in process
   Y *inv; //Inverse mapping (see Create_Inverse func)
   Y tmk;  //Number of words marked.
};

template<class Y,class Z>
CMark<Y,Z>::CMark(const char *namspost) : XPost<Y,Z>(namspost){
   mrk=NULL;
   inv=NULL;
   tx=NULL;
   sx=NULL;
   pdoc=NULL;
}

template<class Y,class Z>
CMark<Y,Z>::CMark(const char *namspost,const char *pnam) : XPost<Y,Z>(namspost,pnam){
   mrk=NULL;
   inv=NULL;
   tx=NULL;
   sx=NULL;
   pdoc=NULL;
}

template<class Y,class Z>
CMark<Y,Z>::~CMark(){
   if(mrk!=NULL)delete [] mrk;
   if(inv!=NULL)delete [] inv;
   if(tx!=NULL)delete [] tx;
   if(sx!=NULL)delete [] sx;
   if(pdoc!=NULL)delete [] pdoc;
}

template<class Y,class Z>
void CMark<Y,Z>::init_cnt(void){
   Y i;
   eps=1.0;

   if(tx==NULL)tx=new Z[this->nwrd];

   if(sx==NULL)sx=new Z[this->nwrd];

   if(pdoc==NULL)pdoc=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)pdoc[i]=eps;
   if(this->pflag)cout<<"end of init\n";
}

template<class Y,class Z>
void CMark<Y,Z>::zerot(void){
   nnx=0.0;
   for(Y i=0;i<this->nwrd;i++)*(tx+i)=0;
}

template<class Y,class Z>
void CMark<Y,Z>::zeros(void){
   nsx=0.0;
   for(Y i=0;i<this->nwrd;i++)*(sx+i)=0;
}

template<class Y,class Z>
void CMark<Y,Z>::freq_tx(void){
   this->gopen_map();
   nnx=this->ndoc;
   for(Y i=0;i<this->nwrd;i++)*(tx+i)=this->freq[i];
}

template<class Y,class Z>
void CMark<Y,Z>::countDoc(Y i){
   Y j,k;
   Z pt;
   pt=pdoc[i];

   this->readp_db(i);
   for(k=0;k<this->nw;k++){
      j=*(this->nwd+k);
      (*(tx+j))+=pt;
   }
}

template<class Y,class Z>
void CMark<Y,Z>::counsDoc(Y i){
   Y j,k;
   Z pt;
   pt=pdoc[i];

   this->readp_db(i);
   for(k=0;k<this->nw;k++){
      j=*(this->nwd+k);
      (*(sx+j))+=pt;
   }
}

template<class Y,class Z>
void CMark<Y,Z>::counbDoc(Y i){
   Y j,k;
   Z pt;
   pt=pdoc[i];

   this->readp_db(i);
   for(k=0;k<this->nw;k++){
      j=*(this->nwd+k);
      (*(sx+j))+=pt;
      (*(tx+j))+=pt;
   }
}

template<class Y,class Z>
void CMark<Y,Z>::countTX(Indx<Y> *cnd){
   Y i,j;

   for(i=0;i<cnd->ix;i++){
      j=cnd->idx[i];
      this->countDoc(j);
      nnx+=pdoc[j];
      this->mark(i,1000,"docs");
   }
}

template<class Y,class Z>
void CMark<Y,Z>::countSX(Indx<Y> *cnd){
   Y i,j;

   for(i=0;i<cnd->ix;i++){
      j=cnd->idx[i];
      this->counsDoc(j);
      nsx+=pdoc[j];
      this->mark(i,1000,"docs");
   }
}

template<class Y,class Z>
void CMark<Y,Z>::countBX(Indx<Y> *cnd){
   Y i,j;

   for(i=0;i<cnd->ix;i++){
      j=cnd->idx[i];
      this->counbDoc(j);
      nnx+=pdoc[j];
      nsx+=pdoc[j];
      this->mark(i,1000,"docs");
   }
}

template<class Y,class Z>
void CMark<Y,Z>::Set_Term_wt(Z cut){
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
void CMark<Y,Z>::Set_Term_alp(Z cut){
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
void CMark<Y,Z>::Set_Term_chi(Z cut){
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
void CMark<Y,Z>::Set_Term_muti(Z cut){
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
void CMark<Y,Z>::Set_Term(Indx<Y> *pTrm){
   Y i;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++) mrk[i]=0;
   
   for(i=0; i<pTrm->ix; i++) mrk[pTrm->idx[i]]=1;
   if(this->pflag)cout << "Number of marked terms= " << pTrm->ix << endl;
}

template<class Y,class Z>
void CMark<Y,Z>::Set_Term_freq(Y nm){
   Y i,nstw=0;
   Z xx=eps*nm;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++){
      if(tx[i]<xx)mrk[i]=0;
      else {
         mrk[i]=1;
         nstw++;
      }
   }
   if(this->pflag)cout << "Number of freq marked terms= " << nstw << endl;
}

template<class Y,class Z>
void CMark<Y,Z>::Set_Term_freq2(Y nm){
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
void CMark<Y,Z>::Set_Term_between(Y nm,Y bg){
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
void CMark<Y,Z>::Set_Term_between2(Y nm,Y bg){
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
void CMark<Y,Z>::Set_All(int n){
   Y i;
   if(mrk==NULL)mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++)mrk[i]=n;
}

template<class Y,class Z>
void CMark<Y,Z>::Set_Char(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strchr(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void CMark<Y,Z>::Set_String(const char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strstr(this->show(i),str))mrk[i]=n;
   }
}

template<class Y,class Z>
void CMark<Y,Z>::Set_NChar(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strchr(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void CMark<Y,Z>::Set_NString(const char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strstr(this->show(i),str))mrk[i]=n;
   }
}

template<class Y,class Z>
void CMark<Y,Z>::Create_Inverse(void){
   Y i,j,k;

   k=0;
   for(i=0;i<this->nwrd;i++)if(mrk[i])k++;
   inv=new Y[k];
   tmk=k;
   k=1;
   for(i=0;i<this->nwrd;i++){
      if(mrk[i]){
         mrk[i]=k;
         inv[k-1]=i;
         k++;
      }
   }
}   

}
#endif
