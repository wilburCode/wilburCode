#ifndef VNAB_H
#define VNAB_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <DataObj.h>
#include <Postg.h>
#include <XPost.h>
#include <Isgrid.h>
using namespace std;
namespace iret {

const double l2=1.0/log(2.0),dmt=20.0,lfac=log(0.65);
const double lfab=log(0.7);
extern double xttrc,xb,xa,xg;

//Local Weighting Functions

float q_const(int lc,long n); //constant 1.0

float d_lc_func(int lc,long n);  //Just returns lc

float d_lc_ratio(int lc,long n); //simple ratio lc/(2+lc)

float d_lc_log(int lc,long n); //simple ln(lc+1)

float s_const(int lc,unsigned int n); //constant 1.0

float sx_const(int lc, long n); //constant 1.0

float s_lc_func(int lc,unsigned int n);  //Just returns lc

float s_lc_ratio(int lc,unsigned int n); //simple ratio lc/(2+lc)

float s_lc_log(int lc,unsigned int n); //simple ln(lc+1)

   //Trec local functions.
float d_len_trec(int lc,long n); //func of len only

float d_inquery_trec(int lc,long n); //Inquery formula (Robertson)

float d_robertson_trec(int lc,long n); //Robertson

float d_wilbur_trec(int lc,long n); //Wilbur Poisson formula

float d_experiment_trec(int lc,long n); //Experimental

   //MED local functions
float d_inquery_med(int lc,long n); //Inquery formula (Robertson)

float d_wilbur_med(int lc,long n); //Wilbur Poisson formula

float d_string(int lc,long n); //Wilbur Poisson formula

float d_prob(int lc,long n); //Simple probability formula

float s_inquery_med(int lc,unsigned int n); //Inquery formula (Robertson)

float s_wilbur_med(int lc,unsigned int n); //Wilbur Poisson formula

float s_string(int lc,unsigned int n); //Wilbur Poisson formula

float s_prob(int lc,unsigned int n); //Simple probability formula

float d_bm25(int lc,long n); //Okapi bm25 formula

//Global Weighting Functions

float global_idf(long n); //IDF 

float global_strict_idf(long n); //IDF but uses exact frequencies

float global_sidf(long n); //IDF but scores even single occuring terms 

float global_const(long n); //Constant 1.0

float global_iqf(long n); //Iqf for Zhiyong

float global_bm25(long n); //Okapi bm25 formula


float sglobal_idf(unsigned int n); //IDF

float sglobal_sidf(unsigned int n); //IDF but scores even single occuring terms

float sglobal_const(unsigned int n); //Constant 1.0

template<class Y,class Z>
class VnbX : public XPost<Y,Z> {
public:
   VnbX(const char *nspost);//name of XPost set
   VnbX(const char *nspost,const char *pnam);//name of XPost set
      //pnam is suffix of path name file.

   ~VnbX(void);

   void create_Weight(float (*global)(long)); //Creates a file
     //.wx with global weights in it
   void create_SWeight(float (*global)(long)); //Creates a file
     //.swx with sqrt of global weights in it
   void gopen_Weight(void); //Maps the wx array
   void gopen_SWeight(void); //Maps the swx array
   void gclose_Weight(void); //Unmaps wx
   void gclose_SWeight(void); //Unmaps swx
   Z    weightG(Ordr<Y,Z> *pOrd,Indx<Y> *pInd); //Looks at intersection of pInd with
     //pOrd and uses size of pOrd as total and gives IDF weight based on it for pInd.

   void CreateLenCos(float (*global)(long)); //Creates a file 
     //.len with document lengths in it. 
     //Cosine form
   void CreateLenCos(Z *wx); //Creates a file 
     //.len with document lengths in it. 
     //Uses wx values for global weights
     //Cosine form
   void CreateLen1(float (*global)(long)); //Creates a file 
     //.len with document lengths in it. 
     //Intended for use with RevisProb1 - harmonic mean
   void CreateLen1(Z *wx); //Creates a file 
     //.len with document lengths in it. 
     //Intended for use with RevisProb1 - harmonic mean
     //Uses wx as replacement for global weights
   void CreateLen2(float (*global)(long)); //Creates a file 
     //.len with document lengths in it. 
     //Intended for use with RevisProb2 - arithmetic mean
   void CreateLen3(float (*global)(long)); //Creates a file 
     //.len with document lengths in it. 
     //Intended for use with RevisProb3 - geometric mean
   void CreateLenProd(void); //Creates a file .len with
     //document lengths as sqrt(N) where N is number of nonzero components
     //Intended for use with RevisProb4
   
   void Load(Y i); //Sets pointers for the ith doc.
     //gopen_db_map must be called first
   void Convert(Doc<Z> *pDoc,float (*d_local)(int,long)); 
     //Creates space, fills, and sets pointers
     //gopen_hash must be called first
   void Convert(Doc<Z> *pDoc); //Assumes local wts in *pDoc 
     //Creates space, fills, and sets pointers
     //gopen_hash must be called first

   void ScoreAll(float (*global)(long)); //Scores all docs
     //Assumes that doc space is filled. 
     //gopen_map must be called first
     //Uses local weights
   void ScoreAllWg(Ordr<Y,Z> *pOrd); //Scores all docs
     //Uses the function weightG to obtain weights for terms loaded
   void ScoreAllmrk(float (*global)(long)); //Scores all docs
     //Assumes that doc space is filled. 
     //gopen_map must be called first
     //Uses local weights. Only includes marked terms in score
   void ScoreAllmrkLoc(float (*global)(long),float (*d_local)(int,long),int *lnt);
     //Assumes that doc space is filled. 
     //gopen_map must be called first
     //Uses local weighting function. Only includes marked terms in score
   void ScoreAllmrkLearned(float (*global)(long),int *lnt,double **wx);
     //Special for an experiment described in LCWT in notebook
   void ScoreAllmrkLearned(Z *wg,int *lnt,double **wx);
     //Special for an experiment described in LCWT in notebook
   void ScoreAllEx(Z *wx); //Scores all docs
     //Assumes that doc space is filled. 
     //gopen_map must be called first
     //Uses local weights
   Y ScoreAllExTop(Z *wx); //Scores all docs
     //Assumes that doc space is filled. 
     //gopen_map must be called first
     //Uses local weights. Returns # of times top score seen
   void ScoreAllExmrk(Z *wx); //Scores all docs
     //Assumes that doc space is filled. 
     //gopen_map must be called first
     //Uses local weights. Only includes marked terms in score
   void ScoreAll(Z *wx); //Scores all docs
     //Assumes that doc space is filled. 
     //gopen_map must be called first
     //Uses wx in place of global function.
     //No local weights
   void ScoreAll(BCount<Z> &Dc); //Scores all docs using terms and
     //numbers in Dc as weights. No local weights
   void ScoreAll(Ordr<Y,Z> *pOrd); //Scores all docs using terms listed
     //in pOrd and their scores in pOrd as global weights. Uses local
     //weights also
   void ScoreAllNlw(Ordr<Y,Z> *pOrd); //Scores all docs using terms listed
     //in pOrd and their scores in pOrd as global weights. No local
     //weights used
   void ScoreAll(Z *sx,Z *wx); //Scores all docs
     //but does not use local weights in docs. wx is global
     //weight array and sx must be an initialized score array
   Z ScoPair(Y i,Y j); //i and j are doc #s in system. Requires that
     //Weights and Lengths be created and opened before running
   Z ScoPairN(Y i,Y j); //i and j are doc #s in system. All weighting must
     //be in the local wts at creation for this function
     //Only marked terms are used in the scoring
   Z ScoPairM(Y i,Y j); //i and j are doc #s in system. All weighting is 
     //in the global wts at creation for this function
     //Only marked terms are used in the scoring

   void gopen_Len(void); //Maps the xlen data produced by CreateLen
     //functions to xlen
     //Call before any Revis functions
   void gclose_Len(void); //Unmaps xlen
   void RevisCos(float (*global)(long)); //Assumes Load 
     //or Convert has been called and ScoreAll has been called. 
     //Produces standard cosine score
     //Assumes CreateLenCos called before
   void RevisCos(Z *wx); //Assumes Load 
     //or Convert has been called and ScoreAllEx has been called. 
     //Produces standard cosine score
     //Assumes CreateLenCos called before
   void RevisX(Z *wx); //Assumes Load 
     //or Convert has been called and ScoreAllEx has been called. 
     //Compares inner product to each len squared and takes smaller
     //Assumes CreateLenCos called before
   void RevisProb1(float (*global)(long)); //Assumes Load 
     //or Convert has been called and ScoreAll has been called. 
     //Converts raw scores to a prob in a symmetric manner. 
     //Just a form of score normalization
     //Create_Len1 must have also been run to create the lengths
   void RevisProb1(Z *wx); //Assumes Load 
     //or Convert has been called and ScoreAll has been called. 
     //Converts raw scores to a prob in a symmetric manner. 
     //Just a form of score normalization
     //Create_Len1 must have also been run to create the lengths
   void RevisProb2(float (*global)(long)); //Assumes Load 
     //or Convert has been called and ScoreAll has been called. 
     //Converts raw scores to a prob in a symmetric manner. 
     //Just a form of score normalization
     //Create_Len2 must have also been run to create the lengths
   void RevisProb3(float (*global)(long)); //Assumes Load 
     //or Convert has been called and ScoreAll has been called. 
     //Converts raw scores to a prob in a symmetric manner. 
     //Just a form of score normalization
     //Create_Len3 must have also been run to create the lengths
   void RevisProd(Z x); //Assuemes scoring has been done
     //Create_LenProd must also have been called
   Ordr<Y,Z> *Skim(Y n); //Skims off the top n scoring docs
   Ordr<Y,Z> *Skimgr(Z s); //Skims off the docs scoring >s
   Ordr<Y,Z> *Skimge(Z s); //Skims off the docs scoring >=s

   //Special mrk array functions
   void Set_All(int n); //Creates and sets mrk array to be n everywhere
   void Set_Char(char c,int n); //Sets mrk n if character c in string
   void Set_NChar(char c,int n); //Sets mrk n if character c not in string
   void Set_String(const char *str,int n); //Sets mrk n if str substring of string
   void Set_NString(const char *str,int n); //Sets mrk n if str not substring of string

   //Data

   int *mrk; //Marks which terms to include in process
   long sflag; //Marks status of local doc space
      //0 unset, 1 from load, 2 from convert
   Z *sco; //score array for whole space
   Z *xlen; //Length data
   Z *wx; //Global weight data
   Z *swx; //sqrt of Global weight data
};

template<class Y,class Z> 
VnbX<Y,Z>::VnbX(const char *nspost) : XPost<Y,Z>(nspost){
   sco=NULL;
   sflag=0;
   mrk=NULL;
}

template<class Y,class Z> 
VnbX<Y,Z>::VnbX(const char *nspost,const char *pnam) : XPost<Y,Z>(nspost,pnam){
   sco=NULL;
   sflag=0;
   mrk=NULL;
}

template<class Y,class Z> 
VnbX<Y,Z>::~VnbX(void){
   if(sco)delete [] sco;
}

template<class Y,class Z> 
void VnbX<Y,Z>::create_Weight(float (*global)(long)){
   Y i,j,k;
   Z xx,yy,zz;

   this->gopen_map();
   xttrc=(double)this->ndoc;
   ofstream *pfout=this->get_Ostr("wx");
   for(i=0;i<this->nwrd;i++){
      xx=(Z)global((long)this->freq[i]);
      pfout->write((char*)&xx,sizeof(Z));
      XPost<Y,Z>::mark(i,10000,"weights");
   }
   this->dst_Ostr(pfout);
}

template<class Y,class Z> 
void VnbX<Y,Z>::create_SWeight(float (*global)(long)){
   Y i,j,k;
   Z xx,yy,zz;

   this->gopen_map();
   xttrc=(double)this->ndoc;
   ofstream *pfout=this->get_Ostr("swx");
   for(i=0;i<this->nwrd;i++){
      xx=(Z)sqrt(global((long)this->freq[i]));
      pfout->write((char*)&xx,sizeof(Z));
      XPost<Y,Z>::mark(i,10000,"weights");
   }
   this->dst_Ostr(pfout);
}

template<class Y,class Z> 
void VnbX<Y,Z>::gopen_Weight(void){
   wx=(Z*)this->get_Mmap("wx");
}

template<class Y,class Z> 
void VnbX<Y,Z>::gopen_SWeight(void){
   swx=(Z*)this->get_Mmap("swx");
}

template<class Y,class Z> 
void VnbX<Y,Z>::gclose_Weight(void){
   this->dst_Mmap("wx",(char*&)wx);
}

template<class Y,class Z> 
void VnbX<Y,Z>::gclose_SWeight(void){
   this->dst_Mmap("swx",(char*&)swx);
}

template<class Y,class Z> 
Z  VnbX<Y,Z>::weightG(Ordr<Y,Z> *pOrd,Indx<Y> *pInd){
   Y i,j,k;
   Indx<Y> *pJnd=pInd->cbool_And((Indx<Y> *)pOrd);
   if(pJnd){
      i=pJnd->ix;
      if(i<20)i=20;
      delete pJnd;
   }
   else i=20;
   return(l2*log(pOrd->ix/((Z)i)));
}

template<class Y,class Z> 
void VnbX<Y,Z>::Load(Y i){
   if(sflag==2){
      delete [] this->nwd;
      delete [] this->lwt;
   }
   this->gopen_db_map();
   this->readp_db(i);
   this->readz_db(i);
   sflag=1;
}

template<class Y,class Z> 
void VnbX<Y,Z>::Convert(Doc<Z> *pDoc,float (*d_local)(int,long)){
   Y i,j,k,u;
   i=pDoc->nw;
   this->gopen_hash();
   if(sflag==2){
      delete [] this->nwd;
      delete [] this->lwt;
   }
   if(i){
      this->nwd=new Y[i];
      this->lwt=new Z[i];
      sflag=2;
   }
   else sflag=0;
   k=0;
   for(j=0;j<i;j++){
      u=this->find(pDoc->word[j]);
      if(u){
         this->nwd[k]=(Y)(u-1);
         this->lwt[k]=(Z)d_local((int)pDoc->lcnt[j],i);
         k++;
      }
   }
   this->nw=(Y)k;
}

template<class Y,class Z>
void VnbX<Y,Z>::Convert(Doc<Z> *pDoc){
   Y i,j,k,u;
   i=pDoc->nw;
   this->gopen_hash();
   if(sflag==2){
      delete [] this->nwd;
      delete [] this->lwt;
   }
   if(i){
      this->nwd=new Y[i];
      this->lwt=new Z[i];
      sflag=2;
   }
   else sflag=0;
   k=0;
   for(j=0;j<i;j++){
      u=this->find(pDoc->word[j]);
      if(u){
         this->nwd[k]=(Y)(u-1);
         this->lwt[k]=(Z)pDoc->lcnt[j];
         k++;
      }
   }
   this->nw=(Y)k;
}

template<class Y,class Z> 
void VnbX<Y,Z>::CreateLenCos(float (*global)(long)){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nwrd;i++){
      pInd=this->readp(i);
      flw=this->readz(i);
      j=pInd->ix;
      xx=(Z)global(j);
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k]*flw[k];
      }
      XPost<Y,Z>::mark(i,1000,"terms");
   }
   for(i=0;i<this->ndoc;i++){
      if(sco[i]>0)sco[i]=1.0/(Z)sqrt((double)sco[i]);
      else sco[i]=1.0;
   }
   bin_Writ("len",this->ndoc*sizeof(Z),(char*)sco);
}

template<class Y,class Z> 
void VnbX<Y,Z>::CreateLenCos(Z *wx){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nwrd;i++){
      pInd=this->readp(i);
      flw=this->readz(i);
      j=pInd->ix;
      xx=(Z)wx[i];
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k]*flw[k];
      }
      XPost<Y,Z>::mark(i,1000,"terms");
   }
   for(i=0;i<this->ndoc;i++){
      if(sco[i]>0)sco[i]=1.0/(Z)sqrt((double)sco[i]);
      else sco[i]=1.0;
   }
   bin_Writ("len",this->ndoc*sizeof(Z),(char*)sco);
}

template<class Y,class Z> 
void VnbX<Y,Z>::CreateLen1(float (*global)(long)){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nwrd;i++){
      pInd=this->readp(i);
      flw=this->readz(i);
      j=pInd->ix;
      xx=global(j);
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k];
      }
      XPost<Y,Z>::mark(i,1000,"terms");
   }
   for(i=0;i<this->ndoc;i++){
      if(sco[i]>0)sco[i]=0.5/sco[i];
      else sco[i]=1.0;
   }
   this->bin_Writ("len",this->ndoc*sizeof(Z),(char*)sco);
}

template<class Y,class Z> 
void VnbX<Y,Z>::CreateLen1(Z *wx){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nwrd;i++){
      pInd=this->readp(i);
      flw=this->readz(i);
      j=pInd->ix;
      xx=(Z)wx[i];
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k];
      }
      XPost<Y,Z>::mark(i,1000,"terms");
   }
   for(i=0;i<this->ndoc;i++){
      if(sco[i]>0)sco[i]=0.5/sco[i];
      else sco[i]=1.0;
   }
   bin_Writ("len",this->ndoc*sizeof(Z),(char*)sco);
}

template<class Y,class Z> 
void VnbX<Y,Z>::CreateLen2(float (*global)(long)){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nwrd;i++){
      pInd=this->readp(i);
      flw=this->readz(i);
      j=pInd->ix;
      xx=(Z)global(j);
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k];
      }
      XPost<Y,Z>::mark(i,1000,"terms");
   }
   for(i=0;i<this->ndoc;i++){
      if(sco[i]>0)sco[i]=0.5*sco[i];
      else sco[i]=1.0;
   }
   bin_Writ("len",this->ndoc*sizeof(Z),(char*)sco);
}

template<class Y,class Z> 
void VnbX<Y,Z>::CreateLen3(float (*global)(long)){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nwrd;i++){
      pInd=this->readp(i);
      flw=this->readz(i);
      j=pInd->ix;
      xx=(Z)global(j);
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k];
      }
      XPost<Y,Z>::mark(i,1000,"terms");
   }
   for(i=0;i<this->ndoc;i++){
      if(sco[i]==0)sco[i]=1.0;
   }
   bin_Writ("len",this->ndoc*sizeof(Z),(char*)sco);
}


template<class Y,class Z> 
void VnbX<Y,Z>::CreateLenProd(void){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nwrd;i++){
      pInd=this->readp(i);
      flw=this->readz(i);
      for(k=0;k<pInd->ix;k++){
         if(flw[k]>0)sco[pInd->idx[k]]++;
      }
      XPost<Y,Z>::mark(i,1000,"terms");
   }
   for(i=0;i<this->ndoc;i++){
      sco[i]=sqrt(sco[i]);
   }
   this->bin_Writ("len",this->ndoc*sizeof(Z),(char*)sco);
}

template<class Y,class Z>
Z VnbX<Y,Z>::ScoPair(Y mi,Y mj){
   Y i,j,k;
   Z *flw,xx,sum=0;
   this->gopen_db_map();
   Load(mi);
   k=this->nw;
   Y *nwk=this->nwd;
   Z *lwk=this->lwt;
   Load(mj);
   i=j=0;
   while((i<k)&&(j<this->nw)){
      if(nwk[i]<this->nwd[j])i++;
      else if(this->nwd[j]<nwk[i])j++;
      else {
         sum+=lwk[i]*(this->lwt[j])*wx[nwk[i]];
         i++;j++;
      }
   }
   return(sum*xlen[mi]*xlen[mj]);
}

template<class Y,class Z>
Z VnbX<Y,Z>::ScoPairN(Y mi,Y mj){
   Y i,j,k;
   Z *flw,xx,sum=0;
   this->gopen_db_map();
   Load(mi);
   k=this->nw;
   Y *nwk=this->nwd;
   Z *lwk=this->lwt;
   Load(mj);
   i=j=0;
   while((i<k)&&(j<this->nw)){
      if(nwk[i]<this->nwd[j])i++;
      else if(this->nwd[j]<nwk[i])j++;
      else {
         if(mrk[nwk[i]])sum+=lwk[i]*(this->lwt[j]);
         i++;j++;
      }
   }
   return(sum);
}

template<class Y,class Z>
Z VnbX<Y,Z>::ScoPairM(Y mi,Y mj){
   Y i,j,k;
   Z *flw,xx,sum=0;
   this->gopen_db_map();
   Load(mi);
   k=this->nw;
   Y *nwk=this->nwd;
   Z *lwk=this->lwt;
   Load(mj);
   i=j=0;
   while((i<k)&&(j<this->nw)){
      if(nwk[i]<this->nwd[j])i++;
      else if(this->nwd[j]<nwk[i])j++;
      else {
         if(mrk[nwk[i]])sum+=wx[nwk[i]];
         i++;j++;
      }
   }
   return(sum);
}

template<class Y,class Z> 
void VnbX<Y,Z>::ScoreAll(float (*global)(long)){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;
  
   for(i=0;i<this->nw;i++){
      pInd=this->readp(this->nwd[i]);
      flw=this->readz(this->nwd[i]);
      j=pInd->ix;
      xx=(Z)global(j)*this->lwt[i];
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k];
      }
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAllWg(Ordr<Y,Z> *pOrd){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   if(sco==NULL)sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nw;i++){
      if(mrk[this->nwd[i]]==1){
         pInd=this->readp(this->nwd[i]);
         flw=this->readz(this->nwd[i]);
         j=pInd->ix;
         xx=weightG(pOrd,pInd)*this->lwt[i];
         for(k=0;k<j;k++){
            sco[pInd->idx[k]]+=xx*flw[k];
         }
      }
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAllmrk(float (*global)(long)){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco==NULL)sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;
cout << "Step 1" << endl;

   for(i=0;i<this->nw;i++){
cout << "Step 2" << endl;
      if(mrk[this->nwd[i]]==1){
cout << "Step 3" << endl;
         pInd=this->readp(this->nwd[i]);
         flw=this->readz(this->nwd[i]);
         j=pInd->ix;
         xx=(Z)global(j)*this->lwt[i];
         for(k=0;k<j;k++){
            sco[pInd->idx[k]]+=xx*flw[k];
         }
      }
cout << "Step 4" << endl;
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAllmrkLoc(float (*global)(long),float (*d_local)(int,long),int *lnt){
   Y i,j,k,mi,jx;
   int u;
   long lmn;
   Z *flw,xx,xmn=0,yy;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nw;i++){
      if(mi=mrk[this->nwd[i]]){
         xmn+=this->lwt[i];
      }
   }
   lmn=rnd((double)xmn);
   for(i=0;i<this->nw;i++){
      if(mi=mrk[this->nwd[i]]){
         u=(int)rnd((double)this->lwt[i]);
         yy=d_local(u,lmn);
         pInd=this->readp(this->nwd[i]);
         flw=this->readz(this->nwd[i]);
         j=pInd->ix;
         xx=(Z)global(j)*yy;
         for(k=0;k<j;k++){
            jx=pInd->idx[k];
            sco[jx]+=xx*d_local((int)rnd((double)flw[k]),(long)lnt[jx]);
         }
      }
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAllmrkLearned(float (*global)(long),int *lnt,double **wx){
   Y i,j,k,ui,uj,uk;
   Z *flw,xx,xmn=0;
   int lmn,u;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nw;i++){
      if(mrk[this->nwd[i]]){
         xmn+=this->lwt[i];
      }
   }
   lmn=(int)rnd((double)xmn);
   for(i=0;i<this->nw;i++){
      if(mrk[this->nwd[i]]){
         u=(int)rnd((double)this->lwt[i]);
         pInd=this->readp(this->nwd[i]);
         flw=this->readz(this->nwd[i]);
         j=pInd->ix;
         xx=(Z)global(j)*wx[lmn][u];
         for(k=0;k<j;k++){
            ui=pInd->idx[k];
            uj=rnd(flw[k]);
            if(uj>9)uj=9;
            uk=lnt[ui];
            if(uk>499)uk=499;
            sco[ui]+=xx*wx[uk][uj];
         }
      }
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAllmrkLearned(Z *wg,int *lnt,double **wx){
   Y i,j,k,ui,uj,uk;
   Z *flw,xx,xmn=0;
   int lmn,u;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nw;i++){
      if(mrk[this->nwd[i]]){
         xmn+=this->lwt[i];
      }
   }
   lmn=(int)rnd((double)xmn);
   for(i=0;i<this->nw;i++){
      if(mrk[this->nwd[i]]){
         u=(int)rnd((double)this->lwt[i]);
         pInd=this->readp(this->nwd[i]);
         flw=this->readz(this->nwd[i]);
         j=pInd->ix;
         xx=wg[this->nwd[i]]*wx[lmn][u];
         for(k=0;k<j;k++){
            ui=pInd->idx[k];
            uj=rnd(flw[k]);
            if(uj>9)uj=9;
            uk=lnt[ui];
            if(uk>499)uk=499;
            sco[ui]+=xx*wx[uk][uj];
         }
      }
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::ScoreAllEx(Z *wx){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nw;i++){
      pInd=this->readp(this->nwd[i]);
      flw=this->readz(this->nwd[i]);
      j=pInd->ix;
      xx=(Z)wx[this->nwd[i]]*this->lwt[i];
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k];
      }
   }
}

template<class Y,class Z>
Y VnbX<Y,Z>::ScoreAllExTop(Z *wx){
   Y i,j,k,ct;
   Z *flw,xx,zz,mx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   mx=0;ct=0;
   for(i=0;i<this->nw;i++){
      pInd=this->readp(this->nwd[i]);
      flw=this->readz(this->nwd[i]);
      j=pInd->ix;
      xx=(Z)wx[this->nwd[i]]*this->lwt[i];
      for(k=0;k<j;k++){
         zz=sco[pInd->idx[k]]+=xx*flw[k];
         if(mx<zz){
            mx=zz;
            ct=1;
         }
         else if(mx==zz)ct++;
      }
   }
   return(ct);
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAllExmrk(Z *wx){
   Y i,j,k;
   Z *flw,xx;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nw;i++){
      if(mrk[this->nwd[i]]){
         pInd=this->readp(this->nwd[i]);
         flw=this->readz(this->nwd[i]);
         j=pInd->ix;
         xx=(Z)wx[this->nwd[i]]*this->lwt[i];
         for(k=0;k<j;k++){
            sco[pInd->idx[k]]+=xx*flw[k];
         }
      }
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::ScoreAll(Z *sx,Z *wx){
   Y i,j,k;
   Z xx;
   Indx<Y> *pInd;
   this->gopen_map();

   for(i=0;i<this->nw;i++){
      pInd=this->readp(this->nwd[i]);
      j=pInd->ix;
      xx=(Z)wx[this->nwd[i]];
      for(k=0;k<j;k++){
         sx[pInd->idx[k]]+=xx;
      }
   }
   sco=sx;
}

template<class Y,class Z> 
void VnbX<Y,Z>::ScoreAll(Z *wx){
   Y i,j,k;
   Z xx;
   Indx<Y> *pInd;
   this->gopen_map();
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;

   for(i=0;i<this->nw;i++){
      pInd=this->readp(this->nwd[i]);
      j=pInd->ix;
      xx=(float)wx[this->nwd[i]];
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx;
      }
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::ScoreAll(BCount<Z> &Dc){
   Y i,j,k;
   Z xx;
   char *pch;
   Indx<Y> *pInd;
   this->gopen_map();
   this->gopen_hash();
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;
  
   Dc.node_first();
   while(Dc.node_next()){
      pch=Dc.show_str();
      if(i=this->find(pch)){
         xx=Dc.count();
         pInd=this->readp(i-1);
         j=pInd->ix;
         for(k=0;k<j;k++){
            sco[pInd->idx[k]]+=xx;
         }
      }
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAll(Ordr<Y,Z> *pOrd){
   Y i,j,k,m;
   Z xx,*flw;
   char *pch;
   Indx<Y> *pInd;
   this->gopen_map();
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;
   for(m=0;m<pOrd->ix;m++){
      i=pOrd->ind(m,xx);
      pInd=this->readp(i);
      flw=this->readz(i);
      j=pInd->ix;
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx*flw[k];
      }
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::ScoreAllNlw(Ordr<Y,Z> *pOrd){
   Y i,j,k,m;
   Z xx;
   char *pch;
   Indx<Y> *pInd;
   this->gopen_map();
   if(sco)delete [] sco;
   sco=new Z[this->ndoc];
   for(i=0;i<this->ndoc;i++)sco[i]=0;
   for(m=0;m<pOrd->ix;m++){
      i=pOrd->ind(m,xx);
      pInd=this->readp(i);
      j=pInd->ix;
      for(k=0;k<j;k++){
         sco[pInd->idx[k]]+=xx;
      }
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::gopen_Len(void){
   xlen=(Z*)this->get_Mmap("len");
}

template<class Y,class Z> 
void VnbX<Y,Z>::gclose_Len(void){
   this->dst_Mmap("len",(char*&)xlen);
}

template<class Y,class Z> 
void VnbX<Y,Z>::RevisCos(float (*global)(long)){
   Y i,j,k;
   Z *flw,xx,sum=0;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;

   for(i=0;i<this->nw;i++){
      j=(long)this->freq[this->nwd[i]];
      xx=(Z)global(j);
      sum+=xx*this->lwt[i]*this->lwt[i];
   }
   if(sum>0)sum=1.0/(Z)sqrt((double)sum);
   else sum=1.0;

   for(i=0;i<this->ndoc;i++){
      sco[i]*=(sum*xlen[i]);
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::RevisCos(Z *wx){
   Y i,j,k;
   Z *flw,xx,sum=0;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;

   for(i=0;i<this->nw;i++){
      xx=(Z)wx[this->nwd[i]];
      sum+=xx*this->lwt[i]*this->lwt[i];
   }
   if(sum>0)sum=1.0/(Z)sqrt((double)sum);
   else sum=1.0;

   for(i=0;i<this->ndoc;i++){
      sco[i]*=(sum*xlen[i]);
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::RevisX(Z *wx){
   Y i,j,k;
   Z *flw,xx,yy,zz,sum=0;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;

   for(i=0;i<this->nw;i++){
      xx=(Z)wx[this->nwd[i]];
      sum+=xx*this->lwt[i]*this->lwt[i];
   }
   if(sum>0)sum=1.0/sum;
   else sum=1.0;

   for(i=0;i<this->ndoc;i++){
      xx=sco[i];
      yy=xx*sum;
      zz=xx*xlen[i]*xlen[i];
      sco[i]=(yy<zz)?yy:zz;
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::RevisProb1(float (*global)(long)){
   Y i,k;
   long j;
   Z *flw,xx,sum=0;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;
   
   for(i=0;i<this->nw;i++){
      j=(long)this->freq[this->nwd[i]];
      xx=(Z)global(j);
      sum+=xx*this->lwt[i];
   }
   if(sum>0)sum=0.5/sum;
   else sum=1.0;

   for(i=0;i<this->ndoc;i++){
      sco[i]*=sum+xlen[i];
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::RevisProb1(Z *wx){
   Y i,j,k;
   Z *flw,xx,sum=0;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;

   for(i=0;i<this->nw;i++){
      xx=(Z)wx[this->nwd[i]];
      sum+=xx*this->lwt[i];
   }
   if(sum>0)sum=0.5/sum;
   else sum=1.0;

   for(i=0;i<this->ndoc;i++){
      sco[i]*=sum+xlen[i];
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::RevisProb2(float (*global)(long)){
   Y i,k;
   long j;
   Z *flw,xx,sum=0;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;

   for(i=0;i<this->nw;i++){
      j=(long)this->freq[this->nwd[i]];
      xx=(Z)global(j);
      sum+=xx*this->lwt[i];
   }
   if(sum>0)sum=0.5*sum;
   else sum=1.0;

   for(i=0;i<this->ndoc;i++){
      sco[i]/=sum+xlen[i];
   }
}

template<class Y,class Z> 
void VnbX<Y,Z>::RevisProb3(float (*global)(long)){
   Y i,k;
   long j;
   Z *flw,xx,sum=0;
   Indx<Y> *pInd;
   this->gopen_map();
   xttrc=(double)this->ndoc;

   for(i=0;i<this->nw;i++){
      j=(long)this->freq[this->nwd[i]];
      xx=(Z)global(j);
      sum+=xx*this->lwt[i];
   }
   if(sum==0)sum=1.0;

   for(i=0;i<this->ndoc;i++){
      sco[i]/=(Z)sqrt((double)sum*xlen[i]);
   }
}


template<class Y,class Z> 
void VnbX<Y,Z>::RevisProd(Z zz){
   Y i,k;
   long j;
   Z *flw,xx,sum=0;

   for(i=0;i<this->ndoc;i++){
      if(xlen[i]>=zz)sco[i]*=xlen[i];
      else sco[i]*=zz;
   }
}

template<class Y,class Z> 
Ordr<Y,Z> *VnbX<Y,Z>::Skim(Y n){
   if(sco==NULL)return(NULL);
   Z sx;
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(n,this->ndoc,sco);

   Y i=0;
   pOrd->ind(i,sx);
   if(sx>0)return(pOrd);
   else return(NULL);
}

template<class Y,class Z> 
Ordr<Y,Z> *VnbX<Y,Z>::Skimgr(Z s){
   if(sco==NULL)return(NULL);
   Y *dx=new Y[this->ndoc];
   Y i,m=0;
   for(i=0;i<this->ndoc;i++){
      if(sco[i]>s){
         dx[m++]=i;
      }
   }
   if(m){
      Indx<Y> *pnd=new Indx<Y>(m,dx,0);
      Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(m,pnd,sco);
      delete pnd;
      delete [] dx;
      return(pOrd);
   } 
   else {
      delete [] dx;
      return(NULL);
   }
}

template<class Y,class Z> 
Ordr<Y,Z> *VnbX<Y,Z>::Skimge(Z s){
   if(sco==NULL)return(NULL);
   Y *dx=new Y[this->ndoc];
   Y i,m=0;
   for(i=0;i<this->ndoc;i++){
      if(sco[i]>=s){
         dx[m++]=i;
      }
   }
   if(m){
      Indx<Y> *pnd=new Indx<Y>(m,dx,0);
      Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(m,pnd,sco);
      delete pnd;
      delete [] dx;
      return(pOrd);
   } 
   else {
      delete [] dx;
      return(NULL);
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::Set_All(int n){
   Y i;
   if(mrk!=NULL)delete [] mrk;
   mrk=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++)mrk[i]=n;
}

template<class Y,class Z>
void VnbX<Y,Z>::Set_Char(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strchr(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::Set_String(const char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(strstr(this->show(i),str))mrk[i]=n;
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::Set_NChar(char c,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strchr(this->show(i),c))mrk[i]=n;
   }
}

template<class Y,class Z>
void VnbX<Y,Z>::Set_NString(const char *str,int n){
   Y i;
   this->gopen_lexos();
   for(i=0;i<this->nwrd;i++){
      if(!strstr(this->show(i),str))mrk[i]=n;
   }
}

}
#endif
