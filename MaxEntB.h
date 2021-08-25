#ifndef MAXB_H
#define MAXB_H

#include <fstream>
#include <iostream>
#include <runn.h>
#include <Vnab.h>
#include <XPost.h>
#include <DataObj.h>
#include <CMark.h>
#include <ap.h>
#include <lbfgs.h>

using namespace std;

void *pMxe;
double fold;
long iter1,iter2;

namespace iret {

template<class Y,class Z>
class MaxEntBX : public CMark<Y,Z> {
   public:
      MaxEntBX(const char *nmspost,const char *nammxe); //Assumes the 
         //XPost is already made 
         //nammxe is the name used for maxent files.
      MaxEntBX(const char *nmspost,const char *pnam,const char *nammxe); 
         //Assumes the XPost is already made 
         //pnam is used in "path_pnam" where to find path, or if 1st char
         //is ':' then remainder is the path
         //nammxe is the name used for maxent files.
      ~MaxEntBX(void);
      //Functions
      void Set_Params(double gm); //Call this function to set parameters
         //gamma is the reciprocal of the variance for the Guassian prior of 
         //the weights for the active features. value=1.0E-4 good on Rebase
      void Init_System(Y m,Indx<Y> **mx); //m is number of classes
         //mx defines the classes. For class i, mx[i] is the index object
         //that contains the numbers of the docs of that class in the training
         //data. The part of the docset not contained in some mx[i] is what is
         //left over as the test set.
      void Set_Probs_EqCls(void); //Sets the total probability of each class
         //at 1/mc. Docs with a class have equal probs
      void Set_Probs_EqDoc(void); //Sets all doc probs to 1/mdoc
      void Create_pto(void); //Creates the 'o' file, mc*nwrd size
         //Also opens it for reading

      Y All_Feature(Y llm); //considers all features that occur for active list
         //llm is lower limit on frequency to be included in the active features
      Y All_FeatureMrk(Y llm); //considers all features that occur for active list
         //llm is lower limit on frequency to be included in the active features
         //only this->mrk[ui]=1 features can be included.
      Y Set_Feature(Indx<Y> *pTrm); //considers all features listed by pTrm as
         //candidates for active list
   
      void Optimize_Active(void); //Optimizes the lambdas over the active list
         //Calls unmap_ptt and Remake_ptt once each time it is called
         //new set of ptt and lambdas
      void Grad(double &fx,double* &gx); //computes function values and gradient
      void Init_ScoringI(Y m); //For internal scoring
         //m is number of classes
         //Sets up all that is needed for scoring. Must call Load function
         //to read in a list of active features that has been previously Saved
         //or must have trained weights by a call to Optimize_Active.
         //Sets the array invs ready for rapid scoring
      double *ScoreAll(void); //Scores all the docs in XPost; uses Value(i,px);
      double *ScoreSet(Indx<Y> *pInd); //Scores the docs listed (must be in XPost).
      void Value(Y i,double *px); //computes for doc number i in XPost set
         //px points at buffer with enough space for mc Zs
         //Uses RawValue(i,px).
      void RawValue(Y i,double *px); //computes for doc number i in XPost set
         //px points at buffer with enough space for mc Zs
         //Produces raw values not exponentiated or normalized.
      double Log_Likelihood(void); //Computes the log likelihood of the training 
         //data based on the computed weights. Only can be run when Init_ScoringI
         //has first been called because must access docs in XPost.
      void Init_ScoringE(Y m); //For external scoring
         //m is number of classes
         //Sets up all that is needed for scoring. Must call Load function
         //to read in a list of active features that has been previously Saved
         //or must have trained weights by a call to Optimize_Active.
         //Sets the array invs ready for rapid scoring
      void Value(Doc<Z> &Doc,double *px); //Uses Doc 
         //and computes and px points at buffer with enough space for mc Zs 
         //which are the p(y|x)
      void RawValue(Doc<Z> &Doc,double *px); //Uses Doc 
         //and computes and px points at buffer with enough space for mc Zs 
         //which are the raw scores
      void Value(List &Lt,double *px); //Uses List 
         //and computes and px points at buffer with enough space for mc Zs 
         //which are the raw scores
      void RawValue(List &Lt,double *px); //Uses List 
         //and computes and px points at buffer with enough space for mc Zs 
         //which are the raw scores

      //Functions to Save and Load active feature set and learned lambdas
      void Save(void);
      void Load(void);

      //Data
      Y mc; //number of classes (y-values)
      Y *markd; //Class marker values, one for each doc in training set,
         //class values 0 to mc-1, mc used for what is left (test set)
         //Begins with first doc in set and ends one beyond last doc.
      Y mdoc; //Number of documents in the training set, marked 0 to mc-1
      double *qdoc; //Probs of docs in classes, mc probabilites
      double *pto; //Observed expectation of each feature in training data, 
         //mc*(#terms). This array just to hold mc numbers for a single term
      ifstream *pfto; //To allow access to 'o' file
      double *ptt; //Current state of solution, mc*(#terms)
      Y ta; //Number of active features
      Y *txa; //Term numbers of the active features
         //In order for efficiency
      Y *tca; //Class numbers of the active features
      Y *invs; //Inverts the txa array to speed scoring
      double *lam; //Array lambdas, one for each active feature
      double *ps; //Working space to compute gradient
      double *pg; //Working space to compute gradient
      double *Zx; //Working space to compute gradient
      double xu; //exponents term linear lam
      double yu; //regularization term quadratic in lam
      double zu; //log term
      double *po; //Working space to compute gradient
      double gam; //Regularization Parameter for normal prior. 
      double *sco; //Score array
      FBase *pFbme; //Pointer at FBase object to carry MaxEnt produced 
         //weights.
};

template<class Y,class Z>
MaxEntBX<Y,Z>::MaxEntBX(const char *nmpost,const char *nammxe) : CMark<Y,Z>(nmpost){
   pFbme=new FBase("maxent",nammxe);
   pfto=NULL;
   qdoc=NULL;
   markd=NULL;
   txa=NULL;
   tca=NULL;
   lam=NULL;
   invs=NULL;
   ps=NULL;
   Zx=NULL;
   po=NULL;
   sco=NULL;
}

template<class Y,class Z>
MaxEntBX<Y,Z>::MaxEntBX(const char *nmpost,const char *pnam,const char *nammxe) : CMark<Y,Z>(nmpost,pnam){
   pFbme=new FBase("maxent",nammxe);
   pfto=NULL;
   qdoc=NULL;
   markd=NULL;
   txa=NULL;
   tca=NULL;
   lam=NULL;
   invs=NULL;
   ps=NULL;
   Zx=NULL;
   po=NULL;
   sco=NULL;
}

template<class Y,class Z>
MaxEntBX<Y,Z>::~MaxEntBX(void){
   if(pfto)pFbme->dst_Istr(pfto);
   delete pFbme;
   if(qdoc)delete [] qdoc;
   if(markd)delete [] markd;
   if(txa)delete [] txa;
   if(tca)delete [] tca;
   if(lam)delete [] lam;
   if(invs)delete [] invs;
   if(sco)delete [] sco;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Init_System(Y m,Indx<Y> **mx){
   Y i,j,u,*pm;

   this->gopen_map();

   mc=m;
   mdoc=0;
   for(i=0;i<m;i++)mdoc+=mx[i]->ix;
   if(this->ndoc<mdoc){cout << "Error in size of training set" << endl;exit(0);}  

   if(markd) delete [] markd;
   markd=new Y[this->ndoc];
   for(i=0;i<this->ndoc;i++)markd[i]=mc;
   for(i=0;i<m;i++){
      u=mx[i]->ix;
      pm=mx[i]->idx;
      for(j=0;j<u;j++)markd[pm[j]]=i;
   }

   pto=new Z[mc];

   ta=0;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Set_Params(double gm){
   gam=gm;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Set_Probs_EqDoc(void){
   Y i;
   double xx=1.0/((double)mdoc);
   if(qdoc) delete [] qdoc;
   qdoc=new double[mc];
   for(i=0;i<mc;i++)qdoc[i]=xx;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Set_Probs_EqCls(void){
   Y i,j,*tt;

   tt=new Y[mc];
   if(qdoc) delete [] qdoc;
   qdoc=new double[mc];
   for(i=0;i<mc;i++)tt[i]=0;
   for(i=0;i<this->ndoc;i++){
      if((j=markd[i])<mc)tt[j]++;
   }
   j=0;
   for(i=0;i<mc;i++)j+=tt[i];
   if(mdoc!=j){cout << "Error in document count!" << endl;exit(0);}
   for(i=0;i<mc;i++)qdoc[i]=1.0/((double)tt[i]*mc);
   delete [] tt;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Create_pto(void){
   Y i,j,k,m,fx;
   double xx,zz,*vp;
   Indx<Y> *pInd;

   if(pfto)pFbme->dst_Istr(pfto);
   vp=new double[mc];
   ofstream *pfout=pFbme->get_Ostr("o",ios::out);
   for(i=0;i<this->nwrd;i++){
      for(j=0;j<mc;j++)vp[j]=0;
      pInd=this->readp(i);
      for(k=0;k<pInd->ix;k++){
         m=markd[pInd->idx[k]];
         if(m<mc)(vp[m])+=qdoc[m];
      }
      pfout->write((char*)vp,mc*sizeof(double));
      this->mark(i,1000,"terms_pto");
   }
   delete [] vp;
   pFbme->dst_Ostr(pfout);
   pfto=pFbme->get_Istr("o",ios::in);
}

template<class Y,class Z>
Y MaxEntBX<Y,Z>::All_Feature(Y llm){
   Y i,j,k,ui,ia;

   pfto->seekg(0,ios::beg);
   ia=0;
   for(ui=0;ui<this->nwrd;ui++){
      pfto->read((char*)pto,mc*sizeof(Z));
      if(this->freq[ui]<llm) continue;
      for(i=0;i<mc;i++){
         if(pto[i])ia++;
      }
      this->mark(ui,100000,"count terms");
   }
   if(tca)delete [] tca;
   if(txa)delete [] txa;
   tca=new Y[ia];
   txa=new Y[ia+1];
   txa[ia]=-1;

   pfto->seekg(0,ios::beg);
   ia=0;
   for(ui=0;ui<this->nwrd;ui++){
      pfto->read((char*)pto,mc*sizeof(Z));
      if(this->freq[ui]<llm) continue;
      for(i=0;i<mc;i++){
         if(pto[i]){
            txa[ia]=ui;
            tca[ia]=i;
            ia++;
         }
      }
      this->mark(ui,100000,"terms");
   }
   ta=ia;
   return(ta);
} 

template<class Y,class Z>
Y MaxEntBX<Y,Z>::All_FeatureMrk(Y llm){
   Y i,j,k,ui,ia;

   pfto->seekg(0,ios::beg);
   ia=0;
   for(ui=0;ui<this->nwrd;ui++){
      pfto->read((char*)pto,mc*sizeof(double));
      if((!this->mrk[ui])||(this->freq[ui]<llm)) continue;
      for(i=0;i<mc;i++){
         if(pto[i])ia++;
      }
      this->mark(ui,100000,"count terms");
   }
   if(tca)delete [] tca;
   if(txa)delete [] txa;
   tca=new Y[ia];
   txa=new Y[ia+1];
   txa[ia]=-1;

   pfto->seekg(0,ios::beg);
   ia=0;
   for(ui=0;ui<this->nwrd;ui++){
      pfto->read((char*)pto,mc*sizeof(double));
      if((!this->mrk[ui])||(this->freq[ui]<llm)) continue;
      for(i=0;i<mc;i++){
         if(pto[i]){
            txa[ia]=ui;
            tca[ia]=i;
            ia++;
         }
      }
      this->mark(ui,100000,"terms");
   }
   ta=ia;
   return(ta);
}

template<class Y,class Z>
Y MaxEntBX<Y,Z>::Set_Feature(Indx<Y> *pTrm){
   Y i,j,k,ui,ia;

   pfto->seekg(0,ios::beg);
   ia=0;
   ui=-1;
   for(j=0;j<pTrm->ix;j++){
      k=pTrm->idx[j];
      while(++ui<=k)pfto->read((char*)pto,mc*sizeof(double));
      for(i=0;i<mc;i++){
         if(pto[i])ia++;
      }
      this->mark(j,10000,"count terms");
   }
   if(tca)delete [] tca;
   if(txa)delete [] txa;
   tca=new Y[ia];
   txa=new Y[ia+1];
   txa[ia]=-1;

   pfto->seekg(0,ios::beg);
   ia=0;
   ui=-1;
   for(j=0;j<pTrm->ix;j++){
      k=pTrm->idx[j];
      while(++ui<=k)pfto->read((char*)pto,mc*sizeof(double));
      for(i=0;i<mc;i++){
         if(pto[i]){
            txa[ia]=ui;
            tca[ia]=i;
            ia++;
         }
      }
      this->mark(ui,100000,"terms");
   }
   ta=ia;
   return(ta);
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Optimize_Active(void){
   Y i,j,k;
   
   ps=new double[mc*this->ndoc];
   Zx=new double[this->ndoc];
   po=new double[ta];
   lam=new double[ta];
   ap::real_1d_array x;
   for(i=0;i<ta;i++){
      pfto->seekg(txa[i]*mc*sizeof(double),ios::beg);
      pfto->read((char*)pto,mc*sizeof(double));
      po[i]=pto[tca[i]];
      lam[i]=0.0;
   }
   x.setspace(1,ta,lam);
   iter1=iter2=0;
   fold=1.0E100;
   const double epsg=1.0E-5;
   const double epsf=0.0001;
   const double epsx=0.0000000000000001;
   int info;
   const int maxits=0;
   pMxe=this;

   lbfgsminimize(ta,5,x,epsg,epsf,epsx,maxits,info);
  
   x.setspace(1,ta,NULL);
   delete [] po;
   delete [] ps;
   delete [] Zx;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Grad(double &fu,double* &g){
   Y i,j,k,u,v,u1;
   double xx,yy,zz,max=1.0;
   double lmm,del;
   Indx<Y> *pInd;
   
   for(i=0;i<mc*this->ndoc;i++) ps[i]=0;
   xu=0.0;
   yu=0.0;
   for(i=0;i<ta;i++){
      yu+=lam[i]*lam[i];
      u=tca[i];
      pInd=readp(txa[i]);
      xu+=lam[i]*po[i];
      for(j=0;j<pInd->ix;j++){
         k=pInd->idx[j];
         ps[k*mc+u]+=lam[i];
      }
   }
   
   zu=0.0;
   for(k=0;k<this->ndoc;k++){
      if((u1=markd[k])<mc){
         Zx[k]=0;
         for(i=0;i<mc;i++){
            Zx[k]+=exp(ps[mc*k+i]);
         }
         zu+=log(Zx[k])*qdoc[u1];
      }
   }
   fu=zu-xu+0.5*gam*yu;
      
   for(i=0;i<ta;i++){
      u=tca[i];
      pInd=this->readp(txa[i]);
      g[i]=-po[i];
      for(j=0;j<pInd->ix;j++){
        k=pInd->idx[j];
        if((u1=markd[k])<mc) g[i]+=(exp(ps[mc*k+u])/Zx[k])*qdoc[u1];
      }
      g[i]+=gam*lam[i];
   }
}

template<class Y,class Z>
double MaxEntBX<Y,Z>::Log_Likelihood(void){
   Y i,j,u;
   double sum=0,xx,*px;

   px=new double[mc];
   for(i=0;i<this->ndoc;i++){
      if((j=markd[i])<mc){
         Value(i,px);
         sum+=log(px[j])*qdoc[j];
      }
   }
   return(sum);
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Init_ScoringE(Y m){
   Y i,j,k;

   mc=m;
   this->gopen_map();
   this->gopen_hash();

   //Set up for rapid scoring
   if(invs)delete [] invs;
   invs=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++)invs[i]=0;

   j=0;
   i=-1;
   while(j<ta){
      if((k=txa[j])!=i){
         invs[k]=j+1;
         i=k;
      }
      j++;
   }
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Init_ScoringI(Y m){
   Y i,j,k;

   mc=m;
   this->gopen_db_map();

   //Set up for rapid scoring
   if(invs)delete [] invs;
   invs=new Y[this->nwrd];
   for(i=0;i<this->nwrd;i++)invs[i]=0;

   j=0;
   i=-1;
   while(j<ta){
      if((k=txa[j])!=i){
         invs[k]=j+1;
         i=k;
      }
      j++;
   }
}

template<class Y,class Z>
double *MaxEntBX<Y,Z>::ScoreAll(void){
   Y i;
   double px[mc];

   if(sco)delete [] sco;
   sco=new double[this->ndoc];

   for(i=0;i<this->ndoc;i++){
      Value(i,px);
      sco[i]=px[0];
      this->mark(i,1000,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
double *MaxEntBX<Y,Z>::ScoreSet(Indx<Y> *pInd){
   Y i,j;
   double px[mc];

   if(sco)delete [] sco;
   sco=new double[this->ndoc];

   for(j=0;j<pInd->ix;j++){
      i=pInd->idx[j];
      Value(i,px);
      sco[i]=px[0];
      this->mark(j,1000,"docs scored");
   }
   return(sco);
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Value(Y n,double *px){
   Y i,j,k,nwt,nt;
   double zz,sum=0;

   RawValue(n,px);
   for(i=0;i<mc;i++){
      zz=exp(px[i]);
      sum+=zz;
      px[i]=zz;
   }
   for(i=0;i<mc;i++)px[i]=px[i]/sum;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::RawValue(Y n,double *px){
   Y i,j,k,nwt,nt;
   double zz,sum=0;

   this->readp_db(n);
   for(i=0;i<mc;i++)px[i]=0;
   for(j=0;j<this->nw;j++){
      nwt=this->nwd[j];
      if((i=invs[nwt])){
         k=i-1;
         do {
            px[tca[k]]+=lam[k];
            k++;
         } while(txa[k]==nwt);
      }
   }
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Value(Doc<Z> &Doc,double *px){
   Y i,j,k,nwt,nt;
   double zz,sum=0;

   RawValue(Doc,px);
   for(i=0;i<mc;i++){
      zz=exp(px[i]);
      sum+=zz;
      px[i]=zz;
   }
   for(i=0;i<mc;i++)px[i]=px[i]/sum;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::RawValue(Doc<Z> &Doc,double *px){
   Y i,j,k,nwt,nt;
   double zz,sum=0;

   for(i=0;i<mc;i++)px[i]=0;
   nt=Doc.num_wrds();
   for(j=0;j<nt;j++){
      if((nwt=find(Doc.word[j]))&&((i=invs[nwt-1]))){
         k=i-1;
         do {
            px[tca[k]]+=lam[k];
            k++;
         } while(txa[k]==nwt-1);
      }
   }
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Value(List &Lt,double *px){
   Y i,j,k,nwt,nt;
   double zz,sum=0;

   RawValue(Lt,px);
   for(i=0;i<mc;i++){
      zz=exp(px[i]);
      sum+=zz;
      px[i]=zz;
   }
   for(i=0;i<mc;i++)px[i]=px[i]/sum;
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::RawValue(List &Lt,double *px){
   Y i,j,k,nwt,nt;
   double zz,sum=0;

   for(i=0;i<mc;i++)px[i]=0;
   Lt.node_first();
   while(Lt.node_next()){
      if((nwt=this->find(Lt.show_str()))&&((i=invs[nwt-1]))){
         k=i-1;
         do {
            px[tca[k]]+=lam[k];
            k++;
         } while(txa[k]==nwt-1);
      }
   }
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Save(void){
   Y i;
   ofstream *pfout=pFbme->get_Ostr("a",ios::out);
   *pfout << ta << endl;
   for(i=0;i<ta;i++){
      *pfout << txa[i] << " " << tca[i] << " " << lam[i] << endl;
   }
   pFbme->dst_Ostr(pfout);
}

template<class Y,class Z>
void MaxEntBX<Y,Z>::Load(void){
   Y i,j,k;
   ifstream *pfin=pFbme->get_Istr("a",ios::in);
   *pfin >> ta;

   if(txa)delete [] txa;
   if(tca)delete [] tca;
   if(lam)delete [] lam;
   txa=new Y[ta+1];
   tca=new Y[ta];
   lam=new double[ta];

   for(i=0;i<ta;i++){
      *pfin >> txa[i] >> tca[i] >> lam[i];
   }
   txa[ta]=-1;
   pFbme->dst_Istr(pfin);
}

}


void funcgrad(ap::real_1d_array qm,double &fz,ap::real_1d_array &qp){
   long i,j,k,u,tt,ip,iq,w;

   double *dp=qp.getcontent();
   ((iret::MaxEntBX<long,double> *)pMxe)->Grad(fz,dp);
   cout << "iter1 " << ++iter1 << endl;
   if(fz<fold){
      fold=fz;
      cout << "iter2 " << ++iter2 <<  " fz " << fz << endl;
   }
}

#endif

