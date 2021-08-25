#ifndef GBOOST_H
#define GBOOST_H

#include <iostream>
#include <fstream>
#include <DataObj.h>
#include <Isgrid.h>
#define EPS 1.0E-5
using namespace std;
namespace iret {

template <class Y,class Z> 
class GBoostX {
public:
  GBoostX(Y trd,Y tsd,Z eps); //trd # training docs
     //tsd # of test docs, eps the limit on p, eps <=p<=1-eps
  ~GBoostX();
  void init(Indx<Y> *gind,Y gr=10000); //Initializes the pdoc array
     //Also sets the mark array for the good using gind
     //Also sets the granularity for Isgrid use.
  void update(Z *trs,Z *tss); //trs is poYer at
     //training score array from some function h & tss is
     //the same poYer at testing score array by same h.
  void update_store(Z *trs,const char *nam); //trs is poYer at
     //at training score array. nam is name under which to store
     //the Isgrid results.

  //Data
  Y tr_doc; //number of training docs
  Y *mark; //Array that marks the good with 1, bad 0
  Z *pdoc; //Probability array for the training docs

  Y ts_doc; //number of testing docs
  Z *ts_sxx; //Final output score for testing docs

  Y grn; //Granularity, default 10,000.
  Y t; //Number of the iteration.
  Z zprod; //Product of Zs
  Z epsilon; //Limit on p.
  Isgrid *pIsg; //PoYer at Isgrid object last used in update.
};

template <class Y,class Z> 
class GBoostX2 {
public:
  GBoostX2(Y nd,Indx<Y> *gd,Indx<Y> *bd,Indx<Y> *ts); //ndoc whole space
     //good docs, bad docs, test docs
  ~GBoostX2();
  void init(Z eps,Y gr=10000); //Initializes the pdoc 
     // and ts_sxx arrays
     //eps the limit on p, eps <=p<=1-eps
  void update(Z *sco); 
     //score array from some function h, over whole space
  void update_store(Z *sco, const char *nam); //sco is poYer at
     //score array. nam is name under which to store
     //the Isgrid results.

  //Data
  Z *pdoc; //Probability array for the training docs

  Y ndoc; //number of docs in whole space
  Y tdoc; //number of docs in training set
  Z *ts_sxx; //Final output score for testing docs
     //array covers all of ndoc, same as sco in update

  Indx<Y> *gdd; //Good set
  Indx<Y> *bdd; //Bad set
  Indx<Y> *tst; //Test set

  Y grn; //Granularity, default 10,000.
  Y t; //Number of the iteration.
  Z zprod; //Product of Zs
  Z epsilon; //Limit on p.
  Isgrid *pIsg; //PoYer at Isgrid object last used in update.
};

template <class Y,class Z> 
class GBoostX3 {
public:
  GBoostX3(Y nd,Y t,Y *tn); //nd whole space
     //t num training docs, tn +1, -1, 0
  ~GBoostX3();
  void init(Z eps,Y gr=10000); //Initializes the pdoc 
     // and ts_sxx arrays
     //eps the limit on p, eps <=p<=1-eps
  void update(Z *sco);
     //score array from some function h, over whole space
  void update_store(Z *sco, const char *nam); //sco is poYer at
     //score array. nam is name under which to store
     //the Isgrid results.
  Z Del_entropy(Y n,Z *sco); //If n is not in training set
     //Computes the expected change in entropy on unseen docs
     //if it is labelled and added to the training set.
     //sco is the score produced by Vnak scoring of doc n
  Z Del_error(Y n,Z *sco); //If n is not in training set
     //Computes the expected change in error rate on unseen docs
     //if it is labelled and added to the training set.
     //sco is the score produced by Vnak scoring of doc n
  Z Gain_Ratio(Y n,Z *sco); //If n is not in training set
     //Computes the expected gain ratio
     //if it is labelled and added to the training set.
     //sco is the score produced by Vnak scoring of doc n

  //Data
  Z *pdoc; //Probability array for the training docs

  Y ndoc; //number of docs in whole space
  Y tdoc; //number of docs in training set
  Z *ts_sxx; //Final output score for testing docs
     //array covers all of ndoc, same as sco in update

  Y *trn; //Stores the training set +1, -1, 0

  Y grn; //Granularity, default 10,000.
  Y t; //Number of the iteration.
  Z lam; //Factor to scale down old scores
  Z zprod; //Product of Zs
  Z epsilon; //Limit on p.
  Isgrid *pIsg; //PoYer at Isgrid object last used in update.
};

template <class Y,class Z> 
class ABoostX {
public:
  ABoostX(Y trd,Y tsd); //trd # training docs
     //tsd # of test docs
  ~ABoostX();
  void init(Indx<Y> *gind); //Initializes the pdoc array
     //Also sets the mark array for the good using gind
  Z update(Z *trs,Z *tss); //trs is poYer at
     //training score array from some function h & tss is
     //the same poYer at testing score array by same h.
     //Returns the optimal alpha value used in update.
  Z update(Z *trs); //trs is poYer at
     //training score array from some function h 
     //Returns the optimal alpha value used in update.
  Z Z_alpha(Z alp); //computes the value for alp.

  //Data
  Y tr_doc; //number of training docs
  Y *mark; //Array that marks the good with 1, bad 0
  Z *tr_sco; //Array of training scores.
  Z *pdoc; //Probability array for the training docs

  Y ts_doc; //number of testing docs
  Z *ts_sxx; //Final output score for testing docs

  Y t; //Number of the iteration.
  Z zprod; //Product of Zs
};

template <class Y,class Z> 
class ABoostX2 {
public:
  ABoostX2(Y nd,Indx<Y> *gd,Indx<Y> *bd,Indx<Y> *ts); //ndoc whole space
     //good docs, bad docs, test docs
  ~ABoostX2();
  void init(void); //Initializes the pdoc
     // and ts_sxx arrays
  Z update(Z *sco);
     //score array from some function h, over whole space
     //Returns the optimal alpha value used in update.
  Z Z_alpha(Z alp); //computes the value for alp.
  Z Z_alphf(Z alp); //computes the value for alp.
     //Uses scf instead of sco

  //Data
  Z *pdoc; //Probability array for the training docs
  Z *sco; //Array of scores
  Z *scf; //Alternative array of scores

  Y ndoc; //number of docs in whole space
  Y tdoc; //number of docs in training set
  Z *ts_sxx; //Final output score for testing docs
     //array covers all of ndoc, same as sco in update

  Indx<Y> *gdd; //Good set
  Indx<Y> *bdd; //Bad set
  Indx<Y> *tst; //Test set

  Y t; //Number of the iteration.
  Z zprod; //Product of Zs
};

template <class Y,class Z> 
class ABoostX3 {
public:
  ABoostX3(Y nd,Indx<Y> *gd,Indx<Y> *bd,Indx<Y> *ts); //ndoc whole space
     //good docs, bad docs, test docs
  ~ABoostX3();
  void init(void); //Initializes the pdoc
     // and ts_sxx arrays

  Z update(Z *sco);
     //score array from some function h, over whole space
     //Returns the optimal alpha value used in update.
  Z find_alpha(Z alp,Z &hnv); //Calls Z_alpha and finds min. poY
  Z find_beta(Z bet,Z &hnv);  //Calls Z_beta and finds min. poY

  Z Z_alpha(Z alp); //computes the value for alp.
  Z Z_alphf(Z alp); //computes the value for alp.
     //Uses scf instead of sco
  Z Z_beta(Z bet); //computes the value for bet.
  Z Z_betf(Z bet); //computes the value for bet.
     //Uses scf instead of sco

  //Data
  Z *pdoc; //Probability array for the training docs
  Z *sco; //Array of scores
  Z *scf; //Alternative array of scores
  Z alpha;
  Z beta;

  Y ndoc; //number of docs in whole space
  Y tdoc; //number of docs in training set
  Z *ts_sxx; //Final output score for testing docs
     //array covers all of ndoc, same as sco in update

  Indx<Y> *gdd; //Good set
  Indx<Y> *bdd; //Bad set
  Indx<Y> *tst; //Test set

  Y t; //Number of the iteration.
  Z zprod; //Product of Zs
};

template <class Y,class Z> 
GBoostX<Y,Z>::GBoostX(Y trd,Y tsd,Z eps) {
   tr_doc=trd;
   ts_doc=tsd;
   epsilon=eps;
   mark=NULL;
   pdoc=NULL;
   ts_sxx=NULL;
   pIsg=NULL;
}

template <class Y,class Z> 
GBoostX<Y,Z>::~GBoostX() {
   if(mark!=NULL)delete [] mark;
   if(pdoc!=NULL)delete [] pdoc;
   if(ts_sxx!=NULL)delete [] ts_sxx;
   if(pIsg)delete pIsg;
}

template <class Y,class Z> 
void GBoostX<Y,Z>::init(Indx<Y> *gind,Y gr){
   Y i,j;
   Z xx=1.0/((Z)tr_doc);

   if(mark!=NULL)delete [] mark;
   mark=new Y[tr_doc];
   if(pdoc!=NULL)delete [] pdoc;
   pdoc=new Z[tr_doc];

   for(i=0;i<tr_doc;i++){
      mark[i]=0;
      pdoc[i]=xx;
   }
   for(i=0;i<gind->ix;i++)mark[*(gind->idx+i)]=1; 

   if(ts_sxx!=NULL)delete [] ts_sxx;
   ts_sxx=new Z[ts_doc];

   for(i=0;i<ts_doc;i++){
      ts_sxx[i]=0.0;
   }
   grn=gr;
   t=0;
   zprod=(Z)tr_doc;
}

template <class Y,class Z> 
void GBoostX<Y,Z>::update(Z *trs,Z *tss){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,sum,deps;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   dmax=dmin=trs[0];
   for(i=1;i<tr_doc;i++){
      xx=trs[i];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<tr_doc;i++){
      if(mark[i])pIsg->add_data(trs[i],pdoc[i],pdoc[i]);
      else pIsg->add_data(trs[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->extend_1df();
   if(pflag)cout << "Average p " << pIsg->avg() << " Information " << pIsg->info() << endl;

   //Update pdoc array
   deps=1.0-epsilon;
   sum=0;
   for(i=0;i<tr_doc;i++){
      xx=pIsg->val_1df(trs[i]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      if(mark[i])pdoc[i]*=sqrt((1.0-xx)/xx);
      else pdoc[i]*=sqrt(xx/(1.0-xx));
      sum+=pdoc[i];
   }
   zprod=zprod*sum;
   for(i=0;i<tr_doc;i++){
      pdoc[i]/=sum;
   }

   //Update tr_sxx
   for(i=0;i<ts_doc;i++){
      xx=pIsg->val_1df(tss[i]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps; 
      ts_sxx[i]+=0.5*log(xx/(1.0-xx));
   }
}
        
template <class Y,class Z> 
void GBoostX<Y,Z>::update_store(Z *trs,const char *nam){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,sum,deps;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   dmax=dmin=trs[0];
   for(i=1;i<tr_doc;i++){
      xx=trs[i];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<tr_doc;i++){
      if(mark[i])pIsg->add_data(trs[i],pdoc[i],pdoc[i]);
      else pIsg->add_data(trs[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->set_name(nam);
   pIsg->write_1df();

   //Update pdoc array
   pIsg->extend_1df();
   if(pflag)cout << "Average p " << pIsg->avg() << " Information " << pIsg->info() << endl;
   deps=1.0-epsilon;
   sum=0;
   for(i=0;i<tr_doc;i++){
      xx=pIsg->val_1df(trs[i]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      if(mark[i])pdoc[i]*=sqrt((1.0-xx)/xx);
      else pdoc[i]*=sqrt(xx/(1.0-xx));
      sum+=pdoc[i];
   }
   zprod=zprod*sum;
   for(i=0;i<tr_doc;i++){
      pdoc[i]/=sum;
   }
}

//GBoostX2

template <class Y,class Z> 
GBoostX2<Y,Z>::GBoostX2(Y nd,Indx<Y> *gd,Indx<Y> *bd,Indx<Y> *ts){
   ndoc=nd;
   pdoc=NULL;
   ts_sxx=NULL;
   gdd=gd;
   bdd=bd;
   tst=ts;
   tdoc=gdd->ix+bdd->ix;
   pIsg=NULL;
}

template <class Y,class Z> 
GBoostX2<Y,Z>::~GBoostX2() {
   if(pdoc!=NULL)delete [] pdoc;
   if(ts_sxx!=NULL)delete [] ts_sxx;
   if(pIsg)delete pIsg;
}

template <class Y,class Z> 
void GBoostX2<Y,Z>::init(Z eps,Y gr){
   Y i,j;
   epsilon=eps;
   grn=gr;

   Z xx=1.0/((Z)tdoc);

   if(pdoc!=NULL)delete [] pdoc;
   pdoc=new Z[ndoc];
   for(i=0;i<ndoc;i++){
      pdoc[i]=xx;
   }

   if(ts_sxx!=NULL)delete [] ts_sxx;
   ts_sxx=new Z[ndoc];

   for(i=0;i<ndoc;i++){
      ts_sxx[i]=0.0;
   }
   t=0;
   zprod=(Z)tdoc;
}

template <class Y,class Z> 
void GBoostX2<Y,Z>::update(Z *sco){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,sum,deps;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   j=gdd->idx[0];
   dmax=dmin=sco[j];
   for(i=1;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx=sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx=sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pIsg->add_data(sco[j],pdoc[j],pdoc[j]);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pIsg->add_data(sco[j],pdoc[j],0.0);
   }
   pIsg->dim1();

   //Update pdoc array
   pIsg->extend_1df();
   if(pflag)cout << "Average p " << pIsg->avg() << \
        " Information " << pIsg->info() << endl;
   deps=1.0-epsilon;
   sum=0;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx=pIsg->val_1df(sco[j]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      pdoc[j]*=sqrt((1.0-xx)/xx);
      sum+=pdoc[j];
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx=pIsg->val_1df(sco[j]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      pdoc[j]*=sqrt(xx/(1.0-xx));
      sum+=pdoc[j];
   }
   zprod=zprod*sum;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pdoc[j]/=sum;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pdoc[j]/=sum;
   }

   //Update tr_sxx
   for(i=0;i<tst->ix;i++){
      j=tst->idx[i];
      xx=pIsg->val_1df(sco[j]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      ts_sxx[j]+=0.5*log(xx/(1.0-xx));
   }
}


template <class Y,class Z> 
void GBoostX2<Y,Z>::update_store(Z *sco,const char *nam){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,sum,deps;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   j=gdd->idx[0];
   dmax=dmin=sco[j];
   for(i=1;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx=sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx=sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pIsg->add_data(sco[j],pdoc[j],pdoc[j]);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pIsg->add_data(sco[j],pdoc[j],0.0);
   }
   pIsg->dim1();
   pIsg->set_name(nam);
   pIsg->write_1df();

   //Update pdoc array
   pIsg->extend_1df();
   if(pflag)cout << "Average p " << pIsg->avg() << \
      " Information " << pIsg->info() << endl;
   deps=1.0-epsilon;
   sum=0;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx=pIsg->val_1df(sco[j]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      pdoc[j]*=sqrt((1.0-xx)/xx);
      sum+=pdoc[j];
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx=pIsg->val_1df(sco[j]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      pdoc[j]*=sqrt(xx/(1.0-xx));
      sum+=pdoc[j];
   }
   zprod=zprod*sum;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pdoc[j]/=sum;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pdoc[j]/=sum;
   }

   //Update tr_sxx
   for(i=0;i<tst->ix;i++){
      j=tst->idx[i];
      xx=pIsg->val_1df(sco[j]);
      if(xx<epsilon)xx=epsilon;
      if(xx>deps)xx=deps;
      ts_sxx[j]+=0.5*log(xx/(1.0-xx));
   }
}

//GBoostX3

template <class Y,class Z> 
GBoostX3<Y,Z>::GBoostX3(Y nd,Y t,Y *tn){
   ndoc=nd;
   pdoc=NULL;
   ts_sxx=NULL;
   tdoc=t;
   trn=tn;
   pIsg=NULL;
}

template <class Y,class Z> 
GBoostX3<Y,Z>::~GBoostX3() {
   if(pdoc!=NULL)delete [] pdoc;
   if(ts_sxx!=NULL)delete [] ts_sxx;
   if(pIsg)delete pIsg;
}

template <class Y,class Z> 
void GBoostX3<Y,Z>::init(Z eps,Y gr){
   Y i,j;
   epsilon=eps;
   grn=gr;

   Z xx=1.0/((Z)tdoc);

   if(pdoc!=NULL)delete [] pdoc;
   pdoc=new Z[ndoc];
   for(i=0;i<ndoc;i++){
      pdoc[i]=xx;
   }

   if(ts_sxx!=NULL)delete [] ts_sxx;
   ts_sxx=new Z[ndoc];

   for(i=0;i<ndoc;i++){
      ts_sxx[i]=0.0;
   }
   t=0;
   zprod=(Z)tdoc;
   lam=0.99;
}

template <class Y,class Z> 
void GBoostX3<Y,Z>::update(Z *sco){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,yy,sum,deps;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   i=0;
   while(!trn[i])i++;
   dmax=dmin=sco[i];
   for(i=0;i<ndoc;i++){
      if(trn[i]){
         xx=sco[i];
         dmax=xx>dmax?xx:dmax;
         dmin=xx<dmin?xx:dmin;
      }
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data(sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data(sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();

   //Update pdoc array
   pIsg->extend_1df();
   if(pflag)cout << "Average p " << pIsg->avg() << \
        " Information " << pIsg->info() << endl;
   deps=1.0-epsilon;
   sum=0;
   for(i=0;i<ndoc;i++){
      if(trn[i]==1){
         xx=pIsg->val_1df(sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         yy=(1.0-xx)/xx;
         ts_sxx[i]-=0.5*log(yy);
         pdoc[i]*=sqrt(yy);
         sum+=pdoc[i];
      }
      else if(trn[i]==-1){
         xx=pIsg->val_1df(sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         yy=xx/(1.0-xx);
         ts_sxx[i]+=0.5*log(yy);
         pdoc[i]*=sqrt(yy);
         sum+=pdoc[i];
      }
      else {
         xx=pIsg->val_1df(sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         ts_sxx[i]+=0.5*log(xx/(1.0-xx));
      }
   }
   zprod=zprod*sum;
   for(i=0;i<ndoc;i++){
      if(trn[i])pdoc[i]/=sum;
   }
}

template <class Y,class Z> 
void GBoostX3<Y,Z>::update_store(Z *sco,const char *nam){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,yy,sum,deps;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   i=0;
   while(!trn[i])i++;
   dmax=dmin=sco[i];
   for(i=0;i<ndoc;i++){
      if(trn[i]){
         xx=sco[i];
         dmax=xx>dmax?xx:dmax;
         dmin=xx<dmin?xx:dmin;
      }
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data(sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data(sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->set_name(nam);
   pIsg->write_1df();

   //Update pdoc array
   pIsg->extend_1df();
   if(pflag)cout << "Average p " << pIsg->avg() << \
        " Information " << pIsg->info() << endl;
   deps=1.0-epsilon;
   sum=0;
   for(i=0;i<ndoc;i++){
      if(trn[i]==1){
         xx=pIsg->val_1df(sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         yy=(1.0-xx)/xx;
         ts_sxx[i]-=0.5*log(yy);
         pdoc[i]*=sqrt(yy);
         sum+=pdoc[i];
      }
      else if(trn[i]==-1){
         xx=pIsg->val_1df(sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         yy=xx/(1.0-xx);
         ts_sxx[i]+=0.5*log(yy);
         pdoc[i]*=sqrt(yy);
         sum+=pdoc[i];
      }
      else {
         xx=pIsg->val_1df(sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         ts_sxx[i]+=0.5*log(xx/(1.0-xx));
      }
   }
   zprod=zprod*sum;
   for(i=0;i<ndoc;i++){
      if(trn[i])pdoc[i]/=sum;
   }
}

template <class Y,class Z> 
Z GBoostX3<Y,Z>::Del_entropy(Y n,Z *sco){
   int pflag=get_qflag();
   Y i,j,k;
   Z Ppos,Pneg,del=0,ps;
   Z dmax,dmin,xx,yy,deps,uu,pre,prf;
   Z sum1=0,sum2=0;

   xx=exp(ts_sxx[n]);
   yy=exp(-ts_sxx[n]);
   Ppos=xx/(xx+yy);
   Pneg=yy/(xx+yy);

   i=0;
   while(!trn[i])i++;
   dmax=dmin=(Z)sco[i];
   for(i=0;i<ndoc;i++){
      if(trn[i]){
         xx=(Z)sco[i];
         dmax=xx>dmax?xx:dmax;
         dmin=xx<dmin?xx:dmin;
      }
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data((Z)sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data((Z)sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->extend_1df();
   deps=1.0-epsilon;

   for(i=0;i<ndoc;i++){
      if(!trn[i]){
         //Prior
         uu=exp(-(Z)(2.0*ts_sxx[i]));
         ps=1.0/(1.0+uu);
         if(ps>0)pre=ps*log(ps);
         else pre=0;
         ps=1.0-ps;
         if(ps>0)pre+=ps*log(ps);
         sum1+=pre;
         //Updated
         xx=pIsg->val_1df((Z)sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         uu*=(1.0-xx)/xx;
         ps=1.0/(1.0+uu);
         if(ps>0)prf=ps*log(ps);
         else prf=0;
         ps=1.0-ps;
         if(ps>0)prf+=ps*log(ps);
         sum2+=prf;
      }
   }
   del+=Ppos*(sum2-sum1);

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(-dmax,-dmin);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data(-(Z)sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data(-(Z)sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->extend_1df();

   sum2=0;
   for(i=0;i<ndoc;i++){
      if(!trn[i]){
         //Prior
         uu=exp(-(Z)(2.0*ts_sxx[i]));
         //Updated
         xx=pIsg->val_1df(-(Z)sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         uu*=(1.0-xx)/xx;
         ps=1.0/(1.0+uu);
         if(ps>0)prf=ps*log(ps);
         else prf=0;
         ps=1.0-ps;
         if(ps>0)prf+=ps*log(ps);
         sum2+=prf;
      }
   }
   del+=Pneg*(sum2-sum1);
   
   return(del);
}

template <class Y,class Z> 
Z GBoostX3<Y,Z>::Del_error(Y n,Z *sco){
   int pflag=get_qflag();
   Y i,j,k;
   Z Ppos,Pneg,del=0,ps;
   Z dmax,dmin,xx,yy,deps,uu,pre,prf;
   Z sum1=0,sum2=0;

   xx=exp(ts_sxx[n]);
   yy=exp(-ts_sxx[n]);
   Ppos=xx/(xx+yy);
   Pneg=yy/(xx+yy);

   i=0;
   while(!trn[i])i++;
   dmax=dmin=(Z)sco[i];
   for(i=0;i<ndoc;i++){
      if(trn[i]){
         xx=(Z)sco[i];
         dmax=xx>dmax?xx:dmax;
         dmin=xx<dmin?xx:dmin;
      }
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data((Z)sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data((Z)sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->extend_1df();
   deps=1.0-epsilon;

   for(i=0;i<ndoc;i++){
      if(!trn[i]){
         //Prior
         uu=exp(-(Z)(2.0*ts_sxx[i]));
         ps=1.0/(1.0+uu);
         if(ps>0.5)pre=1.0-ps;
         else pre=ps;
         sum1+=pre;
         //Updated
         xx=pIsg->val_1df((Z)sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         uu*=(1.0-xx)/xx;
         ps=1.0/(1.0+uu);
         if(ps>0.5)prf=1.0-ps;
         else prf=ps;
         sum2+=prf;
      }
   }
   del+=Ppos*(sum1-sum2);

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(-dmax,-dmin);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data(-(Z)sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data(-(Z)sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->extend_1df();

   sum2=0;
   for(i=0;i<ndoc;i++){
      if(!trn[i]){
         //Prior
         uu=exp(-(Z)(2.0*ts_sxx[i]));
         //Updated
         xx=pIsg->val_1df(-(Z)sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         uu*=(1.0-xx)/xx;
         ps=1.0/(1.0+uu);
         if(ps>0.5)prf=1.0-ps;
         else prf=ps;
         sum2+=prf;
      }
   }
   del+=Pneg*(sum1-sum2);
   
   return(del);
}

template <class Y,class Z> 
Z GBoostX3<Y,Z>::Gain_Ratio(Y n,Z *sco){
   int pflag=get_qflag();
   Y i,j,k;
   Z Ppos,Pneg,del=0,ps;
   Z dmax,dmin,xx,yy,deps,uu,vv,prf;
   Z sum,ddx=0.01;

   xx=exp(ts_sxx[n]);
   yy=exp(-ts_sxx[n]);
   Ppos=xx/(xx+yy);
   Pneg=yy/(xx+yy);

   i=0;
   while(!trn[i])i++;
   dmax=dmin=(Z)sco[i];
   for(i=0;i<ndoc;i++){
      if(trn[i]){
         xx=(Z)sco[i];
         dmax=xx>dmax?xx:dmax;
         dmin=xx<dmin?xx:dmin;
      }
   }
   if(pflag)cout << "Min value: " << dmin << "; Max value: " << dmax << endl;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(dmin,dmax);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data((Z)sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data((Z)sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->extend_1df();
   deps=1.0-epsilon;

   sum=0;
   for(i=0;i<ndoc;i++){
      if(trn[i]==2){
         uu=fabs(ts_sxx[i]);
         xx=pIsg->val_1df((Z)sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         yy=xx/(1.0-xx);
         vv=fabs(0.5*log(yy));
         sum+=vv/(ddx+uu);
      }
   }
   del+=Ppos*sum;

   if(pIsg)delete pIsg;
   pIsg=new Isgrid;
   pIsg->set_xdom(-dmax,-dmin);
   pIsg->set_xgran(grn);
   pIsg->init1();
   for(i=0;i<ndoc;i++){
      if(trn[i]==1)pIsg->add_data(-(Z)sco[i],pdoc[i],pdoc[i]);
      else if(trn[i]==-1)pIsg->add_data(-(Z)sco[i],pdoc[i],0.0);
   }
   pIsg->dim1();
   pIsg->extend_1df();

   sum=0;
   for(i=0;i<ndoc;i++){
      if(trn[i]==2){
         uu=fabs(ts_sxx[i]);
         xx=pIsg->val_1df(-(Z)sco[i]);
         if(xx<epsilon)xx=epsilon;
         if(xx>deps)xx=deps;
         yy=xx/(1.0-xx);
         vv=fabs(0.5*log(yy));
         sum+=vv/(ddx+uu);
      }
   }
   del+=Pneg*sum;
   
   return(del);
}

template <class Y,class Z> 
//ABoostX
ABoostX<Y,Z>::ABoostX(Y trd,Y tsd) {
   tr_doc=trd;
   ts_doc=tsd;
   mark=NULL;
   pdoc=NULL;
   ts_sxx=NULL;
}

template <class Y,class Z> 
ABoostX<Y,Z>::~ABoostX() {
   if(mark!=NULL)delete [] mark;
   if(pdoc!=NULL)delete [] pdoc;
   if(ts_sxx!=NULL)delete [] ts_sxx;
}

template <class Y,class Z> 
void ABoostX<Y,Z>::init(Indx<Y> *gind){
   Y i,j;
   Z xx=1.0/((Z)tr_doc);

   if(mark!=NULL)delete [] mark;
   mark=new Y[tr_doc];
   if(pdoc!=NULL)delete [] pdoc;
   pdoc=new Z[tr_doc];

   for(i=0;i<tr_doc;i++){
      mark[i]=0;
      pdoc[i]=xx;
   }
   for(i=0;i<gind->ix;i++)mark[*(gind->idx+i)]=1;

   if(ts_sxx!=NULL)delete [] ts_sxx;
   ts_sxx=new Z[ts_doc];

   for(i=0;i<ts_doc;i++){
      ts_sxx[i]=0.0;
   }
   t=0;
   zprod=(Z)tr_doc;
}

template <class Y,class Z> 
Z ABoostX<Y,Z>::update(Z *trs,Z *tss){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,sum;
   Z ax,bx,mx;
   tr_sco=trs;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   dmax=dmin=0.0;
   for(i=1;i<tr_doc;i++){
      if(mark[i])xx=trs[i];
      else xx=-trs[i];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   if(dmax==0.0) { cout <<"Max is zero"<<endl; exit(0);}
   ax=300.0/(-dmax);
   if(dmin==0.0) { cout <<"Min is zero"<<endl; exit(0);}
   bx=300.0/(-dmin);

   while(bx-ax>EPS){
      mx=(ax+bx)/2.0; /*midpoint*/
      xx=Z_alpha(mx);
      if(xx>0.0)bx=mx;
      else ax=mx;
   }

   //Update pdoc array
   sum=0;
   for(i=0;i<tr_doc;i++){
      if(mark[i])pdoc[i]*=exp(-ax*trs[i]);
      else pdoc[i]*=exp(ax*trs[i]);
      sum+=pdoc[i];
   }
   zprod=zprod*sum;
   for(i=0;i<tr_doc;i++){
      pdoc[i]/=sum;
   }

   //Update tr_sxx
   for(i=0;i<ts_doc;i++){
      ts_sxx[i]+=ax*tss[i];
   }
   return(ax);
}


template <class Y,class Z> 
Z ABoostX<Y,Z>::Z_alpha(Z alp){
   Y i,j;
   Z xx=0.0;

   for(i=0;i<tr_doc;i++){
      if(mark[i])xx-=pdoc[i]*tr_sco[i]*exp(-alp*tr_sco[i]);
      else xx+=pdoc[i]*tr_sco[i]*exp(alp*tr_sco[i]);
   }
   return(xx);
}

//ABoostX2

template <class Y,class Z> 
ABoostX2<Y,Z>::ABoostX2(Y nd,Indx<Y> *gd,Indx<Y> *bd,Indx<Y> *ts){
   ndoc=nd;
   pdoc=NULL;
   ts_sxx=NULL;
   gdd=gd;
   bdd=bd;
   tst=ts;
   tdoc=gdd->ix+bdd->ix;
}

template <class Y,class Z> 
ABoostX2<Y,Z>::~ABoostX2() {
   if(pdoc!=NULL)delete [] pdoc;
   if(ts_sxx!=NULL)delete [] ts_sxx;
}

template <class Y,class Z> 
void ABoostX2<Y,Z>::init(void){
   Y i,j;

   Z xx=1.0/((Z)tdoc);

   if(pdoc!=NULL)delete [] pdoc;
   pdoc=new Z[ndoc];
   for(i=0;i<ndoc;i++){
      
      pdoc[i]=xx;
   }

   if(ts_sxx!=NULL)delete [] ts_sxx;
   ts_sxx=new Z[ndoc];

   for(i=0;i<ndoc;i++){
      ts_sxx[i]=0.0;
   }
   t=0;
   zprod=(Z)tdoc;
}

template <class Y,class Z> 
Z ABoostX2<Y,Z>::update(Z *scx){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,sum;
   Z ax,bx,mx;
   sco=scx;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   j=gdd->idx[0];
   dmax=dmin=sco[j];
   for(i=1;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx=sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx=-sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }

   if(dmax==0.0) { cout <<"Max is zero"<<endl; exit(0);}
   ax=300.0/(-dmax);
   if(dmin==0.0) { cout <<"Min is zero"<<endl; exit(0);}
   bx=300.0/(-dmin);

   while(bx-ax>EPS){
      mx=(ax+bx)/2.0; /*midpoint*/
      xx=Z_alpha(mx);
      if(xx>0.0)bx=mx;
      else ax=mx;
   }
   //Update pdoc array
   sum=0;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pdoc[j]*=exp(-ax*sco[j]);
      sum+=pdoc[j];
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pdoc[j]*=exp(ax*sco[j]);
      sum+=pdoc[j];
   }
   zprod=zprod*sum;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pdoc[j]/=sum;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pdoc[j]/=sum;
   }

   //Update tr_sxx
   for(i=0;i<tst->ix;i++){
      j=tst->idx[i];
      ts_sxx[j]+=ax*sco[j];
   }
   return(ax);
}


template <class Y,class Z> 
Z ABoostX2<Y,Z>::Z_alpha(Z alp){
   Y i,j;
   Z xx=0.0;

   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx-=pdoc[j]*sco[j]*exp(-alp*sco[j]);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx+=pdoc[j]*sco[j]*exp(alp*sco[j]);
   }
   return(xx);
}

template <class Y,class Z> 
Z ABoostX2<Y,Z>::Z_alphf(Z alp){
   Y i,j;
   Z xx=0.0;

   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx-=pdoc[j]*scf[j]*exp(-alp*scf[j]);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx+=pdoc[j]*scf[j]*exp(alp*scf[j]);
   }
   return(xx);
}

//ABoostX3


template <class Y,class Z> 
ABoostX3<Y,Z>::ABoostX3(Y nd,Indx<Y> *gd,Indx<Y> *bd,Indx<Y> *ts){
   ndoc=nd;
   pdoc=NULL;
   ts_sxx=NULL;
   gdd=gd;
   bdd=bd;
   tst=ts;
   tdoc=gdd->ix+bdd->ix;
}

template <class Y,class Z> 
ABoostX3<Y,Z>::~ABoostX3() {
   if(pdoc!=NULL)delete [] pdoc;
   if(ts_sxx!=NULL)delete [] ts_sxx;
}

template <class Y,class Z> 
void ABoostX3<Y,Z>::init(void){
   Y i,j;

   Z xx=1.0/((Z)tdoc);

   if(pdoc!=NULL)delete [] pdoc;
   pdoc=new Z[ndoc];
   for(i=0;i<ndoc;i++){
      pdoc[i]=xx;
   }

   if(ts_sxx!=NULL)delete [] ts_sxx;
   ts_sxx=new Z[ndoc];

   for(i=0;i<ndoc;i++){
      ts_sxx[i]=0.0;
   }
   t=0;
   zprod=(Z)tdoc;
}

template <class Y,class Z> 
Z ABoostX3<Y,Z>::update(Z *scx){
   int pflag=get_qflag();
   Y i,j,k;
   Z dmax,dmin,xx,sum,ax,bx;
   Z hna,hnb,alphao,betao,lambda=0.0001;
   sco=scx;

   t++;
   if(pflag)cout << t << " Iteration" << endl;

   j=gdd->idx[0];
   dmax=dmin=sco[j];
   for(i=1;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx=sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx=-sco[j];
      dmax=xx>dmax?xx:dmax;
      dmin=xx<dmin?xx:dmin;
   }

   if(dmax==0.0) { cout <<"Max is zero"<<endl; exit(0);}
   ax=300.0/(-dmax);
   if(dmin==0.0) { cout <<"Min is zero"<<endl; exit(0);}
   bx=300.0/(-dmin);

   beta=0;
   hna=bx/5.0;
   alpha=find_alpha(0,hna);
   hnb=30.0;
   beta=find_beta(0,hnb);
   alphao=0;
   betao=0;
   cout << "alpha " << alpha << " beta " << beta << endl;
   while((fabs(alpha-alphao)>lambda)||(fabs(beta-betao)>lambda)){
      alphao=alpha;
      betao=beta;
      alpha=find_alpha(alpha,hna);
      beta=find_beta(beta,hnb);
      cout << "alpha " << alpha << " beta " << beta << endl;
   }
   //Update pdoc array
   sum=0;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pdoc[j]*=exp(-alpha*sco[j]+beta);
      sum+=pdoc[j];
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pdoc[j]*=exp(alpha*sco[j]-beta);
      sum+=pdoc[j];
   }
   zprod=zprod*sum;
   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      pdoc[j]/=sum;
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      pdoc[j]/=sum;
   }

   //Update tr_sxx
   for(i=0;i<tst->ix;i++){
      j=tst->idx[i];
      ts_sxx[j]+=alpha*sco[j]-beta;
   }
   return(alpha);
}

template <class Y,class Z> 
Z ABoostX3<Y,Z>::find_alpha(Z alp,Z &hnv){
   int flag=1,ct=0;
   Z at,bt,mt,xx;
   at=alp-hnv;
   if(Z_alpha(at)>0){
      while(Z_alpha(at-hnv)>0){at-=hnv;ct++;}
      at-=hnv;
      bt=at+hnv;
      flag=0;
   }
   if(flag)bt=alp+hnv;
   if(flag&&(Z_alpha(bt)<0)){
      while(Z_alpha(bt+hnv)<0){bt+=hnv;ct++;}
      bt+=hnv;
      at=bt-hnv;
   }
   if(!ct)hnv=hnv/2.0;
   while(bt-at>EPS){
      mt=(at+bt)/2.0; /*midpoint*/
      xx=Z_alpha(mt);
      if(xx>0.0)bt=mt;
      else at=mt;
   }
   return(at);
}

template <class Y,class Z> 
Z ABoostX3<Y,Z>::find_beta(Z bet,Z &hnv){
   int flag=1,ct=0;
   Z at,bt,mt,xx;
   at=bet-hnv;
   if(Z_beta(at)>0){
      while(Z_beta(at-hnv)>0){at-=hnv;ct++;}
      at-=hnv;
      bt=at+hnv;
      flag=0;
   }
   if(flag)bt=bet+hnv;
   if(flag&&(Z_beta(bt)<0)){
      while(Z_beta(bt+hnv)<0){bt+=hnv;ct++;}
      bt+=hnv;
      at=bt-hnv;
   }
   if(!ct)hnv=hnv/2.0;
   while(bt-at>EPS){
      mt=(at+bt)/2.0; /*midpoint*/
      xx=Z_beta(mt);
      if(xx>0.0)bt=mt;
      else at=mt;
   }
   return(at);
}


template <class Y,class Z> 
Z ABoostX3<Y,Z>::Z_alpha(Z alp){
   Y i,j;
   Z xx=0.0;

   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx-=pdoc[j]*sco[j]*exp(-alp*sco[j]+beta);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx+=pdoc[j]*sco[j]*exp(alp*sco[j]-beta);
   }
   return(xx);
}

template <class Y,class Z> 
Z ABoostX3<Y,Z>::Z_beta(Z bet){
   Y i,j;
   Z xx=0.0;

   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx+=pdoc[j]*exp(-alpha*sco[j]+bet);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx-=pdoc[j]*exp(alpha*sco[j]-bet);
   }
   return(xx);
}

template <class Y,class Z> 
Z ABoostX3<Y,Z>::Z_alphf(Z alp){
   Y i,j;
   Z xx=0.0;

   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx-=pdoc[j]*scf[j]*exp(-alp*scf[j]-beta);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx+=pdoc[j]*scf[j]*exp(alp*scf[j]+beta);
   }
   return(xx);
}

template <class Y,class Z> 
Z ABoostX3<Y,Z>::Z_betf(Z bet){
   Y i,j;
   Z xx=0.0;

   for(i=0;i<gdd->ix;i++){
      j=gdd->idx[i];
      xx+=pdoc[j]*exp(-alpha*scf[j]-bet);
   }
   for(i=0;i<bdd->ix;i++){
      j=bdd->idx[i];
      xx-=pdoc[j]*exp(alpha*scf[j]+bet);
   }
   return(xx);
}

}
#endif
