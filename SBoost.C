#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include "runn.h"
#include "Btree.h"
#include "Word.h"
#include <LUD.h>
#include <Elev.h>
#include <GBoost.h>
#include <Isgrid.h>
#include "SBoost.h"

using namespace std;
namespace iret{
 
 //PAV 
SBoost::SBoost(long dm, long ln, Index *gind, const char *nam): FBase("SBoost", nam){
  dim=dm;
  len=ln;
  gInd=gind; 
  scr=new double*[dim];
  pIs=NULL;
  t=0;
  epsilon=1.0e-12;
  epl=0.0001;
  pflag=1;
}

SBoost::SBoost(long dm,const char *nam): FBase("SBoost", nam){
  dim=dm;
  pflag=1;
  pIs=NULL;
  scr=NULL;
}

SBoost::~SBoost(){
  if(pIs){
     long d;
     for(d=0;d<dim;d++) delete pIs[d];
     delete [] pIs; 
  }
  if(scr!=NULL) delete [] scr;
}

void SBoost::Setup(long i, double *sco){
  scr[i]=sco;
  if(pflag)cout<<"Set up for "<<i<<" completed"<<endl;

}

void SBoost::Learn_Boost(void){
  long l;
  long i,j,k,d, dm, c=0;
  double *dmax,*dmin,xx,sum,deps, **kval, *label, *D, *q, diff=1, q_sum;
  double obj,obo;
  grn=10000;

  kval=new double*[dim]; 
  for(d=0;d<dim;d++) kval[d]=new double[len];
  dmax=new double[dim];
  dmin=new double[dim];
  D=new double[len];
  q=new double[len];
  label=new double[len];
 
  double *pdoc=new double[len];

  for(i=0;i<len;i++){
    pdoc[i]=1;
    for(d=0;d<dim;d++) *(kval[d]+i)=0;
    label[i]=-1;
  }
  for(i=0;i<gInd->ix;i++)label[gInd->idx[i]]=1;

  t++;
  if(pflag)cout << t << " Iteration" << endl;
  for(d=0;d<dim;d++){
    dmax[d]=dmin[d]=*(scr[d]);
    for(i=0;i<len;i++){
       xx=*(scr[d]+i);
       dmax[d]=xx>dmax[d]?xx:dmax[d];
       dmin[d]=xx<dmin[d]?xx:dmin[d];
    }
    if(pflag) cout<<d<<" Min value: "<<dmin[d]<<";Max value: "<<dmax[d]<<endl;
  }
  obj=0;
  while(diff>epl){
    obo=obj;
    if(pflag)cout<<"round="<<c<<endl;
   
    if(pIs){
       for(d=0;d<dim;d++) delete pIs[d];
       delete [] pIs; 
    }
    pIs=(Isgrid**)new long[dim];
    for(d=0;d<dim;d++){
       pIs[d]=new Isgrid;
       pIs[d]->pflag=pflag;
    }
    for(d=0;d<dim;d++){ 
      deps=1.0-epsilon;
      sum=0;
      for(i=0;i<len;i++){
        D[i]=0;
        for(dm=0;dm<dim;dm++) D[i]+=*(kval[dm]+i);
        D[i]-=*(kval[d]+i);
        pdoc[i]=exp(-(label[i]*D[i]));
        sum+=pdoc[i];
      }
      pIs[d]->set_xdom(dmin[d],dmax[d]);
      pIs[d]->set_xgran(grn);
      pIs[d]->init1();
      for(i=0;i<len;i++){
        pdoc[i]=pdoc[i]/sum;
        if(label[i]>0) pIs[d]->add_data(*(scr[d]+i),pdoc[i],pdoc[i]);
        else pIs[d]->add_data(*(scr[d]+i),pdoc[i],0.0);
      } 
      pIs[d]->dim1();
      pIs[d]->extend_1df();
      if(pflag) cout << "Average p " << pIs[d]->avg() << " Information " << pIs[d]->info() << endl;
      for(i=0;i<len;i++){
        xx=pIs[d]->val_1df(*(scr[d]+i));
        if(xx<epsilon)xx=epsilon;
        if(xx>deps)xx=deps;
        *(kval[d]+i)=0.5*log(xx/(1-xx));
      }
    }
    obj=0;
    for(i=0;i<len;i++){ 
      q_sum=0;
      for(d=0;d<dim;d++) q_sum+=*(kval[d]+i);
      q[i]=exp(-q_sum*label[i]);
      obj+=q[i];
    }
    if(c>0) diff=obo-obj;
    else diff=obj;
    c++;
  } 
  for(d=0;d<dim;d++) delete [] kval[d];
  delete [] kval;
  delete [] dmax;
  delete [] dmin;
  delete [] D;
  delete [] q;
  delete [] label;
  delete [] pdoc;
}

void SBoost::Save_Boost(void){
   long d,l;
   char cnam[100];
   for(d=0;d<dim;d++){
      pIs[d]->set_name(add_num(name, d, cnam));
      pIs[d]->write_1df();
   }
}

void SBoost::Load_Boost(void){
  long i,d,l;
  char cnam[100];
  if(pIs){
     for(d=0;d<dim;d++) delete pIs[d];
     delete [] pIs; 
  }
  pIs=(Isgrid**)new long[dim];
  for(i=0;i<dim;i++){
    pIs[i]=new Isgrid(add_num(name, i, cnam));
    pIs[i]->read_1df();
    pIs[i]->extend_1df();
  }
}

//Finds the score by adding up all ts_sxx
double SBoost::Score(double *scc){
  long d;
  double xx,score=0,deps;
  deps=1.0-epsilon;
  for(d=0;d<dim;d++){
    xx=pIs[d]->val_1df(scc[d]);
    if(xx<epsilon) xx=epsilon;
    if(xx>deps) xx=deps; 
    score+=0.5*log(xx/(1.0-xx));
  } 
  return score;
}

  //Optimal Alpha Method

SABoost::SABoost(long dm, long ln, Index *gind, const char *nam): FBase("SABoost", nam){
  dim=dm;
  len=ln;
  epsilon=1.0E-12;
  epl=0.01;
  sum=new double[len];
  gInd=gind; 
  scr=new double*[dim];
  alpha=new double[dim];
  alpha_h=NULL;
  label=NULL;
}

SABoost::SABoost(long dm,const char *nam): FBase("SABoost", nam){
  dim=dm;
  alpha=new double[dim];
  alpha_h=NULL;
  label=NULL;
  scr=NULL;
  sum=NULL;
}

SABoost::~SABoost(){
  long d;
  if(scr!=NULL)delete [] scr;
  if(sum!=NULL)delete [] sum;
  if(alpha_h!=NULL){
     for(d=0;d<dim;d++) if(alpha_h[d]) delete [] alpha_h[d];
     delete [] alpha_h;
  }
  if(alpha!=NULL) delete [] alpha;
  if(label!=NULL) delete [] label;
}

void SABoost::Setup(long i, double *sco){
  scr[i]=sco;
  if(pflag)cout<<"Set up for "<<i<<" completed"<<endl;

}

void SABoost::Learn_Boost(){
  long i,j,k,d, dm, c=0, round, tc;
  double *dmax,*dmin,xx,deps;
  double sum_n, sum_o=1000, diff=1, pow;
  double *ax,*bx,mx;

  alpha_h=new double*[dim]; 
  for(d=0;d<dim;d++) alpha_h[d]=new double[len];
  dmax=new double[dim];
  dmin=new double[dim];
  ax=new double[dim];
  bx=new double[dim];
  label=new double[len];

  for(i=0;i<len;i++){
    for(d=0;d<dim;d++) *(alpha_h[d]+i)=0;
    label[i]=-1;
  }
  for(i=0;i<gInd->ix;i++)label[gInd->idx[i]]=1;

  for(d=0;d<dim;d++){
    dmax[d]=dmin[d]=*(scr[d]);
    for(i=0;i<len;i++){
      if(label[i]>0) xx=*(scr[d]+i);
      else xx=-*(scr[d]+i);
      dmax[d]=xx>dmax[d]?xx:dmax[d];
      dmin[d]=xx<dmin[d]?xx:dmin[d];
    }  
    if(dmax[d]==0.0) { cout <<"Max is zero"<<endl; exit(0);}
    if(dmin[d]==0.0) { cout <<"Min is zero"<<endl; exit(0);}
  }

  round=0;
  while(diff>epl){
    if(pflag)cout<<"round="<<round<<endl; 
    for(d=0;d<dim;d++){
      ax[d]=300.0/(-dmax[d]);
      bx[d]=300.0/(-dmin[d]);
      while(bx[d]-ax[d]>epsilon){
        mx=(ax[d]+bx[d])/2.0; /*midpoint*/
        xx=Z_alpha(mx, d);
        if(xx>0.0)bx[d]=mx;
        else ax[d]=mx;
      }
      alpha[d]=mx;
      for(i=0;i<len;i++) *(alpha_h[d]+i)=*(scr[d]+i)*alpha[d];
     
      sum_n=0;
      for(i=0;i<len;i++){ 
        pow=0;
        for(c=0;c<dim;c++) pow+=*(alpha_h[c]+i);
        if(label[i]>0) sum_n+=exp(-pow);
        else sum_n+=exp(pow);
      }
      diff=sum_o-sum_n;
      sum_o=sum_n;
      if(pflag)cout << "dim " << d << " alpha " << mx << " sum " << sum_o << endl;
    }
    round++;
  }
  delete [] dmax;
  delete [] dmin;
  delete [] ax;
  delete [] bx;
} 
 
void SABoost::Save_Boost(void){
 //Write final thetas or coefficients out in ".cf"-regular and "bcf"-binary files
  long c;
  ofstream *pfout, *ptfout;
  pfout=get_Ostr("cf");
  ptfout=get_Ostr("bcf");
  *pfout<<"Dimention"<<endl<<'\t'<<dim<<endl;
  ptfout->write((char*)&dim, sizeof(long));
  *pfout<<"Optimal Coef. Are:"<<endl;
  for(c=0;c<dim;c++){
    *pfout<<"\t"<<alpha[c]<<endl;
    ptfout->write((char*)(alpha+c), sizeof(double));
  }
  dst_Ostr(pfout);
  dst_Ostr(ptfout);
}

void SABoost::Load_Boost(){
  long dmn;
  ifstream *pfin;
  pfin=get_Istr("bcf");
  pfin->read((char*)&dmn, sizeof(long));
  if(dim==dmn){
    char cnam[10000];
    double val;
    long i;
    //Reads in optimal alpha's
    pfin->read((char*)alpha, (dim+1)*sizeof(double));
    dst_Istr(pfin);
  }
  else{
    cout<<"Error, Dimentions Do Not Match!"<<endl;
    exit(0);
  }
}
 
double SABoost::Z_alpha(double alp, long dm){
   long i,j,d;
   double xx=0.0;
   for(i=0;i<len;i++){
     sum[i]=0;
     for(d=0;d<dm;d++) sum[i]+=*(alpha_h[d]+i);
     for(d=dm+1;d<dim;d++) sum[i]+=*(alpha_h[d]+i);
     if(label[i]>0) xx-=*(scr[dm]+i)*exp(-sum[i])*exp(-*(scr[dm]+i)*alp);
     else xx+=*(scr[dm]+i)*exp(sum[i])*exp(*(scr[dm]+i)*alp);
   }
   return(xx);
}


//Finds the score by adding up weighted scores
double SABoost::Score(double *scc){
  double score=0;
  long d;
  for(d=0;d<dim;d++) score+=scc[d]*alpha[d];
  return score;
}

}


