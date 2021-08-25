#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <runn.h>
#include <Btree.h>
#include <Vnab.h>

using namespace std;
namespace iret {

double xttrc,xb,xa,xg;

   //Local weight functions

float q_const(int lc,long n){
   return(1.0);
}

float d_lc_func(int lc,long n){
   return((float)lc);
}

float d_lc_ratio(int lc,long n){
   return((float)lc/(lc+2.0));
}

float d_lc_log(int lc,long n){
   return((float)log((double)lc+1.0));
}

float s_const(int lc,unsigned int n){
   return(1.0);
}

float sx_const(int lc,long n){
   return(1.0);
}

float s_lc_func(int lc,unsigned int n){
   return((float)lc);
}

float s_lc_ratio(int lc,unsigned int n){
   return((float)lc/(lc+2.0));
}

float s_lc_log(int lc,unsigned int n){
   return((float)log((double)lc+1.0));
}

   //Local TREC functions

float d_len_trec(int lc,long m){
   return(1.0/(1.5+1.5*(float)m/243.128));
}

float d_inquery_trec(int lc,long m){
   return((float)lc/((float)lc+0.5+1.5*(float)m/243.128));
}

float d_robertson_trec(int lc,long m){
   return((float)lc/((float)lc+2.0*(float)m/243.128));
}

float d_wilbur_trec(int lc,long m){
  float u,v;
  double md;
  if(m<dmt)md=dmt;
  else md=(double)m;
  v=(float)exp(md*0.002564+(lc-1.0)*lfac);
  return(1.0/(1.0+v));
}

   //Local Med functions

float d_inquery_med(int lc,long m){
   return((float)lc/((float)lc+0.5+1.5*(float)m/137.924));
}

float d_wilbur_med(int lc,long m){
  float u,v;
  double md;
  if(m<dmt)md=dmt;
  else md=(double)m;
  v=(float)exp(md*0.0044+(lc-1.0)*lfab);
  return(1.0/(1.0+v));
}

float d_string(int lc,long m){
  return(1.0/(5.0+(float)m));
}

float d_prob(int lc,long m){
  return((float)lc/((float)m));
}

//Optimal k1 of 1.9 and b of 1 as in lin & wilbur article
float d_bm25(int lc,long m){
  return((float)lc*2.9/((float)lc+1.9*m/75.3207)); //75.3207 aver doc len in MEDLINE
}

float s_inquery_med(int lc,unsigned int m){
   return((float)lc/((float)lc+0.5+1.5*(float)m/137.924));
}

float s_wilbur_med(int lc,unsigned int m){
  float u,v;
  double md;
  if(m<dmt)md=dmt;
  else md=(double)m;
  v=(float)exp(md*0.0044+(lc-1.0)*lfab);
  return(1.0/(1.0+v));
}

float s_string(int lc,unsigned int m){
  return(1.0/(5.0+(float)m));
}

float s_prob(int lc,unsigned int m){
  return((float)lc/((float)m));
}

   //Global weight functions

float global_idf(long n){
   if(n>dmt)return(l2*log(xttrc/((double)n)));
   else if(n>1)return(l2*log(xttrc/dmt));
   else return(0.0);
}

float global_strict_idf(long n){
   if(n>1)return(l2*log(xttrc/n));
   else return(0.0);
}

float global_sidf(long n){
   if(n>dmt)return(l2*log(xttrc/((double)n)));
   else if(n>0)return(l2*log(xttrc/dmt));
   else return(0.0);
}

float global_const(long n){
   return(1.0);
}

float global_iqf(long n){
   if(n>0){
      float xx=l2*log(xttrc/((double)n));
      return(xx*xx);
   }
   else return(0);
}

float global_bm25(long n){
   float xx=l2*log((xttrc-n+0.5)/(n+0.5));
   return(xx);
}

float sglobal_idf(unsigned int n){
   if(n>dmt)return(l2*log(xttrc/((double)n)));
   else if(n>1)return(l2*log(xttrc/dmt));
   else return(0.0);
}

float sglobal_sidf(unsigned int n){
   if(n>dmt)return(l2*log(xttrc/((double)n)));
   else if(n>0)return(l2*log(xttrc/dmt));
   else return(0.0);
}

float sglobal_const(unsigned int n){
   return(1.0);
}

}
