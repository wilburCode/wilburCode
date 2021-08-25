#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "LexMkv.h"

using namespace std;
namespace iret {

Lex::Lex(void){
   space1=10000;
   str=new char[space1];
   space2=10000;
   addr=new long[space2];
   tx=new char[max_str];
   ht=new long[max_str];
   sn=new long[max_str];
   sm=new long[max_str];
   pflag=1;
}

Lex::~Lex(void){
   delete [] str;
   delete [] addr;
   delete [] tx;
   delete [] ht;
   delete [] sn;
   delete [] sm;
}

char *Lex::show_str(long n){
   if(n<0)return(NULL);
   if(n>=num)return(NULL);
   return(str+addr[n]);
}

void Lex::create_Lex(List &Ls){
   long i=0,j,len;
   long sum1=0,sum2=0;
   char *pch;

   Ls.node_first();
   while(Ls.node_next()){
      len=strlen(Ls.show_str());
      sum1+=len+1;
      sum2++;
      mark(pflag,++i,10000,"strings counted");
   }
   if(sum1>space1){
      space1=sum1;
      delete [] str;
      str=new char[space1];
   }
   if(sum2>space2){
      space2=sum2;
      delete [] addr;
      addr=new long[space2];
   }
   sum1=sum2=i=0;
   Ls.node_first();
   while(Ls.node_next()){
      pch=Ls.show_str();
      addr[sum2]=sum1;
      len=strlen(pch);
      sum1+=len+1;
      sum2++;
      mark(pflag,++i,10000,"strings entered");
   }
   num=sum2;
}

long Lex::find(const char *ssr){
   int j;
   a=b=0;
   if((j=stc_my(ssr,str+addr[0]))<0)return(0);
   else if(j==0)return(1);

   if((j=stc_my(ssr,str+addr[num-1]))>0)return(0);
   else if(j==0)return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_my(ssr,str+addr[i]))==0)return(i+1);
      else if(j<0)y=i;
      else x=i;
   }
   return(0);
}

int Lex::stc_my(const char *ssr,const char *ptr)
   {int i=(a<b) ? a : b;
   const char *p1=ssr+i;
   const char *p2=ptr+i;
   int j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2)return(0);
   else if(*p1<*p2){
      b=i+j;
      return(-1);
   }
   else {
      a=i+j;
      return(1);
   }
}

long Lex::lfind(const char *ssr){
   int i,j,k;
   char *p1;
   const char *p2;

   a=b=0;
   slen=strlen(ssr);
   strcpy(tx,ssr);
   for(i=0;i<slen;i++){
      ht[i]=0;
      sn[i]=0;
      sm[i]=0;
   }

   //Process first string
   p1=tx;
   p2=str+addr[0];
   j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2)return(1);
   else if(*p1<*p2)return(0);
   else {
      if(*p2=='\0'){
         ht[j]=1;
      }
      a=j;
   }
   //Process last string
   p1=tx;
   p2=str+addr[num-1];
   j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2)return(num);
   else if(*p1<*p2){
      b=j;
   }
   else {
      if(*p2=='\0'){
         return(num);
      }
      *p1='\0';
      b=j;
   }

   if(k=ifind(0,num,tx))return(k);
   i=slen-1;
   while(i>0){
      if(ht[i]>0)return(ht[i]);
      else if(ht[i]<0){
         tx[i]='\0';
         if(k=ifind(sn[i],sm[i],tx))return(k);
      }
      i--;
   }
   return(0);
}

long Lex::ifind(long n,long m,const char *ssr){
   int j;
   a=b=0;

   long i,x=n,y=m;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_ly(ssr,i))==0){
         if(a&&!ht[a]){
            ht[a]=-1;
            sn[a]=x;
            sm[a]=i;
         }
         return(i+1);
      }
      else if(j<0)y=i;
      else {
         if((j>1)&&(!ht[j-1])){
            ht[j-1]=-1;
            sn[j-1]=x;
            sm[j-1]=i;
         }
         x=i;
      }
   }
   return(0);
}

int Lex::stc_ly(const char *ssr,long m)
   {int i=(a<b) ? a : b;
   const char *p1=ssr+i;
   const char *p2=str+addr[m]+i;
   int j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2){ht[i+j]=m+1;a=i+j-1;return(0);}
   else if(*p1<*p2){
      b=i+j;
      return(-1);
   }
   else {
      a=i+j;
      if(*p2=='\0'){
         ht[i+j]=m+1;
         return(i+j);
      }
      return(i+j+1);
   }
}

long Lex::tfind(const char *ssr){
   int i,j,k;
   char *p1;
   const char *p2;

   a=b=0;
   slen=strlen(ssr);
   strcpy(tx,ssr);
   for(i=0;i<=slen;i++){
      ht[i]=0;
      sn[i]=0;
      sm[i]=0;
   }

   //Process first string
   p1=tx;
   p2=str+addr[0];
   j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2){ht[j]=1;return(1);}
   else if(*p1<*p2)return(0);
   else {
      if(*p2=='\0'){
         ht[j]=1;
      }
   }
   //Process last string
   p1=tx;
   p2=str+addr[num-1];
   j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2){
      ht[j]=num;
      if(j>1)*(--p1)='\0';
   }
   else if(*p2<*p1){
      if(*p2=='\0'){
         ht[j]=num;
         if(j>1)*(--p1)='\0';
      }
      else *p1='\0';
   }

   ifind(0,num-1,tx);
   i=slen;
   k=0;
   while(i>0){
      if(ht[i]>0)k++;
      else if(ht[i]<0){
         tx[i]='\0';
         if(ifind(sn[i],sm[i],tx))k++;
      }
      i--;
   }
   return(k);
}

long Lex::find_low(const char *ssr){
   long j,k=0;
   a=b=0;
   if((j=stc_low(ssr,str+addr[0]))<0){
      if(j==-2)return(1);
      else return(0);
   }

   if((j=stc_low(ssr,str+addr[num-1]))>0)return(0);
   else if(j==-2)k=j;

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_low(ssr,str+addr[i]))<0){y=i;k=j;}
      else x=i;
   }
   if(k==-2)return(y+1);
   else return(0);
}

int Lex::stc_low(const char *ssr,const char *ptr)
   {int i=(a<b) ? a : b;
   const char *p1=ssr+i;
   const char *p2=ptr+i;
   int j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1<*p2){
      b=i+j;
      if(!*p1)return(-2);
      else return(-1);
   }
   else {
      if(!*p1){b=i+j;return(-2);}
      else {a=i+j;return(1);}
   }
}

long Lex::find_high(const char *ssr){
   long j,k=0;
   a=b=0;
   if((j=stc_high(ssr,str+addr[0]))<0){
      if(j==-1)return(0);
      else if(j==-2)k=j;
   }

   if((j=stc_high(ssr,str+addr[num-1]))>0)return(0);
   else if(j==-2)return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_high(ssr,str+addr[i]))<0){
         if(j==-1)y=i;
         if(j==-2){x=i;k=j;}
      }
      else {x=i;k=j;}
   }
   if(k==-2)return(x+1);
   else return(0);
}

int Lex::stc_high(const char *ssr,const char *ptr)
   {int i=(a<b) ? a : b;
   const char *p1=ssr+i;
   const char *p2=ptr+i;
   int j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1<=*p2){
      if(!*p1){a=i+j;return(-2);}
      else {b=i+j;return(-1);}
   }
   else {
      a=i+j;
      return(1);
   }
}

//Variable order Markov Model

LexMkv::LexMkv(int mn) : Lex(){
   min=mn;
   cnt=new double[space2];
   aug=0;
}

LexMkv::~LexMkv(void){
   delete [] cnt;
}

void LexMkv::set_max(int mx){
   max=mx;
}

void LexMkv::set_aug(double ax){
   aug=ax;
}

void LexMkv::set_fst(char fx,long buf_size){
   buf=new char[buf_size];
   *buf=fx;
   bfc=&(buf[1]);
}

//Suffix Tree Approach

void LexMkv::create_LexMkv(DCount Dc){
   char *pch;
   double xx;
   long i=0,j,len;
   long sum1=0,sum2=0;

   Dc.node_first();
   while(Dc.node_next()){
      pch=Dc.show_str();
      len=strlen(pch)+1;
      if(len>min){
         sum1+=len;
         sum2+=len-min;
      }
      mark(pflag,++i,10000,"strings counted");
   }
   if(sum1>space1){
      space1=sum1;
      delete [] str;
      str=new char[space1];
   }
   if(sum2>space2){
      space2=sum2;
      delete [] addr;
      delete [] cnt;
      addr=new long[space2];
      cnt=new double[space2+1];
   }

   sum1=sum2=i=j=0;
   DCount *pCp=new DCount;
   Dc.node_first();
   while(Dc.node_next()){
      xx=Dc.count();
      pch=Dc.show_str();
      len=strlen(pch);
      if(len>=min){
         strcpy(str+sum1,pch);
         pch=str+sum1;
         pCp->addp_count2(pch,xx+aug);
         for(j=1;j<=len-min;j++){
            pCp->addp_count2(pch+j,xx);
         }
         sum1+=len+1;
      }
      mark(pflag,++i,10000,"strings entered");
   }

   i=(long)str;
   xx=0;
   cnt[0]=xx;
   sum2=0;
   pCp->node_first();
   while(pCp->node_next()){
      j=(long)pCp->show_str()-i;
      xx+=pCp->count();
      addr[sum2]=j;
      cnt[++sum2]=xx;
   }
   tot=pCp->total;
   num=pCp->cnt_key;
   pCp->iclean=1;
   pCp->~DCount();
}

void LexMkv::create_LexMkv(char fx,DCount Dc){
   char *pch;
   double xx;
   long i=0,j,len;
   long sum1=0,sum2=0;

   Dc.node_first();
   while(Dc.node_next()){
      pch=Dc.show_str();
      len=strlen(pch)+2;
      if(len>min+1){
         sum1+=len;
         sum2+=len-min;
      }
      mark(pflag,++i,10000,"strings counted");
   }
   if(sum1>space1){
      space1=sum1;
      delete [] str;
      str=new char[space1];
   }
   if(sum2>space2){
      space2=sum2;
      delete [] addr;
      delete [] cnt;
      addr=new long[space2];
      cnt=new double[space2+1];
   }

   sum1=sum2=i=j=0;
   DCount *pCp=new DCount;
   Dc.node_first();
   while(Dc.node_next()){
      xx=Dc.count();
      pch=Dc.show_str();
      len=strlen(pch);
      if(len>=min){
         str[sum1]=fx;
         strcpy(str+sum1+1,pch);
         pch=str+sum1;
         pCp->addp_count2(pch,xx+aug);
         for(j=1;j<=len-min;j++){
            pCp->addp_count2(pch+j,xx);
         }
         sum1+=len+2;
      }
      mark(pflag,++i,10000,"strings entered");
   }

   i=(long)str;
   xx=0;
   cnt[0]=xx;
   sum2=0;
   pCp->node_first();
   while(pCp->node_next()){
      j=(long)pCp->show_str()-i;
      xx+=pCp->count();
      addr[sum2]=j;
      cnt[++sum2]=xx;
   }
   tot=pCp->total;
   num=pCp->cnt_key;
   pCp->iclean=1;
   pCp->~DCount();
}

double LexMkv::count_Super(const char *stt){
   long i,j;
   double zsum=0;

   if(i=find_low(stt)){
      j=find_high(stt);
      zsum=cnt[j]-cnt[i-1];
   }
   return(zsum);
}

//Probability functions

double LexMkv::prob_ext(const char *stt){
   long i,j,k,m,ok,is;
   char *pch;
   double xx,zz;

   i=strlen(stt);
   pch=new char[i+1];
   strcpy(pch,stt);

   j=1;
   k=1;
   ok=0;
   while(k&&(j<=i)){
      pch[j]='\0';
      k=find_low(pch);
      if(k==0){
         if(j==1){
            delete [] pch;
            return(0);
         }
         else {
            pch[j]=stt[j];
            j--;
            pch[j]='\0';
            m=find_high(pch);
            zz=(cnt[m]-cnt[ok-1])/tot;
            pch[j]=stt[j];
         }
      }
      else {
         if(j<i){
            ok=k;
            pch[j]=stt[j];
            j++;
         }
         else {
            m=find_high(pch);
            zz=(cnt[m]-cnt[k-1])/tot;
            delete [] pch;
            return(zz);
         }
      }
   }
   is=0;
   while(j+is<i){
      is++;
      pch[j+is]='\0';
      while(!(k=find_low(pch+is))){is++;j--;}
      if(j<1){
         delete [] pch;
         return(0);
      }
      pch[j+is]=stt[j+is];
      pch[j-1+is]='\0';
      xx=count_Super(pch+is);
      pch[j-1+is]=stt[j-1+is];

      if(j+is<i){
         ok=k;
         j++;
         while(k&&(j+is<=i)){
            pch[j+is]='\0';
            k=find_low(pch+is);
            if(k==0){
               pch[j+is]=stt[j+is];
               j--;
               pch[j+is]='\0';
               m=find_high(pch+is);
               zz*=(cnt[m]-cnt[ok-1])/xx;
               pch[j+is]=stt[j+is];
            }
            else {
               if(j+is<i){
                  ok=k;
                  pch[j+is]=stt[j+is];
                  j++;
               }
               else {
                  m=find_high(pch+is);
                  zz*=(cnt[m]-cnt[k-1])/xx;
                  delete [] pch;
                  return(zz);
               }
            }
         }
      }
      else {
         m=find_high(pch+is);
         zz*=(cnt[m]-cnt[k-1])/xx;
         delete [] pch;
         return(zz);
      }
   }
   delete [] pch;
   return(zz);
}

double LexMkv::prob_ext_fst(const char *stt){
   strcpy(bfc,stt);
   return(prob_ext(buf));
}

double LexMkv::prob_fxd(const char *stt){
   long i,j,k,m,ok,is,iz;
   char *pch;
   double xx,zz;

   i=strlen(stt);
   pch=new char[i+1];
   strcpy(pch,stt);

   iz=(i<max)?i:max;
   j=1;
   k=1;
   ok=0;
   while(k&&(j<=iz)){
      pch[j]='\0';
      k=find_low(pch);
      if(k==0){
         if(j==1){
            delete [] pch;
            return(0);
         }
         else {
            pch[j]=stt[j];
            j--;
            pch[j]='\0';
            m=find_high(pch);
            zz=(cnt[m]-cnt[ok-1])/tot;
            pch[j]=stt[j];
         }
      }
      else {
         if(j<iz){
            ok=k;
            pch[j]=stt[j];
            j++;
         }
         else {
            m=find_high(pch);
            zz=(cnt[m]-cnt[k-1])/tot;
            if(iz==i){
               delete [] pch;
               return(zz);
            }
            else {
               pch[j]=stt[j];
               k=0;
            }
         }
      }
   }
   is=0;
   while(j+is<i){
      is++;
      pch[j+is]='\0';
      while(!(k=find_low(pch+is))){is++;j--;}
      if(j<1){
         delete [] pch;
         return(0);
      }
      pch[j+is]=stt[j+is];
      pch[j-1+is]='\0';
      xx=count_Super(pch+is);
      pch[j-1+is]=stt[j-1+is];

      iz=is+((max<i-is)?max:(i-is));
      if(j+is<iz){
         ok=k;
         j++;
         while(k&&(j+is<=iz)){
            pch[j+is]='\0';
            k=find_low(pch+is);
            if(k==0){
               pch[j+is]=stt[j+is];
               j--;
               pch[j+is]='\0';
               m=find_high(pch+is);
               zz*=(cnt[m]-cnt[ok-1])/xx;
               pch[j+is]=stt[j+is];
            }
            else {
               if(j+is<iz){
                  ok=k;
                  pch[j+is]=stt[j+is];
                  j++;
               }
               else {
                  m=find_high(pch+is);
                  zz*=(cnt[m]-cnt[k-1])/xx;
                  if(iz==i){
                     delete [] pch;
                     return(zz);
                  }
                  else {
                     pch[j+is]=stt[j+is];
                     k=0;
                  }
               }
            }
         }
      }
      else {
         pch[j+is]='\0';
         m=find_high(pch+is);
         zz*=(cnt[m]-cnt[k-1])/xx;
         if(iz==i){
            delete [] pch;
            return(zz);
         }
         else pch[j+is]=stt[j+is];
      }
   }
   delete [] pch;
   return(zz);
}

double LexMkv::prob_fxd_fst(const char *stt){
   strcpy(bfc,stt);
   return(prob_fxd(buf));
}

void LexMkv::add_Print(DCount &Dc){
   char cx[2],ch;
   int i;
   cx[1]='\0';
   for(i=0;i<128;i++){
      ch=(char)i;
      if(isprint(ch)){
         cx[0]=ch;
         Dc.add_count2(cx,1.0);
      }
   }
}

void LexMkv::add_Char(const char* str,double xx,DCount &Dc){
   char cx[2],ch;
   int i,len;
   cx[1]='\0';
   len=strlen(str);
   for(i=0;i<len;i++){
      ch=str[i];
      if(isprint(ch)){
         cx[0]=ch;
         Dc.add_count2(cx,xx);
      }
   }
}

}
