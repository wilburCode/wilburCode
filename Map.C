#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <cmath>
#include <cstring>
#include <cassert>
#include "runn.h"
#include "Btree.h"
#include "Map.h"

#define alp 0.1
#define mag 500

using namespace std;
namespace iret {

Map::Map(void) : FBase("mapset","null"){
}

Map::Map(const char *nam) : FBase("mapset",nam){
}

Map::~Map(){
}

void Map::create_Map(List &Lst){
   char *wrd,ch[max_str];
   long len,h[max_str],p[max_str];
   long ct,i,j,k,wctu=Lst.cnt_key;
   long n,zt,rt,pc;

   for(ct=0;ct<max_str;ct++)h[ct]=0;
   zt=0;
   Lst.node_first();
   Lst.node_next();
   strcpy(ch,Lst.show_str());
   n=0;
   ct=1;

   while(Lst.node_next()){
      wrd=Lst.show_str();
      len=strlen(wrd);
      rt=0;
      while(ch[rt]==wrd[rt])rt++;
      if(h[rt]==ct-1){
         if(rt>zt)n+=rt-zt+1;
         else n++;
      }
      for(i=rt;i<len+1;i++){
         ch[i]=wrd[i];
         h[i]=ct;
      }
      zt=rt;
      mark(++ct,10000,"term space count");
   }     

   if(pflag)cout << "No. terms " << wctu << " pt space " << n << endl;

   ad=new long[wctu];
   st=new long[wctu];
   pt=new long[n];

   for(ct=0;ct<n;ct++)pt[ct]=0;
   ofstream *pfout=get_Ostr("c",ios::out);
   for(ct=0;ct<max_str;ct++)h[ct]=0;
   zt=0;
   Lst.node_first();
   Lst.node_next();
   strcpy(ch,Lst.show_str());
   len=strlen(ch);
   for(i=0;i<len+1;i++)pfout->put(ch[i]);
   ad[0]=st[0]=0;
   pc=len+1;
   ct=1;

   while(Lst.node_next()){
      ad[ct]=pc;
      wrd=Lst.show_str();     
      len=strlen(wrd);
      rt=0;
      while(ch[rt]==wrd[rt])rt++;
      if(h[rt]<ct-1){
         pt[p[rt]]=ct;
         st[ct]=st[ct-1];
      }
      else {
         k=st[ct-1];
         for(i=zt;i<rt;i++){
            p[i]=k;
            k++;
         }
         st[ct]=k+1;
         pt[k]=ct;
      }
      for(i=rt;i<len+1;i++){
         ch[i]=wrd[i];
         h[i]=ct;
         pfout->put(wrd[i]);
      }
      pc+=len+1-rt;
      zt=rt;
      mark(++ct,10000,"term write count");
   }

   dst_Ostr(pfout);
   bin_Writ("a",sizeof(long)*wctu,(char*)ad);
   bin_Writ("s",sizeof(long)*wctu,(char*)st);
   bin_Writ("p",sizeof(long)*n,(char*)pt);
   delete [] ad;
   delete [] st;
   delete [] pt;
}

void Map::create_Map(long nx,List &Lst){
   char *wrd,ch[max_str];
   long len,h[max_str],p[max_str];
   long ct,i,j,k,wctu=Lst.cnt_key;
   long n,zt,rt,pc;

   for(ct=0;ct<max_str;ct++)h[ct]=0;
   zt=0;
   Lst.node_first();
   Lst.node_next();
   strcpy(ch,Lst.show_str());
   n=0;
   ct=1;

   while(Lst.node_next()){
      wrd=Lst.show_str();
      len=strlen(wrd);
      rt=0;
      while(ch[rt]==wrd[rt])rt++;
      if(h[rt]==ct-1){
         if(rt>zt)n+=rt-zt+1;
         else n++;
      }
      for(i=rt;i<len+1;i++){
         ch[i]=wrd[i];
         h[i]=ct;
      }
      zt=rt;
      mark(++ct,10000,"term space count");
   }

   if(pflag)cout << "No. terms " << wctu << " pt space " << n << endl;

   ad=new long[wctu];
   st=new long[wctu];
   pt=new long[n];

   for(ct=0;ct<n;ct++)pt[ct]=0;
   ofstream *pfout=get_Ostr(nx,"c",ios::out);
   for(ct=0;ct<max_str;ct++)h[ct]=0;
   zt=0;
   Lst.node_first();
   Lst.node_next();
   strcpy(ch,Lst.show_str());
   len=strlen(ch);
   for(i=0;i<len+1;i++)pfout->put(ch[i]);
   ad[0]=st[0]=0;
   pc=len+1;
   ct=1;

   while(Lst.node_next()){
      ad[ct]=pc;
      wrd=Lst.show_str();
      len=strlen(wrd);
      rt=0;
      while(ch[rt]==wrd[rt])rt++;
      if(h[rt]<ct-1){
         pt[p[rt]]=ct;
         st[ct]=st[ct-1];
      }
      else {
         k=st[ct-1];
         for(i=zt;i<rt;i++){
            p[i]=k;
            k++;
         }
         st[ct]=k+1;
         pt[k]=ct;
      }
      for(i=rt;i<len+1;i++){
         ch[i]=wrd[i];
         h[i]=ct;
         pfout->put(wrd[i]);
      }
      pc+=len+1-rt;
      zt=rt;
      mark(++ct,10000,"term write count");
   }

   dst_Ostr(pfout);
   bin_Writ(nx,"a",sizeof(long)*wctu,(char*)ad);
   bin_Writ(nx,"s",sizeof(long)*wctu,(char*)st);
   bin_Writ(nx,"p",sizeof(long)*n,(char*)pt);
   delete [] ad;
   delete [] st;
   delete [] pt;
}

void Map::gopen_map(void){
   char cnam[max_str],*cptr;
   int fld;
   long ct,asize,i;
   
   cf=get_Mmap("c");
   ad=(long*)get_Mmap("a");
   st=(long*)get_Mmap("s");
   pt=(long*)get_Mmap("p");
}

void Map::gopen_map(long nx){
   char cnam[max_str],*cptr;
   int fld;
   long ct,asize,i;
  
   cf=get_Mmap(nx,"c");
   ad=(long*)get_Mmap(nx,"a");
   st=(long*)get_Mmap(nx,"s");
   pt=(long*)get_Mmap(nx,"p");
}

void Map::gopen_Map_copy(Map *pM){
   cf=pM->cf;
   ad=pM->ad;
   st=pM->st;
   pt=pM->pt;
}

long Map::find(const char *wrd){
   long i=1,j,k;
   long n=0;
   const char *xx=cf;
   long len;

   len=strlen(wrd);
   repeat:
      j=0;
      while((i<len)&&(wrd[i]==*(xx+j))){i++;j++;}
      if(wrd[i]>*(xx+j)){
         k=st[n]+j;
         if((k<st[n+1])&&(n=pt[k])){
            xx=cf+ad[n];
            goto repeat;
         }
         else return(0);
      }
      else if(*(xx+j)=='\0')return(n+1);
      else return(0);
}

long Map::bfind(const char *wrd){
   long i=0,j,k;
   long n=0;
   const char *xx=cf;
   long len;

   len=strlen(wrd);
   repeat:
      j=0;
      while((i<len)&&(wrd[i]==*(xx+j))){i++;j++;}
      if(wrd[i]>*(xx+j)){
         k=st[n]+j;
         if((k<st[n+1])&&(n=pt[k])){
            xx=cf+ad[n];
            goto repeat;
         }
         else return(0);
      }
      else if(*(xx+j)=='\0')return(n+1);
      else return(0);
}

long Map::sfind(const char *wrd,State &Sk){
   long i=1,j,k;
   long n=0;
   const char *xx=cf;
   long len;
   Sk.dp=Sk.nm=1;
   Sk.tn[0]=0;
   Sk.mt[0]=0;
   Sk.ch[0]='|';

   len=strlen(wrd);
   repeat:
      j=0;
      while((i<len)&&(wrd[i]==(Sk.ch[i]=*(xx+j)))){i++;j++;}
      if(j){Sk.tn[Sk.dp]=n;Sk.mt[Sk.dp]=j;Sk.dp++;Sk.nm+=j;}
      if(wrd[i]>*(xx+j)){
         k=st[n]+j;
         if((k<st[n+1])&&(n=pt[k])){
            xx=cf+ad[n];
            goto repeat;
         }
         else return(0);
      }
      else if(*(xx+j)=='\0')return(n+1);
      else return(0);
}

long Map::sbfind(const char *wrd,State &Sk){
   long i=0,j,k;
   long n=0;
   const char *xx=cf;
   long len;
   Sk.dp=Sk.nm=0;

   len=strlen(wrd);
   repeat:
      j=0;
      while((i<len)&&(wrd[i]==(Sk.ch[i]=*(xx+j)))){i++;j++;}
      if(j){Sk.tn[Sk.dp]=n;Sk.mt[Sk.dp]=j;Sk.dp++;Sk.nm+=j;}
      if(wrd[i]>*(xx+j)){
         k=st[n]+j;
         if((k<st[n+1])&&(n=pt[k])){
            xx=cf+ad[n];
            goto repeat;
         }
         else return(0);
      }
      else if(*(xx+j)=='\0')return(n+1);
      else return(0);
}

long Map::extend(const char *wrd,long m,State &Sk){
   long i=m,j,k;
   long n=Sk.tn[Sk.dp-1];
   const char *xx=cf+ad[n];
   long len,ip;

   len=strlen(wrd);
   if(i>len)return(0);
   j=Sk.mt[Sk.dp-1];
   ip=Sk.nm;
   while((i<len)&&(wrd[i]==(Sk.ch[ip]=*(xx+j)))){i++;ip++;j++;}
   Sk.mt[Sk.dp-1]=j;
   Sk.nm=ip;
   if(wrd[i]>*(xx+j)){
      k=st[n]+j;
      if((k<st[n+1])&&(n=pt[k])){
         xx=cf+ad[n];
      }
      else return(0);
   }
   else if(*(xx+j)=='\0')return(n+1);
   else return(0);

   repeat:
      j=0;
      while((i<len)&&(wrd[i]==(Sk.ch[ip]=*(xx+j)))){i++;ip++;j++;}
      if(j){Sk.tn[Sk.dp]=n;Sk.mt[Sk.dp]=j;Sk.dp++;Sk.nm+=j;}
      if(wrd[i]>*(xx+j)){
         k=st[n]+j;
         if((k<st[n+1])&&(n=pt[k])){
            xx=cf+ad[n];
            goto repeat;
         }
         else return(0);
      }
      else if(*(xx+j)=='\0')return(n+1);
      else return(0);
}

long Map::extend_s(const char *wrd,long m,State &Sk){
   long i=m,j,k;
   long n=Sk.tn[Sk.dp-1];
   const char *xx=cf+ad[n];
   long len,ip;

   len=strlen(wrd);
   if(i>len)return(-1);
   j=Sk.mt[Sk.dp-1];
   ip=Sk.nm;
   while((i<len)&&(wrd[i]==(Sk.ch[ip]=*(xx+j)))){i++;ip++;j++;}
   Sk.mt[Sk.dp-1]=j;
   Sk.nm=ip;
   if(wrd[i]>*(xx+j)){
      k=st[n]+j;
      if((k<st[n+1])&&(n=pt[k])){
         xx=cf+ad[n];
      }
      else {
         while(Sk.ch[ip]=wrd[i]){ip++;i++;}
         return(0);
      }
   }
   else if(*(xx+j)=='\0'){
      while(Sk.ch[ip]=wrd[i]){ip++;i++;}
      return(n+1);
   }
   else {
      while(Sk.ch[ip]=wrd[i]){ip++;i++;}
      return(0);
   }

   repeat:
      j=0;
      while((i<len)&&(wrd[i]==(Sk.ch[ip]=*(xx+j)))){i++;ip++;j++;}
      if(j){Sk.tn[Sk.dp]=n;Sk.mt[Sk.dp]=j;Sk.dp++;Sk.nm+=j;}
      if(wrd[i]>*(xx+j)){
         k=st[n]+j;
         if((k<st[n+1])&&(n=pt[k])){
            xx=cf+ad[n];
            goto repeat;
         }
         else {
            while(Sk.ch[ip]=wrd[i]){ip++;i++;}
            return(0);
         }
      }
      else if(*(xx+j)=='\0'){
         while(Sk.ch[ip]=wrd[i]){ip++;i++;}
         return(n+1);
      }
      else {
         while(Sk.ch[ip]=wrd[i]){ip++;i++;}
         return(0);
      }
}

long Map::step_Back(State &St){
   switch(St.dp){
      case 1: if(St.mt[0]){
                 St.mt[0]--;
                 St.nm--;
                 return(St.nm+1);
              }
              else return(0);
              break;
      default: 
         if(St.mt[St.dp-1]>1){
            (St.mt[St.dp-1])--;
            St.nm--;
         }
         else {
            St.dp--;
            St.nm--;
         }
         return(St.nm+1);
   }
}

long Map::step_Forward(State &St){
   long id=St.dp-1;
   long n=St.tn[id],j=St.mt[id];
   long k;
   char c=*(cf+ad[n]+j);

   if(c){
      St.mt[id]++;
      St.ch[St.nm]=c;
      St.nm++;
      return(1);
   }
   else {
      k=st[n]+j;
      if((k<st[n+1])&&(n=pt[k])){
         St.tn[St.dp]=n;
         St.ch[St.nm]=*(cf+ad[n]);
         St.nm++;
         St.mt[St.dp]=1;
         St.dp++;
         return(1);
      }
      else return(0);
   }
}

long Map::step_Down(State &St){
   long id=St.dp-1;
   long n=St.tn[id],j=St.mt[id]-1;
   long k=st[n]+j;

   if((k<st[n+1])&&(n=pt[k])){
      if(j){
         St.mt[id]--;
         St.tn[St.dp]=n;
         St.ch[St.nm-1]=*(cf+ad[n]);
         St.mt[St.dp]=1;
         St.dp++;
      }
      else {
         St.tn[id]=n;
         St.ch[St.nm-1]=*(cf+ad[n]);
      }   
      return(1);
   }
   else return(0);
}

void Map::global_one_Edit(long lev,const char* str,State &Sk,Count &Ct){
   lxn=strlen(str);
   ptr=new char[lxn+1];
   strcpy(ptr,str);

   local_one_Edit(str,Sk,Ct);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_one_Edit(str,Sk,Ct);
   }
   delete [] ptr;
}

void Map::local_one_Edit(const char* str,State &Sk,Count &Ct){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(n=extend(ptr,m,Sk)){
      Sk.ch[Sk.nm]='\0';
      Ct.max_count(Sk.ch,n);
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(n=extend(ptr,m,Sk)){
         Sk.ch[Sk.nm]='\0';
         Ct.max_count(Sk.ch,n);
      }
      Sk.Reset(Lx);
      if(n=extend(ptr,m-1,Sk)){
         Sk.ch[Sk.nm]='\0';
         Ct.max_count(Sk.ch,n);
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(n=extend(ptr,m,Sk)){
            Sk.ch[Sk.nm]='\0';
            Ct.max_count(Sk.ch,n);
         }
         Sk.Reset(Lx);
         if(n=extend(ptr,m-1,Sk)){
            Sk.ch[Sk.nm]='\0';
            Ct.max_count(Sk.ch,n);
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(n=extend(ptr,m,Sk)){
               Sk.ch[Sk.nm]='\0'; 
               Ct.max_count(Sk.ch,n); 
            } Sk.Reset(Lx);
            if(n=extend(ptr,m-1,Sk)){
               Sk.ch[Sk.nm]='\0';
               Ct.max_count(Sk.ch,n);
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m<lxn){
      Sk.Reset(Lv);
      ptr[m]=str[m-1];
      ptr[m-1]=str[m];
      if(n=extend(ptr,m-1,Sk)){
         Sk.ch[Sk.nm]='\0';
         Ct.max_count(Sk.ch,n);
      }
      ptr[m]=str[m];
      ptr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

void Map::global_two_Edit(long lev,const char* str,State &Sk,Count &Ct){
   lxn=strlen(str);
   qtr=new char[lxn+1];
   utr=new char[lxn+3];
   strcpy(qtr,str);

   local_two_Edit(str,Sk,Ct);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_two_Edit(str,Sk,Ct);
   }
   delete [] qtr;
   delete [] utr;
}

void Map::local_two_Edit(const char* str,State &Sk,Count &Ct){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(!(n=extend_s(qtr,m,Sk))){
      strcpy(utr,Sk.ch);
      global_one_Edit(m,utr,Sk,Ct);
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(!(n=extend_s(qtr,m,Sk))){
         strcpy(utr,Sk.ch);
         global_one_Edit(m+1,utr,Sk,Ct);
      }
      Sk.Reset(Lx);
      if(!(n=extend_s(qtr,m-1,Sk))){
         strcpy(utr,Sk.ch);
         global_one_Edit(m+1,utr,Sk,Ct);
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(!(n=extend_s(qtr,m,Sk))){
            strcpy(utr,Sk.ch);
            global_one_Edit(m+1,utr,Sk,Ct);
         }
         Sk.Reset(Lx);
         if(!(n=extend_s(qtr,m-1,Sk))){
            strcpy(utr,Sk.ch);
            global_one_Edit(m+1,utr,Sk,Ct);
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(!(n=extend_s(qtr,m,Sk))){
               strcpy(utr,Sk.ch);
               global_one_Edit(m+1,utr,Sk,Ct);
            }
            Sk.Reset(Lx);
            if(!(n=extend_s(qtr,m-1,Sk))){
               strcpy(utr,Sk.ch);
               global_one_Edit(m+1,utr,Sk,Ct);
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m<lxn){
      Sk.Reset(Lv);
      qtr[m]=str[m-1];
      qtr[m-1]=str[m];
      if(!(n=extend_s(qtr,m-1,Sk))){
         strcpy(utr,Sk.ch);
         global_one_Edit(m+1,utr,Sk,Ct);
      }
      qtr[m]=str[m];
      qtr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

//CMap

CMap::CMap(void) : Map(){
}

CMap::CMap(const char *nam,EdtPrb *pEt) : Map(nam){
      change_type("cmapset");
      pEp=pEt;
}

CMap::~CMap(void){
}

void CMap::create_CMap(Count &Ct){
   create_Map(Ct);
   long n,i=0;
   freq=new long[Ct.cnt_key];
   Ct.node_first();
   while(Ct.node_next()){
      freq[i]=Ct.count();
      mark(++i,10000,"count terms");
   }
   bin_Writ("f",Ct.cnt_key*sizeof(long),(char*)freq);
   delete [] freq;
}

void CMap::create_CMap(long nx,Count &Ct){
   create_Map(Ct);
   long n,i=0;
   freq=new long[Ct.cnt_key];
   Ct.node_first();
   while(Ct.node_next()){
      freq[i]=Ct.count();
      mark(++i,10000,"count terms");
   }
   bin_Writ(nx,"f",Ct.cnt_key*sizeof(long),(char*)freq);
   delete [] freq;
}

void CMap::gopen_cmap(void){
   long i;
   gopen_map();
   freq=(long*)get_Mmap("f");

   pEp->gopen_map();
   pSq=new State(max_str);
   n1=pEp->n1;
   n2=pEp->n2;
   n3=pEp->n3;
   cn=pEp->con;
   div_del=pEp->div_del;
   div_ins=pEp->div_ins;
   div_rep=pEp->div_rep;
   div_trn=pEp->div_trn;

   eps=0.075;
   for(i=0;i<81;i++)fac[i]=pow(10.0,eps*(i-80));
}

void CMap::gopen_cmap(long nx){
   long i;
   gopen_map(nx);
   pSq=new State(max_str);
   freq=(long*)get_Mmap(nx,"f");

   pEp->gopen_map();
   n1=pEp->n1;
   n2=pEp->n2;
   n3=pEp->n3;
   cn=pEp->con;
   div_del=pEp->div_del;
   div_ins=pEp->div_ins;
   div_rep=pEp->div_rep;
   div_trn=pEp->div_trn;

   eps=0.075;
   for(i=0;i<81;i++)fac[i]=pow(10.0,eps*(i-80));
}

void CMap::gopen_CMap_copy(CMap *pCM){
   long i;
   gopen_Map_copy((Map *)pCM);
   freq=pCM->freq;

   pSq=new State(max_str);
   n1=pCM->n1;
   n2=pCM->n2;
   n3=pCM->n3;
   cn=pCM->cn;
   div_del=pCM->div_del;
   div_ins=pCM->div_ins;
   div_rep=pCM->div_rep;
   div_trn=pCM->div_trn;

   for(i=0;i<81;i++)fac[i]=pCM->fac[i];
}

long CMap::cfind(const char *wrd,State &Sk,long lm){
   long i=Sk.nm,j,k=Sk.dp,ip,kp;
   long n=Sk.tn[k-1];
   const char *xx=cf+ad[n]+Sk.mt[k-1];
   Sk.mt[k]=0;
   long len=strlen(wrd);
   long cj=0,min;
   ixt=0;

   while(k){
      i-=Sk.mt[k];
      n=Sk.tn[k-1];
      xx=cf+ad[n]+Sk.mt[k-1];
      k--;
      if(*xx)continue;
      min=freq[n];

      ip=i;
      n=0;
      xx=cf;
     repeat:
      j=0;
      while((ip<len)&&(wrd[ip]==*(xx+j))){ip++;j++;}
      if(wrd[ip]>*(xx+j)){
         kp=st[n]+j;
         if((kp<st[n+1])&&(n=pt[kp])){
            xx=cf+ad[n];
            goto repeat;
         }
      }
      else if(*(xx+j)=='\0'){
         min=(min<freq[n])?min:freq[n];
         if(min>cj&&min>=lm){
            ixt=i;
            cj=min;
         }
      }
   }
   if(cj){
      for(j=0;j<ixt;j++)Sk.ch[j]=wrd[j];
      Sk.ch[ixt]='|';
      for(j=ixt;j<len;j++)Sk.ch[j+1]=wrd[j];
      Sk.ch[len+1]='\0';
      return(1);
   }
   else return(0);
}

long CMap::cfind2(const char *wrd,State &Sk,long lm){
   long i=Sk.nm,j,k=Sk.dp,ip,kp;
   long n=Sk.tn[k-1];
   const char *xx=cf+ad[n]+Sk.mt[k-1];
   Sk.mt[k]=0;
   long len=strlen(wrd);
   long cj=0,min;
   ixt=0;

   while(k){
      i-=Sk.mt[k];
      n=Sk.tn[k-1];
      xx=cf+ad[n]+Sk.mt[k-1];
      k--;
      if(*xx)continue;
      min=freq[n];

      ip=i;
      n=0;
      xx=cf;
     repeat:
      j=0;
      while((ip<len)&&(wrd[ip]==*(xx+j))){ip++;j++;}
      if(wrd[ip]>*(xx+j)){
         kp=st[n]+j;
         if((kp<st[n+1])&&(n=pt[kp])){
            xx=cf+ad[n];
            goto repeat;
         }
      }
      else if(*(xx+j)=='\0'){
         min=(min<freq[n])?min:freq[n];
         if(min>cj&&min>=lm){
            ixt=i;
            cj=min;
         }
      }
   }
   if(cj){
      for(j=0;j<ixt;j++)Sk.ch[j]=wrd[j];
      Sk.ch[ixt]='\0';
      return(1);
   }
   else return(0);
}

long CMap::check_sngl(const char *wrd,long l1,long l2){
   long n,mn,min=0;
   char *ppr=new char[l1+l2+2];
   strncpy(ppr,wrd+1,l1);
   ppr[l1]='\0';
   if((n=bfind(ppr))||(l1<=4)){
      if(n)min=freq[n-1];
      else min=1000;
      strncpy(ppr,wrd+l1+2,l2);
      ppr[l2]='\0';
      if(n=bfind(ppr)){
         mn=freq[n-1];
         min=(min<mn)?min:mn;
         return(min);
      }
      else if(l2<=4){
         mn=1000;
         min=(min<mn)?min:mn;
         return(min);
      }
      else return(0);
   }
   else return(0);
}

//XEdit

void CMap::global_one_XEdit(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   ptr=new char[lxn+1];
   strcpy(ptr,str);

   local_one_XEdit(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_one_XEdit(str,Sk);
   }
   delete [] ptr;
}

void CMap::global_ph1_XEdit(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   ptr=new char[lxn+1];
   strcpy(ptr,str);

   if(e1h+1<lev)return;
   while(Sk.nm>e1h)step_Back(Sk);
   local_one_XEdit(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_one_XEdit(str,Sk);
   }
   delete [] ptr;
}

void CMap::local_one_XEdit(const char* str,State &Sk){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(n=extend(ptr,m,Sk)){
      def1=def2*div_ins[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m]]];
      if(def1*freq[n-1]>fmax){
         Sk.ch[Sk.nm]='\0';
         strcpy(cmax,Sk.ch);
         fmax=def1*freq[n-1];
      }
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(n=extend(ptr,m,Sk)){
         def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)ptr[m]]];
         if(def1*freq[n-1]>fmax){
            Sk.ch[Sk.nm]='\0';
            strcpy(cmax,Sk.ch);
            fmax=def1*freq[n-1];
         }
      }
      Sk.Reset(Lx);
      if(n=extend(ptr,m-1,Sk)){
         def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
         if(def1*freq[n-1]>fmax){
            Sk.ch[Sk.nm]='\0';
            strcpy(cmax,Sk.ch);
            fmax=def1*freq[n-1];
         }
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(n=extend(ptr,m,Sk)){
            def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)ptr[m]]];
            if(def1*freq[n-1]>fmax){
               Sk.ch[Sk.nm]='\0';
               strcpy(cmax,Sk.ch);
               fmax=def1*freq[n-1];
            }
         }
         Sk.Reset(Lx);
         if(n=extend(ptr,m-1,Sk)){
            def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
            if(def1*freq[n-1]>fmax){
               Sk.ch[Sk.nm]='\0';
               strcpy(cmax,Sk.ch);
               fmax=def1*freq[n-1];
            }
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(n=extend(ptr,m,Sk)){
               def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)ptr[m]]];
               if(def1*freq[n-1]>fmax){
                  Sk.ch[Sk.nm]='\0';
                  strcpy(cmax,Sk.ch);
                  fmax=def1*freq[n-1];
               }
            } 
            Sk.Reset(Lx);
            if(n=extend(ptr,m-1,Sk)){
               def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
               if(def1*freq[n-1]>fmax){
                  Sk.ch[Sk.nm]='\0';
                  strcpy(cmax,Sk.ch);
                  fmax=def1*freq[n-1];
               }
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m<lxn){
      Sk.Reset(Lv);
      ptr[m]=str[m-1];
      ptr[m-1]=str[m];
      if(n=extend(ptr,m-1,Sk)){
         def1=def2*div_trn[cn[(int)ptr[m-2]]+n1*cn[(int)ptr[m-1]]+n2*cn[(int)ptr[m]]+\
              n3*cn[(int)ptr[m+1]]];
         if(def1*freq[n-1]>fmax){
            Sk.ch[Sk.nm]='\0';
            strcpy(cmax,Sk.ch);
            fmax=def1*freq[n-1];
         }
      }
      ptr[m]=str[m];
      ptr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

void CMap::global_two_XEdit(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   qtr=new char[lxn+1];
   utr=new char[lxn+3];
   strcpy(qtr,str);

   local_two_XEdit(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_two_XEdit(str,Sk);
   }
   delete [] qtr;
   delete [] utr;
}

void CMap::local_two_XEdit(const char* str,State &Sk){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if((n=extend_s(qtr,m,Sk))>-1){
      strcpy(utr,Sk.ch);
      def2=div_ins[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m]]];
      global_one_XEdit(m,utr,Sk);
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if((n=extend_s(qtr,m,Sk))>-1){
         strcpy(utr,Sk.ch);
         def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)qtr[m]]];
         global_one_XEdit(m+1,utr,Sk);
      }
      Sk.Reset(Lx);
      if((n=extend_s(qtr,m-1,Sk))>-1){
         strcpy(utr,Sk.ch);
         def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
         global_one_XEdit(m+1,utr,Sk);
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if((n=extend_s(qtr,m,Sk))>-1){
            strcpy(utr,Sk.ch);
            def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)qtr[m]]];
            global_one_XEdit(m+1,utr,Sk);
         }
         Sk.Reset(Lx);
         if((n=extend_s(qtr,m-1,Sk))>-1){
            strcpy(utr,Sk.ch);
            def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
            global_one_XEdit(m+1,utr,Sk);
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if((n=extend_s(qtr,m,Sk))>-1){
               strcpy(utr,Sk.ch);
               def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)qtr[m]]];
               global_one_XEdit(m+1,utr,Sk);
            }
            Sk.Reset(Lx);
            if((n=extend_s(qtr,m-1,Sk))>-1){
               strcpy(utr,Sk.ch);
               def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
               global_one_XEdit(m+1,utr,Sk);
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m+1<lxn){
      Sk.Reset(Lv);
      qtr[m]=str[m-1];
      qtr[m-1]=str[m];
      if((n=extend_s(qtr,m-1,Sk))>-1){
         strcpy(utr,Sk.ch);
         def2=div_trn[cn[(int)qtr[m-2]]+n1*cn[(int)qtr[m-1]]+n2*cn[(int)qtr[m]]+\
              n3*cn[(int)qtr[m+1]]];
         global_one_XEdit(m+1,utr,Sk);
      }
      qtr[m]=str[m];
      qtr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

//ph2_XEdit

void CMap::global_ph2_XEdit(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   qtr=new char[lxn+1];
   utr=new char[lxn+3];
   strcpy(qtr,str);

   while(Sk.nm>e2h)step_Back(Sk);
   local_ph2_XEdit(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_ph2_XEdit(str,Sk);
   }
   delete [] qtr;
   delete [] utr;
}

void CMap::local_ph2_XEdit(const char* str,State &Sk){
   long i,j,k,n,m,u1,u2;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   u2=(e1w>m+1)?e1w:m+1;
   if(!(n=extend_s(qtr,m,Sk))&&(Sk.nm>=e1w)){
      strcpy(utr,Sk.ch);
      def2=div_ins[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m]]];
      e1h--;
      global_ph1_XEdit(u2-1,utr,Sk);
      e1h++;
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(!(n=extend_s(qtr,m,Sk))&&(Sk.nm>=e1w)){
         strcpy(utr,Sk.ch);
         def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)qtr[m]]];
         global_ph1_XEdit(u2,utr,Sk);
      }
      Sk.Reset(Lx);
      if(!(n=extend_s(qtr,m-1,Sk))&&(Sk.nm>=e1w)){
         strcpy(utr,Sk.ch);
         def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
         e1h++;
         global_ph1_XEdit(u2,utr,Sk);
         e1h--;
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(!(n=extend_s(qtr,m,Sk))&&(Sk.nm>=e1w)){
            strcpy(utr,Sk.ch);
            def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)qtr[m]]];
            global_ph1_XEdit(u2,utr,Sk);
         }
         Sk.Reset(Lx);
         if(!(n=extend_s(qtr,m-1,Sk))&&(Sk.nm>=e1w)){
            strcpy(utr,Sk.ch);
            def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
            e1h++;
            global_ph1_XEdit(u2,utr,Sk);
            e1h--;
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(!(n=extend_s(qtr,m,Sk))&&(Sk.nm>=e1w)){
               strcpy(utr,Sk.ch);
               def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)qtr[m]]];
               global_ph1_XEdit(u2,utr,Sk);
            }
            Sk.Reset(Lx);
            if(!(n=extend_s(qtr,m-1,Sk))&&(Sk.nm>e1w)){
               strcpy(utr,Sk.ch);
               def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
               e1h++;
               global_ph1_XEdit(u2,utr,Sk);
               e1h--;
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m+1<lxn){
      Sk.Reset(Lv);
      qtr[m]=str[m-1];
      qtr[m-1]=str[m];
      if(!(n=extend_s(qtr,m-1,Sk))&&(Sk.nm>e1w)){
         strcpy(utr,Sk.ch);
         def2=div_trn[cn[(int)qtr[m-2]]+n1*cn[(int)qtr[m-1]]+n2*cn[(int)qtr[m]]+\
              n3*cn[(int)qtr[m+1]]];
         global_ph1_XEdit(u2,utr,Sk);
      }
      qtr[m]=str[m];
      qtr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

//SEdit

void CMap::global_one_SEdit(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   ptr=new char[lxn+1];
   strcpy(ptr,str);

   local_one_SEdit(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_one_SEdit(str,Sk);
   }
   delete [] ptr;
}

void CMap::global_ph1_SEdit(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   ptr=new char[lxn+1];
   strcpy(ptr,str);

   while(Sk.nm>e1h)step_Back(Sk);
   local_one_SEdit(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_one_SEdit(str,Sk);
   }
   delete [] ptr;
}

void CMap::local_one_SEdit(const char* str,State &Sk){
   long i,j,k,n,m;
   long ndp;
   Level Lv,Lx;
   double zz;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(n=extend(ptr,m,Sk)){
      def1=def2*div_ins[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m]]];
      if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
      else sumx+=zz=def1*fac[ndp]*ndp;
      if(zz>fmax){
         Sk.ch[Sk.nm]='\0';
         strcpy(cmax,Sk.ch);
         fmax=zz;
      }
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(n=extend(ptr,m,Sk)){
         def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)ptr[m]]];
         if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
         else sumx+=zz=def1*fac[ndp]*ndp;
         if(zz>fmax){
            Sk.ch[Sk.nm]='\0';
            strcpy(cmax,Sk.ch);
            fmax=zz;
         }
      }
      Sk.Reset(Lx);
      if(n=extend(ptr,m-1,Sk)){
         def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
         if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
         else sumx+=zz=def1*fac[ndp]*ndp;
         if(zz>fmax){
            Sk.ch[Sk.nm]='\0';
            strcpy(cmax,Sk.ch);
            fmax=zz;
         }
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(n=extend(ptr,m,Sk)){
            def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)ptr[m]]];
            if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
            else sumx+=zz=def1*fac[ndp]*ndp;
            if(zz>fmax){
               Sk.ch[Sk.nm]='\0';
               strcpy(cmax,Sk.ch);
               fmax=zz;
            }
         }
         Sk.Reset(Lx);
         if(n=extend(ptr,m-1,Sk)){
            def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
            if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
            else sumx+=zz=def1*fac[ndp]*ndp;
            if(zz>fmax){
               Sk.ch[Sk.nm]='\0';
               strcpy(cmax,Sk.ch);
               fmax=zz;
            }
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(n=extend(ptr,m,Sk)){
               def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)ptr[m]]];
               if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
               else sumx+=zz=def1*fac[ndp]*ndp;
               if(zz>fmax){
                  Sk.ch[Sk.nm]='\0';
                  strcpy(cmax,Sk.ch);
                  fmax=zz;
               }
            } 
            Sk.Reset(Lx);
            if(n=extend(ptr,m-1,Sk)){
               def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
               if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
               else sumx+=zz=def1*fac[ndp]*ndp;
               if(zz>fmax){
                  Sk.ch[Sk.nm]='\0';
                  strcpy(cmax,Sk.ch);
                  fmax=zz;
               }
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m<lxn){
      Sk.Reset(Lv);
      ptr[m]=str[m-1];
      ptr[m-1]=str[m];
      if(n=extend(ptr,m-1,Sk)){
         def1=def2*div_trn[cn[(int)ptr[m-2]]+n1*cn[(int)ptr[m-1]]+n2*cn[(int)ptr[m]]+\
              n3*cn[(int)ptr[m+1]]];
         if((ndp=freq[n-1])>80)sumx+=zz=def1*freq[n-1];
         else sumx+=zz=def1*fac[ndp]*ndp;
         if(zz>fmax){
            Sk.ch[Sk.nm]='\0';
            strcpy(cmax,Sk.ch);
            fmax=zz;
         }
      }
      ptr[m]=str[m];
      ptr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

//Exten

void CMap::global_one_Exten(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   ptr=new char[lxn+1];
   strcpy(ptr,str);

   local_one_Exten(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_one_Exten(str,Sk);
   }
   delete [] ptr;
}

void CMap::local_one_Exten(const char* str,State &Sk){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(!(n=extend_s(ptr,m,Sk))){
      def1=def2*div_ins[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m]]];
      if((Sk.nm+lnr2>=lmax)&&((Sk.nm+lnr2>lmax)||(def1>fmax))){
         strcpy(cmax,Sk.ch);
         lmax=Sk.nm+lnr2;
         fmax=def1;
      }
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(!(n=extend_s(ptr,m,Sk))){
         def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)ptr[m]]];
         if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
            strcpy(cmax,Sk.ch);
            lmax=Sk.nm+lnr2-1;
            fmax=def1;
         }
      }
      Sk.Reset(Lx);
      if(!(n=extend_s(ptr,m-1,Sk))){
         def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
         if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
            strcpy(cmax,Sk.ch);
            lmax=Sk.nm+lnr2-1;
            fmax=def1;
         }
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(!(n=extend_s(ptr,m,Sk))){
            def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)ptr[m]]];
            if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
               strcpy(cmax,Sk.ch);
               lmax=Sk.nm+lnr2-1;
               fmax=def1;
            }
         }
         Sk.Reset(Lx);
         if(!(n=extend_s(ptr,m-1,Sk))){
            def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
            if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
               strcpy(cmax,Sk.ch);
               lmax=Sk.nm+lnr2-1;
               fmax=def1;
            }
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(!(n=extend_s(ptr,m,Sk))){
               def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)ptr[m]]];
               if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
                  strcpy(cmax,Sk.ch);
                  lmax=Sk.nm+lnr2-1;
                  fmax=def1;
               }
            }
            Sk.Reset(Lx);
            if(!(n=extend_s(ptr,m-1,Sk))){
               def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
               if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
                  strcpy(cmax,Sk.ch);
                  lmax=Sk.nm+lnr2-1;
                  fmax=def1;
               }
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m<lxn){
      Sk.Reset(Lv);
      ptr[m]=str[m-1];
      ptr[m-1]=str[m];
      if(!(n=extend_s(ptr,m-1,Sk))){
         def1=def2*div_trn[cn[(int)ptr[m-2]]+n1*cn[(int)ptr[m-1]]+n2*cn[(int)ptr[m]]+\
              n3*cn[(int)ptr[m+1]]];
         if((Sk.nm+lnr2>=lmax)&&((Sk.nm+lnr2>lmax)||(def1>fmax))){
            strcpy(cmax,Sk.ch);
            lmax=Sk.nm+lnr2;
            fmax=def1;
         }
      }
      ptr[m]=str[m];
      ptr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

void CMap::global_two_Exten(long lev,const char* str,State &Sk){
   lxn=strlen(str);
   qtr=new char[lxn+1];
   utr=new char[lxn+3];
   strcpy(qtr,str);
   fmax=0;
   lmax=0;

   local_two_Exten(str,Sk);
   while(Sk.nm>=lev){
      step_Back(Sk);
      local_two_Exten(str,Sk);
   }
   delete [] qtr;
   delete [] utr;
}

void CMap::local_two_Exten(const char* str,State &Sk){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(!(n=extend_s(qtr,m,Sk))){
      strcpy(utr,Sk.ch);
      def2=div_ins[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m]]];
      lnr2=0;
      global_one_Exten(m,utr,Sk);
   }
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(!(n=extend_s(qtr,m,Sk))){
         strcpy(utr,Sk.ch);
         def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)qtr[m]]];
         lnr2=-1;
         global_one_Exten(m+1,utr,Sk);
      }
      Sk.Reset(Lx);
      if(!(n=extend_s(qtr,m-1,Sk))){
         strcpy(utr,Sk.ch);
         def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
         lnr2=-1;
         global_one_Exten(m+1,utr,Sk);
      }
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(!(n=extend_s(qtr,m,Sk))){
            strcpy(utr,Sk.ch);
            def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)qtr[m]]];
            lnr2=-1;
            global_one_Exten(m+1,utr,Sk);
         }
         Sk.Reset(Lx);
         if(!(n=extend_s(qtr,m-1,Sk))){
            strcpy(utr,Sk.ch);
            def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
            lnr2=-1;
            global_one_Exten(m+1,utr,Sk);
         }
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(!(n=extend_s(qtr,m,Sk))){
               strcpy(utr,Sk.ch);
               def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)qtr[m]]];
               lnr2=-1;
               global_one_Exten(m+1,utr,Sk);
            }
            Sk.Reset(Lx);
            if(!(n=extend_s(qtr,m-1,Sk))){
               strcpy(utr,Sk.ch);
               def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
               lnr2=-1;
               global_one_Exten(m+1,utr,Sk);
            }
            Sk.Reset(Lx);
         }
      }
   }
   if(m+1<lxn){
      Sk.Reset(Lv);
      qtr[m]=str[m-1];
      qtr[m-1]=str[m];
      if(!(n=extend_s(qtr,m-1,Sk))){
         strcpy(utr,Sk.ch);
         def2=div_trn[cn[(int)qtr[m-2]]+n1*cn[(int)qtr[m-1]]+n2*cn[(int)qtr[m]]+\
              n3*cn[(int)qtr[m+1]]];
         lnr2=0;
         global_one_Exten(m+1,utr,Sk);
      }
      qtr[m]=str[m];
      qtr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
}

//Qxten

long CMap::global_one_Qxten(long lev,const char* str,State &Sk){
   if(lev<7)return(0);
   long n;
   lxn=strlen(str);
   ptr=new char[lxn+2];
   strcpy(ptr,str);
   ptr[lxn+1]='\0';

   n=local_one_Qxten(str,Sk);
   while((!n)&&(Sk.nm>=lev)){
      step_Back(Sk);
      n=local_one_Qxten(str,Sk);
   }
   delete [] ptr;
   return(n);
}

long CMap::local_one_Qxten(const char* str,State &Sk){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(!(n=extend(ptr,m,Sk))){
      def1=def2*div_ins[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m]]];
      if((Sk.nm+lnr2>=lmax)&&((Sk.nm+lnr2>lmax)||(def1>fmax))){
         pSq->Copy(Sk);
         em=en;
         ek=-1;
         lmax=Sk.nm+lnr2;
         fmax=def1;
      }
   }
   else return(n);
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(!(n=extend(ptr,m,Sk))){
         def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)ptr[m]]];
         if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
            pSq->Copy(Sk);
            em=en;
            ek=0;
            lmax=Sk.nm+lnr2-1;
            fmax=def1;
         }
      }
      else return(n);
      Sk.Reset(Lx);
      if(!(n=extend(ptr,m-1,Sk))){
         def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
         if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
            pSq->Copy(Sk);
            em=en;
            ek=1;
            lmax=Sk.nm+lnr2-1;
            fmax=def1;
         }
      }
      else return(n);
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(!(n=extend(ptr,m,Sk))){
            def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)ptr[m]]];
            if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
               pSq->Copy(Sk);
               em=en;
               ek=0;
               lmax=Sk.nm+lnr2-1;
               fmax=def1;
            }
         }
         else return(n);
         Sk.Reset(Lx);
         if(!(n=extend(ptr,m-1,Sk))){
            def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
            if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
               pSq->Copy(Sk);
               em=en;
               ek=1;
               lmax=Sk.nm+lnr2-1;
               fmax=def1;
            }
         }
         else return(n);
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(!(n=extend(ptr,m,Sk))){
               def1=def2*div_rep[cn[(int)ptr[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)ptr[m]]];
               if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
                  pSq->Copy(Sk);
                  em=en;
                  ek=0;
                  lmax=Sk.nm+lnr2-1;
                  fmax=def1;
               }
            }
            else return(n);
            Sk.Reset(Lx);
            if(!(n=extend(ptr,m-1,Sk))){
               def1=def2*div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)ptr[m-2]]+n2*cn[(int)ptr[m-1]]];
               if((Sk.nm+lnr2-1>=lmax)&&((Sk.nm+lnr2-1>lmax)||(def1>fmax))){
                  pSq->Copy(Sk);
                  em=en;
                  ek=1;
                  lmax=Sk.nm+lnr2-1;
                  fmax=def1;
               }
            }
            else return(n);
            Sk.Reset(Lx);
         }
      }
   }
   if(m<lxn){
      Sk.Reset(Lv);
      ptr[m]=str[m-1];
      ptr[m-1]=str[m];
      if(!(n=extend(ptr,m-1,Sk))){
         def1=def2*div_trn[cn[(int)ptr[m-2]]+n1*cn[(int)ptr[m-1]]+n2*cn[(int)ptr[m]]+\
              n3*cn[(int)ptr[m+1]]];
         if((Sk.nm+lnr2>=lmax)&&((Sk.nm+lnr2>lmax)||(def1>fmax))){
            pSq->Copy(Sk);
            em=en;
            ek=0;
            lmax=Sk.nm+lnr2;
            fmax=def1;
         }
      }
      ptr[m]=str[m];
      ptr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
   return(n);
}

long CMap::global_two_Qxten(long lev,const char* str,State &Sk){
   long n;
   lxn=strlen(str);
   qtr=new char[lxn+1];
   utr=new char[lxn+3];
   strcpy(qtr,str);
   fmax=0;
   lmax=0;

   n=local_two_Qxten(str,Sk);
   while((!n)&&(Sk.nm>=lev)){
      step_Back(Sk);
      n=local_two_Qxten(str,Sk);
   }
   delete [] qtr;
   delete [] utr;
   return(n);
}

long CMap::local_two_Qxten(const char* str,State &Sk){
   long i,j,k,n,m;
   long xtn,xmt,xnm,xdp;
   Level Lv,Lx;

   Sk.Save(Lv);
   m=Sk.nm+1;
   if(!(n=extend_s(qtr,m,Sk))){
      strcpy(utr,Sk.ch);
      def2=div_ins[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m]]];
      lnr2=0;
      en=-1;
      if(n=global_one_Qxten(m,utr,Sk))return(n);
   }
   else if(n>0)return(n);
   Sk.Reset(Lv);
   if(step_Forward(Sk)){
      Sk.Save(Lx);
      if(!(n=extend_s(qtr,m,Sk))){
         strcpy(utr,Sk.ch);
         def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
              n3*cn[(int)qtr[m]]];
         lnr2=-1;
         en=0;
         if(n=global_one_Qxten(m+1,utr,Sk))return(n);
      }
      else if(n>0)return(n);
      Sk.Reset(Lx);
      if(!(n=extend_s(qtr,m-1,Sk))){
         strcpy(utr,Sk.ch);
         def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
         lnr2=-1;
         en=1;
         if(n=global_one_Qxten(m+1,utr,Sk))return(n);
      }
      else if(n>0)return(n);
      Sk.Reset(Lx);
      if(step_Down(Sk)){
         Sk.Save(Lx);
         if(!(n=extend_s(qtr,m,Sk))){
            strcpy(utr,Sk.ch);
            def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                 n3*cn[(int)qtr[m]]];
            lnr2=-1;
            en=0;
            if(n=global_one_Qxten(m+1,utr,Sk))return(n);
         }
         else if(n>0)return(n);
         Sk.Reset(Lx);
         if(!(n=extend_s(qtr,m-1,Sk))){
            strcpy(utr,Sk.ch);
            def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
            lnr2=-1;
            en=1;
            if(n=global_one_Qxten(m+1,utr,Sk))return(n);
         }
         else if(n>0)return(n);
         Sk.Reset(Lx);
         while(step_Down(Sk)){
            Lx.tn=Sk.tn[Sk.dp-1];
            if(!(n=extend_s(qtr,m,Sk))){
               strcpy(utr,Sk.ch);
               def2=div_rep[cn[(int)qtr[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)Sk.ch[m-1]]+\
                    n3*cn[(int)qtr[m]]];
               lnr2=-1;
               en=0;
               if(n=global_one_Qxten(m+1,utr,Sk))return(n);
            }
            else if(n>0)return(n);
            Sk.Reset(Lx);
            if(!(n=extend_s(qtr,m-1,Sk))){
               strcpy(utr,Sk.ch);
               def2=div_del[cn[(int)Sk.ch[m-1]]+n1*cn[(int)qtr[m-2]]+n2*cn[(int)qtr[m-1]]];
               lnr2=-1;
               en=1;
               if(n=global_one_Qxten(m+1,utr,Sk))return(n);
            }
            else if(n>0)return(n);
            Sk.Reset(Lx);
         }
      }
   }
   if(m+1<lxn){
      Sk.Reset(Lv);
      qtr[m]=str[m-1];
      qtr[m-1]=str[m];
      if(!(n=extend_s(qtr,m-1,Sk))){
         strcpy(utr,Sk.ch);
         def2=div_trn[cn[(int)qtr[m-2]]+n1*cn[(int)qtr[m-1]]+n2*cn[(int)qtr[m]]+\
              n3*cn[(int)qtr[m+1]]];
         lnr2=0;
         en=0;
         if(n=global_one_Qxten(m+1,utr,Sk))return(n);
      }
      qtr[m]=str[m];
      qtr[m-1]=str[m-1];
   }
   Sk.Reset(Lv);
   if(n>0)return(n);
   else return(0);
}

//Search routines

long CMap::Search_one(const char *str){
   long len=strlen(str),n,fq;
   State Sk(len+3),Sp(len+3);
   double pmax=0;

   if(n=sfind(str,Sk)){
      if((fq=freq[n-1])>1000)return(2);
      else if(fq>80)pmax=(double)fq;
      else pmax=fq*fac[fq]; 
   }

   sumx=pmax;
   def2=1.0;
   fmax=0;
   global_one_SEdit(2,str,Sk);
   if(!n){
      if(fmax)return(1);
      else return(0);
   }
   else if((fmax>0.7*sumx)||(pmax<0.05*sumx))return(1);
   else return(3);
}

long CMap::Search_ph1(const char *str){
   long len=strlen(str),n,fq;
   char c;
   State Sk(len+3),Sp(len+3);
   double pmax=0;

   ln1=0;
   while(((c=str[ln1])!=' ')&&(c!='\0'))ln1++;
   if(c)ln2=strlen(str+ln1+1);
   else ln2=0;
   ln1--;

   if(ln2){
      if(ln1<=2)e2w=ln1+3;
      else e2w=2;
      if(ln2<=2)e1h=ln1+1;
      else e1h=ln1+ln2+3;
   }
   else {
      e2w=2;
      e1h=ln1+1;
   }

   if(n=sfind(str,Sk)){
      if((fq=freq[n-1])>1000)return(2);
      else if(fq>80)pmax=(double)fq;
      else pmax=fq*fac[fq];
   }

   sumx=pmax;
   def2=1.0;
   fmax=0;
   global_ph1_SEdit(e2w,str,Sk);
   if(!n){
      if(fmax)return(1);
      else return(0);
   }
   else if((fmax>0.7*sumx)||(pmax<0.05*sumx))return(1);
   else return(3);
}

long CMap::Search_two(const char *str){
   long len=strlen(str);
   char dnam[max_str];
   State Sk(len+3);

   if(sfind(str,Sk))step_Back(Sk);
   global_two_XEdit(2,str,Sk);
   if(fmax){
      strcpy(dnam,cmax);
      if(Search_one(dnam)==1){
         strcpy(dnam,cmax);
         if(Search_one(dnam)!=1)strcpy(cmax,dnam);
      }
      else strcpy(cmax,dnam);
      return(1);
   }
   else return(0);
}

long CMap::Search_ph2(const char *str){
   long len=strlen(str);
   char dnam[max_str],c;
   State Sk(len+3);
   long i,j,k,u;

   ln1=0;
   while(((c=str[ln1])!=' ')&&(c!='\0'))ln1++;
   if(c)ln2=strlen(str+ln1+1);
   else ln2=0;
   ln1--;

   if(ln1<=2)e2w=e1w=ln1+3;
   else if(ln1<7){e2w=2;e1w=ln1+3;}
   else e2w=e1w=2;

   if(ln2<=2)e2h=e1h=ln1;
   else if(ln2<7){e2h=ln1;e1h=ln1+ln2+1;}
   else e2h=e1h=ln1+ln2+1;

   if(sfind(str,Sk))step_Back(Sk);
   global_ph2_XEdit(e2w,str,Sk);

   if(fmax){
      if((ln1>=7)||(ln2>=7)){
         strcpy(dnam,cmax);
         if(Search_ph1(dnam)==1){
            strcpy(dnam,cmax);
            if(Search_ph1(dnam)!=1)strcpy(cmax,dnam);
         }
         else strcpy(cmax,dnam);
         return(1);
      }
      else return(1);
   }
   else return(0);
}

long CMap::Search_cmb(const char *str){
   char dnam[max_str];
   long n,len,pf,pq;
   double pmax;
 
   fmax=0;
   len=strlen(str);
   if(len<6)return(0);
   n=Search_one(str);
   switch(n){
      case 3: strcpy(cmax,str);
              break;
      case 2: strcpy(cmax,str);
              break;
      case 1: if(len>=6){
                 pmax=fmax;
                 strcpy(dnam,cmax);
                 if(!(Search_one(dnam)==1)){
                    fmax=pmax;
                    strcpy(cmax,dnam);
                 }
              }
              break;
      case 0: if(len>=9){
                 if(Search_spl(str,500))n=4;
                 if(n==0){
                    if(Search_two(str))n=5;
                 }
              }
   }
   if((len>=9)&&((n==3)||(n==1))){
      strcpy(dnam,cmax);
      if(n=find(dnam))pf=freq[n-1];
      else pf=0;
      if((pf<80)&&Search_two(str)){
         pq=freq[find(cmax)-1];
         if((pq>80)&&(pq>10*pf)&&(RBeg(cmax,str)))n=5;
         else strcpy(cmax,dnam);
      }
      else strcpy(cmax,dnam);
   }
   return(n);
}

long CMap::Search_cmg(const char *str){
   char dnam[max_str];
   long n,len;
   double pmax;

   fmax=0;
   len=strlen(str);
   if(len<7)return(0);
   if(irp<3)n=Search_ph1(str);
   else n=Search_one(str);
   switch(n){
      case 3: strcpy(cmax,str);
              break;
      case 2: strcpy(cmax,str);
              break;
      case 1: if(len>=7){
                 pmax=fmax;
                 strcpy(dnam,cmax);
                 if(irp<3){
                    if(!(Search_ph1(dnam)==1)){
                       fmax=pmax;
                       strcpy(cmax,dnam);
                    }
                 }
                 else {
                    if(!(Search_one(dnam)==1)){
                       fmax=pmax;
                       strcpy(cmax,dnam);
                    }
                 }
              }
              break;
      case 0: if(len>=8){
                 if(Search_spg(str,500))n=5;
                 if((len>=15)&&n==0){
                    if(irp==2){
                       if(Search_ph2(str))n=4;
                    }
                    else if(irp==3){
                       if(Search_two(str))n=4;
                    }
                 }
              }
   }
   return(n);
}

long CMap::Search_cm1(const char *str){
   char dnam[max_str];
   long n,len;
   double pmax;

   fmax=0;
   len=strlen(str);
   if(len<7)return(0);
   if(irp<3)n=Search_ph1(str);
   else n=Search_one(str);
   switch(n){
      case 3: strcpy(cmax,str);
              break;
      case 2: strcpy(cmax,str);
              break;
      case 1: if(len>=7){
                 pmax=fmax;
                 strcpy(dnam,cmax);
                 if(irp<3){
                    if((!(Search_ph1(dnam)==1))||(!Rat_check(str,cmax))){
                       fmax=pmax;
                       strcpy(cmax,dnam);
                    }
                 }
                 else {
                    if((!(Search_one(dnam)==1))||(!Rat_check(str,cmax))){
                       fmax=pmax;
                       strcpy(cmax,dnam);
                    }
                 }
              }
              break;
      case 0: if(len>=8){
                 if(irp==2){
                    if(Search_ph2(str))n=4;
                 }
                 else if(irp==3){
                    if(Search_two(str))n=4;
                 }
              }
   }
   return(n);
}

long CMap::Search_cm2(const char *str){
   char dnam[max_str];
   long n,len;
   double pmax;

   fmax=0;
   len=strlen(str);
   n=Search_one(str);
   switch(n){
      case 3: strcpy(cmax,str);
              break;
      case 2: strcpy(cmax,str);
              break;
      case 1: if(len>=8){
                 pmax=fmax;
                 strcpy(dnam,cmax);
                 if(!(Search_one(dnam)==1)){
                    fmax=pmax;
                    strcpy(cmax,dnam);
                 }
              }
              break;
      case 0: if(len>=8){
                 if(Search_two(str))n=4;
                 else n=0;
              }
              else n=0;
   }
   return(n);
}

long CMap::Search_dep(const char *str){
   long lmx0=-1,flag=1;
   long len=strlen(str),i;
   long itr=len/4-1,lsm=0;
   char *ztr=new char[2*len+1];
   strcpy(ztr,str);
   State Sk(3*len);

   lmax=0;
   for(i=0;i<itr;i++){
      lmx0=lmax;
      sfind(ztr,Sk);
      global_two_Exten(2,ztr,Sk);
      lsm+=lmax;
      if(lmx0+1<lmax){
         strcpy(ztr,cmax);
         if(Search_cm2(ztr)){
            delete [] ztr;
            lmax=lsm;
            if(Rat_check(cmax,str))return(i+1);
            else return(0);
         }
      }
      else {
         delete [] ztr;
         return(0);
      }
   }
   delete [] ztr;
   return(0);
}

long CMap::Search_qst(const char *str,State &Sk){
   long lmx0,flag,lev,n;
   long len=strlen(str),i,j,k;
   long itr=len/4-1,iz,sxx,sxo=0;
   char *ztr=new char[3*len+1],*pcx;
   strcpy(ztr,str);
   char *vtr=new char[3*len+1];
   State Sh(3*len+1);

   Sh.Copy(Sk);
   lmx0=lmax=Sk.nm;
   for(i=0;i<itr;i++){
      lev=(2<Sk.nm-8)?Sk.nm-8:2;
      if(n=global_two_Qxten(lev,ztr,Sk))break;
      if(lmx0+4<lmax){
         Sh.Copy(*pSq);
         Sk.Copy(*pSq);
         Sk.ch[Sk.nm]='\0';
         strcpy(vtr,Sk.ch);
         lmx0=Sk.nm;
         k=Sk.nm;
         k-=em;
         k-=ek;
         strcat(vtr,ztr+k);
         strcpy(ztr,vtr);
      }
      else {
         n=0;
         break;
      }
   }
   if(n){
      if(Search_cm2(ztr)){
         delete [] ztr;
         return(1);
      }
   }
   pcx=strchr(ztr+Sh.nm,' ');
   if(!pcx){
      if(Search_cm2(ztr)){
         delete [] ztr;
         return(1);
      }
   }
   else {
      *pcx='\0';
      if(Search_cm2(ztr)){
         j=len;
         sxx=cmstr(strlen(cmax),cmax,str,j);
         if(sxx>4){
            if(str[j]==' ')j++;
            else if(str[j+1]==' ')j+=2;
            ixt=j;
            delete [] ztr;
            return(2);
         }
      }
   }
   j=Sh.nm;
   while(j&&(ztr[j]!=' '))j--;
   if(j<=len)return(0);
   ztr[j]='\0';
   if(Search_cm2(ztr)){
      j=len;
      sxx=cmstr(strlen(cmax),cmax,str,j);
      if(sxx>4){
         if(str[j]==' ')j++;
         else if(str[j+1]==' ')j+=2;
         ixt=j;
         delete [] ztr;
         return(2);
      }
   }

   if(lmax){
      k=Sh.dp;
      i=Sh.nm;
      while(k){
         while(k){
            n=Sh.tn[k-1];
            if(!*(cf+ad[n]+Sh.mt[k-1]))break;
            i-=Sh.mt[k-1];
            k--;
         }
         if(k){
            strncpy(cmax,ztr,i);
            cmax[i]='\0';
            j=len;
            sxx=cmstr(i,cmax,str,j);
            if(sxx>4){
               if(str[j]==' ')j++;
               ixt=j;
               delete [] ztr;
               return(2);
            }
            i-=Sh.mt[k-1];
            k--;
         }
      }
      return(0);
   }
   return(0);
}

long CMap::cmstr(long i,const char *zmax,const char *str,long &j){
   long len=j,jp,j1=0,j2=0,j3=0,k,lim;
   long im=(6<i)?(i-6):0;
   const char *pch=zmax+im;
   long jm=(3<im)?(im-3):0;
   long lch=i-im;
   long *sx=new long[lch];
   long sci,sco=0;

   for(j=0;j<lch;j++)sx[j]=0;
   for(j=jm;j<jm+6;j++){
      sci=0;
      k=0;
      lim=(j+lch<=len)?lch:len-j;
      while(k<lim){
         if(pch[k]==str[j+k]){sx[k]=3;j1=j+k;}
         if(sx[k]){sci++;sx[k]--;}
         k++;
      }
      if(sci>sco){
         sco=sci;
         jp=j1;
         jp=(j2>jp)?j2:jp;
         jp=(j3>jp)?j3:jp;
      }
      j3=j2;
      j2=j1;
      j1=0;    
   }
   j=jp+1;
   delete [] sx;
   return(sco);
}

long CMap::Search_spl(const char *str,long lm){
   long len=strlen(str);
   State Sk(3*len);
   if(len<6)return(0);
   sfind(str,Sk);
   if(cfind(str,Sk,lm)){
      Sk.ch[ixt]='\0';
      if(Search_cm2(Sk.ch))strcpy(rmax,cmax);
      else strcpy(rmax,Sk.ch);
      Sk.ch[ixt]='|';
      if(Search_cm2(Sk.ch+ixt))strcat(rmax,cmax);
      else strcat(rmax,Sk.ch+ixt);
      strcpy(cmax,rmax);
      fmax=100.0;
      return(1);
   }
   else return(0);
}

long CMap::Search_spg(const char *str,long lm){
   long len=strlen(str);
   State Sk(3*len);
   if(len<6)return(0);
   sfind(str,Sk);
   if(cfind2(str,Sk,lm)){
      if(!Search_cm2(Sk.ch))strcpy(cmax,Sk.ch);
      fmax=100.0;
      return(1);
   }
   else return(0);
}

//Extract

void CMap::Extract_sng(const char *str){
   long n;
   irp=0;
   if(n=Search_cmb(str)){
      if(n==1||n>=4)return;
      else {
         strcpy(cmax,str);
         return;
      }
   }
   else if((strlen(str)>11)&&Search_dep(str))return;
   else if(Search_spl(str,1))return;
   else strcpy(cmax,str);
}

void CMap::Extract_par(const char *str,long l1,long l2){
   long n,m,k,f1,f2,fp;
   char *pch=new char[l1+l2+3],cnam[max_str];
   strcpy(pch,str);

   if(l1<3){
      if(l2<3)irp=0;
      else if(l2<7)irp=1;
      else irp=2;
   }
   else if(l1<7){
      if(l2<3)irp=1;
      else irp=2;
   }
   else {
      if(l2<7)irp=2;
      else irp=3;
   }

   pch[l1+1]='\0';
   if((k=find(pch))){
      f1=freq[k-1];
   }
   else if(l1<=4)f1=1000;
   else f1=0;
   pch[l1+1]='|';
   if((k=find(pch+l1+1))){
      f2=freq[k-1];
   }
   else if(l2<=4)f2=1000;
   else f2=0;
   pch[l1+1]=' ';
   f1=(f1<f2)?f1:f2;

   k=find(str);
   if(k)fp=freq[k-1];
   else fp=0;
   if((fp>5)&&(f1>500)){
      strcpy(cmax,str);
      return;
   }
   else if(((l1<=4)||(l2<=4))&&fp&&(f1>50)){
      strcpy(cmax,str);
      return;
   }

   n=Search_cm1(str);
   if(n&&((l1>2)&&(l2>2))){
      k=find(cmax);
      fp=freq[k-1];
      if(fp<f1)n=0;
   }

   if((!n)&&(f1>99)){
      n=2;
      pch[l1+1]='|';
      Extract_sng(pch+l1+1);
      strcpy(cnam,cmax);
      pch[l1+1]='\0';
      Extract_sng(pch);
      strcat(cmax,cnam);
   }
   if(!n){
      n=Search_spl(str,500);
   }
   if(l1+l2>19){
      if(((!irp)||(!f1))&&(!n)){
         pch[l1+1]=' ';
         n=Search_dep(pch);
      }
   }
   if(!n){
      pch[l1+1]='|';
      Extract_sng(pch+l1+1);
      strcpy(cnam,cmax);
      pch[l1+1]='\0';
      Extract_sng(pch);
      strcat(cmax,cnam);
   }
}
 
long CMap::Extract_lng(const char *str,long l1,long l2){
   long n,m,k,f1,f2,fp;
   char *pch=new char[l1+l2+3],cnam[max_str];
   strcpy(pch,str);

   if(l1<3){
      if(l2<3)irp=0;
      else if(l2<7)irp=1;
      else irp=2;
   }
   else if(l1<7){
      if(l2<3)irp=1;
      else irp=2;
   }
   else {
      if(l2<7)irp=2;
      else irp=3;
   }

   pch[l1+1]='\0';
   if((k=find(pch))){
      f1=freq[k-1];
   }
   else if(l1<=4)f1=1000;
   else f1=0;
   pch[l1+1]='|';
   if((k=find(pch+l1+1))){
      f2=freq[k-1];
   }
   else if(l2<=4)f2=1000;
   else f2=0;
   pch[l1+1]=' ';
   f1=(f1<f2)?f1:f2;

   k=find(str);
   if(k)fp=freq[k-1];
   else fp=0;
   if((fp>5)&&(f1>500)){
      strcpy(cmax,str);
      return(1);
   }
   else if(((l1<=4)||(l2<=4))&&fp&&(f1>50)){
      strcpy(cmax,str);
      return(1);
   }

   n=Search_cm1(str);
   if(n&&((l1>2)&&(l2>2))){
      k=find(cmax);
      fp=freq[k-1];
      if(fp<f1)n=0;
   }

   if((l1+l2>19)&&!n){
      n=Search_dep(str);
   }
   return(n);
}

long CMap::Extract_qst(const char *pre,const char *str){
   char cpam[max_str];
   strcpy(cpam,pre);
   long n,i=strlen(cpam);
   State Sk(max_str);

   cpam[i++]=' ';
   cpam[i]='\0';
   strcpy(cpam+i,str);

   if(n=sfind(cpam,Sk)){
      strcpy(cmax,cpam);
      return(1);
   }
   else {
      if(n=Search_qst(cpam,Sk))return(n);
      else return(0);
   }
}

long CMap::Extract_fnl(const char *str,long l1,long l2){
   long n,m,k,f1,f2,fp;
   char *pch=new char[l1+l2+3],cnam[max_str];
   strcpy(pch,str);

   if(l1<3){
      if(l2<3)irp=0;
      else if(l2<7)irp=1;
      else irp=2;
   }
   else if(l1<7){
      if(l2<3)irp=1;
      else irp=2;
   }
   else {
      if(l2<7)irp=2;
      else irp=3;
   }

   pch[l1+1]='\0';
   if((k=find(pch))){
      f1=freq[k-1];
   }
   else if(l1<=4)f1=1000;
   else f1=0;
   pch[l1+1]='|';
   if((k=find(pch+l1+1))){
      f2=freq[k-1];
   }
   else if(l2<=4)f2=1000;
   else f2=0;
   pch[l1+1]=' ';
   f1=(f1<f2)?f1:f2;

   k=find(str);
   if(k)fp=freq[k-1];
   else fp=0;
   if((fp>5)&&(f1>500)){
      strcpy(cmax,str);
      return(1);
   }
   else if(((l1<=4)||(l2<=4))&&fp&&(f1>50)){
      strcpy(cmax,str);
      return(1);
   }

   n=Search_cmg(str);
   if(n&&((l1>2)&&(l2>2))){
      k=find(cmax);
      fp=freq[k-1];
      if(fp<f1)n=0;
   }

   if(n){
      if(n==1||n==4)return(1);
      else if(n==5)return(3);
      else {
         strcpy(cmax,str);
         return(2);
      }
   }
   else if((l1+l2>19)&&(irp==3)&&Search_dep(str))return(1);
   else if((n=check_sngl(str,l1,l2))<250){
      if(n)n=2*n;
      if(n<10)n=10;
      if(Search_spg(str,n))return(3);
      else return(0);
   }
   else return(0);
}

long CMap::Rat_check(const char *str,const char *pch){
   long i,j,k,m;
   long lsn=strlen(str);
   long lbn=strlen(pch);
   i=j=k=m=0;

   while((i<lsn)&&(j<lbn)){
      if(str[i]==pch[j]){i++;j++;}
      else if(str[i+1]==pch[j]){
         if(str[i]==pch[j+1]){i+=2;j+=2;}
         else {k++;i++;}
      }
      else if(str[i]==pch[j+1]){k++;j++;}
      else if(str[i+1]==pch[j+1]){k++;i++;j++;}
      else if(str[i+1]&&(str[i+2]==pch[j])){k++;i++;}
      else if(pch[j+1]&&(str[i]==pch[j+2])){k++;j++;}
      else {k++;i++;j++;}
      if((pch[j]==' ')||(pch[j]=='\0')){
         if((j-m<7)&&(k>1))return(0);
         if(k>2)return(0);
         m=j+1;
         k=0;
      }
   }
   if((j<lbn)&&(i==lsn)){
      while((j<lbn)&&(pch[j]!=' ')&&(pch[j]!='\0')){k++;j++;}
      if((j-m<7)&&(k>1))return(0);
      if(k>2)return(0);
   }
   else if((i<lsn)&&(j==lbn)){
      k+=lsn-i;
      if((j-m<7)&&(k>1))return(0);
      if(k>2)return(0);
   }
   return(1);
}

long CMap::RBeg(const char *str,const char *pch){
   long flag=0;
   if(str[1]!=pch[1])flag++;
   if(str[2]!=pch[2])flag++;
   if(str[3]!=pch[3])flag++;
   if(flag>1)return(0);
   else return(1);
}

//State

State::State(long n){
   dp=nm=0;
   tn=new long[n];
   mt=new long[n];
   ch=new char[n];
}

State::~State(void){
   delete [] tn;
   delete [] mt;
   delete [] ch;
}

void State::view(void){
   long i,j,k=0;

   cout << tn[0] << "|";
   for(j=k;j<k+mt[0]+1;j++){
      if(ch[j]!=' ')cout << ch[j];
      else cout << '_';
   }
   k+=mt[0]+1;
   cout << endl;
   for(i=1;i<dp;i++){
      cout << tn[i] << "|";
      for(j=k;j<k+mt[i];j++){
         if(ch[j]!=' ')cout << ch[j];
         else cout << '_';
      }
      k+=mt[i];
      cout << endl;
   }
}

void State::Save(Level &Lv){
   Lv.tn=tn[dp-1];
   Lv.mt=mt[dp-1];
   Lv.dp=dp;
   Lv.nm=nm;
}

void State::Reset(Level &Lv){
   dp=Lv.dp;
   nm=Lv.nm;
   tn[dp-1]=Lv.tn;
   mt[dp-1]=Lv.mt;
}

void State::Copy(State &St){
   long i;
   dp=St.dp;
   nm=St.nm;
   for(i=0;i<dp;i++){
      tn[i]=St.tn[i];
      mt[i]=St.mt[i];
   }
   for(i=0;i<nm;i++)ch[i]=St.ch[i];
}

//Level
Level::Level(void){
}

Level::~Level(void){
}

//Manager

Manager::Manager(void){
}

Manager::Manager(CMap *pMqS,CMap *pMqT,CMap *pMqP){
   pMS=pMqS;
   pMT=pMqT;
   pMP=pMqP;
}

Manager::~Manager(void){}

void Manager::gopen_map(void){
   pMS->gopen_cmap();
   pMT->gopen_cmap();
   pMP->gopen_cmap();
}

void Manager::gopen_Man_copy(Manager *pMan){
   pMS=new CMap;
   pMT=new CMap;
   pMP=new CMap;

   pMS->gopen_CMap_copy(pMan->pMS);
   pMT->gopen_CMap_copy(pMan->pMT);
   pMP->gopen_CMap_copy(pMan->pMP);
}

void Manager::Reset(void){
   cmax[0]='\0';
   imax=0;
}

void Manager::Extract(const char *str){
   long m,n,i=0,j,k,jz,mz;
   char *pch,*pcc,*ptt,c,cnam[max_str];
   cnam[0]='|';
   pch=cnam+1;

   while((pch[i]=c=str[i])&&(c!=' '))i++;
   if(!c){
      if(i>2){
         pMS->Extract_sng(cnam);
         strcpy(cmax+imax,pMS->cmax);
         imax+=strlen(pMS->cmax);
      }
      else {
         strcpy(cmax+imax,cnam);
         imax+=strlen(cnam);
      }
   }
   else if(c){
      j=i+1;
      while((pch[j]=c=str[j])&&(c!=' '))j++;
      if(!c){
         pMS->Extract_par(cnam,i,j-i-1);
         strcpy(cmax+imax,pMS->cmax);
         imax+=strlen(pMS->cmax);
      }
      else {
         pch[j]='\0';
         n=pMT->Extract_lng(cnam,i,j-i-1);
         if(n){
            pcc=cnam+j+3;
            strcpy(pcc,pMT->cmax);
            m=pMP->Extract_qst(pcc,str+j+1);
            if(m>0&&!Rat_check(&(pMP->cmax[1]),str))m=0;
            if(m){
               strcpy(cmax+imax,pMP->cmax);
               imax+=strlen(pMP->cmax);
               if(m==2)Extract(str+pMP->ixt+j-strlen(pcc));
            }
            else if(!strchr(pcc,' ')){
               m=pMS->Extract_qst(pcc,str+j+1);
               if(m>0&&!Rat_check(&(pMS->cmax[1]),str))m=0;
               if(m){
                  strcpy(cmax+imax,pMS->cmax);
                  imax+=strlen(pMS->cmax);
                  if(m==2)Extract(str+pMS->ixt+j-strlen(pcc));
               }
            }
         }
         if(!n||!m){
            n=pMS->Extract_fnl(cnam,i,j-i-1);
            if((n==1)||(n==2)){
               strcpy(cmax+imax,pMS->cmax);
               imax+=strlen(pMS->cmax);
               Extract(str+j+1);
            }
            else if(n==3){
               strcpy(cmax+imax,pMS->cmax);
               imax+=strlen(pMS->cmax);
               Extract(str+pMS->ixt);
            }
            else {
               pch[i]='\0';
               pMS->Extract_sng(cnam);
               strcpy(cmax+imax,pMS->cmax);
               imax+=strlen(pMS->cmax);
               Extract(str+i+1);
            }
         }
      }
   }
}

long Manager::Rat_check(const char *str,const char *pch){
   long i,j,k,m;
   long lsn=strlen(str);
   long lbn=strlen(pch);
   i=j=k=m=0;

   while((i<lsn)&&(j<lbn)){
      if(str[i]==pch[j]){i++;j++;}
      else if(str[i+1]==pch[j]){
         if(str[i]==pch[j+1]){i+=2;j+=2;}
         else {k++;i++;}
      }
      else if(str[i]==pch[j+1]){k++;j++;}
      else if(str[i+1]==pch[j+1]){k++;i++;j++;}
      else if(str[i+1]&&(str[i+2]==pch[j])){k++;i++;}
      else if(pch[j+1]&&(str[i]==pch[j+2])){k++;j++;}
      else {k++;i++;j++;}
      if((pch[j]==' ')||(pch[j]=='\0')){
         if((j-m<7)&&(k>1))return(0);
         if(k>2)return(0);
         m=j+1;
         k=0;
      }
   }
   if((j<lbn)&&(i==lsn)){
      while((j<lbn)&&(pch[j]!=' ')&&(pch[j]!='\0')){k++;j++;}
      if((j-m<7)&&(k>1))return(0);
      if(k>2)return(0);
   }
   else if((i<lsn)&&(j==lbn)){
      k+=lsn-i;
      if((j-m<7)&&(k>1))return(0);
      if(k>2)return(0);
   }
   return(1);
}

}
