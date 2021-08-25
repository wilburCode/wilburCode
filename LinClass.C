#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "Postg.h"
#include "LinClass.h"
using namespace std;
namespace iret {

LinClass::LinClass(const char *nam) : FBase("linset",nam){
}

LinClass::~LinClass(void){}

Count *LinClass::Count_dbin(const char *post_nam,double *wg){
   char *pch,cnam[10000];
   long i,k;
   FBase Fb("postset",post_nam);
   ifstream *pfin=Fb.get_Istr("s");
   Count *pCt=new Count;

   i=0;
   while(pfin->getline(cnam,10000,'\n')){
      if(wg[i]){
         pCt->add_count2(cnam,i);
      }
      mark(++i,10000,"terms");
   }
   dst_Istr(pfin);

   if(pflag)cout << "Terms in use: " << pCt->cnt_key << endl;
   return(pCt);
}

Count *LinClass::Count_bnby(const char *slice_nam,float *wg){
   char *pch,cnam[10000];
   long i,k,nwrds,nid;
   FBase Fs("slice",slice_nam);
   FBase Fb("blistset",slice_nam);
   Count *pCt=new Count;

   ifstream *pfin=Fb.get_Istr("z");
   *pfin >> nid;
   dst_Istr(pfin); 
   long *inv=(long*)Fb.get_Mmap("inv");
   pfin=Fs.get_Istr("s");
   long *adds=(long*)Fs.get_Mmap("as");

   for(i=0;i<nid;i++){
      if(wg[i]){
         if(k=inv[i]){
            pfin->seekg(adds[k-1],ios::beg);
            pfin->getline(cnam,10000,'\n');
            pCt->add_count2(cnam,i);
         }
         else {cout << "Error in weights!" << endl;exit(0);}
      }
      mark(i,10000,"terms");
   }
   dst_Istr(pfin);

   if(pflag)cout << "Terms in use: " << pCt->cnt_key << endl;
   return(pCt);
}

Count *LinClass::Count_cmls(const char *postg_nam,double *wg){
   char *pch,cnam[10000];
   long i,k,nwrd;

   Postg<char> Psg(postg_nam);

   ifstream *pfin=Psg.get_Istr("n",ios::in);
   *pfin >> nwrd;
   Psg.dst_Istr(pfin);
   char *term=Psg.get_Mmap("s");
   long *sddr=(long*)Psg.get_Mmap("sa");

   Count *pCt=new Count;

   for(i=0;i<nwrd;i++){
      if(wg[i]){
         pCt->add_count2(term+sddr[i],i);
      }
      mark(i,10000,"terms");
   }

   if(pflag)cout << "Terms in use: " << pCt->cnt_key << endl;
   return(pCt);
}

void LinClass::create_Hash(Count *pCt,double *wg,int exc){
   long i,k;
   char *pch,cnam[max_str];

   Lhs.change_type("linset");
   strcpy(cnam,name);
   strcat(cnam,"-hsh");
   Lhs.change_name(cnam);
   Lhs.create_htable(*pCt,exc);

   if(pflag)cout << "Hash table completed" << endl;

   Lhs.gopen_htable_map();
   wt=new float[Lhs.nwrds];
   pCt->node_first();
   while(pCt->node_next()){
      pch=pCt->show_str();
      i=pCt->count();
      if(k=Lhs.find(pch)){
         wt[k-1]=(float)wg[i];
      }
      else {cout << "Error!" << endl;exit(0);}
   }
   bin_Writ("wt",Lhs.nwrds*sizeof(float),(char*)wt);
}

void LinClass::create_Hash(Count *pCt,float *wg,int exc){
   long i,k;
   char *pch,cnam[max_str];

   Lhs.change_type("linset");
   strcpy(cnam,name);
   strcat(cnam,"-hsh");
   Lhs.change_name(cnam);
   Lhs.create_htable(*pCt,exc);

   if(pflag)cout << "Hash table completed" << endl;

   Lhs.gopen_htable_map();
   wt=new float[Lhs.nwrds];
   pCt->node_first();
   while(pCt->node_next()){
      pch=pCt->show_str();
      i=pCt->count();
      if(k=Lhs.find(pch)){
         wt[k-1]=(float)wg[i];
      }
      else {cout << "Error!" << endl;exit(0);}
   }
   bin_Writ("wt",Lhs.nwrds*sizeof(float),(char*)wt);
}

void LinClass::save_Thresh(double th){
   ofstream *pfout=get_Ostr("t");
   *pfout << th << endl;
   dst_Ostr(pfout);
}

void LinClass::save_Thresh(float th){
   ofstream *pfout=get_Ostr("t");
   *pfout << thresh << endl;
   dst_Ostr(pfout);
}

void LinClass::gopen_operate(void){
   char cnam[10000];

   if(Exists("t")){
      ifstream *pfin=get_Istr("t");
      *pfin >> thresh;
      dst_Istr(pfin);
   }
   else {
      thresh=0;
      cout << "Threshhold set to 0!" << endl;
   }
   Lhs.change_type("linset");
   strcpy(cnam,name);
   strcat(cnam,"-hsh");
   Lhs.change_name(cnam); 
   Lhs.gopen_htable_map(); 

   wt=(float*)get_Mmap("wt");
}

float LinClass::weight(const char *str){
   long i;
   if(i=Lhs.find(str))return(wt[i-1]);
   else return(0);
}
 
}
