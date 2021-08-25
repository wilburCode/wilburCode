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
#include <map>
#include <runn.h>
#include <Hash.h>

using namespace std;
namespace iret {

Hash::Hash(void) : FBase("hshset","null"){
   own_hash_mem = false;
}

Hash::Hash(const char *nam) : FBase("hshset",nam){
   own_hash_mem = false;
}

Hash::Hash(int n,const char *nam) : FBase("hshset",n,nam){
   own_hash_mem = false;
}

Hash::~Hash(){
   if(own_hash_mem){
      if(num_file>-1)gclose_htable_map(num_file);
      else gclose_htable_map();
   }
}

void Hash::SetMem(Hash &Hsh){
   this->strmap=Hsh.strmap; //Holds the bit map.
   this->addr  =Hsh.addr; //Holds the offsets to strmap.
   this->nwrds =Hsh.nwrds; //Number of words.
   this->tnum  =Hsh.tnum; //Truncation number, size of har.
   this->harr  =Hsh.harr; //Holds hash array.
   this->farr  =Hsh.farr; //Holds the hash coefficients.
   this->px0   =Hsh.px0;
   this->px1   =Hsh.px1;
   this->px2   =Hsh.px2;
   this->px3   =Hsh.px3;
   this->px4   =Hsh.px4;
   this->px5   =Hsh.px5;
   this->px6   =Hsh.px6;
   this->px7   =Hsh.px7;
   this->px8   =Hsh.px8;
   this->px9   =Hsh.px9;
   this->px10  =Hsh.px10;
   this->px11  =Hsh.px11;
}

void Hash::create_htable(List &Lst,int excess){
   char cnam[max_str],*cptr,*uptr;
   int u,len;
   long ct,i,j,k;
   ofstream *pfout;

   nwrds=Lst.cnt_key;
   ct=nwrds;
   tnum=1;
   u=0;
   while(ct=ct/2){tnum*=2;u++;}
   if(u>30){cout << "Error in size, " << u << endl;exit(0);}
   i=0;
   while((u<32)&&(i<excess)){tnum*=2;u++;i++;}
   tnum--;
   harr=new long[tnum+2];
   for(ct=0;ct<tnum+2;ct++)harr[ct]=0;

   farr=new long[3072];
   ct=1;
   for(i=0;i<3072;i++){
      farr[i]=ct=(ct*331)&tnum;
   }
  
   long *pc0=farr+128,*pc1=farr+384,*pc2=farr+640;
   long *pc3=farr+896,*pc4=farr+1152,*pc5=farr+1408;
   long *pc6=farr+1664,*pc7=farr+1920,*pc8=farr+2176;
   long *pc9=farr+2432,*pc10=farr+2688,*pc11=farr+2944;
   
   Lst.node_first();
   while(Lst.node_next()){
      cptr=Lst.show_str();
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      (harr[ct&tnum])++;
   }

   //Set start points in harr.
   k=0;
   for(i=0;i<tnum+2;i++){
      j=harr[i];
      harr[i]=k;
      k+=j;
   }
   if(k!=nwrds){cout << "Error in summing!" << endl;exit(0);}

   //Write out harr.
   bin_Writ("ha",(tnum+2)*sizeof(long),(char*)harr);

   //Set addresses
   char **addt=new char*[nwrds];
   Lst.node_first();
   while(Lst.node_next()){
      uptr=cptr=Lst.show_str();
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      k=ct&tnum;
      addt[harr[k]]=uptr;
      (harr[k])++;
   }

   //Write out string file
   pfout=get_Ostr("str");
   k=0;
   for(i=0;i<nwrds;i++){
      *pfout << addt[i] << ends;
      len=strlen((char*)addt[i])+1;
      addt[i]=(char*)k;
      k+=len;
   }
   dst_Ostr(pfout);

   //Write out addr file
   bin_Writ("ad",nwrds*sizeof(long),(char*)addt);
   delete [] addt;
   addt=NULL;

   //Write out counts
   pfout=get_Ostr("nm");
   *pfout << nwrds << " " << tnum << " " << k << endl;
   dst_Ostr(pfout);
   delete [] harr;
   delete [] farr;
   harr=NULL;
   farr=NULL;
}

void Hash::create_htable(strMap &Mp,int excess){
   char cnam[max_str];
   const char *cptr,*uptr;
   int u,len;
   long ct,i,j,k;
   ofstream *pfout;

   nwrds=Mp.size();
   ct=nwrds;
   tnum=1;
   u=0;
   while(ct=ct/2){tnum*=2;u++;}
   if(u>30){cout << "Error in size, " << u << endl;exit(0);}
   i=0;
   while((u<32)&&(i<excess)){tnum*=2;u++;i++;}
   tnum--;
   harr=new long[tnum+2];
   for(ct=0;ct<tnum+2;ct++)harr[ct]=0;

   farr=new long[3072];
   ct=1;
   for(i=0;i<3072;i++){
      farr[i]=ct=(ct*331)&tnum;
   }
  
   long *pc0=farr+128,*pc1=farr+384,*pc2=farr+640;
   long *pc3=farr+896,*pc4=farr+1152,*pc5=farr+1408;
   long *pc6=farr+1664,*pc7=farr+1920,*pc8=farr+2176;
   long *pc9=farr+2432,*pc10=farr+2688,*pc11=farr+2944;
   
   typename strMap::iterator p=Mp.begin();
   typename strMap::iterator q=Mp.end();
   while(p!=q){
      cptr=p->first;
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      (harr[ct&tnum])++;
      p++;
   }

   //Set start points in harr.
   k=0;
   for(i=0;i<tnum+2;i++){
      j=harr[i];
      harr[i]=k;
      k+=j;
   }
   if(k!=nwrds){cout << "Error in summing!" << endl;exit(0);}

   //Write out harr.
   bin_Writ("ha",(tnum+2)*sizeof(long),(char*)harr);

   //Set addresses
   const char **addt=new const char*[nwrds];
   p=Mp.begin();
   q=Mp.end();
   while(p!=q){
      uptr=cptr=p->first;
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      k=ct&tnum;
      addt[harr[k]]=uptr;
      (harr[k])++;
      p++;
   }

   //Write out string file
   pfout=get_Ostr("str");
   k=0;
   for(i=0;i<nwrds;i++){
      *pfout << addt[i] << ends;
      len=strlen((char*)addt[i])+1;
      addt[i]=(char*)k;
      k+=len;
   }
   dst_Ostr(pfout);

   //Write out addr file
   bin_Writ("ad",nwrds*sizeof(long),(char*)addt);
   delete [] addt;
   addt=NULL;

   //Write out counts
   pfout=get_Ostr("nm");
   *pfout << nwrds << " " << tnum << " " << k << endl;
   dst_Ostr(pfout);
   delete [] harr;
   delete [] farr;
   harr=NULL;
   farr=NULL;
}

//In memory model intended for small sets
void Hash::create_htableM(List &Lst,int excess){
   char cnam[max_str],*cptr,*uptr;
   int u,len;
   long ct,i,j,k,*barr;
   ofstream *pfout;

   nwrds=Lst.cnt_key;
   ct=nwrds;
   tnum=1;
   u=0;
   while(ct=ct/2){tnum*=2;u++;}
   if(u>30){cout << "Error in size, " << u << endl;exit(0);}
   i=0;
   while((u<32)&&(i<excess)){tnum*=2;u++;i++;}
   tnum--;
   harr=new long[tnum+2];
   barr=new long[tnum+2];
   for(ct=0;ct<tnum+2;ct++)harr[ct]=0;

   farr=new long[3072];
   ct=1;
   for(i=0;i<3072;i++){
      farr[i]=ct=(ct*331)&tnum;
   }
  
   px0=farr+128,px1=farr+384,px2=farr+640;
   px3=farr+896,px4=farr+1152,px5=farr+1408;
   px6=farr+1664,px7=farr+1920,px8=farr+2176;
   px9=farr+2432,px10=farr+2688,px11=farr+2944;
  
   Lst.node_first();
   while(Lst.node_next()){
      cptr=Lst.show_str();
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(px0+u);
                    break;
            case 1: ct+=*(px1+u);
                    break;
            case 2: ct+=*(px2+u);
                    break;
            case 3: ct+=*(px3+u);
                    break;
            case 4: ct+=*(px4+u);
                    break;
            case 5: ct+=*(px5+u);
                    break;
            case 6: ct+=*(px6+u);
                    break;
            case 7: ct+=*(px7+u);
                    break;
            case 8: ct+=*(px8+u);
                    break;
            case 9: ct+=*(px9+u);
                    break;
            case 10: ct+=*(px10+u);
                    break;
            case 11: ct+=*(px11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      (harr[ct&tnum])++;
   }

   //Set start points in harr.
   k=0;
   for(i=0;i<tnum+2;i++){
      j=harr[i];
      barr[i]=harr[i]=k;
      k+=j;
   }
   if(k!=nwrds){cout << "Error in summing!" << endl;exit(0);}

   //Set addresses
   len=0;
   char **addt=new char*[nwrds];
   Lst.node_first();
   while(Lst.node_next()){
      uptr=cptr=Lst.show_str();
      len+=strlen(uptr)+1;
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(px0+u);
                    break;
            case 1: ct+=*(px1+u);
                    break;
            case 2: ct+=*(px2+u);
                    break;
            case 3: ct+=*(px3+u);
                    break;
            case 4: ct+=*(px4+u);
                    break;
            case 5: ct+=*(px5+u);
                    break;
            case 6: ct+=*(px6+u);
                    break;
            case 7: ct+=*(px7+u);
                    break;
            case 8: ct+=*(px8+u);
                    break;
            case 9: ct+=*(px9+u);
                    break;
            case 10: ct+=*(px10+u);
                    break;
            case 11: ct+=*(px11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      k=ct&tnum;
      addt[barr[k]]=uptr;
      (barr[k])++;
   }
   strmap=new char[len];

   //Set up string array
   k=0;
   for(i=0;i<nwrds;i++){
      len=strlen((char*)addt[i])+1;
      strcpy(strmap+k,addt[i]);
      addt[i]=(char*)k;
      k+=len;
   }
   addr=(long*)addt;
   delete [] barr;
}

void Hash::create_htable(int mz,List &Lst,int excess){
   char cnam[max_str],*cptr,*uptr;
   int u,len;
   long ct,i,j,k;
   ofstream *pfout;

   nwrds=Lst.cnt_key;
   ct=nwrds;
   tnum=1;
   u=0;
   while(ct=ct/2){tnum*=2;u++;}
   if(u>30){cout << "Error in size, " << u << endl;exit(0);}
   i=0;
   while((u<32)&&(i<excess)){tnum*=2;u++;i++;}
   tnum--;
   harr=new long[tnum+2];
   for(ct=0;ct<tnum+2;ct++)harr[ct]=0;

   farr=new long[3072];
   ct=1;
   for(i=0;i<3072;i++){
      farr[i]=ct=(ct*331)&tnum;
   }
  
   long *pc0=farr+128,*pc1=farr+384,*pc2=farr+640;
   long *pc3=farr+896,*pc4=farr+1152,*pc5=farr+1408;
   long *pc6=farr+1664,*pc7=farr+1920,*pc8=farr+2176;
   long *pc9=farr+2432,*pc10=farr+2688,*pc11=farr+2944;
   
   Lst.node_first();
   while(Lst.node_next()){
      cptr=Lst.show_str();
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      (harr[ct&tnum])++;
   }

   //Set start points in harr.
   k=0;
   for(i=0;i<tnum+2;i++){
      j=harr[i];
      harr[i]=k;
      k+=j;
   }
   if(k!=nwrds){cout << "Error in summing!" << endl;exit(0);}

   //Write out harr.
   bin_Writ(mz,"ha",(tnum+2)*sizeof(long),(char*)harr);

   //Set addresses
   char **addt=new char*[nwrds];
   Lst.node_first();
   while(Lst.node_next()){
      uptr=cptr=Lst.show_str();
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      k=ct&tnum;
      addt[harr[k]]=uptr;
      (harr[k])++;
   }

   //Write out string file
   pfout=get_Ostr(mz,"str");
   k=0;
   for(i=0;i<nwrds;i++){
      *pfout << addt[i] << ends;
      len=strlen((char*)addt[i])+1;
      addt[i]=(char*)k;
      k+=len;
   }
   dst_Ostr(pfout);

   //Write out addr file
   bin_Writ(mz,"ad",nwrds*sizeof(long),(char*)addt);
   delete [] addt;
   addt=NULL;

   //Write out counts
   pfout=get_Ostr(mz,"nm");
   *pfout << nwrds << " " << tnum << " " << k << endl;
   dst_Ostr(pfout);
   delete [] harr;
   delete [] farr;
   harr=NULL;
   farr=NULL;
}

void Hash::create_htable(int mz,strMap &Mp,int excess){
   char cnam[max_str];
   const char *cptr,*uptr;
   int u,len;
   long ct,i,j,k;
   ofstream *pfout;

   nwrds=Mp.size();
   ct=nwrds;
   tnum=1;
   u=0;
   while(ct=ct/2){tnum*=2;u++;}
   if(u>30){cout << "Error in size, " << u << endl;exit(0);}
   i=0;
   while((u<32)&&(i<excess)){tnum*=2;u++;i++;}
   tnum--;
   harr=new long[tnum+2];
   for(ct=0;ct<tnum+2;ct++)harr[ct]=0;

   farr=new long[3072];
   ct=1;
   for(i=0;i<3072;i++){
      farr[i]=ct=(ct*331)&tnum;
   }
  
   long *pc0=farr+128,*pc1=farr+384,*pc2=farr+640;
   long *pc3=farr+896,*pc4=farr+1152,*pc5=farr+1408;
   long *pc6=farr+1664,*pc7=farr+1920,*pc8=farr+2176;
   long *pc9=farr+2432,*pc10=farr+2688,*pc11=farr+2944;
   
   typename strMap::iterator p=Mp.begin();
   typename strMap::iterator q=Mp.end();
   while(p!=q){
      cptr=p->first;
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      (harr[ct&tnum])++;
      p++;
   }

   //Set start points in harr.
   k=0;
   for(i=0;i<tnum+2;i++){
      j=harr[i];
      harr[i]=k;
      k+=j;
   }
   if(k!=nwrds){cout << "Error in summing!" << endl;exit(0);}

   //Write out harr.
   bin_Writ(mz,"ha",(tnum+2)*sizeof(long),(char*)harr);

   //Set addresses
   const char **addt=new const char*[nwrds];
   p=Mp.begin();
   q=Mp.end();
   while(p!=q){
      uptr=cptr=p->first;
      ct=0;
      i=0;
      while(u=*(cptr++)){
         switch(i){
            case 0: ct+=*(pc0+u);
                    break;
            case 1: ct+=*(pc1+u);
                    break;
            case 2: ct+=*(pc2+u);
                    break;
            case 3: ct+=*(pc3+u);
                    break;
            case 4: ct+=*(pc4+u);
                    break;
            case 5: ct+=*(pc5+u);
                    break;
            case 6: ct+=*(pc6+u);
                    break;
            case 7: ct+=*(pc7+u);
                    break;
            case 8: ct+=*(pc8+u);
                    break;
            case 9: ct+=*(pc9+u);
                    break;
            case 10: ct+=*(pc10+u);
                    break;
            case 11: ct+=*(pc11+u);
                     i-=12;
                    break;
         }
         i++;
      }
      k=ct&tnum;
      addt[harr[k]]=uptr;
      (harr[k])++;
      p++;
   }

   //Write out string file
   pfout=get_Ostr(mz,"str");
   k=0;
   for(i=0;i<nwrds;i++){
      *pfout << addt[i] << ends;
      len=strlen((char*)addt[i])+1;
      addt[i]=(char*)k;
      k+=len;
   }
   dst_Ostr(pfout);

   //Write out addr file
   bin_Writ(mz,"ad",nwrds*sizeof(long),(char*)addt);
   delete [] addt;
   addt=NULL;

   //Write out counts
   pfout=get_Ostr(mz,"nm");
   *pfout << nwrds << " " << tnum << " " << k << endl;
   dst_Ostr(pfout);
   delete [] harr;
   delete [] farr;
   harr=NULL;
   farr=NULL;
}

void Hash::gopen_htable_map(void){
   char cnam[max_str],*cptr;
   int fld;
   long ct,asize,i;
   
   ifstream *pfin=get_Istr("nm");
   *pfin >> nwrds >> tnum >> asize;
   dst_Istr(pfin);

   harr=(long*)get_Mmap("ha");
   addr=(long*)get_Mmap("ad");
   strmap=get_Mmap("str");

   farr=new long[3072];
   ct=1;
   for(i=0;i<3072;i++){
      farr[i]=ct=(ct*331)&tnum;
   }
  
   px0=farr+128,px1=farr+384,px2=farr+640;
   px3=farr+896,px4=farr+1152,px5=farr+1408;
   px6=farr+1664,px7=farr+1920,px8=farr+2176;
   px9=farr+2432,px10=farr+2688,px11=farr+2944;
   own_hash_mem = true;
   num_file = -1;
}

void Hash::gopen_htable_map(int mz){
   char cnam[max_str],*cptr;
   int fld;
   long ct,asize,i;
   
   ifstream *pfin=get_Istr(mz,"nm");
   *pfin >> nwrds >> tnum >> asize;
   dst_Istr(pfin);

   harr=(long*)get_Mmap(mz,"ha");
   addr=(long*)get_Mmap(mz,"ad");
   strmap=get_Mmap(mz,"str");

   farr=new long[3072];
   ct=1;
   for(i=0;i<3072;i++){
      farr[i]=ct=(ct*331)&tnum;
   }
  
   px0=farr+128,px1=farr+384,px2=farr+640;
   px3=farr+896,px4=farr+1152,px5=farr+1408;
   px6=farr+1664,px7=farr+1920,px8=farr+2176;
   px9=farr+2432,px10=farr+2688,px11=farr+2944;
   own_hash_mem = true;
   num_file = mz;
}

long Hash::find(const char *str){
   long ct=0,i=0,k;
   int ic;
   const char *utr=str;
   while(ic=*(utr++)){
      switch(i){
         case 0: ct+=*(px0+ic);
                 break;
         case 1: ct+=*(px1+ic);
                 break;
         case 2: ct+=*(px2+ic);
                 break;
         case 3: ct+=*(px3+ic);
                 break;
         case 4: ct+=*(px4+ic);
                 break;
         case 5: ct+=*(px5+ic);
                 break;
         case 6: ct+=*(px6+ic);
                 break;
         case 7: ct+=*(px7+ic);
                 break;
         case 8: ct+=*(px8+ic);
                 break;
         case 9: ct+=*(px9+ic);
                 break;
         case 10: ct+=*(px10+ic);
                 break;
         case 11: ct+=*(px11+ic);
                  i-=12;
                 break;
      }
      i++;
   }
   k=ct&tnum;
   ct=harr[k+1];
   i=harr[k];
   switch(ct-i){
      case 0: return(0);
              break;
      case 1: if(!strcmp(str,strmap+addr[i]))return(i+1);
              else return(0);
              break;
      case 2: ic=strcmp(str,strmap+addr[i]);
              if(ic>0){
                 if(!strcmp(str,strmap+addr[i+1]))return(i+2);
                 else return(0);
              }
              else if(ic<0)return(0);
              else return(i+1);
              break;
      default: ic=strcmp(str,strmap+addr[i]);
               if(ic<0)return(0);
               else if(!ic)return(i+1);
               ct--;
               ic=strcmp(str,strmap+addr[ct]);
               if(ic>0)return(0);
               else if(!ic)return(ct+1);
               while(ct-i>1){
                  k=(ct+i)/2;
                  ic=strcmp(str,strmap+addr[k]);
                  if(ic>0)i=k;
                  else if(ic<0)ct=k;
                  else return(k+1);
               }
               return(0);
   }
}

void Hash::gclose_htable_map(void){
   dst_Mmap("ha",(char*&)harr);
   dst_Mmap("ad",(char*&)addr);
   dst_Mmap("str",strmap);
   delete [] farr;
   harr=NULL;
   farr=NULL;
   addr=NULL;
   strmap=NULL;
   own_hash_mem = false;
}

void Hash::gclose_htable_map(int mz){
   dst_Mmap(mz,"ha",(char*&)harr);
   dst_Mmap(mz,"ad",(char*&)addr);
   dst_Mmap(mz,"str",strmap);
   delete [] farr;
   harr=NULL;
   farr=NULL;
   addr=NULL;
   strmap=NULL;
   own_hash_mem = false;
}

//Chash code

Chash::Chash() : Hash(){
   change_type("cshset");
   own_chash_mem = false;
}

Chash::Chash(const char *str) : Hash(str){
   change_type("cshset");
   own_chash_mem = false;
}

Chash::Chash(int n,const char *str) : Hash(n,str){
   change_type("cshset");
   own_chash_mem = false;
}


void Chash::SetMem(Chash &Chsh){
   Hash::SetMem((Hash&)Chsh);
   this->cnt   =Chsh.cnt;
}

Chash::~Chash(void){
   if(own_chash_mem){
      if(num_file>-1)gclose_ctable_map(num_file);
      else gclose_ctable_map();
   }
}

void Chash::create_ctable(Count &Ct,int excess){
   create_htable(Ct,excess);
   gopen_htable_map();
   long n,i=0;
   long *pct=new long[Ct.cnt_key];
   Ct.node_first();
   while(Ct.node_next()){
      if(n=find(Ct.show_str())){
         pct[n-1]=Ct.count();
      }        
      else {
         cout << "Error in Count tree!" << endl;exit(0);
      }
      mark(++i,10000,"count terms");
   }
   bin_Writ("ct",Ct.cnt_key*sizeof(long),(char*)pct);
   delete [] pct;
   //cnt=(long*)get_Mmap("ct");
   gclose_htable_map();
}

void Chash::create_ctable(strMap &Mp,int excess){
   create_htable(Mp,excess);
   gopen_htable_map();
   long n,i=0;
   long *pct=new long[Mp.size()];
   typename strMap::iterator p=Mp.begin();
   typename strMap::iterator q=Mp.end();
   while(p!=q){
      if(n=find(p->first)){
         pct[n-1]=p->second;
      }        
      else {
         cout << "Error in Map!" << endl;exit(0);
      }
      p++;
      mark(++i,10000,"count terms");
   }
   bin_Writ("ct",Mp.size()*sizeof(long),(char*)pct);
   delete [] pct;
   gclose_htable_map();
}

void Chash::create_ctable_STerm(strMap &Mp,int excess){
   create_htable(Mp,excess);
   gopen_htable_map();
   long n,i=1;
   long *pct=new long[Mp.size()];
   typename strMap::iterator p=Mp.begin();
   typename strMap::iterator q=Mp.end();
   while(p!=q){
      if(n=find(p->first)){
         pct[n-1]=i;
      }        
      else {
         cout << "Error in Map!" << endl;exit(0);
      }
      p++;
      mark(++i,10000,"count terms");
   }
   bin_Writ("ct",Mp.size()*sizeof(long),(char*)pct);
   delete [] pct;
   gclose_htable_map();
}

void Chash::create_ctable(List &Lt,int excess){
   create_htable(Lt,excess);
   gopen_htable_map();
   long n,i=1;
   long *pct=new long[Lt.cnt_key];
   Lt.node_first();
   while(Lt.node_next()){
      if(n=find(Lt.show_str())){
         pct[n-1]=i;
      }
      else {
         cout << "Error in List tree!" << endl;exit(0);
      }
      mark(++i,10000,"count terms");
   }
   bin_Writ("ct",Lt.cnt_key*sizeof(long),(char*)pct);
   delete [] pct;
   cnt=(long*)get_Mmap("ct");
   gclose_htable_map();
}

void Chash::create_ctable(int mz,Count &Ct,int excess){
   create_htable(mz,Ct,excess);
   gopen_htable_map(mz);
   long n,i=0;
   long *pct=new long[Ct.cnt_key];
   Ct.node_first();
   while(Ct.node_next()){
      if(n=find(Ct.show_str())){
         pct[n-1]=Ct.count();
      }        
      else {
         cout << "Error in Count tree!" << endl;exit(0);
      }
      mark(++i,10000,"count terms");
   }
   bin_Writ(mz,"ct",Ct.cnt_key*sizeof(long),(char*)pct);
   delete [] pct;
   gclose_htable_map(mz);
}

void Chash::create_ctable(int mz,strMap &Mp,int excess){
   create_htable(mz,Mp,excess);
   gopen_htable_map();
   long n,i=0;
   long *pct=new long[Mp.size()];
   typename strMap::iterator p=Mp.begin();
   typename strMap::iterator q=Mp.end();
   while(p!=q){
      if(n=find(p->first)){
         pct[n-1]=p->second;
      }
      else {
         cout << "Error in Map!" << endl;exit(0);
      }
      p++;
      mark(++i,10000,"count terms");
   }
   bin_Writ(mz,"ct",Mp.size()*sizeof(long),(char*)pct);
   delete [] pct;
   gclose_htable_map(mz);
}

void Chash::create_ctable(int mz,List &Lt,int excess){
   create_htable(mz,Lt,excess);
   gopen_htable_map(mz);
   long n,i=1;
   long *pct=new long[Lt.cnt_key];
   Lt.node_first();
   while(Lt.node_next()){
      if(n=find(Lt.show_str())){
         pct[n-1]=i;
      }
      else {
         cout << "Error in List tree!" << endl;exit(0);
      }
      mark(++i,10000,"count terms");
   }
   bin_Writ(mz,"ct",Lt.cnt_key*sizeof(long),(char*)pct);
   delete [] pct;
   gclose_htable_map(mz);
}

void Chash::gopen_ctable_map(void){
   gopen_htable_map();
   cnt=(long*)get_Mmap("ct");
   own_chash_mem = true;
}   

void Chash::gopen_ctable_map(int mz){
   gopen_htable_map(mz);
   cnt=(long*)get_Mmap(mz,"ct");
   own_chash_mem = true;
}   

void Chash::gclose_ctable_map(void){
   gclose_htable_map();
   dst_Mmap("ct",(char*&)cnt);
   cnt=NULL;
   own_chash_mem = false;
}   

void Chash::gclose_ctable_map(int mz){
   gclose_htable_map(mz);
   dst_Mmap(mz,"ct",(char*&)cnt);
   cnt=NULL;
   own_chash_mem = false;
}   

long Chash::count(const char *str){
   long n=find(str);
   if(n)return(cnt[n-1]);
   else return(0);
}

//RelateA code

RelateA::RelateA(void) : FBase("relatea","null"){
   name=NULL;
}

RelateA::RelateA(const char *nam) : FBase("relatea",nam){
}

RelateA::~RelateA(){
}

int RelateA::create_bmatrix(BTList &Btl,int exc1,int exc2){
   char cnam[max_str],*cptr;
   long i,j,k,ct,mapsize;

   nwrd1=Btl.cnt_key;
   nwrd2=Btl.lst->cnt_key;
   cout << nwrd1 << " " << nwrd2 << endl;

   //Create Hash for strings set 1
   strcpy(cnam,name);
   strcat(cnam,"-s1");
   Coord1.change_name(cnam);
   Coord1.create_htable(Btl,exc1);
   //Create Hash for strings set 2
   strcpy(cnam,name);
   strcat(cnam,"-s2");
   Coord2.change_name(cnam);
   Coord2.create_htable(*(Btl.lst),exc2);

   //Write out numbers of strings
   get_pathw(cnam,"relatea",name,"nm");
   ofstream fout(cnam,ios::out);
   fout << nwrd1 << " " << nwrd2 << endl;
   fout.close();

   //Map Hash structures
   Coord1.gopen_htable_map();
   Coord2.gopen_htable_map();

   //Make map
   ad=new long[nwrd1+2];
   ad[0]=0;
   ad[nwrd1+1]=0;
   ct=0;
   mapsize=0;
   Btl.node_first();
   while(Btl.node_next()){
      if(i=Coord1.find(Btl.show_str())){
         mapsize+=ad[i]=Btl.list_size();
      }
      else {cout << "Error in term handling!" << endl;exit(0);}
      mark(++ct,10000,"lex_pair counts");
   }
   map=new long[mapsize];
   k=0;
   for(i=1;i<nwrd1+2;i++){
      j=ad[i];
      ad[i]=k;
      k+=j;
   }
   if(k!=mapsize){cout << "Error in summing!" << endl;exit(0);}

   //Write out "ad" file (offsets into map);
   bin_Writ("ad",(nwrd1+2)*sizeof(long),(char*)ad);
   
   //Set numbers in map.
   ct=0;
   Btl.node_first();
   while(Btl.node_next()){
      if(i=Coord1.find(Btl.show_str())){
         ad1=ad[i];
         k=0;
         Btl.set_ptr();
         while((cptr = Btl.next_ptr()) != NULL) {
            j=Coord2.find(cptr);
            map[ad1+(k++)]=j;
         }
         sSort(k,map+ad1);
      }
      mark(++ct,10000,"lex_entries");
   }
   //Write out map.
   bin_Writ("map",sizeof(long)*mapsize,(char*)map);
   ofstream *pfout=get_Ostr("nm");
   *pfout << nwrd1 << " " << nwrd2 << " " << mapsize << endl;
   dst_Ostr(pfout);
   if(pflag)cout << nwrd1 << " " << nwrd2 << " " << mapsize << endl;

   return(1);
}

void RelateA::gopen_bmatrix(void){
   char cnam[max_str],*cptr;
   long ct,mapsize;

   ifstream *pfin=get_Istr("nm");
   *pfin >> nwrd1 >> nwrd2 >> mapsize;
   if(pflag)cout << nwrd1 << " " << nwrd2 << " " << mapsize << endl;
   dst_Istr(pfin);

   //map Hash structures
   strcpy(cnam,name);
   strcat(cnam,"-s1");
   Coord1.change_name(cnam);
   Coord1.gopen_htable_map();
   strcpy(cnam,name);
   strcat(cnam,"-s2");
   Coord2.change_name(cnam);
   Coord2.gopen_htable_map();

   //map binary files
   ad=(long*)get_Mmap("ad");
   map=(long*)get_Mmap("map");
}

void RelateA::set_wrd(long n){
   ad1=ad[n];
   ad2=ad[n+1]-1;
}

int RelateA::exs_tag(long m){
   long x=ad1;
   long y=ad2;
   long i,j;
   switch(y-x){
      case 0: if(map[x]!=m)return(0);
              else return(1);
              break;
      case 1: if((map[x]!=m)&&(map[y]!=m))return(0);
              else return(1);
              break;
      default: if((j=map[x])>=m){
                  if(j==m)return(1);
                  else return(0);
               }
               else if((j=map[y])<=m){
                  if(j==m)return(1);
                  else return(0);
               }
               else {
                  while(y-x>1){
                     i=(x+y)/2;
                     if((j=map[i])==m)return(1);
                     else if(j<m)x=i;
                     else y=i;
                  }
               }
               return(0);
   }
}
                     
int RelateA::exs_pair(long n,long m){
   set_wrd(n);
   return(exs_tag(m));
}

int RelateA::exs_pair(const char *str1,const char *str2){
   long n,m;
   if(n=Coord1.find(str1)){
      set_wrd(n);
      if(m=Coord2.find(str2)){
         return(exs_tag(m));
      }
      else return(0);
   }
   else return(0);
}

//RelateB code

RelateB::RelateB(void) : FBase("relateb","null"){
   name=NULL;
}

RelateB::RelateB(const char *nam) : FBase("relateb",nam){
}

RelateB::~RelateB(){
}

int RelateB::create_bmatrix(BTList &Btl,int exc1,int exc2){
   char cnam[max_str],*cptr;
   long ct;

   nwrd1=Btl.cnt_key;
   nwrd2=Btl.lst->cnt_key;
   cout << nwrd1 << " " << nwrd2 << endl;

   //Create Hash for strings set 1
   strcpy(cnam,name);
   strcat(cnam,"-s1");
   Coord1.change_name(cnam);
   Coord1.create_htable(Btl,exc1);
   //Create Hash for strings set 2
   strcpy(cnam,name);
   strcat(cnam,"-s2");
   Coord2.change_name(cnam);
   Coord2.create_htable(*(Btl.lst),exc2);

   //Write out number of words and number of tags
   ofstream *pfout=get_Ostr("nm");
   *pfout << nwrd1 << " " << nwrd2 << endl;
   dst_Ostr(pfout);

   //Load into Blist structures
   Coord1.gopen_htable_map();
   Coord2.gopen_htable_map();

   //Make bitmap
   long i,j,k;
   rd=(nwrd2+1)/8+1;
   map=new unsigned char[(nwrd1+1)*rd];
   for(i=0;i<(nwrd1+1)*rd;i++){
      map[i]=0;
      mark(i,100000,"zeros");
   }
   ct=0;
   Btl.node_first();
   while(Btl.node_next()){
      if(i=Coord1.find(Btl.show_str())){
         wrd=map+i*rd;
         Btl.set_ptr();
         while((cptr = Btl.next_ptr()) != NULL) {
            j=Coord2.find(cptr);
            (*(wrd+j/8))+=1<<(j%8);
         }
      }
      else {cout << "Error in term handling!" << endl;exit(0);}
      mark(++ct,10000,"lex_entries");
   }
   //Write out bit map.
   bin_Writ("bit",(nwrd1+1)*rd,(char*)map);

   return(1);
}

void RelateB::gopen_bmatrix(void){
   char cnam[max_str],*cptr;
   long ct;

   ifstream *pfin=get_Istr("nm");
   *pfin >> nwrd1 >> nwrd2;
   if(pflag)cout << nwrd1 << " " << nwrd2 << endl;
   dst_Istr(pfin);

   //map Hash structures
   strcpy(cnam,name);
   strcat(cnam,"-s1");
   Coord1.change_name(cnam);
   Coord1.gopen_htable_map();
   strcpy(cnam,name);
   strcat(cnam,"-s2");
   Coord2.change_name(cnam);
   Coord2.gopen_htable_map();

   //map binary file
   rd=(nwrd2+1)/8+1;
   map=(unsigned char*)get_Mmap("bit");
}

void RelateB::set_wrd(long n){
   wrd=map+n*rd;
}

int RelateB::exs_tag(long m){
   return(*(wrd+m/8)&(1<<(m%8)));
}

int RelateB::exs_pair(long n,long m){
   wrd=map+n*rd;
   return(*(wrd+m/8)&(1<<(m%8)));
}

int RelateB::exs_pair(const char *str1,const char *str2){
   long n,m;
   if(n=Coord1.find(str1)){
      set_wrd(n);
      if(m=Coord2.find(str2)){
         return(exs_tag(m));
      }
      else return(0);
   }
   else return(0);
}

//Lexicon code

Lexicon::Lexicon(const char *nam) : FBase("lexset",nam){
}

Lexicon::~Lexicon(){
}

int Lexicon::create_bmatrix(const char *ppth,int exc1,int exc2){
   char cnam[max_str],*cptr,*uptr;
   long ct;
   BTList Btl;

   ifstream fin(ppth,ios::in);
   if(!fin.is_open()){cout << cnam << " failed to open!" << endl;exit(0);}
   ct=0;
   while(fin.getline(cnam, max_str, '\n')) {
      cptr = strtok(cnam, " ");
      if(*cptr){
         while((uptr = strtok(NULL, " ")) != NULL) {
            if(*uptr)Btl.add_unique(cptr,uptr);
         }
      }
      mark(++ct,10000,"lex_entries");
   }
   fin.close();

   nwrds=Btl.cnt_key;
   ntags=Btl.lst->cnt_key;
   cout << nwrds << " " << ntags << endl;

   //Create Hash for strings set 1
   strcpy(cnam,name);
   strcat(cnam,"-s1");
   Bterms.change_name(cnam);
   Bterms.create_htable(Btl,exc1);
   //Create Hash for strings set 2
   strcpy(cnam,name);
   strcat(cnam,"-s2");
   Btags.change_name(cnam);
   Btags.create_htable(*(Btl.lst),exc2);

   //Write out number of words and number of tags
   ofstream *pfout=get_Ostr("nm");
   *pfout << nwrds << " " << ntags << endl;
   dst_Ostr(pfout);

   //Map Hash structures
   Bterms.gopen_htable_map();
   Btags.gopen_htable_map();

   //Make bitmap 
   long i,j,k;
   rd=(ntags+1)/8+1;
   map=new unsigned char[(nwrds+1)*rd];
   ptag=new long[nwrds+1];
   for(i=0;i<(nwrds+1)*rd;i++){
      map[i]=0; 
      mark(i,100000,"zeros");
   }
   ct=0;
   Btl.node_first();
   while(Btl.node_next()){
      if(i=Bterms.find(Btl.show_str())){
         wrd=map+i*rd;
         Btl.set_ptr();
         cptr = Btl.next_ptr();
         j=Btags.find(cptr);
         ptag[i]=j;
         (*(wrd+j/8))+=1<<(j%8);
         while((cptr = Btl.next_ptr()) != NULL) {
            j=Btags.find(cptr);
            (*(wrd+j/8))+=1<<(j%8);
         }
      }
      else {cout << "Error in term handling!" << endl;exit(0);}
      mark(++ct,10000,"lex_entries");
   }
 
   bin_Writ("bit",(nwrds+1)*rd,(char*)map);
   bin_Writ("pt",(nwrds+1)*sizeof(long),(char*)ptag);
   
   return(1);
}

void Lexicon::gopen_bmatrix(void){
   char cnam[max_str],*cptr;
   long ct;
   
   ifstream *pfin=get_Istr("nm");
   *pfin >> nwrds >> ntags;
   dst_Istr(pfin);
   if(pflag)cout << nwrds << " " << ntags << endl;

   //map Hash structures
   strcpy(cnam,name);
   strcat(cnam,"-s1");
   Bterms.change_name(cnam);
   Bterms.gopen_htable_map();
   strcpy(cnam,name);
   strcat(cnam,"-s2");
   Btags.change_name(cnam);
   Btags.gopen_htable_map();

   //map binary files
   rd=(ntags+1)/8+1;
   map=(unsigned char*)get_Mmap("bit");
   ptag=(long*)get_Mmap("pt");
}

void Lexicon::set_wrd(long n){
   wrd=map+n*rd;
}
   
int Lexicon::exs_tag(long m){
   return(*(wrd+m/8)&(1<<(m%8)));
}

int Lexicon::exs_pair(long n,long m){
   wrd=map+n*rd; 
   return(*(wrd+m/8)&(1<<(m%8)));
}

}

