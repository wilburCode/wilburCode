#ifndef IDSET_H
#define IDSET_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <DataObj.h>
#include <DStor.h>

using namespace std;
namespace iret {

template<class Y>
class Idset : public FBase {
   public:
      Idset(void); //for copying
      Idset(const char *nam,const char *path_nam); //path_nam is pointer at a string
         //sss and reads the path from file path_sss in current directory. But if sss
         //begins with ':' then skips this character and remaining string is the path
         //string itself. 
     ~Idset(void);
      long create_Idset(DStor &Ds); //Returns number in set
      void gopen_Idset(void); //Memory maps data
      //Accessing ids
      void SetMem(Idset &Id); //Sets data and memory to Id data and memory
      long idxf(long m); //Returns mth id if 0<=m<idnum
      Index *idset(Indx<Y> *pInd); //pInd points at index object with
         //indices and largest must be <idnum or returns NULL.
      Index *totid(void); //Returns Index object including all ids, 
         //memory mapped
      //Accessing indices of ids
      long find(long id); //Returns index+1 when id found, else 0
      Indx<Y> *indexset(Index *pId); //Returns a pointer at an index
         //object that lists the indices of *pId that are in the set
         //NULL if the intersection is empty
      Indx<Y> *totindex(void); //Returns list of all indices for set.
      void gclose_Idset(void); //Unmaps the id list
      //DATA
      long idnum;
      long *idx;
};

template<class Y>
Idset<Y>::Idset(void) : FBase("idset","null"){
} 

template<class Y>
Idset<Y>::Idset(const char *nam,const char *path_nam) : FBase("idset",nam){
   if(*path_nam!=':'){
      set_path_name(path_nam);
   }
   else {
      set_path_internal(path_nam+1);
   }
} 

template<class Y>
Idset<Y>::~Idset(){
}  
 
template<class Y>
long Idset<Y>::create_Idset(DStor &Ds){
   long i,j;

   ofstream *pfout=this->get_Ostr("ids",ios::out);
   i=0;
   Ds.gopen_read();
   while(j=Ds.getNext()){
      j--;
      pfout->write((char*)&j,sizeof(long));
      mark(i++,10000,"records");
   }
   this->put_Nnum("num",i);
   this->dst_Ostr(pfout);
   return(i);
}

template<class Y>
void Idset<Y>::gopen_Idset(void){
   this->get_Nnum("num",idnum);
   idx=(long*)this->get_Mmap("ids");
}

template<class Y>
long Idset<Y>::idxf(long m){
   if(m<0)return(-1);
   if(m>=idnum)return(-1);
   return(idx[m]);
}

template<class Y>
Index *Idset<Y>::idset(Indx<Y> *pInd){
   long i;
   if(pInd->idx[pInd->ix-1]>=idnum){
      cout << "Error, index too large!" << endl;return(NULL);
   }
   if(pInd==NULL)return(NULL);
   if(pInd->ix==0)return(NULL);
   Index *pId=new Index(pInd->ix);
   for(i=0;i<pId->ix;i++){
      pId->idx[i]=idx[pInd->idx[i]];
   }
   return(pId);
}

template<class Y>
Index *Idset<Y>::totid(void){
   Index *pId=new Index;
   pId->ix=idnum;
   pId->idx=idx;
   return(pId);
}

template<class Y>
long Idset<Y>::find(long id){
   long i,j,k,x,y;

   if(id<idx[0]){cout << id << " Below range!" << endl;return(0);}
   if(id>idx[idnum-1]){cout << id << " Above range!" << endl;return(0);}
   x=0;
   y=idnum;
   while(y-x>1){
      i=(y+x)/2;
      if(id>=idx[i])x=i;
      else y=i;
   }
   if(id>idx[x])return(0);
   else return(x+1);
}
   
template<class Y>
Indx<Y> *Idset<Y>::indexset(Index *pId){
   long j,k;
   Y i;

   Index *pInd=totid();
   Index *pJnd=pInd->cbool_And(pId);
   if(pJnd)cout << pId->ix-pJnd->ix << " ids lost" << endl;
   else cout << "all lost" << endl;
   if(pJnd==NULL)return(NULL);
   Index *pKnd=pInd->Subvalue(pJnd);
   Indx<Y> *pQnd=new Indx<Y>(pKnd->ix);
   for(i=0;i<pKnd->ix;i++){
      pQnd->idx[i]=pKnd->idx[(long)i];
   }   
   delete pJnd;
   delete pKnd;
   return(pQnd);
}

template<class Y>
Indx<Y> *Idset<Y>::totindex(void){
   return(new Indx<Y>(0,idnum));
}

template<class Y>
void Idset<Y>::gclose_Idset(void){
   this->dst_Mmap("ids",(char*&)idx);
}

template<class Y>
void Idset<Y>::SetMem(Idset &Id){
   idnum=Id.idnum;
   idx  =Id.idx;
}

}
#endif
