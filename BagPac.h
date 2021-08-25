#ifndef BAGPAC_H
#define BAGPAC_H
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <runn.h>
#include <DStor.h>

using namespace std;
namespace iret {

template<class Z>
class BagPac {
   public:
      BagPac(void); 
     ~BagPac(void);
      void add_str(Z m,const char *str); 
      void add_str(Z m,string &str); 
      void pack(DSpan &Ds); //Packs the contents into Ds
      int unpack(DSpan &Ds); //Unpacks the contents of Ds
      void clear(void); //Clears the vectors
      void show_data(void); //Debug function
         //Prints data to stdout
      long id; //Record id
      long date; //Record date
      vector<Z>   si;
      vector<string> ss;
};

template<class Z>
BagPac<Z>::BagPac(void){
} 

template<class Z>
BagPac<Z>::~BagPac(){
}  
 
template<class Z>
void BagPac<Z>::add_str(Z m,const char *str){
   string ptt(str);
   ss.push_back(ptt);
   si.push_back(m);
}

template<class Z>
void BagPac<Z>::add_str(Z m,string &str){
   ss.push_back(str);
   si.push_back(m);
}

template<class Z>
void BagPac<Z>::pack(DSpan &Ds){
   int i=0,j,k=sizeof(Z);
   Ds.id=id;
   Ds.date=date;
   vector<string>::iterator p=ss.begin();
   while(p!=ss.end()){
      i+=p->size()+1;
      p++;
   }
   i+=si.size()*k+sizeof(int);
   if(Ds.bln<i){
      if(Ds.buf)delete [] Ds.buf;
      Ds.buf=new char[i];
      Ds.bln=i;
   }
   Ds.len=i;
   j=si.size();
   char *pc=(char*)&j;
   Ds.buf[0]=pc[0];
   Ds.buf[1]=pc[1];
   Ds.buf[2]=pc[2];
   Ds.buf[3]=pc[3];
   j=sizeof(int);
   typename vector<Z>::iterator q=si.begin();
   while(q!=si.end()){
      pc=(char*)&(*q);
      for(i=0;i<k;i++)Ds.buf[j++]=pc[i];
      q++;
   }
   p=ss.begin();
   string::iterator u;
   while(p!=ss.end()){
      u=p->begin();
      while(u!=p->end()){
         Ds.buf[j++]=*u;
         u++;
      }
      Ds.buf[j++]='\0';
      p++;
   }
}

template<class Z>
int BagPac<Z>::unpack(DSpan &Ds){
   int i=0,j,k=sizeof(Z),*pi;
   Z *pz;
   id=Ds.id;
   date=Ds.date;
   pi=(int*)Ds.buf;
   j=pi[0];
   pz=(Z*)(Ds.buf+sizeof(int));
   for(i=0;i<j;i++){
      si.push_back(pz[i]);
   }
   string stx;
   i=sizeof(Z)*j+sizeof(int);
   while(i<Ds.len){
      stx=Ds.buf+i;
      i+=stx.size()+1;
      ss.push_back(stx);
   }
   return(ss.size());
}

template<class Z>
void BagPac<Z>::clear(void){
   si.clear();
   ss.clear();
}

template<class Z>
void BagPac<Z>::show_data(void){
   cout << "id " << id << endl;
   cout << "date " << date << endl;
   vector<string>::iterator p=ss.begin();
   typename vector<Z>::iterator q=si.begin();
   while(p!=ss.end()){
      cout << *q << "| " << *p << endl;
      p++;
      q++;
   }
}

}
#endif
