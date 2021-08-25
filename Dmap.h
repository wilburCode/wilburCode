#ifndef DMAP_H
#define DMAP_H
#include <fstream>
#include <iostream>
#include <runn.h>
using namespace std;
namespace iret {

template<class Z>
class Dmap : public std::map<const char *,Z,SCmp<const char *> > {
   public:
      Dmap(void); 
     ~Dmap(void);
   void add_count(const char *str,Z n); //Adds the string *str with its count
      //to the tree if not already in list. String is key and count is data.
      //If string is already a key the count is incremented by n.
      //keeps total count, but also counts number of unique keys in count.
   void addp_count(const char *str,Z n); //Adds the string *str with its count
      //just as add_count.
      //Does not make copy of string, but uses the pointer str as key pointer.
   void correct(const char *str,Z n); //If str is in the tree the count is
      //changed to n. Otherwise nothing is done.

   //Functions for maximum calculation
   void max_count(const char *str,Z n); //Adds the string *str with its count
      //to the tree if not already in list. String is key and count is data.
      //If string is already a key the count is max of n and prior value.
      //keeps total count, but also counts number of unique keys in count.
   void maxp_count(const char *str,Z n); //Adds the string *str with its count
      //just as max_count.
      //Does not make copy of string, but uses the pointer str as key pointer.

   //Functions for minium calculation
   void min_count(const char *str,Z n); //Adds the string *str with its count
      //to the tree if not already in list. String is key and count is data.
      //If already a key the count is min of n and prior value.
      //keeps total count, but also counts number of unique keys in count.
   void minp_count(const char *str,Z n); //Adds the string *str with its count
      //just as min_count.
      //Does not make copy of string, but uses the pointer str as key pointer.

   void Merge(Dmap<Z> *pNp); //Adds all content from pNp into map and destroys
      //pNp

   inline void SMclear(void){
      const char *pch;
      Set();
      while(qs!=qz){
         pch=qs->first;
         qn=qs;
         qs++;
         this->erase(qn);
         delete [] pch;
      }
   }

   Z count(const char *str); //Returns the count if a key (in list) otherwise
      //returns 0.
   inline void Set(void){ //Sets iterators ready to iterate through the map
      qs=this->begin();
      qz=this->end();
   }
   void debug(); //Prints to stdout to help debugging
   Z total; //Holds the total of all counts added for all keys.
   long cnt_key; //Holds number of unique keys
   typename map<const char *,Z>::iterator qs,qz,qn;
};

template<class Z>
Dmap<Z>::Dmap() : map<const char *,Z,SCmp<const char *> >() {
   total=0;
   cnt_key=0;
}

template<class Z>
Dmap<Z>::~Dmap() {
  SMclear();
}

template<class Z>
void Dmap<Z>::Merge(Dmap<Z> *pNp){
   const char *pch;
   Z n;
   typename map<const char *,Z>::iterator q;

   pNp->Set();
   while(pNp->qs!=pNp->qz){
      pch=pNp->qs->first;
      n=pNp->qs->second;
      pNp->qn=pNp->qs;
      pNp->qs++;
      pNp->erase(pNp->qn);
      if((q=this->find(pch))!=this->end()){
         q->second+=n;
         delete [] pch;
      }
      else {
         this->insert(make_pair(pch,n));
         cnt_key++;
      }
      total+=n;
   }
   pNp->~Dmap();
}

template<class Z>
void Dmap<Z>::add_count(const char *str,Z n){
   int lxn;
   char *pch;
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      q->second+=n;
   }
   else {
      lxn=strlen(str);
      pch=new char[lxn+1];
      strcpy(pch,str);
      this->insert(make_pair(pch,n));
      cnt_key++;
   }
   total+=n;
}

template<class Z>
void Dmap<Z>::addp_count(const char *str,Z n){
   int lxn;
   char *pch;
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      q->second+=n;
   }
   else {
      this->insert(make_pair(str,n));
      cnt_key++;
   }
   total+=n;
}

template<class Z>
void Dmap<Z>::correct(const char *str,Z n){
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      q->second=n;
   }
}

template<class Z>
void Dmap<Z>::max_count(const char *str,Z n){
   int lxn;
   Z y;
   char *pch;
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      y=q->second;
      if(y<n){
         q->second=n;
         total+=n-y;
      }
   }
   else {
      lxn=strlen(str);
      pch=new char[lxn+1];
      strcpy(pch,str);
      this->insert(make_pair(pch,n));
      cnt_key++;
      total+=n;
   }
}

template<class Z>
void Dmap<Z>::maxp_count(const char *str,Z n){
   int lxn;
   Z y;
   char *pch;
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      y=q->second;
      if(y<n){
         q->second=n;
         total+=n-y;
      }
   }
   else {
      this->insert(make_pair(str,n));
      cnt_key++;
      total+=n;
   }
}

template<class Z>
void Dmap<Z>::min_count(const char *str,Z n){
   int lxn;
   Z y;
   char *pch;
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      y=q->second;
      if(y>n){
         q->second=n;
         total+=n-y;
      }
   }
   else {
      lxn=strlen(str);
      pch=new char[lxn+1];
      strcpy(pch,str);
      this->insert(make_pair(pch,n));
      cnt_key++;
      total+=n;
   }
}

template<class Z>
void Dmap<Z>::minp_count(const char *str,Z n){
   int lxn;
   Z y;
   char *pch;
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      y=q->second;
      if(y>n){
         q->second=n;
         total+=n-y;
      }
   }
   else {
      this->insert(make_pair(str,n));
      cnt_key++;
      total+=n;
   }
}

template<class Z>
Z Dmap<Z>::count(const char *str){
   typename map<const char *,Z>::iterator q,qe=this->end();
   if((q=this->find(str))!=qe){
      return(q->second);
   }
   else return(0);
}

template<class Z>
void Dmap<Z>::debug(void){
   Set();
   while(qs!=qz){
      cout << qs->second << " " << qs->first << endl;
      qs++;
   }
}

typedef Dmap<long> strMap;

//Class to add up values associated with keys
template<class Y,class Z>
class ADmap : public std::map<Y,Z>{
   public:
   ADmap();
   ~ADmap();
   
   void add(Y xm,Z xn);  
   void max(Y xm,Z xn);  
   void min(Y xm,Z xn);  
};

template<class Y,class Z>
ADmap<Y,Z>::ADmap() : map<Y,Z>() {
}

template<class Y,class Z>
ADmap<Y,Z>::~ADmap() {
}

template<class Y,class Z>
void ADmap<Y,Z>::add(Y xm,Z xn){
   typename map<Y,Z>::iterator q,qe=this->end();

   if((q=this->find(xm))!=qe){
      q->second+=xn;
   }
   else {
      this->insert(make_pair(xm,xn));
   }
}

template<class Y,class Z>
void ADmap<Y,Z>::max(Y xm,Z xn){
   typename map<Y,Z>::iterator q,qe=this->end();

   if((q=this->find(xm))!=qe){
      if(q->second<xn)q->second=xn;
   }
   else {
      this->insert(make_pair(xm,xn));
   }
}

template<class Y,class Z>
void ADmap<Y,Z>::min(Y xm,Z xn){
   typename map<Y,Z>::iterator q,qe=this->end();

   if((q=this->find(xm))!=qe){
      if(q->second>xn)q->second=xn;
   }
   else {
      this->insert(make_pair(xm,xn));
   }
}

}
#endif
