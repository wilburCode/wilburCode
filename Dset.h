#ifndef DSET_H
#define DSET_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <set>

using namespace std;
namespace iret {

class Dset : public std::set<const char *,SCmp<const char *> > {
   public:
      Dset(void);
     ~Dset(void);
   void add_key(const char *str); //Adds the string *str
      //to the tree if not already in list. 
      //Number of unique keys in count is kept.
      //Need to call clear to free the memory holding strings 
      //before Dset is deleted to prevent a memory leak
   void addp_key(const char *str); //Adds the pointer to string str
      //Assumes the string is in memory that is safe and will not need
      //to be free'd by this class
      //Number of unique keys in count is kept.
   int ifind(const char *str); //Returns 1 if a key (in list) otherwise
      //returns 0. Use this when Set has not been called since last addition
      //to the set.
   int zfind(const char *str); //Returns 1 if a key (in list) otherwise
      //returns 0. Use this when Set has been called since last addition
   inline void Set(void){ //Sets iterators ready to iterate through the set
      rs=this->begin();
      rz=this->end();
   }
   inline void SMclear(void){
      Set();
      while(rs!=rz){
         delete [] *rs;
         rs++;
      }
      this->clear();
   }
   long cnt_key; //Holds number of unique keys
   typename std::set<const char *,iret::SCmp<const char *> >::iterator rs,rz;
};

} 
#endif

