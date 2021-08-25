#ifndef NSTOR_H
#define NSTOR_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <FBase.h>
#include <sys/mman.h>
#include <set>
#include <string>
#include <map>

using namespace std;
namespace iret {

typedef long* pLong;

class NStor : FBase {
   public:
      NStor(const char *nam,const char *path_nam); //Name of object
         //path_nam is pointer at a string sss and reads the path from file path_sss in 
         //current directory. But if sss begins with ':' then skips this character and 
         //remaining string is the path string itself.
     ~NStor();
      void create_NStor(set<long> &stx);
      void update_NStor(set<long> &stx);
         //files are "n" and "x"

      void gopen_map(void); //Opens map of data
      long find(long m); //If finds number
      //returns its index+1 else returns 0.
      //Does binary search.
      long lfind(long m); //Finds by binary search that i, 0<=i<num,
         //where i is the index of number closest to m from below. If match is
         //perfect returns i+1, else returns -(i+2).
      void gclose_map(void); //Closes map of data

      long num; //Size of array
      long *xn; //Pointer at mapped array
};

}
#endif
