#ifndef PERM_H
#define PERM_H
#include <fstream>
#include <iostream>
#include <runn.h>
using namespace std;
namespace iret {

class Perm {
   public:
      Perm(long nm); //nm is array length
     ~Perm();
      long next(void); //Returns a 1 when a new permutation has been
          //found. Returns 0 when there are no more.
          //After each call to this function that returns a 1 there is
          //a different permutation in mm.
      long num; //array length
      long *mm; //array of indices. Holds the permutations generated
      long *bd; //bounding array
      long *cd; //counting array
      long *rm; //remapping array (work space)
};

class Comb {
   public:
      Comb(long nm); //nm is array length
     ~Comb();
      void set_mm(long t); //Sets the initial configuration of the array mm to
          //have t nonzeroes (left).
          //Here t must be >0.
      long next(void); //Finds next configuration with tt nonzeroes and
          //returns 1 if successful and 0 if unsuccessful.
          //Set up so one is to call first and then mm is available if returns 1
      long num; //array length
      long tt;  //Subset size <=num.
      long ct;  //State counter
      long *mm; //Array indicating to include by 1, not include by 0.
};

class PermPairs {
   public:
      PermPairs(long nm); //nm is array length of perm
     ~PermPairs();
      long Copy(long *tx);    //tm=tx.
      long ConvertF(long *m); //m is perm of length num
         //returns number of 1's.
         //result of this function is representation of m as a tm array of
         //0's and 1's.
      void ConvertB(long *m); //from tm produces corresponding m
      long Consist(void); //Checks tm for consistency and returns 1 if and 0
         //if not.
      long num; //array length of perm
      long tnm; //array length of perm pairs num(num-1)/2
      long *tm; //array of pairs representation, 1 for <, 0 for >.
};

}
#endif
