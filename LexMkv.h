#ifndef LEXMKV_H
#define LEXMKV_H
//In RAM version of variable order Markov string predictor
//Disc based version is in Strset.h & C files
#include <iostream>
#include <fstream>
#include <runn.h>
#include <Btree.h>

using namespace std;
namespace iret {

class Lex {
public:
   Lex(void);
   ~Lex(void);
   void create_Lex(List &Ls);
   char *show_str(long n); //Returns address of nth
      //string.
   long find(const char *str); //If finds string
      //returns its number+1 else returns 0.
      //Does binary search.
   int stc_my(const char *,const char *);
      //Function used to compare two strings.
   long lfind(const char *str); //returns number
      //+1 for longest string that matches an
      //initial segment of the string str.
      //Otherwise returns 0.
   long tfind(const char *str); //returns number of
      //string matches found which begin at beginning
      //of str. Held in ht. If ht[i]=n+1>0 then a match
      //of length i occurs at string n [0,num-1]. Calls ifind.
   long ifind(long,long,const char*); //Called by lfind
      //and tfind
   int stc_ly(const char *,long); //Called in ifind
      //Functions for matches to initial segment
      //list members
   long find_low(const char *str); //Finds first string in lexos
      //that str matches initial segment, returns index+1, or 0
   int stc_low(const char *,const char *);
   long find_high(const char *str); //Finds last string in lexos
      //that str matches initial segment, returns index+1, or 0
   int stc_high(const char *,const char *);
  
   //Data
   long num; //number of strings in the set (file .n)
   long *addr; //array of offsets to strings (file .a)
   char *str;  //pointer at string map (file .s)
   long space1; //Space for str
   long space2; //Space foro addr
   int slen; //Length of string str in tfind call.
   int a; //Variable used for comparison of strings.
   int b; //Variable used for comparison of strings.
   char *tx; //Holds string copy for lfind and tfind.
   long *ht; //Holds the matchs at various depths.
      //At i puts the index+1 of a match of length i.
      //If no match of length i, holds 0.
   long *sn; //Holds the beginning of pair
   long *sm; //Holds the end of pair
   int pflag; //Print flag, 1 gives output, 0 no output
};

class LexMkv : public Lex {
public:
   LexMkv(int mn); //mn is minimum length segment to add.
   ~LexMkv(void);
   void set_max(int mx); //Sets max.
   void set_aug(double ax); //Sets aug.
   void set_fst(char fx,long buf_siz=10000); //Sets fst & buf size

   //Add to DCount, only useful if min is set = 1
   void add_Print(DCount &Dc); //Adds all the printable characters to Dc with
      //counts of 1.0
   void add_Char(const char* str,double xx,DCount &Dc); //Adds the single characters
      //from str with counts of xx to Dc
   
   //Suffix tree approach
   void create_LexMkv(DCount Ct); //Makes the basic string file and a file
      //of counts (.c).
   void create_LexMkv(char fx,DCount Ct); //As above but adds fx as first
      //character in each string.
   double count_Super(const char *stt); //Finds the count of all stings that
      //match this string in their initial part.

   //Probability calculating functions
   double prob_ext(const char *str); //Returns the prob of str
      //Computes using as long a piece as possible at each step.
   double prob_ext_fst(const char *str); //Returns the prob of str
      //Computes using as long a piece as possible at each step.
      //Uses the special first character from set_fst
   double prob_fxd(const char *str); //Returns the prob of str
      //Computes using pieces as long as possible but limited by max.
   double prob_fxd_fst(const char *str); //Returns the prob of str
      //Computes using pieces as long as possible but limited by max.
      //Uses the special first character from set_fst

   //Data
   int min; //Lower limit on length of segments to add.
            //Used in creation.
   int max; //Maximum lookback in computing prob_fxd.
            //Uses back off to obtain probs as needed.
   double aug; //Number to augment counts for first segment of a string
            //Used with add_segs_aug function
            //Default is zero.

   double tot; //Total of all counts in set
   double *cnt; //Holds the counts
private:
   //Used to handle extra character at beginning of string
   char *buf;
   char *bfc;
};

}
#endif

