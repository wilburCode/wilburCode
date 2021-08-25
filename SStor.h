#ifndef SSTOR_H
#define SSTOR_H
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

class Sset : public FBase{
public:
   Sset(const char *nm,const char *path_nm); //nm is name of object
      //path_nam is pointer at a string sss and reads the path from file path_sss in
      //current directory. But if sss begins with ':' then skips this character and 
      //remaining string is the path string itself.
   ~Sset();
   void gopen_write(void); 
   void add_str(const char *pch); //Call to add a string to the set
   void gclose_write(void);
     //Files are .n, .a, and .s

   void gopen_map(void);
   char *show_str(long n); //Returns pointer at the nth string
   void gclose_map(void);
  
   //Data
   long num; //number of strings in the set (file .n)
   long *addr; //array of offsets to strings (file .a)
   char *str;  //pointer at string map (file .s)
   ofstream *pfa; //Used for file creation
   ofstream *pfs; //Used for file creation
};

class SStor : public Sset {
   public:
      SStor(const char *nam,const char *path_nam); //nam is name of object
         //path_nam is pointer at a string sss and reads the path from file path_sss in
         //current directory. But if sss begins with ':' then skips this character and 
         //remaining string is the path string itself.
     ~SStor();
      void create_SStor(set<string> &stx);
      void update_SStor(set<string> &stx);
      void reOrder(void); //Puts "s" set into lexical order
         //Updates will get the strings out of order generally, but not
         //necessary to correct this except occasionally for efficiency's sake

      long find(const char *ssr); //If finds string
      //returns its number+1 else returns 0.
      //Does binary search.
      long lfind(const char *ssr); //Finds by binary search that i, 0<=i<num,
         //where i is the index of string closest to str from below. If match is
         //perfect returns i+1, else returns -(i+1).
      void gopen_append(void); //Opens string object to append in an update operation
      int stc_my(const char *,const char *);
      //Function used to compare two strings.

      int a; //Variable used for comparison of strings.
      int b; //Variable used for comparison of strings.
};

class TStor : public Sset {
   public:
      TStor(const char *nam,int nd,const char *path_nam); //Name of object, nd is dimension of data
         //path_nam is pointer at a string sss and reads the path from file path_sss in
         //current directory. But if sss begins with ':' then skips this character and 
         //remaining string is the path string itself.
         //Usual files for stringset and "nd", "da", "dd" for dimension and data files
     ~TStor();
      void create_TStor(map<string,long *> &mpx);
      void update_TStor(map<string,long *> &mpx);
      void reOrder(void); //Puts "s" set into lexical order
         //Updates will get the strings and data out of order generally, but not
         //necessary to correct this except occasionally for efficiency's sake
      void gopen_TStor(void); //Maps files for reading
      int read(long m,pLong &ar); //Sets data pointer ar to point at data for mth string
         //Returns 0 if m out of range, else 1
      int read(const char *ptr,pLong &ar); //Sets data pointer ar to point at data for string
         //ptr if found. Returns 1 if string pointed at by ptr is found, else 0
      void gclose_TStor(void); //destroys memory maps

      long find(const char *ssr); //If finds string
      //returns its number+1 else returns 0.
      //Does binary search.
      long lfind(const char *ssr); //Finds by binary search that i, 0<=i<num,
         //where i is the index of string closest to str from below. If match is
         //perfect returns i+1, else returns -(i+1).
      void gopen_append(void); //Opens string object to append in an update operation
      int stc_my(const char *,const char *);
      //Function used to compare two strings.

      int a; //Variable used for comparison of strings.
      int b; //Variable used for comparison of strings.
      int ndim; //Dimension of data for each string
      long *da; //Pointer at start of data for each string (an index offset for each string)
      long *dd; //Pointer at data 
};

}
#endif
