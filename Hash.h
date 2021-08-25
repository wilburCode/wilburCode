#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <fstream>
#include <Btree.h>
#include <Thes.h>
#include <FBase.h>
#include <runn.h>
#include <Dmap.h>

namespace iret {

class Hash : public FBase {
public:
  Hash(void);
  Hash(const char *nm);
  Hash(int n,const char *nm); //n gets appended to type if >-1
  ~Hash();

  void SetMem(Hash &Hsh); //sets data and memory to Hsh data & memory
  void create_htable(List &Lst,int excess); //"str" for file of strings, 
      //"ad" for address file, "nm" numbers, 
      //"ha" hash array. Excess is # powers of 2 above size.
  void create_htableM(List &Lst,int excess); //creates in memory ready for use
      //and no need to call gopen or gclose functions
  void create_htable(int mz,List &Lst,int excess); //"str" for file of strings, 
      //Creates a numbered version of above

  void create_htable(strMap &Mp,int excess); 
     //Calls create_htable.
  void create_htable(int mz,strMap &Mp,int excess); 
     //Creates a numbered version of above

  void gopen_htable_map(void); //Creates memory maps
  void gopen_htable_map(int mz); //Creates memory maps
  void gclose_htable_map(void); //Destroys memory maps
     //and deletes memory
  void gclose_htable_map(int mz); //Destroys memory maps
     //and deletes memory
  void gopen_htable_copy(Hash *pH); //Copies memory maps

  long find(const char *str); //Return number+1 if present, else 0.
      //Number is not lexical order but hash order and then lexical
      //within collesion groups.

  //Data
  char *strmap; //Holds the bit map.
  long *addr; //Holds the offsets to strmap.
  long nwrds; //Number of words.
  long tnum; //Truncation number, size of har.
  long *harr; //Holds hash array.
  long *farr; //Holds the hash coefficients.
  long *px0;
  long *px1;
  long *px2;
  long *px3;
  long *px4;
  long *px5;
  long *px6;
  long *px7;
  long *px8;
  long *px9;
  long *px10;
  long *px11;
  bool own_hash_mem; //True if owns the memory maps, else false
  int num_file; //Holds the number used in gopen_map, else -1
};

class Chash : public Hash {
public:
  Chash(void);
  Chash(const char *nm);
  Chash(int n,const char *nm); //n gets appended to type if >-1
  ~Chash(void);

  void SetMem(Chash &Chsh); //sets data and memory to Chsh data & memory
  void create_ctable(Count &Ct,int excess); //Adds "ct" for counts
     //Calls create_htable and then prodoces the array of counts.
  void create_ctable(int mz,Count &Ct,int excess); //Adds "ct" for counts
     //Creates a numbered version of above
  void create_ctable(List &Lt,int excess); //Adds "ct" for term # 
     //and starts the count at 1 and in lexical order. count() will
     //return 0 if term not in list.
  void create_ctable(int mz,List &Lt,int excess); //Adds "ct" for term # 
     //Creates a numbered version of above
 
  void create_ctable(strMap &Mp,int excess); 
     //Adds "ct" for counts
     //Calls create_htable and then prodoces the array of counts.
  void create_ctable_STerm(strMap &Mp,int excess); 
     //Adds "ct" for counts. Counts are index of term + 1
     //Calls create_htable and then prodoces the array of counts.
  void create_ctable(int mz,strMap &Mp,int excess); 
     //Adds "ct" for counts
     //Creates a numbered version of above

  void gopen_ctable_map(void); //Calls gopen_htable_map and also
     //maps "ct" file.
  void gopen_ctable_map(int mz); //Calls gopen_htable_map and also
     //maps "ct" file.
  void gclose_ctable_map(void); //Calls gclose_htable_map and also
     //Unmaps "ct" file.
  void gclose_ctable_map(int mz); //Calls gclose_htable_map and also
     //Unmaps "ct" file.

  long count(const char *str); //Returns count if present, else 0.

  //Data
  long *cnt;
  bool own_chash_mem; //True if owns the memory maps, else false
};

class RelateA : public FBase {
public:
  RelateA(void);
  RelateA(const char *nm);
  ~RelateA();

  int create_bmatrix(BTList &Btl,int exc1,int exc2);
      //BTList must have unique lists (add_unique).
      //exc1 and exc2 are to be between 1 & 3, higher the faster
      //as long as size is not a problem.

  void gopen_bmatrix(void); //maps the data from hash
      //and "arr" files.
  void set_wrd(long n); //n is word number and sets wrd
      //pointer
  int exs_tag(long m); //m is tag number. 1 if present,
      //0 if absent.
  int exs_pair(long n,long m); //m is tag number. 1 if present,
      //0 if absent.
  int exs_pair(const char *str1,const char *str2); //1 if present, else 0.

  //Data
  long *map; //Holds the map array.
  long *ad;  //Holds offsets into map.
  long nw; //Current string number.
  long ad1; //Lower address in map.
  long ad2; //Upper address in map.
  Hash Coord1; //Holds the first coordinate strings.
  Hash Coord2; //Holds the second coordinate strings.
  long nwrd1; //Number of first coordinate strings.
  long nwrd2; //Number of second coordinate strings.
};

class RelateB : public FBase {
public:
  RelateB(void);
  RelateB(const char *nm);
  ~RelateB();

  int create_bmatrix(BTList &Btl,int exc1,int exc2);
      //BTList must have unique lists (add_unique).
      //exc1 and exc2 are to be between 1 & 3, higher the faster
      //as long as size is not a problem.

  void gopen_bmatrix(void); //maps the data from hash
      //and "bit" files.
  void set_wrd(long n); //n is word number and sets wrd
      //pointer
  int exs_tag(long m); //m is tag number. 1 if present,
      //0 if absent.
  int exs_pair(long n,long m); //m is tag number. 1 if present,
      //0 if absent.
  int exs_pair(const char *str1,const char *str2); //1 if present, else 0.

  //Data
  unsigned char *map; //Holds the bit map.
  unsigned char *wrd;  //Points at map for a given word.
  Hash Coord1; //Holds the first coordinate strings.
  Hash Coord2; //Holds the second coordinate strings.
  long nwrd1; //Number of first coordinate strings.
  long nwrd2; //Number of second coordinate strings.
  long rd; //Reduction of second coordinate to bits needs rd char.
};

class Lexicon : public FBase {
public:
  Lexicon(const char *nm);
  ~Lexicon();

  int create_bmatrix(const char *path,int exc1,int exc2);
      //path is the full path and name to lexicon. Data
      //must consist of one word per line and categories
      //on the same line following word and demarcated by spaces.
      //exc1 and exc2 are to be between 1 & 3, higher the faster
      //as long as size is not a problem.

  void gopen_bmatrix(void); //Loads in the data from "s1",
      //"s2" and "bit" files.
  void set_wrd(long n); //n is word number and sets wrd
      //pointer
  int exs_tag(long m); //m is tag number. 1 if present,
      //0 if absent.
  int exs_pair(long n,long m); //m is tag number. 1 if present,
      //0 if absent.

  //Data
  unsigned char *map; //Holds the bit map.
  unsigned char *wrd;  //Points at map for a given word.
  Hash Bterms; //Holds the words.
  Hash Btags; //Holds the tags.
  long nwrds; //Number of words.
  long ntags; //Number of tags.
  long rd; //Reduction of ntags to bits needs rd char.
  long *ptag; //Array of nwrds most probable tag numbers.
};

}
#endif
