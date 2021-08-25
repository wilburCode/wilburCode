#ifndef CLIP_H
#define CLIP_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <Btree.h>
#include <Word.h>
#include <Hash.h>

using namespace std;
namespace iret {

class Clip : public Count {
   public:
      Clip(void); //nm=0
      Clip(long nm); //nm size of longest string to handle
      ~Clip(void);

      void set1(void); //Setup for proc0-3
      void prc0(const char *str); //Produces the strings having alpha or
         //digit and demarcated by other. Other is not included.
      void prc1(const char *str); //Produces the strings demarcated
         //by changes between alpha, digit, and other.
         //The other is always left out of the strings produced.
      void prc2(const char *str,Hash &Hs); //Produces the strings demarcated
         //by changes between alpha, digit, and other and found in Hs.
      void prc2(const char *str,Chash &Cs,long llim); //Produces the strings demarcated
         //by changes between alpha, digit, and other. Only adds strings found
         //in Cs with count less than llim.
      void prc3(const char *str,Hash &Hs); //Removes all chars that are not
         //alpha or digit. From resulting single string adds all contiguous
         //substrings found in Hs.
      void prc3(const char *str,Chash &Cs,long llim); //Removes all chars that are not
         //alpha or digit. From resulting single string adds all contiguous
         //substrings found in Cs and with count less then  llim.

      void spac(const char *str); //Breaks at spaces and takes pieces
      void spac(const char *str,long n); //Breaks at spaces and takes pieces
         //Uses n as the count

      void substr(const char *tok,long ng); //Adds in all substrings with from
         //1 to ng contiguous chars. All get count of 1 on add
      void substr_b(const char *tok,long ng); //Adds in all substrings of the form
         //_abc or abc_ or _abc_ where _ marks token boundary and 
         //number of characters from the token is 1-ng at most. 
      void lsubstr(const char *tok,long ng); //Adds in all substrings with from
         //1 to ng contiguous chars. All get count of 1 on add. Lower cased
      void lsubstr_b(const char *tok,long ng); //Adds in all substrings of the form
         //_abc or abc_ or _abc_ where _ marks token boundary and 
         //number of characters from the token is 1-ng at most. Lower cased.
      void token(const char *tok); //Adds the token into the list
      void ltoken(const char *tok); //Adds the token into the list
         //lower cased.
      void bigram(const char *tk1,const char *tk2); //Adds the tokens as
         //a bigram to the list
      void lbigram(const char *tk1,const char *tk2); //Adds the tokens as
         //a bigram to the list. Lower cased. 
      
      void lincode(double xd,double acc); //acc is a small positive number
         //representing accurracy and xd is the number to be encoded.
      void lincode(double xd,double acc,const char *prf); //acc is a small 
         //positive number representing accurracy and xd is the number 
         //to be encoded. Appends the string prf to the beginning of 
         //the strings that are made. Codes xd as integer multiple of 1/acc.
      void bincode(double xd,long acc); //acc is a small positive integer
         //representing accurracy and xd is the number to be encoded.
      void bincode(double xd,long acc,const char *prf); //acc is a small 
         //positive integer representing accurracy and xd is the number 
         //to be encoded. Appends the string prf to the beginning of the 
         //strings that are made. Does a binary encoding
      void logcode(long n,long amp); //Uses amp & log2 to code to a small integer
      void logcode(long n,long amp,const char *prf); //Uses amp & log2 to code 
         //to a small integer
      void mstring(const char *prf,long m); //converts m to string of 
         //digits and appends prf to the beginning of the string.
      void mstring(const char *prf,long m1,long m2); //converts m1 and 
         //m2 to strings and separates them by '/' and appends prf to 
         //the beginning.
      void mstring(const char *prf,long m1,long m2,long m3); //As previous 
         //but with three numbers

      void addin(Count &Ct); //adds a copy of Ct into current set.
      void addin(Count &Ct,List &Lt); //adds a copy of those members of Ct
         //that are also in Lt into the current set. 
         //All the functions from here down too prod add in strings
         //with a count of 1
      void prefix(const char *prf,Count &Ct); //Takes strings from Ct and appends
         //prf to beginning and adds them and their count to current set.
      void prefix(const char *prf,Count &Ct,List &Lt); //Takes strings from Ct 
         //if they are in Lt and appends prf to beginning and adds them and 
         //their count to current set.
      void suffix(const char *suf,Count &Ct); //Takes strings from Ct and appends
         //suf to end and adds them and their count to current set.
      void suffix(const char *suf,Count &Ct,List &Lt); //Takes strings from Ct 
         //if they are in Lt and appends suf to end and adds them and 
         //their count to current set.
      void ptoken(const char *prf,const char *pch); //appends pch to end of prf
         //without any space inserted and adds to current set.
      void stoken(const char *suf,const char *pch); //appends suf to end of pch
         //without any space inserted and adds to current set.

      void prod(Count &C1,Count &C2); //Adds all pairs of strings st1*st2 to current
         //set.
      void prod(Count &C1,Count &C2,List &Lt); //Adds all pairs of strings 
         //st1*st2 to current set if and only if they are in Lt.

      void fill_grm(const char *txt,int ngram,int aug,Word &Wrd);
           //Fills a Count with ngrams made from the words in the array txt
           //of length ln. Links consecutive words by ngrams of the form "a_bcd"
           //where the a is the first letter of a word and the bcd is the initial
           //part of the next word and total length is ngram. Augments a marked
           //form the first trigram of each word  by entering it aug times and
           //also marks the first letter of each word and adds it in once.
           //Uses the object Wrd of class Word to produce words from which grams
           //made.
      void fill_grmo(const char *txt,int ngram,int aug,Word &Wrd);
           //Like fill_grm, but older version that does not add an initial 2-gram
      void fill_brm(const char *txt,int lgram,Word &Wrd);
           //Fills a Count with 2-lgram ngrams. Takes the ngrams as
           //composed from consecutive letters in the whole string. All
           //ngrams are added in only once.
      void fill_scn(const char *txt,int lgram,Word &Wrd);
           //Fills a Count with lgrams. Takes the lgrams as
           //composed from consecutive chars in the words of the string. All
           //lgrams are added in only once.

      void Set_Up_Word(Word &Wrd); //Sets up word class for standard
           //TexTool processing
      void proc_MrkSng(long len,const char *txt,Word &Wrd,const char *suf); //Adds all
           //single tokens in the text with the suffix suf added to the end
      void proc_title(long len,const char *txt,Word &Wrd);
      void proc_body(long len,const char *txt,Word &Wrd);
           //These functions do the standard document construction
           //as is done in TexTool. 
      void proc_tisng(long len,const char *txt,Word &Wrd); //Like title
           //but only !!t
      void proc_bdsng(long len,const char *txt,Word &Wrd); //Like body 
           //but only !!t
      void proc_Vbarstring(const char *txt); //Processes a string with
           //entities separated by vertical bars. Assumes no spaces at
           //ends of entities. Adds the individual strings
           //Assumes txt is not over 10000 bytes long and ascii
          
           //Applies stemming to words at extraction
      void proc_title_stem(long len,const char *txt,Word &Wrd);
      void proc_body_stem(long len,const char *txt,Word &Wrd);
           //These functions do the standard document construction
           //as is done in TexTool.
      void proc_tisng_stem(long len,const char *txt,Word &Wrd); //Like title
           //but only !!t
      void proc_bdsng_stem(long len,const char *txt,Word &Wrd); //Like body
           //but only !!t

      //Breaks up strings into tokens
      int split(const char *str);
      int split_lower(const char *str); //Lower cases also

   //Data
      long num; //Size of longest string to handle
      char zt[128];
      char *pck;
      long *nck;
      char xmam[10000]; //Work space
      //For splits of strings
      long nx; //Number of tokens
      char **lst; //Holds the tokens
      long word_space; //Space in lst for tokens
         //default 10k
};
}
#endif
