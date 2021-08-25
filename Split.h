#ifndef Split_H
#define Split_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <FBase.h>
#include <DataObj.h>
#include <vector>
#include <map>


using namespace std;
namespace iret {

class Split{
   public:
      Split(void); 
      Split(int wrd_spc); 
     ~Split();

      void token(const char *str); //Produces a list of tokens in order of
         //of occurrence in the string.
      void token_lower(const char *str); //Produces a list of tokens in order of
         //of occurrence in the string. Lowercased.
      void tokenS_lower(const char *str); //Produces a list of tokens in order of
         //of occurrence in the string. Lowercased.
         //Removes punctuation and hyphens.
      void tokenS(const char *str); //Produces a list of tokens in order of
         //of occurrence in the string. 
         //Removes punctuation and hyphens.
      void token_tab(const char *str); //Produces a list of tokens in order of
         //of occurrence in the string. Breaks at tabs
      void token_tab_lower(const char *str); //Produces a list of tokens in order of
         //of occurrence in the string. Lowercased.
      void token_chr(char ch,const char *str); //Produces a list of tokens in order of
         //of occurrence in the string. Breaks at ch.
      void tokenS_chr(char ch,const char *str); //Produces a list of tokens in order of
         //of occurrence in the string. Breaks at ch.
         //Removes punctuation and hyphens.
      void token_chr_lower(char ch, const char *str); //Produces a list of tokens in 
         //order of occurrence in the string. Lowercased.
      void tokenS_chr_lower(char ch, const char *str); //Produces a list of tokens in 
         //order of occurrence in the string. Lowercased.
      void puncChunk(const char *str); //Creates chunks in order with breaks at 
         //punctuation marks as done in the S version of functions. Thus get
         //chunks that essentially take the place of words.
      void loadTrecn(void); //Loads tokens with standard labels for use in
         //tokenQuery function.
      long nNum(char *str); //if string is digits either null terminated or punct
         //terminated it computes the value and returns. If first digit is 0 or no digits
         //returns 0.
      void tokenQuery(const char *str,vector<string> &vx); //Creates chunks in order 
         //with breaks at spaces and ';' and ':' and '"' marks. Also breaks 
         //at '(' and matching ')' at the top level. Returns a vector of strings
         //in 1-1 correspondence with tokens and gives best guess of what the token
         //represents.
      void tokenKN(const char *str); //Splits the string into tokens which
         //alphanumeric. Lower cases all but stop words. For use in processing
         //known item queries and the docs against which they are to retrieve
      void couple(char *buf,char *s1,char *s2); //inserts a space between
         //s1 and s2 and puts result in buf. Assumes that no white space in s1
         //or s2.
      void lower_lst(void); //Lower cases the lst strings produced by 
         //token
      void clear(void); //Clear the lst memory of words
      void clean(const char *pch); //Takes a string as generally produced in
         //Textools with ! and * characters and removes extraneous characters
         //and makes lst of tokens
      char* assembl(void); //Assembles all the strings in lst in order with a
         //space between each and into a single normalized string. Returns a
         //pointer at the string space produced by new
      char* assembl_tab(void); //Assembles all the strings in lst in order with a
         //tab between each and into a single normalized string. Returns a
         //pointer at the string space produced by new
      char* pair(long k); //Assembles the tokens at k and k+1 in lst to make
         //a two word phrase and returns a pointer at the result
      char* segment(long k,long m); //Assembles tokens at k, k+1, ..., k+m-1 to make
         //an m token phrase with spaces between tokens and returns a pointer at result
      void order(void); //Reorders lst into lexical order
      char* leftMod(long k,string st); //adds the string st to the left end of the
         //kth string
      char* rightMod(long k,string st); //adds the string st to the right end of the
         //kth string
      
      //Data
      long word_space; //Space in lst for tokens
         //default 10k
      long num; //Number of tokens
      char **lst; //Holds the tokens
      char *xnam;
      char *ynam; //string space
      map<string,string> Trecn;
};
}
#endif
