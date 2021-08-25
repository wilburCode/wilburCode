#ifndef MAP_H
#define MAP_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <Btree.h>
#include <FBase.h>
#include <EdtPrb.h>

using namespace std;
namespace iret {

class Level {
   public:
      Level(void);
      ~Level(void);
      //Data
      long dp; //Depth
      long nm; //Length
      long tn; //Current term number
      long mt; //Current length of seg
};

class State {
   public:
      State(long n); //Creates a state with memory for a depth n search
      ~State(void);
      void view(void); //For debugging
      void Save(Level &Lv); //Stores current state level
      void Reset(Level &Lv); //Assumes Lv was from a save followed by
        //an extension and resets the state to Lv level.
      void Copy(State &Sk); //Duplicates Sk in current state.

      //Data
      long dp; //Depth of search, initially 0
      long nm; //Length of match held in the state
      long *tn; //Holds term numbers
      long *mt; //Holds seg length at each term number
      char *ch; //Holds the matched portion of string
};

class Map : public FBase {
   public:
      Map(void); 
      Map(const char *nm); 
      ~Map(void);
      void create_Map(List &Ls); //Creates the map structures in several files
        //"c" file for string fragments, "a" for addresses, "s" for space in
        //p file, "p" for place to jump if missmatch.
        //Note that must use add_key_count so that total size of string set
        //is available to the function.
      void create_Map(long n,List &Ls); //Creates the map structures in several files
        //but uses n as a number for the files and this number appears as part of
        //the name.
      void gopen_map(void);
      void gopen_map(long n);
        //Opens files with n as part of the name.
      void gopen_Map_copy(Map *pM); //Copys the addresses of needed files
        //from pM.
      long find(const char* str); //Attempts to find the string str. If
        //successful returns its position+1 where pos. starts at 0.
      long bfind(const char* str); //Attempts to find the string str. If
        //successful returns its position+1 where pos. starts at 0.
        //Like find but should not have an initial | in string.
      long sfind(const char* str,State &Sk); //Attempts to find the 
        //string str. If successful returns its position+1 where pos. 
        //starts at 0. If unsuccessful Sk has the results of the search
      long sbfind(const char* str,State &Sk); //Attempts to find the
        //string str. If successful returns its position+1 where pos.
        //starts at 0. If unsuccessful Sk has the results of the search
        //Like sfind but should not have an | at beginning of string.
      long extend(const char* str,long m,State &Sk); //Attempts to extend Sk 
        //using str starting at position m. If successful returns 
        //its position+1 where pos. starts at 0. If unsuccessful Sk 
      long extend_s(const char* str,long m,State &Sk); //Same as previous 
        //except extends the Sk.ch to full string to end of str if extention
        //fails. Used on local_two_Edit.
      long step_Back(State &St); //Moves back one character (decreases nm
        //by 1. If successful returns 1, else 0.
      long step_Forward(State &St); //Advances nm by 1, starts at top poss.
        //If successful returns 1, else 0.
      long step_Down(State &St); //Moves down by one if possible,
        //If successful returns 1, else 0.
      void local_one_Edit(const char* str,State &Sk,Count &Ct); //Assumes str 
        //is not found. All one edits that occur at level where search halts 
        //are returned in Ct.
      void global_one_Edit(long lev,const char* str,State &Sk,Count &Ct); 
        //Assumes str is not found. All one edits beginning at lev-1 or higher
        //in set are returned in Ct. Of course no one edit is above the state
        //Sk.
      void local_two_Edit(const char* str,State &Sk,Count &Ct); 
      void global_two_Edit(long lev,const char* str,State &Sk,Count &Ct); 
        //Two edits function much like one edits, but only find those
        //two edits whose first edit cannot be extended to a match in the
        //set. 

      //Data
      char *cf; //Points at the "c" file
      long *ad; //Points at the "a" file
      long *st; //Points at the "s" file
      long *pt; //Points at the "p" file
      //for spell checking
      char *ptr; //Work space 
      char *qtr; //Work space 
      char *utr; //Work space 
      long lxn; //Length of string in ptr
};

class CMap : public Map {
   public:
      CMap(void);
      CMap(const char *nm,EdtPrb *pEt);
      ~CMap(void);
      void create_CMap(Count &Ct); //Like Map, but adds an array for counts "f"
      void create_CMap(long n,Count &Ct); //Like Map, but adds an array for counts "f"
        //Uses n as part of file names.
      void gopen_cmap(void);
      void gopen_cmap(long n);
        //Uses n as part of file names.
      void gopen_CMap_copy(CMap *pCM); //Copys the needed files from pCM
      long cfind(const char* str,State &Sk,long lm); //To follow sfind and have
        //the Sk produced by sfind as its argument Sk. It will check if
        //Sk ends on a word and if so will attempt to introduce a space and
        //find the remainder as another term in the database 
        //If successful returns a 1, else a 0. Puts the result in Sk.ch
        //lm is lower limit for both halves frequencies.
      long cfind2(const char* str,State &Sk,long lm); //To follow sfind 
        //Like cfind except puts only first part of split in Sk.ch
        //Called by Search_spg
      long check_sngl(const char *str,long l1,long l2); //Checks the frequency
        //of the individual words of length l1 and l2 and returns the minimum.

      void local_one_XEdit(const char* str,State &Sk); 
      void global_one_XEdit(long lev,const char* str,State &Sk);
      void local_two_XEdit(const char* str,State &Sk);
      void global_two_XEdit(long lev,const char* str,State &Sk);
        //Like the Edit functions but give value to the edits based the 
        //strings involved.
      void global_ph1_XEdit(long lev,const char* str,State &Sk);
      void local_ph2_XEdit(const char* str,State &Sk);
      void global_ph2_XEdit(long lev,const char* str,State &Sk);
        //Like the Edit functions but give value to the edits based the 
        //strings involved. Special for phrases with short words
        //These functions use the limits e1w,e1h,e2w,e2h

        //SEdit is to take a sum of scores and compare with fmax
      void local_one_SEdit(const char* str,State &Sk); 
      void global_one_SEdit(long lev,const char* str,State &Sk);
      void global_ph1_SEdit(long lev,const char* str,State &Sk);
        //Like the preceeding but uses limits based on word length
      void set_lim_ph1(const char *str); //Resets the e2w and e1h limits
        //based on current state of string. Needed after an edit has
        //been done

      long Search_one(const char *str); //Function to search for an edit 
        //solution at depth<=1
      long Search_ph1(const char *str); //Function to search for an edit 
        //solution at depth<=1
        //Uses limits based on word length
      long Search_two(const char *str); //Function to search for an edit 
        //solution at depth=2
      long Search_ph2(const char *str); //Function to search for an edit 
        //solution at depth=2, but limits based on iph array. Invoked in
        //response to flag irp, 0 no special processing, 1 for it.
      long Search_cmb(const char *str); //Function that combines one and
        //two searches to look for a solution. Also has split search.
      long Search_cmg(const char *str); //Function that combines one and
        //two searches to look for a solution. Also has split search.
        //Differs from cmb in putting only first part of split in cmax
      long Search_cm1(const char *str); //Function that combines one and
        //two searches to look for a solution. No split search but does
        //use ph2 search as needed
      long Search_cm2(const char *str); //Function that combines one and
        //two searches to look for a solution. No split search

      long Search_dep(const char *str); //Function to search deep for a
        //solution, uses Search_cm2
      long Rat_check(const char *st1,const char *st2); //Used in Search_dep
        //to check the rationality of the solution.
      long RBeg(const char *st1,const char *st2); //Used to check the
        //first three characters of the strings and returns 1 if two or more
        //of these are the same.

      long Search_qst(const char* str, State &Sk); //Function to extend
        //looking for a solution and split off that solution.
      long cmstr(long i,const char *cmax,const char *str,long &j);
        //Function to find the match point between cmax and str.

      long Search_spl(const char *str,long lm); //Function to split a string if
        //two parts are found in database.
        //lm is lower limit for both halves frequencies.
      long Search_spg(const char *str,long lm); //Function to split a string if
        //two parts are found in database.
        //lm is lower limit for both halves frequencies.
        //Differs from spl in only placing first part of split in cmax
        //Calls cfind2 and called by cmg.

      void local_one_Exten(const char* str,State &Sk); 
      void global_one_Exten(long lev,const char* str,State &Sk);
      void local_two_Exten(const char* str,State &Sk);
      void global_two_Exten(long lev,const char* str,State &Sk);
        //These functions act only when there is no two edit solution
        //They examine the set of two edit matches that are maximal in length
        //and the winner is the highest frequency term
      long local_one_Qxten(const char* str,State &Sk);
      long global_one_Qxten(long lev,const char* str,State &Sk);
      long local_two_Qxten(const char* str,State &Sk);
      long global_two_Qxten(long lev,const char* str,State &Sk);
        //These functions act and return a 1 or 2 edit solution if found
        //If no solution found they examine the set of two edit matches
        //that are maximal in length and the winner is the highest prob term
        //Also keep track of best State and where in query string it reached.


      void Extract_sng(const char *str); //Does like Search_tot, but str one word
      void Extract_par(const char *str,long l1,long l2); //Does a complete search
        //For a solution for a two word str, l1 and l2 the lengths of the terms
      long Extract_lng(const char *str,long l1,long l2); //Does a search for a
        //solution, but only for a single term solution, no splits or separate
        //terms. To use in truncated space
      long Extract_qst(const char *pre,const char *str); //Attempts to extend pre
        //into str to get a longer solution
      long Extract_fnl(const char *str,long l1,long l2); //Attempts to find pair 
        //of words str in database by edit and by split. Returns 1 for edit, 
        //2 for split. l1 and l2 the lengths of the terms.

      //Data
      long *freq; //To map "f" file.
      char cmax[max_str]; //Storage for winning string
      char rmax[max_str]; //Additional Storage for string
      long ixt; //marks where extension stopped in original string
      double sumx; //For the sum of fmax.
      double fmax; //Frequency of winning string
      double def1; //Place to store deficit values for one edit
      double def2; //Place to store deficit values for 1st of 2 edits.
      long lnr2; //Length of match between str and correction in Exten
      long lmax; //Length of match between str and correction in Exten
      State *pSq; //Pointer at best state.
      long en; //Carrier for edit, +1 del, -1 ins
      long em; //Carrier for edit, +1 del, -1 ins, best current
      long ek; //Carrier for 2nd edit, +1 del, -1 ins, best current

      long irp; //1 for ph1 processing, 0 otherwise
                //2 for ph2 & ph1 processing for two edits
      long ln1; //Length of first word in 2 word phrase
      long ln2; //Length of second word in 2 word phrase
      long e1w; //Lower limit of first edit
      long e1h; //Upper limit of first edit
      long e2w; //Lower limit of second edit
      long e2h; //Upper limit of second edit
      double eps; //Used to regulate downweighting for low freq
      double fac[81]; //Factor to decrease a probability of usefulness by
                //depending on the term freq up to 80.

      EdtPrb *pEp; //Used to estimate probs of edits, results stored in def2
      long n1; //Things gotten from EdtPrb
      long n2;
      long n3;
      int *cn;
      double *div_del;
      double *div_ins;
      double *div_rep;
      double *div_trn;
};

class Manager {
   public:
      Manager(void);
      Manager(CMap *pMqS,CMap *pMqT,CMap *pMqP);
      ~Manager(void);

      void gopen_map(void); //Opens all three classes for use.
      void gopen_Man_copy(Manager *pMan); //Copies all three
         //classes from pMan.
      void Reset(void); //Resets imax to 0. cmax to empty string
      void Extract(const char *str); //Used recursively to build cmax
      long Rat_check(const char *st1,const char *st2);
         //Rationality check based on the kth word in st1 and first
         //in st2.

      //Data
      CMap *pMS;
      CMap *pMT;
      CMap *pMP;

      long imax; //marks current end of string in cmax
      char cmax[max_str];
};

}
#endif
