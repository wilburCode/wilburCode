#ifndef EDTPRB_H
#define EDTPRB_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <Btree.h>
#include <FBase.h>

using namespace std;
namespace iret {

//All strings good and bad are assumed to begin with the character '|'
class EdtPrb : public FBase { 
   public:
      EdtPrb(const char *nm);
      ~EdtPrb(void);
      void gopen_map(void);

      void set_mem(void); //Sets memory for counts
      long ct_edits(const char *good,const char *bad);
         //Return -1 if edits too complex
         //Return number of edits otherwise.
      void ct_locat(const char *good,double fg); //Counts the baseline context
         //fg is frequency of the good string in the data
      void ct_edit2(const char *good,const char *bad,double fb); //If ct_edits is good
         //this function used to count the actual edits in context. fb is frequency of
         //bad string in the data


      void convert_prb(void);//Converts the div_ array counts to probs.
      void zero_diagonals(void);//Zeroes the div_rep and div_trn arrays 
         //on the diagonals where the transformation is an identity.
      void write_prb(void);//writes the div_ arrays to files.
         //Files named in usual way with extensions "d", "i", "r", and "t".
      void debug_prb(const char *stt); //Give a word stt and it prints probs for various
         //edits
      //Profiles functions
      void zero_profile(void); //Zeroes out arrays before counting.
      void edit_profile(const char *good,const char *bad,double fb); //If ct_edits is good
         //this funciton counts the edits by position and prints results to screen.
      void see_profile(void); //Writes the edit profiles to screen.

      //Simulation of edits in a string
      char *sim_edit(const char *str); //Introduces a single edit into str and returns a
         //pointer at the altered string in new space.
      char *sim_2edit(const char *str); //Introduces 2 edits into str
      char *sim_3edit(const char *str); //Introduces 3 edits into str

      //Functions to estimate the probability of converting str to btr by mistakes
      double Probc(long lsn,const char *str,long lbn,const char *btr);
      double Probx(long lsn,char *str,long lbn,char *btr); 
      double Probs(long lsn,const char *str,long lbn,const char *btr); 
      long Space(long lsn,const char *str);

      //Data
      long n1;
      long n2;
      long n3;
      long n4;
      int con[128]; //Conversion of characters to ints
      double *bas_del;
      double *div_del; //Counts for deletions
      double *bas_ins;
      double *div_ins;
      double *bas_rep;
      double *div_rep;
      double *bas_trn;
      double *div_trn;
      double *del;
      double *ins;
      double *rep;
      double *trn;
      //String variables
      char *gdd;
      char *bdd;
};

}

#endif
