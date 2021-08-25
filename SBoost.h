#ifndef SBOOST_H
#define SBOOST_H

#include <iostream>
#include <fstream>
#include <DataObj.h>
#include <Isgrid.h>

using namespace std;
namespace iret {

class SBoost  : public FBase {
public:
  SBoost(long dm, long ln, Index *gind,const char *nam);
  //dim is dimention of the data, i.e. how many score arrays are to be 
  //combined.
  //len is the length of score arrays
  //gInd is the pointer to the good elements

  SBoost(long dm,const char *nam);
  //dim is dimention of the data, i.e. how many score arrays are to be 
  //combined.
  
  ~SBoost();
 
  void Setup(long i, double *sco);
  //Sets up scr through dim calls
 
  void Learn_Boost(void);
  void Save_Boost(void);
  void Load_Boost(void);
 
  double Score(double *scc);
  //Predicts the class

  //Data
  long dim;
  long len;//total number of points
  long n;//number of goods
  Index *gInd;
  double **scr;
  Isgrid **pIs;
  long grn; //Granularity, default 10,000
  double zprod; //Product of Zs
  double epsilon; //Limit on probabilities
  long t; //Number of the iteration.
  double epl;//Stopping Criteria for Optimization
  long pflag; //On by default.
};

class SABoost : public FBase {
public:
  SABoost(long dm, long ln, Index *gind, const char *nam);
  //dim is dimention of the data, i.e. how many score arrays are to be 
  //combined.
  //len is the length of score arrays
  //gInd is the pointer to the good elements
  //constructor to be used in train programs

  SABoost(long dm, const char *nam);
  //COnstructor to be used in test programs
  
  ~SABoost();
 
  void Setup(long i, double *sco);
  //Sets up scr through dim calls
 
  void Learn_Boost(void);
  void Save_Boost(void);
  void Load_Boost(void);
 
  double Score(double *scc);
  //Predicts the class

  double Z_alpha(double alp, long dm); //computes the value for alp.

  //Data
  long dim;
  long len;//total number of points
  long n;//number of goods
  Index *gInd;
  double **scr;
  double **alpha_h;
  double *alpha;
  double *label;
  double *sum; //Used for the Z_alpha function.
  double zprod; //Product of Zs
  double epsilon; //Limit on probabilities
  long t; //Number of the iteration.
  double epl;//Stopping Criteria for Optimization
  long pflag; //On by default.
};
}
#endif
