#ifndef HYPER_H
#define HYPER_H

#include <iostream>
#include <fstream>
using namespace std;
namespace iret {

class Hyper {
public:
  Hyper(void); //Creates the shell of the object
  Hyper(long n);
     //The number n is the largest set size that will be used 
     //in calls to functions.
  Hyper(long n,double eps);
     //Like above except allows the approximate function 
     //to be used also.
  ~Hyper();
  double nlog_pval(long n_st,long n_s, long n_t, long N);
     //This returns the negative log10 of the pval that the 
     //overlap will be n_st or greater for the given subsets.
  double nlog_pval_appx(long n_st,long n_s, long n_t, long N);
     //This returns the negative log10 of the pval just
     //as previous function except pval allowed to be less by
     //factor of (1-eps). 
  double nlog_pval_appx2(long n_st,long n_s, long n_t, long N);
  double log_prob(long n_st,long n_s, long n_t, long N);
     //Returns the log10 of the probability that the overlap
     //will be exactly n_st for the two sets.
  double addl(double x,double y);
     //Returns the log of the sum of the two numbers whose
     //logs are input. All logs base 10.
  double log_binomCoeff(long m,long n);
     //Returns the log of binomial coefficient, base 10
  double logOdds(long n_st,long n_s,long n_t,long N);
     //Returns the log of odds ratio that t and s are related
     //Generally pos if n_st larger than expected, neg otherwise
     //log is base 10
  double HlogOdds(long n_st,long n_s,long n_t,long N);
     //Like above, but adjust total to obtain
     //the numerator from HyperG distribution

  long nobj; //total set size;
  double *log_fac; //array for log(factorials).
  double *log_num; //array for log(numbers).
  double eps; //Error factor for speed.
  double epz; //equals epx
};

class SHyper {
public:
  SHyper(void);
     //The number n is the largest set size that will be used 
     //in calls to functions.
  SHyper(double eps);
     //Like above except allows the approximate function 
     //to be used also.
  ~SHyper(void);
  double nlog_pval(long n_st,long n_s, long n_t, long N);
     //This returns the negative log10 of the pval that the 
     //overlap will be n_st or greater for the given subsets.
  double nlog_pval_appx(long n_st,long n_s, long n_t, long N);
     //This returns the negative log10 of the pval just
     //as previous function except pval allowed to be less by
     //factor of (1-eps). 
  void set_log_prob(long n_st,long n_s, long n_t, long N);
     //Sets xi and i1-i4. xi is prob overlap is n_st
  double log_prob(long k);
     //Returns the log10 of the probability that the overlap
     //will be exactly k for the two sets.
  double log_fac(long l,long u); //Returns log10 of u!
  double addl(double x,double y);
     //Returns the log of the sum of the two numbers whose
     //logs are input. All logs base 10.

  double eps; //Error factor for speed.
  double xi;  //Prob of a given overlap
  long i1;
  long i2;
  long i3;
  long i4;
};

}
#endif
