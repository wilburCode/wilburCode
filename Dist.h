#ifndef DIST_H
#define DIST_H

using namespace std;
namespace iret {

double lphi(double);

class LBin{
public:
   LBin(long n); //Allocate n double space 
   ~LBin();
   double combination(long n, long r); 
   //Compute log of the binomial cofficient
   double log_binomial(long n, long r, double p);
   //logarithm of one of term in the binomial expansion
   double log_binomial_pval(long n, long r, double p); 
   //p-value of the above 
   double addl(double x,double y);
   //Given x=log(X) and y=log(Y), compute log(X+Y)

 private:
   double *log_fac; //array for log(factorials).
};

class Binomial{
public:
   Binomial(double errs, double errp,long maxN);
   //errs for searching error eps 
   //errp for p-value error epp
   //maxN is largest integer that will be  used
   ~Binomial(void);
   double fact(long m); //Returns log(m!);
   double combnr(long n, long r); 
   //computes log binomial coefficient 
   double binomial_pval_appx(long n, long r, double p); 
   //fast way compute the p-value
   double upper_limit(long N, long r, double sig);
   //N-number of trials
   //r-number of sucesses
   //sig-significance level
   //errp should be << sig for accuracy in computation
   //errs is error in the upper_limit computed
   //same holds for lower limit
   double lower_limit(long N, long r, double sig);
 private:
  double eps;//error searching
  double epp;//error in p-value 
  double *xfac; //Stores part of log of factorial
    //Uses jumps of 1000
};

//Wilcoxin-Mann_Whitney Test
//k values from one source are mixed with n-k values from a second source
//to make a set of n values. The values are ranked by size and s is the sum
//of the ranks of the k values in the resulting ordering. Their are three
//functions that estimate a pvalue, i.e., the probability that the sum of 
//ranks would be as small as s for a random sample of k values from the 
//ordering. The function wmwcomp(s,k,n) computes an exact value by a
//recursion formula. Its time complexity is of the order of s*k*n or s*s*k
//(the smaller of the two). It is to lengthy when k and s and n are all large.
//Thus wmwappx is the normal approximation based on Larson. It is relatively
//good when k>=10 and n-k>=10. The function wmwpval uses wmwcomp when either
//k or n-k is below 10 and wmwappx otherwise.

class Wmw {
   public:
      Wmw(void);
     ~Wmw(void);

      double wmwpval(long s,long k,long n); //s is rank sum to
         //test, k is number of terms, n the actual space size
         //Combines wmwcomp and napprox when the latter can be
         //used (k>=10 and n-k>=10).
      double wmwcomp(long s,long k,long n); //s is rank sum to
         //test, k is number of terms, n the actual space size
      void setup(long s,long k,long n);  //Sets the reduction
         //and memory for the computation
      void calc(long j); //Computes one step of the algorithm
         //j corresponds to the space size that decreases by one at
         //each step or call
      double wmwappx(long s,long k,long n); //Good when k and
         //n-k both as great as 10 (Larson).
      long sn; //rank sum bound
      long kn; //number of elements in rank sum
      long zn; //Size of space;
      long ct; //Reduced size of space
      double ux; //Numerical value at reduced size
      double wx; //Cumulated value of probability
      double **sx; //origin array
      double **tx; //destination array
      long ka; //Array size for kn
      long sa; //Array size for sn
      double sr2; //Square root of 2
};

}
#endif
