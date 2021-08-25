#ifndef GOLDOPT_H
#define GOLDOPT_H
#include <iostream>
#include <fstream>
namespace iret {
class GoldOpt {
   public:
      GoldOpt(void);
      ~GoldOpt(void);
      void set_init_bracket(double a, double c);
      double give_value();
      void return_value(double xx);
      int next(void);//It will set four points
      void set_limit(double xx);
      double eps; 
      double lower;
      double mid1;
      double mid2;
      double upper;
      double maxf;
      double maxx;
      double f1;
      double f2;
      long m1;
      long flag;
};
}
#endif

