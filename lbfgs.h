#ifndef LBFGS_H
#define LBFGS_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <ap.h>
extern void funcgrad(ap::real_1d_array x, double& f, ap::real_1d_array& g);
using namespace std;
namespace iret {
    void lbfgsminimize(const int& n,
     const int& m,
     ap::real_1d_array& x,
     const double& epsg,
     const double& epsf,
     const double& epsx,
     const int& maxits,
     int& info);
void lbfgslincomb(const int& n,
     const double& da,
     const ap::real_1d_array& dx,
     int sx,
     ap::real_1d_array& dy,
     int sy);
double lbfgsdotproduct(const int& n,
     const ap::real_1d_array& dx,
     int sx,
     const ap::real_1d_array& dy,
     int sy);
void lbfgsmcsrch(const int& n,
     ap::real_1d_array& x,
     double& f,
     ap::real_1d_array& g,
     const ap::real_1d_array& s,
     int sstart,
     double& stp,
     const double& ftol,
     const double& xtol,
     const int& maxfev,
     int& info,
     int& nfev,
     ap::real_1d_array& wa,
     const double& gtol,
     const double& stpmin,
     const double& stpmax);
void lbfgsmcstep(double& stx,
     double& fx,
     double& dx,
     double& sty,
     double& fy,
     double& dy,
     double& stp,
     const double& fp,
     const double& dp,
     bool& brackt,
     const double& stmin,
     const double& stmax,
     int& info);
void lbfgsnewiteration(const ap::real_1d_array& x,
     double f,
     const ap::real_1d_array& g);
}
#endif
