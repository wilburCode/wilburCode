#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <runn.h>
#include "GoldOpt.h"
#define R 0.61803399
#define C (1.0-R)
using namespace std;

namespace iret{

GoldOpt::GoldOpt(){
}

GoldOpt::~GoldOpt(){
}

void GoldOpt::set_init_bracket(double a, double c){
   double xx=C*(c-a);
   lower=a;
   upper=c;
   mid1=a+xx;
   mid2=c-xx;
   flag=2;
   maxf=-1000000000;
}

void GoldOpt::set_limit(double xx){
   eps=xx;
}

double GoldOpt::give_value(){
   switch(flag){
      case 0: return(mid1);
      case 1: return(mid2);
      case 2: return(mid1);
      case 3: return(mid1);
   }
}

void GoldOpt::return_value(double w){
   switch(flag){
      case 0: f1=w;
              if(f1>maxf){
                 maxf=f1;
                 maxx=mid1;
              }
              if(f1<f2){
                 lower=mid1;
                 mid1=mid2;
                 f1=f2;
                 mid2=upper-C*(upper-lower);
                 flag=1;
              }
              else if(f2<f1){
                 upper=mid2;
                 mid2=mid1;
                 f2=f1; 
                 mid1=lower+C*(upper-lower);
                 flag=0;
              }
              else {
                 mid1=0.5*(mid1+mid2);
                 flag=3;
              }
              break;
      case 1: f2=w;
              if(f2>maxf){
                 maxf=f2;
                 maxx=mid2;
              }
              if(f1<f2){
                 lower=mid1;
                 mid1=mid2;
                 f1=f2;
                 mid2=upper-C*(upper-lower);
                 flag=1;
              }
              else if(f2<f1){
                 upper=mid2;
                 mid2=mid1;
                 f2=f1;
                 mid1=lower+C*(upper-lower);
                 flag=0;
              }
              else {
                 mid1=0.5*(mid1+mid2);
                 flag=3;
              }
              break;
      case 2: f1=w;
              if(f1>maxf){
                 maxf=f1;
                 maxx=mid1;
              }
              flag=1;
              break;
      case 3: if(w>maxf){
                 maxf=w;
                 maxx=mid1;
              }
              if(f2>=w){
                 flag=4;
              }
              else {
                mid1=lower+C*(upper-lower);
                lower=mid1;
                upper=mid2;
                mid1=lower+C*(upper-lower);
                mid2=upper-C*(upper-lower);
                flag=2;
              }
              break;
   }
}

int GoldOpt::next(){
   if((upper-lower<eps)||(flag==4))return(0);
   else return(1);
}

}
