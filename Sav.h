#ifndef SAV_H
#define SAV_H
#include <fstream>
#include <iostream>
#include <FBase.h>
#include <string>
#include <runn.h>
using namespace std;
namespace iret {

class Sav : public FBase {
   public:
      Sav(const char *skey,const char *pnam); //created with key and path name
     ~Sav();
      void readf();
      void writef(char *sxx);
      void getval(void); //querys the user for a val
      string *val; //Place to store the value
};
}
#endif
