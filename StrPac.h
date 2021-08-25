#ifndef STRPAC_H
#define STRPAC_H
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <runn.h>
#include <DStor.h>

using namespace std;
namespace iret {

class StrPac {
   public:
      StrPac(void); 
     ~StrPac(void);
      void add_str(const char *str,const char typ); 
      void add_str(string &str,const char typ); 
      void pack(DSpan &Ds); //Packs the contents into Ds
         //cln must be <=size of Ds.buf array for this
         //If same Ds is used on repeated calls this
         //condition will be maintained.
      int unpack(DSpan &Ds); //Unpacks the contents of Ds
      void clear(void); //empties the vectors
      void show_data(void); //Debug function
         //Prints data to stdout
      long id; //Record id
      long date; //Record date
      vector<char>   sc;
      vector<string> ss;
};
}
#endif
