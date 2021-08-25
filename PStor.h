#ifndef PSTOR_H
#define PSTOR_H
#include <fstream>
#include <iostream>
#include <runn.h>
#include <FBase.h>
#include <DataObj.h>
#include <DStor.h>
#include <SStor.h>
#include <NStor.h>
#include <sys/mman.h>
#include <set>
#include <string>
#include <map>

using namespace std;
namespace iret {

typedef long* pLong;

class PStor : public TStor {
public:
   PStor(const char *nm,const char *path_nam); //nm name of object, nd dimension of data
      //path_nam is pointer at a string sss and reads the path from file path_sss in
      //current directory. But if sss begins with ':' then skips this character and 
      //remaining string is the path string itself.
   ~PStor(void);
   void create_Copy(const char *pnam); //Path in path_*pnam where the copy will be written

   //DStor must be of type "med"
   void create_PStor(DStor &Ds); //Takes from all of Ds to populate the data
   void create_PStor(DStor &Ds, Index *pPmid); //Takes from the subset only
   void update_PStor(DStor &Ds); //Takes from all of Ds to populate the data
   void update_PStor(DStor &Ds, Index *pPmid); //Takes from the subset to update
   //PStor records counts for phrases of two or more tokens without stop words or
   //punctuation. Includes only docs with abstracts. The 6 numbers stored for a phrase are
   //n1 = total number of occurrences of phrase in all documents
   //n2 = number of times phrase appears missing left token in same docs with full phrase
   //n3 = number of times phrase appears missing right token in same docs with full phrase
   //n4 = number of times phrase appears as part of longer phrase to left 
   //n5 = number of times phrase appears as part of longer phrase to right
   //n6 = number of docs phrase appears in

   void gopen_AccessN(void); //Opens NStor for access
   bool find(long id); //returns true if id present
   long total(void); //returns the number of documents entered into database
   void gclose_AccessN(void); //Closes NStor

   NStor *pNs; //pointer at NStor object.
};
   
}
#endif
