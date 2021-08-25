#ifndef QSTOR_H
#define QSTOR_H
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

class QStor : public TStor {
public:
   QStor(const char *nm,const char *path_nam); //nm name of object, nd dimension of data
      //path_nam is pointer at a string sss and reads the path from file path_sss in
      //current directory. But if sss begins with ':' then skips this character and 
      //remaining string is the path string itself.
   ~QStor(void);
   void create_Copy(const char *pnam); //Path in path_*pnam where the copy will be written

   //DStor must be of type "doc"
   void create_QStor(DStor &Ds); //Takes from all of Ds to populate the data
   void create_QStor(DStor &Ds, Index *pPmid); //Takes from the subset only
   void update_QStor(DStor &Ds); //Takes from all of Ds to populate the data
   void update_QStor(DStor &Ds, Index *pPmid); //Takes from the subset to update
   //QStor collects counts for the terms used in XPost for MEDLINE. Thus terms have
   //"!!t", "!!p" type endings, etc. Includes MeSH terms in the counts.
   //for a terms records counts from all documents processed even without abstracts
   //n1 = number of documents in which terms occurs
   //n2 = total number of counts for the term in all documents

   void gopen_AccessN(void); //Opens NStor for access
   bool find(long id); //returns true if id is present
   long total(void); //Returns number of documents entered into database
   void gclose_AccessN(void); //Closes NStor

   NStor *pNs; //pointer at NStor object.
};
   
}
#endif
