#ifndef XPOST_H
#define XPOST_H
#include <cmath>
#include <utility>
#include <FBase.h>
#include <Strset.h>
#include <Dmap.h>
#include <Hash.h>
#include <Doc.h>
#include <DStor.h>
#include <BagPac.h>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <queue>

using namespace std;
namespace iret {

extern double xttrc;

class MpTh {
   public:

   MpTh(int);
   ~MpTh(void);
   //Data
   std::mutex g_lockproc;
   std::mutex g_lockque;
   std::condition_variable g_proccheck;
   std::condition_variable g_quecheck;
   bool g_done;
   bool g_filend;
   bool g_loopend;
   long g_cnt;
   std::queue<int> g_indices;
   int nthr;
   int *g_fin;
   long *iw;
   long *ic;
   float *sca;
};

extern class MpTh mt;

void betaProc(int k,DSpan *pDn,strMap *pMp);

void betdProc(int k,string *psxc,strMap *pMp);
 
void mergProc(int k,strMap *pMp,strMap *pNp);

void alphProc(int k,Chash *pCh,DSpan *pDn,int *wrd);

void alphzProc(int k,Chash *pCh,DSpan *pDn,int *wrd,float *xnm);

void alphzfProc(int k,Chash *pCh,DSpan *pDn,int *wrd,float *xnm,float (*d_local)(int,long));

void alpdProc(int k,Chash *pCh,string *psxc,int *addr);

void alpdzProc(int k,Chash *pCh,string *psxc,int *wrd,float *xnm);

void alpdzfProc(int k,Chash *pCh,string *psxc,int *wrd,float *xnm,float (*d_local)(int,long));

class STerm : public FBase {
   public:
      STerm(void); //For copies
      STerm(const char *typ,const char *nam); //typ is type of parent class
      //nam is name of the parent class
      STerm(const char *typ,int n,const char *nam); //typ is type of parent 
      //class, n is appended to end of typ
      //nam is name of the parent class
      STerm(const char *typ,const char *nam,const char *pnam); //typ is type 
      //of parent class nam is name of the parent class
         //pnam is either used in "path_pnam" as place to find path
         //or if begins with ':' is followed by path itself.
     ~STerm();

      void create(strMap &Mp);
      Lexos *gopen_Lexos(void); //Opens the Lexos maps for use
      void gclose_Lexos(void); //Closes the Lexos maps for use
      Chash *gopen_Chash(void); //Opens the Chash maps for use
      void gclose_Chash(void); //Closes the Chash maps for use
      void SetMem(STerm &St); //makes copies of data and memory of St data and memory
      Lexos *pLx; //Points at the lexically order term list object
      Chash *pCh; //Points at the Hash table that can return a term number
};

template<class Y,class Z>
class XPost : public FBase {
   public:
      XPost(void); //For copies
      XPost(const char *nm); //nm is name 
      XPost(const char *nm,const char *pnam); //nm is name 
         //pnam is either used in "path_pnam" as place to find path
         //or if begins with ':' is followed by path itself.
     ~XPost();
      void SetMem(XPost<Y,Z> &Xp); //Sets data and memory to Xp data and memory
      //Create all structures
      void create_Alln(Doc<Z> &SDc); 
         //Creates all files for XPost, postings data only
      void create_Allz(Doc<Z> &SDc);
         //Creates all files for XPost, postings and local Z values 
      void create_AllnM(Doc<Z> &SDc); //Multithreaded version
         //Creates all files for XPost, postings data only
      void create_AllzM(Doc<Z> &SDc); //Multithreaded version
         //Creates all files for XPost, postings and local Z values 
      void create_AllzM(Doc<Z> &SDc,float (*d_local)(int,long)); //Multithreaded version
         //Creates all files for XPost, postings and local Z values 
      void create_AllzG(Doc<Z> &SDc,float (*global)(long));
         //Creates all files for XPost, postings and local Z values 
         //Includes sqrt of global weight in local weight 
      void create_AllzG(Doc<Z> &SDc,float (*global)(long),float (*d_local)(int,long));
         //Creates all files for XPost, postings and local Z values 
         //Includes sqrt of global weight in local weight 
      void create_AllzNorm(Doc<Z> &SDc); //As previous but normalized
      void create_AllzNormG(Doc<Z> &SDc,float (*global)(long),float (*d_local)(int,long)); 
         //As create_AllzG but normalized
         //
      void create_Alln(DStor &Ds); 
         //Creates all files for XPost, postings data only
      void create_Allz(DStor &Ds);
         //Creates all files for XPost, postings and local Z values 
      void create_Allz(DStor &Ds,float (*d_local)(int,long));
         //Creates all files for XPost, postings and local Z values 
      void create_AllnM(DStor &Ds); //Multithreaded version 
         //Creates all files for XPost, postings data only
      void create_AllzM(DStor &Ds); //Multithreaded version
         //Creates all files for XPost, postings and local Z values 
      void create_AllzM(DStor &Ds,float (*d_local)(int,long)); //Multithreaded version
         //Creates all files for XPost, postings and local Z values 
      void create_AllzG(DStor &Ds,float (*global)(long),float (*d_local)(int,long));
         //Creates all files for XPost, postings and local Z values 
         //Includes sqrt of global weight in local weight 
      void create_AllzNorm(DStor &Ds,float (*d_local)(int,long)); 
         //As create_Allz with local Z values, but normalized
      void create_AllzNormG(DStor &Ds,float (*global)(long),float (*d_local)(int,long)); 
         //As create_AllzG but normalized
      //
      //Relative forms. Limits creation to pDnd subset (pDnd is list of indices)
      void create_Alln(Doc<Z> &SDc,Indx<Y> *pDnd); 
         //Creates all files for XPost, postings data only
      void create_Allz(Doc<Z> &SDc,Indx<Y> *pDnd);
         //Creates all files for XPost, postings and local Z values 
      void create_AllzG(Doc<Z> &SDc,Indx<Y> *pDnd,float (*global)(long));
         //Creates all files for XPost, postings and local Z values 
      void create_AllzNorm(Doc<Z> &SDc,Indx<Y> *pDnd); //As previous but normalized
         //
      void create_Alln(DStor &Ds,Index *pDnd); //pDnd is list of pmids
         //Creates all files for XPost, postings data only
      void create_Allz(DStor &Ds,Index *pDnd);
         //Creates all files for XPost, postings and local Z values 
      void create_Allz(DStor &Ds,Index *pDnd,float (*d_local)(int,long));
         //Creates all files for XPost, postings and local Z values 
      void create_AllzG(DStor &Ds,Index *pDnd,float (*global)(long),float (*d_local)(int,long));
         //Creates all files for XPost, postings and local Z values 
         //Includes sqrt of global weight in local weight 
      void create_AllzNorm(DStor &Ds,Index *pDnd,float (*d_local)(int,long)); 
         //As create_Allz with local Z values, but normalized
      void create_AllzNormG(DStor &Ds,Index *pDnd,float (*global)(long),float (*d_local)(int,long)); 
         //As create_AllzG but normalized

      //Term structures
      void create_Terms(Doc<Z> &SDc); //Reads SDc and creates 
         //all the term structures
      void create_Terms(Doc<Z> &SDc,Indx<Y> *pDnd); //Reads Ds and creates 
         //all the term structures. Limits creation to pDnd subset
      void create_TermsM(Doc<Z> &SDc); //Reads SDc and creates 
         //all the term structures. Multithreaded version.
      void create_Terms(DStor &Ds); //Reads SDc and creates 
         //all the term structures
      void create_Terms(DStor &Ds,Index *pDnd); //Reads Ds and creates 
         //all the term structures. Limits creation to pDnd subset
      void create_TermsM(DStor &Ds); //Reads SDc and creates 
         //all the term structures. Multithreaded version
      void record(Index *pDnd); //Here pDnd represents list of pmids
      Index *recover(void); //Reads back in the list of pmids
         //created by record
      void create_SqrtWts(float (*global)(long)); //Creates square roots of global
         //weights as .swx file
      void gopen_hash(void); //Opens pCh for use
      void gclose_hash(void);
      void gclose_lexos(void);
      long find(const char *str); //Returns term num+1 if found, else 0
      void gopen_lexos(void); //Opens pLx for use
      char *show(Y n); //Returns pointer at string numbered n in system
 
      //Binary document representation
      void create_DBinn(Doc<Z> &SDc); //Creates binary doc rep
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinnM(Doc<Z> &SDc); //Creates binary doc rep
         //Also creates a .z file with ndoc, max set size, nwrd. Multithreaded version.
      void create_DBinz(Doc<Z> &SDc); //Creates binary doc rep
         //creates local Z values in .x file
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinzM(Doc<Z> &SDc); //Creates binary doc rep
         //creates local Z values in .x file
         //Also creates a .z file with ndoc, max set size, nwrd 
         //Multithreaded version
      void create_DBinzG(Doc<Z> &SDc); //Creates binary doc rep
         //creates local Z values in .x file
         //Also creates a .z file with ndoc, max set size, nwrd 
         //Must call create_SqrtWts() function first and uses global result
      void create_DBinzM(Doc<Z> &SDc,float (*d_local)(int,long)); 
         //Creates binary doc rep as previous, but also local wts 
         //Also creates a .z file with ndoc, max set size, nwrd 
         //Multithreaded version
      void create_DBinzG(Doc<Z> &SDc,float (*d_local)(int,long)); 
         //Creates binary doc rep as previous, but also local wts times sqrt gl
         //Must call create_SqrtWts() function first and uses global result
      void create_DBinzNorm(Doc<Z> &SDc); //Normalized
      void create_DBinzNormG(Doc<Z> &SDc,float (*d_local)(int,long)); 
         //local wts, like create_DBinzG but normalized
         //
      void create_DBinn(DStor &Ds); //Creates binary doc rep
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinnM(DStor &Ds); //Creates binary doc rep
         //Also creates a .z file with ndoc, max set size, nwrd. Multithreaded version.
      void create_DBinz(DStor &Ds); //Creates binary doc rep
         //creates local Z values in .x file-just local counts
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinzM(DStor &Ds); //Creates binary doc rep
         //creates local Z values in .x file-just local counts
         //Also creates a .z file with ndoc, max set size, nwrd 
         //Multithreaded version
      void create_DBinz(DStor &Ds,float (*d_local)(int,long)); 
         //Creates binary doc rep like previous, but also makes local wts
      void create_DBinzM(DStor &Ds,float (*d_local)(int,long)); 
         //Creates binary doc rep like previous, but also makes local wts
      void create_DBinzG(DStor &Ds,float (*d_local)(int,long)); 
         //Creates binary doc rep as previous, but also local wts times sqrt gl
         //Must call create_SqrtWts() function first and uses global result
      void create_DBinzNorm(DStor &Ds,float (*d_local)(int,long)); 
         //local wts, like create_DBinz but normalized
      void create_DBinzNormG(DStor &Ds,float (*d_local)(int,long)); 
         //local wts, like create_DBinzG but normalized
         //
      //Relative forms. Limits creation to pDnd subset
      void create_DBinn(Doc<Z> &SDc,Indx<Y> *pDnd); //Creates binary doc rep
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinz(Doc<Z> &SDc,Indx<Y> *pDnd); //Creates binary doc rep
         //creates local Z values in .x file
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinzG(Doc<Z> &SDc,Indx<Y> *pDnd); //Creates binary doc rep
         //creates local Z values in .x file
         //Also creates a .z file with ndoc, max set size, nwrd 
         //Must call create_SqrtWts() function first and uses global result
      void create_DBinzNorm(Doc<Z> &SDc,Indx<Y> *pDnd); //Normalized
         //
      void create_DBinn(DStor &Ds,Index *pDnd); //Creates binary doc rep
         //Also creates a .z file with sndoc, max set size, nwrd
      void create_DBinz(DStor &Ds,Index *pDnd); //Creates binary doc rep
         //creates local Z values in .x file - just counts
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinz(DStor &Ds,Index *pDnd,float (*d_local)(int,long)); 
         //Creates binary doc rep
         //creates local Z values in .x file - local weights
         //Also creates a .z file with ndoc, max set size, nwrd 
      void create_DBinzG(DStor &Ds,Index *pDnd,float (*d_local)(int,long)); 
         //Creates binary doc rep
         //creates local Z values in .x file - local weights
         //Also creates a .z file with ndoc, max set size, nwrd 
         //Must call create_SqrtWts() function first and uses global result
      void create_DBinzNorm(DStor &Ds,Index *pDnd,float (*d_local)(int,long)); 
         //like create_DBinz with local weights, but also normalized
      void create_DBinzNormG(DStor &Ds,Index *pDnd,float (*d_local)(int,long)); 
         //like create_DBinzG with local weights, but also normalized
         //
      void gopen_db_map(void); //Maps the addr arrays
         //also maps size arrays.
         //Also reads in .z file with ndoc, nsets and nwrd in it
      void gopen_db(void); //Reads in the data that prev. func. maps
      Indx<Y> *readp_db(Y n); //Sets pointer at doc term indices
         //array for nth term. Sets term count
      void set_Indx(void); //Sets nw and nwd as values for Indx object
      Z *readz_db(Y n); //Sets pointer at array of 
         //local Z values 
      void gclose_db_map(void); //Unmaps the addr arrays and
         //the size arrays.

      //Postings
      void create_Postn(void); //Creates postings files 
         //Also creates a .z file with nwrd, nslcs and ndoc in it
      void create_Postz(void); //Creates files for Z values 
         //Also creates a .z file with nwrd, nslcs and ndoc in it
      void gopen_map(void); //Maps the addr arrays
         //also maps freq array.
         //Also reads in .z file with nwrd, nslcs and ndoc in it
      void gopen_rp(void); //Reads the addr arrays
         //also reads freq array.
         //Also reads in .z file with nwrd, nslcs and ndoc in it
      Indx<Y> *readp(Y n); //Returns pointer at postings
         //index object for nth term.
      Z *readz(Y n); //Returns a pointer at list of 
         //local Z values
      void gclose_map(void); //Unmaps freq and addr


      //GLOBAL DATA
      Y ndoc; //Number of docs in system
      Y nwrd; //Number of terms in system

      //DBIN DATA
      //Single doc read
      Y nw; //Number of terms in doc.
      Y *nwd; //Array of term numbers in doc.
      Z *lwt; //Array of local Z values in doc
      Indx<Y> *pTrm; //represents currently read doc if set_Indx() called
         //Created when gopen_db_map() called.

      //dbin data
      Y *nsiz; //Maps number of terms in each doc
      long *addrd; //Maps the addresses for docs
      Y *term; //Maps the term numbers for all docs
      Z *dw; //Maps local Z values for all docs

      //POSTINGS DATA
      //Single term read
      Indx<Y> *pInd; //frequency and postings
      Z *wz; //local Z values

      //postings data
      Y *freq; //Maps the frequency of terms
      Z *sgw; //Maps square root of global weights
      long *addr; //Maps the addresses for postings 
      Y *post; //Maps the postings file
      Z *wc; //Maps the local weights
      Lexos *pLx; //Points at the lexically order term list object
      Chash *pCh; //Points at the Hash table that can return a term number
      bool own_data_mem;   // true if owns memory maps, else false
      bool inMemPost; //true if postings data read in by gopen_rp.
      bool inMemDocs; //true if document data read in by gopen_db.
      STerm *pSt; //Points at the STerm object
};

template<class Y,class Z>
XPost<Y,Z>::XPost(void) : FBase("xpost","null"){
   open1=open2=open3=open4=0;
   pCh = nullptr;
   pLx = nullptr;
   pSt = new STerm;
   pTrm = nullptr;
   pInd = nullptr;
   own_data_mem = true;
   inMemPost = false;
   inMemDocs = false;
}

template<class Y,class Z>
XPost<Y,Z>::XPost(const char *nm) : FBase("xpost",nm){
   open1=open2=open3=open4=0;
   pCh = nullptr;
   pLx = nullptr;
   pSt=new STerm("xpost",nm);
   pTrm = nullptr;
   pInd = nullptr;
   own_data_mem = true;
   inMemPost = false;
   inMemDocs = false;
} 

template<class Y,class Z>
XPost<Y,Z>::XPost(const char *nm,const char *pnam) : FBase("xpost",nm,pnam){
   open1=open2=open3=open4=0;
   pCh = nullptr;
   pLx = nullptr;
   pSt=new STerm("xpost",nm,pnam);
   pTrm = nullptr;
   pInd = nullptr;
   own_data_mem = true;
   inMemPost = false;
   inMemDocs = false;
} 

template<class Y,class Z>
XPost<Y,Z>::~XPost(){

   if(own_data_mem){
      if(open1)
         gclose_db_map();
      if(open2)
         gclose_map();
   }
   delete pSt;
   if (pInd) {
      pInd->idx = nullptr;
      delete pInd;
   }
   if (pTrm) {
      pTrm->idx = nullptr;
      delete pTrm;
   }
}  

template<class Y,class Z>
void XPost<Y,Z>::SetMem(XPost<Y,Z> &Xp){
   if(Xp.open1){
      this->nwrd=Xp.nwrd;
      this->ndoc=Xp.ndoc;
      this->term=Xp.term;
      this->dw=Xp.dw;
      this->addrd=Xp.addrd;
      this->nsiz=Xp.nsiz;
      this->pTrm=new Indx<Y>;
      this->open1=1;
   }
   if(Xp.open2){
      this->freq=Xp.freq;
      this->addr=Xp.addr;
      this->post=Xp.post;
      this->wc=Xp.wc;
      this->pInd=new Indx<Y>;
      this->open2=1;
   }
   pSt->SetMem(*Xp.pSt);
   if(Xp.open3){
      this->pCh=pSt->pCh;
      this->open3=1;
   }
   if(Xp.open4){
      this->pLx=pSt->pLx;
      this->open4=1;
   }
   own_data_mem = false;
}

template<class Y,class Z>
void XPost<Y,Z>::create_Alln(Doc<Z> &SDc){
   create_Terms(SDc);
   create_DBinn(SDc);
   create_Postn();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Alln(DStor &Ds){
   create_Terms(Ds);
   create_DBinn(Ds);
   create_Postn();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Allz(Doc<Z> &SDc){
   create_Terms(SDc);
   create_DBinz(SDc);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Allz(DStor &Ds){
   create_Terms(Ds);
   create_DBinz(Ds);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Allz(DStor &Ds,float (*d_local)(int,long)){
   create_Terms(Ds);
   create_DBinz(Ds,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllnM(Doc<Z> &SDc){
   create_TermsM(SDc);
   create_DBinnM(SDc);
   create_Postn();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllnM(DStor &Ds){
   create_TermsM(Ds);
   create_DBinnM(Ds);
   create_Postn();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzM(Doc<Z> &SDc){
   create_TermsM(SDc);
   create_DBinzM(SDc);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzM(DStor &Ds){
   create_TermsM(Ds);
   create_DBinzM(Ds);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzM(Doc<Z> &SDc,float (*d_local)(int,long)){
   create_TermsM(SDc);
   create_DBinzM(SDc,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzM(DStor &Ds,float (*d_local)(int,long)){
   create_TermsM(Ds);
   create_DBinzM(Ds,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzG(Doc<Z> &SDc,float (*global)(long)){
   create_Terms(SDc);
   create_SqrtWts(global);
   create_DBinzG(SDc);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzG(Doc<Z> &SDc,float (*global)(long),float (*d_local)(int,long)){
   create_Terms(SDc);
   create_SqrtWts(global);
   create_DBinzG(SDc,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzG(DStor &Ds,float (*global)(long),float (*d_local)(int,long)){
   create_Terms(Ds);
   create_SqrtWts(global);
   create_DBinzG(Ds,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzNormG(DStor &Ds,float (*global)(long),float (*d_local)(int,long)){
   create_Terms(Ds);
   create_SqrtWts(global);
   create_DBinzNormG(Ds,d_local);
   create_Postz();
}
//new_code
template<class Y,class Z>
void XPost<Y,Z>::create_AllzNormG(Doc<Z> &SDc,float (*global)(long),float (*d_local)(int,long)){
   create_Terms(SDc);
   create_SqrtWts(global);
   create_DBinzNormG(SDc,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzNorm(Doc<Z> &SDc){
   create_Terms(SDc);
   create_DBinzNorm(SDc);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzNorm(DStor &Ds,float (*d_local)(int,long)){
   create_Terms(Ds);
   create_DBinzNorm(Ds,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Alln(Doc<Z> &SDc,Indx<Y> *pDnd){
   create_Terms(SDc,pDnd);
   create_DBinn(SDc,pDnd);
   create_Postn();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Alln(DStor &Ds,Index *pDnd){
   create_Terms(Ds,pDnd);
   create_DBinn(Ds,pDnd);
   create_Postn();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Allz(Doc<Z> &SDc,Indx<Y> *pDnd){
   create_Terms(SDc,pDnd);
   create_DBinz(SDc,pDnd);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Allz(DStor &Ds,Index *pDnd){
   create_Terms(Ds,pDnd);
   create_DBinz(Ds,pDnd);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Allz(DStor &Ds,Index *pDnd,float (*d_local)(int,long)){
   create_Terms(Ds,pDnd);
   create_DBinz(Ds,pDnd,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzG(Doc<Z> &SDc,Indx<Y> *pDnd,float (*global)(long)){
   create_Terms(SDc,pDnd);
   create_SqrtWts(global);
   create_DBinzG(SDc,pDnd);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzG(DStor &Ds,Index *pDnd,float (*global)(long),float (*d_local)(int,long)){
   create_Terms(Ds,pDnd);
   create_SqrtWts(global);
   create_DBinzG(Ds,pDnd,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzNormG(DStor &Ds,Index *pDnd,float (*global)(long),float (*d_local)(int,long)){
   create_Terms(Ds,pDnd);
   create_SqrtWts(global);
   create_DBinzNormG(Ds,pDnd,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzNorm(Doc<Z> &SDc,Indx<Y> *pDnd){
   create_Terms(SDc,pDnd);
   create_DBinzNorm(SDc,pDnd);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_AllzNorm(DStor &Ds,Index *pDnd,float (*d_local)(int,long)){
   create_Terms(Ds,pDnd);
   create_DBinzNorm(Ds,pDnd,d_local);
   create_Postz();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Terms(Doc<Z> &SDc){
   Y i,*px;
   long k,lxn;
   Z j;
   const char *pch;
   char *qch;
   strMap Mp;
   typename strMap::iterator q,qe=Mp.end();
   SDc.gopen_map();
   for(i=0;i<SDc.ndoc;i++){
      SDc.read(i);
      while((pch=SDc.show(j))){
         if((q=Mp.find(pch))!=qe){
            q->second+=1;
         }
         else {
            lxn=strlen(pch);
            qch=new char[lxn+1];
            strcpy(qch,pch);
            Mp.insert(make_pair(qch,1));
         }
      }
      SDc.clear();
      mark((long)i,10000,"documents in");
   }

   long n=(long)Mp.size(),m=(long)SDc.ndoc;
   put_Nnum("z",n,m);

   k=0;
   ofstream *pfout=get_Ostr("f",ios::out);
   q=Mp.begin();
   while(q!=qe){
      i=q->second;
      pfout->write((char*)&i,sizeof(Y));
      q++;
      mark(++k,10000,"terms");
   }     
   dst_Ostr(pfout);
   if(this->pflag)cout << Mp.size() << " terms in system" << endl;

   this->map_down((FBase *)pSt);
   pSt->create(Mp);
}

template<class Y,class Z>
void XPost<Y,Z>::create_TermsM(Doc<Z> &SDc){
   long i,sm;
   Y j,k,u,ct,flag,bn[mt.nthr],bx=0;
   int sz;
   strMap *pMp[mt.nthr];
   string sxc[mt.nthr];

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pMp[k]=new strMap;
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(betdProc,k,&sxc[k],pMp[k]));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;

   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfw=get_Ostr("wd",ios::out);
   sm=0;
   
   SDc.gopen_map();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               sxc[k].clear();
               if((sz=SDc.read(sxc[k]))){
                  sz--;
                  pfa->write((char*)&sm,sizeof(long));
                  pfs->write((char*)&sz,sizeof(int));
                  for(j=0;j<sz;j++){
                     pfw->write((char*)&bx,sizeof(Y));
                  }
                  sm+=sz;
                  mt.g_cnt++;
                  mt.iw[k]=1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }

   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfw);
   put_Nnum("sn",sm);
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   ct=2;
   for(k=0;k<mt.nthr;k++)mt.g_fin[k]=k;
   while(ct-1){
      for(k=0;k<mt.nthr;k++)bn[k]=0;
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         ct=0;
         for(k=0;k<mt.nthr;k++){
            if(mt.g_fin[k]>-1){
               bn[mt.g_fin[k]]++;
               ct++;
            }  
         }
      }
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(bn[k]==1){
            if(!flag){u=k;flag=1;}
            else {
               std::thread *pthr = new std::thread(mergProc,k,pMp[u],pMp[k]);
               pthr->detach();
               cout << "Merging " << u << " and " << k << endl;
               mt.g_fin[k]=mt.g_fin[u];
               flag=0;
            }
         }
      }
   }
   cout << "Merges Completed" << endl;
   long n=(long)pMp[0]->size(),m=(long)i;
   put_Nnum("z",n,m);

   k=0;
   ofstream *pfout=get_Ostr("f",ios::out);
   pMp[0]->Set();
   while(pMp[0]->qs!=pMp[0]->qz){
      i=pMp[0]->qs->second;
      pfout->write((char*)&i,sizeof(Y));
      pMp[0]->qs++;
      mark(++k,10000,"terms");
   }
   dst_Ostr(pfout);
   if(this->pflag)cout << pMp[0]->size() << " terms in system" << endl;

   this->map_down((FBase *)pSt);
   pSt->create(*pMp[0]);
   pMp[0]->~Dmap();
   for (k=1;k<mt.nthr;k++) delete pMp[k];
   SDc.gclose_map();
}

template<class Y,class Z>
void XPost<Y,Z>::create_Terms(DStor &Ds){
   Y i,*px;
   long k,lxn;
   Z j;
   const char *pch;
   char *qch;
   strMap Mp;
   typename strMap::iterator q,qe=Mp.end();
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();
   i=0;
   while(Ds.getNext()){
      Ds.read(Dn);
      Bp.unpack(Dn);
      vector<string>::iterator p=Bp.ss.begin();
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         if((q=Mp.find(pch))!=qe){
            q->second+=1;
         }
         else {
            lxn=strlen(pch);
            qch=new char[lxn+1];
            strcpy(qch,pch);
            Mp.insert(make_pair(qch,1));
         }
         p++;
      }
      Bp.clear();
      mark((long)(++i),10000,"documents in");
   }

   long n=(long)Mp.size(),m=(long)i;
   put_Nnum("z",n,m);

   k=0;
   ofstream *pfout=get_Ostr("f",ios::out);
   q=Mp.begin();
   while(q!=qe){
      i=q->second;
      pfout->write((char*)&i,sizeof(Y));
      q++;
      mark(++k,10000,"terms");
   }
   dst_Ostr(pfout);
   if(this->pflag)cout << Mp.size() << " terms in system" << endl;

   this->map_down((FBase *)pSt);
   pSt->create(Mp);
}

template<class Y,class Z>
void XPost<Y,Z>::create_TermsM(DStor &Ds){
   long i,sm;
   Y j,k,u,ct,flag,bn[mt.nthr],bx=0;
   int *pi,sz;
   DSpan Dn[mt.nthr];
   strMap *pMp[mt.nthr];

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pMp[k]=new strMap;
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(betaProc,k,&Dn[k],pMp[k]));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;
   
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfw=get_Ostr("wd",ios::out);
   sm=0;
  
   Ds.gopen_read();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               if(Ds.getNext()){
                  Ds.read(Dn[k]);
                  pi=(int*)Dn[k].buf;
                  sz=pi[0];
                  pfa->write((char*)&sm,sizeof(long));
                  pfs->write((char*)&sz,sizeof(int));
                  for(j=0;j<sz;j++){
                     pfw->write((char*)&bx,sizeof(Y));
                  }
                  sm+=sz;
                  mt.g_cnt++;
                  mt.iw[k]=1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfw);
   put_Nnum("sn",sm);
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   ct=2;
   for(k=0;k<mt.nthr;k++)mt.g_fin[k]=k;
   while(ct-1){
     //     cerr << ct << '\n';
      for(k=0;k<mt.nthr;k++)bn[k]=0;
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         ct=0;
         for(k=0;k<mt.nthr;k++){
            if(mt.g_fin[k]>-1){
               bn[mt.g_fin[k]]++;
               ct++;
            }  
         }
      }
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(bn[k]==1){
            if(!flag){u=k;flag=1;}
            else {
               std::thread *pthr = new std::thread(mergProc,k,pMp[u],pMp[k]);
               pthr->detach();
               cout << "Merging " << u << " and " << k << endl;
               mt.g_fin[k]=mt.g_fin[u];
               flag=0;
            }
         }
      }
   }
   cout << "Merges Completed" << endl;
   long n=(long)pMp[0]->size(),m=(long)i;
   put_Nnum("z",n,m);

   k=0;
   ofstream *pfout=get_Ostr("f",ios::out);
   pMp[0]->Set();
   while(pMp[0]->qs!=pMp[0]->qz){
      i=pMp[0]->qs->second;
      pfout->write((char*)&i,sizeof(Y));
      pMp[0]->qs++;
      mark(++k,10000,"terms");
   }
   dst_Ostr(pfout);
   if(this->pflag)cout << pMp[0]->size() << " terms in system" << endl;

   this->map_down((FBase *)pSt);
   pSt->create(*pMp[0]);
   pMp[0]->~Dmap();
   for (k=0;k<mt.nthr;k++) delete pMp[k];
}

template<class Y,class Z>
void XPost<Y,Z>::create_Terms(Doc<Z> &SDc,Indx<Y> *pDnd){
   Y i,*px;
   long k,lxn;
   Z j;
   const char *pch;
   char *qch;
   strMap Mp;
   typename strMap::iterator q,qe=Mp.end();
   SDc.gopen_map();
   for(i=0;i<pDnd->ix;i++){
      SDc.read(pDnd->idx[i]);
      while((pch=SDc.show(j))){
         if((q=Mp.find(pch))!=qe){
            q->second+=1;
         }
         else {
            lxn=strlen(pch);
            qch=new char[lxn+1];
            strcpy(qch,pch);
            Mp.insert(make_pair(qch,1));
         }
      }
      SDc.clear();
      mark((long)i,10000,"documents in");
   }

   long n=(long)Mp.size(),m=(long)pDnd->ix;
   put_Nnum("z",n,m);

   k=0;
   ofstream *pfout=get_Ostr("f",ios::out);
   q=Mp.begin();
   while(q!=qe){
      i=q->second;
      pfout->write((char*)&i,sizeof(Y));
      q++;
      mark(++k,10000,"terms");
   }     
   dst_Ostr(pfout);
   if(this->pflag)cout << Mp.size() << " terms in system" << endl;

   this->map_down((FBase *)pSt);
   pSt->create(Mp);
}

template<class Y,class Z>
void XPost<Y,Z>::create_Terms(DStor &Ds,Index *pDnd){
   Y i,*px;
   long k,lxn;
   Z j;
   const char *pch;
   char *qch;
   strMap Mp;
   typename strMap::iterator q,qe=Mp.end();
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();
   for(i=0;i<pDnd->ix;i++){
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      vector<string>::iterator p=Bp.ss.begin();
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         if((q=Mp.find(pch))!=qe){
            q->second+=1;
         }
         else {
            lxn=strlen(pch);
            qch=new char[lxn+1];
            strcpy(qch,pch);
            Mp.insert(make_pair(qch,1));
         }
         p++;
      }
      Bp.clear();
      mark((long)i,10000,"documents in");
   }

   long n=(long)Mp.size(),m=(long)pDnd->ix;
   put_Nnum("z",n,m);

   k=0;
   ofstream *pfout=get_Ostr("f",ios::out);
   q=Mp.begin();
   while(q!=qe){
      i=q->second;
      pfout->write((char*)&i,sizeof(Y));
      q++;
      mark(++k,10000,"terms");
   }
   dst_Ostr(pfout);
   if(this->pflag)cout << Mp.size() << " terms in system" << endl;

   this->map_down((FBase *)pSt);
   pSt->create(Mp);
   record(pDnd);
}

template<class Y,class Z>
void XPost<Y,Z>::record(Index *pDnd){
   ofstream *pfmm=get_Ostr("pmid");
   pDnd->write(*pfmm);
   dst_Ostr(pfmm);
}

template<class Y,class Z>
Index *XPost<Y,Z>::recover(void){
   ifstream *pfmm=get_Istr("pmid");
   Index *pDnd = new Index;
   pDnd->read(*pfmm);
   dst_Istr(pfmm);
   return(pDnd);
}

template<class Y,class Z>
void XPost<Y,Z>::create_SqrtWts(float (*global)(long)){
   Y i,j,k;
   Z *zx;
   get_Nnum("z",i,j);
   nwrd=i;
   ndoc=j;
   xttrc=(double)ndoc;
   freq=(Y*)get_Mmap("f");
   zx=new Z[nwrd];
   for(i=0;i<nwrd;i++){
      zx[i]=sqrt((double)global((long)freq[i]));
   }
   bin_Writ("swx",nwrd*sizeof(Z),(char*)zx);
   delete [] zx;
}

//Doc Binary Functions

template<class Y,class Z>
void XPost<Y,Z>::create_DBinn(Doc<Z> &SDc){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)SDc.ndoc;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);

   sm=0;
   for(i=0;i<SDc.ndoc;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(i);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         pfw->write((char*)&k,sizeof(Y));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
}         

template<class Y,class Z>
void XPost<Y,Z>::create_DBinnM(Doc<Z> &SDc){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   int k,flag,u;
   Y *bx,*wrd;
   string sxc[mt.nthr];

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;

   adr=(long*)get_Mmap("ad");
   wrd=(Y*)get_Wmap("wd");
   Chash *pCx[mt.nthr];
   gopen_hash();

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pCx[k]=new Chash;
      pCx[k]->SetMem(*(pCh));
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(alpdProc,k,pCx[k],&sxc[k],wrd));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;

   SDc.gopen_map();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               if(SDc.read(sxc[k])){
                  mt.g_cnt++;
                  mt.ic[k]=adr[i];
                  mt.iw[k]=1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   SDc.gclose_map();
   this->dst_Mmap("ad",(char*&)adr);
   this->mak_Msync("wd",(char*)wrd);
   this->dst_Mmap("wd",(char*&)wrd);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinn(DStor &Ds){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);


   sm=0;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();
   i=0;
   while(Ds.getNext()){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      while(p!=Bp.ss.end()){
         k=(Y)pCh->count((*p).c_str())-1;
         pfw->write((char*)&k,sizeof(Y));
         p++;
      }
      Bp.clear();
      mark(i++,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinnM(DStor &Ds){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   int *pi,flag;
   Y k,*bx,*wrd;
   Z u;
   const char *pch;

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;

   adr=(long*)get_Mmap("ad");
   wrd=(Y*)get_Wmap("wd");
   DSpan Dn[mt.nthr];
   Chash *pCx[mt.nthr];
   gopen_hash();

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pCx[k]=new Chash;
      pCx[k]->SetMem(*pCh);
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(alphProc,k,pCx[k],&Dn[k],wrd));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;
   
   Ds.gopen_read();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               if(Ds.getNext()){
                  Ds.read(Dn[k]);
                  mt.g_cnt++;
                  mt.ic[k]=adr[i];
                  mt.iw[k]=i+1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   this->dst_Mmap("ad",(char*&)adr);
   this->mak_Msync("wd",(char*)wrd);
   this->dst_Mmap("wd",(char*&)wrd);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinn(Doc<Z> &SDc,Indx<Y> *pDnd){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)pDnd->ix;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);

   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(pDnd->idx[i]);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         pfw->write((char*)&k,sizeof(Y));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinn(DStor &Ds,Index *pDnd){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);

   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      while(p!=Bp.ss.end()){
         k=(Y)pCh->count((*p).c_str())-1;
         pfw->write((char*)&k,sizeof(Y));
         p++;
      }
      Bp.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinz(Doc<Z> &SDc){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)SDc.ndoc;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   for(i=0;i<SDc.ndoc;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(i);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}         

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzM(Doc<Z> &SDc){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   int k,flag;
   Y *bx,*wrd;
   Z u,*xnm;
   string sxc[mt.nthr];

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;

   get_Nnum("sn",sm);
   ofstream *pfx=get_Ostr("xd",ios::out);
   u=0;
   for(i=0;i<sm;i++){
      pfx->write((char*)&u,sizeof(Z));
   }
   dst_Ostr(pfx);

   adr=(long*)get_Mmap("ad");
   wrd=(Y*)get_Wmap("wd");
   xnm=(Z*)get_Wmap("xd");
   Chash *pCx[mt.nthr];
   gopen_hash();

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pCx[k]=new Chash;
      pCx[k]->SetMem(*pCh);
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(alpdzProc,k,pCx[k],&sxc[k],wrd,xnm));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;

   SDc.gopen_map();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               if(SDc.read(sxc[k])){
                  mt.g_cnt++;
                  mt.ic[k]=adr[i];
                  mt.iw[k]=1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   SDc.gclose_map();
   this->dst_Mmap("ad",(char*&)adr);
   this->mak_Msync("wd",(char*)wrd);
   this->dst_Mmap("wd",(char*&)wrd);
   this->mak_Msync("xd",(char*)xnm);
   this->dst_Mmap("xd",(char*&)xnm);
   for(i=0;i<mt.nthr;i++){
      delete pCx[i];
   }
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzM(Doc<Z> &SDc,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   int k,flag;
   Y *bx,*wrd;
   Z u,*xnm;
   string sxc[mt.nthr];

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;

   get_Nnum("sn",sm);
   ofstream *pfx=get_Ostr("xd",ios::out);
   u=0;
   for(i=0;i<sm;i++){
      pfx->write((char*)&u,sizeof(Z));
   }
   dst_Ostr(pfx);

   adr=(long*)get_Mmap("ad");
   wrd=(Y*)get_Wmap("wd");
   xnm=(Z*)get_Wmap("xd");
   Chash *pCx[mt.nthr];
   gopen_hash();

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pCx[k]=new Chash;
      pCx[k]->SetMem(*pCh);
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(alpdzfProc,k,pCx[k],&sxc[k],wrd,xnm,d_local));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;

   SDc.gopen_map();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               if(SDc.read(sxc[k])){
                  mt.g_cnt++;
                  mt.ic[k]=adr[i];
                  mt.iw[k]=1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   SDc.gclose_map();
   this->dst_Mmap("ad",(char*&)adr);
   this->mak_Msync("wd",(char*)wrd);
   this->dst_Mmap("wd",(char*&)wrd);
   this->mak_Msync("xd",(char*)xnm);
   this->dst_Mmap("xd",(char*&)xnm);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinz(DStor &Ds){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read(); 

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   i=0;
   while(Ds.getNext()){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         k=(Y)pCh->count((*p).c_str())-1;
         u=(Z)(*q);
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
         p++;
         q++;
      }
      Bp.clear();
      mark(i++,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzM(DStor &Ds){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   int *pi,flag;
   Y k,*bx,*wrd;
   Z u,*xnm;
   const char *pch;

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;

   get_Nnum("sn",sm);
   ofstream *pfx=get_Ostr("xd",ios::out);
   u=0;
   for(i=0;i<sm;i++){
      pfx->write((char*)&u,sizeof(Z));
   }
   dst_Ostr(pfx);

   adr=(long*)get_Mmap("ad");
   wrd=(Y*)get_Wmap("wd");
   xnm=(Z*)get_Wmap("xd");
   DSpan Dn[mt.nthr];
   Chash *pCx[mt.nthr];
   gopen_hash();

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pCx[k]=new Chash;
      pCx[k]->SetMem(*pCh);
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(alphzProc,k,pCx[k],&Dn[k],wrd,xnm));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;
   
   Ds.gopen_read();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               if(Ds.getNext()){
                  Ds.read(Dn[k]);
                  mt.g_cnt++;
                  mt.ic[k]=adr[i];
                  mt.iw[k]=i+1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   this->dst_Mmap("ad",(char*&)adr);
   this->mak_Msync("wd",(char*)wrd);
   this->dst_Mmap("wd",(char*&)wrd);
   this->mak_Msync("xd",(char*)xnm);
   this->dst_Mmap("xd",(char*&)xnm);

   for(k=0;k<mt.nthr;k++){
     // these values must be freed
     // hopefully they handle their values properly, because they are
     // copies
     delete pCh[k];
   }
   gclose_hash();
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzM(DStor &Ds,float (*d_local)(int,long)){
   long sm,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   int *pi,flag;
   Y k,*bx,*wrd;
   Z u,*xnm;
   const char *pch;

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;

   get_Nnum("sn",sm);
   ofstream *pfx=get_Ostr("xd",ios::out);
   u=0;
   for(i=0;i<sm;i++){
      pfx->write((char*)&u,sizeof(Z));
   }
   dst_Ostr(pfx);

   adr=(long*)get_Mmap("ad");
   wrd=(Y*)get_Wmap("wd");
   xnm=(Z*)get_Wmap("xd");
   DSpan Dn[mt.nthr];
   Chash *pCx[mt.nthr];
   gopen_hash();

   mt.g_done=false;
   mt.g_cnt=mt.nthr;
   std::vector<std::thread> threads;
   for(k=0;k<mt.nthr;k++){
      pCx[k]=new Chash;
      pCx[k]->SetMem(*pCh);
      mt.iw[k]=0;
      mt.g_fin[k]=0;
      threads.push_back(std::thread(alphzfProc,k,pCx[k],&Dn[k],wrd,xnm,d_local));
   }
   for(k=0;k<mt.nthr;k++)mt.g_indices.push(k);
   mt.g_loopend=false;
   
   Ds.gopen_read();
   i=0;

   while((!mt.g_loopend)||mt.g_cnt){
      {
         std::unique_lock<std::mutex> lck(mt.g_lockque); //Acquire a lock on mt.g_lockque, but unlocks mt.g_lockque on destruction
         while(!mt.g_indices.empty()){
            k=mt.g_indices.front();
            mt.g_indices.pop();
            mt.g_cnt--;
            flag=0;
            while((!mt.g_loopend)&&(!flag)){
               if(Ds.getNext()){
                  Ds.read(Dn[k]);
                  mt.g_cnt++;
                  mt.ic[k]=adr[i];
                  mt.iw[k]=i+1;
                  flag=1;
                  mark(++i,10000,"documents in");
               }
               else mt.g_loopend=true;
            }
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
      if(mt.g_cnt){
         std::unique_lock<std::mutex> lck(mt.g_lockque);
         mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
         mt.g_quecheck.wait(lck);//lock created in previous line blocks this thread and then is unlocked
      }
   }
   mt.g_done=true;
   mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   flag=1;
   while(flag){
      flag=0;
      for(k=0;k<mt.nthr;k++){
         if(!mt.g_fin[k]){
            flag=1;
            cout << "Waiting on thread " << k << endl;
            break;
         }
      }
      mt.g_proccheck.notify_all(); //Notify any thread waiting on mt.g_proccheck
   }
   for(auto& t : threads)t.join();
   cout << "Joined successfully" << endl;
   this->dst_Mmap("ad",(char*&)adr);
   this->mak_Msync("wd",(char*)wrd);
   this->dst_Mmap("wd",(char*&)wrd);
   this->mak_Msync("xd",(char*)xnm);
   this->dst_Mmap("xd",(char*&)xnm);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinz(DStor &Ds,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,vi;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,sum;
   const char *pch;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   i=0;
   while(Ds.getNext()){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      sum=0;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         p++;
         q++;
      }
      vi=rnd(sum);
      q=Bp.si.begin();
      while(q!=Bp.si.end()){
         u=(Z)d_local(rnd(*q),sum);
         pfx->write((char*)&u,sizeof(Z));
         q++;
      }
      Bp.clear();
      mark(i++,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzNorm(DStor &Ds,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,sum;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,ssq;
   const char *pch;
   vector<Z> vz;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   i=0;
   while(Ds.getNext()){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      sum=0;
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         p++;
         q++;
      }
      q=Bp.si.begin();
      ssq=0;
      while(q!=Bp.si.end()){
         u=(Z)d_local(*q,sum);
         vz.push_back(u);
         ssq+=u*u;
         q++;
      }
      ssq=(Z)sqrt((double)ssq);
      typename std::vector<Z>::iterator s=vz.begin();
      while(s!=vz.end()){
         u=*s/ssq;
         pfx->write((char*)&u,sizeof(Z));
         s++;
      }
      Bp.clear();
      vz.clear();
      mark(i++,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzG(DStor &Ds,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,sum;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;
   vector<long> vb;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   i=0;
   while(Ds.getNext()){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      sum=0;
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         vb.push_back(k);
         p++;
         q++;
      }
      vector<long>::iterator r=vb.begin();
      q=Bp.si.begin();
      while(q!=Bp.si.end()){
         u=(Z)d_local(*q,sum);
         u*=sgw[*r];
         pfx->write((char*)&u,sizeof(Z));
         q++;
         r++;
      }
      Bp.clear();
      vb.clear();
      mark(i++,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzNormG(DStor &Ds,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,sum;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,ssq;
   const char *pch;
   vector<long> vb;
   vector<Z> vz;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   i=0;
   while(Ds.getNext()){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      sum=0;
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         vb.push_back(k);
         p++;
         q++;
      }
      vector<long>::iterator r=vb.begin();
      q=Bp.si.begin();
      ssq=0;
      while(q!=Bp.si.end()){
         u=(Z)d_local(*q,sum);
         u*=sgw[*r];
         vz.push_back(u);
         ssq+=u*u;
         q++;
         r++;
      }
      ssq=(Z)sqrt((double)ssq);
      typename std::vector<Z>::iterator s=vz.begin();
      while(s!=vz.end()){
         u=*s/ssq;
         pfx->write((char*)&u,sizeof(Z));
         s++;
      }
      Bp.clear();
      vb.clear();
      vz.clear();
      mark(i++,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzG(Doc<Z> &SDc){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)SDc.ndoc;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   for(i=0;i<SDc.ndoc;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(i);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         u*=sgw[k];
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzG(Doc<Z> &SDc,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,sum;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)SDc.ndoc;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   for(i=0;i<SDc.ndoc;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(i);
      sum=0;
      while((pch=SDc.show(u))){
         if(strstr(pch,"!!t"))sum+=u;
      }
      n1=rnd(sum);
      k=(Y)SDc.nw;
      sm+=k;
      pfs->write((char*)&k,sizeof(Y));
      SDc.reset();
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         u=sgw[k]*d_local((int)rnd(u),n1);
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}
//new_code
template<class Y,class Z>
void XPost<Y,Z>::create_DBinzNormG(Doc<Z> &SDc,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,sum,ssq;
   const char *pch;
   vector<Y> vb;
   vector<Z> vz;
   typename std::vector<Y>::iterator r,se;
   typename std::vector<Z>::iterator s;

   SDc.gopen_map();
   ndoc=(Y)SDc.ndoc;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   for(i=0;i<SDc.ndoc;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(i);
      sum=0;
      while((pch=SDc.show(u))){
         if(strstr(pch,"!!t"))sum+=u;
      }
      n1=rnd(sum);
      k=(Y)SDc.nw;
      sm+=k;
      pfs->write((char*)&k,sizeof(Y));
      ssq=0;
      SDc.reset();
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         vb.push_back(k);
         u=sgw[k]*d_local((int)rnd(u),n1);
         vz.push_back(u);
         ssq+=u*u;
      }
      ssq=(Z)sqrt((double)ssq);
      r=vb.begin(); se=vb.end();
      s=vz.begin();
      while(r!=se){
         k=*r;
         u=*s/ssq;
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
         r++;s++;
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinz(Doc<Z> &SDc,Indx<Y> *pDnd){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)pDnd->ix;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(pDnd->idx[i]);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzG(Doc<Z> &SDc,Indx<Y> *pDnd){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)pDnd->ix;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(pDnd->idx[i]);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         u*=sgw[k];
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzNorm(Doc<Z> &SDc){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,sum;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)SDc.ndoc;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   for(i=0;i<SDc.ndoc;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(i);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      sum=0;
      while((pch=SDc.show(u))){
         sum+=u*u;
      }
      sum=(Z)sqrt((double)sum);
      SDc.reset();
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         pfw->write((char*)&k,sizeof(Y));
         u/=sum;
         pfx->write((char*)&u,sizeof(Z));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzNorm(Doc<Z> &SDc,Indx<Y> *pDnd){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,sum;
   const char *pch;

   SDc.gopen_map();
   ndoc=(Y)pDnd->ix;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      SDc.read(pDnd->idx[i]);
      k=(Y)SDc.nw;
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      sum=0;
      while((pch=SDc.show(u))){
         sum+=u*u;
      }
      sum=(Z)sqrt((double)sum);
      SDc.reset();
      while((pch=SDc.show(u))){
         k=(Y)pCh->count(pch)-1;
         pfw->write((char*)&k,sizeof(Y));
         u/=sum;
         pfx->write((char*)&u,sizeof(Z));
      }
      SDc.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinz(DStor &Ds,Index *pDnd){
   long sm=0,i,j,m,n1,n2;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read(); 

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         k=(Y)pCh->count((*p).c_str())-1;
         u=(Z)(*q);
         pfw->write((char*)&k,sizeof(Y));
         pfx->write((char*)&u,sizeof(Z));
         p++;
         q++;
      }
      Bp.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinz(DStor &Ds,Index *pDnd,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,sum;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      sum=0;
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         p++;
         q++;
      }
      q=Bp.si.begin();
      while(q!=Bp.si.end()){
         u=(Z)d_local(*q,sum);
         pfx->write((char*)&u,sizeof(Z));
         q++;
      }
      Bp.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzNorm(DStor &Ds,Index *pDnd,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,sum;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,ssq;
   const char *pch;
   vector<Z> vz;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      sum=0;
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         p++;
         q++;
      }
      q=Bp.si.begin();
      ssq=0;
      while(q!=Bp.si.end()){
         u=(Z)d_local(*q,sum);
         vz.push_back(u);
         ssq+=u*u;
         q++;
      }
      ssq=(Z)sqrt((double)ssq);
      typename std::vector<Z>::iterator s=vz.begin();
      while(s!=vz.end()){
         u=*s/ssq;
         pfx->write((char*)&u,sizeof(Z));
         s++;
      }
      Bp.clear();
      vz.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzG(DStor &Ds,Index *pDnd,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,sum;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u;
   const char *pch;
   vector<long> vb;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      sum=0;
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         vb.push_back(k);
         p++;
         q++;
      }
      vector<long>::iterator r=vb.begin();
      q=Bp.si.begin();
      while(q!=Bp.si.end()){
         u=(Z)d_local(*q,sum);
         u*=sgw[*r];
         pfx->write((char*)&u,sizeof(Z));
         q++;
         r++;
      }
      Bp.clear();
      vb.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::create_DBinzNormG(DStor &Ds,Index *pDnd,float (*d_local)(int,long)){
   long sm=0,i,j,m,n1,n2,sum;
   long *sn,*ss,sp,*px=NULL,*adr;
   Y k;
   Z u,ssq;
   const char *pch;
   vector<long> vb;
   vector<Z> vz;
   DSpan Dn;
   BagPac<Z> Bp;
   Ds.gopen_read();

   get_Nnum("z",i,j);
   ndoc=(Y)j;
   nwrd=(Y)i;
   gopen_hash();
   ofstream *pfw=get_Ostr("wd",ios::out);
   ofstream *pfa=get_Ostr("ad",ios::out);
   ofstream *pfs=get_Ostr("sd",ios::out);
   ofstream *pfx=get_Ostr("xd",ios::out);

   sgw=(Z*)get_Mmap("swx");
   sm=0;
   for(i=0;i<pDnd->ix;i++){
      pfa->write((char*)&sm,sizeof(long));
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      k=(Y)Bp.ss.size();
      pfs->write((char*)&k,sizeof(Y));
      sm+=k;
      vector<string>::iterator p=Bp.ss.begin();
      typename vector<Z>::iterator q=Bp.si.begin();
      sum=0;
      while(p!=Bp.ss.end()){
         pch=(*p).c_str();
         k=(Y)pCh->count(pch)-1;
         if(strstr(pch,"!!t"))sum+=*q;
         pfw->write((char*)&k,sizeof(Y));
         vb.push_back(k);
         p++;
         q++;
      }
      vector<long>::iterator r=vb.begin();
      q=Bp.si.begin();
      ssq=0;
      while(q!=Bp.si.end()){
         u=(Z)d_local(*q,sum);
         u*=sgw[*r];
         vz.push_back(u);
         ssq+=u*u;
         q++;
         r++;
      }
      ssq=(Z)sqrt((double)ssq);
      typename std::vector<Z>::iterator s=vz.begin();
      while(s!=vz.end()){
         u=*s/ssq;
         pfx->write((char*)&u,sizeof(Z));
         s++;
      }
      Bp.clear();
      vb.clear();
      vz.clear();
      mark(i,10000,"docs");
   }
   dst_Ostr(pfw);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
   dst_Ostr(pfx);
}

template<class Y,class Z>
void XPost<Y,Z>::gopen_db_map(void){
   long i,j;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   ndoc=(Y)j;
   if(!open1){
      term=(Y*)get_Mmap("wd");
      if(Exists("xd"))dw=(Z*)get_Mmap("xd");
      else dw=NULL;
      addrd=(long*)get_Mmap("ad");
      nsiz=(Y*)get_Mmap("sd");
      pTrm=new Indx<Y>;
      open1=1;
   }
}

template<class Y,class Z>
void XPost<Y,Z>::gopen_db(void){
   long i,j;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   ndoc=(Y)j;
   if(!open1){
      term=(Y*)get_Read("wd");
      if(Exists("xd"))dw=(Z*)get_Read("xd");
      else dw=NULL;
      addrd=(long*)get_Read("ad");
      nsiz=(Y*)get_Read("sd");
      pTrm=new Indx<Y>;
      open1=1;
      inMemDocs = true;
cout << "Doc Read Data Read In!" << endl;
   }
}

template<class Y,class Z>
Indx<Y> *XPost<Y,Z>::readp_db(Y n){
   pTrm->ix=nw=nsiz[n];
   pTrm->idx=nwd=term+addrd[n];
   return(pTrm);
}

template<class Y,class Z>
void XPost<Y,Z>::set_Indx(void){
   pTrm->ix=nw;
   pTrm->idx=nwd;
}

template<class Y,class Z>
Z *XPost<Y,Z>::readz_db(Y n){
   nw=nsiz[n];
   lwt=dw+addrd[n];
   return(lwt);
}

template<class Y,class Z>
void XPost<Y,Z>::gclose_db_map(void){
   long i;
   if(open1){
      if(!inMemDocs){
         dst_Mmap("wd",(char*&)term);
         if(Exists("xd"))dst_Mmap("xd",(char*&)dw);
         dst_Mmap("ad",(char*&)addrd);
         dst_Mmap("sd",(char*&)nsiz);
      }
      else {
         delete [] term;
         if(Exists("xd"))delete [] dw;
         delete [] addrd;
         delete [] nsiz;
      }

      pTrm->idx=NULL;
      delete pTrm;
      pTrm = nullptr;
      open1=0;
   }
}

template<class Y,class Z>
void XPost<Y,Z>::create_Postn(void){
   long sm=0,i,j,*adr;
   Y *px,k,m,n,*nc;

   get_Nnum("z",i,j);
   nwrd=(Y)i;
   ndoc=(Y)j;
   freq=(Y*)get_Mmap("f");
   sm=0;
   adr=new long[i];
   nc=new Y[i];
   for(m=0;m<i;m++){
      adr[m]=sm;
      sm+=freq[m];
      nc[m]=0;
   }
   px=new Y[sm];

   gopen_db_map();
   for(n=0;n<ndoc;n++){
      readp_db(n);
      for(m=0;m<nw;m++){
         k=nwd[m];
         px[adr[k]+nc[k]]=n;
         nc[k]++;
      }
      mark((long)n,10000,"docs");
   }
   bin_Writ("a",nwrd*sizeof(long),(char*)adr);
   bin_Writ("p",sm*sizeof(Y),(char*)px);
   delete [] px;
   delete [] adr;
   delete [] nc;
}         

template<class Y,class Z>
void XPost<Y,Z>::create_Postz(void){
   long sm=0,i,j,*adr;
   Y *px,k,m,n,*nc;
   Z *cx,u;

   get_Nnum("z",i,j);
   nwrd=(Y)i;
   ndoc=(Y)j;
   freq=(Y*)get_Mmap("f");
   sm=0;
   adr=new long[i];
   nc=new Y[i];
   for(m=0;m<i;m++){
      adr[m]=sm;
      sm+=freq[m];
      nc[m]=0;
   }
   px=new Y[sm];
   cx=new Z[sm];

   gopen_db_map();
   for(n=0;n<ndoc;n++){
      readp_db(n);
      readz_db(n);
      for(m=0;m<nw;m++){
         k=nwd[m];
         u=lwt[m];
         px[adr[k]+nc[k]]=n;
         cx[adr[k]+nc[k]]=u;
         nc[k]++;
      }
      mark((long)n,10000,"docs");
   }
   bin_Writ("a",nwrd*sizeof(long),(char*)adr);
   bin_Writ("p",sm*sizeof(Y),(char*)px);
   bin_Writ("x",sm*sizeof(Z),(char*)cx);
   delete [] px;
   delete [] cx;
   delete [] adr;
   delete [] nc;
   dst_Mmap("f",(char*&)freq);
}         

template<class Y,class Z>
void XPost<Y,Z>::gopen_map(void){
   long i,j;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   ndoc=(Y)j;
   if(!open2){
      freq=(Y*)get_Mmap("f");
      addr=(long*)get_Mmap("a");
      post=(Y*)get_Mmap("p");
      if(Exists("x"))wc=(Z*)get_Mmap("x");
      pInd=new Indx<Y>;
      open2=1;
   }
}

template<class Y,class Z>
void XPost<Y,Z>::gopen_rp(void){
   long i,j;
   get_Nnum("z",i,j);
   nwrd=(Y)i;
   ndoc=(Y)j;
   if(!open2){
      freq=(Y*)get_Read("f");
      addr=(long*)get_Read("a");
      post=(Y*)get_Read("p");
      if(Exists("x"))wc=(Z*)get_Read("x");
      pInd=new Indx<Y>;
      inMemPost = true;
      open2=1;
   }
cout << "Postings Read Data Read In!" << endl;
}

template<class Y,class Z>
Indx<Y> *XPost<Y,Z>::readp(Y n){
   pInd->ix=freq[n];
   pInd->idx=post+addr[n];
   return(pInd);
}

template<class Y,class Z>
Z *XPost<Y,Z>::readz(Y n){
   wz=wc+addr[n];
   return(wz);
}

template<class Y,class Z>
void XPost<Y,Z>::gclose_map(void){
   if(open2){
      if(!inMemPost){
         dst_Mmap("f",(char*&)freq);
         dst_Mmap("a",(char*&)addr);
         dst_Mmap("p",(char*&)post);
         if(Exists("x"))dst_Mmap("x",(char*&)wc);
      }
      else {
         delete [] freq;
         delete [] addr;
         delete [] post;
         if(Exists("x"))delete [] wc;
      }
      pInd->idx=NULL;
      delete pInd;
      pInd = nullptr;
      open2=0;
   }
}

template<class Y,class Z>
void XPost<Y,Z>::gopen_hash(void){
   if(!open3){
      this->map_down((FBase *)pSt);
      pCh=pSt->gopen_Chash();
      open3=1;
   }
}

template<class Y,class Z>
long XPost<Y,Z>::find(const char *str){
   return(pCh->count(str));
}

template<class Y,class Z>
void XPost<Y,Z>::gopen_lexos(void){
   if(!open4){
      this->map_down((FBase *)pSt);
      pLx=pSt->gopen_Lexos();
      open4=1;
   }
}

template<class Y,class Z>
char *XPost<Y,Z>::show(Y n){
   return(pLx->show_str((long)n));
}

template<class Y,class Z>
void XPost<Y,Z>::gclose_hash(void){
   if(open3){
      pSt->gclose_Chash();
      pCh = nullptr;
      open3=0;
   }
}

template<class Y,class Z>
void XPost<Y,Z>::gclose_lexos(void){
   if(open4){
      pSt->gclose_Lexos();
      pLx = nullptr;
      open4=0;
   }
}

}
#endif
