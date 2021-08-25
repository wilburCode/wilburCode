#ifndef DOC_H
#define DOC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <Btree.h>
#include <FBase.h>
#include <Word.h>
#include <Clip.h>
#include <Dmap.h>
#include <DataObj.h>
#include <DStor.h>
#include <BagPac.h>
#include <Dset.h>
#include <vector>
#include <map>
#include <string>

using namespace std;
namespace iret {

class Post; //forward declaration

template<class Z>
class Doc : public FBase {

   public:
      Doc(); //Creates the space for a document to be stored.
      Doc(const char *nm); //Creates the space for a document to be stored. 
           //Stores the string *nm in *name for use in reading or writing.
      Doc(const char *nm,const char *pnm); //Creates the space for a document to 
           //be stored. Name stored. *pnm is path file suffix
      ~Doc();

   //Reading functions.
      void gopen_map(void); //Mmaps data and sets up w file reading depending
           //Sets gct=-1.
      bool read(void); //Fills a document. To be used only after a call to 
           //gopen_map(). Reads at the current file pointer position.
           //Increments gct by 1.
      bool read(Doc<Z> &Dc); //Fills a document, Dc. To be used as previous function 
      void read(long n); //Fills a document. To be used only after a
           //call to gopen_map(). Reads in document numbered n (starts at 0). 
           //Sets gct=n.
      int  size(long n); //Reads only the number of terms in doc n
      int read(string &sx); //Reads the document into a string. To be used for sequential 
           //access to the whole docset, one document at a time. Returns # of terms if doc read
           //else zero if end of file.
      void extract(string sx,Dset &Ds); //sx is string form of Doc produced by the previous
           //function. Simply puts the text strings into Ds. Does not begin with Ds.SMclear().
      void extract(string sx,Dmap<Z> &Dm); //sx is string form of Doc produced by the previous
           //function. Puts the text strings as keys and paired numbers into Dm.
           //Does not begin with Dm.SMclear()
      void gclose_map(void); //munmaps files 

   //Writing functions.
      void gopen_write(void); //Sets up for writing
      void write(); //Writes out a document. Assumes Doc is closed 
           //Increments the global counter.
      void write(Dmap<Z> *pMp); //Writes as a document
      void write(string sx); //Writes as a document
      void gclose_write(); //Writes out document number and closes associated
           //files.

   //Filling functions.
      void open(); //Creates the map at pMp and enters the document data into it. 
      void open(Dmap<Z> &Mp); //Assumes that the current document has been cleared
           //Takes Mp as the map pointed to by pMp.
      void add(const char *str,Z n); //Adds the string into the data and adds its 
           //count n also. During this addition the map at pMp is used.
      void add(Dmap<Z> &Mp); //Adds the data from Mp to the map at pMp.
           //Mp is unchanged by this operation.
      void add(Doc<Z> *pDoc); //Adds the data at pDoc to the map at pMp.
      void add_OneDels(int n,const char *); //adds 1-grams of string if len>=n 
      void add_Pair(const char *,const char *,Z); //Adds the pair of strings, no space
      void add_Trip(const char *,const char *,const char *,Z); 
           //Adds the three strings. No space is added.
      void add_Quad(const char *,const char *,const char *,const char *,Z); 
           //Adds the four strings. No space is added.
      //add max functions
      void add_max(const char *str,Z n); //Adds the string into the data and uses max 
           //count n. During this addition the map at pMp is used.
      void add_Pair_max(const char *,const char *,Z); //Adds the pair of strings, no space
      void add_Trip_max(const char *,const char *,const char *,Z); 
           //Adds the three strings. No space is added. Uses max as count
           //This is max of current number and the number already present if present
      void add_Quad_max(const char *,const char *,const char *,const char *,Z); 
           //Adds the four strings. No space is added.
           //This is max of current number and the number already present if present
      void close(); //Transfers the data into the proper document structures 
           //and deletes *pMp. Close means close for additions to the document. 
           //It can always be opened again.
      void build(Dmap<Z> &Mp); //Skips open and close. Assumes already closed and
           //is ready to write when finished.
      void buildLCW(Clip *pCl,float (*d_local)(int,long));
           //like build, but creates local weights, designed to work with standard
           //Clip functions to create MEDLINE terms.
      void clear(); //Clears the data from a document so that it can be used to 
           //build a new document. Sets sizes to zero.
      void copy(Doc<Z> *pDoc); //Assumes a document is in memory in pDoc and copies
           //to fill the local document. 
           //Local document should be closed before this is used.

   //Access functions.
      inline void reset(){ct=0;}; //A call to the show function advances the 
            //counter in the document. This function resets.
      inline const char *show(Z &j){ //Gives the address to a term and its 
                 //local value is placed in j.
         if(ct<nw){              
            j=lcnt[ct];
            ct++;
            return(word[ct-1]);
         }
         else {
            j=0;
            return(NULL);
         }
      }

   //DStor docset creation functions
      void create_Docset(DStor &Ds);
      void create_Docset(DStor &Ds,Index *pDnd); //pDnd index into PMIDs

   //Debugging functions.
      void debug(void);
           //Like write except writes to cout.
      void debug_map(void);
           //Writes data in pMp to cout

   //Global data.
      long ndoc; //Number of documents in set.
      long *addr; //Array for addresses of documents.
   //File pointers used in writing and accessing data
      ifstream *pfid; //Points at the file stream object for .w file.
      ofstream *pfod; //Points at the file stream object for .w file.
      ofstream *pfa; //For writing document address for .a file
      long gct; //Global document counter used in writing documents.

   //Local data.
      int nw; //number of words in document
      int ct; //current word in current document
      vector<const char *> word; //Holds the words in the document.    
      vector<Z> lcnt;//Holds the local terms frequencies in the document.
      Dmap<Z> *pMp; //Used with the open, add, and close functions.
};

template<class Z>
Doc<Z>::Doc(void) : FBase("docset","null"){
   nw=0;
   open1=0;
   open2=0;
   pMp=NULL;
}

template<class Z>
Doc<Z>::Doc(const char *nam) : FBase("docset",nam){
   nw=0;
   open1=0;
   open2=0;
   pMp=NULL;
}

template<class Z>
Doc<Z>::Doc(const char *nam,const char *pnam) : FBase("docset",nam,pnam){
   nw=0;
   open1=0;
   open2=0;
   pMp=NULL;
}

template<class Z>
Doc<Z>::~Doc(){
   if(open2)close();
   for(int i=0;i<nw;i++){
      delete [] word[i];
   }
   if(pMp)delete pMp;
}

template<class Z>
void Doc<Z>::gopen_map(void){
   char cnam[max_str];
   ifstream *pfin;

   pfin=get_Istr("n");
   *pfin >> ndoc;
   dst_Istr(pfin);

   if(!open1){
      addr=(long*)get_Mmap("a");
      pfid=get_Istr("w",ios::in);

      nw=0;
      gct=0;
      open1=1;
   }
   else { //may have problem when fwd is already being read
      pfid->seekg(0, ios::beg);
      gct=0;
   }
}

template<class Z>
bool Doc<Z>::read(void){
   Z xx;
   string sx;
   char *px;
   if(gct>=ndoc)return(false);
   *pfid >> nw;
   pfid->get();
   for(int i=0;i<nw;i++){
      *pfid >> xx;
      lcnt.push_back(xx);
      pfid->get();
      getline(*pfid,sx);
      px=new char[sx.size()+1];
      strcpy(px,sx.c_str());
      word.push_back(px);
   }
   ct=0;
   gct++;
   return(true);
}

template<class Z>
bool Doc<Z>::read(Doc<Z> &Dc){
   Z xx;
   string sx;
   char *px;
   if(gct>=ndoc)return(false);
   *pfid >> Dc.nw;
   pfid->get();
   for(int i=0;i<Dc.nw;i++){
      *pfid >> xx;
      Dc.lcnt.push_back(xx);
      pfid->get();
      getline(*pfid,sx);
      px=new char[sx.size()+1];
      strcpy(px,sx.c_str());
      Dc.word.push_back(px);
   }
   Dc.ct=0;
   gct++;
   return(true);
}

template<class Z>
void Doc<Z>::read(long n){
   pfid->seekg(*(addr+n));
   read();
}

template<class Z>
int Doc<Z>::size(long n){
   int u;
   pfid->seekg(*(addr+n));
   *pfid >> u;
   return(u);
}

template<class Z>
void Doc<Z>::gclose_map(void){
   if(open1){
      dst_Mmap("a",(char*&)addr);
      dst_Istr(pfid);
      open1=0;
   }
}

template<class Z>
void Doc<Z>::gopen_write(void){
   if(!open1){
      pfod=get_Ostr("w",ios::out);
      pfa= get_Ostr("a",ios::out);
      gct=0;
      open1=1;
      nw=ct=0;
      pMp=new Dmap<Z>;
   }
}

template<class Z>
void Doc<Z>::write(void){
   long u=pfod->tellp();
   pfa->write((char*)&u,sizeof(long));

   *pfod << nw << endl;
   for(ct=0;ct<nw;ct++)*pfod << lcnt[ct] << "\t" << word[ct] << endl;
   *pfod << endl;

   gct++;
}

template<class Z>
void Doc<Z>::write(Dmap<Z> *pMp){
   long u=pfod->tellp();
   pfa->write((char*)&u,sizeof(long));

   *pfod << pMp->size() << endl;
   pMp->Set();
   while(pMp->qs!=pMp->qz){
      *pfod << pMp->qs->second << "\t" << pMp->qs->first << endl;
      pMp->qs++;
   }
   *pfod << endl;

   gct++;
}

template<class Z>
void Doc<Z>::write(string sx){
   long u=pfod->tellp();
   pfa->write((char*)&u,sizeof(long));

   *pfod << sx << endl;
   gct++;
}

template<class Z>
void Doc<Z>::gclose_write(void){
   if(open1){
      ndoc=gct;
      put_Nnum("n",ndoc);

      dst_Ostr(pfod);
      pfod=NULL;
      dst_Ostr(pfa);
      pfa=NULL;
      open1=0;
   }
}

template<class Z>
void Doc<Z>::open(void){
   int i;
   if(pMp==NULL)pMp=new Dmap<Z>;
   if(!open2){
      for(i=0;i<nw;i++){
         pMp->insert(make_pair(word[i],lcnt[i]));
      }
      nw=ct=0;
      open2=1;
   }
}

template<class Z>
void Doc<Z>::open(Dmap<Z> &Mp){
   if(open2){
      close();
      clear();
   }
   if(pMp)delete pMp;
   pMp=&Mp;
   open2=1;
}

template<class Z>
void Doc<Z>::add(const char *str,Z n){
   pMp->add_count(str,n);
}

template<class Z>
void Doc<Z>::add_max(const char *str,Z n){
   pMp->max_count(str,n);
}

template<class Z>
void Doc<Z>::add(Dmap<Z> &Mp){
   Mp.Set();
   while(Mp.qs!=Mp.qz){
      pMp->add_count(Mp.qs->first,Mp.qs->second);
      Mp.qs++;
   }
}

template<class Z>
void Doc<Z>::add(Doc<Z> *pDoc){
   long i;
   for(i=0;i<pDoc->nw;i++){
      pMp->add_count(pDoc->word[i],pDoc->lcnt[i]);
   }
}

template<class Z>
void Doc<Z>::add_OneDels(int n,const char *str){
   int lxn,i,j,k;
   char cnam[1000];
   lxn=strlen(str);
   if(lxn<n)return;
   lxn=(lxn<1000)?lxn:1000;
   for(i=0;i<lxn;i++){
      strncpy(cnam,str,i);
      strcpy(cnam+i,str+i+1);
      add(cnam,1.0);
   }
}

template<class Z>
void Doc<Z>::add_Pair(const char *str,const char *ptr,Z r){
   string sx;
   sx=str;
   sx+=ptr;
   add(sx.c_str(),r);
}

template<class Z>
void Doc<Z>::add_Pair_max(const char *str,const char *ptr,Z r){
   string sx;
   sx=str;
   sx+=ptr;
   add_max(sx.c_str(),r);
}

template<class Z>
void Doc<Z>::add_Trip(const char *str,const char *ptr,const char *utr,Z r){
   string sx;
   sx=str;
   sx+=ptr;
   sx+=utr;
   add(sx.c_str(),r);
}

template<class Z>
void Doc<Z>::add_Trip_max(const char *str,const char *ptr,const char *utr,Z r){
   string sx;
   sx=str;
   sx+=ptr;
   sx+=utr;
   add_max(sx.c_str(),r);
}

template<class Z>
void Doc<Z>::add_Quad(const char *str,const char *ptr,const char *utr,const char *ztr,Z r){
   string sx;
   sx=str;
   sx+=ptr;
   sx+=utr;
   sx+=ztr;
   add(sx.c_str(),r);
}

template<class Z>
void Doc<Z>::add_Quad_max(const char *str,const char *ptr,const char *utr,const char *ztr,Z r){
   string sx;
   sx=str;
   sx+=ptr;
   sx+=utr;
   sx+=ztr;
   add_max(sx.c_str(),r);
}

template<class Z>
void Doc<Z>::close(void){
   int lxn;
   char *pch;
   if(open2){
      nw=0;
      pMp->Set();
      while(pMp->qs!=pMp->qz){
         lcnt.push_back(pMp->qs->second);
         word.push_back(pMp->qs->first);
         nw++;pMp->qs++;
      }
      pMp->clear();
      open2=0;
   }
}

template<class Z>
void Doc<Z>::build(Dmap<Z> &Mp){
   int lxn;
   char *pch;
   if(open2){
      close();
      clear();
      open2=0;
   }
   else clear();

   nw=0;
   Mp.Set();
   while(Mp.qs!=Mp.qz){
      lcnt.push_back(Mp.qs->second);
      word.push_back(Mp.qs->first);
      nw++;Mp.qs++;
   }
   Mp.clear();
}

template<class Z>
void Doc<Z>::buildLCW(Clip *pCl,float (*d_local)(int,long)){
   int lxn;
   long len=pCl->total;
   char *pch;
   if(open2){
      close();
      clear();
      open2=0;
   }
   else clear();

   nw=0;
   pCl->node_first();
   while(pCl->node_next()){
      lcnt.push_back(d_local((int)pCl->count(),len));
      word.push_back(pCl->show_str());
      nw++;
   }
}

template<class Z>
void Doc<Z>::clear(void){
   for(int i=0;i<nw;i++){
      delete [] word[i];
   }
   lcnt.clear();
   word.clear();
   ct=nw=0;
}

template<class Z>
void Doc<Z>::copy(Doc<Z> *pDoc){
   int i;
   if(open2){
      close();
      clear();
      open2=0;
   }
   else clear();
   nw=pDoc->nw;
   for(i=0;i<nw;i++){
      word.push_back(pDoc->word[i]);
      lcnt.push_back(pDoc->lcnt[i]);
   }
   pDoc->lcnt.clear();
   pDoc->word.clear();
   pDoc->nw=0;
   pDoc->ct=0;
   ct=0;
}

template<class Z>
void Doc<Z>::debug(void){
   cout << nw << endl;
   for(ct=0;ct<nw;ct++)cout << lcnt[ct] << "\t" << word[ct] << endl;
}

template<class Z>
void Doc<Z>::debug_map(void){
   if(open2){
      cout << pMp->size() << endl;
      pMp->Set();
      while(pMp->qs!=pMp->qz){
         cout << pMp->qs->second << "\t" << pMp->qs->first << endl;
         pMp->qs++;
      }
   }
}

template<class Z>
void Doc<Z>::create_Docset(DStor &Ds){
   long i;
   vector<string>::iterator p;
   typename vector<Z>::iterator q;
   DSpan Dn;
   BagPac<Z> Bp;
   gopen_write();
   Ds.gopen_read();
   i=0;
   while(Ds.getNext()){
      Ds.read(Dn);
      Bp.unpack(Dn);
      Dmap<Z> Mp;
      p=Bp.ss.begin();
      q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         Mp.add_count((*p).c_str(),(Z)(*q));
         p++;
         q++;
      }
      Bp.clear();
      build(Mp);
      write();
      clear();
      mark((long)(++i),10000,"documents in");
   }
   gclose_write();
}

template<class Z>
void Doc<Z>::create_Docset(DStor &Ds,Index *pDnd){
   long i;
   vector<string>::iterator p;
   typename vector<Z>::iterator q;
   DSpan Dn;
   BagPac<Z> Bp;
   gopen_write();
   Ds.gopen_read();
   for(i=0;i<pDnd->ix;i++){
      Ds.read(pDnd->idx[i],Dn);
      Bp.unpack(Dn);
      Dmap<Z> Mp;
      p=Bp.ss.begin();
      q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         Mp.add_count((*p).c_str(),(Z)(*q));
         p++;
         q++;
      }
      Bp.clear();
      build(Mp);
      write();
      clear();
      mark((long)i,10000,"documents in");
   }  
   gclose_write();
}

template<class Z>
int Doc<Z>::read(string &sx){
   string sz;
   int nx=0;

   sx.clear();
   if(getline(*pfid,sz)){
      sx.append(sz);
      sx.push_back('\n');
      while(sz.size()){
         getline(*pfid,sz);
         sx.append(sz);
         sx.push_back('\n');
         nx++;
      }
      return(nx);
   }    
   else return(0);
}

template<class Z>
void Doc<Z>::extract(string sx,Dset &Ds){
   int i,j;
   string sz;
   istringstream istm(sx);
   istm >> i;
   for(j=0;j<i;j++){
      getline(istm,sz,'\t');
      getline(istm,sz);
      Ds.add_key(sz.c_str());
   }
}

template<class Z>
void Doc<Z>::extract(string sx,Dmap<Z> &Dm){
   int i,j;
   Z zz;
   string sz;
   istringstream istm(sx);
   istm >> i;
   for(j=0;j<i;j++){
      istm >> zz;
      istm.get();
      getline(istm,sz);
      Dm.add_count(sz.c_str(),zz);
   }
}

}
#endif
