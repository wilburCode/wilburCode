#include <vector>
#include <Clip.h>
#include <Word.h>
#include <BagPac.h>
#include <QStor.h>
namespace iret {

QStor::QStor(const char *nam,const char *path_nam) : TStor(nam,2,path_nam) {
   change_type("tstor");
   pNs=new NStor(name,path_nam);
   open2=0;
} 

QStor::~QStor(){
}  
 
void QStor::create_Copy(const char *pnaz){
   FBase Fb("tstor",name,pnaz);
   char anam[max_str],bnam[max_str];
  
   get_pathx(anam,"a");
   Fb.get_pathx(bnam,"a");
   string sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());
   get_pathx(anam,"da");
   Fb.get_pathx(bnam,"da");
   sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());
   get_pathx(anam,"dd");
   Fb.get_pathx(bnam,"dd");
   sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());
   get_pathx(anam,"n");
   Fb.get_pathx(bnam,"n");
   sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());
   get_pathx(anam,"nd");
   Fb.get_pathx(bnam,"nd");
   sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());
   get_pathx(anam,"s");
   Fb.get_pathx(bnam,"s");
   sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());

   string szz;
   if(eflag==0){
      szz=":";
      szz+=path;
   }
   else if(eflag==2)szz=pnam;
   FBase Fn("nstor",name,szz.c_str());
   Fb.change_type("nstor");

   Fn.get_pathx(anam,"n");
   Fb.get_pathx(bnam,"n");
   sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());
   Fn.get_pathx(anam,"x");
   Fb.get_pathx(bnam,"x");
   sz="cp ";
   sz+=anam;
   sz+=" ";
   sz+=bnam;
   cout << "Running: " << sz << endl;
   system(sz.c_str());
}

void QStor::create_QStor(DStor &Dst){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,bi;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   BagPac<int> Bp;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long*>::iterator si;
   vector<string>::iterator p;
   vector<int>::iterator q;

   i=0;
   while(Dst.getNext()){
      Dst.read(Ds);
      Bp.unpack(Ds);
      nxx.insert(Bp.id);
      p=Bp.ss.begin();
      q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         si=mxx.find(*p);
         u=*q;
         if(si==mxx.end()){
            ax=new long[2];
            ax[0]=1;
            ax[1]=u;
            mxx.insert(pair<string,long*>(*p,ax));
         }
         else {
            ax=si->second;
            ax[0]++;
            ax[1]+=u;
         }
         p++;
         q++;
      }
      Bp.clear();
      mark(++i,100000,"docs");
   }
   pNs->create_NStor(nxx);
   create_TStor(mxx);
}

void QStor::update_QStor(DStor &Dst){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,pmid;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   BagPac<int> Bp;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long*>::iterator si;
   vector<string>::iterator p;
   vector<int>::iterator q;
   pNs->gopen_map();

   i=0;
   while(pmid=Dst.getNext()){
      if(pNs->find(pmid-1)){mark(++i,100000,"docs");continue;}
      Dst.read(Ds);
      Bp.unpack(Ds);
      nxx.insert(Bp.id);
      p=Bp.ss.begin();
      q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         si=mxx.find(*p);
         u=*q;
         if(si==mxx.end()){
            ax=new long[2];
            ax[0]=1;
            ax[1]=u;
            mxx.insert(pair<string,long*>(*p,ax));
         }
         else {
            ax=si->second;
            ax[0]++;
            ax[1]+=u;
         }
         p++;
         q++;
      }
      Bp.clear();
      mark(++i,100000,"docs");
   }
   pNs->gclose_map();
   pNs->update_NStor(nxx);
   update_TStor(mxx);
}

void QStor::create_QStor(DStor &Dst,Index *pMid){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,bi;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   BagPac<int> Bp;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long*>::iterator si;
   vector<string>::iterator p;
   vector<int>::iterator q;

   i=0;
   for(i=0;i<pMid->ix;i++){
      if(!Dst.read(pMid->idx[i],Ds)){mark(i+1,100000,"docs");continue;}
      Bp.unpack(Ds);
      nxx.insert(Bp.id);
      p=Bp.ss.begin();
      q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         si=mxx.find(*p);
         u=*q;
         if(si==mxx.end()){
            ax=new long[2];
            ax[0]=1;
            ax[1]=u;
            mxx.insert(pair<string,long*>(*p,ax));
         }
         else {
            ax=si->second;
            ax[0]++;
            ax[1]+=u;
         }
         p++;
         q++;
      }
      Bp.clear();
      mark(++i,100000,"docs");
   }
   pNs->create_NStor(nxx);
   create_TStor(mxx);
}

void QStor::update_QStor(DStor &Dst,Index *pMid){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,pmid;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   BagPac<int> Bp;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long*>::iterator si;
   vector<string>::iterator p;
   vector<int>::iterator q;
   pNs->gopen_map();

   i=0;
   for(i=0;i<pMid->ix;i++){
      if(pNs->find(pMid->idx[i])){mark(i+1,100000,"docs");continue;}
      if(!Dst.read(pMid->idx[i],Ds)){mark(i+1,100000,"docs");continue;}
      Bp.unpack(Ds);
      nxx.insert(Bp.id);
      p=Bp.ss.begin();
      q=Bp.si.begin();
      while(p!=Bp.ss.end()){
         si=mxx.find(*p);
         u=*q;
         if(si==mxx.end()){
            ax=new long[2];
            ax[0]=1;
            ax[1]=u;
            mxx.insert(pair<string,long*>(*p,ax));
         }
         else {
            ax=si->second;
            ax[0]++;
            ax[1]+=u;
         }
         p++;
         q++;
      }
      Bp.clear();
      mark(i+1,100000,"docs");
   }
   pNs->gclose_map();
   pNs->update_NStor(nxx);
   update_TStor(mxx);
}

void QStor::gopen_AccessN(void){
   pNs->gopen_map();
}

bool QStor::find(long m){
   if(pNs->find(m))return(true);
   else return(false);
}

long QStor::total(void){
   return(pNs->num);
}

void QStor::gclose_AccessN(void){
   pNs->gclose_map();
}

} 
