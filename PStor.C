#include <vector>
#include <Clip.h>
#include <Word.h>
#include <StrPac.h>
#include <PStor.h>
namespace iret {

PStor::PStor(const char *nam,const char *path_nam) : TStor(nam,6,path_nam) {
   change_type("tstor");
   pNs=new NStor(name,path_nam);
   open2=0;
} 

PStor::~PStor(){
}  
 
void PStor::create_Copy(const char *pnaz){
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

void PStor::create_PStor(DStor &Dst){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,bi;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   StrPac Sp;
   Word Wd,W2;
   Clip Cl;
   Cl.Set_Up_Word(Wd);
   Wd.stop=1;
   Cl.Set_Up_Word(W2);
   W2.stop=1;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long>::iterator si,sj;
   map<string,long*>::iterator sk,sm;
   string stz;

   i=0;
   while(Dst.getNext()){
      Dst.read(Ds);
      Sp.unpack(Ds);
      flag=1;
      j=0;
      while(j<Sp.sc.size()){
         if(Sp.sc[j]=='t')u=j;
         else if(Sp.sc[j]=='a'){
            flag=0;
            break;
         }
         j++;
      }
      if(flag){Sp.clear();mark(++i,100000,"docs");continue;}
      nxx.insert(Sp.id);
      map<string,long> mx;
      {  
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
        j=u;
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
     }
     for(si=mx.begin();si!=mx.end();si++){
        strcpy(cnam,si->first.c_str());
        if(pch=strchr(cnam,' ')){
           sk=mxx.find(stz.assign(cnam));
           if(sk==mxx.end()){
              ax=new long[6];
              ax[0]=0;ax[1]=0;ax[2]=0;ax[3]=0;ax[4]=0;ax[5]=0;
              pair<map<string,long*>::iterator,bool> su=mxx.insert(make_pair(cnam,ax));
              if(su.second)sk=su.first;
              else {cout << "Error in insertion!" << endl;exit(1);}
           }
           j=si->second;
           pch++;
           sj=mx.find(stz.assign(pch));
           k=sj->second;
           pch=strrchr(cnam,' ');
           *pch='\0';
           sj=mx.find(stz.assign(cnam));
           m=sj->second;
           ax=sk->second;
           ax[0]+=j;
           ax[1]+=k-j;
           ax[2]+=m-j;
           ax[5]++;
        }
     }
     Sp.clear();
     mark(++i,100000,"docs");
   }
   i=0;
   for(sk=mxx.begin();sk!=mxx.end();sk++){
      strcpy(cnam,sk->first.c_str());
      j=0;
      pch=cnam;
      while(ptr=strchr(pch,' ')){
         j++;
         pch=++ptr;
      }
      if(j>1){
         ptr=strchr(cnam,' ');
         k=(sk->second)[0];
         sm=mxx.find(stz.assign(++ptr));
         (sm->second)[3]+=k;
         ptr=strrchr(cnam,' ');
         *ptr='\0';
         sm=mxx.find(stz.assign(cnam));
         (sm->second)[4]+=k;
      }
      mark(++i,100000,"terms");
   }
   pNs->create_NStor(nxx);
   create_TStor(mxx);
}

void PStor::update_PStor(DStor &Dst){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,pmid;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   StrPac Sp;
   Word Wd,W2;
   Clip Cl;
   Cl.Set_Up_Word(Wd);
   Wd.stop=1;
   Cl.Set_Up_Word(W2);
   W2.stop=1;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long>::iterator si,sj;
   map<string,long*>::iterator sk,sm;
   string stz;
   pNs->gopen_map();

   i=0;
   while(pmid=Dst.getNext()){
      if(pNs->find(pmid-1)){mark(++i,100000,"docs");continue;}
      Dst.read(Ds);
      Sp.unpack(Ds);
      flag=1;
      j=0;
      while(j<Sp.sc.size()){
         if(Sp.sc[j]=='t')u=j;
         else if(Sp.sc[j]=='a'){
            flag=0;
            break;
         }
         j++;
      }
      if(flag){Sp.clear();mark(++i,100000,"docs");continue;}
      nxx.insert(Sp.id);
      map<string,long> mx;
      {  
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
        j=u;
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
     }
     for(si=mx.begin();si!=mx.end();si++){
        strcpy(cnam,si->first.c_str());
        if(pch=strchr(cnam,' ')){
           sk=mxx.find(stz.assign(cnam));
           if(sk==mxx.end()){
              ax=new long[6];
              ax[0]=0;ax[1]=0;ax[2]=0;ax[3]=0;ax[4]=0;ax[5]=0;
              pair<map<string,long*>::iterator,bool> su=mxx.insert(make_pair(cnam,ax));
              if(su.second)sk=su.first;
              else {cout << "Error in insertion!" << endl;exit(1);}
           }
           j=si->second;
           pch++;
           sj=mx.find(stz.assign(pch));
           k=sj->second;
           pch=strrchr(cnam,' ');
           *pch='\0';
           sj=mx.find(stz.assign(cnam));
           m=sj->second;
           ax=sk->second;
           ax[0]+=j;
           ax[1]+=k-j;
           ax[2]+=m-j;
           ax[5]++;
        }
     }
     Sp.clear();
     mark(++i,100000,"docs");
   }
   pNs->gclose_map();
   i=0;
   for(sk=mxx.begin();sk!=mxx.end();sk++){
      strcpy(cnam,sk->first.c_str());
      j=0;
      pch=cnam;
      while(ptr=strchr(pch,' ')){
         j++;
         pch=++ptr;
      }
      if(j>1){
         ptr=strchr(cnam,' ');
         k=(sk->second)[0];
         sm=mxx.find(stz.assign(++ptr));
         (sm->second)[3]+=k;
         ptr=strrchr(cnam,' ');
         *ptr='\0';
         sm=mxx.find(stz.assign(cnam));
         (sm->second)[4]+=k;
      }
      mark(++i,100000,"terms");
   }
   pNs->update_NStor(nxx);
   update_TStor(mxx);
}

void PStor::create_PStor(DStor &Dst,Index *pMid){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,bi;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   StrPac Sp;
   Word Wd,W2;
   Clip Cl;
   Cl.Set_Up_Word(Wd);
   Wd.stop=1;
   Cl.Set_Up_Word(W2);
   W2.stop=1;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long>::iterator si,sj;
   map<string,long*>::iterator sk,sm;
   string stz;

   i=0;
   for(i=0;i<pMid->ix;i++){
      if(!Dst.read(pMid->idx[i],Ds)){mark(i+1,100000,"docs");continue;}
      Sp.unpack(Ds);
      flag=1;
      j=0;
      while(j<Sp.sc.size()){
         if(Sp.sc[j]=='t')u=j;
         else if(Sp.sc[j]=='a'){
            flag=0;
            break;
         }
         j++;
      }
      if(flag){Sp.clear();mark(i+1,100000,"docs");continue;}
      nxx.insert(Sp.id);
      map<string,long> mx;
      {  
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
        j=u;
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
     }
     for(si=mx.begin();si!=mx.end();si++){
        strcpy(cnam,si->first.c_str());
        if(pch=strchr(cnam,' ')){
           sk=mxx.find(stz.assign(cnam));
           if(sk==mxx.end()){
              ax=new long[6];
              ax[0]=0;ax[1]=0;ax[2]=0;ax[3]=0;ax[4]=0;ax[5]=0;
              pair<map<string,long*>::iterator,bool> su=mxx.insert(make_pair(cnam,ax));
              if(su.second)sk=su.first;
              else {cout << "Error in insertion!" << endl;exit(1);}
           }
           j=si->second;
           pch++;
           sj=mx.find(stz.assign(pch));
           k=sj->second;
           pch=strrchr(cnam,' ');
           *pch='\0';
           sj=mx.find(stz.assign(cnam));
           m=sj->second;
           ax=sk->second;
           ax[0]+=j;
           ax[1]+=k-j;
           ax[2]+=m-j;
           ax[5]++;
        }
     }
     Sp.clear();
     mark(++i,100000,"docs");
   }
   i=0;
   for(sk=mxx.begin();sk!=mxx.end();sk++){
      strcpy(cnam,sk->first.c_str());
      j=0;
      pch=cnam;
      while(ptr=strchr(pch,' ')){
         j++;
         pch=++ptr;
      }
      if(j>1){
         ptr=strchr(cnam,' ');
         k=(sk->second)[0];
         sm=mxx.find(stz.assign(++ptr));
         (sm->second)[3]+=k;
         ptr=strrchr(cnam,' ');
         *ptr='\0';
         sm=mxx.find(stz.assign(cnam));
         (sm->second)[4]+=k;
      }
      mark(i+1,100000,"terms");
   }
   pNs->create_NStor(nxx);
   create_TStor(mxx);
}

void PStor::update_PStor(DStor &Dst,Index *pMid){
   long i=0,*ax,cn=1;
   long j,k,m,n,flag,u,pmid;
   char *ptr,*pch,bnam[10000],cnam[100000],xnam[10000];
   DSpan Ds;
   StrPac Sp;
   Word Wd,W2;
   Clip Cl;
   Cl.Set_Up_Word(Wd);
   Wd.stop=1;
   Cl.Set_Up_Word(W2);
   W2.stop=1;
   Dst.gopen_read();
   Dst.reset();
   map<string,long*> mxx;
   set<long> nxx;
   map<string,long>::iterator si,sj;
   map<string,long*>::iterator sk,sm;
   string stz;
   pNs->gopen_map();

   i=0;
   for(i=0;i<pMid->ix;i++){
      if(pNs->find(pMid->idx[i])){mark(i+1,100000,"docs");continue;}
      if(!Dst.read(pMid->idx[i],Ds)){mark(i+1,100000,"docs");continue;}
      Sp.unpack(Ds);
      flag=1;
      j=0;
      while(j<Sp.sc.size()){
         if(Sp.sc[j]=='t')u=j;
         else if(Sp.sc[j]=='a'){
            flag=0;
            break;
         }
         j++;
      }
      if(flag){Sp.clear();mark(i+1,100000,"docs");continue;}
      nxx.insert(Sp.id);
      map<string,long> mx;
      {  
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
        j=u;
        {
           strcpy(cnam,Sp.ss[j].c_str());
           Wd.convert_protected(cnam,strlen(cnam));
           Wd.modify('\'');
           Wd.phrase();
           for(k=0;k<Wd.cnt;k++){
              W2.convert(Wd.list[k],strlen(Wd.list[k]));
              W2.single();
              for(m=0;m<W2.cnt;m++){
                 strcpy(bnam,W2.list[m]);
                 si=mx.find(stz.assign(bnam));
                 if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                 else (si->second)++;
                 for(n=m+1;n<W2.cnt;n++){
                    strcat(bnam," ");
                    strcat(bnam,W2.list[n]);
                    si=mx.find(stz.assign(bnam));
                    if(si==mx.end())mx.insert(make_pair(stz.assign(bnam),cn));
                    else (si->second)++;
                 }
              }
              W2.clear_list();
           }
           Wd.clear_list();
        }
     }
     for(si=mx.begin();si!=mx.end();si++){
        strcpy(cnam,si->first.c_str());
        if(pch=strchr(cnam,' ')){
           sk=mxx.find(stz.assign(cnam));
           if(sk==mxx.end()){
              ax=new long[6];
              ax[0]=0;ax[1]=0;ax[2]=0;ax[3]=0;ax[4]=0;ax[5]=0;
              pair<map<string,long*>::iterator,bool> su=mxx.insert(make_pair(cnam,ax));
              if(su.second)sk=su.first;
              else {cout << "Error in insertion!" << endl;exit(1);}
           }
           j=si->second;
           pch++;
           sj=mx.find(stz.assign(pch));
           k=sj->second;
           pch=strrchr(cnam,' ');
           *pch='\0';
           sj=mx.find(stz.assign(cnam));
           m=sj->second;
           ax=sk->second;
           ax[0]+=j;
           ax[1]+=k-j;
           ax[2]+=m-j;
           ax[5]++;
        }
     }
     Sp.clear();
     mark(i+1,100000,"docs");
   }
   pNs->gclose_map();
   i=0;
   for(sk=mxx.begin();sk!=mxx.end();sk++){
      strcpy(cnam,sk->first.c_str());
      j=0;
      pch=cnam;
      while(ptr=strchr(pch,' ')){
         j++;
         pch=++ptr;
      }
      if(j>1){
         ptr=strchr(cnam,' ');
         k=(sk->second)[0];
         sm=mxx.find(stz.assign(++ptr));
         (sm->second)[3]+=k;
         ptr=strrchr(cnam,' ');
         *ptr='\0';
         sm=mxx.find(stz.assign(cnam));
         (sm->second)[4]+=k;
      }
      mark(++i,100000,"terms");
   }
   pNs->update_NStor(nxx);
   update_TStor(mxx);
}

void PStor::gopen_AccessN(void){
   pNs->gopen_map();
}

bool PStor::find(long m){
   if(pNs->find(m))return(true);
   else return(false);
}

long PStor::total(void){
   return(pNs->num);
}

void PStor::gclose_AccessN(void){
   pNs->gclose_map();
}

} 
