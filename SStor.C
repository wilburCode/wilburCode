#include <vector>
#include <SStor.h>
namespace iret {

Sset::Sset(const char *nam,const char *path_nam) : FBase("sset",nam,path_nam){
   open1=0;
}

Sset::~Sset(void){
}

void Sset::gopen_write(void){
   num=0;
   pfa=get_Ostr("a");
   pfs=get_Ostr("s");
}

void Sset::add_str(const char *pch){
   long i=pfs->tellp();   
   pfa->write((const char*)&i,sizeof(long));
   *pfs << pch << ends;
   num++;
}

void Sset::gclose_write(void){
   ofstream *pfn=get_Ostr("n");
   *pfn << num << endl;
   dst_Ostr(pfn);
   dst_Ostr(pfa);
   dst_Ostr(pfs);
}

void Sset::gopen_map(void){
   if(!open1){
      ifstream *pfn=get_Istr("n");
      *pfn >> num;
      dst_Istr(pfn);

      addr=(long*)get_Mmap("a");
      str=get_Mmap("s");
      open1=1;
   }
}

char *Sset::show_str(long n){
   if(n<0)return(NULL);
   if(n>=num)return(NULL);
   return(str+addr[n]);
}

void Sset::gclose_map(void){
   if(open1){
      dst_Mmap("a",(char*&)addr);
      dst_Mmap("s",str);
      open1=0;
   }
}

//SStor updateable string store

SStor::SStor(const char *nam,const char *path_nam) : Sset(nam,path_nam) {
   change_type("sstor");
} 

SStor::~SStor(){
}  
 
void SStor::create_SStor(set<string> &stx){
   long i=0;
   gopen_write();
   set<string>::iterator si;
   for(si=stx.begin();si!=stx.end();si++){
      add_str(si->c_str());
      mark(++i,1000,"string");
   }
   gclose_write();
}

void SStor::update_SStor(set<string> &stx){
   long ct=0,i,j,k,m;
   gopen_map();
   vector<long> vx;
   vector<set<string>::iterator> vy;
   set<string>::iterator si;
   i=0;
   for(si=stx.begin();si!=stx.end();si++){
      while((i<num)&&((m=strcmp(str+addr[i],si->c_str()))<0))i++;
      if(i<num){
         if(!m)j=i;
         else j=-i-1;
      }
      else j=-num-1;
      if(j<0){
         vy.push_back(si);
         vx.push_back(-j);
      }
      mark(++ct,1000,"strings one");
   }
   gclose_map();
   k=vx.size();
   if(!k)return;
   gopen_append();
   for(i=0;i<k;i++){
      si=vy[i];
      add_str(si->c_str());
   }
   gclose_write();
   vy.clear();
   get_Nnum("n",num);
   long *ta = (long*)get_Wmap("a");
   long *tb=new long[k];
   if(tb==NULL){cout << "Out of memory!!" << endl;exit(0);}
   for(i=0;i<k;i++)tb[i]=ta[num-k+i];
   m=k-1;
   i=num-1;
   while(vx[m]==num-k+1){
      ta[i--]=tb[m--];
   } 
   j=num-k-1;
   while(j>=0){
      ta[i--]=ta[j];
      while((m>=0)&&(vx[m]==j+1))ta[i--]=tb[m--];
      j--;
   }
   mak_Msync("a",(char *)ta);
   dst_Mmap("a",(char *&)ta);
   delete [] tb;
}   

void SStor::reOrder(void){
   long i;
   char bnam[max_str],cnam[max_str];
   gopen_map();
   string szz;
   if(eflag==0){
      szz=":";
      szz+=path;
   }
   else if(eflag==2)szz=pnam;
   Sset Nw("reorder",szz.c_str());
   Nw.change_type("sstor");
   Nw.gopen_write();
   for(i=0;i<num;i++){
      Nw.add_str(show_str(i));
      mark(i,10000,"strings");
   }
   Nw.gclose_write();
   gclose_map();
   //replace old files
   strcpy(bnam,"mv ");
   Nw.get_pathx(bnam+3,"a");
   strcat(bnam," ");
   i=strlen(bnam);
   get_pathx(bnam+i,"a");
   system(bnam);
   Nw.get_pathx(bnam+3,"s");
   strcat(bnam," ");
   i=strlen(bnam);
   get_pathx(bnam+i,"s");
   system(bnam);
   strcpy(bnam,"rm -f ");
   Nw.get_pathx(bnam+6,"n");
   system(bnam);
}

long SStor::find(const char *ssr){
   int j;
   a=b=0;
   if((j=stc_my(ssr,str+addr[0]))<0)return(0);
   else if(j==0)return(1);

   if((j=stc_my(ssr,str+addr[num-1]))>0)return(0);
   else if(j==0)return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_my(ssr,str+addr[i]))==0)return(i+1);
      else if(j<0)y=i;
      else x=i;
   }
   return(0);
}

long SStor::lfind(const char *ssr){
   int j;
   a=b=0;
   if((j=stc_my(ssr,str+addr[0]))<0)return(-1);
   else if(j==0)return(1);

   if((j=stc_my(ssr,str+addr[num-1]))>0)return(-num-1);
   else if(j==0)return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_my(ssr,str+addr[i]))==0)return(i+1);
      else if(j<0)y=i;
      else x=i;
   }
   return(-x-2);
}

int SStor::stc_my(const char *ssr,const char *ptr)
   {int i=(a<b) ? a : b;
   const char *p1=ssr+i;
   const char *p2=ptr+i;
   int j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2)return(0);
   else if(*p1<*p2){
      b=i+j;
      return(-1);
   }
   else {
      a=i+j;
      return(1);
   }
}

void SStor::gopen_append(void){
   get_Nnum("n",num);
   pfa=get_Ostr("a",ios::app);
   pfs=get_Ostr("s",ios::app);
}   

//TStor for strings and data 

TStor::TStor(const char *nam,int nd,const char *path_nam) : Sset(nam,path_nam) {
   ndim=nd;
   change_type("tstor");
   open2=0;
} 

TStor::~TStor(){
}  
 
void TStor::create_TStor(map<string,long*> &mpx){
   long i=0,*mx,cn=0;
   gopen_write();
   ofstream *pfda=get_Ostr("da",ios::out);
   ofstream *pfdd=get_Ostr("dd",ios::out);
   map<string,long*>::iterator si;
   for(si=mpx.begin();si!=mpx.end();si++){
      add_str(si->first.c_str());
      mx=si->second;
      pfdd->write((char*)mx,sizeof(long)*ndim);
      delete [] mx;
      pfda->write((char*)&cn,sizeof(long));
      cn+=ndim;
      mark(++i,1000,"string");
   }
   dst_Ostr(pfdd);
   dst_Ostr(pfda);
   put_Nnum("nd",ndim);
   gclose_write();
   time_Stamp();
}

void TStor::update_TStor(map<string,long*> &mpx){
   long i,j,k,m,*mx,*mz,cn;
   long ct=0;
   gopen_map();
   dd=(long*)get_Wmap("dd");
   da=(long*)get_Mmap("da");
   vector<long> vx;
   map<string,long*>::iterator si;
   vector<map<string,long*>::iterator> vy;
   i=0;
   for(si=mpx.begin();si!=mpx.end();si++){
      while((i<num)&&((m=strcmp(str+addr[i],si->first.c_str()))<0))i++;
      if((i<num)&&(!m)){
         mx=si->second;
         mz=dd+da[i];
         for(j=0;j<ndim;j++)mz[j]+=mx[j];
         delete [] mx;
      }
      else {
         vy.push_back(si);
         vx.push_back(i+1);
      }
      mark(++ct,10000,"strings one");
   }
   mak_Msync("dd",(char*)dd);
   dst_Mmap("dd",(char*&)dd);
   dst_Mmap("da",(char*&)da);
   gclose_map();
   k=vx.size();
   if(!k)return;
   ofstream *pfda=get_Ostr("da",ios::app);
   ofstream *pfdd=get_Ostr("dd",ios::app);
   gopen_append();
   get_Nnum("nd",ndim);
   cn=num*ndim;
   for(i=0;i<k;i++){
      si=vy[i];
      add_str(si->first.c_str());
      mx=si->second;
      pfdd->write((char*)mx,sizeof(long)*ndim);
      delete [] mx;
      pfda->write((char*)&cn,sizeof(long));
      cn+=ndim;
      mark(i,1000,"string");
   }
   vy.clear();
   dst_Ostr(pfdd);
   dst_Ostr(pfda);
   gclose_write();
   get_Nnum("n",num);
   long *ta = (long*)get_Wmap("a");
   long *ma = (long*)get_Wmap("da");
   long *tb=new long[k];
   long *mb=new long[k];
   if((tb==NULL)||(mb==NULL)){cout << "Out of memory!!" << endl;exit(0);}
   for(i=0;i<k;i++){
      tb[i]=ta[num-k+i];
      mb[i]=ma[num-k+i];
   }
   m=k-1;
   i=num-1;
   while(vx[m]==num-k+1){
      ta[i]=tb[m];
      ma[i--]=mb[m--];
   } 
   j=num-k-1;
   while(j>=0){
      ta[i]=ta[j];
      ma[i--]=ma[j];
      while((m>=0)&&(vx[m]==j+1)){
         ta[i]=tb[m];
         ma[i--]=mb[m--];
      }
      j--;
   }
/*
   //Begin test
   str=get_Mmap("s");
   for(i=1;i<num;i++){
      if(strcmp(str+ta[i-1],str+ta[i])>=0){cout << "Error in order of strings at " << i << " " << str+ta[i-1] << " | " << str+ta[i] << endl;exit(0);}
      mark(i,1000000,"test");
   }
   dst_Mmap("s",str);
   //End test
*/
   mak_Msync("a",(char*)ta);
   dst_Mmap("a",(char *&)ta);
   mak_Msync("da",(char*)ma);
   dst_Mmap("da",(char *&)ma);
   delete [] tb;
   delete [] mb;
   time_Stamp();
}   

void TStor::reOrder(void){
   long i,cn,*mx;
   char bnam[max_str],cnam[max_str];
   gopen_map();
   string szz;
   if(eflag==0){
      szz=":";
      szz+=path;
   }
   else if(eflag==2)szz=pnam;
   Sset Nw("reorder",szz.c_str());
   Nw.change_type("tstor");
   Nw.gopen_write();
   get_Nnum("nd",ndim);
   ofstream *pfdd=Nw.get_Ostr("dd");
   ofstream *pfda=Nw.get_Ostr("da");
   dd=(long*)get_Mmap("dd");
   da=(long*)get_Mmap("da");
   cn=0;
   for(i=0;i<num;i++){
      Nw.add_str(show_str(i));
      mx=dd+da[i];
      pfdd->write((char*)mx,sizeof(long)*ndim);
      pfda->write((char*)&cn,sizeof(long));
      cn+=ndim;
      mark(i,1000,"strings");
   }
   Nw.gclose_write();
   gclose_map();
   dst_Ostr(pfdd);
   dst_Ostr(pfda);
   dst_Mmap("dd",(char *&)dd);
   dst_Mmap("da",(char *&)da);
   //replace old files
   strcpy(bnam,"mv ");
   Nw.get_pathx(bnam+3,"a");
   strcat(bnam," ");
   i=strlen(bnam);
   get_pathx(bnam+i,"a");
   system(bnam);
   Nw.get_pathx(bnam+3,"s");
   strcat(bnam," ");
   i=strlen(bnam);
   get_pathx(bnam+i,"s");
   system(bnam);
   Nw.get_pathx(bnam+3,"dd");
   strcat(bnam," ");
   i=strlen(bnam);
   get_pathx(bnam+i,"dd");
   system(bnam);
   Nw.get_pathx(bnam+3,"da");
   strcat(bnam," ");
   i=strlen(bnam);
   get_pathx(bnam+i,"da");
   system(bnam);
   strcpy(bnam,"rm -f ");
   Nw.get_pathx(bnam+6,"n");
   system(bnam);
   time_Stamp();
}

void TStor::gopen_TStor(void){
   if(!open2){
      gopen_map();
      dd=(long*)get_Mmap("dd");
      da=(long*)get_Mmap("da");
      get_Nnum("nd",ndim);
      open2=1;
   }
}

int TStor::read(long m,pLong &ar){
   if((m<0)||(m>=num))return(0);
   else {
      ar=dd+da[m];
      return(1);
   }
}

int TStor::read(const char *ptr,pLong &ar){
   long j=find(ptr);
   if(j){
      ar=dd+da[j-1];
      return(1);
   }
   else return(0);
}
   
void TStor::gclose_TStor(void){
   if(open2){
      gclose_map();
      dst_Mmap("dd",(char *&)dd);
      dst_Mmap("da",(char *&)da);
      open2=0;
   }
}

long TStor::find(const char *ssr){
   int j;
   a=b=0;
   if((j=stc_my(ssr,str+addr[0]))<0)return(0);
   else if(j==0)return(1);

   if((j=stc_my(ssr,str+addr[num-1]))>0)return(0);
   else if(j==0)return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_my(ssr,str+addr[i]))==0)return(i+1);
      else if(j<0)y=i;
      else x=i;
   }
   return(0);
}

long TStor::lfind(const char *ssr){
   int j;
   a=b=0;
   if((j=stc_my(ssr,str+addr[0]))<0)return(-1);
   else if(j==0)return(1);

   if((j=stc_my(ssr,str+addr[num-1]))>0)return(-num-1);
   else if(j==0)return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if((j=stc_my(ssr,str+addr[i]))==0)return(i+1);
      else if(j<0)y=i;
      else x=i;
   }
   return(-x-2);
}

int TStor::stc_my(const char *ssr,const char *ptr)
   {int i=(a<b) ? a : b;
   const char *p1=ssr+i;
   const char *p2=ptr+i;
   int j=0;
   while((*p1==*p2)&&(*p1!='\0')){
      j++;
      p1++;
      p2++;
   }
   if(*p1==*p2)return(0);
   else if(*p1<*p2){
      b=i+j;
      return(-1);
   }
   else {
      a=i+j;
      return(1);
   }
}

void TStor::gopen_append(void){
   get_Nnum("n",num);
   pfa=get_Ostr("a",ios::app);
   pfs=get_Ostr("s",ios::app);
}   

} 
