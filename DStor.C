#include <iostream>
#include <iostream>
#include <runn.h>
#include <DStor.h>
namespace iret {

void DSpan::showData(void){
   int i,j,k;
   cout << "id:   " << id << endl;
   cout << "date: " << date << endl;
   cout << "len:  " << len << endl;
   cout << "bln:  " << bln << endl;
   for(i=0;i<len;i++){
      cout << i << " " << (int)buf[i] << " (" << buf[i] << ")" << endl;
   }
}

DSlice::DSlice(const char *nam,int nm) : FBase("dstor",nam){
   nsl  =nm;
   pout_slice=NULL;
   pin_slice =NULL;
   ids  =NULL;
   dat  =NULL;
   adr  =NULL;
   sts  =NULL;
   ord  =NULL;
   siz  =0;
   open1=0; //gopen_add marker
   open2=0; //gopen_read marker
   open3=0;
} 

DSlice::~DSlice(){
   if(open1)gclose_add();
   if(open2==1)gclose_read();
   if(open3)gclose_markStatus();
}  
 
bool DSlice::gopen_add(void){
   bool bo;
   if(open1)return(true);
   pout_slice=get_Ostr(nsl,"slice",ios::app);
   bo=pout_slice->is_open();
   if(bo)open1=1;
   return(bo);
}

bool DSlice::add(DSpan &Ds){
   pout_slice->write((char*)&(Ds.id),sizeof(long));
   pout_slice->write((char*)&(Ds.date),sizeof(long));
   pout_slice->write((char*)&(Ds.len),sizeof(int));
   pout_slice->write(Ds.buf,Ds.len);
   pout_slice->write((char*)&(Ds.len),sizeof(int));
   return(pout_slice->good());
}

void DSlice::gclose_add(void){
   if(!open1)return;
   dst_Ostr(pout_slice);
   open1=0;
}

bool DSlice::createIndex(void){
   long i,j,k,ki,kd,ct=0,flag;
   int xp,yp;
   bool bo,co;

   if(open1)gclose_add();
   if(get_Fsiz(nsl,"ids")){
      bo=updateIndex();
      co=createStatus();
      return(bo&co);
   }
   pin_slice=get_Istr(nsl,"slice",ios::in);
   ofstream *pout_ids=get_Ostr(nsl,"ids",ios::out);
   ofstream *pout_date=get_Ostr(nsl,"dat",ios::out);
   ofstream *pout_addr=get_Ostr(nsl,"adr",ios::out);
   flag=1;
   while(flag){
      j=pin_slice->tellg();
      pin_slice->read((char*)&ki,sizeof(long));
      if(pin_slice->gcount()==0)break;
      pout_ids->write((char*)&ki,sizeof(long));
      pin_slice->read((char*)&kd,sizeof(long));
      pout_date->write((char*)&kd,sizeof(long));
      pin_slice->read((char*)&xp,sizeof(int));
      for(k=0;k<xp;k++)pin_slice->get();
      pin_slice->read((char*)&yp,sizeof(int));
      if(xp!=yp){
         cout << "Index Creation Failure, broken record id " << ki << " date " << kd << " slice " << nsl << endl;
         return(false);
      }
      pout_addr->write((char*)&j,sizeof(long));
      mark(++ct,10000,"create index");
   }
   put_Nnum(nsl,"siz",ct);
   dst_Istr(pin_slice);
   pout_ids->close();
   pout_ids->clear();
   pout_date->close();
   pout_date->clear();
   pout_addr->close();
   pout_addr->clear();
   if(ct)return(createStatus());
   else return(false);
}
   
bool DSlice::updateIndex(void){
   long i,j,k,ki,kd,ct,flag;
   int xp,yp;
   //Get span count
   if(!get_Fsiz(nsl,"ids")){cout << "File \"ids\" does not exist" << endl; return(false);}
   get_Nnum(nsl,"siz",ct);
   //Read through last span
   adr=(long*)get_Mmap(nsl,"adr");
   pin_slice=get_Istr(nsl,"slice",ios::in);
   pin_slice->seekg(adr[ct-1]);
   pin_slice->read((char*)&k,sizeof(long));
   if(pin_slice->gcount()==0){cout << "Error in read, no bytes" << endl; return(false);}
   pin_slice->read((char*)&k,sizeof(long));
   pin_slice->read((char*)&xp,sizeof(int));
   for(k=0;k<xp;k++)pin_slice->get();
   pin_slice->read((char*)&yp,sizeof(int));
   if(xp!=yp)return(false);
   //Open files to append to
   ofstream *pout_ids=get_Ostr(nsl,"ids",ios::app);
   ofstream *pout_date=get_Ostr(nsl,"dat",ios::app);
   ofstream *pout_addr=get_Ostr(nsl,"adr",ios::app);
   flag=1;
   while(flag){
      j=pin_slice->tellg();
      pin_slice->read((char*)&ki,sizeof(long));
      if(pin_slice->gcount()==0)break;
      pout_ids->write((char*)&ki,sizeof(long));
      pin_slice->read((char*)&kd,sizeof(long));
      pout_date->write((char*)&kd,sizeof(long));
      pin_slice->read((char*)&xp,sizeof(int));
      for(k=0;k<xp;k++)pin_slice->get();
      pin_slice->read((char*)&yp,sizeof(int));
      if(xp!=yp){
         cout << "Index Update Failure, broken record id " << ki << " date " << kd << " slice " << nsl << endl;
         return(false);
      }
      pout_addr->write((char*)&j,sizeof(long));
      mark(++ct,10000,"update index");
   }
   put_Nnum(nsl,"siz",ct);
   dst_Istr(pin_slice);
   pout_ids->close();
   pout_ids->clear();
   pout_date->close();
   pout_date->clear();
   pout_addr->close();
   pout_addr->clear();
   return(true);
}

bool DSlice::createStatus(void){
   long i,j,k,m,ct,flag;
   //Get span count
   if(!get_Fsiz(nsl,"ids")){cout << "File \"ids\" does not exist" << endl; return(false);}
   get_Nnum(nsl,"siz",ct);
   ids=(long*)get_Read(nsl,"ids");
   dat=(long*)get_Read(nsl,"dat");
   ord=new long[ct];
   for(i=0;i<ct;i++)ord[i]=i;
   hSort(ct,ids,dat,ord);
   sts=new char[ct];
   for(i=0;i<ct;i++)sts[i]='o';
   k=ids[0];
   j=dat[0];
   i=m=0;
   while(i<ct){
      while((i+1<ct)&&(ids[i+1]==k)){
         i++;  
         if(j<dat[i]){j=dat[i];m=i;}
      }    
      sts[ord[m]]='a';
      i++;
      if(i<ct){
         k=ids[i];
         j=dat[i];
         m=i;
      }
   }
   bin_Writ(nsl,"sts",ct,sts);
   return(true);
}
   
bool DSlice::checkStatus(long &cur,long &obs){
   long i,ct;
   //Get span count
   if(!get_Fsiz(nsl,"sts")){cout << "File \"sts\" does not exist" << endl; return(false);}
   get_Nnum(nsl,"siz",ct);
   sts=get_Mmap(nsl,"sts");
   cur=obs=0;
   for(i=0;i<ct;i++){
      if(sts[i]=='o')obs++;
      else cur++;
   }
   dst_Mmap(nsl,"sts",(char*&)sts);
   cout << "Slice " << nsl << " of DStor " << this->name << endl;
   cout << "Current  records = " << cur << endl;
   cout << "Obsolete records = " << obs << endl;
   cout << endl;
   return(true);
}
   
bool DSlice::gopen_read(void){
   long i,j,k,m,ct,flag;
   char ch;
   if(open2)return(true);
   //Get span count
   if(!get_Fsiz(nsl,"ids")){ cerr << "ids file " << nsl << " empty" << endl;return(false);}
   get_Nnum(nsl,"siz",ct);
   //Open index files to read
   ifstream *pin_ids=get_Istr(nsl,"ids",ios::in);
   if(!(pin_ids->is_open()))return(false);
   ifstream *pin_sts=get_Istr(nsl,"sts",ios::in);
   if(!(pin_sts->is_open()))return(false);
   ids=new long[ct];
   ord=new long[ct];
   j=0;
   for(i=0;i<ct;i++){
      pin_sts->read(&ch,sizeof(char));
      pin_ids->read((char*)&m,sizeof(long));
      if(ch=='a'){
         ids[j]=m;
         ord[j++]=i;
      }
   }
   pin_ids->close();
   pin_ids->clear();
   delete pin_ids;
   pin_sts->close();
   pin_sts->clear();
   delete pin_sts;
   csiz=j;
   hSort(csiz,ids,ord);
   adr=(long*)get_Mmap(nsl,"adr");
   dat=(long*)get_Mmap(nsl,"dat");
   pin_slice=get_Istr(nsl,"slice",ios::in);
   if(!(pin_slice->is_open()))return(false);
   open2=1;
   return(true);
}
   
bool DSlice::gopen_read(long date_thr){
   long i,j,k,m,ct,flag;
   if(open2)return(true);
   //Get span count
   if(!get_Fsiz(nsl,"ids")){ cerr << "ids file " << nsl << " empty" << endl;return(false);}
   get_Nnum(nsl,"siz",ct);
   //Open index files to read
   ifstream *pin_ids=get_Istr(nsl,"ids",ios::in);
   if(!(pin_ids->is_open()))return(false);
   ifstream *pin_dat=get_Istr(nsl,"dat",ios::in);
   if(!(pin_dat->is_open()))return(false);
   ids=new long[ct];
   dat=new long[ct];
   ord=new long[ct];
   j=0;
   for(i=0;i<ct;i++){
      pin_dat->read((char*)&k,sizeof(long));
      pin_ids->read((char*)&m,sizeof(long));
      if(k<=date_thr){
         dat[j]=k;
         ids[j]=m;
         ord[j++]=i;
      }
   }
   pin_ids->close();
   pin_ids->clear();
   delete pin_ids;
   pin_dat->close();
   pin_dat->clear();
   delete pin_dat;
   csiz=j;
   hSort(csiz,ids,dat,ord);
   //Mark just latest as active
   sts=new char[ct];
   for(i=0;i<ct;i++)sts[i]='o';
   k=ids[0];
   j=dat[0];
   i=m=0;
   while(i<csiz){
      while((i+1<ct)&&(ids[i+1]==k)){
         i++;  
         if(j<dat[i]){j=dat[i];m=i;}
      }    
      sts[ord[m]]='a';
      i++;
      if(i<csiz){
         k=ids[i];
         j=dat[i];
         m=i;
      }
   }
   //Remove 'o' spans
   j=0;
   for(i=0;i<csiz;i++){
      if(sts[ord[i]]=='a'){
         ids[j]=ids[i];
         ord[j++]=ord[i];
      }
   }
   csiz=j;
   delete [] dat;
   delete [] sts;

   adr=(long*)get_Mmap(nsl,"adr");
   dat=(long*)get_Mmap(nsl,"dat");
   pin_slice=get_Istr(nsl,"slice",ios::in);
   if(!(pin_slice->is_open()))return(false);
   open2=1;
   return(true);
}
   
   
bool DSlice::create_map(const char *prvt){
   long i,j,k,m,ct,flag;
   char ch;
   if(open2)return(true);
   //Get span count
   if(!get_Fsiz(nsl,"ids")){ cerr << "ids file " << nsl << " empty" << endl;return(false);}
   get_Nnum(nsl,"siz",ct);
   //Open index files to read
   ifstream *pin_ids=get_Istr(nsl,"ids",ios::in);
   if(!(pin_ids->is_open()))return(false);
   ifstream *pin_sts=get_Istr(nsl,"sts",ios::in);
   if(!(pin_sts->is_open()))return(false);
   ids=new long[ct];
   ord=new long[ct];
   j=0;
   for(i=0;i<ct;i++){
      pin_sts->read(&ch,sizeof(char));
      pin_ids->read((char*)&m,sizeof(long));
      if(ch=='a'){
         ids[j]=m;
         ord[j++]=i;
      }
   }
   pin_ids->close();
   pin_ids->clear();
   pin_sts->close();
   pin_sts->clear();
   csiz=j;
   hSort(csiz,ids,ord);
   FBase Fb("dstor",name,prvt);
   Fb.put_Nnum(nsl,"csiz",csiz);
   Fb.bin_Writ(nsl,"ids",sizeof(long)*csiz,(char*)ids);
   Fb.bin_Writ(nsl,"ord",sizeof(long)*csiz,(char*)ord);
   delete [] ids;
   delete [] ord;
   return(true);
}
   
   
bool DSlice::create_map(long date_thr,const char *prvt){
   long i,j,k,m,ct,flag;
   if(open2)return(true);
   //Get span count
   if(!get_Fsiz(nsl,"ids")){ cerr << "ids file " << nsl << " empty" << endl;return(false);}
   get_Nnum(nsl,"siz",ct);
   //Open index files to read
   ifstream *pin_ids=get_Istr(nsl,"ids",ios::in);
   if(!(pin_ids->is_open()))return(false);
   ifstream *pin_dat=get_Istr(nsl,"dat",ios::in);
   if(!(pin_dat->is_open()))return(false);
   ids=new long[ct];
   dat=new long[ct];
   ord=new long[ct];
   j=0;
   for(i=0;i<ct;i++){
      pin_dat->read((char*)&k,sizeof(long));
      pin_ids->read((char*)&m,sizeof(long));
      if(k<=date_thr){
         dat[j]=k;
         ids[j]=m;
         ord[j++]=i;
      }
   }
   pin_ids->close();
   pin_ids->clear();
   pin_dat->close();
   pin_dat->clear();
   csiz=j;
   hSort(csiz,ids,dat,ord);
   //Mark just latest as active
   sts=new char[ct];
   for(i=0;i<ct;i++)sts[i]='o';
   k=ids[0];
   j=dat[0];
   i=m=0;
   while(i<csiz){
      while((i+1<ct)&&(ids[i+1]==k)){
         i++;  
         if(j<dat[i]){j=dat[i];m=i;}
      }    
      sts[ord[m]]='a';
      i++;
      if(i<csiz){
         k=ids[i];
         j=dat[i];
         m=i;
      }
   }
   //Remove 'o' spans
   j=0;
   for(i=0;i<csiz;i++){
      if(sts[ord[i]]=='a'){
         ids[j]=ids[i];
         ord[j++]=ord[i];
      }
   }
   csiz=j;
   delete [] dat;
   delete [] sts;
   FBase Fb("dstor",name,prvt);
   Fb.put_Nnum(nsl,"csiz",csiz);
   Fb.bin_Writ(nsl,"ids",sizeof(long)*csiz,(char*)ids);
   Fb.bin_Writ(nsl,"ord",sizeof(long)*csiz,(char*)ord);
   delete [] ids;
   delete [] ord;

   return(true);
}
   
bool DSlice::gopen_map(const char *prvt){
   long i,j,k,m,ct,flag;
   char ch;
   if(open2)return(true);
   FBase Fb("dstor",name,prvt);
   //Get span count
   if(!Fb.get_Fsiz(nsl,"ids")){ cerr << "ids file " << nsl << " empty" << endl;return(false);}
   Fb.get_Nnum(nsl,"csiz",csiz);
   ids=(long*)Fb.get_Mmap(nsl,"ids");
   ord=(long*)Fb.get_Mmap(nsl,"ord");
   adr=(long*)get_Mmap(nsl,"adr");
   dat=(long*)get_Mmap(nsl,"dat");
   pin_slice=get_Istr(nsl,"slice",ios::in);
   if(!(pin_slice->is_open()))return(false);
   open2=2;
   return(true);
}
   
bool DSlice::read(long id,DSpan &Ds){
   long i,j,k,m,n;
   int yp;

   if(!open2)return(false);
   if(id<ids[0]){cout << "Below lower limit" << endl;return(false);}
   i=0;
   if(id>ids[csiz-1]){cout << "Above upper limit" << endl;return(false);}
   j=csiz;
   while((j-i)>1){
      k=(i+j)/2;
      if(id>=ids[k])i=k;
      else if(id<ids[k])j=k;
   }
   if(id>ids[i])return(false);
   //Read span
   pin_slice->seekg(adr[ord[i]]);
   pin_slice->read((char*)&(Ds.id),sizeof(long));
   pin_slice->read((char*)&(Ds.date),sizeof(long));
   pin_slice->read((char*)&(Ds.len),sizeof(int));
   if(Ds.len>Ds.bln){delete [] Ds.buf;Ds.buf=new char[Ds.len];Ds.bln=Ds.len;}
   pin_slice->read(Ds.buf,Ds.len);
   pin_slice->read((char*)&yp,sizeof(int));
   if(Ds.len!=yp)return(false);
   else return(true);
}

bool DSlice::exists(long id){
   long i,j,k,m,n;
   int yp;

   if(!open2)return(false);
   if(id<ids[0]){cout << "Below lower limit" << endl;return(false);}
   i=0;
   if(id>ids[csiz-1]){cout << "Above upper limit" << endl;return(false);}
   j=csiz;
   while((j-i)>1){
      k=(i+j)/2;
      if(id>=ids[k])i=k;
      else if(id<ids[k])j=k;
   }
   if(id>ids[i])return(false);
   else return(true);
}

long DSlice::date(long id){
   long i,j,k,m,n;
   int yp;

   if(!open2)return(0);
   if(id<ids[0]){cout << "Below lower limit" << endl;return(0);}
   i=0;
   if(id>ids[csiz-1]){cout << "Above upper limit" << endl;return(0);}
   j=csiz;
   while((j-i)>1){
      k=(i+j)/2;
      if(id>=ids[k])i=k;
      else if(id<ids[k])j=k;
   }
   if(id>ids[i])return(0);
   else return(dat[ord[i]]);
}

void DSlice::reset(void){
   cind=-1;
}
 
long DSlice::getNext(void){
   if(!open2)return(0);
   cind++;
   if(cind<csiz)return(ids[cind]+1);
   else return(0);
}

long DSlice::date(void){
   if(!open2)return(0);
   if(cind<csiz)return(dat[ord[cind]]);
   else return(0);
}
   
bool DSlice::read(DSpan &Ds){
   long i,j,k,m,n;
   int yp;

   if(!open2)return(false);
   //Read span
   pin_slice->seekg(adr[ord[cind]]);
   pin_slice->read((char*)&(Ds.id),sizeof(long));
   pin_slice->read((char*)&(Ds.date),sizeof(long));
   pin_slice->read((char*)&(Ds.len),sizeof(int));
   if(Ds.len>Ds.bln){delete [] Ds.buf;Ds.buf=new char[Ds.len];Ds.bln=Ds.len;}
   pin_slice->read(Ds.buf,Ds.len);
   pin_slice->read((char*)&yp,sizeof(int));
   if(Ds.len!=yp)return(false);
   else return(true);
}
   
void DSlice::gclose_read(void){
   if(!open2)return;
   dst_Mmap(nsl,"adr",(char*&)adr);
   dst_Mmap(nsl,"dat",(char*&)dat);
   delete [] ids;
   delete [] ord;
   dst_Istr(pin_slice);
   open2=0;
}

void DSlice::gclose_map(const char *prvt){
   if(!open2)return;
   dst_Mmap(nsl,"adr",(char*&)adr);
   dst_Mmap(nsl,"dat",(char*&)dat);
   FBase Fb("dstor",name,prvt);
   Fb.dst_Mmap(nsl,"ids",(char*&)ids);
   Fb.dst_Mmap(nsl,"ord",(char*&)ord);
   dst_Istr(pin_slice);
   open2=0;
}

bool DSlice::gopen_markStatus(void){
   if(open3)return(true);
   if(!get_Fsiz(nsl,"sts")){ cerr << "sts file " << nsl << " empty" << endl;return(false);}
   sts=get_Wmap(nsl,"sts");
   open3=1;
   return(true);
}
   
bool DSlice::markStatus(long id,char status){
   long i,j,k,m,n;
   if(!open3)return(false);
   if(id<ids[0]){cout << "Below lower limit" << endl;return(false);}
   i=0;
   if(id>ids[csiz-1]){cout << "Above upper limit" << endl;return(false);}
   j=csiz;
   while((j-i)>1){
      k=(i+j)/2;
      if(id>=ids[k])i=k;
      else if(id<ids[k])j=k;
   }
   if(id>ids[i])return(0);
   sts[ord[i]]=status;
   return(true);
}
        
bool DSlice::markStatus(char status){   
   if(!open3)return(false);
   sts[ord[cind]]=status;
   return(true);
}
   
void DSlice::gclose_markStatus(void){
   if(!open3)return;
   dst_Mmap(nsl,"sts",sts);
   open3=0;
}
   
bool DSlice::cleanSlice(void){
   char bnam[max_str],cnam[max_str];
   long i,ct=0;
   if(!get_Fsiz(nsl,"ids")){ cerr << "ids file " << nsl << " empty" << endl;return(false);}
   ofstream *pout_temp=get_Ostr(nsl,"temp");
   if(!pout_temp->is_open())return(false);
   ofstream *pout_ids=get_Ostr(nsl,"ids");
   if(!pout_ids->is_open())return(false);
   ofstream *pout_dat=get_Ostr(nsl,"dat");
   if(!pout_dat->is_open())return(false);
   ofstream *pout_adr=get_Ostr(nsl,"temp_adr");
   if(!pout_adr->is_open())return(false);
   
   DSpan Ds;
   reset();
   while(getNext()){
      if(!read(Ds))return(false);
      i=pout_temp->tellp();
      pout_adr->write((char*)&i,sizeof(long));
      pout_temp->write((char*)&(Ds.id),sizeof(long));
      pout_ids->write((char*)&(Ds.id),sizeof(long));
      pout_temp->write((char*)&(Ds.date),sizeof(long));
      pout_dat->write((char*)&(Ds.date),sizeof(long));
      pout_temp->write((char*)&(Ds.len),sizeof(int));
      pout_temp->write(Ds.buf,Ds.len);
      pout_temp->write((char*)&(Ds.len),sizeof(int));
      ct++;
   }
   put_Nnum(nsl,"siz",ct);
   bool bo=pout_temp->good();
   dst_Ostr(pout_temp);
   dst_Ostr(pout_ids);
   dst_Ostr(pout_dat);
   dst_Ostr(pout_adr);
   //replace "slice" file
   strcpy(bnam,"mv ");
   get_pathx(cnam,nsl,"temp");
   strcat(bnam,cnam);
   strcat(bnam," ");
   get_pathx(cnam,nsl,"slice");
   strcat(bnam,cnam);
   if(system(bnam))return(false);
   //replace "adr" file
   dst_Mmap(nsl,"adr",(char*&)adr);
   delete [] ids;
   delete [] ord;
   dst_Istr(pin_slice);
   open2=0;
   strcpy(bnam,"mv ");
   get_pathx(cnam,nsl,"temp_adr");
   strcat(bnam,cnam);
   strcat(bnam," ");
   get_pathx(cnam,nsl,"adr");
   strcat(bnam,cnam);
   if(system(bnam))return(false);
   createStatus();
   return(bo);
}
         
//DStor
DStor::DStor(const char *nam,const char *path_nam) : FBase("dstor",nam){
   long i,j,k;
   if(*path_nam!=':'){
      set_path_name(path_nam);
      if(!readlims()){cout << "Must have \"ctp\" file!" << endl;return;}
      pSl=new DSlice*[num];
      for(i=0;i<num;i++){
         pSl[i]=new DSlice(name,i);
         pSl[i]->set_path_name(path_nam);
      }
   }
   else {
      set_path_internal(path_nam+1);
      if(!readlims()){cout << "Must have \"ctp\" file!" << endl;return;}
      pSl=new DSlice*[num];
      for(i=0;i<num;i++){
         pSl[i]=new DSlice(name,i);
         pSl[i]->set_path_internal(path_nam+1);
      }
   }
   begg=0;
   endd=num;
   open1=0;
   open2=0;
   open3=0;
}

DStor::DStor(const char *nam,const char *path_nam,long bg,long ed) : FBase("dstor",nam){
   long i,j,k;
   begg=bg;
   endd=ed;
   if(*path_nam!=':'){
      set_path_name(path_nam);
      if(!readlims()){cout << "Must have \"ctp\" file!" << endl;return;}
      pSl=new DSlice*[num];
      for(i=begg;i<endd;i++){
         pSl[i]=new DSlice(name,i);
         pSl[i]->set_path_name(path_nam);
      }
   }
   else {
      set_path_internal(path_nam+1);
      if(!readlims()){cout << "Must have \"ctp\" file!" << endl;return;}
      pSl=new DSlice*[num];
      for(i=begg;i<endd;i++){
         pSl[i]=new DSlice(name,i);
         pSl[i]->set_path_internal(path_nam+1);
      }
   }
   open1=0;
   open2=0;
   open3=0;
}

DStor::~DStor(void){
   long i,j,k;
   for(i=begg;i<endd;i++){delete pSl[i];}
   delete [] pSl;
   delete pBn;
}
 
bool DStor::readlims(void){
   long i,j,k;
   ifstream *pfin=get_Istr("ctp");
   if(!pfin->is_open())return(false);
   num=0;
   while(*pfin >> i)num++;
   dst_Istr(pfin);
   pBn=new Bnum(num);
   mm=pBn->mm;
   i=0;
   pfin=get_Istr("ctp");
   while(*pfin >> k){
      mm[i]=k;
      i++;
   }
   num--;
   dst_Istr(pfin);
   return(true);
}

bool DStor::gopen_add(void){
   long i;

   if(open1)return(true);
   for(i=begg;i<endd;i++){
      if(!pSl[i]->gopen_add())return(false);
   }
   open1=1;
   return(true);
}

bool DStor::add(DSpan &Ds){
   long i,j,k;

   if(Ds.id<mm[begg])return(false);
   if(Ds.id>=mm[endd])return(false);
   i=pBn->index(Ds.id);
   return(pSl[i]->add(Ds));
}

void DStor::gclose_add(void){
   long i;

   if(!open1)return;
   for(i=begg;i<endd;i++)pSl[i]->gclose_add();
   open1=0;
}

bool DStor::createIndices(void){
   long i;
   bool bo=false;

   for(i=begg;i<endd;i++){
      if(pSl[i]->createIndex())bo=true;
   }
   return(bo);
}

void DStor::checkStatus(void){
   long i,tcur=0,tobs=0,j,k;
   for(i=begg;i<endd;i++){
      if(!pSl[i]->checkStatus(j,k))cout << "No index for slice " << i << " of DStor " << this->name << endl << endl;
      tcur+=j;
      tobs+=k;
   }
   cout << endl << "Total current  records " << tcur << endl;
   cout << "Total obsolete records " << tobs << endl;
}

void DStor::gopen_read(void){
   long i;

   if(open2){reset();return;}
   for(i=begg;i<endd;i++){
      if(!pSl[i]->gopen_read()){cerr << "Files for " << i << " failed to open for read" << endl;}
   }
   reset();
   open2=1;
}

void DStor::gopen_read(long date_thr){
   long i;

   if(open2){reset();return;}
   for(i=begg;i<endd;i++){
      if(!pSl[i]->gopen_read(date_thr)){cerr << "Files for " << i << " failed to open for read" << endl;}
   }
   reset();
   open2=1;
}
   
void DStor::create_map(const char *prvt){
   long i;

   for(i=begg;i<endd;i++){
      if(!pSl[i]->create_map(prvt)){cerr << "Files for " << i << " failed to create map" << endl;}
   }
}

void DStor::create_map(long date_thr,const char *prvt){
   long i;

   for(i=begg;i<endd;i++){
      if(!pSl[i]->create_map(date_thr,prvt)){cerr << "Files for " << i << " failed to create map" << endl;}
   }
}
   
void DStor::gopen_map(const char *prvt){
   long i;

   if(open2){reset();return;}
   for(i=begg;i<endd;i++){
      if(!pSl[i]->gopen_map(prvt)){cerr << "Maps for " << i << " failed to open for access" << endl;}
   }
   reset();
   open2=2;
}

bool DStor::read(long id,DSpan &Ds){
   long i,j,k;

   if(id<mm[begg])return(false);
   if(id>=mm[endd])return(false);
   i=pBn->index(id);
   return(pSl[i]->read(id,Ds));
}
   
bool DStor::exists(long id){
   long i,j,k;

   if(id<mm[begg])return(false);
   if(id>=mm[endd])return(false);
   i=pBn->index(id);
   return(pSl[i]->exists(id));
}

long DStor::date(long id){
   long i,j,k;

   if(id<mm[begg])return(0);
   if(id>=mm[endd])return(0);
   i=pBn->index(id);
   return(pSl[i]->date(id));
}

void DStor::reset(void){
   long i;

   for(i=begg;i<endd;i++)pSl[i]->reset();
   curf=begg;
}

long DStor::getNext(void){
   long i,j,k;

   while((curf<endd)&&(!(i=pSl[curf]->getNext())))curf++;
   if(curf<endd)return(i);
   else return(0);
}
      
long DStor::date(void){
   long i,j,k;

   if(curf<endd)return(pSl[curf]->date());
   else return(0);
}
      
bool DStor::read(DSpan &Ds){
   long i,j,k;

   return(pSl[curf]->read(Ds));
}

void DStor::gclose_read(void){
   long i;

   if(!open2)return;
   for(i=begg;i<endd;i++)pSl[i]->gclose_read();
   open2=0;
}

void DStor::gclose_map(const char *prvt){
   long i;

   if(!open2)return;
   for(i=begg;i<endd;i++)pSl[i]->gclose_map(prvt);
   open2=0;
}

void DStor::gopen_markStatus(void){
   long i;

   if(open3)return;
   for(i=begg;i<endd;i++){
      pSl[i]->gopen_markStatus();
   }
   open3=1;
}

bool DStor::markStatus(long id,char status){
   long i,j,k;

   if(id<mm[begg])return(false);
   if(id>=mm[endd])return(false);
   i=pBn->index(id);
   return(pSl[i]->markStatus(id,status));
}

bool DStor::markStatus(char status){
   long i,j,k;

   return(pSl[curf]->markStatus(status));
}

void DStor::gclose_markStatus(void){
   long i;

   if(!open3)return;
   for(i=begg;i<endd;i++)pSl[i]->gclose_markStatus();
   open3=0;
}
   
bool DStor::cleanAll(void){
   long i;

   for(i=begg;i<endd;i++){
      if(!pSl[i]->cleanSlice()){cout << "Clean failed for " << i << endl;}
   }
   open2=0;
   return(true);
}

} 
