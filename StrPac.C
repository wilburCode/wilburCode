#include <StrPac.h>
namespace iret {

StrPac::StrPac(void){
} 

StrPac::~StrPac(){
}  
 
void StrPac::add_str(const char *str,const char typ){
   string ptt(str);
   ss.push_back(ptt);
   sc.push_back(typ);
}

void StrPac::add_str(string &str,const char typ){
   ss.push_back(str);
   sc.push_back(typ);
}

void StrPac::pack(DSpan &Ds){
   int i=(int)ss.size();
   vector<string>::iterator p=ss.begin();
   Ds.id=id;
   Ds.date=date;
   while(p!=ss.end()){
      i+=p->size()+1;
      p++;
   }
   if(Ds.bln<i){
      if(Ds.buf)delete [] Ds.buf;
      Ds.buf=new char[i];
      Ds.bln=i;
   }
   Ds.len=i;
   p=ss.begin();
   vector<char>::iterator q=sc.begin();
   string::iterator u;
   int j=0;
   while(p!=ss.end()){
      Ds.buf[j++]=*q;
      u=p->begin();
      while(u!=p->end()){
         Ds.buf[j++]=*u;
         u++;
      }
      Ds.buf[j++]='\0';
      p++;
      q++;
   }
}

int StrPac::unpack(DSpan &Ds){
   int i=0;
   id=Ds.id;
   date=Ds.date;
   char *ptt=Ds.buf;
   string stx;
   while(i<Ds.len){
      sc.push_back(ptt[i++]);
      stx=ptt+i;
      i+=stx.size()+1;
      ss.push_back(stx);
   }
   return(ss.size());
}

void StrPac::clear(void){
   sc.clear();
   ss.clear();
}

void StrPac::show_data(void){
   cout << "id " << id << endl;
   cout << "date " << date << endl;
   vector<string>::iterator p=ss.begin();
   vector<char>::iterator q=sc.begin();
   while(p!=ss.end()){
      cout << *q << "| " << *p << endl;
      p++;
      q++;
   }
}

} 
