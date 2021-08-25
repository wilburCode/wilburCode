#include "Dset.h"
namespace iret {

Dset::Dset() : set<const char *,SCmp<const char *> >(){
   cnt_key=0;
} 

Dset::~Dset(){
}  

void Dset::add_key(const char *str){
   int lxn;
   char *pch;
   if((rs=this->find(str))==this->end()){
      lxn=strlen(str);
      pch=new char[lxn+1];
      strcpy(pch,str);
      this->insert(pch);
      cnt_key++;
   }
}

void Dset::addp_key(const char *str){
   int lxn;
   char *pch;
   if((rs=this->find(str))==this->end()){
      this->insert(str);
      cnt_key++;
   }
}
 
int Dset::ifind(const char *str){
   rz=this->end();
   if((rs=this->find(str))!=rz){
      return(1);
   }
   else return(0);
}

int Dset::zfind(const char *str){
   if((rs=this->find(str))!=rz){
      return(1);
   }
   else return(0);
}

}
   
