#include "Sav.h"

namespace iret {

Sav::Sav(const char *skey,const char *pnam) : FBase("sav",skey,pnam){
   readf();
} 

Sav::~Sav(){
   delete val;
}  
 
void Sav::readf(void){
   char anam[10000];
   int len;
   if(Exists("v")){
      ifstream *pfin=get_Istr("v");
      pfin->getline(anam,10000);
      val=new string(anam);
      dst_Istr(pfin);
   }
}

void Sav::writef(char *sxx){
   int len;
   ofstream *pfout=get_Ostr("v");
   *pfout << sxx << endl;
   dst_Ostr(pfout);
   val=new string(sxx);
}

void Sav::getval(void){
   string ss="string value for key ";
   ss+=name;
   char *ptt=cstra(0,NULL,NULL,ss.c_str());
   writef(ptt);
}
} 
