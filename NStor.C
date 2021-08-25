#include <vector>
#include <NStor.h>
namespace iret {

NStor::NStor(const char *nam,const char *path_nam) : FBase("nstor",nam,path_nam) {
   open1=0;
} 

NStor::~NStor(){
}  
 
void NStor::create_NStor(set<long> &stx){
   long i=0,k;
   ofstream *pfx=get_Ostr("x",ios::out);
   set<long>::iterator si;
   for(si=stx.begin();si!=stx.end();si++){
      k=*si;
      pfx->write((char*)&k,sizeof(long));
      mark(++i,100000,"numbers");
   }
   dst_Ostr(pfx);
   k=stx.size();
   put_Nnum("n",k);
}

void NStor::gopen_map(void){
   if(!open1){
      get_Nnum("n",num);
      xn=(long*)get_Mmap("x");
      open1=1;
   }
}

long NStor::find(long m){
   if(m<xn[0])return(0);
   else if(m==xn[0])return(1);

   if(m>xn[num-1])return(0);
   else if(m==xn[num-1])return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if(m==xn[i])return(i+1);
      else if(m<xn[i])y=i;
      else x=i;
   }
   return(0);
}

long NStor::lfind(long m){
   if(m<xn[0])return(-1);
   else if(m==xn[0])return(1);

   if(m>xn[num-1])return(-num-1);
   else if(m==xn[num-1])return(num);

   long i,x=0,y=num-1;
   while(y-x>1){
      i=(y+x)/2;
      if(m==xn[i])return(i+1);
      else if(m<xn[i])y=i;
      else x=i;
   }
   return(-x-2);
}

void NStor::gclose_map(void){
   if(open1){
      dst_Mmap("x",(char*&)xn);
      open1=0;
   }
}

void NStor::update_NStor(set<long> &stx){
   long i=0,j,k,m;
   gopen_map();
   vector<long> vx;
   set<long>::iterator si;
   set<long> sty;
   for(si=stx.begin();si!=stx.end();si++){
      j=lfind(*si);
      if(j<0){
         sty.insert(*si);
         vx.push_back(-j);
      }
      mark(++i,1000,"strings one");
   }
   gclose_map();
   k=vx.size();
   if(!k)return;

   ofstream *pfx=get_Ostr("x",ios::app);
   for(si=sty.begin();si!=sty.end();si++){
      j=*si;
      pfx->write((char*)&j,sizeof(long));
      num++;
   }
   dst_Ostr(pfx);
   put_Nnum("n",num);
   long *ta = (long*)get_Wmap("x");
   long *tb=new long[k];
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
   mak_Msync("x",(char *)ta);
   dst_Mmap("x",(char *&)ta);
   delete [] tb;
}   

}
