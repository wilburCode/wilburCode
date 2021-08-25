#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <cmath>
#include <cstring>
#include <cassert>
#include "runn.h"
#include "Btree.h"
#include "EdtPrb.h"

using namespace std;
namespace iret {

EdtPrb::EdtPrb(const char *nam) : FBase("edtprb",nam){
   long i;
   con[0]=0;
   for(i=1;i<97;i++)con[i]=1;
   for(i=97;i<123;i++)con[i]=i-94;
   for(i=123;i<128;i++)con[i]=1;
   con[124]=2;

   n1=29;
   n2=29*n1;
   n3=29*n2;
   n4=29*n3;
   cflag=oflag=0;
   gdd=new char[max_str];
   bdd=new char[max_str];
   gdd[0]='|';
   bdd[0]='|';
}

EdtPrb::~EdtPrb(){
   delete [] gdd;
   delete [] bdd;
}

void EdtPrb::gopen_map(void){
   if(!oflag){
      div_del=(double*)get_Mmap("d");
      div_ins=(double*)get_Mmap("i");
      div_rep=(double*)get_Mmap("r");
      div_trn=(double*)get_Mmap("t");
      oflag=1;
   }
}   

void EdtPrb::set_mem(void){
   long i;
   bas_del=new double[n3];
   div_del=new double[n3];
   bas_ins=new double[n2];
   div_ins=new double[n3];
   bas_rep=new double[n3];
   div_rep=new double[n4];
   bas_trn=new double[n4];
   div_trn=new double[n4];
   for(i=0;i<n2;i++){
      bas_ins[i]=0;
   }
   for(i=0;i<n3;i++){
      bas_del[i]=0;
      div_del[i]=0;
      div_ins[i]=0;
      bas_rep[i]=0;
   }
   for(i=0;i<n4;i++){
      div_rep[i]=0;
      bas_trn[i]=0;
      div_trn[i]=0;
   }
}

long EdtPrb::ct_edits(const char *gd,const char *bd){
   char *xx=(char *)gd;
   char *yy=(char *)bd;
   long ct=0;

   while(*xx||*yy){
      while(*xx&&(*xx==*yy)){xx++;yy++;}
      if(*xx){
          if(*yy){
             if((*(xx+1)==*yy)&&((*xx==*(yy+1))&&(*(xx+2)==*(yy+2)))){
                xx+=2;yy+=2;ct++;
             }
             else if(*xx==*(yy+1)){yy++;ct++;}
             else if(*(xx+1)==*yy){xx++;ct++;}
             else if(*(xx+1)==*(yy+1)){xx++;yy++;ct++;}
             else return(-1);
          }
          else {
             if(!*(xx+1)){
                ct++;
                return(ct);
             }
             else return(-1);
          }
      }
      else if(*yy){
         if(!*(yy+1)){
            ct++;
            return(ct);
         }
         else return(-1);
      }
   }
   return(ct);
}

void EdtPrb::ct_locat(const char *gd,double fg){
   char *xx=(char *)gd;
   long len=strlen(xx),i;

   for(i=0;i<len-2;i++){
      bas_ins[con[(int)xx[i]]+n1*con[(int)xx[i+1]]]+=fg;
      bas_del[con[(int)xx[i+1]]+n1*con[(int)xx[i]]+n2*con[(int)xx[i+2]]]+=fg;
      bas_rep[con[(int)xx[i]]+n1*con[(int)xx[i+1]]+n2*con[(int)xx[i+2]]]+=fg;
      bas_trn[con[(int)xx[i]]+n1*con[(int)xx[i+1]]+n2*con[(int)xx[i+2]]+\
              n3*con[(int)xx[i+3]]]+=fg;
   }
   i=len-2;
   bas_ins[con[(int)xx[i]]+n1*con[(int)xx[i+1]]]+=fg;
   bas_del[con[(int)xx[i+1]]+n1*con[(int)xx[i]]+n2*con[(int)xx[i+2]]]+=fg;
   bas_rep[con[(int)xx[i]]+n1*con[(int)xx[i+1]]+n2*con[(int)xx[i+2]]]+=fg;
   i=len-1;
   bas_ins[con[(int)xx[i]]+n1*con[(int)xx[i+1]]]+=fg;
}

void EdtPrb::ct_edit2(const char *gd,const char *bd,double fb){
   char *xx=(char *)gd;
   char *yy=(char *)bd;
   long ct=0;
   while(*xx||*yy){
      while(*xx&&(*xx==*yy)){xx++;yy++;}
      if(*xx){
          if(*yy){
             if((*(xx+1)==*yy)&&((*xx==*(yy+1))&&(*(xx+2)==*(yy+2)))){
                div_trn[con[(int)*(xx-1)]+n1*con[(int)xx[0]]+n2*con[(int)xx[1]]+\
                        n3*con[(int)xx[2]]]+=fb;
                xx+=2;yy+=2;
             }
             else if(*xx==*(yy+1)){
                div_ins[con[(int)yy[0]]+n1*con[(int)*(xx-1)]+n2*con[(int)xx[0]]]+=fb;
                yy++;
             }
             else if(*(xx+1)==*yy){
                div_del[con[(int)xx[0]]+n1*con[(int)*(xx-1)]+n2*con[(int)xx[1]]]+=fb;
                xx++;
             }
             else if(*(xx+1)==*(yy+1)){
                div_rep[con[(int)yy[0]]+n1*con[(int)*(xx-1)]+n2*con[(int)xx[0]]+\
                            n3*con[(int)xx[1]]]+=fb;
                xx++;yy++;
             }
          }
          else {
             if(!*(xx+1)){
                div_del[con[(int)xx[0]]+n1*con[(int)*(xx-1)]+n2*con[(int)xx[1]]]+=fb;
                xx++;
             }
          }
      }
      else if(*yy){
         if(!*(yy+1)){
            div_ins[con[(int)yy[0]]+n1*con[(int)*(xx-1)]+n2*con[(int)xx[0]]]+=fb;
            yy++;
         }
      }
   }
}

void EdtPrb::zero_profile(void){
   long i;
   del=new double[max_str];
   ins=new double[max_str];
   rep=new double[max_str];
   trn=new double[max_str];

   for(i=0;i<max_str;i++){
      del[i]=ins[i]=rep[i]=trn[i]=0;
   }
}

void EdtPrb::edit_profile(const char *gd,const char *bd,double fb){
   char *xx=(char *)gd;
   char *yy=(char *)bd;
   long ct=0,i,j;

   i=j=0;
   while(xx[i]||yy[j]){
      while(xx[i]&&(xx[i]==yy[j])){i++;j++;}
      if(xx[i]){
          if(yy[j]){
             if(xx[i+1]==yy[j+1]){
                rep[i]+=fb;
                i++;j++;
             }
             else if((xx[i+1]==yy[j])&&((xx[i]==yy[j+1])&&(xx[i+2]==yy[j+2]))){
                trn[i]+=fb;
                i+=2;j+=2;
             }
             else if(xx[i]==yy[j+1]){
                ins[i]+=fb;
                j++;
             }
             else if(xx[i+1]==yy[j]){
                del[i]+=fb;
                i++;
             }
          }
          else {
             if(!xx[i+1]){
                del[i]+=fb;
                i++;
             }
          }
      }
      else if(yy[j]){
         if(!yy[j+1]){
            ins[i]+=fb;
            j++;
         }
      }
   }
}

void EdtPrb::see_profile(void){
   long i;
   for(i=0;i<max_str;i++){
      if(del[i])cout << "del " << i << del[i] << endl;
   }
   for(i=0;i<max_str;i++){
      if(ins[i])cout << "ins " << i << ins[i] << endl;
   }
   for(i=0;i<max_str;i++){
      if(rep[i])cout << "rep " << i << rep[i] << endl;
   }
   for(i=0;i<max_str;i++){
      if(trn[i])cout << "trn " << i << trn[i] << endl;
   }
}

void EdtPrb::convert_prb(void){
   long i,j;
   double sum1,sum2,zz;

   sum1=sum2=0;
   for(i=0;i<n3;i++)sum1+=bas_del[i];
   for(i=0;i<n3;i++)sum2+=div_del[i];
   zz=sum2/sum1;
   for(i=0;i<n3;i++){
      if(bas_del[i])div_del[i]=(div_del[i]+1.0)/(bas_del[i]+2.0);
      else div_del[i]=zz;
   }
   if(pflag)cout << "Aver prob of del: " << zz << endl;

   sum1=sum2=0;
   for(i=0;i<n2;i++)sum1+=bas_ins[i];
   for(i=0;i<n3;i++)sum2+=div_ins[i];
   zz=sum2/(sum1*n1);
   for(i=0;i<n2;i++){
      if(bas_ins[i]){
         for(j=0;j<n1;j++){
            div_ins[j+n1*i]=(div_ins[j+n1*i]+1.0)/(bas_ins[i]+2.0);
         }
      }
      else {
         for(j=0;j<n1;j++){
            div_ins[j+n1*i]=zz;
         }
      }
   }
   if(pflag)cout << "Aver prob of ins: " << zz << endl;

   sum1=sum2=0;
   for(i=0;i<n3;i++)sum1+=bas_rep[i];
   for(i=0;i<n4;i++)sum2+=div_rep[i];
   zz=sum2/(sum1*n1);
   for(i=0;i<n3;i++){
      if(bas_rep[i]){
         for(j=0;j<n1;j++){
            div_rep[j+n1*i]=(div_rep[j+n1*i]+1.0)/(bas_rep[i]+2.0);
         }
      }
      else {
         for(j=0;j<n1;j++){
            div_rep[j+n1*i]=zz;
         }
      }
   }
   if(pflag)cout << "Aver prob of rep: " << zz << endl;

   sum1=sum2=0;
   for(i=0;i<n4;i++)sum1+=bas_trn[i];
   for(i=0;i<n4;i++)sum2+=div_trn[i];
   zz=sum2/sum1;
   for(i=0;i<n4;i++){
      if(bas_trn[i])div_trn[i]=(div_trn[i]+1.0)/(bas_trn[i]+2.0);
      else div_trn[i]=zz;
   }
   if(pflag)cout << "Aver prob of trn: " << zz << endl;
}


void EdtPrb::zero_diagonals(void){
   long i,j,k;

   for(i=0;i<n1;i++){
      for(j=0;j<n1;j++){
         for(k=0;k<n1;k++){
            div_rep[i+j*n1+i*n2+k*n3]=0;
         }
      }
   }
   for(i=0;i<n1;i++){
      for(j=0;j<n1;j++){
         for(k=0;k<n1;k++){
            div_trn[j+i*n1+i*n2+k*n3]=0;
         }
      }
   }
}

void EdtPrb::write_prb(void){
   bin_Writ("d",n3*sizeof(double),(char*)div_del);
   bin_Writ("i",n3*sizeof(double),(char*)div_ins);
   bin_Writ("r",n4*sizeof(double),(char*)div_rep);
   bin_Writ("t",n4*sizeof(double),(char*)div_trn);
}
   
void EdtPrb::debug_prb(const char *stt){
   long i,len=strlen(stt);
   long ik,ip;
   char *pch;
   
   start:
   ik=clnga(0,0,"-ik","1 for d, 2 for i, 3 for r, 4 for t");
   switch(ik){
      case 1: cout << "Probs for deletes:" << endl;
         for(i=1;i<len;i++){
            cout << i << " " << stt[i] << " " << div_del[con[(int)stt[i]]+\
                 n1*con[(int)stt[i-1]]+n2*con[(int)stt[i+1]]] << endl;
         }
         break;
      case 2: pch=cstra(0,0,"-pch","letter to insert");
         cout << "Probs for inserting '" << pch[0] << "' after:" << endl;
         for(i=0;i<len;i++){
            cout << i << " " << stt[i] << " " << div_ins[con[(int)pch[0]]+n1*con[(int)stt[i]]+\
                 n2*con[(int)stt[i+1]]] << endl;
         }
         break;
      case 3: pch=cstra(0,0,"-pch","letter to substitute(replace with)");
         cout << "Probs for replacing with '" << pch[0]<< "' at:" << endl;
         for(i=1;i<len;i++){
            cout << i << " " << stt[i] << " " << div_rep[con[(int)pch[0]]+n1*con[(int)stt[i-1]]+\
                 n2*con[(int)stt[i]]+n3*con[(int)stt[i+1]]] << endl;
         }
         break;
      case 4: cout << "Probs for switching pair starting at:" << endl;
         for(i=1;i<len-1;i++){
            cout << i << " " << stt[i] << " " << div_trn[con[(int)stt[i-1]]+\
                 n1*con[(int)stt[i]]+n2*con[(int)stt[i+1]]+n3*con[(int)stt[i+2]]] << endl;
         }
   }
   goto start;
}

char *EdtPrb::sim_edit(const char *str){
   char *gnam,*bnam,c;
   long freq,lng,i,j,flag,seed;
   double de[200],in[200],rp[200],tr[200];
   double sum,sum2,rx,rx2,zz;

   lng=strlen(str)+1;
   gnam=new char[lng+1];
   gnam[0]='|';
   strcpy(&gnam[1],str);

   bnam=new char[lng+2];
   for(i=1;i<lng;i++){
      de[i]=0;
      in[i]=0;
      rp[i]=0;
      tr[i]=0;
   }
   sum=0;
   for(i=1;i<lng;i++){
      de[i]=div_del[con[(int)gnam[i-1]]+n1*con[(int)gnam[i]]+\
              n2*con[(int)gnam[i+1]]];
      sum+=de[i];
   }
   for(i=1;i<=lng;i++){
      for(j=3;j<n1;j++)in[i]+=div_ins[j+n1*con[(int)gnam[i-1]]+\
                              +n2*con[(int)gnam[i]]];
      sum+=in[i];
   }
   for(i=1;i<lng;i++){
      for(j=3;j<n1;j++)rp[i]+=div_rep[j+n1*con[(int)gnam[i-1]]+\
              n2*con[(int)gnam[i]]+n3*con[(int)gnam[i+1]]];
      sum+=rp[i];
   }
   for(i=1;i<lng-1;i++){
       tr[i]=div_trn[con[(int)gnam[i-1]]+n1*con[(int)gnam[i]]+\
              n2*con[(int)gnam[i+1]]+n3*con[(int)gnam[i+2]]];
      sum+=tr[i];
   }
 
   rx=drand48();
   rx=sum*rx;
   sum=0;
   i=1;
   flag=1;
   while((i<lng)&&(sum+de[i]<rx)){sum+=de[i];i++;}
   if(i<lng){
      c=gnam[i];
      gnam[i]='\0';
      strcpy(bnam,&gnam[1]);
      strcat(bnam,&gnam[i+1]);
      flag=0;
   }
   if(flag){
      i=1;
      while((i<lng+1)&&(sum+in[i]<rx)){sum+=in[i];i++;}
      if(i<lng+1){
         sum2=0;
         rx2=in[i]*drand48();
         j=3;
         while((j<n1)&&(sum2+(zz=div_ins[j+n1*con[(int)gnam[i-1]]+\
                                        n2*con[(int)gnam[i]]])<rx2)){
            sum2+=zz;j++;}
         if(j==n1)j--;
         c=gnam[i];
         gnam[i]='\0';
         strcpy(bnam,&gnam[1]);
         bnam[i-1]=j+94;
         bnam[i]='\0';
         gnam[i]=c;
         strcat(bnam,&gnam[i]);
         flag=0;
      }
   }
   if(flag){
      i=1;
      while((i<lng)&&(sum+rp[i]<rx)){sum+=rp[i];i++;}
      if(i<lng){
         sum2=0;
         rx2=rp[i]*drand48();
         j=3;
         while((j<n1)&&(sum2+(zz=div_rep[j+n1*con[(int)gnam[i-1]]+\
                                        n2*con[(int)gnam[i]]+n3*\
                                        con[(int)gnam[i+1]]])<rx2)){
            sum2+=zz;j++;}
         if(j==n1)j--;
         strcpy(bnam,&gnam[1]);
         bnam[i-1]=j+94;
         flag=0;
      }
   }
   if(flag){
      i=1;
      while((i<lng-1)&&(sum+tr[i]<rx)){sum+=tr[i];i++;}
      if(i==lng-1)i--;
      strcpy(bnam,&gnam[1]);
      bnam[i-1]=gnam[i+1];
      bnam[i]=gnam[i];
   }
   return(bnam);
}

char *EdtPrb::sim_2edit(const char *str){
   char *gnam,*bnam;

   gnam=sim_edit(str);
   bnam=sim_edit(gnam);
   delete [] gnam;
   return(bnam);
}

char *EdtPrb::sim_3edit(const char *str){
   char *gnam,*bnam;

   gnam=sim_2edit(str);
   bnam=sim_edit(gnam);
   delete [] gnam;
   return(bnam);
}

double EdtPrb::Probc(long lsn,const char *str,long lbn,const char *btr){
   strcpy(&gdd[1],str);
   strcpy(&bdd[1],btr);

   return(Probx(lsn+1,gdd,lbn+1,bdd));
}

double EdtPrb::Probx(long lsn,char *str,long lbn,char *btr){
   long i,j,k,m;
   double dx=1.0,ux=1.0;

   i=j=1;

   while((i<lsn)&&(j<lbn)){
      if(str[i]==btr[j]){i++;j++;}
      else if(str[i]==btr[j+1]){
         if(str[i+1]==btr[j]){dx=dx*div_trn[con[(int)str[i-1]]+n1*con[(int)str[i]]\
            +n2*con[(int)str[i+1]]+n3*con[(int)str[i+2]]];i+=2;j+=2;}
         else {dx=dx*div_ins[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
               n2*con[(int)str[i]]];j++;}
      }
      else if(str[i+1]==btr[j]){dx=dx*div_del[con[(int)str[i]]+n1*con[(int)str[i-1]]+\
              n2*con[(int)str[i+1]]];i++;}
      else if(str[i+1]==btr[j+1]){dx=dx*div_rep[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
            n2*con[(int)str[i]]+n3*con[(int)str[i+1]]];i++;j++;}
      else if(btr[j+1]&&(str[i]==btr[j+2])){dx=dx*div_ins[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
               n2*con[(int)str[i]]];j++;}
      else if(str[i+1]&&(str[i+2]==btr[j])){dx=dx*div_del[con[(int)str[i]]+n1*con[(int)str[i-1]]+\
              n2*con[(int)str[i+1]]];i++;}
      else {dx=dx*div_rep[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
            n2*con[(int)str[i]]+n3*con[(int)str[i+1]]];i++;j++;}
   }
   if(i<lsn){
      m=lsn-1;
      while(m>=i){dx=dx*div_del[con[(int)str[m]]+n1*con[(int)str[m-1]]]\
              ;m--;}
   }
   if(j<lbn){
      m=lbn-1;
      while(m>=j){dx=dx*div_ins[con[(int)btr[m]]+n1*con[(int)btr[m-1]]]\
              ;m--;}
   }

   if((i=Space(lsn,str))&&(j=Space(lbn,btr))){
      str[i]='\0';
      btr[j]='\0';
      ux=Probs(i,str,j,btr);
      i++;j++;
      ux=ux*Probx(lsn-i,str+i,lbn-j,btr+j);
      dx=(dx>ux)?dx:ux;
      str[i-1]=' ';
      btr[j-1]=' ';
   }

   return(dx);
} 

double EdtPrb::Probs(long lsn,const char *str,long lbn,const char *btr){
   long i,j,k,m;
   char *utr,*vtr;
   double dx=1.0;

   i=j=1;

   while((i<lsn)&&(j<lbn)){
      if(str[i]==btr[j]){i++;j++;}
      else if(str[i]==btr[j+1]){
         if(str[i+1]==btr[j]){dx=dx*div_trn[con[(int)str[i-1]]+n1*con[(int)str[i]]\
            +n2*con[(int)str[i+1]]+n3*con[(int)str[i+2]]];i+=2;j+=2;}
         else {dx=dx*div_ins[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
               n2*con[(int)str[i]]];j++;}
      }
      else if(str[i+1]==btr[j]){dx=dx*div_del[con[(int)str[i]]+n1*con[(int)str[i-1]]+\
              n2*con[(int)str[i+1]]];i++;}
      else if(str[i+1]==btr[j+1]){dx=dx*div_rep[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
            n2*con[(int)str[i]]+n3*con[(int)str[i+1]]];i++;j++;}
      else if(btr[j+1]&&(str[i]==btr[j+2])){dx=dx*div_ins[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
               n2*con[(int)str[i]]];j++;}
      else if(str[i+1]&&(str[i+2]==btr[j])){dx=dx*div_del[con[(int)str[i]]+n1*con[(int)str[i-1]]+\
              n2*con[(int)str[i+1]]];i++;}
      else {dx=dx*div_rep[con[(int)btr[j]]+n1*con[(int)str[i-1]]+\
            n2*con[(int)str[i]]+n3*con[(int)str[i+1]]];i++;j++;}
   }
   if(i<lsn){
      m=lsn-1;
      while(m>=i){dx=dx*div_del[con[(int)str[m]]+n1*con[(int)str[m-1]]]\
              ;m--;}
   }
   if(j<lbn){
      m=lbn-1;
      while(m>=j){dx=dx*div_ins[con[(int)btr[m]]+n1*con[(int)btr[m-1]]]\
              ;m--;}
   }

   return(dx);
}

long EdtPrb::Space(long lsn,const char *str){
   long i=0;
   int j;

   while((j=str[i])&&(j!=32))i++;
   if(i<lsn)return(i);
   else return(0);
}

}
