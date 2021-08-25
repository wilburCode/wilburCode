#include "Perm.h"
namespace iret {

Perm::Perm(long nm){
   long i;
   num=nm;
   mm=new long[num];
   bd=new long[num];
   cd=new long[num];
   rm=new long[num];
   for(i=0;i<num;i++){
      mm[i]=i;
      bd[i]=i+1;
      cd[i]=0;
   }
   mm[0]=1;
   mm[1]=0;
   cd[1]=-1;
} 

Perm::~Perm(){
   delete [] mm;
   delete [] bd;
   delete [] cd;
   delete [] rm;
}  
 
long Perm::next(void){
   long i,j,k,x,y;

   i=1;
   while((i<num-1)&&(cd[i]==bd[i]))i++;
   while((i<num)&&(cd[i]<bd[i])){
      for(j=0;j<=i;j++)rm[j]=mm[(j+1)%(i+1)];
      for(j=0;j<=i;j++)mm[j]=rm[j];
      cd[i]++;
      if(cd[i]<bd[i])break;
      else i++;
   }
   for(j=0;j<i;j++)cd[j]=0;
   if(i==num)return(0);
   else return(1);
}

//Combinations

Comb::Comb(long nm){
   num=nm;
   mm=new long[num];
}

Comb::~Comb(){
   delete [] mm;
}

void Comb::set_mm(long t){
   long j;
   tt=t;
   for(j=0;j<num;j++){
      if(j<tt)*(mm+j)=1;
      else *(mm+j)=0;
   }
   ct=0;
}

long Comb::next(){
   int i,k=num-1,flag=1;
   if(!ct){ct++;return(1);}
   else ct++;
   while(mm[k]==0)k--;
   if(k<num-1){
      mm[k]=0;
      mm[k+1]=1;
   }
   else {
      long u=0;
      while((k>=0)&&mm[k]){k--;u++;}
      if(u==tt)flag=0;
      else {
         while(mm[k]==0)k--;
         mm[k]=0;
         mm[k+1]=1;
         for(i=k+2;i<k+2+u;i++)mm[i]=1;
         for(i=k+2+u;i<num;i++)mm[i]=0;
      }
   }
return(flag);
}

//Permutation Pairs Representation

PermPairs::PermPairs(long nm){
   num=nm;
   tnm=num*num;
   tm=new long[tnm];
}

PermPairs::~PermPairs(){
   delete [] tm;
}

long PermPairs::Copy(long *tx){
   long i;
   long j=0;
   for(i=0;i<tnm;i++)j+=tm[i]=tx[i];
   return(j);
}

long PermPairs::ConvertF(long *m){
   long i,j,k,ct=0;

   for(i=0;i<num-1;i++){
      for(j=i+1;j<num;j++){
         if(m[i]<m[j]){tm[m[i]*num+m[j]]=1;ct++;}
         else tm[m[j]*num+m[i]]=0;
      }
   }
   return(ct);
}

void PermPairs::ConvertB(long *m){
   long i,j,k;

   long *ix=new long[num];
   for(i=0;i<num;i++)ix[i]=0;
   for(i=0;i<num-1;i++){
      for(j=i+1;j<num;j++){
         if(tm[i*num+j])ix[j]++;
         else ix[i]++;
      }
   }
   for(i=0;i<num;i++)m[ix[i]]=i;
   delete [] ix;
}

long PermPairs::Consist(void){
   long i,j,k,n,flag;

   flag=1;
   for(i=0;i<num-2;i++){
      for(j=i+1;j<num-1;j++){
         for(k=j+1;k<num;k++){
            if(((n=tm[i*num+j])==tm[j*num+k])&&(n!=tm[i*num+k])){flag=0;break;}
         }
         if(!flag)break;
      }
      if(!flag)break;
   }
   return(flag);
}

} 
