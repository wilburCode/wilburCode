#ifndef APST_H
#define APST_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <DataObj.h>
#include <XPost.h>
#include <runn.h>
#include <StrPac.h>
#include <DStor.h>
#include <Split.h>
#include <Idset.h>
#include <string>

using namespace std;
namespace iret {

//template class APost 
template<class Y,class Z>
class APost : public XPost<Y,Z> {
public:
   APost(const char *nspost,const char *pnam);//nspost is name of XPost set
      //pnam is either used in "path_pnam" as place to find path
      //or if begins with ':' is followed by path itself.
   ~APost(void);
   void createAuthLists(DStor &Dst,Idset<Y> &Id);
   void gopen_AuthLists(void);
   int  findAuth(Y dn,long na); //na is number of an author name form
      //function returns i if na appears as a name form for ith author
      //else returns 0
   int  compAuths(long xof,Y dn); //Returns i if author at xof is name
      //compatible with ith author in doc # dn, else returns 0
   int  cxAuth(long xof); //returns 1 if slots 3 or 4 are nonzero (complex name)
   Z    scoAuths(Y d1,Y d2,Y u,float *weg); //Scores all compatible author pairs in the two documents
      //weg is the weight array that are not in the name space defined by u term number.
   Z    DscoAuths(Y d1,Y d2,Y u,float *weg); //Scores all compatible author pairs in the two documents
      //weg is the weight array that are not in the name space defined by u term number. For debugging
   long ReMrk(int nx,long nb,Y *mrk); //nb is number of a base name form
      //Sets mrk array to nx for all names that have nb as initial segment
   long TotalNames(void); //Finds the total of name counts in play
   void seeConsis(long nb,Y dn); //Writes out the name belonging to doc
      //dn and having nb as base form. For debugging
   
   //Data pointers
   long *anm; //Number of names for each document
   long *atd; //term numbers for name forms, in groups of 6/author
   long *ofs; //offset for each documents group of names

};

template<class Y,class Z>
APost<Y,Z>::APost(const char *namspost,const char *pnam) : XPost<Y,Z>(namspost,pnam){
}

template<class Y,class Z>
APost<Y,Z>::~APost(){
}

template<class Y,class Z>
void APost<Y,Z>::createAuthLists(DStor &Dst,Idset<Y> &Id){
  long i,j,k,ix,kx,u,u1,v,flag,tr[6],zr=0,v1,v2,tkc;
  char *ptr,*pch,*qch,ch,hnam[500];
  Split Sp,Sq;
  string s1,s2,s3;

  this->gopen_hash();

  StrPac Pk;
  DSpan Ds;
  Dst.gopen_read();
  vector<string>::iterator p;
  vector<char>::iterator q;
  Id.gopen_Idset();
  Index *pPmid=Id.totid();

  ofstream *pFs=this->get_Ostr("autAdd",ios::out);
  ofstream *pFd=this->get_Ostr("autDat",ios::out);

  for(j=0;j<pPmid->ix;j++){
     Dst.read(pPmid->idx[j],Ds);
     Pk.unpack(Ds);
     p=Pk.ss.begin();
     q=Pk.sc.begin();
     flag=1;
     while(q!=Pk.sc.end()){
        if(*q=='u'){
           if((*p).size()>400){
              strncpy(hnam,(*p).c_str(),400);
              i=399;
              while((i>0)&&(hnam[i]!=';'))i--;
              if(i>0){
                 hnam[i]='\0';
                 Sp.token_chr_lower(';',hnam);
              }
           }
           else Sp.token_chr_lower(';',(*p).c_str());
           tkc=0;
           for(i=0;i<Sp.num;i++){
              tr[0]=tr[1]=tr[2]=tr[3]=tr[4]=tr[5]=0;
              Sq.tokenS_chr(',',Sp.lst[i]);
              if(Sq.num>1){ //more than one token in the name
                 s1=Sq.lst[0];
                 s1+=" ";
                 s1+=Sq.lst[1];
                 s1+="!!a";
                 u=this->find(s1.c_str()); //The full name as given in the record
                 pch=Sq.lst[1];
                 v1=v2=0;
                 while((*pch!=' ')&&(*pch)){pch++;v1++;}
                 if(*pch==' '){ //Has at least a middle initial
                    if(v1>1){ //Has a first name given, not just first initial
                       *pch='\0';
                       s1=Sq.lst[0];
                       s1+=" ";
                       s1+=Sq.lst[1];
                       s1+="!!a";
                       u1=this->find(s1.c_str()); //First name and last name
                    }
                    pch++;
                    qch=pch;
                    while(*pch){pch++;v2++;}
                 }
                 if(v2>0){ //Has middle initial
                    if(v1>1){ //Has first name, more than just initial
                       tr[5]=u; //Full name
                       tr[3]=u1; //First name greater than just initial 
                       if(v2>1){ //Has middle name, more than just initial
                          s1=Sq.lst[0];
                          s1+=" ";
                          Sq.lst[1][1]='\0';
                          s1+=Sq.lst[1];
                          s1+=" ";
                          s1+=qch;
                          s1+="!!a";
                          tr[4]=this->find(s1.c_str()); //First initial and middle name more than initial
                       }
                    }
                    else if(v2>1)tr[4]=u; //First initial and middle name more than initial
                 }
                 else {
                    if(v1>1)tr[3]=u; //First name greater than just initial
                    else tr[1]=u; //Just first initial and last name
                 }
              }
              else if(Sq.num){
                 s1=Sq.lst[0];
                 s1+="!!a";
                 u=this->find(s1.c_str());
                 tr[0]=u; //Single last name (could be more than one token)
              }
              if(Sq.num>2){
                 s1=Sq.lst[0];
                 s1+=" ";
                 s1+=Sq.lst[2];
                 s1+="!!a";
                 u=this->find(s1.c_str());
                 if(*(Sq.lst[2]+1)!='\0'){
                    tr[2]=u; //Two initials and last name
                    *(Sq.lst[2]+1)='\0';
                    s1=Sq.lst[0];
                    s1+=" ";
                    s1+=Sq.lst[2];
                    s1+="!!a";
                    u=this->find(s1.c_str());
                    tr[1]=u; //first initial and last name
                 }
                 else tr[1]=u; //first initial and last name
              }
              if(Sq.num){
                 pFd->write((char*)tr,sizeof(long)*6);
                 tkc++;
              }
              Sq.clear();
           }
           pFs->write((char*)&tkc,sizeof(long));
           Sp.clear();
           flag=0;
        }
        p++;q++;
     }
     Pk.clear();
     if(flag)pFs->write((char*)&zr,sizeof(long));
     mark(1,j,10000,"documents in");
  }
  this->dst_Ostr(pFd);
  this->dst_Ostr(pFs);
}   

template<class Y,class Z>
void APost<Y,Z>::gopen_AuthLists(void){
   long i,j,k;
   atd=(long*)this->get_Mmap("autDat");
   anm=(long*)this->get_Mmap("autAdd");
   ofs=new long[this->ndoc+1];
   k=0;
   for(i=0;i<this->ndoc;i++){
      ofs[i]=k;
      k+=6*anm[i];
   }
   ofs[i]=k;
}

template<class Y,class Z>
int APost<Y,Z>::findAuth(Y dn,long na){
   long i,j,k,np=na+1,flag=0;
   k=ofs[dn];
   for(i=0;i<anm[dn];i++){
      j=k+6*i;
      if(atd[j]==np){flag=i+1;break;}
      if(atd[j+1]==np){flag=i+1;break;}
      if(atd[j+2]==np){flag=i+1;break;}
      if(atd[j+3]==np){flag=i+1;break;}
      if(atd[j+4]==np){flag=i+1;break;}
      if(atd[j+5]==np){flag=i+1;break;}
   }
   return(flag);
}

template<class Y,class Z>
int APost<Y,Z>::compAuths(long xof,Y dn){
   long n,m,i,j,k=ofs[dn],flag;
   for(i=0;i<anm[dn];i++){
      flag=1;
      for(j=0;j<6;j++){
         if((n=atd[xof+j])&&(m=atd[k+j])&&(n!=m)){
            flag=0;break;
         }
      }
      if(flag)return(i+1);
      k+=6;
   }
   return(0);
}

template<class Y,class Z>
int APost<Y,Z>::cxAuth(long xof){
   long n,m,i,j,flag;

   if((n=atd[xof+3])||(m=atd[xof+4]))return(1);
   else return(0);
}

template<class Y,class Z>
Z APost<Y,Z>::scoAuths(Y d1,Y d2,Y u,float *weg){
   long n,m,i,j,flag,p,q;
   long k1=ofs[d1],k2;
   float sco=0,xx,yy;

   if((anm[d1]==1)||(anm[d2]==1))return(0);
   for(i=0;i<anm[d1];i++){
      if(atd[k1+1]==u+1){k1+=6;continue;}
      xx=0;
      k2=ofs[d2];
      for(j=0;j<anm[d2];j++){
         if(atd[k1+1]==atd[k2+1]){
            yy=0;
            for(m=0;m<6;m++){
               if((p=atd[k1+m])&&(q=atd[k2+m])){
                  if(p==q){
                     if(yy<weg[p-1])yy=weg[p-1];
                  }
                  else {
                     yy=0;
                     break;
                  }
               }
            }
            if(xx<yy)xx=yy;
         }
         k2+=6;
      }
      sco+=xx;
      k1+=6;
   }
   sco/=(anm[d1]-1)*(anm[d2]-1);
   return(sco);
}

template<class Y,class Z>
Z APost<Y,Z>::DscoAuths(Y d1,Y d2,Y u,float *weg){
   long n,m,i,j,flag,p,q,r;
   long k1=ofs[d1],k2;
   float sco=0,xx,yy;

   if((anm[d1]==1)||(anm[d2]==1))return(0);
   for(i=0;i<anm[d1];i++){
      if(atd[k1+1]==u+1){k1+=6;continue;}
      xx=0;
      k2=ofs[d2];
      for(j=0;j<anm[d2];j++){
         if(atd[k1+1]==atd[k2+1]){
            yy=0;
            for(m=0;m<6;m++){
               if((atd[k1+m])&&(atd[k2+m])){
                  p=atd[k1+m];
                  q=atd[k2+m];
                  if(p==q){
                     if(yy<weg[p-1])yy=weg[p-1];
                  }
                  else {
                     yy=0;
                     break;
                  }
               }
            }
            if(xx<yy){xx=yy;r=p-1;}
         }
         k2+=6;
      }
      if(xx>0)cout << xx << " " << r << " " << this->show(r) << endl;
      sco+=xx;
      k1+=6;
   }
   sco/=(anm[d1]-1)*(anm[d2]-1);
   cout << "Returned value " << sco << endl;
   return(sco);
}

template<class Y,class Z>
void APost<Y,Z>::seeConsis(long nb,Y dn){
   long n,m,i,j,k=ofs[dn],flag;

   m=ofs[dn];
   for(i=0;i<anm[dn];i++){
      n=m+6*i;
      if(atd[n+1]==(nb+1)){
         cout << "1 " << this->show(nb) << endl;
         if(atd[n+2])cout << "2 " << this->show(atd[n+2]-1) << endl;
         if(atd[n+3])cout << "3 " << this->show(atd[n+3]-1) << endl;
         if(atd[n+4])cout << "4 " << this->show(atd[n+4]-1) << endl;
         if(atd[n+5])cout << "5 " << this->show(atd[n+5]-1) << endl;
      }
   }
}

template<class Y,class Z>
long APost<Y,Z>::ReMrk(int nx,long nb,Y *mrk){
   Y i,j,k,m,n,sb=0;
   long *fqq,cf=0,ttm=0,tless=0;
   const char *ptr,*pch;
   string s1,s2;
   size_t u,v;

   mrk[nb]=nx;
   s2=this->show(nb);
   u=s2.size()-3;
   s1.assign(s2,0,u);
   ptr=s1.c_str();
   for(i=0;i<this->nwrd;i++){
      if(!strncmp(ptr,this->pLx->show_str(i),u)){
          tless+=this->freq[i];
          mrk[i]=0;
      }
      mark(1,i,10000000,"auths");
   }
   return(tless);
}

template<class Y,class Z>
long APost<Y,Z>::TotalNames(void){
   Y i,j,k,m,sb=0;
   long *fqq,cf=0,ttm=0,tNam=0;
   char *ptr,*pch,ch;
   string s1,s2;
   size_t u,v;

   for(i=0;i<this->ndoc;i++)ttm+=anm[i];
   return(ttm);
}

}
#endif
