#ifndef MCSVMG_H
#define MCSVMG_H

#include <iostream>
#include <fstream>
#include <DataObj.h>
#include <XPost.h>
#include <Hyper.h>
#include <Elev.h>
#include <CMark.h>

using namespace std;
namespace iret {

//template class MCSvmG
template<class Y,class Z>
class MCSvmG : public CMark<Y,Z> {
public:
   MCSvmG(const char *nspost);//name of XPost set
   MCSvmG(const char *nspost,const char *pnam);//name of SPost set
      //pnam is either used in "path_pnam" as place to find path
      //or if begins with ':' is followed by path itself.
   ~MCSvmG(void);

   //Learning functions
   void Init_System(Y m,Indx<Y> **mx); //m is number of classes
      //mx defines the classes. For class i, mx[i] is the index object
      //that contains the numbers of the docs of that class in the training
      //data. The part of the docset not contained in some mx[i] is what is
      //left over as the test set.

   void LearnSGD(Y iter1,Y iter2); //Simple SGD like update. Optimize with lam=0 and early stopping
     //iter2 is number of learning cycles. Find iter=8 is good for large MEDLINE collections, 
     //iter=9 for smaller training sets. iter1 is number of constraint addiition cycles.

   Y *ClassAll(void); //Puts the class number in the cls array and returns pointer to it 
   Y *ClassSet(Indx<Y> *ind);

   //Relate to the documents
   Y tdoc; //Number of training documents 
   bool **ca; //Training selection array
   Y *markd; //Marks the class with 0 to mc-1, mc is not training
   Y mc; //Number of classes

   //Relate to the terms processing
   Z lambda; //lambda
   Z *wt; //weights of the process
   Z th; //Threshhold
   Y *sr;  //Maps the tdoc training set
 
   Y *cls; //Holds the class labels as integers 0 to mc-1 inclusive
};

template<class Y,class Z>
MCSvmG<Y,Z>::MCSvmG(const char *namspost) : CMark<Y,Z>(namspost){
   wt=NULL;
   ca=NULL;
   sr=NULL;
   cls=NULL;
   markd=NULL;
}

template<class Y,class Z>
MCSvmG<Y,Z>::MCSvmG(const char *namspost,const char *pnam) : CMark<Y,Z>(namspost,pnam){
   wt=NULL;
   ca=NULL;
   sr=NULL;
   cls=NULL;
   markd=NULL;
}

template<class Y,class Z>
MCSvmG<Y,Z>::~MCSvmG(){
   if(wt!=NULL)this->dst_Mmap("wt",(char*&)wt);
   if(ca!=NULL){
      for(i=0;i<this->ndoc;i++){
         if(markd[i]<mc)delete [] ca[i];
      }
      delete [] ca;
   }
   if(sr!=NULL)delete [] sr;
   if(cls!=NULL)delete [] cls;
   if(markd!=NULL)delete [] markd;
}

template<class Y,class Z>
void MCSvmG<Y,Z>::Init_System(Y m,Indx<Y> **mx){
   Y i,j,u,*pm;
   bool *ux;

   this->gopen_db_map();

   mc=m;
   tdoc=0;
   for(i=0;i<m;i++)tdoc+=mx[i]->ix;
   if(this->ndoc<tdoc){cout << "Error in size of training set" << endl;exit(0);}  

   if(markd) delete [] markd;
   markd=new Y[this->ndoc];
   for(i=0;i<this->ndoc;i++)markd[i]=mc;
   for(i=0;i<m;i++){
      u=mx[i]->ix;
      pm=mx[i]->idx;
      for(j=0;j<u;j++)markd[pm[j]]=i;
   }
   ca=new bool*[this->ndoc];
   for(i=0;i<this->ndoc;i++){
      if(markd[i]<mc){
         ux=ca[i]=new bool[mc];
         for(j=0;j<mc;j++)ux[j]=false;
         while((j=zrand(mc))==markd[i]);
         ux[j]=true;
      }
   }
}

template<class Y,class Z>
void MCSvmG<Y,Z>::LearnSGD(Y iter1,Y iter2){
   Y i,j,k,n,sc,iz,iu,rp,ix;
   Y flag,mode,ct,gt,cwt,i1,i2,cwn;
   Z xx,yy,cv,*sumx,*w1,*w2,losx;
   Z sum1,sum2,zz,uu,eta=0.002;
   bool *ux;

   sumx=new Z[mc];
   this->gopen_db_map();
   if(sr!=NULL)delete [] sr;

   cwn=this->nwrd;
   cwt=cwn+1;
   wt=new Z[cwt];
   for(i=0;i<cwt;i++)wt[i]=0;
   ofstream *pfout=this->get_Ostr("wt",ios::out);
   for(i=0;i<mc;i++){
      pfout->write((char*)wt,cwt*sizeof(Z));
   }
   this->dst_Ostr(pfout);
   delete [] wt;
   wt=(Z)this->get_Wmap("wt");

   sr=new Y[mc*tdoc];

   for(rp=0;rp<iter1;rp++){ //Outer loop iter over adding constraints
      //Set the list of training examples
      n=0;
      for(i=0;i<this->ndoc;i++){
         if(markd[i]<mc){
            ux=ca[i];
            k=mc*i;
            for(j=0;j<mc;j++){
               if(ux[j]){
                  sr[n++]=k+j;
               }
            }
         }
      }
      //Zero the wt array
      for(i=0;i<cwt*mc;i++)wt[i]=0;

      for(gt=0;gt<iter2;gt++){ //Inner loop for SGD learning
         xshuffle(n,sr);
         losx;
         for(i=0;i<n;i++){
            j=sr[i]/mc;
            i1=markd[j];
            i2=sr[i]%mc;
            w1=wt+i1*cwt;
            w2=wt+i2+cwt;
            readp_db(j);
            sum1=*(w1+cwn);
            sum2=*(w2+cwn);
            for(iz=0;iz<this->nw;iz++){
               if(this->mrk[iu=this->nwd[iz]]){
                  sum1+=w1[iu];
                  sum2-=w2[iu];
               }
            }
            //Loss zz calculation from p=sum
            if(sum1<2.0+sum2){
               losx+=2.0+sum2-sum1;
               *(w1+cwn)+=eta;
               *(w2+cwn)+=eta;
               for(iz=0;iz<this->nw;iz++){
                  if(this->mrk[iu=this->nwd[iz]]){
                     w1[iu]+=eta;
                     w2[iu]-=eta;
                  }
               }
            }
         }
         if(this->pflag)cout << "E(tr loss) " << losx/n << " iter2 " << gt+1 << endl;
      }
      //Adding constraints
      for(i=0;i<this->ndoc;i++){
         if((i1=markd[i])<mc){
            readp_db(i);
            for(j=0;j<mc;j++)sumx[j]=*(wt+j*cwt+cwn);
            for(iz=0;iz<this->nw;iz++){
               if(this->mrk[iu=this->nwd[iz]]){
                  for(j=0;j<mc;j++){
                     sumx[j]+=wt[j*cwt+iu];
                  }
               }
            }
            sum2=sumx[i1]-2.0;
            i2=-1;
            for(j=0;j<i1;j++){
               if(sum2<sumx[j]){
                  sum2=sumx[j];
                  i2=j;
               }
            }   
            for(j=i1+1;j<mc;j++){
               if(sum2<sumx[j]){
                  sum2=sumx[j];
                  i2=j;
               }
            }   
            //Loss zz calculation from p=sum
            if((i2>-1)&&(sumx[i1]<2.0+sum2)){
               ca[i][i2]=true;
            }
         }
      }
      if(this->pflag)cout << "iter1 " << rp+1 << endl;
   }
   delete [] sr;
   sr=NULL;
   delete [] sumx;
   this->mak_Msync("wt",(char*)wt);
}

template<class Y,class Z>
Y *MCSvmG<Y,Z>::ClassAll(void){
   Y i,j,n,iz,iu,i2,cwn,cwt;
   Z *sumx,sum2;
   
   this->gopen_db_map();
   cwn=this->nwrd;
   cwt=cwn+1;
   sumx=new Z[mc];
   if(cls!=NULL)delete [] cls;
   cls=new Y[this->ndoc];
      
   for(i=0;i<this->ndoc;i++){
      readp_db(i);
      for(j=0;j<mc;j++)sumx[j]=*(wt+j*cwt+cwn);
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            for(j=0;j<mc;j++){
               sumx[j]+=wt[j*cwt+iu];
            }
         }
      }
      sum2=sumx[0];
      i2=0;
      for(j=1;j<mc;j++){
         if(sum2<sumx[j]){
            sum2=sumx[j];
            i2=j;
         }
      }
      cls[i]=i2;
      this->mark(j+1,100,"docs scored");
   }
   return(cls);
}  

template<class Y,class Z>
Y *MCSvmG<Y,Z>::ClassSet(Indx<Y> *ind){
   Y i,j,n,iz,iu,i2,cwn,cwt;
   Z *sumx,sum2;

   this->gopen_db_map();
   cwn=this->nwrd;
   cwt=cwn+1;
   sumx=new Z[mc];
   if(cls!=NULL)delete [] cls;
   cls=new Y[this->ndoc];

   for(n=0;n<ind->ix;n++){
      i=ind->idx[n];
      readp_db(i);
      for(j=0;j<mc;j++)sumx[j]=*(wt+j*cwt+cwn);
      for(iz=0;iz<this->nw;iz++){
         if(this->mrk[iu=this->nwd[iz]]){
            for(j=0;j<mc;j++){
               sumx[j]+=wt[j*cwt+iu];
            }
         }
      }
      sum2=sumx[0];
      i2=0;
      for(j=1;j<mc;j++){
         if(sum2<sumx[j]){
            sum2=sumx[j];
            i2=j;
         }
      }
      cls[i]=i2;
      this->mark(n+1,100,"docs scored");
   }
   return(cls);
}

}
#endif
