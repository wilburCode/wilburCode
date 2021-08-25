#ifndef ELEV_H
#define ELEV_H

#include <iostream>
#include <fstream>
#include <Vnab.h>
#include <DataObj.h>

using namespace std;
namespace iret {

class Resamp {
   public:
      Resamp(long num); //num is the size of the space
      ~Resamp(void);
      void zero(void); //Sets all counts to zero.
      void sample(long m,Index *pind); //Samples m times from set defined by
         //pind. Records counts in zsamp.
      void origin(Index *pind); //Sets the pind points a count of 1.
      long num; //size of space (should agree with Elev.num).
      long *zsamp; //Array of size num to hold counts.
};

class Elev {
   public:
      Elev(long nm); //nm is the number of entities that are ranked.
         //Usually this is the number of objects in the database. 
      ~Elev();

      void load(ifstream &fin); //Loads in the relevance data for a single 
         //query. Uses a standard format. Begins with query# and # of items.
         //This is then followed by that many items in form of entity index
         //number followed by relevance (between 0 and 1). 
      void load(Index *ind); //Loads in the relevance data where the
         //entity numbers in ind->idx[i] are the relevanant ones and all
         //have a relevance of 1.0.

      void process_ranks(long *idx); //idx is the array of indices of the
         //ranked set of num entities. The arrays pc and sc and tm are updated.
      void process_ranks(long ix,long *idx); //idx is the array of indices of the
         //ranked set of ix<=num entities. The arrays pc and sc and tm are updated.
      void process_ranks(Order *pord); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.
      void recall_prec_graph(double rc1,double rc2,Order *pord,ofstream &fout); 
         //Prints out the recall-precision pair (plus the score) for those 
         //recalls that fall in the range between rc1 and rc2. Starts at the 
         //top of the order and puts out a pair for each relevant data point 
         //encountered.
      void roc_graph(double rc1,double rc2,Order *pord,ofstream &fout); 
         //Prints out the recall-false drop pair (plus the score) for those 
         //recalls that fall in the range between rc1 and rc2. Starts at the 
         //top of the order and puts out a pair for each relevant data point 
         //encountered.
      double roc_scor(Order *pord); 
         //Returns the roc score for the data
      double ave_prec(Order *pord); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.
      double prcrec_BreakEven(Order *pord); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.
      double precision(long rnk,Order *pord); //Returns the precision over 
         //the ranks 1 down to rnk.
      double recall(long rnk,Order *pord); //Returns the recall over 
         //the ranks 1 down to rnk.
      double current(void); //Uses pc to compute spl (after call to 
         //process_ranks(). Returns spl.
      double summary(void); //Uses tm and sc to compute tc. Also sets spc
         //and returns it as a value.
      
      //Processing for paired bootstrap testing
      void rev_load(Resamp &Rsp); //Recalculates pt based on zsamp and rev.
      void process_ranks(long *idx,Resamp &Rsp); //idx is the array of indices of the
         //ranked set of num entities. The arrays pc and sc and tm are updated.
      void process_ranks(long ix,long *idx,Resamp &Rsp); //idx is the array of indices of the
         //ranked set of ix<=num entities. The arrays pc and sc and tm are updated.
      void process_ranks(Order *pord,Resamp &Rsp); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.

      //Error rate function
      double Error(double *sx); //After the load function has been called
          //this function will return the error rate based on Sign(sx) for
          //each point and num points.
      double Thresh(Order *pord); //Computes the optimal threshhold to use
          //with the scoring array.
      double Error(double thr,double *sx); //Computes the error rate based on
          //sx>thr is counted success and sx<=thr is counted failure.

      long num; //Number of documents in the database.
      double *rev; //Array of relevances of the entities in the database.
      double pt; //Holds the sum of relevances from rev.
      double *px; //Array of probabilities of relevance (isotonic).
      
      double rc[11]; //The eleven recall points.
      double pc[11]; //The eleven precisions. For one query.
      long tm; //Total number of queries thus far entered.
      double sc[11]; //The running total of eleven precisions 
         //for current set of queries that have been entered.
      double tc[11]; //Produced by summary() as the current summary
         //eleven precisions.
      double spl; //Local 11-point average precision just for current query
         //produced from pc by call to current().
      double spc; //Summary 11-point average precision 
};

template<class Y,class Z>
class EvalX {
   public:
      EvalX(Y nm); //nm is the number of entities that are ranked.
         //Usually this is the number of objects in the database. 
      ~EvalX(void);

      void load(Indx<Y> *ind); //Loads in the relevance data where the
         //entity numbers in ind->idx[i] are the relevanant ones and all
         //have a relevance of 1.0.
      void unload(Indx<Y> *ind); //Unloads the relevance data where the
         //entity numbers in ind->idx[i] were the relevant ones and all
         //are reset to 0.

      void recall_prec_graph(Z rc1,Z rc2,Ordr<Y,Z> *pOrd,ofstream &fout); 
         //Prints out the recall-precision pair (plus the score) for those 
         //recalls that fall in the range between rc1 and rc2. Starts at the 
         //top of the order and puts out a pair for each relevant data point 
         //encountered.
      void roc_graph(Z rc1,Z rc2,Ordr<Y,Z> *pOrd,ofstream &fout); 
         //Prints out the recall-false drop pair (plus the score) for those 
         //recalls that fall in the range between rc1 and rc2. Starts at the 
         //top of the order and puts out a pair for each relevant data point 
         //encountered.
      Z roc_scor(Ordr<Y,Z> *pOrd); 
         //Returns the roc score for the data
      Z ave_prec(Ordr<Y,Z> *pOrd); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.
      Z ave_prec(Ordr<Y,Z> *pOrd,Z px); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.
         //px is relevant material that does not appear in pOrd
      Z prcrec_BreakEven(Ordr<Y,Z> *pOrd); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.
      Z prcrec_BreakEven(Ordr<Y,Z> *pOrd,Z px); //Must have pord->pInd->ix<=num.
         //This function avoids the problem that may arise with duplicate scores
         //that causes ambiguous ordering and unreliable measurements.
         //px is relevant material that does not appear in pOrd
      Z precision(Y rnk,Ordr<Y,Z> *pOrd); //Returns the precision over 
         //the ranks 1 down to rnk.
      Z recall(Y rnk,Ordr<Y,Z> *pOrd); //Returns the recall over 
         //the ranks 1 down to rnk.
      //Special for case of non-integral relevance numbers
      Z DCG(Ordr<Y,Z> *pOrd); //Must have pord->pInd->ix<=num.
         //Discounted cumulative gain with variable relevance values allowed
      Z nDCG(Ordr<Y,Z> *pOrd,Z &xDCG,Z &xIDCG); //Must have pord->pInd->ix<=num.
         //Normalized Discounted cumulative gain with variable relevance values allowed
         //Also computes DCG and Ideal DCG and returns in parameter list.
      //Special  in that rx <= pOrd->ix is necessary to get values for higher ranks 
      //than full data size
      Z DCG(Y rx,Ordr<Y,Z> *pOrd); //Must have pord->pInd->ix<=num.
         //Discounted cumulative gain with variable relevance values allowed
      Z nDCG(Y rx,Ordr<Y,Z> *pOrd,Z &xDCG,Z &xIDCG); //Must have pord->pInd->ix<=num.
         //Normalized Discounted cumulative gain with variable relevance values allowed
         //Also computes DCG and Ideal DCG and returns in parameter list.

      //Error rate function
      Z Error(Z *sx); //After the load function has been called
          //this function will return the error rate based on Sign(sx) for
          //each point and num points.
      Z Thresh(Ordr<Y,Z> *pOrd); //Computes the optimal threshhold to use
          //with the scoring array.
      Z Error(Z thr,Z *sx); //Computes the error rate based on
          //sx>thr is counted success and sx<=thr is counted failure.
      //Correlation functions
      Z Pearson(Y n,Z *sx,Z *sy); //standard method, n number of points
      Z Spearman(Y n,Z *sx,Z *sy); //standard method, handls equal values in
        //the score arrays

      Y num; //Number of documents in the database.
      Z *rev; //Array of relevances of the entities in the database.
      Z pt; //Holds the sum of relevances from rev.
      Z *px; //Array of probabilities of relevance (isotonic).
};

template<class Y,class Z>
EvalX<Y,Z>::EvalX(Y nm){
   num=nm;
   rev=new Z[num];
}

template<class Y,class Z>
EvalX<Y,Z>::~EvalX(void){
   delete [] rev;
}

template<class Y,class Z>
void EvalX<Y,Z>::load(Indx<Y> *ind){
   long i;

   pt=0;
   for(i=0;i<num;i++)*(rev+i)=0;
   if(ind){
      for(i=0;i<ind->ix;i++){
         pt+=*(rev+ind->idx[i])=1.0;
      }
   }
}

template<class Y,class Z>
void EvalX<Y,Z>::unload(Indx<Y> *ind){
   long i;

   if(ind){
      for(i=0;i<ind->ix;i++){
         *(rev+ind->idx[i])=0;
      }
   }
}

template<class Y,class Z>
Z EvalX<Y,Z>::ave_prec(Ordr<Y,Z> *pord){
   Y i,j,k,u,ix;
   Z xx,tt,pz,rz;
   Z rtt=0,sum,ptt,ftt=0;
   Z so,si;

   if(pt<=0)return(0);
   if(!pord->ix)return(0);

   ix=pord->ix;
   u=pord->ind(0,so);
   sum=rev[u];
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         if(sum>0){
            ptt=0;
            if(sum>1.0)rz=(Z)(sum-1.0)/(xx-1.0);
            else rz=0;
            for(k=j;k<i;k++)ptt+=(rtt+rz*(k-j)+1.0)/((Z)k+1.0);
            ftt+=sum*ptt/xx;
            rtt+=sum;
         }
         sum=0;
         j=i;
         so=si;
      }
      sum+=rev[u];
   }
   xx=(Z)(ix-j);
   if(sum>0){
      ptt=0;
      if(sum>1.0)rz=(Z)(sum-1.0)/(xx-1.0);
      else rz=0;
      for(k=j;k<ix;k++)ptt+=(rtt+rz*(k-j)+1.0)/((Z)k+1.0);
      ftt+=sum*ptt/xx;
      rtt+=sum;
   }
   return(ftt/pt);
}

template<class Y,class Z>
Z EvalX<Y,Z>::ave_prec(Ordr<Y,Z> *pord,Z px){
   Y i,j,k,u,ix;
   Z xx,tt,pz,rz;
   Z rtt=0,sum,ptt,ftt=0;
   Z so,si;

   if(pt<=0)return(0);
   if(!pord->ix)return(0);

   ix=pord->ix;
   u=pord->ind(0,so);
   sum=rev[u];
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         if(sum>0){
            ptt=0;
            if(sum>1.0)rz=(Z)(sum-1.0)/(xx-1.0);
            else rz=0;
            for(k=j;k<i;k++)ptt+=(rtt+rz*(k-j)+1.0)/((Z)k+1.0);
            ftt+=sum*ptt/xx;
            rtt+=sum;
         }
         sum=0;
         j=i;
         so=si;
      }
      sum+=rev[u];
   }
   xx=(Z)(ix-j);
   if(sum>0){
      ptt=0;
      if(sum>1.0)rz=(Z)(sum-1.0)/(xx-1.0);
      else rz=0;
      for(k=j;k<ix;k++)ptt+=(rtt+rz*(k-j)+1.0)/((Z)k+1.0);
      ftt+=sum*ptt/xx;
      rtt+=sum;
   }
   return((ftt/rtt)*(pt/(pt+px)));
}

template<class Y,class Z>
Z EvalX<Y,Z>::DCG(Ordr<Y,Z> *pord){
   Y i,j,k,u,v,ix;
   Z xx,tt,pz,*rz;
   Z rtt=0,sum,ptt,ftt=0;
   Z so,si;

   if(pt<=0)return(0);
   if(!pord->ix)return(0);

   ix=pord->ix;
   rz=new Z[ix];
   for(i=0;i<ix;i++){
      u=pord->ind(i,si);
      rz[i]=rev[u];
   }
   sum=rz[0];
   u=pord->ind(0,so);
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         sum/=xx;
         for(k=j;k<i;k++){
            rz[k]=sum;
         }
         sum=0;
         j=i;
         so=si;
      }
      sum+=rz[i];
   }
   if(i>j){
      xx=(Z)(i-j);
      sum/=xx;
      for(k=j;k<i;k++){
         rz[k]=sum;
      }
   }
      
   sum=0;
   for(i=0;i<ix;i++){
      sum+=rz[i]/(l2*log(i+2));
   }
   delete [] rz;
   return(sum);
}

template<class Y,class Z>
Z EvalX<Y,Z>::DCG(Y rx,Ordr<Y,Z> *pord){ //rx <=pord->ix is necessary
   Y i,j,k,u,v,ix;
   Z xx,tt,pz,*rz;
   Z rtt=0,sum,ptt,ftt=0;
   Z so,si;

   if(pt<=0)return(0);
   if(!pord->ix)return(0);

   ix=pord->ix;
   rz=new Z[ix];
   for(i=0;i<ix;i++){
      u=pord->ind(i,si);
      rz[i]=rev[u];
   }
   sum=rz[0];
   u=pord->ind(0,so);
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         sum/=xx;
         for(k=j;k<i;k++){
            rz[k]=sum;
         }
         sum=0;
         j=i;
         so=si;
      }
      sum+=rz[i];
   }
   if(i>j){
      xx=(Z)(i-j);
      sum/=xx;
      for(k=j;k<i;k++){
         rz[k]=sum;
      }
   }

   sum=0;
   for(i=0;i<rx;i++){
      sum+=rz[i]/(l2*log(i+2));
   }
   delete [] rz;
   return(sum);
}

template<class Y,class Z>
Z EvalX<Y,Z>::nDCG(Ordr<Y,Z> *pord,Z &xDCG,Z &xIDCG){
   Y i,j,k,u,v,ix;
   Z xx,tt,pz,*rz;
   Z rtt=0,sum,ptt,ftt=0;
   Z so,si;

   if(pt<=0){
      xDCG=xIDCG=0;
      return(0);
   }
   if(!pord->ix){
      xDCG=xIDCG=0;
      return(0);
   }

   ix=pord->ix;
   Ordr<Y,Z> *pRrd=new Ordr<Y,Z>(ix,num,rev);
   cout << endl;
   xDCG=DCG(pord);
   xIDCG=DCG(pRrd);
   delete pRrd;
   if(xIDCG>0)return(xDCG/xIDCG);
   else return(0);
}

template<class Y,class Z>
Z EvalX<Y,Z>::nDCG(Y rx,Ordr<Y,Z> *pord,Z &xDCG,Z &xIDCG){ //rx <=pord->ix is necessary
   Y i,j,k,u,v,ix;
   Z xx,tt,pz,*rz;
   Z rtt=0,sum,ptt,ftt=0;
   Z so,si;

   if(pt<=0){
      xDCG=xIDCG=0;
      return(0);
   }

   ix=pord->ix;
   Ordr<Y,Z> *pRrd=new Ordr<Y,Z>(ix,num,rev);
   cout << endl;
   xDCG=DCG(rx,pord);
   xIDCG=DCG(rx,pRrd);
   delete pRrd;
   if(xIDCG>0)return(xDCG/xIDCG);
   else return(0);
}

template<class Y,class Z>
Z EvalX<Y,Z>:: prcrec_BreakEven(Ordr<Y,Z> *pord){
   Y i,j,k,u,ix;
   Z tt,xx;
   Z *rnw,sum;
   Z so,si;
   Z yo,y;

   if(pt<=0)return(0);
   if(!pord->ix)return(0);
   ix=pord->ix;
   rnw=new Z[ix];
   u=pord->ind(0,so);
   sum=rev[u];
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         for(k=j;k<i;k++)rnw[k]=sum/xx;
         sum=0;
         j=i;
         so=si;
      }
      sum+=rev[u];
   }
   xx=(Z)(i-j);
   for(k=j;k<i;k++)rnw[k]=sum/xx;
   yo=0;
   u=(ix<pt)?ix:rnd(pt);
   for(i=0;i<u;i++){
       yo+=rnw[i];
   }
   xx=yo/pt;
   delete [] rnw;
   return (xx);
}

template<class Y,class Z>
Z EvalX<Y,Z>:: prcrec_BreakEven(Ordr<Y,Z> *pord,Z px){
   Y i,j,k,u,ix;
   Z tt,xx;
   Z *rnw,sum;
   Z so,si;
   Z yo,y;

   if(pt<=0)return(0);
   if(!pord->ix)return(0);
   ix=pord->ix;
   rnw=new Z[ix];
   u=pord->ind(0,so);
   sum=rev[u];
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         for(k=j;k<i;k++)rnw[k]=sum/xx;
         sum=0;
         j=i;
         so=si;
      }
      sum+=rev[u];
   }
   xx=(Z)(i-j);
   for(k=j;k<i;k++)rnw[k]=sum/xx;
   yo=0;
   u=(ix<(pt+px))?ix:rnd(pt+px);
   for(i=0;i<u;i++){
       yo+=rnw[i];
   }
   xx=yo/(pt+px);
   delete [] rnw;
   return (xx);
}

template<class Y,class Z>
Z EvalX<Y,Z>::precision(Y rnk,Ordr<Y,Z> *pord){
   Y i,j,k,u,ix;
   Z xx,tt,pz,rz;
   Z *rnw,sum,ptt;
   Z so,si;

   ix=pord->ix;
   rnw=new Z[ix];
   u=pord->ind(0,so);
   sum=rev[u];
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         for(k=j;k<i;k++)rnw[k]=sum/xx;
         sum=0;
         j=i;
         so=si;
      }
      sum+=rev[u];
   }
   xx=(Z)(i-j);
   for(k=j;k<i;k++)rnw[k]=sum/xx;

   u=(rnk<pord->ix)?rnk:pord->ix;
   tt=0;
   for(i=0;i<u;i++){
      tt+=*(rnw+i);
   }
   ptt=tt/((Z)rnk);
   delete [] rnw;
   return(ptt);
}

template<class Y,class Z>
Z EvalX<Y,Z>::recall(Y rnk,Ordr<Y,Z> *pord){
   Y i,j,k,u,ix;
   Z xx,tt,pz,rz;
   Z *rnw,sum,ptt;
   Z so,si;

   if(pt==0.0)return(0);

   ix=pord->ix;
   rnw=new Z[ix];
   u=pord->ind(0,so);
   sum=rev[u];
   j=0;
   for(i=1;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         xx=(Z)(i-j);
         for(k=j;k<i;k++)rnw[k]=sum/xx;
         sum=0;
         j=i;
         so=si;
      }
      sum+=rev[u];
   }
   xx=(Z)(i-j);
   for(k=j;k<i;k++)rnw[k]=sum/xx;

   u=(rnk<pord->ix)?rnk:pord->ix;
   tt=0;
   for(i=0;i<u;i++){
      tt+=*(rnw+i);
   }
   ptt=tt/pt;
   delete [] rnw;
   return(ptt);
}

template<class Y,class Z>
void EvalX<Y,Z>::recall_prec_graph(Z rc1,Z rc2,Ordr<Y,Z> *pord,ofstream &fout){
   Y i,j,k,u,ix=pord->ix;
   Z sum=0,rss;
   Z si;

   for(i=0;i<ix;i++){
      u=pord->ind(i,si);
      if(rev[u]>0){
         sum+=rev[u];
         rss=sum/pt;
         if((rc1<=rss)&&(rss<=rc2)){
            fout << rss << " " << sum/((double)(i+1)) << " " << si << endl;
         }
      }
   }
}

template<class Y,class Z>
void EvalX<Y,Z>::roc_graph(Z rc1,Z rc2,Ordr<Y,Z> *pord,ofstream &fout){
   Y i,j,k,u,ix=pord->ix;
   Z sum=0,rss=0,xss;
   Z si;

   xss=ix-pt;
   for(i=0;i<ix;i++){
      u=pord->ind(i,si);
      if(rev[u]>0){
         sum+=rev[u];
         rss=sum/pt;
      }
      else {
         if((rc1<=rss)&&(rss<=rc2)){
            fout << (i+1-sum)/xss << " " << rss << " " << si << endl;
         }
      }
   }
}

template<class Y,class Z>
Z EvalX<Y,Z>::roc_scor(Ordr<Y,Z> *pord){
   Y i,j,k,u,ix=pord->ix;
   Y ip,ig;
   Z xx,tt,pz,rz;
   Z sum=0,xss,smx=0,rss=0;
   Z so,si;

   if(pt==0.0)return(0);

   u=pord->ind(0,so);
   so++;
   xss=ix-pt;
   ip=ig=0;
   for(i=0;i<ix;i++){
      u=pord->ind(i,si);
      if(so>si){
         if(ig){
            smx+=ig*(0.5*ip+sum)/pt;
         }
         sum+=ip;
         so=si;
         ip=ig=0;
      }
      if(rev[u]>0)ip++;
      else ig++;
   }
   if(ig){
      smx+=ig*(0.5*ip+sum)/pt;
   }
   return(smx/xss);
}

template<class Y,class Z>
Z EvalX<Y,Z>::Error(Z *sx){
   Y i;
   Z sum=0;

   for(i=0;i<num;i++){
      if((rev[i]>0.5)&&(sx[i]>0))sum++;
      else if((rev[i]<0.5)&&(sx[i]<=0))sum++;
   }
   return((num-sum)/num);
}

template<class Y,class Z>
Z EvalX<Y,Z>::Thresh(Ordr<Y,Z> *pord){
   Y i,j,k=-1;
   Z xx,zz,er,rr;
   Z xc,th;

   pord->ind(0,xc);
   xc++;
   th=xc;
   er=rr=pt;
   i=0;
   while(i<num){
      j=pord->ind(i,zz);
      if(rev[j]>0.5)rr--;
      else rr++;
      if(zz<xc){
         if(rr<er){
            th=zz;
            er=rr;
            k=i;
         }
         xc=zz;
      }
      i++;
   }
   if((k>-1)&&(k<num-1)){
      pord->ind(k+1,zz);
      th=0.5*(th+zz);
   }
   else if(k==num-1)th--;
   return(th);
}

template<class Y,class Z>
Z EvalX<Y,Z>::Error(Z thr,Z *sx){
   Y i;
   Z sum=0;

   for(i=0;i<num;i++){
      if((rev[i]>0.5)&&(sx[i]>thr))sum++;
      else if((rev[i]<0.5)&&(sx[i]<=thr))sum++;
   }
   return((num-sum)/num);
}

template<class Y,class Z>
Z EvalX<Y,Z>::Pearson(Y n,Z *sx,Z *sy){
   Y i,j,k;
   Z uu,vv,xx,yy,zz;
  
   xx=yy=0;
   for(i=0;i<n;i++){
      xx+=sx[i];
      yy+=sy[i];
   }
   xx/=n;
   yy/=n;
   uu=vv=zz=0;
   for(i=0;i<n;i++){
      uu+=(sx[i]-xx)*(sx[i]-xx);
      vv+=(sy[i]-yy)*(sy[i]-yy);
      zz+=(sx[i]-xx)*(sy[i]-yy);
   }
   return(zz/(sqrt(uu)*sqrt(vv)));
}

template<class Y,class Z>
Z EvalX<Y,Z>::Spearman(Y n,Z *su,Z *sv){
   Y i,j,k;
   Z uu,vv,xx,yy,zz;
   Z *rmk=new Z[n];
   Z *rnk=new Z[n];
   Z *sx=new Z[n];
   Z *sy=new Z[n];
   Y *orm=new Y[n];
   Y *orn=new Y[n];
   Y *rmv=new Y[n];
   Y *rnv=new Y[n];
   for(i=0;i<n;i++){
      sx[i]=su[i];
      sy[i]=sv[i];
      orm[i]=orn[i]=i;
   }
   hSort(n,sx,orm);
   hSort(n,sy,orn);
   for(i=0;i<n;i++){
     rmv[orm[i]]=i;
     rnv[orn[i]]=i;
   }
   //sx computation
   xx=sx[0];
   zz=0;
   k=0;
   for(i=1;i<n;i++){
      if(sx[i]==xx){
         zz+=i;
      }
      else {
         uu=zz/(i-k);
         for(j=k;j<i;j++)rmk[j]=uu;
         xx=sx[i];
         zz=i;
         k=i;
      }
   } 
   i=n;
   uu=zz/(i-k);
   for(j=k;j<i;j++)rmk[j]=uu;
   //sy computation
   xx=sy[0];
   zz=0;
   k=0;
   for(i=1;i<n;i++){
      if(sy[i]==xx){
         zz+=i;
      }
      else {
         uu=zz/(i-k);
         for(j=k;j<i;j++)rnk[j]=uu;
         xx=sy[i];
         zz=i;
         k=i;
      }
   } 
   i=n;
   uu=zz/(i-k);
   for(j=k;j<i;j++)rnk[j]=uu;
   for(i=0;i<n;i++){
      sx[i]=rmk[rmv[i]];
      sy[i]=rnk[rnv[i]];
   }
   delete [] rmk;
   delete [] rnk;
   delete [] orm;
   delete [] orn;
   delete [] rmv;
   delete [] rnv;
   zz=Pearson(n,sx,sy);
   delete [] sx;
   delete [] sy;
   return(zz);
}

}
#endif
