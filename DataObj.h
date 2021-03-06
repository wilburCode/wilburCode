#ifndef DATAOBJ_H
#define DATAOBJ_H

#include <runn.h>
#include <FBase.h>
#include <set>

using namespace std;
namespace iret {

class Order; //Forward declaration.
class Indes; //Forward declaration.

class Index {
   public:
      Index(void); //Simple constructor;
      Index(long num); //Sets up memory;
      Index(long begin,long end_plus_one);
      Index(long num,long *imd,long oflag); //Sets up memory and populates. Orders
         //if oflag=1. Orders and uniques if oflag=2.
      Index(Index *ind); //copy constructor
      Index(Indes *ins); //copy constructor
      ~Index(void); //Frees memory if pointer not null.
      void sSort(void); //Performs a simple heapsort on idx.
      long unique(void); //Checks for duplicates and removes.
        //Must be ordered first.
      Index *cbool_And(Index *jnd);//Ands with the argument
      Index *cbool_Or(Index *jnd);//Ors with the argument
      Index *cbool_Butnot(Index *jnd);//Removes the argument

        //Identifies a subset by indices
      long Subvalue(long j); //on success returns value k+1 such that idx[k]=j.
          //Otherwise returns 0.
      Index *Subvalue(Index *jnd); //jnd must represent a subset and the 
          //returned Index object will give the values k such that idx[k]
          //is among jnd->idx values.
      long Lbound(long n); //Returns the smallest value j such that
          //n<=idx[j]. If there is no such j it returns ix.
      long Ubound(long n); //Returns the largest value j such that
          //idx[j]<n. If there is no such j it returns 0
      void Compose(Index *pKnd); //Assumes that the current objects
          //idx values are less than pKnd->ix and does function
          //composition.

      Index *Subinterval(long n,long m); //Makes the index from
         //interval n<=x<m (idx[n] to idx[m-1]).
      Index *Subsample(long n,long seed); //Generates a random sample of
         //idx of size n if ix>n, else returns this.
      Index *Subsample(long n); //Generates a random sample of
         //idx of size n if ix>n, else returns this.
         //Uses current state of random number generator

        //Subset by threshold sxx[n] and idx[n] defined.
      Index *Greater(float *sxx,float thresh); //All values such
          //that *(sxx+n)>thresh.
      Index *Greateq(float *sxx,float thresh); //All values such
          //that *(sxx+n)>=thresh.
      Index *Greater(double *sxx,double thresh); //All values such
      Index *Greateq(double *sxx,double thresh); //All values such
      Index *Lesser(float *sxx,float thresh); //All values such
          //that *(sxx+n)<thresh.
      Index *Lesseq(float *sxx,float thresh); //All values such
          //that *(sxx+n)<=thresh.
      Index *Lesser(double *sxx,double thresh); //All values such
      Index *Lesseq(double *sxx,double thresh); //All values such
          //that *(sxx+n)<=thresh.

        //Subset by threshold sxx[idx] defined.
      Index *Greater(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)>thresh.
      Index *Greateq(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)>=thresh.
      Index *Greater(double thresh,double *sxx); //Subset of idx values 
      Index *Greateq(double thresh,double *sxx); //Subset of idx values 
      Index *Lesser(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)<thresh.
      Index *Lesseq(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)<=thresh.
      Index *Lesser(double thresh,double *sxx); //Subset of idx values 
      Index *Lesseq(double thresh,double *sxx); //Subset of idx values 
        //Orders by threshold sxx[idx] defined.
      Order *oGreater(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)>thresh made into an Order.
      Order *oGreateq(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)>=thresh made into an Order.
      Order *oGreater(double thresh,double *sxx); //Subset of idx values
      Order *oGreateq(double thresh,double *sxx); //Subset of idx values
      Order *oLesser(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)<thresh made into an Order.
      Order *oLesseq(float thresh,float *sxx); //Subset of idx values such
          //that *(sxx+idx)<=thresh made into an Order.
      Order *oLesser(double thresh,double *sxx); //Subset of idx values
      Order *oLesseq(double thresh,double *sxx); //Subset of idx values
      void Save(const char *nam); //Saves the current object in "index"_nam.x file
      void Load(const char *nam); //Loads the object "index"_nam.x into current
          //object
      void Save(const char *nam,const char *pnam); //Saves the current object in "index"_nam.x file
          //pnam is name associated with path_pnam where path is found
      void Load(const char *nam,const char *pnam); //Loads the object "index"_nam.x into current
          //object
          //pnam is name associated with path_pnam where path is found
      void write(ofstream &fout); //Writes the object into a file
      void read(ifstream &fin); //Reads the object in from a file
      void debug(void); //Writes object to standard out.

      long ix; //Number of elements.
      long *idx; //Elements in order.
};

class Indes {
   public:
      Indes(void); //Simple constructor;
      Indes(unsigned int num); //Sets up memory;
      Indes(unsigned int begin,unsigned int end_plus_one);
      Indes(Indes *ind); //copy constructor
      Indes(Index *ind); //copy constructor
      ~Indes(void); //Frees memory if pointer not null.

      unsigned int ix; //Number of elements.
      unsigned int *idx; //Elements in order.
};

class CValid {
public:
   CValid(Index *gdd, long n);
   CValid(Index *gdd, Index *bdd);
   ~CValid();
   void cross_valid(long m, long seed);
   Index *ind;
   Index *gnd;
   Index *bnd;
   Index **pGTS;
   Index **pBTS;
   Index **pWTS;
   Index **pGTR;
   Index **pBTR;
   Index **pWTR;
   long setn;

};

class Order {
   public:
      Order(void); //Simple constructor.
      Order(long n,long m,float *sco); //sco[i] is score for index i.
         //m is database & sco array size. Produces an Order of size n.
      Order(long n,long m,double *sco); //sco[i] is score for index i.
         //m is database & sco array size. Produces an Order of size n.
      Order(long n,Index *ind,float *sco); //sco[ind[i]] is interpretation.
         //Gives the top n.
      Order(long n,Index *ind,double *sco); //sco[ind[i]] is interpretation.
         //Gives the top n.
      Order(long n,float *sco,Index *ind); //idx[i] is index of
         //sco[i]. Produces an Order size n.
      Order(long n,double *sco,Index *ind); //idx[i] is index of
         //sco[i]. Produces an Order size n.
      Order(Order *pord); //copy constructor
      ~Order(void); //Removes data.
      Order *cbool_And(Index *jnd); //Ands with the argument but produces
         //an Order as output with the same scores and order restricted
         //to the subset.
      Order *cbool_Butnot(Index *jnd); //Removes the argument set and
         //produces an Order as output with the restricted parental
         //scores and order.
      Index *Topn(long n); //Returns an index object that contains the n
         //top scoring elements if n<=pInd->ix, else returns copy of pInd
      float Precision(long n,Index *ind); //Assumes that ind is a subset
         //of the index object of this order and returns the precision over
         //the top n scores. Deals with problem of identical scores.
      long CtGreateq(float thresh); //Returns the number of elements with
         //scores greater than or equal to thresh.
      double SVar(Order *qOrd); //Produces a score which is ratio of sum
         //of squares of scores of qOrd objects in this to sum of squares
         //of scores in all of qOrd
      double SVar(Order *qOrd,long m); //Produces a score which is ratio of sum
         //of squares of scores of top m qOrd objects in this to sum of squares
         //of scores of same in qOrd
      void write(ofstream &fout); //Writes the object into a file
      void read(ifstream &fin); //Reads the object in from a file
      void debug(void); //Prints data to standard out.
      void debugW(void); //Prints data to standard out.

      //Data Access
      long num(void);  //Returns number of objects in set.
      long ind(long i,float &sco); //Returns doc number (index) and
         //score pair that correlate and are in decreasing score order as i
         //increases from 0 to num-1.
      long inv(long i); //Returns the number j with i=order[j], else pInd->ix
      long *seq(void); //Returns pointer at new array that contains the array
         //pInd->idx in decreasing score order. 

      Index *pInd; //Holds the list in standard form of index set.
      long *order; //Holds the order of decreasing score.
      float *score; //Holds the scores in the ord order.
};

//This class can be used as convenient way to process a theme without any special terms
//Terms can be here by number as in Thln.h&C or by ID as in Focus.h&C
class WOrder : public Order {
  public:
   WOrder(void);
   WOrder(Order *pOrd,float *weg); //Uses weg[pOrd->pInd->idx[order[i]]]
      //to obtain the weights for terms. ith weight corresponds to ind(i,.)
   WOrder(Order *pOrd,double *weg);
      //Both these constructors destroy the Order object pointed
      //to by pOrd. A case of catabolism.
   ~WOrder(void);

   float ave_sco(long n); //Avererages the n top scores
   float dice_sco(WOrder *wpord); //Produces dice similarity score.
   int subset_equal(int n, WOrder *wpord); //Returns 1 if the two objects have the same terms
      //in the top n spots, order not important
   int positive(int n); //returns 1 if first n weights are positive, else 0
   //Data
   float *weight;
};

template<class Y,class Z> class Ordr;

//Template Index class
template<class Z>
class Indx {
   public:
      Indx(void); //Simple constructor;
      Indx(Z num); //Sets up memory;
      Indx(Z begin,Z end_plus_one); //Sets memory and fills with values in
         //the interval [begin,end_plus_one)
      Indx(Z num,Z *imd,long oflag); //Sets up memory and populates. Orders
         //if oflag=1. Orders and uniques if oflag=2.
      Indx(Indx<Z> *ind); //copy constructor
      Indx(set<Z> Sx); //constructs from set of numbers Sx
      ~Indx(void); //Frees memory if pointer not null.
      void sSort(void); //Performs a simple heapsort on idx.
      long unique(void); //Checks for duplicates and removes.
        //Must be ordered first.
      Indx<Z> *cbool_And(Indx<Z> *jnd);//Ands with the argument
      Indx<Z> *cbool_Or(Indx<Z> *jnd);//Ors with the argument
      Indx<Z> *cbool_Butnot(Indx<Z> *jnd);//Removes the argument
      Z CAnd(Indx<Z> *jnd);//Returns number of elements in common
      void Compose(Indx<Z> *pKnd); //Assumes that the current objects
          //idx values are less than pKnd->ix and does function
          //composition.
      inline void Subs(Indx<Z> *pKnd){ //Substitution
          ix=pKnd->ix;
          idx=pKnd->idx;
      }

        //Identifies a subset by indices
      Z Subvalue(Z j); //on success returns value k+1 such that idx[k]=j.
          //Otherwise returns 0.
      Indx<Z> *Subvalue(Indx<Z> *jnd); //jnd must represent a subset and the 
          //returned Indx object will give the values k such that idx[k]
          //is among jnd->idx values.
      Z Lbound(Z n); //Returns the smallest value j such that
          //n<=idx[j]. If there is no such j it returns ix.
      Z Ubound(Z n); //Returns the largest value j such that
          //idx[j]<n. If there is no such j it returns 0

      Indx<Z> *Subinterval(Z n,Z m); //Makes the index from
         //interval n<=x<m (idx[n] to idx[m-1]).
      Indx<Z> *Subsample(Z n,long seed); //Generates a random sample of
         //idx of size n if ix>n, else returns this.
      Indx<Z> *Subsample(Z n); //Generates a random sample of
         //idx of size n if ix>n, else returns this.
         //Uses current state of random number generator

        //Subset by threshold sxx[n] and idx[n] defined.
      template<typename U>
      Indx<Z> *Greater(U *sxx,U thresh); //All values such
          //that *(sxx+n)>thresh.
      template<typename U>
      Indx<Z> *Greateq(U *sxx,U thresh); //All values such
          //that *(sxx+n)>=thresh.
      template<typename U>
      Indx<Z> *Lesser(U *sxx,U thresh); //All values such
          //that *(sxx+n)<thresh.
      template<typename U>
      Indx<Z> *Lesseq(U *sxx,U thresh); //All values such
          //that *(sxx+n)<=thresh.

        //Subset by threshold sxx[idx] defined.
      template<typename U>
      Indx<Z> *Greater(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)>thresh.
      template<typename U>
      Indx<Z> *Greateq(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)>=thresh.
      template<typename U>
      Indx<Z> *Lesser(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)<thresh.
      template<typename U>
      Indx<Z> *Lesseq(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)<=thresh.

        //Orders by threshold sxx[idx] defined.
      template<typename U>
      Ordr<Z,U> *oGreater(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)>thresh made into an Order.
      template<typename U>
      Ordr<Z,U> *oGreateq(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)>=thresh made into an Order.
      template<typename U>
      Ordr<Z,U> *oLesser(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)<thresh made into an Order.
      template<typename U>
      Ordr<Z,U> *oLesseq(U thresh,U *sxx); //Subset of idx values such
          //that *(sxx+idx)<=thresh made into an Order.

      void Save(const char *nam); //Saves the current object in "index"_nam.x file
      void Load(const char *nam); //Loads the object "index"_nam.x into current
         //object
      void Save(const char *nam,const char *p_nam); //Saves the current object 
         //in "index"_nam.x file and uses p_nam as name for setting path
      void Load(const char *nam,const char *p_nam); //Loads the object 
         //"index"_nam.x into current object and uses p_nam as name for setting path
      void write(ofstream &fout); //Writes the object into a file
      void read(ifstream &fin); //Reads the object in from a file
      void debug(void); //Writes object to standard out.

      Z ix; //Number of elements.
      Z *idx; //Elements in order.
};

template<class Z>
Indx<Z>::Indx(void){
   ix=0;
   idx=NULL;
}  
   
template<class Z>
Indx<Z>::Indx(Z num){
   ix=num;
   idx=new Z[ix];
}  

template<class Z>
Indx<Z>::Indx(Z begin,Z end_plus_one){
   ix=end_plus_one-begin;
   idx=new Z[ix];
   Z i;
   for(i=0;i<ix;i++)idx[i]=begin+i;
}  
   
template<class Z>
Indx<Z>::Indx(Z num,Z *imd,long oflag){
   Z i;
      
   ix=num;
   idx=new Z[ix];
   for(i=0;i<num;i++)*(idx+i)=*(imd+i);
   if(oflag)sSort();
   if(oflag>1)unique();
}     
      
template<class Z>
Indx<Z>::Indx(Indx<Z> *ind){
   Z i;
      
   ix=ind->ix;
   idx=new Z[ix];
   for(i=0;i<ix;i++)idx[i]=ind->idx[i];
}

template<class Z>
Indx<Z>::Indx(set<Z> Sx){
   long i=0;

   ix=Sx.size();
   idx=new Z[ix];
   typename set<Z>::iterator bg,en;
   bg=Sx.begin();
   en=Sx.end();
   while(bg!=en){
      idx[i++]=*bg;
      bg++;
   }
   sSort();
}

template<class Z>
Indx<Z>::~Indx(void){
   if(idx!=NULL)delete [] idx;
}

template<class Z>
void Indx<Z>::sSort(void){
  Z k, j, ir, i;
  Z rra;

  if(ix<2)return;

  k=(ix>>1);
  ir=ix-1;
  for(;;) {
    if(k>0) {
      rra=idx[--k];
    }
    else {
      rra=idx[ir];
      idx[ir] = idx[0];
      if(--ir ==0) {
        idx[0]=rra;
        return;
      }
    }
    i=k;
    j=((k+1)<<1)-1;
    while(j<=ir) {
      if(j<ir && (idx[j]<idx[j+1])) ++j;
      if(rra<idx[j]) {
        idx[i]=idx[j];
        j +=(i=j)+1;
      }
      else j=ir+1;
    }
    idx[i]=rra;
  }
}

template<class Z>
long Indx<Z>::unique(void){
   Z i,k;
   Z j;

   if(ix<2)return(ix);
   k=1;
   j=*idx;
   for(i=1;i<ix;i++){
      if(j<*(idx+i))*(idx+(k++))=j=*(idx+i);
   }
   if(k<ix){
      Z *jdx=new Z[k];
      for(i=0;i<k;i++)*(jdx+i)=*(idx+i);
      delete [] idx;
      idx=jdx;
      ix=k;
   }
   return(k);
}

template<class Z>
Indx<Z> *Indx<Z>::cbool_And(Indx<Z> *jnd){
   Z i,j,k,su,w,p,m,bu;
   Z *pdx,*bdx,*sdx,bx,sx;

   if(jnd==NULL){
      return(NULL);
   }
   if(ix>jnd->ix){
      bdx=idx;
      bx=ix;
      sdx=jnd->idx;
      sx=jnd->ix;
   }
   else {
      bdx=jnd->idx;
      bx=jnd->ix;
      sdx=idx;
      sx=ix;
   }
   pdx=new Z[sx];
   i=j=k=0;
   while((i<sx)&&(j<bx)){
      bu=*(bdx+j);
      while((i<sx)&&(*(sdx+i)<bu))i++;
      if(i<sx){
         su=*(sdx+i);
         if(su==bu){
            *(pdx+k)=bu;
            k++;j++;i++;
         }
         else {
            if(bx-j>sx-i)w=(bx-j)/(sx-i);
            else w=1;
            while((j+w<bx)&&(su>*(bdx+j+w)))j+=w;
            if(j+w>=bx){
               w=bx-j-1;
               if(su>*(bdx+j+w))i=sx;
            }
            if(i<sx){
               if(su==*(bdx+j+w)){
                  *(pdx+k)=su;
                  k++;i++;j+=w+1;
               }
               else {
                  p=j+w;
                  while(p-j>1){
                     m=(j+p)/2;
                     if(su<*(bdx+m))p=m;
                     else j=m;
                  }
                  if(su==*(bdx+j)){
                     *(pdx+k)=su;
                     k++;
                  }
                  i++;j++;
               }
            }
         }
      }
   }
   if(k==0){
      delete [] pdx;
      return(NULL);
   }
   else {
      Indx<Z> *pnd=new Indx<Z>(k);
      for(i=0;i<k;i++)*(pnd->idx+i)=*(pdx+i);
      delete [] pdx;
      return(pnd);
   }
}

template<class Z>
Indx<Z> *Indx<Z>::cbool_Or(Indx<Z> *jnd){
   Z i,j,k,bx,ii,jj;
   Z *iix,*jjx,iu,ju,*pdx;

   if(jnd==NULL){
      Indx<Z> *pnd=new Indx<Z>(ix,idx,0);
      return(pnd);
   }

   ii=ix;
   iix=idx;
   jj=jnd->ix;
   jjx=jnd->idx;
   bx=ii+jj;
   pdx=new Z[bx];
   i=j=k=0;
   while((i<ii)&&(j<jj)){
      ju=*(jjx+j);
      while((i<ii)&&((iu=*(iix+i))<ju)){
         *(pdx+k)=iu;
         k++;i++;
      }
      if(i<ii){
         if(iu==ju){
            *(pdx+k)=iu;
            k++;i++;j++;
         }
         else {
            while((j<jj)&&(iu>(ju=*(jjx+j)))){
               *(pdx+k)=ju;
               k++;j++;
            }
            if(j<jj){
               if(iu==ju){
                  *(pdx+k)=iu;
                  k++;i++;j++;
               }
            }
         }
      }
   }
   while(i<ii){
      *(pdx+k)=*(iix+i);
      k++;i++;
   }
   while(j<jj){
      *(pdx+k)=*(jjx+j);
      k++;j++;
   }

   if(k==0){
      delete [] pdx;
      return(NULL);
   }
   else {
      Indx<Z> *pnd=new Indx<Z>(k);
      for(i=0;i<k;i++)*(pnd->idx+i)=*(pdx+i);
      delete [] pdx;
      return(pnd);
   }
}

template<class Z>
Indx<Z> *Indx<Z>::cbool_Butnot(Indx<Z> *jnd){
   Z i,j,k,w,p,m,flab,bx,sx,*pdx;
   Z *bdx,*sdx,bu,su;

   if(jnd==NULL){
      Indx<Z> *pnd=new Indx<Z>(ix,idx,0);
      return(pnd);
   }

   if(ix>jnd->ix){
      bdx=idx;
      bx=ix;
      sdx=jnd->idx;
      sx=jnd->ix;
      flab=1;
   }
   else {
      bdx=jnd->idx;
      bx=jnd->ix;
      sdx=idx;
      sx=ix;
      flab=0;
   }
   pdx=new Z[ix];
   for(i=0;i<ix;i++)*(pdx+i)=1; //Initialize as marker.
   if(flab){ //Case ind is big.
   i=j=k=0;
   while((i<sx)&&(j<bx)){
      bu=*(bdx+j);
      while((i<sx)&&(*(sdx+i)<bu))i++;
      if(i<sx){
         su=*(sdx+i);
         if(su==bu){
            *(pdx+j)=0;
            k++;j++;i++;
         }
         else {
            if(bx-j>sx-i)w=(bx-j)/(sx-i);
            else w=1;
            while((j+w<bx)&&(su>*(bdx+j+w)))j+=w;
            if(j+w>=bx){
               w=bx-j-1;
               if(su>*(bdx+j+w))i=sx;
            }
            if(i<sx){
               if(su==*(bdx+j+w)){
                  *(pdx+j+w)=0;
                  k++;i++;j+=w+1;
               }
               else {
                  p=j+w;
                  while(p-j>1){
                     m=(j+p)/2;
                     if(su<*(bdx+m))p=m;
                     else j=m;
                  }
                  if(su==*(bdx+j)){
                     *(pdx+j)=0;
                     k++;
                  }
                  i++;j++;
               }
            }
         }
      }
   }
   } //End of case ind is big.
   else { //Case ind is small.
   i=j=k=0;
   while((i<sx)&&(j<bx)){
      bu=*(bdx+j);
      while((i<sx)&&(*(sdx+i)<bu))i++;
      if(i<sx){
         su=*(sdx+i);
         if(su==bu){
            *(pdx+i)=0;
            k++;j++;i++;
         }
         else {
            if(bx-j>sx-i)w=(bx-j)/(sx-i);
            else w=1;
            while((j+w<bx)&&(su>*(bdx+j+w)))j+=w;
            if(j+w>=bx){
               w=bx-j-1;
               if(su>*(bdx+j+w))i=sx;
            }
            if(i<sx){
               if(su==*(bdx+j+w)){
                  *(pdx+i)=0;
                  k++;i++;j+=w+1;
               }
               else {
                  p=j+w;
                  while(p-j>1){
                     m=(j+p)/2;
                     if(su<*(bdx+m))p=m;
                     else j=m;
                  }
                  if(su==*(bdx+j)){
                     *(pdx+i)=0;
                     k++;
                  }
                  i++;j++;
               }
            }
         }
      }
   }
   } //End of case ind is small.

   j=ix-k;
   if(k==0){
      delete [] pdx;
      Indx<Z> *pnd=new Indx<Z>(ix,idx,0);
      return(pnd);
   }
   else if(j==0){
      delete [] pdx;
      return(NULL);
   }
   else {
      Indx<Z> *pnd=new Indx<Z>(j);
      j=0;
      for(i=0;i<ix;i++){
         if(*(pdx+i)){
            *(pnd->idx+j)=*(idx+i);
            j++;
         }
      }
      delete [] pdx;
      return(pnd);
   }
}

template<class Z>
Z Indx<Z>::CAnd(Indx<Z> *jnd){
   Z i,j,k,su,w,p,m,bu;
   Z *bdx,*sdx,bx,sx;

   if(jnd==NULL){
      return(0);
   }
   if(ix>jnd->ix){
      bdx=idx;
      bx=ix;
      sdx=jnd->idx;
      sx=jnd->ix;
   }
   else {
      bdx=jnd->idx;
      bx=jnd->ix;
      sdx=idx;
      sx=ix;
   }
   i=j=k=0;
   while((i<sx)&&(j<bx)){
      bu=*(bdx+j);
      while((i<sx)&&(*(sdx+i)<bu))i++;
      if(i<sx){
         su=*(sdx+i);
         if(su==bu){
            k++;j++;i++;
         }
         else {
            if(bx-j>sx-i)w=(bx-j)/(sx-i);
            else w=1;
            while((j+w<bx)&&(su>*(bdx+j+w)))j+=w;
            if(j+w>=bx){
               w=bx-j-1;
               if(su>*(bdx+j+w))i=sx;
            }
            if(i<sx){
               if(su==*(bdx+j+w)){
                  k++;i++;j+=w+1;
               }
               else {
                  p=j+w;
                  while(p-j>1){
                     m=(j+p)/2;
                     if(su<*(bdx+m))p=m;
                     else j=m;
                  }
                  if(su==*(bdx+j)){
                     k++;
                  }
                  i++;j++;
               }
            }
         }
      }
   }
   return(k);
}

template<class Z>
void Indx<Z>::Compose(Indx<Z> *pInd){
   Z i;
   
   for(i=0;i<ix;i++){
      idx[i]=pInd->idx[idx[i]];
   }
}

template<class Z>
Z Indx<Z>::Subvalue(Z j){
   Z x,y,m;
   Z i,k;

   if(j<=(k=idx[0])){
      if(j!=k)return(0);
      else return(1);
   }
   if(j>=(k=idx[ix-1])){
      if(j!=k)return(0);
      else return(ix);
   }
   x=0;
   y=ix-1;
   if(y==1)return(0);

   while(y-x>1){
      m=(y+x)/2;
      if(j>(k=idx[m]))x=m;
      else if(j<k)y=m;
      else return(m+1);
   }
   return(0);
}

template<class Z>
Indx<Z> *Indx<Z>::Subvalue(Indx *jnd){
   if(jnd==NULL)return(NULL);
   else if(jnd->ix==0)return(NULL);

   Indx<Z> *knd=new Indx<Z>(jnd->ix);

   Z i=0,j=0;
   Z k;
   while(j<jnd->ix){
      k=jnd->idx[j];
      while(idx[i]<k)i++;
      knd->idx[(j++)]=i;
   }

   return(knd);
}

template<class Z>
Indx<Z> *Indx<Z>::Subinterval(Z n,Z m){
   if(m<=n)return(NULL);
   if(n<0)return(NULL);
   if(ix<m)return(NULL);
   Indx<Z> *pind=new Indx<Z>(m-n);
   Z i,*ptr;
   ptr=pind->idx;
   for(i=n;i<m;i++)*(ptr++)=*(idx+i);
   return(pind);
}

template<class Z>
Indx<Z> *Indx<Z>::Subsample(Z n,long seed){
   Z i,j,*ptr;

   if(n>=ix){
      Indx<Z> *pnd=new Indx<Z>(ix,idx,0);
      return(pnd);
   }

   srandom((unsigned int)seed);

   Z *udx=new Z[ix];
   ptr=udx;
   for(i=0;i<n;i++)*(ptr++)=1;
   for(i=n;i<ix;i++)*(ptr++)=0;
   xshuffle(ix,udx);
   Indx<Z> *rind=new Indx<Z>(n);
   j=0;
   ptr=udx;
   for(i=0;i<ix;i++){
      if(*(ptr++))*(rind->idx+(j++))=*(idx+i);
   }
   delete [] udx;
   return(rind);
}

template<class Z>
Indx<Z> *Indx<Z>::Subsample(Z n){
   Z i,j,*ptr;

   if(n>=ix){
      Indx<Z> *pnd=new Indx<Z>(ix,idx,0);
      return(pnd);
   }

   Z *udx=new Z[ix];
   ptr=udx;
   for(i=0;i<n;i++)*(ptr++)=1;
   for(i=n;i<ix;i++)*(ptr++)=0;
   xshuffle(ix,udx);
   Indx<Z> *rind=new Indx<Z>(n);
   j=0;
   ptr=udx;
   for(i=0;i<ix;i++){
      if(*(ptr++))*(rind->idx+(j++))=*(idx+i);
   }
   delete [] udx;
   return(rind);
}

template<class Z>
Z Indx<Z>::Lbound(Z n){
   if(idx[0]>=n)return(0);
   if(idx[ix-1]<n)return(ix);
   Z i=0,j=ix-1;
   Z k;

   while(j-i>1){
      k=(i+j)/2;
      if(idx[k]<n)i=k;
      else j=k;
   }
   return(j);
}

template<class Z>
Z Indx<Z>::Ubound(Z n){
   if(idx[ix-1]<n)return(ix-1);
   if(idx[0]>=n)return(-1);
   Z i=0,j=ix-1;
   Z k;

   while(j-i>1){
      k=(i+j)/2;
      if(idx[k]<n)i=k;
      else j=k;
   }
   return(i);
}

template<class Z> template<typename U>
Indx<Z> *Indx<Z>::Greater(U *sxx,U thresh){
   Z i,k,ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)>thresh)ct++;
   }
   if(!ct)return(NULL);
   Indx<Z> *ind=new Indx<Z>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)>thresh)ind->idx[ct++]=idx[i];
   }
   return(ind);
}

template<class Z> template<typename U>
Indx<Z> *Indx<Z>::Greateq(U *sxx,U thresh){
   Z i,k,ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)>=thresh)ct++;
   }
   if(!ct)return(NULL);
   Indx<Z> *ind=new Indx<Z>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)>=thresh)ind->idx[ct++]=idx[i];
   }
   return(ind);
}

template<class Z> template<typename U>
Indx<Z> *Indx<Z>::Lesser(U *sxx,U thresh){
   Z i,k,ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)<thresh)ct++;
   }
   if(!ct)return(NULL);
   Indx<Z> *ind=new Indx<Z>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)<thresh)ind->idx[ct++]=idx[i];
   }
   return(ind);
}

template<class Z> template<typename U>
Indx<Z> *Indx<Z>::Lesseq(U *sxx,U thresh){
   Z i,k,ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)<=thresh)ct++;
   }
   if(!ct)return(NULL);
   Indx<Z> *ind=new Indx<Z>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+i)<=thresh)ind->idx[ct++]=idx[i];
   }
   return(ind);
}

//sxx[idx[]]
template<class Z> template<typename U>
Indx<Z> *Indx<Z>::Greater(U thresh,U *sxx){
   Z i,k,ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+idx[i])>thresh)ct++;
   }
   if(!ct)return(NULL);
   Indx<Z> *ind=new Indx<Z>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+(k=idx[i]))>thresh)ind->idx[ct++]=k;
   }
   return(ind);
}

template<class Z> template<typename U>
Indx<Z> *Indx<Z>::Greateq(U thresh,U *sxx){
   Z i,k,ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+idx[i])>=thresh)ct++;
   }
   if(!ct)return(NULL);
   Indx<Z> *ind=new Indx<Z>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if(*(sxx+(k=idx[i]))>=thresh)ind->idx[ct++]=k;
   }
   return(ind);
}

template<class Z> template<typename U>
Ordr<Z,U> *Indx<Z>::oGreater(U thresh,U *sxx){
   Z i,k,ct=0;
   U ss;

   for(i=0;i<ix;i++){
      if(*(sxx+idx[i])>thresh)ct++;
   }
   if(!ct)return(NULL);
   Ordr<Z,U> *jOrd=new Ordr<Z,U>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if((ss=*(sxx+(k=idx[i])))>thresh){
         jOrd->score[ct]=ss;
         jOrd->order[ct]=ct;
         jOrd->idx[ct++]=k;
      }
   }
   if(ct>1)hRort(ct,jOrd->score,jOrd->order);
   return(jOrd);
}

template<class Z> template<typename U>
Ordr<Z,U> *Indx<Z>::oGreateq(U thresh,U *sxx){
   Z i,k,ct=0;
   U ss;

   for(i=0;i<ix;i++){
      if(*(sxx+idx[i])>=thresh)ct++;
   }
   if(!ct)return(NULL);
   Ordr<Z,U> *jOrd=new Ordr<Z,U>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if((ss=*(sxx+(k=idx[i])))>=thresh){
         jOrd->score[ct]=ss;
         jOrd->order[ct]=ct;
         jOrd->idx[ct++]=k;
      }
   }
   if(ct>1)hRort(ct,jOrd->score,jOrd->order);
   return(jOrd);
}

template<class Z> template<typename U>
Ordr<Z,U> *Indx<Z>::oLesser(U thresh,U *sxx){
   Z i,k,ct=0;
   U ss;

   for(i=0;i<ix;i++){
      if(*(sxx+idx[i])<thresh)ct++;
   }
   if(!ct)return(NULL);
   Ordr<Z,U> *jOrd=new Ordr<Z,U>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if((ss=*(sxx+(k=idx[i])))<thresh){
         jOrd->score[ct]=ss;
         jOrd->order[ct]=ct;
         jOrd->idx[ct++]=k;
      }
   }
   if(ct>1)hRort(ct,jOrd->score,jOrd->order);
   return(jOrd);
}

template<class Z> template<typename U>
Ordr<Z,U> *Indx<Z>::oLesseq(U thresh,U *sxx){
   Z i,k,ct=0;
   U ss;

   for(i=0;i<ix;i++){
      if(*(sxx+idx[i])<=thresh)ct++;
   }
   if(!ct)return(NULL);
   Ordr<Z,U> *jOrd=new Ordr<Z,U>(ct);
   ct=0;
   for(i=0;i<ix;i++){
      if((ss=*(sxx+(k=idx[i])))<=thresh){
         jOrd->score[ct]=ss;
         jOrd->order[ct]=ct;
         jOrd->idx[ct++]=k;
      }
   }
   if(ct>1)hRort(ct,jOrd->score,jOrd->order);
   return(jOrd);
}

template<class Z>
void Indx<Z>::Save(const char *nam){
   FBase Fb("index",nam);
   ofstream *pfout=Fb.get_Ostr("x",ios::out);
   write(*pfout);
   Fb.dst_Ostr(pfout);
}

template<class Z>
void Indx<Z>::Load(const char *nam){
   FBase Fb("index",nam);
   ifstream *pfin=Fb.get_Istr("x",ios::in);
   read(*pfin);
   Fb.dst_Istr(pfin);
}

template<class Z>
void Indx<Z>::Save(const char *nam,const char *p_nam){
   FBase Fb("index",nam,p_nam);
   ofstream *pfout=Fb.get_Ostr("x",ios::out);
   write(*pfout);
   Fb.dst_Ostr(pfout);
}

template<class Z>
void Indx<Z>::Load(const char *nam,const char *p_nam){
   FBase Fb("index",nam,p_nam);
   ifstream *pfin=Fb.get_Istr("x",ios::in);
   read(*pfin);
   Fb.dst_Istr(pfin);
}

template<class Z>
void Indx<Z>::write(ofstream &fout){
   Z i;

   fout << ix << endl;
   for(i=0;i<ix;i++){
      fout << " " << *(idx+i) << endl;
   }
}

template<class Z>
void Indx<Z>::read(ifstream &fin){
   Z i;

   if(idx!=NULL)delete [] idx;

   fin >> ix;
   idx=new Z[ix];
   for(i=0;i<ix;i++){
      fin >> idx[i];
   }
}

template<class Z>
void Indx<Z>::debug(void){
   Z i;

   for(i=0;i<ix;i++)cout << i << " " << idx[i] << endl;
}

//Template Order class
template<class Y,class Z>
class Ordr : public Indx<Y> {
   public:
      Ordr(void); //Simple constructor.
      Ordr(Y n); //Constructer creates arrays
      Ordr(Y n,Y m,Z *sco); //sco[i] is score for index i.
         //m is database & sco array size. Produces an Order of size n.
      Ordr(Y n,Indx<Y> *ind,Z *sco); //sco[ind[i]] is interpretation.
         //Gives the top n.
      Ordr(Y n,Z *sco,Indx<Y> *ind); //idx[i] is index of
         //sco[i]. Produces an Order size n.
      Ordr(map<Y,Z> Mp);  //Produces from map Mp
      Ordr(Ordr<Y,Z> *pord); //copy constructor
      ~Ordr(void); //Removes data.
      Ordr<Y,Z> *cbool_And(Indx<Y> *jnd); //Ands with the argument but produces
         //an Order as output with the same scores and order restricted
         //to the subset.
      Ordr<Y,Z> *cbool_Butnot(Indx<Y> *jnd); //Removes the argument set and
         //produces an Order as output with the restricted parental
         //scores and order.
      Indx<Y> *Topn(Y n); //Returns an index object that contains the n
         //top scoring elements if n<=pInd->ix, else returns copy of pInd
      Z Precision(Y n,Indx<Y> *ind); //Assumes that ind is a subset
         //of the index object of this order and returns the precision over
         //the top n scores. Deals with problem of identical scores.
      Y CtGreateq(Z thresh); //Returns the number of elements with
         //scores greater than or equal to thresh.
      Z SVar(Ordr<Y,Z> *qOrd); //Produces a score which is ratio of sum
         //of squares of scores of qOrd objects in this to sum of squares
         //of scores in all of qOrd
      Z SVar(Ordr<Y,Z> *qOrd,Y m); //Produces a score which is ratio of sum
         //of squares of scores of top m qOrd objects in this to sum of squares
         //of scores of same in qOrd
      Z Proj(void); //Returns sum of squares of scores.
      Order *Convert(void); //Converts to standard Order object 
      void Save_O(const char *nam,const char *p_nam); //Saves the current object 
         //in "index"_nam.x file and uses p_nam as name for setting path
      void Load_O(const char *nam,const char *p_nam); //Loads the object 
         //"index"_nam.x into current object and uses p_nam as name for setting path
      void write_o(ofstream &fout); //Writes the object into a file
      void read_o(ifstream &fin); //Reads the object in from a file
      void debug(void); //Prints data to standard out.
      void debugW(void); //Prints data to standard out.

      //Data Access
      Y ind(Y i,Z &sco); //Returns doc number (index) and
         //score pair that correlate and are in decreasing score order as i
         //increases from 0 to num-1.
      Y inv(Y i); //Returns the number j with i=order[j], else pInd->ix
      Y *seq(void); //Returns pointer at new array that contains the array
         //pInd->idx in decreasing score order. 
      Y *invert(void); //Returns pointer at new array ax that satisfies 
         //ax[order[i]]=i. Thus order[ax[j]]=j and idx[j]=ind(ax[j],sx).

      Y *order; //Holds the order of decreasing score.
      Z *score; //Holds the scores in the ord order.
};

template<class Y,class Z>
Ordr<Y,Z>::Ordr(void) : Indx<Y>(){
   order=NULL;
   score=NULL;
}

template<class Y,class Z>
Ordr<Y,Z>::Ordr(Y n) : Indx<Y>(n){
   order=new Y[n];
   score=new Z[n];
}

template<class Y,class Z>
Ordr<Y,Z>::~Ordr(void){
   if(order)delete [] order;
   if(score)delete [] score;
}

template<class Y,class Z>
Ordr<Y,Z>::Ordr(Y n,Y m,Z *sco) : Indx<Y>(){
   Y i,j,k,*pt,ir,ii;
   Y u,*ord,*inv;
   Z *scx,ss,xx,*bt,*dt;
   n=(m<n)?m:n;

   if(n<2){
      if(n<1){
         this->idx=NULL;
         order=NULL;
         score=NULL;
         return;
      }
      ss=*sco;
      i=1;
      j=0;
      while(i<m){
         xx=sco[i];
         if(ss<xx){
            j=i;
            ss=xx;
         }
         i++;
      }
      this->idx=new Y[1];
      *this->idx=j;
      this->ix=n;
      order=new Y[1];
      score=new Z[1];
      *order=0;
      *score=ss;
      return;
   }

   this->idx=new Y[n];
   scx=new Z[n];

   pt=this->idx;
   bt=scx;
   dt=sco;
   for(i=0;i<n;i++){
      *(pt++)=i;
      *(bt++)=*(dt++);
   }

   //Build the initial heap
   k=(n>>1);
   ir=n-1;
   while(k){
      ss=scx[(--k)];
      ii=this->idx[k];

      i=k;
      j=((k+1)<<1)-1;
      while(j<=ir){
         if(j<ir && scx[j]>scx[j+1])++j;
         if(ss>scx[j]){
            scx[i]=scx[j];
            this->idx[i]=this->idx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      scx[i]=ss;
      this->idx[i]=ii;
   }
   //Filter the remaining points into heap
   xx=*scx;
   for(u=n;u<m;u++){
      if((ss=*(dt++))>xx){
         ii=u;
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && scx[j]>scx[j+1])++j;
            if(ss>scx[j]){
               scx[i]=scx[j];
               this->idx[i]=this->idx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         scx[i]=ss;
         this->idx[i]=ii;
         xx=*scx;
      }
   }
   //Order the heap by promotion & filtering
   for(;;){
      ss=scx[ir];
      ii=this->idx[ir];
      scx[ir]=scx[0];
      this->idx[ir]=this->idx[0];
      if((--ir)==0){
         scx[0]=ss;
         this->idx[0]=ii;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && scx[j]>scx[j+1])++j;
         if(ss>scx[j]){
            scx[i]=scx[j];
            this->idx[i]=this->idx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      scx[i]=ss;
      this->idx[i]=ii;
   }

   ord=new Y[n];
   inv=new Y[n];
   pt=inv;
   for(i=0;i<n;i++)*(pt++)=i;
   hSort(n,this->idx,inv);
   pt=inv;
   for(i=0;i<n;i++)*(ord+*(pt++))=i;
   delete [] inv;

   this->ix=n;
   order=ord;
   score=scx;
}

template<class Y,class Z>
Ordr<Y,Z>::Ordr(Y n,Indx<Y> *ind,Z *sco) : Indx<Y>(){
   Y i,j,k,*pt,ir,ii;
   Y u,*ord,*inv;
   Z *scx,ss,xx,*bt,*dt;
   Y m=ind->ix,*udx=ind->idx;
   n=(m<n)?m:n;

   if(n<2){
      if(n<1){
         this->idx=NULL;
         order=NULL;
         score=NULL;
         return;
      }
      ss=sco[udx[0]];
      i=1;
      j=udx[0];
      while(i<m){
         xx=sco[udx[i]];
         if(ss<xx){
            j=udx[i];
            ss=xx;
         }
         i++;
      }
      this->idx=new Y[1];
      *this->idx=j;
      this->ix=n;
      order=new Y[1];
      score=new Z[1];
      *order=0;
      *score=ss;
      return;
   }

   this->idx=new Y[n];
   scx=new Z[n];

   pt=this->idx;
   bt=scx;
   for(i=0;i<n;i++){
      *(pt++)=udx[i];
      *(bt++)=sco[udx[i]];
   }

   //Build the initial heap
   k=(n>>1);
   ir=n-1;
   while(k){
      ss=scx[(--k)];
      ii=this->idx[k];

      i=k;
      j=((k+1)<<1)-1;
      while(j<=ir){
         if(j<ir && scx[j]>scx[j+1])++j;
         if(ss>scx[j]){
            scx[i]=scx[j];
            this->idx[i]=this->idx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      scx[i]=ss;
      this->idx[i]=ii;
   }
   //Filter the remaining points into heap
   xx=*scx;
   for(u=n;u<m;u++){
      if((ss=sco[udx[u]])>xx){
         ii=udx[u];
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && scx[j]>scx[j+1])++j;
            if(ss>scx[j]){
               scx[i]=scx[j];
               this->idx[i]=this->idx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         scx[i]=ss;
         this->idx[i]=ii;
         xx=*scx;
      }
   }
   //Order the heap by promotion & filtering
   for(;;){
      ss=scx[ir];
      ii=this->idx[ir];
      scx[ir]=scx[0];
      this->idx[ir]=this->idx[0];
      if((--ir)==0){
         scx[0]=ss;
         this->idx[0]=ii;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && scx[j]>scx[j+1])++j;
         if(ss>scx[j]){
            scx[i]=scx[j];
            this->idx[i]=this->idx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      scx[i]=ss;
      this->idx[i]=ii;
   }

   ord=new Y[n];
   inv=new Y[n];
   pt=inv;
   for(i=0;i<n;i++)*(pt++)=i;
   hSort(n,this->idx,inv);
   pt=inv;
   for(i=0;i<n;i++)*(ord+*(pt++))=i;
   delete [] inv;

   this->ix=n;
   order=ord;
   score=scx;
}

template<class Y,class Z>
Ordr<Y,Z>::Ordr(Y n,Z *sco,Indx<Y>*ind) : Indx<Y>(){
   Y i,j,k,*pt,*qt,ir,ii;
   Y u,*ord,*inv;
   Z *scx,ss,xx,*bt,*dt;
   n=((ind->ix)<n)?(ind->ix):n;

   if(n<2){
      if(n<1){
         this->idx=NULL;
         order=NULL;
         score=NULL;
         return;
      }
      ss=*sco;
      i=1;
      j=0;
      while(i<ind->ix){
         xx=sco[i];
         if(ss<xx){
            j=i;
            ss=xx;
         }
         i++;
      }
      this->idx=new Y[1];
      *this->idx=ind->idx[0];
      this->ix=n;
      order=new Y[1];
      score=new Z[1];
      *order=0;
      *score=ss;
      return;
   }

   this->idx=new Y[n];
   scx=new Z[n];

   pt=this->idx;
   qt=ind->idx;
   bt=scx;
   dt=sco;
   for(i=0;i<n;i++){
      *(pt++)=*(qt++);
      *(bt++)=*(dt++);
   }

   //Build the initial heap
   k=(n>>1);
   ir=n-1;
   while(k){
      ss=scx[(--k)];
      ii=this->idx[k];

      i=k;
      j=((k+1)<<1)-1;
      while(j<=ir){
         if(j<ir && scx[j]>scx[j+1])++j;
         if(ss>scx[j]){
            scx[i]=scx[j];
            this->idx[i]=this->idx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      scx[i]=ss;
      this->idx[i]=ii;
   }
   //Filter the remaining points into heap
   xx=*scx;
   for(u=n;u<ind->ix;u++){
      if((ss=*(dt++))>xx){
         ii=ind->idx[u];
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && scx[j]>scx[j+1])++j;
            if(ss>scx[j]){
               scx[i]=scx[j];
               this->idx[i]=this->idx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         scx[i]=ss;
         this->idx[i]=ii;
         xx=*scx;
      }
   }
   //Order the heap by promotion & filtering
   for(;;){
      ss=scx[ir];
      ii=this->idx[ir];
      scx[ir]=scx[0];
      this->idx[ir]=this->idx[0];
      if((--ir)==0){
         scx[0]=ss;
         this->idx[0]=ii;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && scx[j]>scx[j+1])++j;
         if(ss>scx[j]){
            scx[i]=scx[j];
            this->idx[i]=this->idx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      scx[i]=ss;
      this->idx[i]=ii;
   }

   ord=new Y[n];
   inv=new Y[n];
   pt=inv;
   for(i=0;i<n;i++)*(pt++)=i;
   hSort(n,this->idx,inv);
   pt=inv;
   for(i=0;i<n;i++)*(ord+*(pt++))=i;
   delete [] inv;

   this->ix=n;
   order=ord;
   score=scx;
}

template<class Y,class Z>
Ordr<Y,Z>::Ordr(Ordr<Y,Z> *pOrd) : Indx<Y>(){
   Y i;

   this->ix=pOrd->ix;
   this->idx=new Y[pOrd->ix];
   order=new Y[pOrd->ix];
   score=new Z[pOrd->ix];
   for(i=0;i<pOrd->ix;i++){
      this->idx[i]=pOrd->idx[i];
      order[i]=pOrd->order[i];
      score[i]=pOrd->score[i];
   }
}

template<class Y,class Z>
Ordr<Y,Z>::Ordr(map<Y,Z> Mp) : Indx<Y>(){
   long i;

   typename map<Y,Z>::iterator bg,en;
   i=Mp.size();
   Y *nx=new Y[i];
   Z *sx=new Z[i];
   bg=Mp.begin();
   en=Mp.end();
   i=0;
   while(bg!=en){
      nx[i]=bg->first;
      sx[i]=bg->second;
      i++;bg++;
   }
   hSort(i,nx,sx);
   Indx<Y> *pInd=new Indx<Y>();
   pInd->ix=i;
   pInd->idx=nx;
   Ordr<Y,Z> *pOrd=new Ordr<Y,Z>(i,sx,pInd);
   this->ix=pOrd->ix;
   this->idx=pOrd->idx;
   order=pOrd->order;
   score=pOrd->score;
   delete [] nx;
   delete [] sx;
   pInd->idx=NULL;
   delete pInd;
   pOrd->idx=NULL;
   pOrd->order=NULL;
   pOrd->score=NULL;
   delete pOrd;
}

template<class Y,class Z>
Y Ordr<Y,Z>::ind(Y i,Z &sco){
   sco=*(score+i);
   return(*(this->idx+*(order+i)));
}

template<class Y,class Z>
Y Ordr<Y,Z>::inv(Y i){
   Y j=0;
   while((j<this->ix)&&(i!=order[j]))j++;
   return(j);
}

template<class Y,class Z>
Y *Ordr<Y,Z>::seq(void){
   Y *ssq=new Y[this->ix];
   Y i;
   for(i=0;i<this->ix;i++)ssq[i]=*(this->idx+*(order+i));
   return(ssq);
}

template<class Y,class Z>
Y *Ordr<Y,Z>::invert(void){
   Y *ssq=new Y[this->ix];
   Y i;
   for(i=0;i<this->ix;i++)ssq[*(order+i)]=i;
   return(ssq);
}

template<class Y,class Z>
Ordr<Y,Z> *Ordr<Y,Z>::cbool_And(Indx<Y> *jnd){
   Y i,j,k,*sub;

   Indx<Y> *pUnd=(Indx<Y>*)this;
   Indx<Y> *pind=pUnd->cbool_And(jnd);
   if(pind==NULL)return(NULL);
   if(pind->ix<1){delete pind;return(NULL);}

   sub=new Y[this->ix];
   for(i=0;i<this->ix;i++)sub[i]=0;
   Y *pi=this->idx;
   Y *pj=pind->idx;
   i=j=0;
   while(j<pind->ix){
      while(pi[i]<pj[j])i++;
      sub[i]=j+1;
      j++;
   }
   Ordr<Y,Z> *psub=new Ordr<Y,Z>;
   psub->idx=pind->idx;
   psub->ix=pind->ix;
   psub->order=new Y[pind->ix];
   psub->score=new Z[pind->ix];
   j=0;
   for(i=0;i<this->ix;i++){
      if((k=sub[order[i]])>0){
         psub->order[j]=k-1;
         psub->score[j]=score[i];
         j++;
      }
   }
   delete [] sub;
   pind->idx=NULL;
   delete pind;
   return(psub);
}

template<class Y,class Z>
Ordr<Y,Z> *Ordr<Y,Z>::cbool_Butnot(Indx<Y> *jnd){
   Y i,j,k,*sub;

   Indx<Y> *pUnd=(Indx<Y>*)this;
   Indx<Y> *pind=pUnd->cbool_Butnot(jnd);
   if(pind==NULL)return(NULL);
   if(pind->ix<1)return(NULL);

   sub=new Y[this->ix];
   for(i=0;i<this->ix;i++)sub[i]=0;
   Y *pi=this->idx;
   Y *pj=pind->idx;
   i=j=0;
   while(j<pind->ix){
      while(pi[i]<pj[j])i++;
      sub[i]=j+1;
      j++;
   }
   Ordr<Y,Z> *psub=new Ordr<Y,Z>;
   psub->idx=pind->idx;
   psub->ix=pind->ix;
   psub->order=new Y[pind->ix];
   psub->score=new Z[pind->ix];
   j=0;
   for(i=0;i<this->ix;i++){
      if((k=sub[order[i]])>0){
         psub->order[j]=k-1;
         psub->score[j]=score[i];
         j++;
      }
   }
   delete [] sub;
   return(psub);
}

template<class Y,class Z>
Z Ordr<Y,Z>::Precision(Y n,Indx<Y> *ind){
   if(n>this->ix)return(0);
   else if(!n)return(0);
   Indx<Y> *jnd=this->Subvalue(ind);
   Z cx=0.0,sx,ss,vx;
   Y i,j,k=0;
   Y *bz=new Y[this->ix];
   for(i=0;i<this->ix;i++)bz[i]=0;
   for(i=0;i<jnd->ix;i++)bz[jnd->idx[i]]=1;
   sx=score[0];
   while(k<n){
      i=0;
      vx=0.0;
      while((k+i<this->ix)&&(score[k+i]==sx)){
         vx+=bz[order[k+i]];
         i++;
      }
      if(k+i<n){
         cx+=vx;
         k+=i;
         sx=score[k];
      }
      else {
         cx+=vx*(n-k)/((Z)i);
         k+=i;
      }
   }
   delete jnd;
   delete [] bz;
   return(cx/((Z)n));
}

template<class Y,class Z>
Y Ordr<Y,Z>::CtGreateq(Z thresh){
   Z sss,si,sj;
   Y i,j,k;

   i=0;
   ind(i,sss);
   if(sss<thresh)return(0);
   j=this->ix-1;
   ind(j,sss);
   if(sss>=thresh)return(j+1);
   while(j-i>1){
      k=(i+j)/2;
      ind(k,sss);
      if(sss>=thresh)i=k;
      else j=k;
   }
   return(i+1);
}

template<class Y,class Z>
void Ordr<Y,Z>::debug(void){
   Y i,j,k;
   Z ss;
   cout << this->ix << endl;
   for(i=0;i<this->ix;i++){
      j=ind(i,ss);
      cout << i << " " << j << " " << ss << endl;
   }
}

template<class Y,class Z>
void Ordr<Y,Z>::debugW(void){
   Y i;

   cout << this->ix << endl;
   for(i=0;i<this->ix;i++){
      cout << " " << *(this->idx+i) << " " << *(order+i) << " " \
         << *(score+i) << endl;
   }
}

template<class Y,class Z>
Order *Ordr<Y,Z>::Convert(void){
   long i;
   Index *pnd=new Index((long)this->ix);
   for(i=0;i<this->ix;i++)pnd->idx[i]=this->idx[i];
   Order *pOrd=new Order;
   pOrd->pInd=pnd;
   pOrd->order=new long[this->ix];
   pOrd->score=new float[this->ix];
   for(i=0;i<this->ix;i++){
      pOrd->order[i]=(long)order[i];
      pOrd->score[i]=(float)score[i];
   }
   return(pOrd);
}

template<class Y,class Z>
void Ordr<Y,Z>::Save_O(const char *nam,const char *p_nam){
   FBase Fb("order",nam);
   Fb.set_path_name(p_nam);
   ofstream *pfout=Fb.get_Ostr("x",ios::out);
   write_o(*pfout);
   Fb.dst_Ostr(pfout);
}

template<class Y,class Z>
void Ordr<Y,Z>::Load_O(const char *nam,const char *p_nam){
   FBase Fb("order",nam);
   Fb.set_path_name(p_nam);
   ifstream *pfin=Fb.get_Istr("x",ios::in);
   read_o(*pfin);
   Fb.dst_Istr(pfin);
}

template<class Y,class Z>
void Ordr<Y,Z>::write_o(ofstream &fout){
   Y i;

   fout << this->ix << endl;
   for(i=0;i<this->ix;i++){
      fout << " " << *(this->idx+i) << " " << *(order+i) << " " \
         << *(score+i) << endl;
   }
}

template<class Y,class Z>
void Ordr<Y,Z>::read_o(ifstream &fin){
   Y i;

   if(this->idx!=NULL)delete [] this->idx;
   if(order!=NULL)delete [] order;
   if(score!=NULL)delete [] score;

   fin >> i;
   this->idx=new Y[i];
   this->ix=i;
   order=new Y[i];
   score=new Z[i];
   for(i=0;i<this->ix;i++){
      fin >> this->idx[i] >> order[i] >> score[i];
   }
}

template<class Y,class Z>
Indx<Y> *Ordr<Y,Z>::Topn(Y n){
   if(n>=this->ix){
      Indx<Y> *pnd=new Indx<Y>(this->ix,this->idx,0);
      return(pnd);
   }
   else {
      Y i,j;
      Z xx;
      Indx<Y> *pnd=new Indx<Y>(n);
      Y *md=new Y[this->ix];
      for(j=0;j<this->ix;j++)md[j]=0;
      for(i=0;i<n;i++){
         j=order[i];
         md[j]=1;
      }
      i=0;
      for(j=0;j<this->ix;j++){
         if(md[j]){
            pnd->idx[i++]=this->idx[j];
         }
      }
      delete [] md;    
      return(pnd);
   }
}

template<class Y,class Z>
Z Ordr<Y,Z>::SVar(Ordr<Y,Z> *qOrd){
   Y i,j,k;
   Z sum1=0,sum2=0;
   Z xx,yy;
   if(!qOrd)return(0);
   for(i=0;i<this->ix;i++){
      j=this->ind(i,xx);
      if(qOrd->Subvalue(j))sum1+=xx*xx;
   }
   k=qOrd->ix;
   for(i=0;i<k;i++){
      sum2+=(qOrd->score[i])*(qOrd->score[i]);
   }
   if(sum2>0)return(sum1/sum2);
   else return(0);
}

template<class Y,class Z>
Z Ordr<Y,Z>::SVar(Ordr<Y,Z> *qOrd,Y m){
   Y i,j,k;
   Z sum1=0,sum2=0;
   Z xx,yy;
   if(!qOrd)return(0);
   Indx<Y> *sInd=new Indx<Y>(m);
   for(i=0;i<m;i++){
      j=qOrd->ind(i,xx);
      sInd->idx[i]=j;
      sum2+=xx*xx;
   }
   sInd->sSort();
   for(i=0;i<this->ix;i++){
      j=this->ind(i,xx);
      if(sInd->Subvalue(j)){
         sum1+=xx*xx;
      }
   }
   delete sInd;
   if(sum2>0)return(sum1/sum2);
   else return(0);
}

template<class Y,class Z>
Z Ordr<Y,Z>::Proj(void){
   Y i,j,k;
   Z xx,sum1=0;
   for(i=0;i<this->ix;i++){
      j=this->ind(i,xx);
      sum1+=xx*xx;
   }
   return(sum1);
}

template<class Y>
class CValdX {
public:
   CValdX(Indx<Y> *gdd, long n);
   CValdX(Indx<Y> *gdd, Indx<Y> *bdd);
   ~CValdX();
   void cross_valid(long m, long seed);
   Indx<Y> *ind;
   Indx<Y> *gnd;
   Indx<Y> *bnd;
   Indx<Y> **pGTS;
   Indx<Y> **pBTS;
   Indx<Y> **pWTS;
   Indx<Y> **pGTR;
   Indx<Y> **pBTR;
   Indx<Y> **pWTR;
   long setn;
};

template<class Y>
CValdX<Y>::CValdX(Indx<Y> *gdd, long n){
   ind=new Indx<Y>(0,n);
   gnd=new Indx<Y>(gdd);
   bnd=ind->cbool_Butnot(gnd);
}

template<class Y>
CValdX<Y>::CValdX(Indx<Y> *gdd, Indx<Y> *bdd){
   gnd=new Indx<Y>(gdd);
   bnd=new Indx<Y>(bdd);
   ind=gnd->cbool_Or(bnd);
}

template<class Y>
void CValdX<Y>::cross_valid(long m, long seed){

   long i,j,k,size,blk,rem;
   setn=m;

   pBTS=new Indx<Y>*[setn];
   pGTS=new Indx<Y>*[setn];
   pWTS=new Indx<Y>*[setn];
   pBTR=new Indx<Y>*[setn];
   pGTR=new Indx<Y>*[setn];
   pWTR=new Indx<Y>*[setn];

   size=gnd->ix;
   if(size<setn) {
      cout <<"Size of Relevant set is smaller than the the number of cross validation"<<endl;
      exit(0);
   }

   long *sizg, *sizb;

   sizg=new long[setn];
   blk=size/setn;
   rem=size%setn;

   for(i=0;i<setn;i++){
      if(i<rem){sizg[i]=blk+1;}
      else {sizg[i]=blk;}
   }

   sizb=new long[setn];
   size=bnd->ix;

   if(size<setn) {
      cout <<"Size of Non-Relevant set is smaller than the number of cross validations"<<endl;
      exit(0);
   }

   blk=size/setn;
   rem=size%setn;
   k=0;
   for(i=0;i<setn;i++){
      if(i<rem){sizb[i]=blk+1;}
      else {sizb[i]=blk;}
   }

   Indx<Y> *tmp1=new Indx<Y>(gnd);
   Indx<Y> *tmp2=new Indx<Y>(bnd);
   Indx<Y> *tmp11,*tmp22;
   for(i=0;i<setn-1;i++){

      pGTS[i]=tmp1->Subsample(sizg[i],seed);
      pBTS[i]=tmp2->Subsample(sizb[i],seed);
      pWTS[i]=pGTS[i]->cbool_Or(pBTS[i]);

      pGTR[i]=gnd->cbool_Butnot(pGTS[i]);
      pBTR[i]=bnd->cbool_Butnot(pBTS[i]);
      pWTR[i]=pGTR[i]->cbool_Or(pBTR[i]);

      tmp11=tmp1->cbool_Butnot(pGTS[i]);
      delete tmp1;
      tmp1=tmp11;
      tmp22=tmp2->cbool_Butnot(pBTS[i]);
      delete tmp2;
      tmp2=tmp22;
   }
   pGTS[i]=tmp1;
   pBTS[i]=tmp2;
   pWTS[i]=pGTS[i]->cbool_Or(pBTS[i]);

   pGTR[i]=gnd->cbool_Butnot(pGTS[i]);
   pBTR[i]=bnd->cbool_Butnot(pBTS[i]);
   pWTR[i]=pGTR[i]->cbool_Or(pBTR[i]);
   
   delete [] sizg;
   delete [] sizb;
}

template<class Y>
CValdX<Y>::~CValdX(void){
   for(long i=0;i<setn;i++){
      delete pGTS[i];
      delete pBTS[i];
      delete pWTS[i];
      delete pGTR[i];
      delete pBTR[i];
      delete pWTR[i];
   }
   delete [] pGTS;
   delete [] pBTS;
   delete [] pWTS;
   delete [] pGTR;
   delete [] pBTR;
   delete [] pWTR;  

   if (ind) delete ind;
   if (gnd) delete gnd;
   if (bnd) delete bnd;
}

template<class Y>
class Store {
public:
   Store(long n); //size of array
      //Creates array wx
   ~Store(void); 
   void CopyIn(Y *u); //Makes a copy of u
   void CopyOut(Y *u); //Copys from memory to u
   long num; //Size of array wx
   Y *wx; //Memory for storage
};

template<class Y>
Store<Y>::Store(long n){
   num=n;
   wx=new Y[n];
}

template<class Y>
Store<Y>::~Store(void){
   delete [] wx;
}

template<class Y>
void Store<Y>::CopyIn(Y *u){
   long i;
   for(i=0;i<num;i++)wx[i]=u[i];
}

template<class Y>
void Store<Y>::CopyOut(Y *u){
   long i;
   for(i=0;i<num;i++)u[i]=wx[i];
}

//This class can be used as convenient way to process a theme without any special terms
//Terms can be here by number as in Thln.h&C or by ID as in Focus.h&C
template<class Y,class Z>
class WOrdr : public Ordr<Y,Z> {
  public:
   WOrdr(void);
   WOrdr(Ordr<Y,Z> *pOrd,Z *weg); //Uses weg[pOrd->idx[order[i]]]
      //to obtain the weights for terms. ith weight corresponds to ind(i,.)
      //constructor destroys the Order object pointed
      //to by pOrd. A case of catabolism.
   ~WOrdr(void);

   Z ave_sco(Y n); //Avererages the n top scores
   Z dice_sco(WOrdr<Y,Z> *wpord); //Produces dice similarity score.
   Z sub_sco(WOrdr<Y,Z> *wpord); //Produces subset similarity score. What frac of wpord sco is in this
   int subset_equal(Y n, WOrdr *wpord); //Returns 1 if the two objects have the same terms
      //in the top n spots, order not important
   int positive(Y n); //returns 1 if first n weights are positive, else 0
   void write_w(ofstream &fout); //Writes the object into a file
   void read_w(ifstream &fin); //Reads the object in from a file
   void debug(void); //Prints data to standard out.
   //Data
   Z *weight;
};

template<class Y,class Z>
WOrdr<Y,Z>::WOrdr(void) : Ordr<Y,Z>(){
weight=NULL;
}

template<class Y,class Z>
WOrdr<Y,Z>::WOrdr(Ordr<Y,Z> *pOrd,Z *weg){
   if(!pOrd){cout << "Error, pOrd is NULL!" << endl;exit(0);}
   this->idx=pOrd->idx;
   this->ix=pOrd->ix;
   pOrd->idx=NULL;
   this->order=pOrd->order;
   pOrd->order=NULL;
   this->score=pOrd->score;
   pOrd->score=NULL;

   Y i;
   weight=new Z[this->ix];
   for(i=0;i<this->ix;i++)weight[i]=weg[this->idx[this->order[i]]];
   delete pOrd;
}

template<class Y,class Z>
WOrdr<Y,Z>::~WOrdr(void){
   if(weight!=NULL)delete [] weight;
}

template<class Y,class Z>
Z WOrdr<Y,Z>::ave_sco(Y n){
   Y i;
   Z xx,sum=0;

   for(i=0;i<n;i++){
      this->ind(i,xx);
      sum+=xx;
   }
   return(sum/((Z)n));
}

template<class Y,class Z>
Z WOrdr<Y,Z>::dice_sco(WOrdr<Y,Z> *wpord){

   Y i,j,k;
   Z s1=0,s2=0,sc1=0,sc2=0;
   Z xi,yj;
   Y m=this->ix;
   Y n=wpord->ix;

   Indx<Y> *pJnd=this->cbool_And((Indx<Y> *)wpord);
   if(!pJnd)return(0);
   for(i=0;i<m;i++){
      s1+=this->score[i];
      k=this->ind(i,xi);
      if(pJnd->Subvalue(k))sc1+=xi;
   }
   for(j=0;j<n;j++){
      s2+=wpord->score[j];
      k=wpord->ind(j,yj);
      if(pJnd->Subvalue(k))sc2+=yj;
   }
   delete pJnd;
   return((sc1+sc2)/(s1+s2));
}

template<class Y,class Z>
Z WOrdr<Y,Z>::sub_sco(WOrdr<Y,Z> *wpord){

   Y i,j,k;
   Z s1=0,s2=0,sc1=0,sc2=0;
   Z xi,yj;
   Y n=wpord->ix;

   Indx<Y> *pJnd=this->cbool_And((Indx<Y> *)wpord);
   if(!pJnd)return(0);
   for(j=0;j<n;j++){
      s2+=wpord->score[j];
      k=wpord->ind(j,yj);
      if(pJnd->Subvalue(k))sc2+=yj;
   }
   delete pJnd;
   return((sc2)/(s2));
}

template<class Y,class Z>
int WOrdr<Y,Z>::subset_equal(Y n, WOrdr<Y,Z> *wpord){

   Y i,j,k,v;
   Y flag;
   Z xi,yj;
   Y m=this->ix;
   Y r=wpord->ix;
   if((m<n)||(r<n))return(0);

   for(i=0;i<n;i++){
      k=this->ind(i,xi);
      flag=1;
      for(v=0;v<n;v++){
         if(k==(wpord->ind(v,yj))){flag=0;break;}
      }
      if(flag)return(0);
   }
   return(1);
}

template<class Y,class Z>
int WOrdr<Y,Z>::positive(Y n){
   Y i=0,flag=1;
   while(i<n){
      if(weight[i]<=0)return(0);
      i++;
   }
   return(1);
}

template<class Y,class Z>
void WOrdr<Y,Z>::write_w(ofstream &fout){
   Y i;

   fout << this->ix << endl;
   for(i=0;i<this->ix;i++){
      fout << " " << *(this->idx+i) << " " << *(this->order+i) << " " \
         << *(this->score+i) << " " << *(this->weight+i) << endl;
   }
}

template<class Y,class Z>
void WOrdr<Y,Z>::read_w(ifstream &fin){
   Y i;

   if(this->idx!=NULL)delete [] this->idx;
   if(this->order!=NULL)delete [] this->order;
   if(this->score!=NULL)delete [] this->score;
   if(this->weight!=NULL)delete [] this->weight;

   fin >> i;
   this->idx=new Y[i];
   this->ix=i;
   this->order=new Y[i];
   this->score=new Z[i];
   this->weight=new Z[i];
   for(i=0;i<this->ix;i++){
      fin >> this->idx[i] >> this->order[i] >> this->score[i] >> this->weight[i];
   }
}

template<class Y,class Z>
void WOrdr<Y,Z>::debug(void){
   Y i,j,k;
   Z ss;
   cout << this->ix << endl;
   for(i=0;i<this->ix;i++){
      j=this->ind(i,ss);
      cout << i << " " << j << " " << ss << " " << this->weight[i] << endl;
   }
}

}
#endif
