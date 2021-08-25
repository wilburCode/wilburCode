#ifndef HEAP_H
#define HEAP_H
#include <fstream>
#include <iostream>
#include <runn.h>
using namespace std;
namespace iret {

template <class X>
class Heap {
   public:
      Heap(long nm,X *ra); 
     ~Heap();
      //F stands for forward order and captures the largest values
      void add_pointF(X rx); //Adds a single point to 
      long orderF(void); //Orders the heap post additions
         //points are then in decreasing order in rx
         //Returns the number of points captured
      //R stands for reverse order and captures the smallest values
      void add_pointR(X rx); //Adds a single point to 
      long orderR(void); //Orders the heap post additions
         //points are then in increasing order in rx
         //Returns the number of points captured

      long num; //Number of points allowed
      long cnm; //Current number of points
      X *mx;
      long ir; //Used in calculations
         //is num-1
};

template <class X>
Heap<X>::Heap(long nm,X *ra){
   num=nm;
   mx=ra;
   ir=num-1;
   cnm=0;
   if(num<2){cout << "Error in heap constuctor!" << endl;exit(1);}
} 

template <class X>
Heap<X>::~Heap(){
}  
 
template <class X>
void Heap<X>::add_pointF(X xx){
   long i,j,k,n;
   X ss;

   if(num>cnm){
      mx[cnm++]=xx;
      if(cnm==num){
         //Build the initial heap
         k=(num>>1);
         while(k){
            ss=mx[(--k)];

            i=k;
            j=((k+1)<<1)-1;
            while(j<=ir){
               if(j<ir && mx[j]>mx[j+1])++j;
               if(ss>mx[j]){
                  mx[i]=mx[j];
                  j+=(i=j)+1;
               }
               else j=ir+1;
            }
            mx[i]=ss;
         }
      }
   }
   else {
      ss=*mx;
      if(xx>ss){
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && mx[j]>mx[j+1])++j;
            if(xx>mx[j]){
               mx[i]=mx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=xx;
      }
   }
}

template <class X>
long Heap<X>::orderF(void){
   long i,j,k;
   X ss;
   //If cnm<num build initial heap
   if(cnm<2)return(cnm);
   if(cnm<num){
      ir=cnm-1;
      k=(cnm>>1);
      while(k){
         ss=mx[(--k)];

         i=k;
         j=((k+1)<<1)-1;
         while(j<=ir){
            if(j<ir && mx[j]>mx[j+1])++j;
            if(ss>mx[j]){
               mx[i]=mx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=ss;
      }
   }
   //Order the heap by promotion & filtering
   ir=cnm-1;
   for(;;){
      ss=mx[ir];
      mx[ir]=mx[0];
      if((--ir)==0){
         mx[0]=ss;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && mx[j]>mx[j+1])++j;
         if(ss>mx[j]){
            mx[i]=mx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      mx[i]=ss;
   }
   return(cnm);
}
   

template <class X>
void Heap<X>::add_pointR(X xx){
   long i,j,k,n;
   X ss;

   if(num>cnm){
      mx[cnm++]=xx;
      if(cnm==num){
         //Build the initial heap
         k=(num>>1);
         while(k){
            ss=mx[(--k)];

            i=k;
            j=((k+1)<<1)-1;
            while(j<=ir){
               if(j<ir && mx[j]<mx[j+1])++j;
               if(ss<mx[j]){
                  mx[i]=mx[j];
                  j+=(i=j)+1;
               }
               else j=ir+1;
            }
            mx[i]=ss;
         }
      }
   }
   else {
      ss=*mx;
      if(xx<ss){
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && mx[j]<mx[j+1])++j;
            if(xx<mx[j]){
               mx[i]=mx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=xx;
      }
   }
}

template <class X>
long Heap<X>::orderR(void){
   long i,j,k;
   X ss;
   //If cnm<num build initial heap
   if(cnm<2)return(cnm);
   if(cnm<num){
      ir=cnm-1;
      k=(cnm>>1);
      while(k){
         ss=mx[(--k)];

         i=k;
         j=((k+1)<<1)-1;
         while(j<=ir){
            if(j<ir && mx[j]<mx[j+1])++j;
            if(ss<mx[j]){
               mx[i]=mx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=ss;
      }
   }
   //Order the heap by promotion & filtering
   ir=cnm-1;
   for(;;){
      ss=mx[ir];
      mx[ir]=mx[0];
      if((--ir)==0){
         mx[0]=ss;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && mx[j]<mx[j+1])++j;
         if(ss<mx[j]){
            mx[i]=mx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      mx[i]=ss;
   }
   return(cnm);
}

template <class X,class Y>
class Heap2 {
   public:
      Heap2(long nm,X *ra,Y *rb); 
     ~Heap2();
      //F stands for forward order and captures the largest values
      void add_pointF(X rx,Y ry); //Adds a single point to 
      long orderF(void); //Orders the heap post additions
         //points are then in decreasing order in rx
         //Returns the number of points captured
      //R stands for reverse order and captures the smallest values
      void add_pointR(X rx,Y ry); //Adds a single point to 
      long orderR(void); //Orders the heap post additions
         //points are then in increasing order in rx
         //Returns the number of points captured

      long num; //Number of points allowed
      long cnm; //Current number of points
      X *mx;
      Y *nx;
      long ir; //Used in calculations
         //is num-1
};

template <class X,class Y>
Heap2<X,Y>::Heap2(long nm,X *ra,Y *rb){
   num=nm;
   mx=ra;
   nx=rb;
   ir=num-1;
   cnm=0;
   if(num<2){cout << "Error in heap constructor!" << endl;exit(1);}
} 

template <class X,class Y>
Heap2<X,Y>::~Heap2(){
}  
 
template <class X,class Y>
void Heap2<X,Y>::add_pointF(X xx,Y yy){
   long i,j,k,n;
   X ss;
   Y tt;

   if(num>cnm){
      mx[cnm]=xx;
      nx[cnm++]=yy;
      if(cnm==num){
         //Build the initial heap
         k=(num>>1);
         while(k){
            ss=mx[(--k)];
            tt=nx[k];

            i=k;
            j=((k+1)<<1)-1;
            while(j<=ir){
               if(j<ir && mx[j]>mx[j+1])++j;
               if(ss>mx[j]){
                  mx[i]=mx[j];
                  nx[i]=nx[j];
                  j+=(i=j)+1;
               }
               else j=ir+1;
            }
            mx[i]=ss;
            nx[i]=tt;
         }
      }
   }
   else {
      ss=*mx;
      if(xx>ss){
         tt=*nx;
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && mx[j]>mx[j+1])++j;
            if(xx>mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=xx;
         nx[i]=yy;
      }
   }
}

template <class X,class Y>
long Heap2<X,Y>::orderF(void){
   long i,j,k;
   X ss;
   Y tt;
   //Build the initial heap if cnm<num
   if(cnm<2)return(cnm);
   if(cnm<num){
      ir=cnm-1;
      k=(cnm>>1);
      while(k){
         ss=mx[(--k)];
         tt=nx[k];

         i=k;
         j=((k+1)<<1)-1;
         while(j<=ir){
            if(j<ir && mx[j]>mx[j+1])++j;
            if(ss>mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=ss;
         nx[i]=tt;
      }
   }
   //Order the heap by promotion & filtering
   ir=cnm-1;
   for(;;){
      ss=mx[ir];
      tt=nx[ir];
      mx[ir]=mx[0];
      nx[ir]=nx[0];
      if((--ir)==0){
         mx[0]=ss;
         nx[0]=tt;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && mx[j]>mx[j+1])++j;
         if(ss>mx[j]){
            mx[i]=mx[j];
            nx[i]=nx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      mx[i]=ss;
      nx[i]=tt;
   }
   return(cnm);
}
   
template <class X,class Y>
void Heap2<X,Y>::add_pointR(X xx,Y yy){
   long i,j,k,n;
   X ss;
   Y tt;

   if(num>cnm){
      mx[cnm]=xx;
      nx[cnm++]=yy;
      if(cnm==num){
         //Build the initial heap
         k=(num>>1);
         while(k){
            ss=mx[(--k)];
            tt=nx[k];

            i=k;
            j=((k+1)<<1)-1;
            while(j<=ir){
               if(j<ir && mx[j]<mx[j+1])++j;
               if(ss<mx[j]){
                  mx[i]=mx[j];
                  nx[i]=nx[j];
                  j+=(i=j)+1;
               }
               else j=ir+1;
            }
            mx[i]=ss;
            nx[i]=tt;
         }
      }
   }
   else {
      ss=*mx;
      if(xx<ss){
         tt=*nx;
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && mx[j]<mx[j+1])++j;
            if(xx<mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=xx;
         nx[i]=yy;
      }
   }
}

template <class X,class Y>
long Heap2<X,Y>::orderR(void){
   long i,j,k;
   X ss;
   Y tt;
   //Build the initial heap if cnm<num
   if(cnm<2)return(cnm);
   if(cnm<num){
      ir=cnm-1;
      k=(cnm>>1);
      while(k){
         ss=mx[(--k)];
         tt=nx[k];

         i=k;
         j=((k+1)<<1)-1;
         while(j<=ir){
            if(j<ir && mx[j]<mx[j+1])++j;
            if(ss<mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=ss;
         nx[i]=tt;
      }
   }
   //Order the heap by promotion & filtering
   ir=cnm-1;
   for(;;){
      ss=mx[ir];
      tt=nx[ir];
      mx[ir]=mx[0];
      nx[ir]=nx[0];
      if((--ir)==0){
         mx[0]=ss;
         nx[0]=tt;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && mx[j]<mx[j+1])++j;
         if(ss<mx[j]){
            mx[i]=mx[j];
            nx[i]=nx[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      mx[i]=ss;
      nx[i]=tt;
   }
   return(cnm);
}

template <class X,class Y,class Z>
class Heap3 {
   public:
      Heap3(long nm,X *ra,Y *rb,Z *rc); 
     ~Heap3();
      //F stands for forward order and captures the largest values
      void add_pointF(X rx,Y ry,Z rz); //Adds a single point to 
      long orderF(void); //Orders the heap post additions
      //R stands for reverse order and captures the smallest values
      //Returns the number of points captured
      void add_pointR(X rx,Y ry,Z rz); //Adds a single point to 
      long orderR(void); //Orders the heap post additions
      //Returns the number of points captured

      long num; //Number of points allowed
      long cnm; //Current number of points
      X *mx;
      Y *nx;
      Z *ox;
      long ir; //Used in calculations
         //is num-1
};

template <class X,class Y,class Z>
Heap3<X,Y,Z>::Heap3(long nm,X *ra,Y *rb,Z *rc){
   num=nm;
   mx=ra;
   nx=rb;
   ox=rc;
   ir=num-1;
   cnm=0;
   if(num<2){cout << "Error in heap constructor!" << endl;exit(1);}
} 

template <class X,class Y,class Z>
Heap3<X,Y,Z>::~Heap3(){
}  
 
template <class X,class Y,class Z>
void Heap3<X,Y,Z>::add_pointF(X xx,Y yy,Z zz){
   long i,j,k,n;
   X ss;
   Y tt;
   Z uu;

   if(num>cnm){
      mx[cnm]=xx;
      nx[cnm]=yy;
      ox[cnm++]=zz;
      if(cnm==num){
         //Build the initial heap
         k=(num>>1);
         while(k){
            ss=mx[(--k)];
            tt=nx[k];
            uu=ox[k];

            i=k;
            j=((k+1)<<1)-1;
            while(j<=ir){
               if(j<ir && mx[j]>mx[j+1])++j;
               if(ss>mx[j]){
                  mx[i]=mx[j];
                  nx[i]=nx[j];
                  ox[i]=ox[j];
                  j+=(i=j)+1;
               }
               else j=ir+1;
            }
            mx[i]=ss;
            nx[i]=tt;
            ox[i]=uu;
         }
      }
   }
   else {
      ss=*mx;
      if(xx>ss){
         tt=*nx;
         uu=*ox;
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && mx[j]>mx[j+1])++j;
            if(xx>mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               ox[i]=ox[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=xx;
         nx[i]=yy;
         ox[i]=zz;
      }
   }
}

template <class X,class Y,class Z>
long Heap3<X,Y,Z>::orderF(void){
   long i,j,k;
   X ss;
   Y tt;
   Z uu;
   //Build the initial heap if cnm<num
   if(cnm<2)return(cnm);
   if(cnm<num){
      ir=cnm-1;
      k=(cnm>>1);
      while(k){
         ss=mx[(--k)];
         tt=nx[k];
         uu=ox[k];

         i=k;
         j=((k+1)<<1)-1;
         while(j<=ir){
            if(j<ir && mx[j]>mx[j+1])++j;
            if(ss>mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               ox[i]=ox[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=ss;
         nx[i]=tt;
         ox[i]=uu;
      }
   }
   //Order the heap by promotion & filtering
   ir=cnm-1;
   for(;;){
      ss=mx[ir];
      tt=nx[ir];
      uu=ox[ir];
      mx[ir]=mx[0];
      nx[ir]=nx[0];
      ox[ir]=ox[0];
      if((--ir)==0){
         mx[0]=ss;
         nx[0]=tt;
         ox[0]=uu;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && mx[j]>mx[j+1])++j;
         if(ss>mx[j]){
            mx[i]=mx[j];
            nx[i]=nx[j];
            ox[i]=ox[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      mx[i]=ss;
      nx[i]=tt;
      ox[i]=uu;
   }
   return(cnm);
}
   
template <class X,class Y,class Z>
void Heap3<X,Y,Z>::add_pointR(X xx,Y yy,Z zz){
   long i,j,k,n;
   X ss;
   Y tt;
   Z uu;

   if(num>cnm){
      mx[cnm]=xx;
      nx[cnm]=yy;
      ox[cnm++]=zz;
      if(cnm==num){
         //Build the initial heap
         k=(num>>1);
         while(k){
            ss=mx[(--k)];
            tt=nx[k];
            uu=ox[k];

            i=k;
            j=((k+1)<<1)-1;
            while(j<=ir){
               if(j<ir && mx[j]<mx[j+1])++j;
               if(ss<mx[j]){
                  mx[i]=mx[j];
                  nx[i]=nx[j];
                  ox[i]=ox[j];
                  j+=(i=j)+1;
               }
               else j=ir+1;
            }
            mx[i]=ss;
            nx[i]=tt;
            ox[i]=uu;
         }
      }
   }
   else {
      ss=*mx;
      if(xx<ss){
         tt=*nx;
         uu=*ox;
         i=0;
         j=1;
         while(j<=ir){
            if(j<ir && mx[j]<mx[j+1])++j;
            if(xx<mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               ox[i]=ox[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=xx;
         nx[i]=yy;
         ox[i]=zz;
      }
   }
}

template <class X,class Y,class Z>
long Heap3<X,Y,Z>::orderR(void){
   long i,j,k;
   X ss;
   Y tt;
   Z uu;
   //Build the initial heap if cnm<num
   if(cnm<2)return(cnm);
   if(cnm<num){
      ir=cnm-1;
      k=(cnm>>1);
      while(k){
         ss=mx[(--k)];
         tt=nx[k];
         uu=ox[k];

         i=k;
         j=((k+1)<<1)-1;
         while(j<=ir){
            if(j<ir && mx[j]<mx[j+1])++j;
            if(ss<mx[j]){
               mx[i]=mx[j];
               nx[i]=nx[j];
               ox[i]=ox[j];
               j+=(i=j)+1;
            }
            else j=ir+1;
         }
         mx[i]=ss;
         nx[i]=tt;
         ox[i]=uu;
      }
   }
   //Order the heap by promotion & filtering
   ir=cnm-1;
   for(;;){
      ss=mx[ir];
      tt=nx[ir];
      uu=ox[ir];
      mx[ir]=mx[0];
      nx[ir]=nx[0];
      ox[ir]=ox[0];
      if((--ir)==0){
         mx[0]=ss;
         nx[0]=tt;
         ox[0]=uu;
         break;
      }
      i=0;
      j=1;
      while(j<=ir){
         if(j<ir && mx[j]<mx[j+1])++j;
         if(ss<mx[j]){
            mx[i]=mx[j];
            nx[i]=nx[j];
            ox[i]=ox[j];
            j+=(i=j)+1;
         }
         else j=ir+1;
      }
      mx[i]=ss;
      nx[i]=tt;
      ox[i]=uu;
   }
   return(cnm);
}

} 
#endif
