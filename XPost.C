#include <Hash.h>
#include <Doc.h>
#include <DStor.h>
#include <Dmap.h>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <queue>
#include <XPost.h>
namespace iret {

MpTh mt(10);

MpTh::MpTh(int m){
   nthr=m;
   g_fin=new int[m];
   iw=new long[m];
   ic=new long[m];
   sca=new float[m];
}

MpTh::~MpTh(void){
   delete [] g_fin;
   delete [] iw;
   delete [] ic;
   delete [] sca;
}

void betaProc(int k,DSpan *pDn,strMap *pMp){ //k is the number of the thread, beta is original
   int i,j,n,m,u,ia,ib,ix,iu;
   const char *pch;
   bool pushed=true;
   BagPac<int> Bp;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         Bp.unpack(*pDn);
         vector<string>::iterator p=Bp.ss.begin();
         while(p!=Bp.ss.end()){
            pch=(*p).c_str();
            pMp->add_count(pch,1);
            p++;
         }
         Bp.clear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
         //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

void betdProc(int k,string *psxc,strMap *pMp){ //k is the number of the thread, beta is original
   int i,j,n,m,u,ia,ib,ix,iu;
   const char *pch;
   float sx;
   bool pushed=true;
   Doc<float> Dc;
   Dset Ds;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         Dc.extract(*psxc,Ds);
         Ds.Set();
         while(Ds.rs!=Ds.rz){
            pMp->add_count(*Ds.rs,1);
            Ds.rs++;
         }
         Ds.SMclear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
         //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

void mergProc(int k,strMap *pMp,strMap *pNp){ //merges pMp and pNp
   pMp->Merge(pNp);
   {
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      mt.g_fin[k]=-1;
   }
}

void alphProc(int k,Chash *pCh,DSpan *pDn,int *wrd){
   int i,j,*ladr;
   bool pushed=true;
   BagPac<float> Bp;
   vector<string>::iterator p;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         ladr=wrd+mt.ic[k];
         Bp.unpack(*pDn);
         i=0;
         p=Bp.ss.begin();
         while(p!=Bp.ss.end()){
            j=(int)pCh->count((*p).c_str())-1;
            *(ladr+i)=j;
            p++;i++;
         }
         Bp.clear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
      //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

void alphzProc(int k,Chash *pCh,DSpan *pDn,int *wrd,float *xnm){
   int i,j,*ladr;
   float *zadr;
   bool pushed=true;
   BagPac<float> Bp;
   vector<string>::iterator p;
   vector<float>::iterator q;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         ladr=wrd+mt.ic[k];
         zadr=xnm+mt.ic[k];
         Bp.unpack(*pDn);
         i=0;
         p=Bp.ss.begin();
         q=Bp.si.begin();
         while(p!=Bp.ss.end()){
            j=(int)pCh->count((*p).c_str())-1;
            *(ladr+i)=j;
            *(zadr+i)=*q;
            p++;q++;i++;
         }
         Bp.clear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
      //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

void alphzfProc(int k,Chash *pCh,DSpan *pDn,int *wrd,float *xnm,float (*d_local)(int,long)){
   int i,j,*ladr;
   long vi;
   float *zadr,u,sum;
   bool pushed=true;
   BagPac<float> Bp;
   vector<string>::iterator p;
   vector<float>::iterator q;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         ladr=wrd+mt.ic[k];
         zadr=xnm+mt.ic[k];
         Bp.unpack(*pDn);
         i=0;
         sum=0;
         p=Bp.ss.begin();
         q=Bp.si.begin();
         while(p!=Bp.ss.end()){
            j=(int)pCh->count((*p).c_str())-1;
            *(ladr+i)=j;
            if(strstr((*p).c_str(),"!!t"))sum+=*q;
            p++;q++;i++;
         }
         vi=rnd(sum);
         i=0;
         q=Bp.si.begin();
         while(q!=Bp.si.end()){
            j=(int)rnd(*q);
            u=d_local(j,vi);
            *(zadr+i)=u;
            q++;i++;
         }
         Bp.clear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
      //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

void alpdProc(int k,Chash *pCh,string *psxc,int *addr){
   int i,j,*ladr;
   float xx;
   const char *pch;
   bool pushed=true;
   Doc<float> Dc;
   Dset Ds;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         ladr=addr+mt.ic[k];
         i=0;
         Dc.extract(*psxc,Ds);
         Ds.Set();
         while(Ds.rs!=Ds.rz){
            j=(int)pCh->count(*(Ds.rs))-1;
            *(ladr+i)=j;
            Ds.rs++;
            i++;
         }
         Ds.SMclear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
      //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

void alpdzProc(int k,Chash *pCh,string *psxc,int *wrd,float *xnm){
   int i,j,*ladr;
   float xx,*zadr;
   const char *pch;
   bool pushed=true;
   Doc<float> Dc;
   Dmap<float> Dm;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         ladr=wrd+mt.ic[k];
         zadr=xnm+mt.ic[k];
         i=0;
         Dc.extract(*psxc,Dm);
         Dm.Set();
         while(Dm.qs!=Dm.qz){
            j=(int)pCh->count(Dm.qs->first)-1;
            *(ladr+i)=j;
            *(zadr+i)=Dm.qs->second;
            Dm.qs++;
            i++;
         }
         Dm.SMclear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
      //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

void alpdzfProc(int k,Chash *pCh,string *psxc,int *wrd,float *xnm,float (*d_local)(int,long)){
   int i,j,*ladr,vi;
   float xx,*zadr,sum;
   const char *pch;
   bool pushed=true;
   Doc<float> Dc;
   Dmap<float> Dm;

   while((!mt.g_done)||mt.g_cnt){
      if(mt.iw[k]){
         ladr=wrd+mt.ic[k];
         zadr=xnm+mt.ic[k];
         sum=0;
         i=0;
         Dc.extract(*psxc,Dm);
         Dm.Set();
         while(Dm.qs!=Dm.qz){
            j=(int)pCh->count(Dm.qs->first)-1;
            *(ladr+i)=j;
            if(strstr(Dm.qs->first,"!!t"))sum+=Dm.qs->second;
            Dm.qs++;
            i++;
         }
         vi=rnd(sum);
         i=0;
         Dm.Set();
         while(Dm.qs!=Dm.qz){
            *(zadr+i)=d_local(rnd(Dm.qs->second),vi);
            Dm.qs++;
            i++;
         }
         Dm.SMclear();
         pushed=false;
      }
      std::unique_lock<std::mutex> lckque(mt.g_lockque); //obtain a lock of mt.g_lockque
      if(!pushed){
         mt.iw[k]=0;
         mt.g_indices.push(k);
         pushed=true;
      }
      mt.g_quecheck.notify_one(); //Notify main that its mt.g_quecheck wait is over/unblocked
      mt.g_proccheck.wait(lckque); //Wait blocks thread and releases the lock on mt.g_lockque allowing other threads
      //access to mt.g_lockque and the que
   }
   mt.g_fin[k]=1;
}

//STerm

STerm::STerm(void) : FBase("null","null"){
   pLx=NULL;
   pCh=NULL;
   open1=open2=0;
}

STerm::STerm(const char *par,const char *nam) : FBase(par,nam){
   pLx=NULL;
   pCh=NULL;
   open1=open2=0;
}

STerm::STerm(const char *par,int n,const char *nam) : FBase(par,n,nam){
   pLx=NULL;
   pCh=NULL;
   open1=open2=0;
}

STerm::STerm(const char *par,const char *nam,const char *pnam) : FBase(par,nam,pnam){
   pLx=NULL;
   pCh=NULL;
   open1=open2=0;
}

STerm::~STerm(void){
   if(pLx)delete pLx;
   if(pCh)delete pCh;
}

void STerm::SetMem(STerm &St){
   if(St.open2){
      this->pCh=new Chash;
      this->pCh->SetMem(*St.pCh);
      this->open2=1;
   }
   if(St.open1){
      this->pLx=new Lexos;
      this->pLx->SetMem(*St.pLx);
      this->open1=1;
   }
}

void STerm::create(strMap &Mp){
   pLx=new Lexos(name);
   map_down_sub((FBase *)pLx,"lexos");
   pLx->create_Lexos(Mp);
   pCh=new Chash(name);
   map_down_sub((FBase *)pCh,"chash");
   pCh->create_ctable_STerm(Mp,3);
}

Lexos *STerm::gopen_Lexos(void){
   if(!open1){
      pLx=new Lexos(name);
      map_down_sub((FBase *)pLx,"lexos");
      pLx->gopen_map();
      open1=1;
   }
   return(pLx);
}

void STerm::gclose_Lexos(void){
   if(open1){
      pLx->gclose_map();
      delete pLx;
      pLx=NULL;
      open1=0;
   }
}

Chash *STerm::gopen_Chash(void){
   if(!open2){
      delete pCh;
      pCh=new Chash(name);
      map_down_sub((FBase *)pCh,"chash");
      pCh->gopen_ctable_map();
      open2=1;
   }
   return(pCh);
}

void STerm::gclose_Chash(void){
   if(open2){
      pCh->gclose_ctable_map();
      delete pCh;
      pCh=NULL;
      open2=0;
   }
}

} 
