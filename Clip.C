#include <Clip.h>
namespace iret {

Clip::Clip(void) : Count(){
   num=0;
   pck=NULL;
   nck=NULL;
   word_space=10000;
   lst=new char*[word_space];
} 

Clip::Clip(long nm) : Count(){
   num=nm;
   pck=new char[num];
   nck=new long[num];
   word_space=10000;
   lst=new char*[word_space];
} 

Clip::~Clip(){
    if(pck!=NULL)delete [] pck;
    if(nck!=NULL)delete [] nck;
    delete [] lst;
}  
 
void Clip::set1(void){
   long i,j,k;

   for(i=0;i<32;i++)zt[i]=0;
   for(i=32;i<48;i++)zt[i]=1;
   for(i=48;i<58;i++)zt[i]=2;
   for(i=58;i<65;i++)zt[i]=1;
   for(i=65;i<91;i++)zt[i]=3;
   for(i=91;i<97;i++)zt[i]=1;
   for(i=97;i<123;i++)zt[i]=3;
   for(i=123;i<127;i++)zt[i]=1;
   zt[127]=0;
}

void Clip::prc0(const char *pch){
   long i,j,k,len,jo,j2;
   char c;

   strcpy(pck,pch);
   len=strlen(pck);
   for(j=0;j<len;j++)if(zt[pck[j]]==0)pck[j]=32;
   j=0;
   while(j<len){
      while((j<len)&&(zt[pck[j]]<2))j++; //find beginning of string
      if(j<len){
         k=0;
         while(zt[pck[j+k]]>=2)k++; //Find end of string
         c=pck[j+k];
         pck[j+k]='\0';
         add_count2(pck+j,1);
         pck[j+k]=c;
         j+=k;
      }
   }
}

void Clip::prc1(const char *pch){
   long i,j,k,len,jo,j2;

   strcpy(pck,pch);
   i=0;
   while(pck[i]){
      if(!zt[pck[i]])pck[i]=32;
      i++;
   }

   i=j=0;
   while(pck[j]){ 
      while(zt[pck[i]]==1)i++;
      j=i;
      while(zt[pck[j]]>1)j++; 
      if(pck[j])j2=j+1;//Mark end of processing
      else j2=j;//so process ends at end of string
      while(j>i){
         pck[j]='\0';
         k=i;
         while(k<j){
            add_count2(pck+k,1);
            jo=zt[pck[k]];
            while((k<j)&&(zt[pck[k]]==jo))k++;
         }
         jo=zt[pck[j-1]];
         k=j-1;
         while((k>=i)&&(zt[pck[k]]==jo))k--;
         j=k+1;
      }
      i=j=j2;
   }
} 

void Clip::prc2(const char *pch,Hash &Hs){
   long i,j,k,len,jo,j2;

   len=strlen(pch);
   strcpy(pck,pch);
   for(j=0;j<len;j++)if(zt[pck[j]]==0)pck[j]=32;
   i=0;
   k=1;
   jo=zt[pck[0]];
   for(j=0;j<len;j++){
      if(zt[pck[j]]!=jo){k++;jo=zt[pck[j]];}
      if(zt[pck[j]]!=1){
         pck[i]=pck[j];
         nck[i]=k;
         i++;
      }
   }
   len=i;
   jo=nck[0];
   j=0;
   while(j<len){
      k=1;
      while((j+k<len)&&(nck[j+k]==jo)){nck[j+k]=0;k++;}
      nck[j]=1;
      j+=k;
      jo=nck[j];
   }
   j=len;
   while(j>0){
      pck[j]='\0';
      k=j-1;
      while(k>=0){
         while(!nck[k])k--;
         if(Hs.find(pck+k))add_count2(pck+k,1);
         k--;
      }
      j--;
      while(!nck[j])j--;
   }
}

void Clip::prc2(const char *pch,Chash &Cs,long llim){
   long i,j,k,len,jo,j2;

   len=strlen(pch);
   strcpy(pck,pch);
   for(j=0;j<len;j++)if(zt[pck[j]]==0)pck[j]=32;
   i=0;
   k=1;
   jo=zt[pck[0]];
   for(j=0;j<len;j++){
      if(zt[pck[j]]!=jo){k++;jo=zt[pck[j]];}
      if(zt[pck[j]]!=1){
         pck[i]=pck[j];
         nck[i]=k;
         i++;
      }
   }
   len=i;
   jo=nck[0];
   j=0;
   while(j<len){
      k=1;
      while((j+k<len)&&(nck[j+k]==jo)){nck[j+k]=0;k++;}
      nck[j]=1;
      j+=k;
      jo=nck[j];
   }
   j=len;
   while(j>0){
      pck[j]='\0';
      k=j-1;
      while(k>=0){
         while(!nck[k])k--;
         if((j2=Cs.count(pck+k))&&(j2<llim))add_count2(pck+k,1);
         k--;
      }
      j--;
      while(!nck[j])j--;
   }
}

void Clip::prc3(const char *pch,Hash &Hs){
   long i,j,k,len,jo,j2;

   len=strlen(pch);
   strcpy(pck,pch);
   for(j=0;j<len;j++)if(zt[pck[j]]==0)pck[j]=32;
   i=0;
   for(j=0;j<len;j++){
      if(zt[pck[j]]>1){
         pck[i]=pck[j];
         i++;
      }
   }
   len=i;
   j=len;
   while(j>0){
      pck[j]='\0';
      k=j-1;
      while(k>=0){
         if(Hs.find(pck+k))add_count2(pck+k,1);
         k--;
      }
      j--;
   }
}

void Clip::prc3(const char *pch,Chash &Cs,long llim){
   long i,j,k,len,jo,j2;

   len=strlen(pch);
   strcpy(pck,pch);
   for(j=0;j<len;j++)if(zt[pck[j]]==0)pck[j]=32;
   i=0;
   for(j=0;j<len;j++){
      if(zt[pck[j]]>1){
         pck[i]=pck[j];
         i++;
      }
   }
   len=i;
   j=len;
   while(j>0){
      pck[j]='\0';
      k=j-1;
      while(k>=0){
         if((jo=Cs.count(pck+k))&&(jo<llim))add_count2(pck+k,1);
         k--;
      }
      j--;
   }
}

//Simple tokenization
void Clip::spac(const char *str){
   long i,j;
   char *pch,c;
   pch=new char[strlen(str)+1];
   strcpy(pch,str);

   i=0;
   while(pch[i]){
      while(pch[i]==' ')i++;
      j=i;
      while((c=pch[i])&&(c!=' '))i++;
      if(i>j){
         pch[i]='\0';
         add_count2(pch+j,1);
      }
      if(c)i++;
   }
   delete [] pch;
}

void Clip::spac(const char *str,long n){
   long i,j;
   char *pch,c;
   pch=new char[strlen(str)+1];
   strcpy(pch,str);

   i=0;
   while(pch[i]){
      while(pch[i]==' ')i++;
      j=i;
      while((c=pch[i])&&(c!=' '))i++;
      if(i>j){
         pch[i]='\0';
         add_count2(pch+j,n);
      }
      if(c)i++;
   }
   delete [] pch;
}

//fill_grm functions
void Clip::fill_grm(const char *str,int gram_size,int aug,Word &Wrd){
   long j,k,n,lxn,ln=strlen(str);
   int *pint,gp,grmm2;
   char ch,srr[256],*pch,brd[256],frs[3],scc[4];

   grmm2=gram_size-2;
   brd[0]='\0';
   brd[1]=' ';
   scc[2]=frs[1]='$';
   scc[3]=frs[2]='\0';
   n=0;
   j=0;
   Wrd.convert_protected(str,ln);
   while(lxn=Wrd.wordf(j,ln)){
         k=0;
         gp=(lxn<grmm2)?lxn:grmm2;
         while(k<gp){brd[2+k]=Wrd.wrd[k];k++;}
         brd[k+2]='\0';
         if(brd[0]!='\0')this->add_count2(brd,1);

         scc[0]=frs[0]=brd[0]=Wrd.wrd[0];
         scc[1]=Wrd.wrd[1];
         if(lxn>gram_size){
            //Handle first character
            this->add_count2(frs,1);
            //Add the 2gram
            this->add_count2(scc,1);
            //Define first ngram
            ch=Wrd.wrd[gram_size];
            Wrd.wrd[gram_size]='\0';
            //Enter the first ngram into the tree
            strcpy(srr,Wrd.wrd);
            this->add_count2(srr,1);
            //Produce the marked ngram and enter into tree
            strcat(srr,"!");
            if(aug>0)this->add_count2(srr,aug);
            Wrd.wrd[gram_size]=ch;
            //Enter the remainder of the ngrams into tree
            for(k=1;k<lxn-gram_size+1;k++){
               ch=Wrd.wrd[gram_size+k];
               Wrd.wrd[gram_size+k]='\0';
               strcpy(srr,&Wrd.wrd[k]);
               this->add_count2(srr,1);
               Wrd.wrd[gram_size+k]=ch;
            }
         }
         else {
            //Handle first character
            this->add_count2(frs,1);
            if(lxn>1)this->add_count2(scc,1);
            //Enter the first ngram into tree
            strcpy(srr,Wrd.wrd);
            this->add_count2(srr,1);
            //Mark the ngram as first and enter into tree
            strcat(srr,"!");
            if(aug>0)this->add_count2(srr,aug);
         }
   }
}

void Clip::fill_grmo(const char *str,int gram_size,int aug,Word &Wrd){
   long j,k,n,lxn,ln=strlen(str);
   int *pint,gp,grmm2;
   char ch,srr[256],*pch,brd[256],frs[3],scc[4];

   grmm2=gram_size-2;
   brd[0]='\0';
   brd[1]=' ';
   frs[1]='$';
   frs[2]='\0';
   n=0;
   j=0;
   Wrd.convert_protected(str,ln);
   while(lxn=Wrd.wordf(j,ln)){
         k=0;
         gp=(lxn<grmm2)?lxn:grmm2;
         while(k<gp){brd[2+k]=Wrd.wrd[k];k++;}
         brd[k+2]='\0';
         if(brd[0]!='\0')this->add_count2(brd,1);

         frs[0]=brd[0]=Wrd.wrd[0];
         if(lxn>gram_size){
            //Handle first character
            this->add_count2(frs,1);
            //Define first ngram
            ch=Wrd.wrd[gram_size];
            Wrd.wrd[gram_size]='\0';
            //Enter the first ngram into the tree
            strcpy(srr,Wrd.wrd);
            this->add_count2(srr,1);
            //Produce the marked ngram and enter into tree
            strcat(srr,"!");
            if(aug>0)this->add_count2(srr,aug);
            Wrd.wrd[gram_size]=ch;
            //Enter the remainder of the ngrams into tree
            for(k=1;k<lxn-gram_size+1;k++){
               ch=Wrd.wrd[gram_size+k];
               Wrd.wrd[gram_size+k]='\0';
               strcpy(srr,&Wrd.wrd[k]);
               this->add_count2(srr,1);
               Wrd.wrd[gram_size+k]=ch;
            }
         }
         else {
            //Handle first character
            this->add_count2(frs,1);
            //Enter the first ngram into tree
            strcpy(srr,Wrd.wrd);
            this->add_count2(srr,1);
            //Mark the ngram as first and enter into tree
            strcat(srr,"!");
            if(aug>0)this->add_count2(srr,aug);
         }
   }
}

void Clip::fill_brm(const char *str,int lim_gram_size,Word &Wrd){
   long i,j,k,n,lxn,ln=strlen(str);
   int *pint,gp,grmm2;
   char *chd,*chz,*czz,*pch,cho;
   char *zlt;

   Wrd.convert_protected(str,ln);
   zlt=Wrd.acc_conv();
   pch=new char[ln+1];
   j=0;
   for(i=0;i<ln;i++){
      if(gp=zlt[i])pch[j++]=gp;
   }
   chd=pch;
   if(j>1)chz=pch+2;
   for(i=1;i<j;i++){
      czz=chd;
      cho=*chz;
      *chz='\0';
      while((czz+1)!=chz){
         this->add_count2(czz,1);
         czz++;
      }
      *chz=cho;
      chz++;
      if(i+2>lim_gram_size)chd++;
   }
   delete [] pch;
}

void Clip::fill_scn(const char *str,int lim_gram_size,Word &Wrd){
   long i,j,lxn,ln=strlen(str);
   char *chd;

   j=0;
   Wrd.convert_protected(str,ln);
   while(lxn=Wrd.wordf(j,ln)){
      if(lxn>=3){
         chd=new char[lxn+1];
         strcpy(chd,Wrd.wrd);
         i=3;
         while(i<=lxn){
            chd[i]='\0';
            add_count2(chd+i-3,1);
            chd[i]=Wrd.wrd[i];
            i++;
         }
      }
   }
}

void Clip::Set_Up_Word(Word &Wrd){
   Wrd.byte_lim=64;
   Wrd.back_lim=48;
   Wrd.set_map(".,:;!?",'\024',LOWERCASE);
   Wrd.restore_map("\'");
}

void Clip::proc_MrkSng(long lg,const char *txt,Word &Wrd,const char *suf){
   long k,len;
   char chp[max_str];

   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      strcat(chp,suf);
      add_count2(chp,1);
   }
   Wrd.clear_list();
}

void Clip::proc_title(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];

   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='T';
      *(chp+len+3)='\0';
      add_countz(chp,1);
      *(chp+len+2)='t';
      add_count2(chp,2);
   }
   Wrd.clear_list();
   Wrd.multiple(2);
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='P';
      *(chp+len+3)='\0';
      add_countz(chp,1);
      *(chp+len+2)='p';
      add_countz(chp,1);
    }
    Wrd.clear_list();
}


void Clip::proc_body(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];

   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='t';
      *(chp+len+3)='\0';
      add_count2(chp,1);
   }
   Wrd.clear_list();
   Wrd.multiple(2);
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='p';
      *(chp+len+3)='\0';
      add_countz(chp,1);
    }
    Wrd.clear_list();
}

void Clip::proc_title_stem(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];

   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single_stem();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='T';
      *(chp+len+3)='\0';
      add_count2(chp,1);
      *(chp+len+2)='t';
      add_count2(chp,2);
   }
   Wrd.clear_list();
   Wrd.multiple_stem(2);
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='P';
      *(chp+len+3)='\0';
      add_count2(chp,1);
      *(chp+len+2)='p';
      add_count2(chp,1);
    }
    Wrd.clear_list();
}


void Clip::proc_body_stem(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];

   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single_stem();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='t';
      *(chp+len+3)='\0';
      add_count2(chp,1);
   }
   Wrd.clear_list();
   Wrd.multiple_stem(2);
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='p';
      *(chp+len+3)='\0';
      add_count2(chp,1);
    }
    Wrd.clear_list();
}

void Clip::proc_tisng(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];
   
   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='T';
      *(chp+len+3)='\0';
      add_countz(chp,1);
      *(chp+len+2)='t';
      add_count2(chp,2);
   } 
   Wrd.clear_list();
}     
      
      
void Clip::proc_bdsng(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];
            
   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='t';
      *(chp+len+3)='\0';
      add_count2(chp,1);
   }
   Wrd.clear_list();
}

void Clip::proc_tisng_stem(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];

   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single_stem();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='T';
      *(chp+len+3)='\0';
      add_count2(chp,1);
      *(chp+len+2)='t';
      add_count2(chp,2);
   }
   Wrd.clear_list();
}


void Clip::proc_bdsng_stem(long lg,const char *txt,Word &Wrd){
   long k,len;
   char chp[max_str];

   Wrd.convert_protected(txt,lg);
   Wrd.modify('\'');
   Wrd.single_stem();
   for(k=0;k<Wrd.cnt;k++){
      strcpy(chp,*(Wrd.list+k));
      len=strlen(chp);
      *(chp+len)='!';
      *(chp+len+1)='!';
      *(chp+len+2)='t';
      *(chp+len+3)='\0';
      add_count2(chp,1);
   }
   Wrd.clear_list();
}

void Clip::proc_Vbarstring(const char *txt){
   char *ptr,*pch;
   strcpy(xmam,txt);
   ptr=xmam;
   while(pch=strchr(ptr,'|')){
      *pch='\0';
      if(*ptr){
         if(strchr(ptr,'*'))add_countz(ptr,1);
         else add_count2(ptr,1);
      }
      ptr=pch+1;
  }
  if(*ptr){
     if(strchr(ptr,'*'))add_countz(ptr,1);
     else add_count2(ptr,1);
  }
}

void Clip::substr(const char *tok,long ng){
   long i,j,k,n=strlen(tok);
   char *ptr=new char[n+1];
   for(i=0;i<n;i++){
      k=(ng<n-i)?ng:n-i;
      for(j=0;j<k;j++){
         ptr[j]=tok[i+j];
         ptr[j+1]='\0';
         add_count2(ptr,1);
      }
   }
   delete [] ptr;
}

void Clip::substr_b(const char *tok,long ng){
   long i,j,k,n=strlen(tok);
   char *ptr=new char[n+3];
   ptr[0]='_';
   k=(ng<n)?ng:n;
   for(j=0;j<k;j++){
      ptr[j+1]=tok[j];
      ptr[j+2]='\0';
      add_count2(ptr,1);
   }
   if(n<=ng){
      ptr[j+1]='_';
      ptr[j+2]='\0';
      add_count2(ptr,1);
   }
   strcpy(ptr,tok);
   ptr[n]='_';
   ptr[n+1]='\0';
   for(j=n-k;j<n;j++){
      add_count2(ptr+j,1);
   }
   delete [] ptr;
}

void Clip::lsubstr(const char *tok,long ng){
   long i,j,k,n=strlen(tok);
   char *ptr=new char[n+1];
   for(i=0;i<n;i++){
      k=(ng<n-i)?ng:n-i;
      for(j=0;j<k;j++){
         ptr[j]=tolower(tok[i+j]);
         ptr[j+1]='\0';
         add_count2(ptr,1);
      }
   }
   delete [] ptr;
}

void Clip::lsubstr_b(const char *tok,long ng){
   long i,j,k,n=strlen(tok);
   char *ptr=new char[n+3];
   ptr[0]='_';
   k=(ng<n)?ng:n;
   for(j=0;j<k;j++){
      ptr[j+1]=tolower(tok[j]);
      ptr[j+2]='\0';
      add_count2(ptr,1);
   }
   if(n<=ng){
      ptr[j+1]='_';
      ptr[j+2]='\0';
      add_count2(ptr,1);
   }
   for(i=0;i<n;i++)ptr[i]=tolower(tok[i]);
   ptr[n]='_';
   ptr[n+1]='\0';
   for(j=n-k;j<n;j++){
      add_count2(ptr+j,1);
   }
   delete [] ptr;
}

void Clip::token(const char *tok){
   add_count2(tok,1);
}

void Clip::ltoken(const char *tok){
   long i,n=strlen(tok);
   char *ptr=new char[n+1];
   for(i=0;i<n;i++)ptr[i]=tolower(tok[i]);
   ptr[n]='\0';
   add_count2(ptr,1);
   delete [] ptr;
}

void Clip::bigram(const char *tk1,const char *tk2){
   long i,n=strlen(tk1),m=strlen(tk2);
   char *ptr=new char[n+m+2];

   strcpy(ptr,tk1);
   ptr[n]=' ';
   ptr[n+1]='\0';
   strcat(ptr,tk2);
   add_count2(ptr,1);
   delete [] ptr;
}

void Clip::lbigram(const char *tk1,const char *tk2){
   long i,n=strlen(tk1),m=strlen(tk2);
   char *ptr=new char[n+m+2];

   for(i=0;i<n;i++)ptr[i]=tolower(tk1[i]);
   ptr[n]=' ';
   for(i=n+1;i<n+1+m;i++)ptr[i]=tolower(tk2[i-n-1]);
   ptr[n+1+m]='\0';
   add_count2(ptr,1);
   delete [] ptr;
}

void Clip::lincode(double xd,double acc){
   long i,j,k;
   char cnam[max_str];

   if(xd==0)return;
   if(xd<0){cnam[0]='n';xd=-xd;}
   else cnam[0]='p';
   i=1;
   while(i*acc<xd){
      long_str(cnam+1,i);
      add_count2(cnam,1);
      i++;
   }
}

void Clip::lincode(double xd,double acc,const char *prf){
   long i,j,k;
   char cnam[max_str];
   strcpy(cnam,prf);
   long len=strlen(cnam)+1;

   if(xd==0)return;
   if(xd<0){cnam[len-1]='n';xd=-xd;}
   else cnam[len-1]='p';
   i=1;
   while(i*acc<xd){
      long_str(cnam+len,i);
      add_count2(cnam,1);
      i++;
   }
}

void Clip::bincode(double xd,long acc){
   long i,j,k;
   double xx,yy,zz;
   char cnam[max_str];

   if(xd==0)return;
   if(xd<0){cnam[0]='n';xd=-xd;}
   else cnam[0]='p';
   i=0;
   xx=1.0;
   while((i<100)&&(xx<xd)){xx*=2.0;i++;}
   if(i==100)xd=xx*0.99999;
   j=i;
   while(j>=-acc){
      if((xx>xd)&&(xx<=2.0*xd)){
         long_str(cnam+1,j);
         add_count2(cnam,1);
         xd-=0.5*xx;
      }
      xx*=0.5;
      j--;
   }   
}

void Clip::bincode(double xd,long acc,const char *prf){
   long i,j,k,len;
   double xx,yy,zz;
   char cnam[max_str];
   strcpy(cnam,prf);
   len=strlen(cnam)+1;

   if(xd==0)return;
   if(xd<0){cnam[len-1]='n';xd=-xd;}
   else cnam[len-1]='p';
   i=0;
   xx=1.0;
   while((i<100)&&(xx<xd)){xx*=2.0;i++;}
   if(i==100)xd=xx*0.99999;
   j=i;
   while(j>=-acc){
      if((xx>xd)&&(xx<=2.0*xd)){
         long_str(cnam+len,j);
         add_count2(cnam,1);
         xd-=0.5*xx;
      }
      xx*=0.5;
      j--;
   }
}

void Clip::logcode(long n,long amp){
   long i,j,k;
   char cnam[max_str];
   
   if(n==0)return;
   if(n<0){cnam[0]='n';n=-n;}
   else cnam[0]='p';

   i=j=0;
   k=1;
   while(i<n){
      i+=k;
      j++;
      long_str(cnam+1,j);
      add_count2(cnam,1);
      if(j%amp==0)k*=2;
   }
}

void Clip::logcode(long n,long amp,const char *prf){
   long i,j,k,len;
   char cnam[max_str];
   strcpy(cnam,prf);
   len=strlen(cnam)+1;

   if(n==0)return;
   if(n<0){cnam[len-1]='n';n=-n;}
   else cnam[len-1]='p';

   i=j=0;
   k=1;
   while(i<n){
      i+=k;
      j++;
      long_str(cnam+len,j);
      add_count2(cnam,1);
      if(j%amp==0)k*=2;
   }
}

void Clip::mstring(const char *prf,long m){
   char cnam[30];
   strcpy(cnam,prf);
   long len=strlen(cnam);
   long_str(cnam+len,m);
   add_count2(cnam,1);
}

void Clip::mstring(const char *prf,long m1,long m2){
   char cnam[50],bnam[30];
   strcpy(cnam,prf);
   long len=strlen(cnam);
   long_str(cnam+len,m1);
   long_str(bnam,m2);
   strcat(cnam,"/");
   strcat(cnam,bnam);
   add_count2(cnam,1);
}

void Clip::mstring(const char *prf,long m1,long m2,long m3){
   char cnam[70],bnam[20],anam[20];
   strcpy(cnam,prf);
   long len=strlen(cnam);
   long_str(cnam+len,m1);
   long_str(bnam,m2);
   long_str(anam,m3);
   strcat(cnam,"/");
   strcat(cnam,bnam);
   strcat(cnam,"/");
   strcat(cnam,anam);
   add_count2(cnam,1);
}

void Clip::prefix(const char *prf,Count &Ct){
   char cnam[max_str];
   long len=strlen(prf);
   strcpy(cnam,prf);
   Ct.node_first();
   while(Ct.node_next()){
      strcpy(cnam+len,Ct.show_str());
      add_count2(cnam,1);
   }
}

void Clip::prefix(const char *prf,Count &Ct,List &Lt){
   char cnam[max_str];
   long len=strlen(prf);
   strcpy(cnam,prf);
   Ct.node_first();
   while(Ct.node_next()){
      strcpy(cnam+len,Ct.show_str());
      if(Lt.search(cnam)){
         add_count2(cnam,1);
      }
   }
}

void Clip::suffix(const char *suf,Count &Ct){
   char cnam[max_str];
   Ct.node_first();
   while(Ct.node_next()){
      strcpy(cnam,Ct.show_str());
      strcat(cnam,suf);
      add_count2(cnam,1);
   }
}

void Clip::suffix(const char *suf,Count &Ct,List &Lt){
   char cnam[max_str];
   Ct.node_first();
   while(Ct.node_next()){
      strcpy(cnam,Ct.show_str());
      strcat(cnam,suf);
      if(Lt.search(cnam)){
         add_count2(cnam,1);
      }
   }
}

void Clip::addin(Count &Ct){
   Ct.node_first();
   while(Ct.node_next()){
      add_count2(Ct.show_str(),Ct.count());
   }
}

void Clip::addin(Count &Ct,List &Lt){
   Ct.node_first();
   while(Ct.node_next()){
      if(Lt.search(Ct.show_str())){
         add_count2(Ct.show_str(),Ct.count());
      }
   }
}

void Clip::ptoken(const char *prf,const char *pch){
   char cnam[max_str];
   strcpy(cnam,prf);
   strcat(cnam,pch);
   add_count2(cnam,1);
}

void Clip::stoken(const char *suf,const char *pch){
   char cnam[max_str];
   strcpy(cnam,pch);
   strcat(cnam,suf);
   add_count2(cnam,1);
}

void Clip::prod(Count &C1,Count &C2){
   long i,j,n=C1.cnt_key;
   char **cns=(char**)new long[n];
   char cnam[max_str];

   C1.node_first();
   i=0;
   while(C1.node_next()){
      cns[i++]=C1.show_str();
   }
   C2.node_first();
   while(C2.node_next()){
      strcpy(cnam,C2.show_str());
      j=strlen(cnam);
      strcat(cnam,"*");
      j++;
      for(i=0;i<n;i++){
         strcpy(cnam+j,cns[i]);
         add_count2(cnam,1);
      }
   }
   delete [] cns;
}
      
void Clip::prod(Count &C1,Count &C2,List &Lt){
   long i,j,n=C1.cnt_key;
   char **cns=(char**)new long[n];
   char cnam[max_str];

   C1.node_first();
   i=0;
   while(C1.node_next()){
      cns[i++]=C1.show_str();
   }
   C2.node_first();
   while(C2.node_next()){
      strcpy(cnam,C2.show_str());
      j=strlen(cnam);
      strcat(cnam,"*");
      j++;
      for(i=0;i<n;i++){
         strcpy(cnam+j,cns[i]);
         if(Lt.search(cnam))add_count2(cnam,1);
      }
   }
   delete [] cns;
}

int Clip::split(const char *str){
   int i,j=0,k=0;
   while(str[j]){
      while(isblank(str[j]))j++;
      i=j;
      while((str[j])&&(!isblank(str[j])))j++;
      if(j>i){
         lst[k]=new char[j-i+1];
         strncpy(lst[k],str+i,j-i);
         lst[k][j-i]='\0';
         k++;
      }
   }
   nx=k;
   return(nx);
}

int Clip::split_lower(const char *str){
   char *pt;
   strcpy(xmam,str);
   pt=xmam;
   while(*pt){
      *pt=tolower(*pt);
      pt++;
   }
   split(xmam);
   return(nx);
}

}
