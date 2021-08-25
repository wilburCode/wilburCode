#include "AbbrvE.h"
namespace iret {

  AbbrvE::AbbrvE(long ta,long wrd_spc){
    tta=ta;
    word_space=wrd_spc;
    abbl=new char*[tta];
    abbs=new char*[tta];
    rt=new int[tta];
    nt=new int[tta];
    lst=new char*[word_space];
    numa=num=0;
    pMt=new MPtok;
    setup_Test();
  }

  AbbrvE::~AbbrvE(){
    if(numa)cleara();
    clear();
    delete [] abbl;
    delete [] abbs;
    delete [] rt;
    delete [] nt;
    delete [] lst;
    delete pMt;
  }
  
void AbbrvE::Extract(char*pch){
   long i,j,k,u,flag;
   int ix;

   if ( strlen(pch) <= 0 )	// no text to look at
     return;

   token(pch);

   i=j=k=0;
   flag=0;
   while(i<num){
      if(!strcmp("(",lst[i])){
         if(flag)k=j+1;
         if((i>k)&&(strcmp(")",lst[i-1]))){
            j=i;
            flag=1;
         }
      }
      if(!strcmp(")",lst[i])){
         if(!flag){j=k=i+1;}
         else {
            if(((j>k)&&(i<j+12))&&(i>j+1)){
               if(k<j-10)k=j-10;
               strcpy(cnam,lst[k]);
               for(u=k+1;u<j;u++){
                  strcat(cnam," ");
                  strcat(cnam,lst[u]);
               }
               ix=strlen(cnam);
               abbl[numa]=new char[ix+1];
               strcpy(abbl[numa],cnam);

               strcpy(cnam,lst[j+1]);
               for(u=j+2;u<i;u++){
                  strcat(cnam," ");
                  strcat(cnam,lst[u]);
               }
               nt[numa]=i-j-1;
               ix=strlen(cnam);
               abbs[numa]=new char[ix+1];
               strcpy(abbs[numa],cnam);
               if(Test(abbs[numa]))numa++;
	       else{ //if test done earlier would not need to allocate memory
		     //until known to be needed
		 delete [] abbs[numa];
		 delete [] abbl[numa];
	       }
               flag=0;
            }
            else {
               flag=0;
               k=i+1;
            }
         }
      }
      i++;
   }
}


//modified  Jan-9-2008
//extract SF in [], parse until ';' or ',' in () or []
void AbbrvE::Extract2(const char*pch){
   long i,j,k,u,ii,jj,kk,flag;
   int ix;
   char openCh[2], closeCh[2];

   token2(pch); // alpha beta (AB) -> alpha beta ( AB )
 
   for(jj=0; jj<2; jj++) {//deal with both () & []
   i=j=k=0;
   flag=0;
   
   if(jj==0) { strcpy(openCh,"("); strcpy(closeCh,")"); }
   else if(jj==1) { strcpy(openCh,"["); strcpy(closeCh,"]"); }

   while(i<num){
      if(!strcmp(openCh,lst[i])){
         if(flag)k=j+1; //increment after seeing both '(' and ')'
         if((i>k)&&(strcmp(closeCh,lst[i-1]))){
            j=i; //index of '(' 
            flag=1;
         }
      }
      if(!strcmp(closeCh,lst[i])){
         if(!flag){j=k=i+1;} //next token
         else {
            if(((j>k)&&(i<j+12))&&(i>j+1)){ 
               if(k<j-10)k=j-10; 
               strcpy(cnam,lst[k]);
               for(u=k+1;u<j;u++){
                  strcat(cnam," ");
                  strcat(cnam,lst[u]);
               }
               ix=strlen(cnam);
               abbl[numa]=new char[ix+1];
               strcpy(abbl[numa],cnam);

               strcpy(cnam,lst[j+1]);
               for(u=j+2;u<i;u++){
                  strcat(cnam," ");
                  strcat(cnam,lst[u]);
               }
               nt[numa]=i-j-1; //# abbr tokens
               ix=strlen(cnam);
	       
	       //---- parse until ';' or ',' 
	       ii=0;
	       while(ii<ix) {
		 if( ((cnam[ii]==';')&&(cnam[ii+1]==' ')) || 
                     ((cnam[ii]==',')&&(cnam[ii+1]==' ')) ) {
                    ix=ii+1;
		    cnam[ii]='\0';
		    break;
		  }
		  ii++;
	       }
	       //----
	       
               abbs[numa]=new char[ix+1];
               strcpy(abbs[numa],cnam);
                 if(Test(abbs[numa]))numa++;
		 else{ //if test done earlier would not need to allocate memory
		   //until known to be needed
		   delete [] abbs[numa];
		   delete [] abbl[numa];
		 }
               flag=0;
            }
            else {
               flag=0;
               k=i+1;
            }
         }
      }
      i++;
   }
   }
}


void AbbrvE::token(const char *pch){
   long i=1,j=0,k=0;
   long u=1,flag=0;
   char c,*str=cnam;
   int size=strlen(pch);
   if(size>cnam_size) {
     cerr<<"Scratch space "<<cnam_size<<", needed "<<size<<endl;
     exit(1);
   }
   clear();                     // ready space for tokens
   cnam[0]=pch[0];
   while(c=pch[i]){
      switch(c){
         case '(': if(isblank(str[u-1])){
                      str[u++]=pch[i++];
                      if(!isblank(pch[i])){
                         str[u++]=' ';
                      }
                      flag=1;
                   }
                   else str[u++]=pch[i++];
                   break;
         case ')': if(flag){
                      if(!isblank(str[u-1])){
                         str[u++]=' ';
                         str[u++]=pch[i++];
                      }
                      if(!isblank(pch[i]))str[u++]=' ';
                      flag=0;
                   }
                   else str[u++]=pch[i++];
                   break;
         default: str[u++]=pch[i++];
      }
   }
   while((u>0)&&(isblank(str[u-1])))u--;
   str[u]='\0';    
  
   while(str[j]){
      while(isblank(str[j]))j++;
      i=j;
      while((str[j])&&(!isblank(str[j])))j++;
      lst[k]=new char[j-i+1];
      strncpy(lst[k],str+i,j-i);
      lst[k][j-i]='\0';
      if(str[j]){
         k++;
         j++;
      }
   }
   num=k+1;
}


//both () & [] Jan-9-2008
//(G(1)) -> ( G(1) ) Jan-28-2008
void AbbrvE::token2(const char *pch){
   long i=1,j=0,k=0;
   long u=1;
   vector<bool> openChFlag1,openChFlag2;
   long cflag;
   long ii, jj, kk, sz;
   char c,*str=cnam;
   clear();                     // ready space for tokens
   cnam[0]=pch[0];
   while(c=pch[i]){
      switch(c){
        case '(':   	   
                   //--- (h)alpha -> (h)alpha, (h)-alpha -> ( h ) -alpha
                   ii=kk=i;
		   cflag=0; 
		   while(pch[ii] && !isblank(pch[ii])) { //pch[ii] can be '\0'
		     if(pch[ii]=='(') cflag -= 1;
		     else if(pch[ii]==')') { cflag += 1; kk=ii; }
		     ii++; 
		   }

		   if(!cflag && isalnum(pch[kk+1])) { //if alnum right after ')'
		     while(i<ii) str[u++]=pch[i++];
		     break;
		   }
		   //---
		   
	           if(isblank(str[u-1])){
                      str[u++]=pch[i++];
                      if(!isblank(pch[i])){ 
                         str[u++]=' ';
                      }
                      openChFlag1.push_back(true);
                   }
                   else {
		      str[u++]=pch[i++];
                      openChFlag1.push_back(false);
		   }

                   break;

        case ')':  sz = openChFlag1.size();
	           if(sz>0 && openChFlag1[sz-1]){ //modified Jan-28-08
		      if(!isblank(str[u-1])){ 
			 str[u++]=' ';
                         str[u++]=pch[i++]; //pch[i++] is ')'
                      }
		      //---added (Jan-11-08): (BIV; ), -> ( BIV; ) ,
		      else if(!isblank(pch[i+1])){ 
			 str[u++]=pch[i++]; //pch[i++] is ')'
		      }
	              //---

                      if(!isblank(pch[i]))str[u++]=' '; //pch[i] must be after ')'
                   }
                   else str[u++]=pch[i++];

		   if(sz>0) openChFlag1.pop_back();

                   break;

         case '[': 
                   //--- [h]alpha -> [h]alpha
                   ii=kk=i;
		   cflag=0; 
		   while(pch[ii] && !isblank(pch[ii])) { //pch[ii] can be '\0'
		     if(pch[ii]=='[') cflag -= 1;
		     else if(pch[ii]==']') { cflag += 1; kk=ii; }
		     ii++; 
		   }

		   if(!cflag && isalnum(pch[kk+1])) { //if alnum right after ')'
		     while(i<ii) str[u++]=pch[i++];
		     break;
		   }
		   //---

                   if(isblank(str[u-1])){
                      str[u++]=pch[i++];
                      if(!isblank(pch[i])){
                         str[u++]=' ';
                      }
                      openChFlag2.push_back(true);
                   }
                   else {
		      str[u++]=pch[i++];
                      openChFlag2.push_back(false);
		   }

                   break;

        case ']':  sz=openChFlag2.size(); 
	           if(sz>0 && openChFlag2[sz-1]){ //modified Jan-28-08
                      if(!isblank(str[u-1])){ 
                         str[u++]=' ';
                         str[u++]=pch[i++];
                      }
		      //---added (Jan-11-08): [BIV; ], -> [ BIV; ] ,
		      else if(!isblank(pch[i+1])){ 
			 str[u++]=pch[i++];
		      }
	              //---
                      if(!isblank(pch[i]))str[u++]=' ';
                   }
                   else str[u++]=pch[i++];

		   if(sz>0) openChFlag2.pop_back();
		   
                   break;
         default: str[u++]=pch[i++];
      }
   }
   while((u>0)&&(isblank(str[u-1])))u--;
   str[u]='\0';    
      
   while(str[j]){
      while(isblank(str[j]))j++;
      i=j;
      while((str[j])&&(!isblank(str[j])))j++;
      lst[k]=new char[j-i+1];
      strncpy(lst[k],str+i,j-i);
      lst[k][j-i]='\0';
      if(str[j]){
         k++;
         j++;
      }
   }
   num=k+1;
}


  void AbbrvE::clear(void){
    for ( int i=0; i<num; i++ ) {
      delete [] lst[i];
    }
    num=0;
  }

void AbbrvE::cleara(void){
   long i;
   for(i=0;i<numa;i++){
      delete [] abbl[i];
      delete [] abbs[i];
   }
   numa=0;
}

#if 0

//no space before and after abbs[] (because of using token)
int AbbrvE::Test(const char *str){
   long i,j,k;
   char b,c;

   if(strchr(str,'='))return(0);
   if(!strcmp(str,"author's transl"))return(0);
   if(!strcmp(str,"proceedings"))return(0);
   //---added (Jan-11-08) & (Apr 08)
   if(!strcmp(str,"see"))return(0);
   if(!strcmp(str,"and"))return(0);
   if(!strcmp(str,"comment"))return(0);
   if(!strcmp(str,"letter"))return(0);
   //---
   if((str[0]=='e')&&(str[1]=='g')){
      if(!(c=str[2])||(c=='.')||(c==','))return(0);
   }
   if((str[0]=='s')&&(str[1]=='e')&&(str[2]=='e')&&(((b=str[3])==' ')||(b==',')))return(0);
   if('p'==tolower(str[0])){
      if(strchr(str+1,'<'))return(0);
   }
   i=j=k=0;
   while((c=str[i])&&(c!=' ')){
      i++;
      if(isdigit(c))j++;
      if(isalpha(c))k++;
      if((i==j)&&(i==3))return(0);
   }
   if((i==j)||(k==0))return(0);
   else return(1);
} 

#endif

  bool AbbrvE::prefix_match( const char *str ) {
    size_t size = strlen(str);
    for ( int i = 0; i < prefix.size(); ++i ) {
      string& pre = prefix[i];
      if ( size > pre.size()  and
           0 == pre.compare( 0, pre.size(), str, pre.size() ) )
        return true;
    }
    return false;
  }

    
//no space before and after abbs[] (because of using token)
bool AbbrvE::Test(const char *str){

   if ( match.find(str) != match.end() ) return false;
   if ( prefix_match(str) ) return false;

   size_t length, letters, digits;
   length = letters = digits = 0;

   char c;
   while((c=str[length])&&(c!=' ')){
     length++;
     if ( isdigit(c) ) digits++;
     if ( isalpha(c) ) letters++;
     
     if( length==digits  and  length>=3 ) return false;
   }
   if ( digits == length ) return false;
   if ( letters <= 0 ) return false;

   return true;
} 

  void AbbrvE::setup_Test( void ) {

    match.insert("author's transl");
    match.insert("proceedings");
    match.insert("see");
    match.insert("and");
    match.insert("comment");
    match.insert("letter");
    match.insert("eg");

    prefix.push_back("=");
    prefix.push_back("eg.");
    prefix.push_back("eg,");
    prefix.push_back("see ");
    prefix.push_back("see,");
    prefix.push_back("p<");
    prefix.push_back("P<");

    // rules added in 2010
    match.insert("e.g.");
    match.insert("ie");
    match.insert("i.e.");
    match.insert("mean");
    match.insert("age");
    match.insert("std");
    match.insert("range");
    match.insert("young");
    match.insert("old");
    match.insert("male");
    match.insert("female");

  }

void AbbrvE::Rate(void){
   long i,j,k1,k2,k3,k4;
   char c,c1,c2,c3,c4;

   for(i=0;i<numa;i++)rt[i]=1;

   i=0;
   k1=k2=k3=k4=0;
   while(i<numa){
      if(k1&&(abbs[i][1]=='\0')){
         switch(c1){
            case 'a': if((c=abbs[i][0])=='b'){
                         rt[k1-1]=0;
                         rt[i]=0;
                         k1=i+1;
                         c1='b';
                      }
                      else if('a'==c){
                         rt[i]=0;
                      }
                      break;
            case 'b': if((c=abbs[i][0])=='c'){
                         rt[i]=0;
                         k1=i+1;
                         c1='c';
                      }
                      else if(('a'<=c)&&(c<=c1)){
                         rt[i]=0;
                      }
                      break;
            case 'c': if((c=abbs[i][0])=='d'){
                         rt[i]=0;
                         k1=i+1;
                         c1='d';
                      }
                      else if(('a'<=c)&&(c<=c1)){
                         rt[i]=0;
                      }
                      break;
            case 'd': 
            case 'e': if((c=abbs[i][0])=='e'){
                         rt[i]=0;
                         k1=i+1;
                         c1='e';
                      }
                      else if(('a'<=c)&&(c<=c1)){
                         rt[i]=0;
                      }
                      break;
         }
      }
      if(k2&&(abbs[i][1]=='\0')){
         switch(c2){
            case 'A': if((c=abbs[i][0])=='B'){
                         rt[k2-1]=0;
                         rt[i]=0;
                         k2=i+1;
                         c2='B';
                      }
                      else if(('A'<=c)&&(c<=c2)){
                         rt[i]=0;
                      }
                      break;
            case 'B': if((c=abbs[i][0])=='C'){
                         rt[i]=0;
                         k2=i+1;
                         c2='B';
                      }
                      else if(('A'<=c)&&(c<=c2)){
                         rt[i]=0;
                      }
                      break;
            case 'C': if((c=abbs[i][0])=='D'){
                         rt[i]=0;
                         k2=i+1;
                         c2='D';
                      }
                      else if(('A'<=c)&&(c<=c2)){
                         rt[i]=0;
                      }
                      break;
            case 'D': 
            case 'E': if((c=abbs[i][0])=='E'){
                         rt[i]=0;
                         k2=i+1;
                         c2='E';
                      }
                      else if(('A'<=c)&&(c<=c2)){
                         rt[i]=0;
                      }
                      break;
         }
      }
      if(k3){
         switch(c3){
            case 'i': if(!strcmp(abbs[i],"ii")){
                         rt[k3-1]=0;
                         rt[i]=0;
                         k3=i+1;
                         c3='2';
                      }
                      else if(!strcmp(abbs[i],"i")){
                         rt[i]=0;
                      }
                      break;
            case '2': if(!strcmp(abbs[i],"iii")){
                         rt[i]=0;
                         k3=i+1;
                         c3='3';
                      }
                      else if(!strcmp(abbs[i],"i")||!strcmp(abbs[i],"ii")){
                         rt[i]=0;
                      }
                      break;
            case '3': if(!strcmp(abbs[i],"iv")){
                         rt[i]=0;
                         k3=i+1;
                         c3='4';
                      }
                      else if(!strcmp(abbs[i],"i")||!strcmp(abbs[i],"ii")){
                         rt[i]=0;
                      }
                      else if(!strcmp(abbs[i],"iii")){
                         rt[i]=0;
                      }
                      break;
            case '4': 
            case '5': if(!strcmp(abbs[i],"v")){
                         rt[i]=0;
                         k3=i+1;
                         c3='5';
                      }
                      else if(!strcmp(abbs[i],"i")||!strcmp(abbs[i],"ii")){
                         rt[i]=0;
                      }
                      else if(!strcmp(abbs[i],"iii")||!strcmp(abbs[i],"iv")){
                         rt[i]=0;
                      }
                      else if(!strcmp(abbs[i],"v")){
                         rt[i]=0;
                      }
                      break;
         }
      }
      if(k4){
         switch(c4){
            case 'I': if(!strcmp(abbs[i],"II")){
                         rt[k4-1]=0;
                         rt[i]=0;
                         k4=i+1;
                         c4='2';
                      }
                      else if(!strcmp(abbs[i],"I")){
                         rt[i]=0;
                      }
                      break;
            case '2': if(!strcmp(abbs[i],"III")){
                         rt[i]=0;
                         k4=i+1;
                         c4='3';
                      }
                      else if(!strcmp(abbs[i],"I")||!strcmp(abbs[i],"II")){
                         rt[i]=0;
                      }
                      break;
            case '3': if(!strcmp(abbs[i],"IV")){
                         rt[i]=0;
                         k4=i+1;
                         c4='4';
                      }
                      else if(!strcmp(abbs[i],"I")||!strcmp(abbs[i],"II")){
                         rt[i]=0;
                      }
                      else if(!strcmp(abbs[i],"III")){
                         rt[i]=0;
                      }
                      break;
            case '4': 
            case '5': if(!strcmp(abbs[i],"V")){
                         rt[i]=0;
                         k4=i+1;
                         c4='5';
                      }
                      else if(!strcmp(abbs[i],"I")||!strcmp(abbs[i],"II")){
                         rt[i]=0;
                      }
                      else if(!strcmp(abbs[i],"III")||!strcmp(abbs[i],"IV")){
                         rt[i]=0;
                      }
                      else if(!strcmp(abbs[i],"V")){
                         rt[i]=0;
                      }
                      break;
         }
      }
      if(abbs[i][1]=='\0'){
         if(!k1&&((c=abbs[i][0])=='a')){k1=i+1;c1='a';}
         else if(!k2&&(c=='A')){k2=i+1;c2='A';}
         else if(!k3&&(c=='i')){k3=i+1;c3='i';}
         else if(!k4&&(c=='I')){k4=i+1;c4='I';}
      }
      i++;
   }
}

void AbbrvE::Proc(char *pxh){
   long i,j;
   char *pch,*ptr;
   pMt->segment(pxh);
   for(i=0;i<pMt->sent.size();i++){
     Extract2( (pMt->sent[i]).c_str() );
   }
   
   Rate();
   j=0;
   for(i=0;i<numa;i++){
      if(rt[i]){
         if(j<i){
            pch=abbl[i];
            if(ptr=strchr(pch,'|')){
               *ptr='/';
               ptr++;
               while(ptr=strchr(pch,'|')){
                  *ptr='/';
                  ptr++;
               }
            }
            abbl[j]=pch;
            pch=abbs[i];
            if(ptr=strchr(pch,'|')){
               *ptr='/';
               ptr++;
               while(ptr=strchr(pch,'|')){
                  *ptr='/';
                  ptr++;
               }
            }
            abbs[j]=pch;
            nt[j]=nt[i];
         }
         j++;
      }
      else {
         delete [] abbl[i];
         delete [] abbs[i];
      }
   }

   numa=j;
}

}
