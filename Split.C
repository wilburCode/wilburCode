#include <Btree.h>
#include <Split.h>
#include <vector>
#include <string>

namespace iret {

Split::Split(){
   word_space=10000;
   lst=new char*[word_space];
   xnam=new char[10*word_space];
   ynam=new char[10*word_space];
} 

Split::Split(int wrd_spc){
   word_space=wrd_spc;
   lst=new char*[word_space];
   xnam=new char[10*word_space];
   ynam=new char[10*word_space];
} 
   
Split::~Split(){
   if(lst)delete [] lst;
   delete [] xnam;
   delete [] ynam;
}  

void Split::token(const char *str){
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
   num=k;
}

void Split::token_tab(const char *str){
  int i,j=0,k=0;
  while(str[j]){
    while(isblank(str[j]))j++;
    i=j;
    while((str[j])&&(str[j]!='\t'))j++;
    if(j>i){
      lst[k]=new char[j-i+1];
      strncpy(lst[k],str+i,j-i);
      lst[k][j-i]='\0';
      k++;
    }
  }
  num=k;
}

void Split::token_chr(char ch,const char *str){
  int i,j=0,k=0;
  while(str[j]){
    while(isblank(str[j])||(str[j]==ch))j++;
    i=j;
    while((str[j])&&(str[j]!=ch))j++;
    if(j>i){
      lst[k]=new char[j-i+1];
      strncpy(lst[k],str+i,j-i);
      lst[k][j-i]='\0';
      k++;
      if(str[j])j++;
    }
  }
  num=k;
}

void Split::tokenS_chr(char ca,const char *str){
   char *pt,ch;
   int s1=0,s2=0,s3=0,s4=0;
   strcpy(xnam,str);
   pt=xnam;
   switch(*pt){
      case '\n': *pt=' '; break;
      case '-': *pt=' '; break;
      case '.': *pt=' '; break;
      case ',': *pt=' '; break;
      case ':': *pt=' '; break;
      case '?': *pt=' '; break;
      case '!': *pt=' '; break;
      case ';': *pt=' '; break;
      case '"': *pt=' '; break;
      case '(': *pt=' '; s1++;break;
      case '<': *pt=' '; s2++;break;
      case '{': *pt=' '; s3++;break;
      case '[': *pt=' '; s4++;break;
   }
   pt++;
   while(*pt){
      switch(*pt){
         case '\n': *pt=' '; break;
         case '-': *pt=' '; break;
         case '"': *pt=' '; break;
         case '.': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ',': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ':': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '?': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '!': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ';': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ')': if(s1){*pt=' ';s1--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '>': if(s2){*pt=' ';s2--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '}': if(s3){*pt=' ';s3--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ']': if(s4){*pt=' ';s4--;} 
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '(': if((*(pt-1)=' ')||s1){s1++;*pt=' ';} break;
         case '<': if((*(pt-1)=' ')||s2){s2++;*pt=' ';} break;
         case '{': if((*(pt-1)=' ')||s3){s3++;*pt=' ';} break;
         case '[': if((*(pt-1)=' ')||s4){s4++;*pt=' ';} break;
      }
      pt++;
   }
   token_chr(ca,xnam);
}

void Split::puncChunk(const char *str){
   char ch;
   int i,j;
   strcpy(xnam,str);
   i=j=0;
   switch(xnam[0]){
      case '.': xnam[0]=' '; break;
      case ',': xnam[0]=' '; break;
      case ':': xnam[0]=' '; break;
      case '?': xnam[0]=' '; break;
      case '!': xnam[0]=' '; break;
      case ';': xnam[0]=' '; break;
   }
   j++;
   while(xnam[j]){
      switch(xnam[j]){
         case '.': if(((ch=xnam[j+1])==' ')||(ch=='\0')){
                      if(j>i+2){
                         xnam[j]='\0'; 
                         lst[num]=new char[j-i+1];
                         strcpy(lst[num++],xnam+i);
                         if(ch)j+=2;
                         else j++;
                         i=j;
                      }
                      else if(j>i){
                         if(ch)j+=2;
                         else j++;
                      }
                      else {i++;j++;}
                   }
                   else j++;
                   break;
         case ',': if(((ch=xnam[j+1])==' ')||(ch=='\0')){
                      if(j>i){
                         xnam[j]='\0'; 
                         lst[num]=new char[j-i+1];
                         strcpy(lst[num++],xnam+i);
                         if(ch)j+=2;
                         else j++;
                         i=j;
                      }
                      else {i++;j++;}
                   }
                   else j++;
                   break;
         case ':': if(((ch=xnam[j+1])==' ')||(ch=='\0')){
                      if(j>i){
                         xnam[j]='\0'; 
                         lst[num]=new char[j-i+1];
                         strcpy(lst[num++],xnam+i);
                         if(ch)j+=2;
                         else j++;
                         i=j;
                      }
                      else {i++;j++;}
                   }
                   else j++;
                   break;
         case '?': if(((ch=xnam[j+1])==' ')||(ch=='\0')){
                      if(j>i){
                         xnam[j]='\0'; 
                         lst[num]=new char[j-i+1];
                         strcpy(lst[num++],xnam+i);
                         if(ch)j+=2;
                         else j++;
                         i=j;
                      }
                      else {i++;j++;}
                   }
                   else j++;
                   break;
         case '!': if(((ch=xnam[j+1])==' ')||(ch=='\0')){
                      if(j>i){
                         xnam[j]='\0'; 
                         lst[num]=new char[j-i+1];
                         strcpy(lst[num++],xnam+i);
                         if(ch)j+=2;
                         else j++;
                         i=j;
                      }
                      else {i++;j++;}
                   }
                   else j++;
                   break;
         case ';': if(((ch=xnam[j+1])==' ')||(ch=='\0')){
                      if(j>i){
                         xnam[j]='\0'; 
                         lst[num]=new char[j-i+1];
                         strcpy(lst[num++],xnam+i);
                         if(ch)j+=2;
                         else j++;
                         i=j;
                      }
                      else {i++;j++;}
                   }
                   else j++;
                   break;
         default: j++;
      }
   }
   if(j>i){
      lst[num]=new char[j-i+1];
      strcpy(lst[num++],xnam+i);
   }
}

void Split::loadTrecn(void){
   string sx,sy;
   ifstream fin("/panfs/pan1.be-md.ncbi.nlm.nih.gov/cbb-wilbur/wilbur/CentosCPP64/lib/FIXED_DATA/Token.recognition",ios::in);
   Trecn.clear();
   while(fin >> sx >> sy){
      Trecn.insert(std::pair<string,string>(sx,sy));
   }
}

long Split::nNum(char *str){
   long u,j,k;
   char anam[10000];

   j=0;
   while(isdigit(str[j]))j++;
   if((!j)||(str[0]=='0'))return(0);
   else if(!str[j]){
      str_long(str,u);
      return(u);
   }
   else if(ispunct(str[j])){
      strncpy(anam,str,j-1);
      anam[j]='\0';
      str_long(str,u);
      return(u);
   }
   else return(0);
}
     
void Split::tokenQuery(const char *str,vector<string> &vx){
   char ch,*pch,*ptr,anam[10000];
   int i=0,j=0,k=0,ip,flag;
   long u;
   strcpy(xnam,str);
   string sx,sy,sz;
   vx.clear();

   while(xnam[j]){
      while(xnam[j]&&(isblank(xnam[j])||(strchr(":;,",xnam[j]))))j++;
      if(xnam[j]=='('){
         i=j++;ip=1;
         while((xnam[j])&&((xnam[j]!=')')||(ip>1))){
            if(xnam[j]=='(')ip++;
            else if(xnam[j]==')')ip--;
            j++;
         }
         if(xnam[j]&&(j>i+1)){
            lst[k]=new char[j-i];
            strncpy(lst[k],xnam+i+1,j-i-1);
            lst[k][j-i-1]='\0';
            sx="parenthetical";vx.push_back(sx);
            j++;k++;
         }
         else if(j>i+1){
            lst[k]=new char[j-i];
            strcpy(lst[k],xnam+i+1);
            sx="parenthetical_left";vx.push_back(sx);
            k++;
         }
         else if(xnam[j])j++;
      }
      else if(xnam[j]=='['){
         i=j++;
         while((xnam[j])&&(xnam[j]!=']')){
            j++;
         }
         if(xnam[j]&&(j>i+1)){
            lst[k]=new char[j-i+2];
            strncpy(lst[k],xnam+i,j-i+1);
            lst[k][j-i+1]='\0';
            sx="Tag";vx.push_back(sx);
            j++;k++;
         }
         else if(j>i+1){
            lst[k]=new char[j-i+1];
            strcpy(lst[k],xnam+i);
            sx="Tag_left";vx.push_back(sx);
            k++;
         }
         else if(xnam[j])j++;
      }
      else if(xnam[j]=='"'){
         i=j++;
         while((xnam[j])&&(xnam[j]!='"')){
            j++;
         }
         if(xnam[j]&&(j>i+1)){
            lst[k]=new char[j-i];
            strncpy(lst[k],xnam+i+1,j-i-1);
            lst[k][j-i-1]='\0';
            sx="quote";vx.push_back(sx);
            j++;k++;
         }
         else if(j>i+1){
            lst[k]=new char[j-i];
            strcpy(lst[k],xnam+i+1);
            sx="quote_left";vx.push_back(sx);
            k++;
         }
         else if(xnam[j])j++;
      }
      else {
         i=j;
         flag=0;
         while((xnam[j])&&((!isblank(xnam[j]))&&(!strchr("([",xnam[j])))){
            switch(xnam[j]){
               case '.': if((j>i)&&(!isdigit(xnam[j-1]))){
                            if(isdigit(xnam[j+1]))flag=1;
                            else if(xnam[j+1]==',')flag=2;
                         }
                         break;
               case ':': if(j>i)flag=1;
                         break;
               case ';': if(j>i)flag=1;
                         break;
               case ',': if(j>i)flag=1;
                         break;
            }
            if(flag==2)xnam[j]=';';
            j++;
            if(flag)break;
         }
         if(xnam[j]=='['){
            while(xnam[j]&&(xnam[j]!=']'))j++;
            if(xnam[j]==']'){
               lst[k]=new char[j-i+2];
               strncpy(lst[k],xnam+i,j-i+1);
               lst[k][j-i+1]='\0';
               sx="Tagged";vx.push_back(sx);
               j++;k++;
            }
            else {
               lst[k]=new char[j-i+1];
               strcpy(lst[k],xnam+i);
               sx="Tagged_left";vx.push_back(sx);
               k++;
            }
         }      
         else if(xnam[j]&&((j>i+1)||((j>i)&&(!strchr(".,;:",xnam[i]))))){
            lst[k]=new char[j-i+1];
            strncpy(lst[k],xnam+i,j-i);
            lst[k][j-i]='\0';
            sx="NULL";vx.push_back(sx);
            if(isblank(xnam[j]))j++;k++;
         }
         else if((j>i+1)||((j>i)&&(!strchr(".,;:",xnam[i])))){
            lst[k]=new char[j-i+1];
            strcpy(lst[k],xnam+i);
            sx="NULL";vx.push_back(sx);
            k++;
         }
      }
   }
   num=k;
   //typing the tokens
   map<string,string>::iterator mb,me;
   me=Trecn.end();
   for(i=0;i<num;i++){
      sx=vx[i];
      if(sx!="NULL")continue;
      sy=lst[i];
      if((mb=Trecn.find(sy))!=me){
         vx[i]=mb->second;
         if(vx[i]=="PageIndicator"){
            if((i+1<num)&&(isdigit(lst[i+1][0])||isdigit(lst[i+1][1])))continue;
         }
         else if(vx[i]=="VolumeIndicator"){
            if((i+1<num)&&isdigit(lst[i+1][0]))continue;
         }
         else continue;
      }
      if(ispunct(sy.back())){
         strcpy(anam,lst[i]);
         anam[sy.size()-1]='\0';
         sz=anam;
         if(((mb=Trecn.find(sz))!=me)&&(mb->second=="Month")){
            vx[i]="MonthP";
            continue;
         }
      }
      if(isdigit(lst[i][0])||(isdigit(lst[i][1])&&(ptr=strchr(lst[i],'-'))&&((*(ptr+1)==lst[i][0])||isdigit(*(ptr+1))))){
         if(strchr(lst[i],'['))vx[i]="NumericT";
         else vx[i]="Numeric";
      }
      else {
         j=0;
         while(lst[i][j]&&isalpha(lst[i][j]))j++;
         if(!lst[i][j]){
            k=0;
            while(lst[i][k]&&isupper(lst[i][k]))k++;
            if(!lst[i][k])vx[i]="Abbr";
            else if(k>0)vx[i]="W";
            else vx[i]="w";
         } 
         else if(lst[i][j]&&strchr(":;,.",lst[i][j])&&(!lst[i][j+1])){
            k=0;
            while(isupper(lst[i][k]))k++;
            switch(lst[i][j]){
               case '.': if(k==j)vx[i]="Abbr.";
                         else if(k>0)vx[i]="W.";
                         else vx[i]="w.";
                         break;
               case ':': if(k==j)vx[i]="Abbr:";
                         else if(k>0)vx[i]="W:";
                         else vx[i]="w:";
                         break;
               case ';': if(k==j)vx[i]="Abbr;";
                         else if(k>0)vx[i]="W;";
                         else vx[i]="w;";
                         break;
               case ',': if(k==j)vx[i]="Abbr,";
                         else if(k>0)vx[i]="W,";
                         else vx[i]="w,";
                         break;
            }
         }
         else if(j>0){
            k=0;
            while(lst[i][k]&&isascii(lst[i][k])&&(lst[i][k]!='['))k++;
            if(!lst[i][k])vx[i]="ASCII";
            else if((k>0)&&(lst[i][k]=='['))vx[i]="Tagged";
            else if(lst[i][k]=='[')vx[i]="Tag";
            else vx[i]="NonASCII";
         }   
         else if(lst[i][j]=='[')vx[i]="Tag";
         else vx[i]="NonAlpha";
      }
   }
   //type logic
   int iy,iv,id,ii;
   iy=id=iv=ii=ip=-9;
   for(i=0;i<num;i++){
      sx=vx[i];
      pch=lst[i];
      if(sx=="PageIndicator"){
         if((i+1<num)&&(vx[i+1]=="Numeric")){
            vx[i+1]="Page";
            if(i&&(vx[i-1]=="Numeric")){
               if((i>1)&&(vx[i-2]=="Numeric")){
                  vx[i-2]="Volume";iv=i-2;
                  vx[i-1]="Issue";ii=i-1;
               }
               else {
                  vx[i-1]="Volume";
                  iv=i-1;
               }
            }
            i++;
         }
         else if(i+1<num){
            vx[i+1]="Page";
            ip=i+1;
            i++;
         }
      }
      else if(sx=="VolumeIndicator"){
         if((i+1<num)&&(vx[i+1]=="Numeric")){
            vx[i+1]="Volume";
            i++;
         }
         else if(i+1<num){
            vx[i+1]="Volume";
            iv=i+1;
            i++;
         }
      }
      else if(sx=="Numeric"){
         j=0;
         while(isdigit(pch[j]))j++;
         if(!j){
            vx[i]="Page";
            ip=i;
         }
         else if((!pch[j])&&(j<5)){
            str_long(pch,u);
            if((u>1900)&&(u<2018)&&((iy==-9)||(i-iy>3))){
               vx[i]="Year";
               iy=i;
            }
            else if(i&&(vx[i-1]=="Month")){
               if(u&&(u<32)){
                  vx[i]="Day";
                  id=i;
                  if((i+1<num)&&(vx[i+1]=="Numeric")){vx[i+1]="Volume";iv=i+1;i++;}
               }
               else {
                  vx[i]="Volume";
                  iv=i;
               }
            }
            else if(iy==i-1){
               vx[i]="Volume";
               iv=i;
            }
            else if((iv==i-1)||(ii==i-1)){
               vx[i]="Page";
               ip=i;
            }
         }
         else if(ispunct(pch[j])){
            strncpy(anam,pch,j);
            anam[j]='\0';
            str_long(anam,u);
            if((u>1900)&&(u<2018)&&((iy==-9)||(i-iy>3))){
               if(pch[j]==':'){
                  vx[i]="YearR";
               }
               else if(pch[j]=='/'){
                  vx[i]="Date";
               }
               else {
                  vx[i]="Year";
                  iy=i;
               }
            }
            else if(i&&(vx[i-1]=="Month")){
               if(u&&(u<32)){
                  vx[i]="Day";
                  id=i;
                  if((i+1<num)&&(vx[i+1]=="Numeric")&&(u=nNum(lst[i+1]))){
                     if(u<1900){vx[i+1]="Volume";iv=i+1;}
                     else if(u<2018){vx[i+1]="Year";iy=i+1;}
                     i++;
                  }
               }
               else {
                  vx[i]="Volume";
                  iv=i;
               }
            }
            else if(pch[j]==':'){
               if(iv==i-1)vx[i]="Issue";
               else vx[i]="Volume";
               if((i+1<num)&&(vx[i+1]=="Numeric")){
                  vx[i+1]="Page";
                  ip=i+1;
                  i++;
               }
            }
            else if(pch[j]=='-'){
               if(strchr(pch,'['))vx[i]="Tagged";
               else if((pch[0]!='0')&&isdigit(pch[j+1])){
                  if(i&&(ip==i-1)){vx[i-1]="Issue";ii=i-1;}
                  else if(i&&(vx[i-1]=="Numeric")){vx[i-1]="Volume";iv=i-1;}
                  vx[i]="Page";
                  ip=i;
               }
            }
            else if(iy==i-1){
               vx[i]="Volume";
               iv=i;
            }
            else if((iv==i-1)||(ii==i-1)){
               vx[i]="Page";
               ip=i;
            }
         }
      }
      else if(sx=="parenthetical"){
         j=0;
         while(isdigit(pch[j]))j++;
         if((!pch[j])&&(j<5)){
            str_long(pch,u);
            if((u>1900)&&(u<2018)&&((iy==-9)||(i-iy>3))){
               vx[i]="Year";
               iy=i;
            }
            else {
               if(i&&(vx[i-1]=="Numeric")){
                  vx[i-1]="Volume";
                  iv=i-1;
                  vx[i]="Issue";
                  ii=i;
               }
               else if(i&&((vx[i-1]=="Volume")||(vx[i-1]=="Year")||(vx[i-1]=="Day"))){
                  vx[i]="Issue";
                  ii=i;
               }
            }
         }
      }
   }
   
}

void Split::token_lower(const char *str){
   char *pt;
   strcpy(xnam,str);
   pt=xnam;
   while(*pt){
      *pt=tolower(*pt);
      pt++;
   }
   token(xnam);
}

void Split::tokenS_lower(const char *str){
   char *pt,ch;
   int s1=0,s2=0,s3=0,s4=0;
   strcpy(xnam,str);
   pt=xnam;
   switch(*pt){
      case '\n': *pt=' '; break;
      case '-': *pt=' '; break;
      case '"': *pt=' '; break;
      case '.': *pt=' '; break;
      case ',': *pt=' '; break;
      case ':': *pt=' '; break;
      case '?': *pt=' '; break;
      case '!': *pt=' '; break;
      case ';': *pt=' '; break;
      case '(': *pt=' '; s1++;break;
      case '<': *pt=' '; s2++;break;
      case '{': *pt=' '; s3++;break;
      case '[': *pt=' '; s4++;break;
      default: *pt=tolower(*pt);
   }
   pt++;
   while(*pt){
      switch(*pt){
         case '\n': *pt=' '; break;
         case '-': *pt=' '; break;
         case '"': *pt=' '; break;
         case '.': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0'; 
                   break;
         case ',': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0'; 
                   break;
         case ':': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0'; 
                   break;
         case '?': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0'; 
                   break;
         case '!': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0'; 
                   break;
         case ';': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0'; 
                   break;
         case ')': if(s1){*pt=' ';s1--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '>': if(s2){*pt=' ';s2--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '}': if(s3){*pt=' ';s3--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ']': if(s4){*pt=' ';s4--;} 
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '(': if((*(pt-1)==' ')||s1){s1++;*pt=' ';} break;
         case '<': if((*(pt-1)==' ')||s2){s2++;*pt=' ';} break;
         case '{': if((*(pt-1)==' ')||s3){s3++;*pt=' ';} break;
         case '[': if((*(pt-1)==' ')||s4){s4++;*pt=' ';} break;
         default: *pt=tolower(*pt);
      }
      pt++;
   }
   token(xnam);
}

void Split::tokenS(const char *str){
   char *pt,ch;
   int s1=0,s2=0,s3=0,s4=0;
   strcpy(xnam,str);
   pt=xnam;
   switch(*pt){
      case '\n': *pt=' '; break;
      case '-': *pt=' '; break;
      case '"': *pt=' '; break;
      case '.': *pt=' '; break;
      case ',': *pt=' '; break;
      case ':': *pt=' '; break;
      case '?': *pt=' '; break;
      case '!': *pt=' '; break;
      case ';': *pt=' '; break;
      case '(': *pt=' '; s1++;break;
      case '<': *pt=' '; s2++;break;
      case '{': *pt=' '; s3++;break;
      case '[': *pt=' '; s4++;break;
   }
   pt++;
   while(*pt){
      switch(*pt){
         case '\n': *pt=' '; break;
         case '-': *pt=' '; break;
         case '"': *pt=' '; break;
         case '.': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ',': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ':': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '?': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '!': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ';': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ')': if(s1){*pt=' ';s1--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '>': if(s2){*pt=' ';s2--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '}': if(s3){*pt=' ';s3--;}
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ']': if(s4){*pt=' ';s4--;} 
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '(': if((*(pt-1)==' ')||s1){s1++;*pt=' ';} break;
         case '<': if((*(pt-1)==' ')||s2){s2++;*pt=' ';} break;
         case '{': if((*(pt-1)==' ')||s3){s3++;*pt=' ';} break;
         case '[': if((*(pt-1)==' ')||s4){s4++;*pt=' ';} break;
      }
      pt++;
   }
   token(xnam);
}

void Split::token_tab_lower(const char *str){
   char *pt;
   strcpy(xnam,str);
   pt=xnam;
   while(*pt){
      *pt=tolower(*pt);
      pt++;
   }
   token_tab(xnam);
}

void Split::token_chr_lower(char ch,const char *str){
   char *pt,ca=tolower(ch);
   strcpy(xnam,str);
   pt=xnam;
   while(*pt){
      *pt=tolower(*pt);
      pt++;
   }
   token_chr(ca,xnam);
}

void Split::tokenS_chr_lower(char ca,const char *str){
   char *pt,ch,cb;
   int s1=0,s2=0,s3=0,s4=0;
   strcpy(xnam,str);
   pt=xnam;
   switch(*pt){
      case '\n': *pt=' '; break;
      case '-': *pt=' '; break;
      case '"': *pt=' '; break;
      case '.': *pt=' '; break;
      case ',': *pt=' '; break;
      case ':': *pt=' '; break;
      case '?': *pt=' '; break;
      case '!': *pt=' '; break;
      case ';': *pt=' '; break;
      case '(': *pt=' '; s1++;break;
      case '<': *pt=' '; s2++;break;
      case '{': *pt=' '; s3++;break;
      case '[': *pt=' '; s4++;break;
      default: *pt=tolower(*pt);
   }
   pt++;
   while(*pt){
      switch(*pt){
         case '\n': *pt=' '; break;
         case '-': *pt=' '; break;
         case '"': *pt=' '; break;
         case '.': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ',': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ':': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '?': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '!': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ';': if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ')': if(s1){*pt=' ';s1--;} 
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '>': if(s2){*pt=' ';s2--;} 
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '}': if(s3){*pt=' ';s3--;} 
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case ']': if(s4){*pt=' ';s4--;} 
                   else if((ch=*(pt+1))==' ')*pt=' ';
                   else if(ch=='\0')*pt='\0';
                   break;
         case '(': if((*(pt-1)==' ')||s1){s1++;*pt=' ';} break;
         case '<': if((*(pt-1)==' ')||s2){s2++;*pt=' ';} break;
         case '{': if((*(pt-1)==' ')||s3){s3++;*pt=' ';} break;
         case '[': if((*(pt-1)==' ')||s4){s4++;*pt=' ';} break;
         default: *pt=tolower(*pt);
      }
      pt++;
   }
   cb=tolower(ca);
   token_chr(cb,xnam);
}

void Split::couple(char *buf,char *s1,char *s2){
   strcpy(buf,s1);
   strcat(buf," ");
   strcat(buf,s2);
}

void Split::lower_lst(void){
   int i;
   char *pt;
   for(i=0;i<num;i++){
      pt=lst[i];
      while(*pt){
         *pt=tolower(*pt);
         pt++;
      }
   }
}

void Split::clear(void){
   int i;
   for(i=0;i<num;i++){
      delete [] lst[i];
   }
   num=0;
}

void Split::clean(const char *pch){
   long i,j,k;
   char ptr;
   strcpy(xnam,pch);
   i=strlen(xnam);
   if((xnam[i-2]=='!')&&(xnam[i-3]=='!')){
      xnam[i-3]='\0';
      i-=3;
   }
   else if((xnam[i-1]=='!')||(xnam[i-1]=='*')){
      xnam[i-1]='\0';
      i--;
   }
   j=0;
   for(j=0;j<i;j++){
      if((xnam[j]=='!')||(xnam[j]=='*'))xnam[j]=' ';
   }
   token(xnam);
}
   
char* Split::assembl(void){
   long i,j,k;

   if(num>0)strcpy(xnam,lst[0]);
   else return(NULL);
   for(i=1;i<num;i++){
      strcat(xnam," ");
      strcat(xnam,lst[i]);
   }
   i=strlen(xnam);
   char *pzz=new char[i+1];
   strcpy(pzz,xnam);
   return(pzz);
}

char* Split::assembl_tab(void){
   long i,j,k;

   if(num>0)strcpy(xnam,lst[0]);
   else return(NULL);
   for(i=1;i<num;i++){
      strcat(xnam,"\t");
      strcat(xnam,lst[i]);
   }
   i=strlen(xnam);
   char *pzz=new char[i+1];
   strcpy(pzz,xnam);
   return(pzz);
}

char* Split::pair(long m){
   long i;

   if(num>m+1){
      strcpy(xnam,lst[m]);
      strcat(xnam," ");
      strcat(xnam,lst[m+1]);
      i=strlen(xnam);
      char *pzz=new char[i+1];
      strcpy(pzz,xnam);
      return(pzz);
   }
   else return(NULL);
}

char* Split::segment(long k,long m){
   long len=0,i,j;

   for(i=k;i<k+m;i++){
      len+=strlen(lst[i])+1;
   }
   char *ux=new char[len];
   strcpy(ux,lst[k]);
   len=strlen(lst[k]);
   for(i=k+1;i<k+m;i++){
      ux[len++]=' ';
      strcpy(ux+len,lst[i]);
      len+=strlen(lst[i]);
   }
   ux[len]='\0';
   return(ux);
}
   
void Split::order(void){
  long k, j, ir, i;
  char* rra;

  if(num<=1)return;

  k=(num>>1);
  ir=num-1;
  for(;;) {
    if(k>0) {
      rra=lst[--k];
    }
    else {
      rra=lst[ir];
      lst[ir] = lst[0];
      if(--ir ==0) {
        lst[0]=rra;
        return;
      }
    }
    i=k;
    j=((k+1)<<1)-1;
    while(j<=ir) {
      if(j<ir && (strcmp(lst[j],lst[j+1])<0)) ++j;
      if(strcmp(rra,lst[j])<0) {
        lst[i]=lst[j];
        j +=(i=j)+1;
      }
      else j=ir+1;
    }
    lst[i]=rra;
  }
}

char* Split::leftMod(long k,string st){
   if((k<0)||(k>=num))return(NULL);
   string sx=st+lst[k];
   delete [] lst[k];
   lst[k]=new char[sx.size()+1];
   strcpy(lst[k],sx.c_str());
   return(lst[k]);
}
   
char* Split::rightMod(long k,string st){
   if((k<0)||(k>=num))return(NULL);
   string sx=lst[k]+st;
   delete [] lst[k];
   lst[k]=new char[sx.size()+1];
   strcpy(lst[k],sx.c_str());
   return(lst[k]);
}

}
