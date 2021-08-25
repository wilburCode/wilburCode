#ifndef DSTOR_H
#define DSTOR_H
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <FBase.h>
#include <Bnum.h>
#include <runn.h>
using namespace std;
namespace iret {

class DSpan {
public:
  DSpan(void){bln=10000;len=0;buf=new char[10000];}
  ~DSpan(void){if(buf!=NULL)delete [] buf;}
  void showData(void); //for debugging
  long id;
  long date;
  int len;
  char * buf;
  int bln; //Current length of buf
};

class DSlice : public FBase {
   public:
      DSlice(const char *nam,int nm); 
     ~DSlice();
      //Add to the store
      bool gopen_add(void); //Opens files ready to add records to the store
         //No index files are used and no check is done when adding records.
         //They are simply appended to the appropriate file.
         //returns true of opened, false otherwise
      bool add(DSpan &Ds); //Writes the data to
         //appropriate file. Returns true if success. If fails returns false
      void gclose_add(void); //Close file.
      //Create indices
      //create index files ".ids", ".dat", ".adr" all longs in binary 
      //Also creates file ".siz" for total number of records in slice
      bool createIndex(void); //returns 1 for success
         //otherwise 0. 
      bool createStatus(void); //writes ".sts" file of type char
         //reads the ".ids" and ".dat" files and bases on them
         //for any id all but most recent date are marked obsolete ('o')
         //Most recent record marked with 'a' for active
         //Returns number of obsolete records in the file 
      bool updateIndex(void); //Returns true if successful
         //Uses current index files to access last record before update. Then reads 
         //through new records appending to all index files. Finally calls createStatus
         //if records have been added. If no records added does nothing
      bool checkStatus(long &cur,long &obs); //Reads the status array and reports status to stdout
         //If no sts array exists returns false, else true
         //cur contains number of current records in slice and obs number of obsolete
      //Read from store
      bool gopen_read(void); //Prepares for reading
         //reads in the id files omitting id for any record marked as obsolete ('o')
         //in ".sts" file. Re-orders these arrays so ids are in increasing order for each 
         //file. Uses an additional array to keep the original index for each record and 
         //memory maps the ".adr" files so can access records based on this.
      bool gopen_read(long date_thr); //As gopen_read() except uses the record with the 
         //latest date that is <=date_thr for a given id. If  none leaves that id out
      bool exists(long id); //Checks if pmid in this slice and returns true if is, else false
      long date(long id); //Checks if pmid in this slice and returns date if is, else zero
      bool read(long id,DSpan &Ds);
         //On call len must be the length in chars of buf and buf must come from new.
         //on successful return length will contain number of bytes read in. If len is
         //too small, buf will be deleted and reallocated to allow the read.
         //Does binary search to find correct file and then binary search in the id array
         //for that file to find record. If not found returns false or if read operation fails
         //returns false. At end of read compares the tail count with head count and if not 
         //equal it returns a false. If tail and head count agree then returns a true 
         //and this represents a successful read.
      void reset(void); //Sets cind to -1 
      long getNext(void); //Returns the next listed record id+1. Or 0 if fails (no more records).
      long date(void); //returns date for pmid last returned by getNext()
      bool read(DSpan &Ds);//read after a successful call to getNext().
         //Returns false if read fails
         //Returns false if head and tail counts not same, true if are the same.
      void gclose_read(void); //Deletes arrays and closes files
      bool gopen_markStatus(void); //Opens the "sts" file for random writes
         //Mark functions require the gopen_read function be called before they can be used
      bool markStatus(long id,char status); //Uses the state with gopen_read() to find 
         //active record if it exists in the store and mark it with status char. 
         //Does this by random access write in
         // "sts" file and changing the value. Results will be lost if records 
         //added and indices rebuilt before a clean operation. If finds record and marks 
         //returns true, else falsee 
      bool markStatus(char status); //Uses the state with gopen_read() and current active record 
         //performs random access write in "sts" file to mark status  
      void gclose_markStatus(void); //Closes the "sts" file
      //Clean slice by removing all obsolete records
      bool cleanSlice(void);//cleans all files. Assumes some form of call to gopen_read
         //Runs through the array of ids and reads and writes in id order.  
         //Writing is into a temporary file and when done deletes original file and 
         //renames the new file to the correct name. Also deletes index files so 
         //indices can be recreated.
      //Special functions for more rapid access to the data
      bool create_map(const char *prvt); //In path_*prvt is the path where these 
         //files will be stored
      bool create_map(long date_thr,const char *prvt); //In path_*prvt is the path where these 
         //files will be stored. Here date_thr has same meaning as in gopen_read function
      bool gopen_map(const char *prvt); //Takes place of gopen_ read and sets up for reading with 
         //standard functions. *prvt has some use here as in create_map
      void gclose_map(const char *prvt); //Takes place of gclose_read 
         // *prvt has some use here as in create_map

      int nsl; //slice number
      ofstream *pout_slice; //Pointer at output file streams. Opened in append mode
      ifstream *pin_slice;  //Pointer at input file streams for reading 
      ofstream *pout_sts;  //Pointer at output file for marking status
      long *ids; //Pointer at array of ids
      long *dat; //Pointer at array of dates for records
      long *adr; //Pointer at array of record addresses
      char *sts; //Pointer at array of status chars
      long *ord; //Pointer at array of indices
      long siz; //Total number of records in slice
      long csiz; //Current number of records indexed
      long cind; //Current record index 
};

class DStor : public FBase {
   public:
      DStor(const char *nam,const char *path_nam); //path_nam is pointer at a string
         //sss and reads the path from file path_sss in current directory. But if sss
         //begins with ':' then skips this character and remaining string is the path
         //string itself.
      DStor(const char *nam,const char *path_nam,long bg,long ed); //Like previous
         //constructor except that bg and ed are beginning slice and end slice+1, resp.
     ~DStor();
      //Read in basic id ranges and store
      bool readlims(void); //This reads the file of cut points and returns the
         //number found or zero if the numbers are not in increasing order.
         //The file extension is ".ctp" and the file is humanly edited and not
         //changed by any functions in the class.
         //There will be num+1 cut points. The smallest is a lower bound for ids 
         //and the largest is 1 larger than the upper bound. These values are 
         //stored in the array mm and a binary search is done to find which file 
         //an id belongs to.
      //Add to the store
      bool gopen_add(void); //calls gopen_add in each slice
         //No index files are used and no check is done when adding records.
         //They are simply appended to the appropriate file.
      bool add(DSpan &Ds); //Writes the data to
         //appropriate file. Returns true if success. Fails if id out of range and 
         //returns false Uses binary search on the mm array to find correct file where to append
      void gclose_add(void); //Closes all files.
      //Create indices
      //create index files ".ids", ".dat", ".adr" all longs in binary and ".sts" and ".siz"
      bool createIndices(void); //For all files, calls createIndex, true for success, 
         //otherwise false. Uses slice functions to perform
      void checkStatus(void); //Reads the status file for each slice and reports how
         //many records and how many obsolete for each slice
      //Read from store
      void gopen_read(void); //Prepares for reading. Calls gopen_read() for each slice
      void gopen_read(long date_thr); //Prepares for reading. Calls same for each slice
      bool exists(long id); //Checks if pmid in DStor and returns true if is, else false
      long date(long id); //Checks if pmid in DStor and returns date if is, else returns zero
      bool read(long id,DSpan &Ds);
         //Calls read in appropriate slice after a binary search. Returns true if success
         //otherwise false 
      void reset(void); //Sets the state to 0 slice and calls reset in all slices
      long getNext(void); //Returns the next 'a' record id+1. Or 0 if fails (no more records).
      long date(void); //Returns date for pmid last returned by getNext() 
      bool read(DSpan &Ds);//read after a successful call t  getNext().
         //Returns false if read fails
         //Returns false if head and tail counts not same, true if are the same.
      void gclose_read(void); //Deletes arrays and closes files
      //Clean files by removing all obsolete records
      void gopen_markStatus(void); //Opens the "sts" files for random writes
         //Mark functions require the gopen_read function be called before they can be used
      bool markStatus(long id,char status); //Uses the state with gopen_read() to find 
         //active record if it exists in the store and mark it with status. Calls 
         //function in slice
      bool markStatus(char status); //Uses the state with gopen_read() calls 
         //corresponding slice function marks current record 
      void gclose_markStatus(void); //Closes the "sts" files
      //Clean slice by removing all obsolete records
      bool cleanAll(void);//cleans all files. This must follow a call to some 
         //form of gopen_read. Calls cleanSlice() functions in slices.
      //Special functions for more rapid access to the data
      void create_map(const char *prvt); //In path_*prvt is the path where these 
         //files will be stored
      void create_map(long date_thr,const char *prvt); //In path_*prvt is the path where these 
         //files will be stored. Here date_thr has same meaning as in gopen_read function
      void gopen_map(const char *prvt); //Takes place of gopen_read and sets up for reading with 
         //standard functions. *prvt has some use here as in create_map
      void gclose_map(const char *prvt); //Takes place of gclose_read 
         // *prvt has some use here as in create_map
        
      long num; //Number of  slices
      long begg; //Beginning slice number
      long endd; //One larger than largest slice number to use
      long *mm; //Array of boundary values for id ranges for slices
         //read in by readlims() in constructor
      DSlice **pSl; //Array of pointers at slice objects
      Bnum *pBn;  //Binary search object for random access
      long curf; //Current file number
};

}
#endif
