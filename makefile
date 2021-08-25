LIBPATH=/panfs/pan1.be-md.ncbi.nlm.nih.gov/cbb-wilbur/wilbur/CentosCPP64/lib
LIBCOM=/home/comeau/lib
INCCOM=/home/comeau/include
ALT_LIBPATH=/panfs/pan1.be-md.ncbi.nlm.nih.gov/cbb-wilbur/wilbur/CentosCPP64/altlib
INCPATH=/panfs/pan1.be-md.ncbi.nlm.nih.gov/cbb-wilbur/wilbur/CentosCPP64/include
TRASHFILES = $(LIBDIR)/*.o *~ *.bak core	

.KEEP_STATE:
$(N): $(N).o /panfs/pan1.be-md.ncbi.nlm.nih.gov/cbb-wilbur/wilbur/CentosCPP64/lib/libops.a
	g++ $(OS) -o $(N) $(N).o -L$(LIBPATH) -lops -lpcre -L$(LIBCOM) -lDCCiret -D_POSIX_C_SOURCE=199506L -lpthread

.SUFFIXES:
%.o: %.C 
	g++ -std=c++11 -c $(OS) -I$(INCPATH) -I$(INCCOM) $< -o $@

clean: rm -f $(TRASHFILES)

