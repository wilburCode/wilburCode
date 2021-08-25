SRC_DIR=/panfs/pan1.be-md.ncbi.nlm.nih.gov/cbb-wilbur/wilbur/CentosCPP64/lib
TRASHFILES = *.o *~ *.bak core	
LIB_INC=-I/panfs/pan1.be-md.ncbi.nlm.nih.gov/cbb-wilbur/wilbur/CentosCPP64/include
#.KEEP_STATE:
libops.a: runn.o  Btree.o Word.o   \
	FBase.o Patt.o Hyper.o GBoost.o	QStor.o\
	Istreg.o Thes.o Split.o\
	Isgrid.o Istgen.o Bnum.o LUD.o Elev.o\
	Sock.o	\
	Bootstrap.o Hash.o DataObj.o Brill.o Filter.o \
	Repeat.o Vnab.o PostBrill_FN.o PostBrill_FP.o \
	LinClass.o Strset.o Mahal.o Dist.o DStor.o \
	Lglin.o SBoost.o SCmls.o PStor.o SStor.o NStor.o\
	Island.o huge.o mHMM.o Sav.o \
	Nnls.o Repeats.o Clip.o \
	Perm.o lbfgs.o XPost.o Dset.o	\
	Map.o EdtPrb.o MPlex.o MPtok.o MPtag.o \
	GoldOpt.o \
	MPpar.o LexMkv.o StrPac.o\
	MPcandc.o \
	Ab3P.o AbbrStra.o AbbrvE.o \

#	rm -f libops.a
	ar rus $@ $?
#	$(CCC) -xar -o $@ $?
#rm -f *.o
#

%.o: $(SRC_DIR)/%.C
	g++ -std=c++11 -c -gdwarf-2 $< -o $@ $(LIB_INC)

#.C~.o:
#	$(GET) $(GFLAG) -p $< > $*.C

clean: rm -f $(TRASHFILES)

