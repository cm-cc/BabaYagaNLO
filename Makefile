EXP= DEFAULT

DEXP=
ifneq ($(EXP),DEFAULT)
DEXP=-D$(EXP)
endif
ifeq ($(EXP),DEFAULT)
EXP=
endif

RECOLA=
RECOLA=-DRECOLA

COLLIER=
COLLIER=-DCOLLIER

LT =
LT =-DLOOPTOOLS

QUAD=
QUAD=-quad

LIBFILES=
LTVER=2.16
LTDIR=LoopTools-$(LTVER)
ifeq ($(LT),-DLOOPTOOLS)
  LIBFILES += $(LTDIR)/lib64/libooptools$(QUAD).a
  LTINC=-I$(LTDIR)/include
  LTSTRING =-looptools
  ifdef QUAD
    LTSTRING =-looptools$(QUAD)
    LTQUAD   = --quad
    QUADTYPE =-DRealType=real*16 -DComplexType=complex*32
  endif
endif

FRED=alphaQEDc19/libfred.a
FRED=

CLLDIR = collier/COLLIER-1.2.9
ifeq ($(COLLIER),-DCOLLIER)
  LIBFILES += $(CLLDIR)/lib/libcollier.a
  CLLMOD = -I$(CLLDIR)/include/
  CLLLIB = -L$(CLLDIR) -lcollier
endif

ifeq ($(RECOLA),-DRECOLA)
  RCLVERSION = 1.5.0
  RCLDIR = recolas/recola-$(RCLVERSION)
  RCLMOD = -I$(RCLDIR)/include/
  RCLLIB = -L$(RCLDIR)/lib/ -lrecola
  LIBFILES += $(RCLDIR)/lib/librecola.a
endif

EXE = babayaga$(EXP)

PAWyn =
PAW =
CERNLIB =


F77 = gfortran
CC  = gcc
# tune for your processor
FFLAGS = 

FFLAGS = -O3 -ffast-math -funroll-loops -fPIC
ifeq ($(OPT),YES)
 FFLAGS = -O3 -march=native -mtune=native  -ffast-math -funroll-loops -fPIC
endif


CFLAGS=$(FFLAGS)

# FFLAGS = -Wall -Wextra -fcheck=all -finit-real=snan -finit-integer=-999 -finit-logical=true\
# 	-fbacktrace -g -O0

## C ranlux optimizations, AVX2 faster but only on recent CPUS, SSE2 slower but present on
## any neowulf node. If both are specified, -DSSE2 is ignored (and won't run on old neowulf nodes).
RLXOPT=-DSSE2
ifeq ($(OPT),YES)
 RLXOPT=-DSSE2 -DAVX2
endif

FPCHECK = trapfpe.c
FPCHECK =

EXTRADEPS = #Makefile

VPHLMNT=vp_hlmnt_v2_1
VPHLMNT=vp_hlmnt_v2_1_1
VPHLMNT=vp_hlmnt_v2_2

OBJECTS = gen_events.o initcloseby.o cuts.o sv.o matrix_model.o mapmomenta.o loops.o ffpi.o\
          routines.o sampling.o phasespace.o distributions.o $(VPHLMNT).o\
          hadr5n16.o hadr5n09.o hadr5x23.o userinterface.o intpl.o Rteubner.o\
          hadr5n17.o hadr5n12.o c_rnlx_interface.o ranlux_common.o ranlxd.o ranlxs.o recola_int.o\
          hard_ampl.o pent.o storage.o muemuegg.o ALPHA.o

F77 += $(FFLAGS)

default: $(EXE)

SAVEDIR = release
RELEASEDIR = BabaYaga
pack: # use only to release BABAYAGA
	mkdir -p $(RELEASEDIR)/form &&\
	mkdir -p $(RELEASEDIR)/ffpi-models/ &&\
	mkdir -p $(RELEASEDIR)/c_ranlux &&\
	mkdir -p $(RELEASEDIR)/collier/ &&\
	mkdir -p $(RELEASEDIR)/$(RCLDIR)/ &&\
	cp -ra oneloop/ $(RELEASEDIR)/oneloop &&\
	cp -ra LoopTools-$(LTVER)-clean/ $(RELEASEDIR)/LoopTools-$(LTVER) &&\
	cp -ra LoopTools-$(LTVER)-clean/ $(RELEASEDIR)/LoopTools-$(LTVER)-clean &&\
	cp -ra $(CLLDIR)-clean/ $(RELEASEDIR)/$(CLLDIR) &&\
	cp -ra $(RCLDIR)-clean/* $(RELEASEDIR)/$(RCLDIR)/ &&\
	cp -ra $(CLLDIR)-clean/ $(RELEASEDIR)/$(CLLDIR)-clean/ &&\
	cp -ra $(RCLDIR)-clean/ $(RELEASEDIR)/$(RCLDIR)-clean/ &&\
	cp -ra pipirad/ mumurad/ eerad/ pions/ $(RELEASEDIR)/ &&\
        cp -r input_rad Makefile README shared.F\
        main.F loops.F gen_events.F commonmain.F cuts.F sv.F matrix_model.F mapmomenta.F vpol_novosibirsk.dat\
        routines.F invariants.h sampling.F phasespace.F distributions.F $(VPHLMNT).F\
        hadr5n16.F hadr5n09.F userinterface.F intpl.F Rteubner.F strong2020common.F driver_gen_events.F\
        hadr5n17.F hadr5n12.F hadr5x23.F hard_ampl.F interface.F mtx_eeenudbb.F ffpi.F initcloseby.F\
        recola_int.F storage.F dalhadshigh17.F dalhadslow17.F dalhadt17.F dalhadthigh17.F prm.f\
        vpol_all_bare_sum_v2.9.dat vpol_bare_lept_v2.9.dat muemuegg.F muemuegg-plus.f pent.F ALPHA.F\
        $(RELEASEDIR) &&\
	cp form/*.[fF] $(RELEASEDIR)/form/ &&\
	cp c_ranlux/* $(RELEASEDIR)/c_ranlux/ &&\
	mkdir $(SAVEDIR) ;\
	tar cjvf $(SAVEDIR)/babayaga.tar.bz2 $(RELEASEDIR)/ &&\
	rm -rf $(RELEASEDIR)
clean:
	rm -f $(OBJECTS) *.o *.a
deepclean:
	rm -rf $(OBJECTS) *.o *.a *.so $(EXE) $(EXE)full *~ form/*~ $(LTDIR)/lib64/ $(LTDIR)/bin/ $(LTDIR)/build/ $(LTDIR)/makefile $(LTDIR)/include/ $(LTDIR)/build$(QUAD)/ $(CLLDIR)/build/ $(CLLDIR)/modules/* $(CLLDIR)/include/ $(CLLDIR)/collierC*.cmake $(CLLDIR)/lib/ $(CLLDIR)/lib/libcollier.a $(RCLDIR)/lib/ $(RCLDIR)/build/ $(CLLDIR)/libcollier.a $(RCLDIR)/librecola.a

# C version of ranlux by Martin Luscher
# http://luscher.web.cern.ch/luscher/ranlux/index.html
c_rnlx_interface.o: c_ranlux/c_rnlx_interface.c $(EXTRADEPS) Makefile
	$(CC) $(CFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/c_rnlx_interface.c
ranlxd.o: c_ranlux/ranlxd.c $(EXTRADEPS) Makefile 
	$(CC) $(CFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlxd.c
ranlxs.o: c_ranlux/ranlxs.c $(EXTRADEPS) Makefile
	$(CC) $(CFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlxs.c
ranlux_common.o: c_ranlux/ranlux_common.c $(EXTRADEPS) Makefile
	$(CC) $(CFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlux_common.c
####

# source files
main.o: main.F commonmain.F $(EXTRADEPS) Makefile
	$(F77) -c main.F
gen_events.o: gen_events.F commonmain.F $(EXTRADEPS) collier
	$(F77) -c $(COLLIER) $(CLLMOD) gen_events.F
initcloseby.o: initcloseby.F commonmain.F $(EXTRADEPS) collier
	$(F77) -c $(COLLIER) $(CLLMOD) initcloseby.F
cuts.o: cuts.F $(EXTRADEPS) strong2020common.F
	$(F77) $(DEXP) -c cuts.F
matrix_model.o: matrix_model.F  form/formme.F form/ee3gexact.F form/borngg.F form/formmemm.F pipirad/hard/*.F $(EXTRADEPS)
	$(F77) -c matrix_model.F
sv.o: sv.F $(EXTRADEPS) Makefile  collier
	$(F77) -c $(COLLIER) $(CLLMOD) $(QUADTYPE) sv.F
2loop.o: 2loop.F $(EXTRADEPS)
	$(F77) -c 2loop.F
hplog.o: 2loop/hplog.F $(EXTRADEPS)
	$(F77) -c 2loop/hplog.F
userinterface.o: userinterface.F $(EXTRADEPS) collier
	$(F77) $(DEXP) $(RECOLA) -c userinterface.F
phasespace.o: phasespace.F $(EXTRADEPS) strong2020common.F
	$(F77) -c $(DEXP) phasespace.F
storage.o: storage.F $(EXTRADEPS)
	$(F77) -c storage.F
hadr5n16.o: hadr5n16.F $(EXTRADEPS)
	$(F77) -c hadr5n16.F
hadr5n17.o: hadr5n17.F $(EXTRADEPS)
	$(F77) -c hadr5n17.F
ffpi.o: ffpi.F $(EXTRADEPS)
	$(F77) -c ffpi.F
muemuegg.o: muemuegg.F $(EXTRADEPS)
	$(F77) -c muemuegg.F
hadr5x23.o: hadr5x23.F $(EXTRADEPS)
	$(F77) -c hadr5x23.F
hadr5n12.o: hadr5n12.F $(EXTRADEPS)
	$(F77) -c hadr5n12.F
hadr5n09.o: hadr5n09.F $(EXTRADEPS)
	$(F77) -c hadr5n09.F
intpl.o: intpl.F $(EXTRADEPS)
	$(F77) -c intpl.F
$(VPHLMNT).o: $(VPHLMNT).F $(EXTRADEPS)
	$(F77) -c $(VPHLMNT).F
Rteubner.o: Rteubner.F $(EXTRADEPS)
	$(F77) -c Rteubner.F
mapmomenta.o: mapmomenta.F $(EXTRADEPS)
	$(F77) -c mapmomenta.F
sampling.o: sampling.F $(EXTRADEPS) strong2020common.F
	$(F77) -c $(DEXP) sampling.F
loops.o: loops.F $(EXTRADEPS) collier
	$(F77) -c $(COLLIER) $(CLLMOD) $(QUADTYPE) loops.F
routines.o: routines.F $(EXTRADEPS) collier
	$(F77) $(COLLIER) $(CLLMOD) -c routines.F
distributions.o: distributions.F shared.F $(EXTRADEPS) Makefile strong2020common.F
	$(F77) $(DEXP) -c distributions.F
#ranlux.o: ranlux.F $(EXTRADEPS)
#	$(F77) -c ranlux.F
recola_int.o: recola_int.F $(EXTRADEPS) Makefile $(RCLDIR)/lib/librecola.a
	$(F77) $(RECOLA) $(RCLMOD) -c recola_int.F	
hard_ampl.o: hard_ampl.F $(EXTRADEPS)
	$(F77) -c hard_ampl.F
pent.o: pent.F $(EXTRADEPS) $(LTDIR)/lib64/libooptools$(QUAD).a
	$(F77) -c $(LT) $(LTINC) $(QUADTYPE) pent.F

# lightweight ALPHA version, by Mauro Moretti
ALPHA.o: ALPHA.F prm.f
	$(F77) -ffixed-line-length-132 -c ALPHA.F
############

looptools: $(LTDIR)/lib64/libooptools$(QUAD).a
$(LTDIR)/lib64/libooptools$(QUAD).a:
	@echo "Building LoopTools"
	@echo " "
	cd $(LTDIR) && ./configure $(LTQUAD) --prefix=. && make -j`nproc` && make install

collier: $(CLLDIR)/lib/libcollier.a
$(CLLDIR)/lib/libcollier.a:
	@echo "Building Collier"
	@echo " "
	mkdir -p $(CLLDIR)/build/
	cd $(CLLDIR)/build/ && cmake -Dstatic=On .. -DCMAKE_INSTALL_PREFIX=..  && make && make install

recola: $(RCLDIR)/lib/librecola.a collier
$(RCLDIR)/lib/librecola.a: $(CLLDIR)/lib/libcollier.a
	@echo "Building RECOLA"
	@echo " "
	mkdir -p $(RCLDIR)/build/
	cd $(RCLDIR)/build/ && cmake .. -Dcollier_path=../../../$(CLLDIR) -DCMAKE_INSTALL_PREFIX=.. -Dstatic=yes && make && make install

# babayaga library
libbabayaga.a: $(OBJECTS)
	ar cr libbabayaga.a $(OBJECTS)

# babayaga library + all external libraries
libbabayagafull.a: $(OBJECTS) $(LIBFILES)
	@rm -f libbabayagafull.a libbabayagafull.so
#	@echo "Creating BabaYaga static and shared libraries, including all external libraries"
#	@echo "Creating BabaYaga library, including all external libraries"
#	@echo " "
	@mkdir -p .arlibs/
	@cp -a $(LIBFILES) .arlibs/
	@cp -a $(OBJECTS) .arlibs/
	@cd .arlibs && ar x libooptools$(QUAD).a && ar x libcollier.a && ar x librecola.a
	@ar cr libbabayagafull.a .arlibs/*.o
	@ld -shared -o libbabayagafull.so .arlibs/*.o
	@rm -rf .arlibs

# executable
$(EXE): $(LIBFILES) libbabayaga.a main.o libbabayagafull.a
	$(F77) main.o -L. -lbabayaga $(FPCHECK) -o $(EXE) $(FRED) $(CLLLIB) $(RCLLIB) -L$(LTDIR)/lib64 $(LTSTRING)

$(EXE)full: $(LIBFILES) libbabayagafull.a main.o
	$(F77) main.o -L. -lbabayagafull -o $(EXE)full -Wl,-rpath=.
