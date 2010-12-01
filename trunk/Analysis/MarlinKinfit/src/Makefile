#========================================================
#	a simple makefile to compile our first programs
#========================================================
.SUFFIXES : .o .c .C .h

CC = g++ -g -Wall -pg
CFLAGS= -I`root-config --incdir` -I ../..
CLIBS= -L/cern/pro/lib -lmathlib -lkernlib -lpacklib -L /usr/lib -lshift -lnsl -lcrypt -ldl
LDFLAGS= $(CFLAGS)  `root-config --libs` -lgsl
FLIBS= -lg2c -lm 
F77=g77
FC=g77
# CLIBS=`cernlib`
GSLINCDIR=/opt/products/gsl/1.9/include
INCDIR=../include
LIBDIR=./

all:  .rootversion lib

clean:
	rm -f Makefile.dep *.o *.a

.rootversion: rootversion

rootversion:
	@if [[ -f .rootversion ]]; \
        then \
          echo ".rootversion exists: " ; \
          cat .rootversion ; \
          if [[ `cat .rootversion` == `root-config --version` ]] ; \
          then echo "root version has not changed." ; \
          else echo "root version has changed! => will recompile everything."; \
            root-config --version > .rootversion ; \
            make clean;  \
          fi \
        else echo ".rootversion will be created." \
          root-config --version > .rootversion ; \
          make clean ;\
        fi

Makefile.dep: ${INCDIR}/*.h *.cc .rootversion
	touch Makefile.dep
	makedepend -fMakefile.dep -I${INCDIR} -I`root-config --incdir` *.cc

tags:
	ctags ./*.[hC]
        
depend: Makefile.dep
        
include Makefile.dep


LIB_OBJS= BaseFitObject.o \
             ParticleFitObject.o JetFitObject.o NeutrinoFitObject.o \
             PhotonFitObject.o PhotonFitObjectPxyg.o \
             BaseConstraint.o \
             ParticleConstraint.o  \
             PConstraint.o MomentumConstraint.o  \
             MassConstraint.o \
             BaseFitter.o \
             NewtonFitter.o OPALFitter.o  \
             OPALFitterGSL.o NewtonFitterGSL.o  \
             BaseTracer.o \
             TextTracer.o RootTracer.o \
             TopEventILC.o TopEventILCSoft.o FourVector.o \
             TwoB4JPairing.o  FourJetPairing.o \
             SoftGaussParticleConstraint.o SoftGaussMassConstraint.o \
             SoftGaussMomentumConstraint.o \
             SoftBWParticleConstraint.o SoftBWMassConstraint.o \
             dgesvd.o
             
lib:libkinfit.a
        
libkinfit.a: $(LIB_OBJS)
	ar rv libkinfit.a $(LIB_OBJS)


.cc.o:
	$(CC) $< $(CFLAGS) -I ${INCDIR} -I ${GSLINCDIR} -c -o $@

