###########################################################
#
# Example GNUmakefile to build an Application with Marlin
#
# @author Frank Gaede
# @version $Id: GNUmakefile,v 1.1 2006-02-24 11:01:06 gaede Exp $
#     last modification by A.Zhelezov -- include dependencies
#                                     and HBOOK if needed
###########################################################

# requires path to MARLIN and Marlin v00-09 or higher
# can be called from top level makefile


#----------------------------------------------------------------------------
#                user variables 
#----------------------------------------------------------------------------


# name of the program/library
PROGNAME = AbsCalibration


#----  additional user include/libraries
#-- optionally can be specified in the to level GNUmakefile of $MARLIN
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


USERINCLUDES = -I ../../MarlinCommon/include
USERINCLUDES += -I../../../../CEDViewer/include
USERINCLUDES += -I/opt/products/gsl/1.6/include
USERINCLUDES += -I /opt/products/cernlib/2003/include/cfortran

ifndef USERLIBS
 USERLIBS =  -L ../../MarlinCommon/lib -lMarlinCommon
 USERLIBS += -L /opt/products/gsl/1.6/lib -lgsl -lgslcblas
 USERLIBS += -L /opt/products/lib/ -lstdc++ -lg2c
endif


#----------------------------------------------------------------------------
#        do not modify anything below in order to use top level makefile
#----------------------------------------------------------------------------

# the base directory of your application
# this is where you  would normaly have the lib, bin and source directories
BASEDIR = ..
INSTALL = $(BASEDIR)/lib
BINDIR = $(BASEDIR)/bin

ifdef MARLINDEBUG
 CPPFLAGS = -g -Wall -ansi
else
 CPPFLAGS = -O3 -Wall -ansi
endif


# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

AR = ar

#--------------------------------------------------------------------------------------
LIBS +=   -Wl,-whole-archive  -L $(INSTALL)  -l$(PROGNAME) -Wl,-no-whole-archive
#--------------------------------------------------------------------------------------

INCLUDE +=  -I $(BASEDIR)/include
INCLUDE += $(USERINCLUDES)
CPPFLAGS = -g -ansi
CPPFLAGS += $(INCLUDE)
CPPFLAGS += `. $(MARLIN)/bin/marlin_includes.sh`

LIBS += $(USERLIBS)
LIBS += `. $(MARLIN)/bin/marlin_libs.sh`

all:	lib

bin:  lib $(BINDIR)/$(PROGNAME)

$(BINDIR)/$(PROGNAME): $(INSTALL)/lib$(PROGNAME).a
	mkdir -p $(BINDIR)
	$(CXX) -o $(BINDIR)/$(PROGNAME) $(CPPFLAGS) $(MARLIN)/src/Marlin.cc $(LIBS) 

objects = $(patsubst %.cc,%.o, $(wildcard *.cc)  )

$(INSTALL)/lib$(PROGNAME).a: $(objects)
	mkdir -p $(INSTALL)
	rm -f $@
	$(AR) qc $@ $+
	ranlib $@

lib: $(INSTALL)/lib$(PROGNAME).a


clean:
	rm -f *.o *.P $(INSTALL)/lib$(PROGNAME).a $(BINDIR)/$(PROGNAME)

doc:
	cd $(BASEDIR)/doc ; doxygen ;

MAKEDEPEND=$(CXX) -M $(CPPFLAGS) -o $*.d $<


%.P : %.cc
	$(MAKEDEPEND)
	@sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' < $*.d > $@; \
	rm -f $*.d; [ -s $@ ] || rm -f $@

include $(objects:.o=.P)

