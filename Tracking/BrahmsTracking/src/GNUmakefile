###########################################################
#
# Example GNUmakefile to build an Application with Marlin
#
# @author Frank Gaede
# @version $Id: GNUmakefile,v 1.3 2005-11-03 15:16:14 aplin Exp $
#
##########################################################

# requires path to MARLIN and Marlin v00-09 or higher
# can be called from top level makefile


#----------------------------------------------------------------------------
#                user variables 
#----------------------------------------------------------------------------


# name of the program/library
PROGNAME = BrahmsTracking


#----  additional user include/libraries
#-- optionally can be specified in the to level GNUmakefile of $MARLIN

ifndef FORTRANINCLUDES
FORTRANINCLUDES = -I$(MARLIN)/packages/MarlinUtil/include/
endif

ifndef USERINCLUDES
 USERINCLUDES  =  -I $(MARLIN)/packages/MarlinUtil/include
 USERINCLUDES +=  -Df2cFortran
 #USERINCLUDES = `root-config --cflags`
endif

ifndef USERLIBS
 USERLIBS = -L -I $(MARLIN)/packages/MarlinUtil/lib -lMarlinUtil
 USERLIBS += -lg2c
 USERLIBS +=  -L /opt/products/cernlib/pro/lib/ -lmathlib -lkernlib 
 #USERLIBS = `root-config --libs`
endif


#----------  BrahmsTracking needs to compile some f77 files ----------
#brahms object files. Note: "-Wno-globals" used to stop endless warnings when using cernlib funcitons remove this if you feel uneasy about ignoring these warnings
./f77/%.o: ./f77/%.F
	$(FC) -g -Wno-globals -o $@ -ff2c -c $? $(FORTRANINCLUDES)

FOBJS = $(patsubst %.F,%.o, $(wildcard f77/*.F) )
# FOBJS also added to objects below

./C/%.o: ./C/%.C
	$(CC) -g -ansi -pedantic -o $@ -c $?
COBJS = $(patsubst %.c,%.o, $(wildcard C/*.c) )

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
CPPFLAGS += $(INCLUDE)
CPPFLAGS += `. $(MARLIN)/bin/marlin_includes.sh`

LIBS += $(USERLIBS)
LIBS += `. $(MARLIN)/bin/marlin_libs.sh`

bin:  lib $(BINDIR)/$(PROGNAME)

$(BINDIR)/$(PROGNAME): $(INSTALL)/lib$(PROGNAME).a
	mkdir -p $(BINDIR)
	$(CXX) -o $(BINDIR)/$(PROGNAME) $(CPPFLAGS) $(MARLIN)/src/Marlin.cc $(LIBS) 

objects = $(patsubst %.cc,%.o, $(wildcard *.cc)  )
objects += $(FOBJS)
objects += $(COBJS)




$(INSTALL)/lib$(PROGNAME).a: $(objects)
	mkdir -p $(INSTALL)
	$(AR) cr $(INSTALL)/lib$(PROGNAME).a  $?

lib: $(INSTALL)/lib$(PROGNAME).a


clean:
	rm -f *.o $(INSTALL)/lib$(PROGNAME).a $(BINDIR)/$(PROGNAME) $(FOBJS) $(COBJS)

doc:
	cd $(BASEDIR)/doc ; doxygen ;

