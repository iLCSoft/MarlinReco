###########################################################
#
# Example GNUmakefile to build an Application with Marlin
#
# @author Frank Gaede
# @author Jan Engels
# @version $Id: GNUmakefile,v 1.4 2007-02-16 16:52:43 engels Exp $
#
###########################################################
                                                                                                                                                             

#----------------------------------------------------------------------------
#        user variables
#----------------------------------------------------------------------------

# name of the program
PROGNAME = TPCDigi

# source file for the binary
# leave this commented out if your Application does not produce a binary file
#BINSRC = $(PROGNAME).cc

# marlin working directory
ifndef MARLINWORKDIR
 MARLINWORKDIR=$(MARLIN)
 export MARLINWORKDIR
endif

# additional user include/libraries
ifndef USERINCLUDES
 USERINCLUDES += -I$(MARLINWORKDIR)/packages/MarlinUtil/include
 USERINCLUDES += -I$(GSL_HOME)/include
endif

ifndef USERLIBS
 USERLIBS +=
endif

#----------------------------------------------------------------------------
#        do not modify anything below in order to use top level makefile
#----------------------------------------------------------------------------

# directories
INCDIR = ../include
BINDIR = $(MARLINWORKDIR)/bin
DOCDIR = $(MARLINWORKDIR)/doc/MarlinReco/$(PROGNAME)
TMPDIR = $(MARLINWORKDIR)/tmp/MarlinReco
LIBDIR = $(TMPDIR)/lib
OBJDIR = $(TMPDIR)/obj/$(PROGNAME)
DEPDIR = $(TMPDIR)/dep/$(PROGNAME)

ifdef BINSRC
# binary file with path
BIN = $(BINDIR)/$(patsubst %.cc,%, $(BINSRC))

# object file for the binary without path
BINOBJS = $(patsubst %.cc,%.o,$(BINSRC))

# object file for the binary with path
BINOBJ = $(OBJDIR)/$(patsubst %.cc,%.o,$(BINSRC))
endif

# library file without path
LIBNAME = lib$(PROGNAME).a

# library file with path
LIB = $(LIBDIR)/$(LIBNAME)

# source files wothout path
SRCS = $(filter-out $(BINSRC), $(wildcard *.cc))

# object files without path
OBJS = $(patsubst %.cc,%.o, $(SRCS))

# object files with path
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJS))

# dependencie files without path
DEPS = $(patsubst %.o,%.d,$(OBJS) $(BINOBJS))

# dependencie files with path
DEP = $(patsubst %,$(DEPDIR)/%,$(DEPS))

# c++ compiler
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

# compiler flags
ifdef MARLINDEBUG
 CXXFLAGS = -g -Wall -ansi
else
 CXXFLAGS = -O3 -Wall -ansi
endif

# pre-processor flags
CPPFLAGS = -I $(INCDIR)
CPPFLAGS += $(USERINCLUDES)
CPPFLAGS += `. $(MARLINWORKDIR)/bin/marlin_includes.sh`

# libraries
LIBS = -Wl,-whole-archive -L $(LIBDIR) -l$(PROGNAME) -Wl,-no-whole-archive
LIBS += $(USERLIBS)
LIBS += `. $(MARLINWORKDIR)/bin/marlin_libs.sh`

#compile command
COMPILE = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c

# linker
LD = $(CXX)

# linker flags
LDFLAGS = -Wl,--demangle

# archiever
AR = ar cru

# check if a directory exists
DIR_EXISTS = test -d

.PHONY: all lib doc clean distclean

all: lib

ifdef BINSRC
.PHONY: bin
bin: lib $(BINOBJ) $(BIN)
	rm -f $(BINOBJ)
endif

lib: $(LIB)

doc:
	@$(DIR_EXISTS) $(DOCDIR) || mkdir -p $(DOCDIR) ; \
	export DOX_PROGNAME=$(PROGNAME) ; \
	export DOX_DOCDIR=$(DOCDIR) ; \
	cd ../doc ; doxygen ;

clean:
	rm -f $(OBJ) $(BINOBJ) $(BIN) $(LIB)

distclean:
	rm -rf $(BIN) $(LIB) $(TMPDIR)/*/$(PROGNAME) $(DOCDIR)

# tell make to include dependencies
-include $(DEP)

# binary rule
$(BIN) : $(OBJ) $(BINOBJ)
	@$(DIR_EXISTS) $(BINDIR) || mkdir -p $(BINDIR)
	$(LD) $(BINOBJ) $(LDFLAGS) $(LIBS) -o $@

# objects rule
$(OBJDIR)/%.o: %.cc
	@$(DIR_EXISTS) $(OBJDIR) || mkdir -p $(OBJDIR)
	$(COMPILE) $< -o $@

# library rule
$(LIB) : $(OBJ)
	@$(DIR_EXISTS) $(LIBDIR) || mkdir -p $(LIBDIR)
	$(AR) $@ $?

# dependencies rule
$(DEPDIR)/%.d : %.cc
	@$(DIR_EXISTS) $(DEPDIR) || mkdir -p $(DEPDIR)
	@echo "Creating dependencies for $<" ; \
	set -e; $(CXX) -MM $(CPPFLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/$(subst /,\/,$(OBJDIR))\/\1.o $(subst /,\/,$(DEPDIR))\/$*.d : /g'\
	> $@; [ -s $@ ] || rm -f $@

# end
