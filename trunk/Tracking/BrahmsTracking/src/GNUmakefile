###########################################################
#
# Example GNUmakefile to build an Application with Marlin
#
# @author Frank Gaede
# @author Jan Engels
# @version $Id: GNUmakefile,v 1.9 2007-03-01 15:54:34 engels Exp $
#
###########################################################


#----------------------------------------------------------------------------
#        user variables
#----------------------------------------------------------------------------

# name of the program
PROGNAME = BrahmsTracking

# source files for the binary
# leave this commented out if your Application does not produce binary files
#BINSRCS = $(PROGNAME).cc test.cc

# marlin working directory
ifndef MARLINWORKDIR
 MARLINWORKDIR=$(MARLIN)
 export MARLINWORKDIR
endif

# additional user include/libraries
ifndef USERINCLUDES
 USERINCLUDES += -I$(MARLINWORKDIR)/packages/MarlinUtil/include
 USERINCLUDES += -Df2cFortran
endif

ifndef USERLIBS
 USERLIBS +=
endif

F77INCLUDES = -I$(MARLINWORKDIR)/packages/MarlinUtil/include/

#----------------------------------------------------------------------------
#        do not modify anything below in order to use top level makefile
#----------------------------------------------------------------------------

# directories
INCDIR = ../include
LIBDIR = ../lib
BINDIR = $(MARLINWORKDIR)/bin
DOCDIR = $(MARLINWORKDIR)/doc/MarlinReco/$(PROGNAME)
TMPDIR = $(MARLINWORKDIR)/tmp/MarlinReco
OBJDIR = $(TMPDIR)/obj/$(PROGNAME)
DEPDIR = $(TMPDIR)/dep/$(PROGNAME)
F77DIR = ./f77
CDIR = ./C

# library file without path
LIBNAME = lib$(PROGNAME).a

# library file with path
LIB = $(LIBDIR)/$(LIBNAME)

ifdef BINSRCS
# binary files with path
BIN = $(patsubst %.cc,$(BINDIR)/%,$(BINSRCS))

# object files for the binaries without path
BINOBJS = $(patsubst %.cc,%.o,$(BINSRCS))

# object files for the binaries with path
BINOBJ = $(patsubst %,$(OBJDIR)/%,$(BINOBJS))
endif

# source files without path
SRCS = $(filter-out $(BINSRCS),$(wildcard *.cc))
F77SRCS = $(patsubst $(F77DIR)/%,%,$(wildcard $(F77DIR)/*.F))
CSRCS = $(patsubst $(CDIR)/%,%,$(wildcard $(CDIR)/*.c))

# object files without path
OBJS = $(patsubst %.cc,%.o,$(SRCS))
F77OBJS = $(patsubst %.F,%.o,$(F77SRCS))
COBJS = $(patsubst %.c,%.o,$(CSRCS))

# object files with path
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJS))
F77OBJ = $(patsubst %,$(OBJDIR)/%,$(F77OBJS))
COBJ = $(patsubst %,$(OBJDIR)/%,$(COBJS))

# dependency files without path
DEPS = $(patsubst %.o,%.d,$(OBJS) $(BINOBJS) $(F77OBJS) $(COBJS))

# dependency files with path
DEP = $(patsubst %,$(DEPDIR)/%,$(DEPS))

# c++ compiler
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

# fortran compiler
ifdef MY_FC
  FC = $(MY_FC)
else
  FC = g77
endif

# c compiler
ifdef MY_CC
  CC = $(MY_CC)
else
  CC = gcc
endif

# compiler flags
ifdef MARLINDEBUG
 CXXFLAGS = -g -Wall -ansi
 FCFLAGS = -g -Wno-globals -ff2c
 CCFLAGS = -g -ansi -pedantic
else
 CXXFLAGS = -O3 -Wall -ansi
 FCFLAGS = -O3 -Wall -ansi -ff2c
 CCFLAGS = -O3 -Wall -ansi -pedantic
endif

# pre-processor flags
CPPFLAGS = -I$(INCDIR)
CPPFLAGS += $(USERINCLUDES)
CPPFLAGS += `. $(MARLINWORKDIR)/bin/marlin_includes.sh`

# libraries
#LIBS = -Wl,-whole-archive -L $(LIBDIR) -l$(PROGNAME) -Wl,-no-whole-archive
LIBS += $(USERLIBS)
LIBS += `. $(MARLINWORKDIR)/bin/marlin_libs.sh`

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

ifdef BINSRCS
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
	rm -f $(OBJ) $(F77OBJ) $(COBJ) $(BINOBJ) $(BIN) $(LIB)

distclean:
	rm -rf $(BIN) $(LIB) $(TMPDIR)/*/$(PROGNAME) $(DOCDIR)

# tell make to include dependencies
-include $(DEP)

# binary rule
$(BIN) : $(OBJ) $(F77OBJ) $(COBJ) $(BINOBJ)
	@$(DIR_EXISTS) $(BINDIR) || mkdir -p $(BINDIR)
	$(LD) $(OBJDIR)/$(notdir $@).o $(LDFLAGS) $(LIBS) -o $@

# objects rule
$(OBJDIR)/%.o: %.cc
	@$(DIR_EXISTS) $(OBJDIR) || mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# objects rule for fortran
$(OBJDIR)/%.o: $(F77DIR)/%.F
	$(FC) $(FCFLAGS) $(F77INCLUDES) -c $< -o $@

# objects rule for c
$(OBJDIR)/%.o: $(CDIR)/%.c
	$(CC) $(CCFLAGS) -c $< -o $@

# library rule
$(LIB) : $(OBJ) $(F77OBJ) $(COBJ)
	@$(DIR_EXISTS) $(LIBDIR) || mkdir -p $(LIBDIR)
	$(AR) $@ $?

# dependency rules for c++
$(DEPDIR)/%.d : %.cc
	@$(DIR_EXISTS) $(DEPDIR) || mkdir -p $(DEPDIR)
	@echo "Creating dependencies for $<" ; \
	set -e; $(CXX) -MM $(CPPFLAGS) $< \
	| sed 's|\($*\)\.o[ :]*|$(OBJDIR)/\1.o $(DEPDIR)/$*.d : |g'\
	> $@; [ -s $@ ] || rm -f $@

# dependency rules for fortran
$(DEPDIR)/%.d : $(F77DIR)/%.F
	@$(DIR_EXISTS) $(DEPDIR) || mkdir -p $(DEPDIR)
	@echo "Creating dependencies for $<" ; \
	set -e; $(FC) -MM $(FCFLAGS) $< \
	| sed 's|\($*\)\.o[ :]*|$(OBJDIR)/\1.o $(DEPDIR)/$*.d : |g'\
	> $@; [ -s $@ ] || rm -f $@

# dependency rules for c
$(DEPDIR)/%.d : $(CDIR)/%.c
	@$(DIR_EXISTS) $(DEPDIR) || mkdir -p $(DEPDIR)
	@echo "Creating dependencies for $<" ; \
	set -e; $(CC) -MM $(CCFLAGS) $< \
	| sed 's|\($*\)\.o[ :]*|$(OBJDIR)/\1.o $(DEPDIR)/$*.d : |g'\
	> $@; [ -s $@ ] || rm -f $@

# end
