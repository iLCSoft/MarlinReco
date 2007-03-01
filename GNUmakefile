######################################################
#
# Toplevel Makefile to build MarlinReco
#
# @author Frank Gaede
# @author Jan Engels
#
######################################################

ifndef MARLINWORKDIR
 MARLINWORKDIR=$(MARLIN)
 export MARLINWORKDIR
endif

PROGNAME = MarlinReco
LIBNAME = lib$(PROGNAME).a

LIBDIR = ./lib
TMPDIR = $(MARLINWORKDIR)/tmp/$(PROGNAME)
DOCDIR = $(MARLINWORKDIR)/doc/$(PROGNAME)

LIB = $(LIBDIR)/$(LIBNAME)

subdirs = $(patsubst %/src/GNUmakefile,%, $(wildcard ./*/*/src/GNUmakefile))
objs = $(wildcard $(TMPDIR)/obj/*/*.o)

.PHONY: all lib doc clean distclean packages

all: lib

lib: packages $(LIB)

packages:
	@for dir in $(subdirs); do \
	echo "*************************************" ; \
	echo "*   Building $(PROGNAME) Library $$dir ..." ; \
	echo "*************************************" ; \
	$(MAKE) -C $$dir/src lib ; done;

$(LIB): $(objs)
	@echo "*************************************" ; \
	echo "*    Linking $(PROGNAME) library ..." ; \
	echo "*************************************" ; \
	test -d $(LIBDIR) || mkdir -p $(LIBDIR) ;
	ar cru $@ $(TMPDIR)/obj/*/*.o

doc:
	@for dir in $(subdirs); do \
	echo "*************************************"; \
	echo "*    Creating Documentation for $$dir..."; \
	echo "*************************************"; \
	$(MAKE) -C $$dir/src doc; done;

clean:
	@for dir in $(subdirs); do \
	echo "*************************************"; \
	echo "*    Cleaning MarlinReco Package $$dir..."; \
	echo "*************************************"; \
	$(MAKE) -C $$dir/src clean; done; \
	rm -f $(LIB)

distclean:
	@for dir in $(subdirs); do \
	echo "*************************************"; \
	echo "*    Cleaning MarlinReco Package $$dir..."; \
	echo "*************************************"; \
	$(MAKE) -C $$dir/src distclean; done; \
	rm -rf $(LIB) $(TMPDIR) $(DOCDIR)
