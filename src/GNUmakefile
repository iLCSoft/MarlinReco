ifndef MARLINWORKDIR
 MARLINWORKDIR=$(MARLIN)
 export MARLINWORKDIR
endif

PROGNAME = MarlinReco

subdirs := $(patsubst %/src/GNUmakefile,%, $(wildcard $(MARLINWORKDIR)/packages/MarlinReco/*/*/src/GNUmakefile))

packages := $(foreach file, $(subdirs), $(shell basename $(file)))

libs := $(patsubst %,$(MARLINWORKDIR)/tmp/$(PROGNAME)/lib/lib%.a, $(packages))

DOCDIR = $(MARLINWORKDIR)/doc/$(PROGNAME)

.PHONY: all lib doc clean distclean $(libs)

all: lib

lib: $(MARLINWORKDIR)/lib/lib$(PROGNAME).a

$(MARLINWORKDIR)/lib/lib$(PROGNAME).a : $(libs)
	@echo "*************************************"
	@echo "*    linking $(PROGNAME) library     *"
	@echo "*************************************"
	ar cru $@ $(MARLINWORKDIR)/tmp/$(PROGNAME)/obj/*/*.o

$(libs):
	@echo "*************************************"
	@echo "*   building library $@ ..."
	@echo "*************************************"
	$(MAKE) -C $(shell find ../ -name $(patsubst lib%.a,%,$(notdir $@)) -type d 2>/dev/null)/src lib

doc:
	@for dir in $(subdirs); do \
	echo "Creating documentation for $$dir..."; \
	$(MAKE) -C $$dir/src doc; done;

clean:
	@for dir in $(subdirs); do \
	echo "Clearing $$dir..."; \
	$(MAKE) -C $$dir/src clean; done; \
	rm -f $(MARLINWORKDIR)/lib/lib$(PROGNAME).a

distclean:
	rm -rf $(MARLINWORKDIR)/lib/lib$(PROGNAME).a $(MARLINWORKDIR)/tmp/$(PROGNAME) $(DOCDIR)

# end
