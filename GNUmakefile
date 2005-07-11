# Collect subdirectories with packages
subdirs := $(patsubst %/src/GNUmakefile,%, $(wildcard */*/src/GNUmakefile) )
libdirs := $(patsubst %,%/lib, $(subdirs))



.PHONY: all lib clean $(subdirs) 


all: lib

#fixme: avoid to rebuild the library always?
lib: $(subdirs)
	@echo "*************************************"
	@echo "*    linking MarlinReco library     *"
	@echo "*************************************"
	@echo $(libdirs)
	mkdir -p lib
	mkdir -p build
	cd build && ls *.a| xargs -n1 ar x
	rm -f lib/libMarlinReco.a
	ar cr lib/libMarlinReco.a build/*.o
	rm -rf build
	ranlib lib/libMarlinReco.a





$(subdirs):
	@echo "*************************************"
	@echo "*   building subdir $@ ..."
	@echo "*************************************"
	$(MAKE) -C $@/src lib
	@mkdir -p build
	@cp $@/lib/lib*.a build


clean:
	for dir in $(subdirs); do \
	   $(MAKE) -C $$dir/src clean; \
	done

