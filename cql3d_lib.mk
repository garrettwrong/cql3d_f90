MODDIR=$(OBJ)/mod
LIBDIR=$(OBJ)/lib
DYLDIR=$(OBJ)/dylib
EXEDIR=$(OBJ)/exe

## TRANSP build stuff, should set FC90 and FFLAGS.
SHR = $(CODESYSDIR)/source/misc/makeflags.mk
include $(SHR)
export

libcql3d_lib: install
	@echo CQL3D Post Install Cleanup
	$(MAKE) -f makefile_gfortran64.CentOS7 clean

pkg: clean
	$(MAKE) -f makefile_gfortran64.CentOS7 -j all

install: pkg
	@test -d $(DYLDIR) || mkdir -p $(DYLDIR)
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@test -d $(MODDIR) || mkdir -p $(MODDIR)
	@test -d $(EXEDIR) || mkdir -p $(EXEDIR)
	@cp libxcql3d.so $(DYLDIR)/libcql3d_lib.so
	@cp libxcql3d.a $(LIBDIR)/cql3d_lib.a
	@cp *.mod $(MODDIR)
	@cp xcql3d $(EXEDIR)/xcql3d

clean:
	$(MAKE) -f makefile_gfortran64.CentOS7 clean
realclean: clean

.PHONY: realclean
.PHONY: clean
.PHONY: install
.PHONY: pkg
.PHONY: libcql3d_lib
