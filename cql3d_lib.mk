MODDIR=$(OBJ)/mod
LIBDIR=$(OBJ)/lib
DYLDIR=$(OBJ)/dylib
EXEDIR=$(OBJ)/exe
BUILDDIR=$(OBJ)/obj/cql3d_lib

## TRANSP build stuff, should set FC90 and FFLAGS.
SHR = $(CODESYSDIR)/source/misc/makeflags.mk
include $(SHR)
export

libcql3d_lib: install

pkg:
	@test -d $(DYLDIR) || mkdir -p $(DYLDIR)
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@test -d $(MODDIR) || mkdir -p $(MODDIR)
	@test -d $(EXEDIR) || mkdir -p $(EXEDIR)
	@test -d $(BUILDDIR) || mkdir -p $(BUILDDIR)
	$(MAKE) -f makefile_gfortran64.CentOS7 -j all

install: pkg
	@cp $(BUILDDIR)/libxcql3d.so $(DYLDIR)/libcql3d_lib.so
	@cp $(BUILDDIR)/libxcql3d.a $(LIBDIR)/cql3d_lib.a
	@cp $(BUILDDIR)/*.mod $(MODDIR)
	@cp $(BUILDDIR)/xcql3d $(EXEDIR)/xcql3d

clean:
	$(MAKE) -f makefile_gfortran64.CentOS7 clean
realclean: clean

.PHONY: realclean
.PHONY: clean
.PHONY: install
.PHONY: pkg
.PHONY: libcql3d_lib
