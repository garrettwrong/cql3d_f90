MODDIR=$(OBJ)/mod
LIBDIR=$(OBJ)/lib
DYLDIR=$(OBJ)/dylib
EXEDIR=$(OBJ)/exe

# These are needed to invoke cql3d from library as done in a_cqlp
# One day these can migrate to some sort of single api module.
# For now this is not so bad. I think it even works.
PUBMODS=$(wildcard *.mod)

## TRANSP build stuff, should set FC90 and FFLAGS.
SHR = $(CODESYSDIR)/source/misc/makeflags.mk
include $(SHR)
export

libcql3d_lib: pkg
	@test -d $(DYLDIR) || mkdir -p $(DYLDIR)
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@test -d $(MODDIR) || mkdir -p $(MODDIR)
	@test -d $(EXEDIR) || mkdir -p $(EXEDIR)
	@cp libxcql3d.so $(DYLDIR)/libcql3d_lib.so
	@cp libxcql3d.a $(LIBDIR)/cql3d_lib.a
	@cp $(PUBMODS) $(MODDIR)
	@cp xcql3d $(EXEDIR)/xcql3d

pkg:
	$(MAKE) -f makefile_gfortran64.CentOS7 -j all

clean:
	$(MAKE) -f makefile_gfortran64.CentOS7 clean
realclean: clean

.PHONY: realclean
.PHONY: clean
.PHONY: pkg
.PHONY: libcql3d_lib


