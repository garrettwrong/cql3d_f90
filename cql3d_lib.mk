MODDIR=$(OBJ)/mod
LIBDIR=$(OBJ)/lib
DYLDIR=$(OBJ)/dylib
EXEDIR=$(OBJ)/exe

## TRANSP build stuff, should set FC90 and FFLAGS.
SHR = $(CODESYSDIR)/source/misc/makeflags.mk
include $(SHR)
export

all: libcql3d_lib libmpi_cql3d_lib

libcql3d_lib: pkg
	@test -d $(DYLDIR) || mkdir -p $(DYLDIR)
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@test -d $(MODDIR) || mkdir -p $(MODDIR)
	@test -d $(EXEDIR) || mkdir -p $(EXEDIR)
	@cp libxcql3d.so $(DYLDIR)/libcql3d_lib.so
	@cp libxcql3d.a $(LIBDIR)/cql3d_lib.a
	@cp *.mod $(MODDIR)
	@cp xcql3d $(EXEDIR)/xcql3d

libmpi_cql3d_lib: libcql3d_lib
	$(MAKE) -f makefile_gfortran64.CentOS7 -j mpi
	@cp libmpi_xcql3d.so $(DYLDIR)/libmpi_cql3d_lib.so
	@cp libmpi_xcql3d.a $(LIBDIR)/mpi_cql3d_lib.a
	@cp mpi_xcql3d $(EXEDIR)/mpi_xcql3d

pkg:
	$(MAKE) -f makefile_gfortran64.CentOS7 -j all

clean:
	$(MAKE) -f makefile_gfortran64.CentOS7 clean
realclean: clean

.PHONY: realclean
.PHONY: clean
.PHONY: pkg
.PHONY: libcql3d_lib
.PHONY: all
