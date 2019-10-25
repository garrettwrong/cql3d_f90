MODDIR=$(OBJ)/mod
MPIMODDIR=$(OBJ)/mpi_mod
LIBDIR=$(OBJ)/lib
DYLDIR=$(OBJ)/dylib
EXEDIR=$(OBJ)/exe

## TRANSP build stuff, should set FC90 and FFLAGS.
SHR = $(CODESYSDIR)/source/misc/makeflags.mk
include $(SHR)
export

libmpi_cql3d_lib: | clean install clean

install: pkg
	@cp libmpi_xcql3d.so $(DYLDIR)/libmpi_cql3d_lib.so
	@cp libmpi_xcql3d.a $(LIBDIR)/mpi_cql3d_lib.a
	@cp cqlmpilib_mod.mod $(MPIMODDIR)
	@cp mpi_xcql3d $(EXEDIR)/mpi_xcql3d

pkg:
	$(MAKE) -f makefile_gfortran64.CentOS7 -j mpi

clean:
	$(MAKE) -f makefile_gfortran64.CentOS7 clean
realclean: clean

.PHONY: realclean
.PHONY: clean
.PHONY: install
.PHONY: pkg
.PHONY: libmpi_cql3d_lib
.PHONY: all
