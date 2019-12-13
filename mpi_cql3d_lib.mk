MODDIR=$(OBJ)/mod
MPIMODDIR=$(OBJ)/mpi_mod
LIBDIR=$(OBJ)/lib
DYLDIR=$(OBJ)/dylib
EXEDIR=$(OBJ)/exe
BUILDDIR=$(OBJ)/obj/mpi_cql3d_lib

## TRANSP build stuff, should set FC90 and FFLAGS.
SHR = $(CODESYSDIR)/source/misc/makeflags.mk
include $(SHR)
export

libmpi_cql3d_lib: install

pkg:
	@test -d $(DYLDIR) || mkdir -p $(DYLDIR)
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@test -d $(MODDIR) || mkdir -p $(MODDIR)
	@test -d $(MPIMODDIR) || mkdir -p $(MPIMODDIR)
	@test -d $(EXEDIR) || mkdir -p $(EXEDIR)
	@test -d $(BUILDDIR) || mkdir -p $(BUILDDIR)
	$(MAKE) -f makefile_gfortran64.CentOS7 -j mpi

install: pkg
	$(MAKE) -f makefile_gfortran64.CentOS7 -j mpi
	@cp $(BUILDDIR)/libmpi_xcql3d.so $(DYLDIR)/libmpi_cql3d_lib.so
	@cp $(BUILDDIR)/libmpi_xcql3d.a $(LIBDIR)/mpi_cql3d_lib.a
	@cp $(BUILDDIR)/cqlmpilib_mod.mod $(MPIMODDIR)
	@cp $(BUILDDIR)/mpi_xcql3d $(EXEDIR)/mpi_xcql3d

clean:
	$(MAKE) -f makefile_gfortran64.CentOS7 clean
realclean: clean

.PHONY: realclean
.PHONY: clean
.PHONY: install
.PHONY: pkg
.PHONY: libmpi_cql3d_lib
