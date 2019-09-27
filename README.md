# cql3d_f90

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/garrettwrong/cql3d_f90?include_prereleases)
![master pipeline](https://gitlab.com/geebdubya/cql3d_f90/badges/master/pipeline.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3463591.svg)](https://doi.org/10.5281/zenodo.3463591)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

## About

cql3d_f90 is a derived work of the FORTRAN IV/F77 CompX (compxco) CQL3D,
    a relativistic Collisional/QuasiLinear Fokker-Planck code,
        by R. W. Harvey et al.

Per Harvey, compx\_archive/compx\_LICENSE:

```
....
    The primary reference for the code is:
    ``The CQL3D Fokker-Planck Code'',  R.W. Harvey and M.G. McCoy,
    Proc. of IAEA Technical Committee Meeting on Advances in Simulation
    and Modeling of Thermonuclear Plasmas, Montreal, 1992, p. 489-526,
    IAEA, Vienna (1993), available through NTIS/DOC (National Technical
    Information Service, U.S. Dept. of Commerce), Order No. DE93002962.
    See also, http://www.compxco.com/cql3d.html, CQL3D Manual.

    The CQL3D code was built starting from the 2D single flux surface
    Coulomb collision solver in the CQL code: G.D. Kerbel and M.G. McCoy,
    "Kinetic theory and simulation of multispecies plasmas in tokamaks
    excited with electromagnetic waves in the ion-cyclotron range of
    frequencies", Phys. Fluids 28, 3629 (1985).
...
                 GNU Public License Distribution Only
...
```

Harvey is an excellent scientific author, and I suggest you read his summaries of
the physics behind this code here `https://www.compxco.com/cql3d.html`.

As in the original code, the standalone executable expects to read the `cqlinput`
file.  You can find exhausting notes on the options in cqlinput_help and further
in the compx_archive, particularly in a_change, and the many READMEs.

## Why does cql3d_f90 exist?

Unfortunately, the original code required a little more attention than could be
delivered by CompX.  Expecting the current and next generation of CQL users will
not have the intimate familiarity with the code required, the TRANSP group at
PPPL, initially working with CompX, have begun this cql3d_f90 project with the
main goals of a safer and easier to maintain code that can more easily integrate
with other projects. So far we have accomplished much in that direction:

     * Bulk rewrite to F90 using modules and derived types.
     * Removal of common blocks (spare freya).
     * Improved numerical resolution, and bug fixes including bitwise
          reproducibility, culminating in fully automated testing.
     * MPI Re-Implementation and fixes MPI.
     * Static and shared library based API implementation supporting
          further integration with other softwares.
     * Consolidation of makefiles.
     * CentOS7 distribution via Docker container.
     * Makefile option to compile and run without pgplot library.
     * ...and more... see cql3d_f90 git...

Currently the library integration is under test via TRANSP at PPPL.  All of the
test cases in the `ci` directory have been reviewed with CompX to be physically
equivalent.. Between cql3d_f90 and legacy CQL3D, the cql3d_f90 results are
generally believed to be numerically superiour, mainly based on demonstrable
improved stability and resolution.  Although we have not tested all options
exhaustively, our reviews with CompX suggest we have decent coverage.

At this time there has been no changes or additions targeting underlying physics.

There is still much work to be done to improve these codes, and I
hope that having a more collaboritive open project with some automation inspires
the community to do so.  I am also hopeful that we will soon have the complete
support of CompX now that the bulk derivation has been produced.

## Getting Started

Standalone cql3d_f90 should build and run natively on most versions of Linux, 
and also on various OSX systems.  Additionally, I have provided a Docker container
that should allow you to easily run or develop this code in _any_ 
environment supporting Docker.  The Docker image includes all required dependent
packages prebuilt for you, so it is totally turn-key.

Get a copy of the code:

```
git clone git@github.com:garrettwrong/cql3d_f90.git
cd cql3d_f90
```

Then pick the section that is most relevant to you below.  If you have trouble,
with the dependencies, or are on windows, try using the distributed docker
release, which should give you a consistent/known environment.

### Linux

The `makefile_gfortran64.CentOS7` is the general starting place for linux.  As
the name implies it was developed on CentOS7 for gfortran, but easily changed.

Assuming you can build codes using gfortran, along with having Netcdf and pgplot:

```
make -f makefile_gfortran64.CentOS7
```

You may optionally add a -j flag to the make line.

#### MPI

The mpi version expects you have MPI installed, resolving to `mpifort`, and with
required paths/libraries configured.

```
make -f makefile_gfortran64.CentOS7 mpi
mpirun -n 3 ./mpi_xcql3d
```

Obviously you might want to use more than 3 processors.

#### Fancy compilers (Intel)

For ifort 2015-2019 I was simply able to change the makefile line:

```
#from
COMPILER  = gfortran
#to
COMPILER  = ifort
```

Easy peasy.  You may wish to change some flags on your own.

There is a known optimization bug in some iforts,
including 2015, that manifests as stack corruption induced segfault.
 If you think that is affecting your run, it is avoidable with the flag
 `-heap-arrays`, which disables that intel optimization.

#### Don't like PGPLOT?

Try the following or setting equivalent in your makefile:

`export NOPGPLOT=TRUE`

### OSX

Much of this code was actually developed on a mac.  Using `homebrew` I installed
 dependencies:

```
brew cask install gfortran
brew install netcdf pgplot
brew install open-mpi
```

then, when you beleive you can build F codes:
```
make -f makefile_gfortran64.OSX -j
./xcql3d
```

### Using Docker

You can run a known image, the same as used for testing, by using Docker.

This image comes with all required software installed.

First pull the image (this may take a few minutes the first time):

```
docker pull registry.gitlab.com/geebdubya/cql3d_f90/centos-7:latest
```

Then you can make and run the software using the container:


```
docker run -v `pwd`:`pwd` -w `pwd` make -f makefile_gfortran64.CentOS7
docker run -v `pwd`:`pwd` -w `pwd`   ./xcql3d
```

similarly for mpi:

```
docker run -v `pwd`:`pwd` -w `pwd` make -f makefile_gfortran64.CentOS7 mpi
docker run -v `pwd`:`pwd` -w `pwd`  mpirun -n 3 ./mpi_xcql3d
```

If you want to alter the environment, you can start with the
provided file `ci/Dockerfile.centos7` used in this distribution, building and
 tagging as you would any docker image.


## Support

This project was supported by TRANSP development at the
Princeton Plasma Physics Laboratory via
`U.S. Department of Energy (DE-AC02-09CH11466)`.

