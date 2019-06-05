FROM centos:7

#setup the base sysem
RUN yum upgrade -y && \
    yum group install "Development Tools" "X Window System" -y && \
    yum install gcc-gfortran netcdf netcdf-fortran netcdf-devel netcdf-fortran-devel emacs libX11-devel -y

# setup enough PGPLOT to work with cql3d
RUN curl -o pgplot522.tar.gz ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot522.tar.gz && \
    tar -xvzf pgplot522.tar.gz && \
    cd pgplot && \
    mv drivers.list drivers.list.bak && \
    echo """  NUDRIV 0 /NULL      Null device (no output)                           Std F77 \
  PSDRIV 1 /PS        PostScript printers, monochrome, landscape        Std F77 \
  PSDRIV 2 /VPS       Postscript printers, monochrome, portrait         Std F77 \
  PSDRIV 3 /CPS       PostScript printers, color, landscape             Std F77 \
  PSDRIV 4 /VCPS      PostScript printers, color, portrait              Std F77 \
  XWDRIV 1 /XWINDOW   Workstations running X Window System               C """ > drivers.list && \
    ./makemake $PWD linux g77_gcc && \
    sed -i 's/=g77/=gfortran/g'  makefile && \
    sed -i 's#/usr/X11R6/lib#/usr/lib64/#g' makefile && \
    sed -i 's#/usr/X11R6/include#/usr/include#g' makefile && \
    make && \
    export PGPLOT_DIR=$HOME/pgplot && \
    export LD_LIBRARY_PATH=$PGPLOT:$LD_LIBRARY_PATH

# setup Remik's nccmp tool
RUN curl -L -o nccmp-1.8.3.1.tar.gz https://gitlab.com/remikz/nccmp/-/archive/1.8.3.1/nccmp-1.8.3.1.tar.gz && \
    tar -xvzf nccmp-1.8.3.1.tar.gz && \
    cd nccmp-1.8.3.1 && \
    chmod +x ./configure # wut && \
    ./configure && \
    make && \
    make install
