# I have used h77 to denote f77 formatted header
# while .h is f90.
# as the need for f77 goes away, then will call everything .h again...
# YuP[2019-05-31]: Also needed mpilib.h for MPI compilations/runs:
#      see mpins_par.f, mpilib.f90 and doparallel.py

INCLUDES  = frcomm.h77 frname.h77 frname_decl.h77 \
	frname.h frname_decl.h trans.h wpadvnc.h
