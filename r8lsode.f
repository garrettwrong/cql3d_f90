c    Obtained from http://www.netlib.no/netlib/odepack/
c    Contains lsode and related subroutines, according to
c    file doc in odepack.
c    BobH,  990303.
c    A DOUBLE PRECISION VERSION.
c
c
c YuP[2017] Modifications: func. vnorm() is renamed into vnorm_lsode(),
c to avoid conflicts with a possible common-name function.
c Also, the common blocks /ls0001/ and /eh0001/ are renamed into
c /ls_lsode/ and /eh_lsode/


      subroutine lsode (f, neq, y, t, tout, itol, rtol, atol, itask,
     1            istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      use iso_c_binding, only : c_double
      use r8subs_mod, only : d1mach
      save
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      real(c_double) y, t, tout, rtol, atol, rwork
      dimension neq(*), y(*), rtol(*), atol(*), rwork(lrw), iwork(liw)
c-----------------------------------------------------------------------
c this is the march 30, 1987 version of
c lsode.. livermore solver for ordinary differential equations.
c this version is in double precision.
c
c lsode solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c lsode is a package based on the gear and gearb packages, and on the
c october 23, 1978 version of the tentative odepack user interface
c standard, with minor modifications.
c-----------------------------------------------------------------------
c reference..
c     alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c-----------------------------------------------------------------------
c author and contact.. alan c. hindmarsh,
c                      computing and mathematics research div., l-316
c                      lawrence livermore national laboratory
c                      livermore, ca 94550.
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsode package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including optional communication, nonstandard options,
c and instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first provide a subroutine of the form..
c               subroutine f (neq, t, y, ydot)
c               dimension y(neq), ydot(neq)
c which supplies the vector function f by loading ydot(i) with f(i).
c
c b. next determine (or guess) whether or not the problem is stiff.
c stiffness occurs when the jacobian matrix df/dy has an eigenvalue
c whose real part is negative and large in magnitude, compared to the
c reciprocal of the t span of interest.  if the problem is nonstiff,
c use a method flag mf = 10.  if it is stiff, there are four standard
c choices for mf, and lsode requires the jacobian matrix in some form.
c this matrix is regarded either as full (mf = 21 or 22),
c or banded (mf = 24 or 25).  in the banded case, lsode requires two
c half-bandwidth parameters ml and mu.  these are, respectively, the
c widths of the lower and upper parts of the band, excluding the main
c diagonal.  thus the band consists of the locations (i,j) with
c i-ml .le. j .le. i+mu, and the full bandwidth is ml+mu+1.
c
c c. if the problem is stiff, you are encouraged to supply the jacobian
c directly (mf = 21 or 24), but if this is not feasible, lsode will
c compute it internally by difference quotients (mf = 22 or 25).
c if you are supplying the jacobian, provide a subroutine of the form..
c               subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
c               dimension y(neq), pd(nrowpd,neq)
c which supplies df/dy by loading pd as follows..
c     for a full jacobian (mf = 21), load pd(i,j) with df(i)/dy(j),
c the partial derivative of f(i) with respect to y(j).  (ignore the
c ml and mu arguments in this case.)
c     for a banded jacobian (mf = 24), load pd(i-j+mu+1,j) with
c df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
c pd from the top down.
c     in either case, only nonzero elements need be loaded.
c
c d. write a main program which calls subroutine lsode once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsode.  on the first call to lsode, supply arguments as follows..
c f      = name of subroutine for right-hand side vector f.
c          this name must be declared external in calling program.
c neq    = number of first order ode-s.
c y      = array of initial values, of length neq.
c t      = the initial value of the independent variable.
c tout   = first point where output is desired (.ne. t).
c itol   = 1 or 2 according as atol (below) is a scalar or array.
c rtol   = relative tolerance parameter (scalar).
c atol   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be roughly less (in magnitude) than
c             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
c             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
c          thus the local error test passes if, in each component,
c          either the absolute error is less than atol (or atol(i)),
c          or the relative error is less than rtol.
c          use rtol = 0.0 for pure absolute error control, and
c          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
c          control.  caution.. actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c itask  = 1 for normal computation of output values of y at t = tout.
c istate = integer flag (input and output).  set istate = 1.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             20 + 16*neq                    for mf = 10,
c             22 +  9*neq + neq**2           for mf = 21 or 22,
c             22 + 10*neq + (2*ml + mu)*neq  for mf = 24 or 25.
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least..
c             20        for mf = 10,
c             20 + neq  for mf = 21, 22, 24, or 25.
c          if mf = 24 or 25, input in iwork(1),iwork(2) the lower
c          and upper half-bandwidths ml,mu.
c liw    = declared length of iwork (in user-s dimension).
c jac    = name of subroutine for jacobian matrix (mf = 21 or 24).
c          if used, this name must be declared external in calling
c          program.  if not used, pass a dummy name.
c mf     = method flag.  standard values are..
c          10 for nonstiff (adams) method, no jacobian used.
c          21 for stiff (bdf) method, user-supplied full jacobian.
c          22 for stiff method, internally generated full jacobian.
c          24 for stiff method, user-supplied banded jacobian.
c          25 for stiff method, internally generated banded jacobian.
c note that the main program must declare arrays y, rwork, iwork,
c and possibly atol.
c
c e. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsode was successful, negative otherwise.
c          -1 means excess work done on this call (perhaps wrong mf).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of mf or tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c
c f. to continue the integration after a successful return, simply
c reset tout and call lsode again.  no other parameters need be reset.
c
c-----------------------------------------------------------------------
c example problem.
c
c the following is a simple example problem, with the coding
c needed for its solution by lsode.  the problem is from chemical
c kinetics, and consists of the following three rate equations..
c     dy1/dt = -.04*y1 + 1.e4*y2*y3
c     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c     dy3/dt = 3.e7*y2**2
c on the interval from t = 0.0 to t = 4.e10, with initial conditions
c y1 = 1.0, y2 = y3 = 0.  the problem is stiff.
c
c the following coding solves this problem with lsode, using mf = 21
c and printing results at t = .4, 4., ..., 4.e10.  it uses
c itol = 2 and atol much smaller for y2 than y1 or y3 because
c y2 has much smaller values.
c at the end of the run, statistical quantities of interest are
c printed (see optional outputs in the full description below).
c
c     external fex, jex
c     double precision atol, rtol, rwork, t, tout, y
c     dimension y(3), atol(3), rwork(58), iwork(23)
c     neq = 3
c     y(1) = 1.d0
c     y(2) = 0.d0
c     y(3) = 0.d0
c     t = 0.d0
c     tout = .4d0
c     itol = 2
c     rtol = 1.d-4
c     atol(1) = 1.d-6
c     atol(2) = 1.d-10
c     atol(3) = 1.d-6
c     itask = 1
c     istate = 1
c     iopt = 0
c     lrw = 58
c     liw = 23
c     mf = 21
c     do 40 iout = 1,12
c       call lsode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
c    1     iopt,rwork,lrw,iwork,liw,jex,mf)
c       write(6,20)t,y(1),y(2),y(3)
c 20    format(7h at t =,e12.4,6h   y =,3e14.6)
c       if (istate .lt. 0) go to 80
c 40    tout = tout*10.d0
c     write(6,60)iwork(11),iwork(12),iwork(13)
c 60  format(/12h no. steps =,i4,11h  no. f-s =,i4,11h  no. j-s =,i4)
c     stop
c 80  write(6,90)istate
c 90  format(///22h error halt.. istate =,i3)
c     stop
c     end
c
c     subroutine fex (neq, t, y, ydot)
c     double precision t, y, ydot
c     dimension y(3), ydot(3)
c     ydot(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
c     ydot(3) = 3.d7*y(2)*y(2)
c     ydot(2) = -ydot(1) - ydot(3)
c     return
c     end
c
c     subroutine jex (neq, t, y, ml, mu, pd, nrpd)
c     double precision pd, t, y
c     dimension y(3), pd(nrpd,3)
c     pd(1,1) = -.04d0
c     pd(1,2) = 1.d4*y(3)
c     pd(1,3) = 1.d4*y(2)
c     pd(2,1) = .04d0
c     pd(2,3) = -pd(1,3)
c     pd(3,2) = 6.d7*y(2)
c     pd(2,2) = -pd(1,2) - pd(3,2)
c     return
c     end
c
c the output of this program (on a cdc-7600 in single precision)
c is as follows..
c
c   at t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
c   at t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
c   at t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
c   at t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
c   at t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
c   at t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
c   at t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
c   at t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
c   at t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
c   at t =  4.0000e+08   y =  5.494529e-06  2.197824e-11  9.999945e-01
c   at t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01
c   at t =  4.0000e+10   y = -7.170586e-08 -2.868234e-13  1.000000e+00
c
c   no. steps = 330  no. f-s = 405  no. j-s =  69
c-----------------------------------------------------------------------
c full description of user interface to lsode.
c
c the user interface to lsode consists of the following parts.
c
c i.   the call sequence to subroutine lsode, which is a driver
c      routine for the solver.  this includes descriptions of both
c      the call sequence arguments and of user-supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  descriptions of other routines in the lsode package that may be
c      (optionally) called by the user.  these provide the ability to
c      alter error message handling, save and restore the internal
c      common, and obtain specified derivatives of the solution y(t).
c
c iii. descriptions of common blocks to be declared in overlay
c      or similar environments, or to be saved when doing an interrupt
c      of the problem and continued solution later.
c
c iv.  description of two routines in the lsode package, either of
c      which the user may replace with his own version, if desired.
c      these relate to the measurement of errors.
c
c-----------------------------------------------------------------------
c part i.  call sequence.
c
c the call sequence parameters used for input only are
c     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
c and those used for both input and output are
c     y, t, istate.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine lsode to the user-s calling program.)
c
c the legality of input parameters will be thoroughly checked on the
c initial call for the problem, but not checked thereafter unless a
c change in input parameters is flagged by istate = 3 on input.
c
c the descriptions of the call arguments are as follows.
c
c f      = the name of the user-supplied subroutine defining the
c          ode system.  the system must be put in the first-order
c          form dy/dt = f(t,y), where f is a vector-valued function
c          of the scalar t and the vector y.  subroutine f is to
c          compute the function f.  it is to have the form
c               subroutine f (neq, t, y, ydot)
c               dimension y(1), ydot(1)
c          where neq, t, and y are input, and the array ydot = f(t,y)
c          is output.  y and ydot are arrays of length neq.
c          (in the dimension statement above, 1 is a dummy
c          dimension.. it can be replaced by any value.)
c          subroutine f should not alter y(1),...,y(neq).
c          f must be declared external in the calling program.
c
c          subroutine f may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in f) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y below.
c
c          if quantities computed in the f routine are needed
c          externally to lsode, an extra call to f should be made
c          for this purpose, for consistent and accurate results.
c          if only the derivative dy/dt is needed, use intdy instead.
c
c neq    = the size of the ode system (number of first order
c          ordinary differential equations).  used only for input.
c          neq may be decreased, but not increased, during the problem.
c          if neq is decreased (with istate = 3 on input), the
c          remaining components of y should be left undisturbed, if
c          these are to be accessed in f and/or jac.
c
c          normally, neq is a scalar, and it is generally referred to
c          as a scalar in this user interface description.  however,
c          neq may be an array, with neq(1) set to the system size.
c          (the lsode package accesses only neq(1).)  in either case,
c          this parameter is passed as the neq argument in all calls
c          to f and jac.  hence, if it is an array, locations
c          neq(2),... may be used to store other integer data and pass
c          it to f and/or jac.  subroutines f and/or jac must include
c          neq in a dimension statement in that case.
c
c y      = a real array for the vector of dependent variables, of
c          length neq or more.  used for both input and output on the
c          first call (istate = 1), and only for output on other calls.
c          on the first call, y must contain the vector of initial
c          values.  on output, y contains the computed solution vector,
c          evaluated at t.  if desired, the y array may be used
c          for other purposes between calls to the solver.
c
c          this array is passed as the y argument in all calls to
c          f and jac.  hence its length may exceed neq, and locations
c          y(neq+1),... may be used to store other real data and
c          pass it to f and/or jac.  (the lsode package accesses only
c          y(1),...,y(neq).)
c
c t      = the independent variable.  on input, t is used only on the
c          first call, as the initial point of the integration.
c          on output, after each call, t is the value at which a
c          computed solution y is evaluated (usually the same as tout).
c          on an error return, t is the farthest point reached.
c
c tout   = the next value of t at which a computed solution is desired.
c          used only for input.
c
c          when starting the problem (istate = 1), tout may be equal
c          to t for one call, then should .ne. t for the next call.
c          for the initial t, an input value of tout .ne. t is used
c          in order to determine the direction of the integration
c          (i.e. the algebraic sign of the step sizes) and the rough
c          scale of the problem.  integration in either direction
c          (forward or backward in t) is permitted.
c
c          if itask = 2 or 5 (one-step modes), tout is ignored after
c          the first call (i.e. the first call with tout .ne. t).
c          otherwise, tout is required on every call.
c
c          if itask = 1, 3, or 4, the values of tout need not be
c          monotone, but a value of tout which backs up is limited
c          to the current internal t interval, whose endpoints are
c          tcur - hu and tcur (see optional outputs, below, for
c          tcur and hu).
c
c itol   = an indicator for the type of error control.  see
c          description below under atol.  used only for input.
c
c rtol   = a relative error tolerance parameter, either a scalar or
c          an array of length neq.  see description below under atol.
c          input only.
c
c atol   = an absolute error tolerance parameter, either a scalar or
c          an array of length neq.  input only.
c
c             the input parameters itol, rtol, and atol determine
c          the error control performed by the solver.  the solver will
c          control the vector e = (e(i)) of estimated local errors
c          in y, according to an inequality of the form
c                      rms-norm of ( e(i)/ewt(i) )   .le.   1,
c          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
c          and the rms-norm (root-mean-square norm) here is
c          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
c          is a vector of weights which must always be positive, and
c          the values of rtol and atol should all be non-negative.
c          the following table gives the types (scalar/array) of
c          rtol and atol, and the corresponding form of ewt(i).
c
c             itol    rtol       atol          ewt(i)
c              1     scalar     scalar     rtol*abs(y(i)) + atol
c              2     scalar     array      rtol*abs(y(i)) + atol(i)
c              3     array      scalar     rtol(i)*abs(y(i)) + atol
c              4     array      array      rtol(i)*abs(y(i)) + atol(i)
c
c          when either of these parameters is a scalar, it need not
c          be dimensioned in the user-s calling program.
c
c          if none of the above choices (with itol, rtol, and atol
c          fixed throughout the problem) is suitable, more general
c          error controls can be obtained by substituting
c          user-supplied routines for the setting of ewt and/or for
c          the norm calculation.  see part iv below.
c
c          if global errors are to be estimated by making a repeated
c          run on the same problem with smaller tolerances, then all
c          components of rtol and atol (i.e. of ewt) should be scaled
c          down uniformly.
c
c itask  = an index specifying the task to be performed.
c          input only.  itask has the following values and meanings.
c          1  means normal computation of output values of y(t) at
c             t = tout (by overshooting and interpolating).
c          2  means take one step only and return.
c          3  means stop at the first internal mesh point at or
c             beyond t = tout and return.
c          4  means normal computation of output values of y(t) at
c             t = tout but without overshooting t = tcrit.
c             tcrit must be input as rwork(1).  tcrit may be equal to
c             or beyond tout, but not behind it in the direction of
c             integration.  this option is useful if the problem
c             has a singularity at or beyond t = tcrit.
c          5  means take one step, without passing tcrit, and return.
c             tcrit must be input as rwork(1).
c
c          note..  if itask = 4 or 5 and the solver reaches tcrit
c          (within roundoff), it will return t = tcrit (exactly) to
c          indicate this (unless itask = 4 and tout comes before tcrit,
c          in which case answers at t = tout are returned first).
c
c istate = an index used for input and output to specify the
c          the state of the calculation.
c
c          on input, the values of istate are as follows.
c          1  means this is the first call for the problem
c             (initializations will be done).  see note below.
c          2  means this is not the first call, and the calculation
c             is to continue normally, with no change in any input
c             parameters except possibly tout and itask.
c             (if itol, rtol, and/or atol are changed between calls
c             with istate = 2, the new values will be used but not
c             tested for legality.)
c          3  means this is not the first call, and the
c             calculation is to continue normally, but with
c             a change in input parameters other than
c             tout and itask.  changes are allowed in
c             neq, itol, rtol, atol, iopt, lrw, liw, mf, ml, mu,
c             and any of the optional inputs except h0.
c             (see iwork description for ml and mu.)
c          note..  a preliminary call with tout = t is not counted
c          as a first call here, as no initialization or checking of
c          input is done.  (such a call is sometimes useful for the
c          purpose of outputting the initial conditions.)
c          thus the first call for which tout .ne. t requires
c          istate = 1 on input.
c
c          on output, istate has the following values and meanings.
c           1  means nothing was done, as tout was equal to t with
c              istate = 1 on input.  (however, an internal counter was
c              set to detect and prevent repeated calls of this type.)
c           2  means the integration was performed successfully.
c          -1  means an excessive amount of work (more than mxstep
c              steps) was done on this call, before completing the
c              requested task, but the integration was otherwise
c              successful as far as t.  (mxstep is an optional input
c              and is normally 500.)  to continue, the user may
c              simply reset istate to a value .gt. 1 and call again
c              (the excess work step counter will be reset to 0).
c              in addition, the user may increase mxstep to avoid
c              this error return (see below on optional inputs).
c          -2  means too much accuracy was requested for the precision
c              of the machine being used.  this was detected before
c              completing the requested task, but the integration
c              was successful as far as t.  to continue, the tolerance
c              parameters must be reset, and istate must be set
c              to 3.  the optional output tolsf may be used for this
c              purpose.  (note.. if this condition is detected before
c              taking any steps, then an illegal input return
c              (istate = -3) occurs instead.)
c          -3  means illegal input was detected, before taking any
c              integration steps.  see written message for details.
c              note..  if the solver detects an infinite loop of calls
c              to the solver with illegal input, it will cause
c              the run to stop.
c          -4  means there were repeated error test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              the problem may have a singularity, or the input
c              may be inappropriate.
c          -5  means there were repeated convergence test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              this may be caused by an inaccurate jacobian matrix,
c              if one is being used.
c          -6  means ewt(i) became zero for some i during the
c              integration.  pure relative error control (atol(i)=0.0)
c              was requested on a variable which has now vanished.
c              the integration was successful as far as t.
c
c          note..  since the normal output value of istate is 2,
c          it does not need to be reset for normal continuation.
c          also, since a negative input value of istate will be
c          regarded as illegal, a negative output value requires the
c          user to change it, and possibly other inputs, before
c          calling the solver again.
c
c iopt   = an integer flag to specify whether or not any optional
c          inputs are being used on this call.  input only.
c          the optional inputs are listed separately below.
c          iopt = 0 means no optional inputs are being used.
c                   default values will be used in all cases.
c          iopt = 1 means one or more optional inputs are being used.
c
c rwork  = a real working array (double precision).
c          the length of rwork must be at least
c             20 + nyh*(maxord + 1) + 3*neq + lwm    where
c          nyh    = the initial value of neq,
c          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
c                   smaller value is given as an optional input),
c          lwm   = 0             if miter = 0,
c          lwm   = neq**2 + 2    if miter is 1 or 2,
c          lwm   = neq + 2       if miter = 3, and
c          lwm   = (2*ml+mu+1)*neq + 2 if miter is 4 or 5.
c          (see the mf description for meth and miter.)
c          thus if maxord has its default value and neq is constant,
c          this length is..
c             20 + 16*neq                  for mf = 10,
c             22 + 16*neq + neq**2         for mf = 11 or 12,
c             22 + 17*neq                  for mf = 13,
c             22 + 17*neq + (2*ml+mu)*neq  for mf = 14 or 15,
c             20 +  9*neq                  for mf = 20,
c             22 +  9*neq + neq**2         for mf = 21 or 22,
c             22 + 10*neq                  for mf = 23,
c             22 + 10*neq + (2*ml+mu)*neq  for mf = 24 or 25.
c          the first 20 words of rwork are reserved for conditional
c          and optional inputs and optional outputs.
c
c          the following word in rwork is a conditional input..
c            rwork(1) = tcrit = critical value of t which the solver
c                       is not to overshoot.  required if itask is
c                       4 or 5, and ignored otherwise.  (see itask.)
c
c lrw    = the length of the array rwork, as declared by the user.
c          (this will be checked by the solver.)
c
c iwork  = an integer work array.  the length of iwork must be at least
c             20        if miter = 0 or 3 (mf = 10, 13, 20, 23), or
c             20 + neq  otherwise (mf = 11, 12, 14, 15, 21, 22, 24, 25).
c          the first few words of iwork are used for conditional and
c          optional inputs and optional outputs.
c
c          the following 2 words in iwork are conditional inputs..
c            iwork(1) = ml     these are the lower and upper
c            iwork(2) = mu     half-bandwidths, respectively, of the
c                       banded jacobian, excluding the main diagonal.
c                       the band is defined by the matrix locations
c                       (i,j) with i-ml .le. j .le. i+mu.  ml and mu
c                       must satisfy  0 .le.  ml,mu  .le. neq-1.
c                       these are required if miter is 4 or 5, and
c                       ignored otherwise.  ml and mu may in fact be
c                       the band parameters for a matrix to which
c                       df/dy is only approximately equal.
c
c liw    = the length of the array iwork, as declared by the user.
c          (this will be checked by the solver.)
c
c note..  the work arrays must not be altered between calls to lsode
c for the same problem, except possibly for the conditional and
c optional inputs, and except for the last 3*neq words of rwork.
c the latter space is used for internal scratch space, and so is
c available for use by the user outside lsode between calls, if
c desired (but not for use by f or jac).
c
c jac    = the name of the user-supplied routine (miter = 1 or 4) to
c          compute the jacobian matrix, df/dy, as a function of
c          the scalar t and the vector y.  it is to have the form
c               subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
c               dimension y(1), pd(nrowpd,1)
c          where neq, t, y, ml, mu, and nrowpd are input and the array
c          pd is to be loaded with partial derivatives (elements of
c          the jacobian matrix) on output.  pd must be given a first
c          dimension of nrowpd.  t and y have the same meaning as in
c          subroutine f.  (in the dimension statement above, 1 is a
c          dummy dimension.. it can be replaced by any value.)
c               in the full matrix case (miter = 1), ml and mu are
c          ignored, and the jacobian is to be loaded into pd in
c          columnwise manner, with df(i)/dy(j) loaded into pd(i,j).
c               in the band matrix case (miter = 4), the elements
c          within the band are to be loaded into pd in columnwise
c          manner, with diagonal lines of df/dy loaded into the rows
c          of pd.  thus df(i)/dy(j) is to be loaded into pd(i-j+mu+1,j).
c          ml and mu are the half-bandwidth parameters (see iwork).
c          the locations in pd in the two triangular areas which
c          correspond to nonexistent matrix elements can be ignored
c          or loaded arbitrarily, as they are overwritten by lsode.
c               jac need not provide df/dy exactly.  a crude
c          approximation (possibly with a smaller bandwidth) will do.
c               in either case, pd is preset to zero by the solver,
c          so that only the nonzero elements need be loaded by jac.
c          each call to jac is preceded by a call to f with the same
c          arguments neq, t, and y.  thus to gain some efficiency,
c          intermediate quantities shared by both calculations may be
c          saved in a user common block by f and not recomputed by jac,
c          if desired.  also, jac may alter the y array, if desired.
c          jac must be declared external in the calling program.
c               subroutine jac may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in jac) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y above.
c
c mf     = the method flag.  used only for input.  the legal values of
c          mf are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, and 25.
c          mf has decimal digits meth and miter.. mf = 10*meth + miter.
c          meth indicates the basic linear multistep method..
c            meth = 1 means the implicit adams method.
c            meth = 2 means the method based on backward
c                     differentiation formulas (bdf-s).
c          miter indicates the corrector iteration method..
c            miter = 0 means functional iteration (no jacobian matrix
c                      is involved).
c            miter = 1 means chord iteration with a user-supplied
c                      full (neq by neq) jacobian.
c            miter = 2 means chord iteration with an internally
c                      generated (difference quotient) full jacobian
c                      (using neq extra calls to f per df/dy value).
c            miter = 3 means chord iteration with an internally
c                      generated diagonal jacobian approximation.
c                      (using 1 extra call to f per df/dy evaluation).
c            miter = 4 means chord iteration with a user-supplied
c                      banded jacobian.
c            miter = 5 means chord iteration with an internally
c                      generated banded jacobian (using ml+mu+1 extra
c                      calls to f per df/dy evaluation).
c          if miter = 1 or 4, the user must supply a subroutine jac
c          (the name is arbitrary) as described above under jac.
c          for other values of miter, a dummy argument can be used.
c-----------------------------------------------------------------------
c optional inputs.
c
c the following is a list of the optional inputs provided for in the
c call sequence.  (see also part ii.)  for each such input variable,
c this table lists its name as used in this documentation, its
c location in the call sequence, its meaning, and the default value.
c the use of any of these inputs requires iopt = 1, and in that
c case all of these inputs are examined.  a value of zero for any
c of these optional inputs will cause the default value to be used.
c thus to use a subset of the optional inputs, simply preload
c locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
c then set those of interest to nonzero values.
c
c name    location      meaning and default value
c
c h0      rwork(5)  the step size to be attempted on the first step.
c                   the default value is determined by the solver.
c
c hmax    rwork(6)  the maximum absolute step size allowed.
c                   the default value is infinite.
c
c hmin    rwork(7)  the minimum absolute step size allowed.
c                   the default value is 0.  (this lower bound is not
c                   enforced on the final step before reaching tcrit
c                   when itask = 4 or 5.)
c
c maxord  iwork(5)  the maximum order to be allowed.  the default
c                   value is 12 if meth = 1, and 5 if meth = 2.
c                   if maxord exceeds the default value, it will
c                   be reduced to the default value.
c                   if maxord is changed during the problem, it may
c                   cause the current order to be reduced.
c
c mxstep  iwork(6)  maximum number of (internally defined) steps
c                   allowed during one call to the solver.
c                   the default value is 500.
c
c mxhnil  iwork(7)  maximum number of messages printed (per problem)
c                   warning that t + h = t on a step (h = step size).
c                   this must be positive to result in a non-default
c                   value.  the default value is 10.
c-----------------------------------------------------------------------
c optional outputs.
c
c as optional additional output from lsode, the variables listed
c below are quantities related to the performance of lsode
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on any successful return from lsode, and on any return with
c istate = -1, -2, -4, -5, or -6.  on an illegal input return
c (istate = -3), they will be unchanged from their existing values
c (if any), except possibly for tolsf, lenrw, and leniw.
c on any error return, outputs relevant to the error will be defined,
c as noted below.
c
c name    location      meaning
c
c hu      rwork(11) the step size in t last used (successfully).
c
c hcur    rwork(12) the step size to be attempted on the next step.
c
c tcur    rwork(13) the current value of the independent variable
c                   which the solver has actually reached, i.e. the
c                   current internal mesh point in t.  on output, tcur
c                   will always be at least as far as the argument
c                   t, but may be farther (if interpolation was done).
c
c tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
c                   computed when a request for too much accuracy was
c                   detected (istate = -3 if detected at the start of
c                   the problem, istate = -2 otherwise).  if itol is
c                   left unaltered but rtol and atol are uniformly
c                   scaled up by a factor of tolsf for the next call,
c                   then the solver is deemed likely to succeed.
c                   (the user may also ignore tolsf and alter the
c                   tolerance parameters in any other way appropriate.)
c
c nst     iwork(11) the number of steps taken for the problem so far.
c
c nfe     iwork(12) the number of f evaluations for the problem so far.
c
c nje     iwork(13) the number of jacobian evaluations (and of matrix
c                   lu decompositions) for the problem so far.
c
c nqu     iwork(14) the method order last used (successfully).
c
c nqcur   iwork(15) the order to be attempted on the next step.
c
c imxer   iwork(16) the index of the component of largest magnitude in
c                   the weighted local error vector ( e(i)/ewt(i) ),
c                   on an error return with istate = -4 or -5.
c
c lenrw   iwork(17) the length of rwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c leniw   iwork(18) the length of iwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c the following two arrays are segments of the rwork array which
c may also be of interest to the user as optional outputs.
c for each array, the table below gives its internal name,
c its base address in rwork, and its description.
c
c name    base address      description
c
c yh      21             the nordsieck history array, of size nyh by
c                        (nqcur + 1), where nyh is the initial value
c                        of neq.  for j = 0,1,...,nqcur, column j+1
c                        of yh contains hcur**j/factorial(j) times
c                        the j-th derivative of the interpolating
c                        polynomial currently representing the solution,
c                        evaluated at t = tcur.
c
c acor     lenrw-neq+1   array of size neq used for the accumulated
c                        corrections on each step, scaled on output
c                        to represent the estimated local error in y
c                        on the last step.  this is the vector e in
c                        the description of the error control.  it is
c                        defined only on a successful return from lsode.
c
c-----------------------------------------------------------------------
c part ii.  other routines callable.
c
c the following are optional calls which the user may make to
c gain additional capabilities in conjunction with lsode.
c (the routines xsetun and xsetf are designed to conform to the
c slatec error handling package.)
c
c     form of call                  function
c   call xsetun(lun)          set the logical unit number, lun, for
c                             output of messages from lsode, if
c                             the default is not desired.
c                             the default value of lun is 6.
c
c   call xsetf(mflag)         set a flag to control the printing of
c                             messages by lsode.
c                             mflag = 0 means do not print. (danger..
c                             this risks losing valuable information.)
c                             mflag = 1 means print (the default).
c
c                             either of the above calls may be made at
c                             any time and will take effect immediately.
c
c   call srcom(rsav,isav,job) saves and restores the contents of
c                             the internal common blocks used by
c                             lsode (see part iii below).
c                             rsav must be a real array of length 218
c                             or more, and isav must be an integer
c                             array of length 41 or more.
c                             job=1 means save common into rsav/isav.
c                             job=2 means restore common from rsav/isav.
c                                srcom is useful if one is
c                             interrupting a run and restarting
c                             later, or alternating between two or
c                             more problems solved with lsode.
c
c   call intdy(,,,,,)         provide derivatives of y, of various
c        (see below)          orders, at a specified point t, if
c                             desired.  it may be called only after
c                             a successful return from lsode.
c
c the detailed instructions for using intdy are as follows.
c the form of the call is..
c
c   call intdy (t, k, rwork(21), nyh, dky, iflag)
c
c the input parameters are..
c
c t         = value of independent variable where answers are desired
c             (normally the same as the t last returned by lsode).
c             for valid results, t must lie between tcur - hu and tcur.
c             (see optional outputs for tcur and hu.)
c k         = integer order of the derivative desired.  k must satisfy
c             0 .le. k .le. nqcur, where nqcur is the current order
c             (see optional outputs).  the capability corresponding
c             to k = 0, i.e. computing y(t), is already provided
c             by lsode directly.  since nqcur .ge. 1, the first
c             derivative dy/dt is always available with intdy.
c rwork(21) = the base address of the history array yh.
c nyh       = column length of yh, equal to the initial value of neq.
c
c the output parameters are..
c
c dky       = a real array of length neq containing the computed value
c             of the k-th derivative of y(t).
c iflag     = integer flag, returned as 0 if k and t were legal,
c             -1 if k was illegal, and -2 if t was illegal.
c             on an error return, a message is also written.
c-----------------------------------------------------------------------
c part iii.  common blocks.
c
c if lsode is to be used in an overlay situation, the user
c must declare, in the primary overlay, the variables in..
c   (1) the call sequence to lsode,
c   (2) the two internal common blocks
c         /ls_lsode/  of length  257  (218 double precision words
c                         followed by 39 integer words),
c         /eh_lsode/  of length  2 (integer words).
c
c if lsode is used on a system in which the contents of internal
c common blocks are not preserved between calls, the user should
c declare the above two common blocks in his main program to insure
c that their contents are preserved.
c
c if the solution of a given problem by lsode is to be interrupted
c and then later continued, such as when restarting an interrupted run
c or alternating between two or more problems, the user should save,
c following the return from the last lsode call prior to the
c interruption, the contents of the call sequence variables and the
c internal common blocks, and later restore these values before the
c next lsode call for that problem.  to save and restore the common
c blocks, use subroutine srcom (see part ii above).
c
c-----------------------------------------------------------------------
c part iv.  optionally replaceable solver routines.
c
c below are descriptions of two routines in the lsode package which
c relate to the measurement of errors.  either routine can be
c replaced by a user-supplied version, if desired.  however, since such
c a replacement may have a major impact on performance, it should be
c done only when absolutely necessary, and only with great caution.
c (note.. the means by which the package version of a routine is
c superseded by the user-s version may be system-dependent.)
c
c (a) ewset.
c the following subroutine is called just before each internal
c integration step, and sets the array of error weights, ewt, as
c described under itol/rtol/atol above..
c     subroutine ewset (neq, itol, rtol, atol, ycur, ewt)
c where neq, itol, rtol, and atol are as in the lsode call sequence,
c ycur contains the current dependent variable vector, and
c ewt is the array of weights set by ewset.
c
c if the user supplies this subroutine, it must return in ewt(i)
c (i = 1,...,neq) a positive quantity suitable for comparing errors
c in y(i) to.  the ewt array returned by ewset is passed to the
c vnorm_lsode routine (see below), and also used by lsode in the computation
c of the optional output imxer, the diagonal jacobian approximation,
c and the increments for difference quotient jacobians.
c
c in the user-supplied version of ewset, it may be desirable to use
c the current values of derivatives of y.  derivatives up to order nq
c are available from the history array yh, described above under
c optional outputs.  in ewset, yh is identical to the ycur array,
c extended to nq + 1 columns with a column length of nyh and scale
c factors of h**j/factorial(j).  on the first call for the problem,
c given by nst = 0, nq is 1 and h is temporarily set to 1.0.
c the quantities nq, nyh, h, and nst can be obtained by including
c in ewset the statements..
c     double precision h, rls
c     common /ls_lsode/ rls(218),ils(39)
c     nq = ils(35)
c     nyh = ils(14)
c     nst = ils(36)
c     h = rls(212)
c thus, for example, the current value of dy/dt can be obtained as
c ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
c unnecessary when nst = 0).
c
c (b) vnorm_lsode.
c the following is a real function routine which computes the weighted
c root-mean-square norm of a vector v..
c     d = vnorm_lsode(n, v, w)
c where..
c   n = the length of the vector,
c   v = real array of length n containing the vector,
c   w = real array of length n containing weights,
c   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
c vnorm_lsode is called with n = neq and with w(i) = 1.0/ewt(i), where
c ewt is as set by subroutine ewset.
c
c if the user supplies this function, it should return a non-negative
c value of vnorm_lsode suitable for use in the error control in lsode.
c none of the arguments should be altered by vnorm_lsode.
c for example, a user-supplied vnorm_lsode routine might..
c   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
c   -ignore some components of v in the norm, with the effect of
c    suppressing the error control on those components of y.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c other routines in the lsode package.
c
c in addition to subroutine lsode, the lsode package includes the
c following subroutines and function routines..
c  intdy    computes an interpolated value of the y vector at t = tout.
c  stode    is the core integrator, which does one step of the
c           integration and the associated error control.
c  cfode    sets all method coefficients and test constants.
c  prepj    computes and preprocesses the jacobian matrix j = df/dy
c           and the newton iteration matrix p = i - h*l0*j.
c  solsy    manages solution of linear system in chord iteration.
c  ewset    sets the error weight vector ewt before each step.
c  vnorm_lsode    computes the weighted r.m.s. norm of a vector.
c  srcom    is a user-callable routine to save and restore
c           the contents of the internal common blocks.
c  dgefa and dgesl   are routines from linpack for solving full
c           systems of linear algebraic equations.
c  dgbfa and dgbsl   are routines from linpack for solving banded
c           linear systems.
c  daxpy, dscal, idamax, and ddot   are basic linear algebra modules
c           (blas) used by the above linpack routines.
c  d1mach   computes the unit roundoff in a machine-independent manner.
c  xerrwv, xsetun, and xsetf   handle the printing of all error
c           messages and warnings.  xerrwv is machine-dependent.
c note..  vnorm_lsode, idamax, ddot, and d1mach are function routines.
c all the others are subroutines.
c
c the intrinsic and external routines used by lsode are..
c dabs, dmax1, dmin1, dfloat, max0, min0, mod, dsign, dsqrt, and write.
c
c a block data subprogram is also included with the package,
c for loading some of the variables in internal common.
c
c-----------------------------------------------------------------------
c the following card is for optimized compilation on llnl compilers.
clll. optimize
c-----------------------------------------------------------------------
      external prepj, solsy
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, iflag, imxer, kgo, lf0,
     1   leniw, lenrw, lenwm, ml, mord, mu, mxhnl0, mxstp0
      real(c_double) rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real(c_double) atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0, vnorm_lsode
CXXX     2   d1mach, vnorm_lsode
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following internal common block contains
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of the block is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsode first,
c then those local to subroutine stode, and finally those used
c for communication.  the block is declared in subroutines
c lsode, intdy, stode, prepj, and solsy.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls_lsode/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      meth = mf/10
      miter = mf - 10*meth
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 5) go to 608
      if (miter .le. 3) go to 30
      ml = iwork(1)
      mu = iwork(2)
      if (ml .lt. 0 .or. ml .ge. n) go to 609
      if (mu .lt. 0 .or. mu .ge. n) go to 610
 30   continue
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0d0
      hmxi = 0.0d0
      hmin = 0.0d0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0d0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0d0) go to 615
      hmxi = 0.0d0
      if (hmax .gt. 0.0d0) hmxi = 1.0d0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0d0) go to 616
c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
c-----------------------------------------------------------------------
 60   lyh = 21
      if (istate .eq. 1) nyh = n
      lwm = lyh + (maxord + 1)*nyh
      if (miter .eq. 0) lenwm = 0
      if (miter .eq. 1 .or. miter .eq. 2) lenwm = n*n + 2
      if (miter .eq. 3) lenwm = n + 2
      if (miter .ge. 4) lenwm = (2*ml + mu + 1)*n + 2
      lewt = lwm + lenwm
      lsavf = lewt + n
      lacor = lsavf + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      if (miter .eq. 0 .or. miter .eq. 3) leniw = 20
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      do 70 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0d0) go to 619
        if (atoli .lt. 0.0d0) go to 620
 70     continue
      if (istate .eq. 1) go to 100
c if istate = 3, set flag to signal parameter changes to stode. --------
      jstart = -1
      if (nq .le. maxord) go to 90
c maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------
      do 80 i = 1,n
 80     rwork(i+lsavf-1) = rwork(i+lwm-1)
c reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
 90   if (miter .gt. 0) rwork(lwm) = dsqrt(uround)
      if (n .eq. nyh) go to 200
c neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0d0
      go to 200
c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  uround = d1mach(4)
      tn = t
      if (itask .ne. 4 .and. itask .ne. 5) go to 110
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0d0) go to 625
      if (h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0)
     1   h0 = tcrit - t
 110  jstart = 0
      if (miter .gt. 0) rwork(lwm) = dsqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0d0
      nqu = 0
      ccmax = 0.3d0
      maxcor = 3
      msbp = 20
      mxncf = 10
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0))
      nfe = 1
c load the initial value vector in yh. ---------------------------------
      do 115 i = 1,n
 115    rwork(i+lyh-1) = y(i)
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0d0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 120 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 621
 120    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                       1
c where   w0     = max ( abs(t), abs(tout) ),
c         f(i)   = i-th component of initial value of f,
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      if (h0 .ne. 0.0d0) go to 180
      tdist = dabs(tout - t)
      w0 = dmax1(dabs(t),dabs(tout))
      if (tdist .lt. 2.0d0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = dmax1(tol,rtol(i))
 140  if (tol .gt. 0.0d0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = dabs(y(i))
        if (ayi .ne. 0.0d0) tol = dmax1(tol,atoli/ayi)
 150    continue
 160  tol = dmax1(tol,100.0d0*uround)
      tol = dmin1(tol,0.001d0)
      sum = vnorm_lsode(n, rwork(lf0), rwork(lewt))
      sum = 1.0d0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0d0/dsqrt(sum)
      h0 = dmin1(h0,tdist)
      h0 = dsign(h0,tout-t)
c adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = dabs(h0)*hmxi
      if (rh .gt. 1.0d0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270
c-----------------------------------------------------------------------
c block d.
c the next code block is for continuation calls only (istate = 2 or 3)
c and is to check stop conditions before taking a step.
c-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0d0 + 100.0d0*uround)
      if ((tp - tout)*h .gt. 0.0d0) go to 623
      if ((tn - tout)*h .lt. 0.0d0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
      if ((tcrit - tout)*h .lt. 0.0d0) go to 625
      if ((tn - tout)*h .lt. 0.0d0) go to 245
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
 245  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      if (istate .eq. 2) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stode.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 510
 260    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm_lsode(n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0d0) go to 280
      tolsf = tolsf*2.0d0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xerrwv(50hlsode--  warning..internal t (=r1) and h (=r2) are,
     1   50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h      such that in the machine, t + h = t on the next step  ,
     1   60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      (h = step size). solver will continue anyway,
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv(50hlsode--  above warning has been issued i1 times.  ,
     1   50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      it will not be issued again for this problem,
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
 290  continue
c-----------------------------------------------------------------------
c     call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prepj,solsy)
c-----------------------------------------------------------------------
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),
     2   f, jac, prepj, solsy)
      kgo = 1 - kflag
      go to (300, 530, 540), kgo
c-----------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).  test for stop conditions.
c-----------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 310  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
c itask = 3.  jump to exit if tout was reached. ------------------------
 330  if ((tn - tout)*h .ge. 0.0d0) go to 400
      go to 250
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  if ((tn - tout)*h .lt. 0.0d0) go to 345
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
c-----------------------------------------------------------------------
c block g.
c the following block handles all successful returns from lsode.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.
c if istate = 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c-----------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xerrwv(
     1  60hlsode--  repeated calls with istate = 1 and tout = t (=r1)  ,
     1   60, 301, 0, 0, 0, 0, 1, t, 0.0d0)
      go to 800
c-----------------------------------------------------------------------
c block h.
c the following block handles all unsuccessful returns other than
c those for illegal input.  first the error message routine is called.
c if there was an error test or convergence test failure, imxer is set.
c then y is loaded from yh, t is set to tn, and the illegal input
c counter illin is set to 0.  the optional outputs are loaded into
c the work arrays before returning.
c-----------------------------------------------------------------------
c the maximum number of steps was taken before reaching tout. ----------
 500  call xerrwv(50hlsode--  at current t (=r1), mxstep (=i1) steps   ,
     1   50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      taken on this call before reaching tout     ,
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xerrwv(50hlsode--  at t (=r1), ewt(i1) has become r2 .le. 0.,
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  call xerrwv(50hlsode--  at t (=r1), too much accuracy requested  ,
     1   50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      for precision of machine..  see tolsf (=r2) ,
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xerrwv(50hlsode--  at t(=r1) and step size h(=r2), the error,
     1   50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      test failed repeatedly or with abs(h) = hmin,
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xerrwv(50hlsode--  at t (=r1) and step size h (=r2), the    ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(50h      corrector convergence failed repeatedly     ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(30h      or with abs(h) = hmin   ,
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
c compute imxer if relevant. -------------------------------------------
 560  big = 0.0d0
      imxer = 1
      do 570 i = 1,n
        size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer
c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xerrwv(30hlsode--  istate (=i1) illegal ,
     1   30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
      go to 700
 602  call xerrwv(30hlsode--  itask (=i1) illegal  ,
     1   30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
      go to 700
 603  call xerrwv(50hlsode--  istate .gt. 1 but lsode not initialized  ,
     1   50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      go to 700
 604  call xerrwv(30hlsode--  neq (=i1) .lt. 1     ,
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
      go to 700
 605  call xerrwv(50hlsode--  istate = 3 and neq increased (i1 to i2)  ,
     1   50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 606  call xerrwv(30hlsode--  itol (=i1) illegal   ,
     1   30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
      go to 700
 607  call xerrwv(30hlsode--  iopt (=i1) illegal   ,
     1   30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 608  call xerrwv(30hlsode--  mf (=i1) illegal     ,
     1   30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
      go to 700
 609  call xerrwv(50hlsode--  ml (=i1) illegal.. .lt.0 or .ge.neq (=i2),
     1   50, 9, 0, 2, ml, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 610  call xerrwv(50hlsode--  mu (=i1) illegal.. .lt.0 or .ge.neq (=i2),
     1   50, 10, 0, 2, mu, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 611  call xerrwv(30hlsode--  maxord (=i1) .lt. 0  ,
     1   30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)
      go to 700
 612  call xerrwv(30hlsode--  mxstep (=i1) .lt. 0  ,
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
      go to 700
 613  call xerrwv(30hlsode--  mxhnil (=i1) .lt. 0  ,
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
      go to 700
 614  call xerrwv(40hlsode--  tout (=r1) behind t (=r2)      ,
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv(50h      integration direction is given by h0 (=r1)  ,
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
      go to 700
 615  call xerrwv(30hlsode--  hmax (=r1) .lt. 0.0  ,
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
      go to 700
 616  call xerrwv(30hlsode--  hmin (=r1) .lt. 0.0  ,
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
      go to 700
 617  call xerrwv(
     1  60hlsode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2),
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 618  call xerrwv(
     1  60hlsode--  iwork length needed, leniw (=i1), exceeds liw (=i2),
     1   60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
      go to 700
 619  call xerrwv(40hlsode--  rtol(i1) is r1 .lt. 0.0        ,
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
      go to 700
 620  call xerrwv(40hlsode--  atol(i1) is r1 .lt. 0.0        ,
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv(40hlsode--  ewt(i1) is r1 .le. 0.0         ,
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
      go to 700
 622  call xerrwv(
     1  60hlsode--  tout (=r1) too close to t(=r2) to start integration,
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(
     1  60hlsode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  ,
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(
     1  60hlsode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ,
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(
     1  60hlsode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ,
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv(50hlsode--  at start of problem, too much accuracy   ,
     1   50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  60h      requested for precision of machine..  see tolsf (=r1) ,
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv(50hlsode--  trouble from intdy. itask = i1, tout = r1,
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xerrwv(50hlsode--  repeated occurrences of illegal input    ,
     1   50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
c
 800  call xerrwv(50hlsode--  run aborted.. apparent infinite loop     ,
     1   50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
      return
c----------------------- end of subroutine lsode -----------------------
      end
      block data for_lsode
c-----------------------------------------------------------------------
c this data subprogram loads variables into the internal common
c blocks used by the odepack solvers.  the variables are
c defined as follows..
c   illin  = counter for the number of consecutive times the package
c            was called with illegal input.  the run is stopped when
c            illin reaches 5.
c   ntrep  = counter for the number of consecutive times the package
c            was called with istate = 1 and tout = t.  the run is
c            stopped when ntrep reaches 5.
c   mesflg = flag to control printing of error messages.  1 means print,
c            0 means no printing.
c   lunit  = default value of logical unit number for printing of error
c            messages.
c-----------------------------------------------------------------------
      integer illin, iduma, ntrep, idumb, iowns, icomm, mesflg, lunit
      real*8 rowns, rcomm
      common /ls_lsode/ rowns(209), rcomm(9),
     1   illin, iduma(10), ntrep, idumb(2), iowns(6), icomm(19)
      common /eh_lsode/ mesflg, lunit
      data illin/0/, ntrep/0/
      data mesflg/1/, lunit/6/
c
c----------------------- end of block data -----------------------------
      end block data

      subroutine cfode (meth, elco, tesco)
      use iso_c_binding, only : c_double
clll. optimize
      integer meth
      integer i, ib, nq, nqm1, nqp1
      real(c_double) elco, tesco
      real(c_double)  agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
c-----------------------------------------------------------------------
c cfode is called by the integrator routine to set coefficients
c needed there.  the coefficients for the current method, as
c given by the value of meth, are set for all orders and saved.
c the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
c (a smaller value of the maximum order is also allowed.)
c cfode is called once at the beginning of the problem,
c and is not called again unless and until meth is changed.
c
c the elco array contains the basic method coefficients.
c the coefficients el(i), 1 .le. i .le. nq+1, for the method of
c order nq are stored in elco(i,nq).  they are given by a genetrating
c polynomial, i.e.,
c     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
c for the implicit adams methods, l(x) is given by
c     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
c for the bdf methods, l(x) is given by
c     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
c where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c the tesco array contains test constants used for the
c local error test and the selection of step size and/or order.
c at order nq, tesco(k,nq) is used for the selection of step
c size at order nq - 1 if k = 1, at order nq if k = 2, and at order
c nq + 1 if k = 3.
c-----------------------------------------------------------------------
      dimension pc(12)
c
      go to (100, 200), meth
c
 100  elco(1,1) = 1.0d0
      elco(2,1) = 1.0d0
      tesco(1,1) = 0.0d0
      tesco(2,1) = 2.0d0
      tesco(1,2) = 1.0d0
      tesco(3,12) = 0.0d0
      pc(1) = 1.0d0
      rqfac = 1.0d0
      do 140 nq = 2,12
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/dfloat(nq)
        nqm1 = nq - 1
        fnqm1 = dfloat(nqm1)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0d0
        do 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0d0
        tsign = 1.0d0
        do 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/dfloat(i)
 120      xpin = xpin + tsign*pc(i)/dfloat(i+1)
c store coefficients in elco and tesco. --------------------------------
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0d0
        do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/dfloat(i)
        agamq = rqfac*xpin
        ragq = 1.0d0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/dfloat(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
      return
c
 200  pc(1) = 1.0d0
      rq1fac = 1.0d0
      do 230 nq = 1,5
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        fnq = dfloat(nq)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0d0
        do 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
        do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0d0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = dfloat(nqp1)/elco(1,nq)
        tesco(3,nq) = dfloat(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
      return
c----------------------- end of subroutine cfode -----------------------
      end
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
      use iso_c_binding, only : c_double

clll. optimize
c-----------------------------------------------------------------------
c this subroutine sets the error weight vector ewt according to
c     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
c with the subscript on rtol and/or atol possibly replaced by 1 above,
c depending on the value of itol.
c-----------------------------------------------------------------------
      integer n, itol
      integer i
      real(c_double) rtol, atol, ycur, ewt
      dimension rtol(*), atol(*), ycur(n), ewt(n)
c
      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(i)
      return
c----------------------- end of subroutine ewset -----------------------
      end
      subroutine intdy (t, k, yh, nyh, dky, iflag)
      use iso_c_binding, only : c_double
clll. optimize
      integer k, nyh, iflag
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real(c_double) t, yh, dky
      real(c_double) rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real(c_double) c, r, s, tp
      dimension yh(nyh,*), dky(*)
      common /ls_lsode/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c intdy computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.  this routine
c is called within the package with k = 0 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in the usage documentation.)
c-----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c-----------------------------------------------------------------------
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0d0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0d0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = dfloat(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = dfloat(ic)
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   call xerrwv(30hintdy--  k (=i1) illegal      ,
     1   30, 51, 0, 1, k, 0, 0, 0.0d0, 0.0d0)
      iflag = -1
      return
 90   call xerrwv(30hintdy--  t (=r1) illegal      ,
     1   30, 52, 0, 0, 0, 0, 1, t, 0.0d0)
      call xerrwv(
     1  60h      t not in interval tcur - hu (= r1) to tcur (=r2)      ,
     1   60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
c----------------------- end of subroutine intdy -----------------------
      end
      subroutine prepj (neq, y, yh, nyh, ewt, ftem, savf, wm, iwm,
     1   f, jac)
      use iso_c_binding, only : c_double
clll. optimize
      external f, jac
      integer neq, nyh, iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, ier, ii, j, j1, jj, lenp,
     1   mba, mband, meb1, meband, ml, ml3, mu, np1
      real(c_double) y, yh, ewt, ftem, savf, wm
      real(c_double) rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real(c_double) con, di, fac, hl0, r, r0, srur, yi, yj, yjj,
     1   vnorm_lsode
      dimension neq(*), y(*), yh(nyh,*), ewt(*), ftem(*), savf(*),
     1   wm(*), iwm(*)
      common /ls_lsode/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c prepj is called by stode to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c here j is computed by the user-supplied routine jac if
c miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
c if miter = 3, a diagonal approximation to j is used.
c j is stored in wm and replaced by p.  if miter .ne. 3, p is then
c subjected to lu decomposition in preparation for later solution
c of linear systems with p as coefficient matrix. this is done
c by dgefa if miter = 1 or 2, and by dgbfa if miter = 4 or 5.
c
c in addition to variables described previously, communication
c with prepj uses the following..
c y     = array containing predicted values on entry.
c ftem  = work array of length n (acor in stode).
c savf  = array containing f evaluated at predicted y.
c wm    = real work space for matrices.  on output it contains the
c         inverse diagonal matrix if miter = 3 and the lu decomposition
c         of p if miter is 1, 2 , 4, or 5.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround), used in numerical jacobian increments.
c         wm(2) = h*el0, saved for later use if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c el0   = el(1) (input).
c ierpj = output error flag,  = 0 if no trouble, .gt. 0 if
c         p matrix found to be singular.
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         (or approximation) is now current.
c this routine also uses the common variables el0, h, tn, uround,
c miter, n, nfe, and nje.
c-----------------------------------------------------------------------
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
      go to (100, 200, 300, 400, 500), miter
c if miter = 1, call jac and multiply by scalar. -----------------------
 100  lenp = n*n
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0d0
      call jac (neq, tn, y, 0, 0, wm(3), n)
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 240
c if miter = 2, make n calls to f to approximate j. --------------------
 200  fac = vnorm_lsode(n, savf, ewt)
      r0 = 1000.0d0*dabs(h)*uround*dfloat(n)*fac
      if (r0 .eq. 0.0d0) r0 = 1.0d0
      srur = wm(1)
      j1 = 2
      do 230 j = 1,n
        yj = y(j)
        r = dmax1(srur*dabs(yj),r0/ewt(j))
        y(j) = y(j) + r
        fac = -hl0/r
        call f (neq, tn, y, ftem)
        do 220 i = 1,n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
        y(j) = yj
        j1 = j1 + n
 230    continue
      nfe = nfe + n
c add identity matrix. -------------------------------------------------
 240  j = 3
      np1 = n + 1
      do 250 i = 1,n
        wm(j) = wm(j) + 1.0d0
 250    j = j + np1
c do lu decomposition on p. --------------------------------------------
      call dgefa (wm(3), n, n, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  wm(2) = hl0
      r = el0*0.1d0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (neq, tn, y, wm(3))
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1d0*r0 - h*(wm(i+2) - savf(i))
        wm(i+2) = 1.0d0
        if (dabs(r0) .lt. uround/ewt(i)) go to 320
        if (dabs(di) .eq. 0.0d0) go to 330
        wm(i+2) = 0.1d0*r0/di
 320    continue
      return
 330  ierpj = 1
      return
c if miter = 4, call jac and multiply by scalar. -----------------------
 400  ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      do 410 i = 1,lenp
 410    wm(i+2) = 0.0d0
      call jac (neq, tn, y, ml, mu, wm(ml3), meband)
      con = -hl0
      do 420 i = 1,lenp
 420    wm(i+2) = wm(i+2)*con
      go to 570
c if miter = 5, make mband calls to f to approximate j. ----------------
 500  ml = iwm(1)
      mu = iwm(2)
      mband = ml + mu + 1
      mba = min0(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = vnorm_lsode(n, savf, ewt)
      r0 = 1000.0d0*dabs(h)*uround*dfloat(n)*fac
      if (r0 .eq. 0.0d0) r0 = 1.0d0
      do 560 j = 1,mba
        do 530 i = j,n,mband
          yi = y(i)
          r = dmax1(srur*dabs(yi),r0/ewt(i))
 530      y(i) = y(i) + r
        call f (neq, tn, y, ftem)
        do 550 jj = j,n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = dmax1(srur*dabs(yjj),r0/ewt(jj))
          fac = -hl0/r
          i1 = max0(jj-mu,1)
          i2 = min0(jj+ml,n)
          ii = jj*meb1 - ml + 2
          do 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      continue
 560    continue
      nfe = nfe + mba
c add identity matrix. -------------------------------------------------
 570  ii = mband + 2
      do 580 i = 1,n
        wm(ii) = wm(ii) + 1.0d0
 580    ii = ii + meband
c do lu decomposition of p. --------------------------------------------
      call dgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c----------------------- end of subroutine prepj -----------------------
      end
      subroutine solsy (wm, iwm, x, tem)
      use iso_c_binding, only : c_double
clll. optimize
      integer iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, meband, ml, mu
      real(c_double) wm, x, tem
      real(c_double) rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real(c_double) di, hl0, phl0, r
      dimension wm(*), iwm(*), x(*), tem(*)
      common /ls_lsode/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord iteration.  it is called if miter .ne. 0.
c if miter is 1 or 2, it calls dgesl to accomplish this.
c if miter = 3 it updates the coefficient h*el0 in the diagonal
c matrix, and then computes the solution.
c if miter is 4 or 5, it calls dgbsl.
c communication with solsy uses the following variables..
c wm    = real work space containing the inverse diagonal matrix if
c         miter = 3 and the lu decomposition of the matrix otherwise.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround) (not used here),
c         wm(2) = hl0, the previous value of h*el0, used if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c tem   = vector of work space of length n, not used in this version.
c iersl = output flag (in common).  iersl = 0 if no trouble occurred.
c         iersl = 1 if a singular matrix arose with miter = 3.
c this routine also uses the common variables el0, h, miter, and n.
c-----------------------------------------------------------------------
      iersl = 0
      go to (100, 100, 300, 400, 400), miter
 100  call dgesl (wm(3), n, n, iwm(21), x, 0)
      return
c
 300  phl0 = wm(2)
      hl0 = h*el0
      wm(2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        di = 1.0d0 - r*(1.0d0 - 1.0d0/wm(i+2))
        if (dabs(di) .eq. 0.0d0) go to 390
 320    wm(i+2) = 1.0d0/di
 330  do 340 i = 1,n
 340    x(i) = wm(i+2)*x(i)
      return
 390  iersl = 1
      return
c
 400  ml = iwm(1)
      mu = iwm(2)
      meband = 2*ml + mu + 1
      call dgbsl (wm(3), meband, n, ml, mu, iwm(21), x, 0)
      return
c----------------------- end of subroutine solsy -----------------------
      end
      subroutine srcom (rsav, isav, job)
      use iso_c_binding, only : c_double
      save
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001 and eh0001, which are used internally
c by one or more odepack solvers.
c
c rsav = real array of length 218 or more.
c isav = integer array of length 41 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ieh, ils
      integer i, lenils, lenrls
      real(c_double) rsav,   rls
      dimension rsav(*), isav(*)
      common /ls_lsode/ rls(218), ils(39)
      common /eh_lsode/ ieh(2)
      data lenrls/218/, lenils/39/
c
      if (job .eq. 2) go to 100
c
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      isav(lenils+1) = ieh(1)
      isav(lenils+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      ieh(1) = isav(lenils+1)
      ieh(2) = isav(lenils+2)
      return
c----------------------- end of subroutine srcom -----------------------
      end
      subroutine stode (neq, y, yh, nyh, yh1, ewt, savf, acor,
     1   wm, iwm, f, jac, pjac, slvs)
      use iso_c_binding, only : c_double
clll. optimize
      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      real(c_double) y, yh, yh1, ewt, savf, acor, wm
      real(c_double) conit, crate, el, elco, hold, rmax, tesco,
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real(c_double) dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     1   r, rh, rhdn, rhsm, rhup, told, vnorm_lsode
      dimension neq(*), y(*), yh(nyh,*), yh1(*), ewt(*), savf(*),
     1   acor(*), wm(*), iwm(*)
      common /ls_lsode/ conit, crate, el(13), elco(13,12),
     1   hold, rmax, tesco(3,12),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),
     3   ialth, ipup, lmax, meo, nqnyh, nslp,
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c stode performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stode is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used, or the jacobian structure.
c communication with stode is done with the following variables..
c
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of length n containing multiplicative weights
c          for local error measurements.  local errors in y(i) are
c          compared to 1.0/ewt(i) in various error tests.
c savf   = an array of working storage, of length n.
c          also used for input of yh(*,maxord+2) when jstart = -1
c          and maxord .lt. the current order nq.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          and p = i - h*el0*jac, if a chord method is being used.
c slvs   = name of routine to solve linear system in chord iteration.
c ccmax  = maximum relative change in h*el0 before pjac is called.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h, maxord,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was succesful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c              -3  fatal error in pjac or slvs.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c maxcor = the maximum number of corrector iterations allowed.
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).
c mxncf  = maximum number of convergence failures allowed.
c meth/miter = the method flags.  see description in driver.
c n      = the number of first-order differential equations.
c-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0d0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0d0
      rc = 0.0d0
      el0 = 1.0d0
      crate = 0.7d0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c if an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c if the caller has changed meth, cfode is called to reset
c the coefficients of the method.
c if the caller has changed maxord to a value less than the current
c order nq, nq is reduced to maxord, and a new h chosen accordingly.
c if h is to be changed, yh must be rescaled.
c if h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfode (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dfloat(nq+2)
      ddn = vnorm_lsode(n, savf, ewt)/tesco(1,l)
      exdn = 1.0d0/dfloat(l)
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
      rh = dmin1(rhdn,1.0d0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = dmin1(rh,dabs(h/hold))
      h = hold
      go to 175
c-----------------------------------------------------------------------
c cfode is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
 140  call cfode (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dfloat(nq+2)
      go to (160, 170, 200), iret
c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = dmax1(rh,hmin/dabs(h))
 175  rh = dmin1(rh,rmax)
      rh = rh/dmax1(1.0d0,dabs(h)*hmxi*rh)
      r = 1.0d0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every msbp steps.
c-----------------------------------------------------------------------
 200  if (dabs(rc-1.0d0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0d0
      nslp = nst
      crate = 0.7d0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0d0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm_lsode(n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnorm_lsode(n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  if (m .ne. 0) crate = dmax1(0.2d0*crate,del/delp)
      dcon = del*dmin1(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0d0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
      delp = del
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or mxncf failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0d0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (dabs(h) .le. hmin*1.00001d0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25d0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm_lsode(n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0d0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0d0
      if (dabs(h) .le. hmin*1.00001d0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0d0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0d0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm_lsode(n, savf, ewt)/tesco(3,nq)
      exup = 1.0d0/dfloat(l+1)
      rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
 540  exsm = 1.0d0/dfloat(l)
      rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
      rhdn = 0.0d0
      if (nq .eq. 1) go to 560
      ddn = vnorm_lsode(n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0d0/dfloat(nq)
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1d0) go to 610
      r = el(l)/dfloat(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1d0)) go to 610
      if (kflag .le. -2) rh = dmin1(rh,0.2d0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1d0
      rh = dmax1(hmin/dabs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0d0
 700  r = 1.0d0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return
c----------------------- end of subroutine stode -----------------------
      end
      real(c_double) function vnorm_lsode(n, v, w)
      use iso_c_binding, only : c_double
clll. optimize
c-----------------------------------------------------------------------
c this function routine computes the weighted root-mean-square norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vnorm_lsode = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
c-----------------------------------------------------------------------
      integer n,   i
      real(c_double) v, w,   sum
      dimension v(n), w(n)
      sum = 0.0d0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      vnorm_lsode = dsqrt(sum/dfloat(n))
      return
c----------------------- end of function vnorm_lsode -------------------------
      end
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      use iso_c_binding, only : c_double
      integer msg, nmes, nerr, level, ni, i1, i2, nr,
     1   i, lun, lunit, mesflg, ncpw, nch, nwds
      real(c_double) r1, r2
      dimension msg(nmes)
c-----------------------------------------------------------------------
c subroutines xerrwv, xsetf, and xsetun, as given here, constitute
c a simplified version of the slatec error handling package.
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.
c this version is in real(c_double).
c
c all arguments are input arguments.
c
c msg    = the message (hollerith literal or integer array).
c nmes   = the length of msg (number of characters).
c nerr   = the error number (not used).
c level  = the error level..
c          0 or 1 means recoverable (control returns to caller).
c          2 means fatal (run is aborted--see note below).
c ni     = number of integers (0, 1, or 2) to be printed with message.
c i1,i2  = integers to be printed, depending on ni.
c nr     = number of reals (0, 1, or 2) to be printed with message.
c r1,r2  = reals to be printed, depending on nr.
c
c note..  this routine is machine-dependent and specialized for use
c in limited context, in the following ways..
c 1. the number of hollerith characters stored per word, denoted
c    by ncpw below, is a data-loaded constant.
c 2. the value of nmes is assumed to be at most 60.
c    (multi-line messages are generated by repeated calls.)
c 3. if level = 2, control passes to the statement   stop
c    to abort the run.  this statement may be machine-dependent.
c 4. r1 and r2 are assumed to be in real(c_double) and are printed
c    in d21.13 format.
c 5. the common block /eh_lsode/ below is data-loaded (a machine-
c    dependent feature) with default values.
c    this block is needed for proper retention of parameters used by
c    this routine which the user can reset by calling xsetf or xsetun.
c    the variables in this block are as follows..
c       mesflg = print control flag..
c                1 means print all messages (the default).
c                0 means no printing.
c       lunit  = logical unit number for messages.
c                the default is 6 (machine-dependent).
c-----------------------------------------------------------------------
c the following are instructions for installing this routine
c in different machine environments.
c
c to change the default output unit, change the data statement
c in the block data subprogram below.
c
c for a different number of characters per word, change the
c data statement setting ncpw below, and format 10.  alternatives for
c various computers are shown in comment cards.
c
c for a different run-abort command, change the statement following
c statement 100 at the end.
c-----------------------------------------------------------------------
      common /eh_lsode/ mesflg, lunit
c-----------------------------------------------------------------------
c the following data-loaded value of ncpw is valid for the cdc-6600
c and cdc-7600 computers.
c     data ncpw/10/
c the following is valid for the cray-1 computer.
c     data ncpw/8/
c the following is valid for the burroughs 6700 and 7800 computers.
c     data ncpw/6/
c the following is valid for the pdp-10 computer.
c     data ncpw/5/
c the following is valid for the vax computer with 4 bytes per integer,
c and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.
      data ncpw/4/
c the following is valid for the pdp-11, or vax with 2-byte integers.
c     data ncpw/2/
c-----------------------------------------------------------------------
      if (mesflg .eq. 0) go to 100
c get logical unit number. ---------------------------------------------
      lun = lunit
c get number of words in message. --------------------------------------
      nch = min0(nmes,60)
      nwds = nch/ncpw
      if (nch .ne. nwds*ncpw) nwds = nwds + 1
c write the message. ---------------------------------------------------
      write (lun, 10) (msg(i),i=1,nwds)
c-----------------------------------------------------------------------
c the following format statement is to have the form
c 10  format(1x,mmann)
c where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.
c the following is valid for ncpw = 10.
c 10  format(1x,6a10)
c the following is valid for ncpw = 8.
c 10  format(1x,8a8)
c the following is valid for ncpw = 6.
c 10  format(1x,10a6)
c the following is valid for ncpw = 5.
c 10  format(1x,12a5)
c the following is valid for ncpw = 4.
  10  format(1x,15a4)
c the following is valid for ncpw = 2.
c 10  format(1x,30a2)
c-----------------------------------------------------------------------
      if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,23hin above message,  i1 =,i10)
      if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
      if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,23hin above message,  r1 =,d21.13)
      if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,15hin above,  r1 =,d21.13,3x,4hr2 =,d21.13)
c abort the run if level = 2. ------------------------------------------
 100  if (level .ne. 2) return
      stop 'r8lsode:'
c----------------------- end of subroutine xerrwv ----------------------
      end
      subroutine xsetf (mflag)
      use iso_c_binding, only : c_double

c
c this routine resets the print control flag mflag.
c
      integer mflag, mesflg, lunit
      common /eh_lsode/ mesflg, lunit
c
      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return
c----------------------- end of subroutine xsetf -----------------------
      end
      subroutine xsetun (lun)
      use iso_c_binding, only : c_double

c
c this routine resets the logical unit number for messages.
c
      integer lun, mesflg, lunit
      common /eh_lsode/ mesflg, lunit
c
      if (lun .gt. 0) lunit = lun
      return
c----------------------- end of subroutine xsetun ----------------------
      end



c=====================================================================
c----------STUB-ROUTINES: Full versions need to be brought in from 
c----------LINPACK, if they are to be called.  [BH, 100622]

      subroutine dgesl(dum,idum1,idum2,idum3,ddum,idum4)
      use iso_c_binding, only : c_double
      real(c_double) dum,ddum
      dimension ddum(*)
      integer idum1,idum2,idum3,idum4
      write(*,*)'LINPACK subroutine dgesl not presently installed'
      STOP
      return
      end

      subroutine dgbsl(dum,idum1,idum2,idum3,idum4,idum5,
     +     ddum,idum6)
      use iso_c_binding, only : c_double
      real(c_double) dum,ddum
      dimension ddum(*)
      integer idum1,idum2,idum3,idum4,idum5,idum6
      write(*,*)'LINPACK subroutine dgbsl not presently installed'
      STOP
      return
      end

      subroutine dgefa(dum,idum1,idum2,idum3,idum4)
      use iso_c_binding, only : c_double

      real(c_double) dum
      integer idum1,idum2,idum3,idum4
      write(*,*)'LINPACK subroutine dgefa not presently installed'
      STOP
      return
      end

      subroutine dgbfa(dum,idum1,idum2,idum3,idum4,
     +     idum5,idum6)
      use iso_c_binding, only : c_double

      real(c_double) dum
      integer idum1,idum2,idum3,idum4,idum5,idum6
      write(*,*)'LINPACK subroutine dgbfa not presently installed'
      STOP
      return
      end
c=====================================================================
