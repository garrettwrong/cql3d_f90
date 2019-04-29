
!..................................................................
!
!   Include file for dskin.f, 
!     which reads diskf distribution function file from CQL3D.
!
!   parameters iya,jxa,lrza,ngena need to be set 
!     in accord with iy,jx,lrz,ngen in diskf file.
!
!   Note: You might need to know that storage is slightly 
!         different here than in CQL3D.
!
!..................................................................



!      parameter(iya=80)
      parameter(iya=120)
!      parameter(jxa=80)
      parameter(jxa=120)
!      parameter(lrza=15)
      parameter(lrza=20)
      parameter(ngena=1)


      real*8 x(jxa),y(iya,lrza),rovera(lrza),elecfld(lrza), &
            bthr(lrza),btoru(lrza),bthr0(lrza),btor0(lrza), &
            reden(lrza,ngena),temp(lrza,ngena)
      real*8 radmin,vnorm,vmaxdvt,eovedd
      real*8 bnumb(ngena),fmass(ngena)
      integer*4 iy(lrza),itl(lrza),itu(lrza)
      real*8 f(iya,jxa,lrza,ngena)

    ! should be in comm.f90 now
      !COMMON /dskincomm/x,y,rovera,elecfld,bthr,btoru,bthr0,btor0
      !COMMON /dskincomm/reden,temp,radmin,vnorm,vmaxdvt,eovedd
      !COMMON /dskincomm/bnumb,fmass,itl,itu
      !COMMON /dskincomm/f

!      IDL related test of pointers:
!      common/dptr/fptr,roveraptr,radminptr
!      pointer(fptr,f_idl(1:iya,1:jxa,lrza,ngena))
!      pointer(roveraptr,rovera_idl(1:lrza))
!      pointer(radminptr,radmin_idl)
