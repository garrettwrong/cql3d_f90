module ntdstore_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use restvty_mod, only : restvty

  !---END USE

!
!

contains

      subroutine ntdstore
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!

      if (n .ne. nchec*(n/nchec)) return
      if (nch(l_) .eq. nonch) nch(l_) =0
      call restvty
      nch(l_)=nch(l_)+1
      if (n .eq. 0) nch(l_)=1
      ptime(nch(l_),l_)=timet
      do 30 k=1,ngen
        if (cqlpmod .ne. "enabled") then
          pdens(nch(l_),k,l_)=xlndn(k,lr_)/zmaxpsi(lr_) ! <n>_FSA ZOW only

          pdenm(nch(l_),k,l_)=reden(k,lr_)   !  n_midplane
          pengy(nch(l_),k,l_)=energy(k,lr_)   ! <energy>_FSA
          pengym(nch(l_),k,l_)=energym(k,lr_) ! energy at midplane
          if (efswtchn.eq."neo_hh".and.k.eq.kelecg) then
             pcurr(nch(l_),k,l_)=currpar(lr_)
          else
             pcurr(nch(l_),k,l_)=curr(k,lr_)/3.e9 !A/cm2, <jpar>_FSA
          endif
          pcurrm(nch(l_),k,l_)=currm(k,lr_)/3.e9 !A/cm2, jpar at midplane
          pefld(nch(l_),l_)=elecfld(lr_)
        else
          pdens(nch(l_),k,l_)=denpar(k,ls_)
          pengy(nch(l_),k,l_)=enrgypa(k,ls_)
          pcurr(nch(l_),k,l_)=currm(k,l_)/3.e9
          pefld(nch(l_),l_)=elecfld(lr_)
        endif
        do 20 is=-1,15
          pentr(nch(l_),k,is,l_)=entr(k,is,l_)
 20     continue
 30   continue

      if (cqlpmod .ne. "enabled") then
         k=1
         pdenra(nch(l_),l_)=denra(k,lr_)
         pcurra(nch(l_),l_)=curra(k,lr_)/3.e9
         pfdenra(nch(l_),l_)=fdenra(k,lr_)
         pfcurra(nch(l_),l_)=fcurra(k,lr_)
         pucrit(nch(l_),l_)=ucrit(k,lr_)
         peoe0(nch(l_),l_)=eoe0(k,lr_)
         psrc(nch(l_),l_)=srckotot(lr_)
         peoed(nch(l_),l_)=elecfld(lr_)/elecr(lr_)
      endif

      consnp(nch(l_),l_)=consn(l_)
      constp(nch(l_),l_)=conserv
 60   continue
      if (kelecg.eq.0) go to 70
      if (n.lt.nonel .or. n.gt.noffel) go to 70
      if (l_ .eq. lmdpln_) then
        restp(nch(l_),lr_)=resist
        restnp(nch(l_),lr_)=resistn
        sptzrp(nch(l_),l_)=sptzr(l_)
        if (cqlpmod .ne. "enabled") rovsp(nch(l_),lr_)=rovs(lr_)
        vpov(nch(l_),lr_)=vpovth
      endif
      if (cqlpmod .eq. "enabled") rovsp(nch(l_),ls_)=rovsloc(l_)
 70   continue

      return
      end subroutine ntdstore
      
end module ntdstore_mod
