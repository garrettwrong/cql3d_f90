c
      subroutine bavdens(k)
      use param_mod
      use cqcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..........................................................
c     This routine computes the bounce average of certain
c     density related quantities.
c.............................................................


      do 1 i=1,iy
        bavdn(i,lr_)=reden(k,lr_)
        bavpd(i,lr_)=batot(i,lr_)*sinn(i,lmdpln_)*reden(k,lr_)
 1    continue

c..................................................................
c     Exit if locquas.ne. "enabled"
c..................................................................

      if (kelecg.gt.0.or.locquas.eq."disabled".or.kelecm.ne.k) return
      call bcast(temc3,zero,iy)
      call bcast(temc4,zero,iy)

      
      lrange=lz
      if (numclas .eq. 1) lrange=lz/2+1

      
      do 30 l=1,lrange
        do 31 i=1,imax(l,lr_)        

          ! 1. All orbits that reach/pass a given poloidal node l:
          ! passing (i<itl), or trapped that could reach/pass this l;
          ! also includes Last Trapped, with tip at l=lz_bmax
          !(for such particle, lmax(itl)=lz_bmax; see micxinil)
            ax=dtau(i,l,lr_)/tau(i,lr_)
            y1=densz(l,ngen+1,negyrg,lr_)
            temc3(i)=temc3(i)+y1*ax
            temc4(i)=temc4(i)+(y1/bbpsi(l,lr_))*ax
                    
          ! 2a. Trapped, with tips between l and l+1 (ABOVE midplane):
          if (l.eq.lmax(i,lr_) .and. l.lt.lz_bmax(lr_)) then 
            ! Add contribution from orbit tip:
            ax=dtau(i,l+1,lr_)/tau(i,lr_)
            y1=densz(l,  ngen+1,negyrg,lr_)
            y2=densz(l+1,ngen+1,negyrg,lr_)
            xx=zboun(i,lr_)
            qq=exlin(y1,y2,z(l,lr_),z(l+1,lr_),xx)
            temc3(i)=temc3(i)+qq*ax
            temc4(i)=temc4(i)+(qq/psif(zboun(i,lr_)))*ax  
          endif
                  
          ! 2b. Trapped, with tips between l and l-1 (BELOW midplane):
          if (l.eq.lmax(i+iyh,lr_) .and. l.gt.lz_bmax(lr_)) then 
            ! Add contribution from orbit tip:
            ax=dtau(i,l-1,lr_)/tau(i,lr_) 
            y1=densz(l,  ngen+1,negyrg,lr_)
            y2=densz(l-1,ngen+1,negyrg,lr_)
            xx=zboun(i,lr_)
            qq=exlin(y1,y2,z(l,lr_),z(l-1,lr_),xx)
            temc3(i)=temc3(i)+qq*ax
            temc4(i)=temc4(i)+(qq/psif(zboun(i,lr_)))*ax
          endif
          
 31     continue
 30   continue
 
 
      ! Symmetrize around pitch=pi/2  
      do 40 i=1,iyh
        iii=iy+1-i
        temc3(iii)=temc3(i)
        temc4(iii)=temc4(i)
 40   continue
      call dcopy(iy,temc3,1,bavdn,1)
      
      
      ! Symmetrize around pitch=pi/2  
      do 50 i=2,iyh
        iii=iy+1-i
        bavpd(i,lr_)=(temc4(i)/sinn(i,lmdpln_)**2
     1    -temc3(i))*tann(i,lmdpln_)**2*sinn(i,lmdpln_)
        bavpd(iii,lr_)=+bavpd(i,lr_)
 50   continue
      bavpd(1,lr_)=0.
      bavpd(iy,lr_)=0.
            
      
      return
      end
