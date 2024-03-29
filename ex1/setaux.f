c     ============================================
      subroutine setaux(mbc,mx,xlower,dx,maux,aux)
c     ============================================
c     
c     # set auxiliary arrays 
c     # variable coefficient acoustics
c     #  aux(i,1) = impedance Z in i'th cell
c     #  aux(i,2) = sound speed c in i'th cell
c     
c     # Piecewise constant medium with single interface at x=0
c     # Density and sound speed to left and right are set in setprob.f
c

      use amr_module, only : NEEDS_TO_BE_SET
      use adjoint_module, only: innerprod_index, adjoint_flagging
      implicit none

      integer, intent(in) :: mbc, mx, maux
      double precision, intent(in) :: xlower, dx
      double precision, intent(out) :: aux(maux, 1-mbc:mx+mbc)

      common /comaux/ Zl, cl, Zr, cr
      double precision Zl, cl, Zr, cr

      integer i
      double precision xcell

      if (adjoint_flagging) then
          do i=1-mbc,mx+mbc

             if (aux(1,i) .eq. NEEDS_TO_BE_SET) then
                aux(innerprod_index,i) = 0.d0
             endif

            enddo
      endif

      do i=1-mbc,mx+mbc

         xcell = xlower + (i-0.5d0)*dx
         if (xcell .lt. 0.0d0) then
            aux(1,i) = Zl
            aux(2,i) = cl
         else
            aux(1,i) = Zr
            aux(2,i) = cr
         endif
      enddo

      return
      end
