c     
c     
c=========================================================
      subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c=========================================================
c     
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c     
c     
      implicit none

      integer, intent(in) :: meqn, mbc, mx, maux
      double precision, intent(in) :: xlower, dx, aux
      double precision, intent(out) :: q
      dimension q(meqn, 1-mbc:mx+mbc)
      dimension aux(maux, 1-mbc:mx+mbc)

      common /cqinit/ betal, betar, freql, freqr
      double precision betal, betar, freql, freqr

      integer i
      double precision xcell, al,ar,x1,x2,beta
c     
c     
      al = 1.d0
      ar = 0.d0
      x1 = -10.d0
      x2 = -1.d0
      beta = 5.d0

      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx


c     # wave packet with left-going and right-going wave packets:
c        q(1,i) = ar*dexp(-betar*(xcell-5.0d0)**2) * dsin(freqr*xcell)
c    .            + al*dexp(-betal*(xcell+5.0d0)**2) * dsin(freql*xcell)
c        q(2,i) = -ar*dexp(-betar*(xcell-5.0d0)**2) * dsin(freqr*xcell)
c    .            + al*dexp(-betal*(xcell+5.0d0)**2) * dsin(freql*xcell)

         q(1,i) = al/(1.d0 + exp(-beta*(xcell-x1)) + 
     .                       exp(beta*(xcell-x2))) * 
     .            sin(0.25d0*abs(xcell)**2.5d0)
         q(2,i) = q(1,i)

 150  continue
c     
      return
      end
