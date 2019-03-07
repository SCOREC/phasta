      module e3_dc_func_m
c------------------------------------------------------------------------
c  move the calculation of yiA0DCyj and the DC factor from e3dc.f to this
c  module
c------------------------------------------------------------------------
        use propar_m, only: npro
        use conpar_m, only: nflow
        use intpt_m, only: ngauss, intp
        use number_def_m
        use genpar_m, only: ipord
        implicit none
c
c
      contains
        subroutine calc_e3_dc_factor(DC,   gAgyi, 
     &                               g1yi, g2yi, g3yi, A0,
     &                               raLS, rtLS, giju, A0DC,
     &                               epsM, iDC)
c.......................................................................
c  calculate yiA0DCyj and the DC factor(numerical viscousity)        
c.......................................................................
          implicit none
c
          real*8, dimension(npro,nflow), intent(in):: g1yi, g2yi, g3yi
          real*8, dimension(npro,5,5), intent(in):: A0
          real*8, dimension(npro), intent(in):: raLS, rtLS
          real*8, dimension(npro,6), intent(in):: giju
          real*8, dimension(npro,4), intent(in):: A0DC
          real*8 :: epsM
          integer :: iDC
c
          real*8, dimension(npro,ngauss), intent(out):: DC
          real*8, dimension(npro,15), intent(out):: gAgyi
c
          real*8, dimension(npro):: gnorm
          real*8, dimension(npro,15):: A0gyi
          real*8, dimension(npro,6):: yiA0DCyj
          real*8 :: fact
          
c ... -----------------------> initialize <----------------------------
c
            A0gyi    = zero
            gAgyi    = zero
            yiA0DCyj = zero
            DC       = zero
c.... ----------------------->  global gradient  <----------------------
c
c.... calculate (A0 y_,j) --> A0gyi
c
c  A0 Y_{,1}
c
            A0gyi( :,1) = A0(:,1,1)*g1yi(:,1)
     &                  + A0(:,1,2)*g1yi(:,2)
     &                  + A0(:,1,3)*g1yi(:,3)
     &                  + A0(:,1,4)*g1yi(:,4)
     &                  + A0(:,1,5)*g1yi(:,5)      
            A0gyi( :,2) = A0(:,2,1)*g1yi(:,1)
     &                  + A0(:,2,2)*g1yi(:,2)
     &                  + A0(:,2,3)*g1yi(:,3)
     &                  + A0(:,2,4)*g1yi(:,4)
     &                  + A0(:,2,5)*g1yi(:,5)
            A0gyi( :,3) = A0(:,3,1)*g1yi(:,1)
     &                  + A0(:,3,2)*g1yi(:,2)
     &                  + A0(:,3,3)*g1yi(:,3)
     &                  + A0(:,3,4)*g1yi(:,4)
     &                  + A0(:,3,5)*g1yi(:,5)
            A0gyi( :,4) = A0(:,4,1)*g1yi(:,1)
     &                  + A0(:,4,2)*g1yi(:,2)
     &                  + A0(:,4,3)*g1yi(:,3)
     &                  + A0(:,4,4)*g1yi(:,4)
     &                  + A0(:,4,5)*g1yi(:,5)
            A0gyi( :,5) = A0(:,5,1)*g1yi(:,1)
     &                  + A0(:,5,2)*g1yi(:,2)
     &                  + A0(:,5,3)*g1yi(:,3)
     &                  + A0(:,5,4)*g1yi(:,4)
     &                  + A0(:,5,5)*g1yi(:,5)
c
c...A0 Y_{,2}
c
            A0gyi( :,6) = A0(:,1,1)*g2yi(:,1)
     &                  + A0(:,1,2)*g2yi(:,2)
     &                  + A0(:,1,3)*g2yi(:,3)
     &                  + A0(:,1,4)*g2yi(:,4)
     &                  + A0(:,1,5)*g2yi(:,5)
            A0gyi( :,7) = A0(:,2,1)*g2yi(:,1)
     &                  + A0(:,2,2)*g2yi(:,2)
     &                  + A0(:,2,3)*g2yi(:,3)
     &                  + A0(:,2,4)*g2yi(:,4)
     &                  + A0(:,2,5)*g2yi(:,5)
            A0gyi( :,8) = A0(:,3,1)*g2yi(:,1)
     &                  + A0(:,3,2)*g2yi(:,2)
     &                  + A0(:,3,3)*g2yi(:,3)
     &                  + A0(:,3,4)*g2yi(:,4)
     &                  + A0(:,3,5)*g2yi(:,5)
            A0gyi( :,9) = A0(:,4,1)*g2yi(:,1)
     &                  + A0(:,4,2)*g2yi(:,2)
     &                  + A0(:,4,3)*g2yi(:,3)
     &                  + A0(:,4,4)*g2yi(:,4)
     &                  + A0(:,4,5)*g2yi(:,5)
            A0gyi(:,10) = A0(:,5,1)*g2yi(:,1)
     &                  + A0(:,5,2)*g2yi(:,2)
     &                  + A0(:,5,3)*g2yi(:,3)
     &                  + A0(:,5,4)*g2yi(:,4)
     &                  + A0(:,5,5)*g2yi(:,5)
c
c  A0 Y_{,3}
c
            A0gyi(:,11) = A0(:,1,1)*g3yi(:,1)
     &                  + A0(:,1,2)*g3yi(:,2)
     &                  + A0(:,1,3)*g3yi(:,3)
     &                  + A0(:,1,4)*g3yi(:,4)
     &                  + A0(:,1,5)*g3yi(:,5)
            A0gyi(:,12) = A0(:,2,1)*g3yi(:,1)
     &                  + A0(:,2,2)*g3yi(:,2)
     &                  + A0(:,2,3)*g3yi(:,3)
     &                  + A0(:,2,4)*g3yi(:,4)
     &                  + A0(:,2,5)*g3yi(:,5)
            A0gyi(:,13) = A0(:,3,1)*g3yi(:,1)
     &                  + A0(:,3,2)*g3yi(:,2)
     &                  + A0(:,3,3)*g3yi(:,3)
     &                  + A0(:,3,4)*g3yi(:,4)
     &                  + A0(:,3,5)*g3yi(:,5)
            A0gyi(:,14) = A0(:,4,1)*g3yi(:,1)
     &                  + A0(:,4,2)*g3yi(:,2)
     &                  + A0(:,4,3)*g3yi(:,3)
     &                  + A0(:,4,4)*g3yi(:,4)
     &                  + A0(:,4,5)*g3yi(:,5)
            A0gyi(:,15) = A0(:,5,1)*g3yi(:,1)
     &                  + A0(:,5,2)*g3yi(:,2)
     &                  + A0(:,5,3)*g3yi(:,3)
     &                  + A0(:,5,4)*g3yi(:,4)
     &                  + A0(:,5,5)*g3yi(:,5)
c
c.... calculate (giju A0 y_,j) --> gAgyi
c

            gAgyi( :,1) = giju(:,1)*A0gyi( :,1)
     &                  + giju(:,4)*A0gyi( :,6)
     &                  + giju(:,5)*A0gyi(:,11)
c
            gAgyi( :,2) = giju(:,1)*A0gyi( :,2)
     &                  + giju(:,4)*A0gyi( :,7)
     &                  + giju(:,5)*A0gyi(:,12)
c
            gAgyi( :,3) = giju(:,1)*A0gyi( :,3)
     &                + giju(:,4)*A0gyi( :,8)
     &                + giju(:,5)*A0gyi(:,13)
c
            gAgyi( :,4) = giju(:,1)*A0gyi( :,4)
     &                  + giju(:,4)*A0gyi( :,9)
     &                  + giju(:,5)*A0gyi(:,14)
c
            gAgyi( :,5) = giju(:,1)*A0gyi( :,5)
     &                  + giju(:,4)*A0gyi(:,10)
     &                  + giju(:,5)*A0gyi(:,15)
c
            gAgyi( :,6) = giju(:,4)*A0gyi( :,1)
     &                  + giju(:,2)*A0gyi( :,6)
     &                  + giju(:,6)*A0gyi(:,11)
c
            gAgyi( :,7) = giju(:,4)*A0gyi( :,2)
     &                  + giju(:,2)*A0gyi( :,7)
     &                  + giju(:,6)*A0gyi(:,12)
c
            gAgyi( :,8) = giju(:,4)*A0gyi( :,3)
     &                  + giju(:,2)*A0gyi( :,8)
     &                  + giju(:,6)*A0gyi(:,13)
c
            gAgyi( :,9) = giju(:,4)*A0gyi( :,4)
     &                  + giju(:,2)*A0gyi( :,9)
     &                  + giju(:,6)*A0gyi(:,14)
c
            gAgyi(:,10) = giju(:,4)*A0gyi( :,5)
     &                  + giju(:,2)*A0gyi(:,10)
     &                  + giju(:,6)*A0gyi(:,15)
c
            gAgyi(:,11) = giju(:,5)*A0gyi( :,1)
     &                  + giju(:,6)*A0gyi( :,6)
     &                  + giju(:,3)*A0gyi(:,11)
c
            gAgyi(:,12) = giju(:,5)*A0gyi( :,2)
     &                  + giju(:,6)*A0gyi( :,7)
     &                  + giju(:,3)*A0gyi(:,12)
c
            gAgyi(:,13) = giju(:,5)*A0gyi( :,3)
     &                  + giju(:,6)*A0gyi( :,8)
     &                  + giju(:,3)*A0gyi(:,13)
c
            gAgyi(:,14) = giju(:,5)*A0gyi( :,4)
     &                  + giju(:,6)*A0gyi( :,9)
     &                  + giju(:,3)*A0gyi(:,14)
c
            gAgyi(:,15) = giju(:,5)*A0gyi( :,5)
     &                  + giju(:,6)*A0gyi(:,10)
     &                  + giju(:,3)*A0gyi(:,15)
c	
c... the denominator term of the DC factor
c... evaluation of the term  Y,i.A0DC Y,j 
c
            yiA0DCyj(:,1) = A0DC(:,1)*g1yi(:,1)**2
     &                    + two*g1yi(:,1)*A0DC(:,2)*g1yi(:,5)
     &                    + A0DC(:,3)*g1yi(:,2)**2
     &                    + A0DC(:,3)*g1yi(:,3)**2
     &                    + A0DC(:,3)*g1yi(:,4)**2
     &                    + A0DC(:,4)*g1yi(:,5)**2
c
            yiA0DCyj(:,2) = A0DC(:,1)*g2yi(:,1)**2
     &                    + two*g2yi(:,1)*A0DC(:,2)*g2yi(:,5)
     &                    + A0DC(:,3)*g2yi(:,2)**2
     &                    + A0DC(:,3)*g2yi(:,3)**2
     &                    + A0DC(:,3)*g2yi(:,4)**2
     &                    + A0DC(:,4)*g2yi(:,5)**2
c
            yiA0DCyj(:,3) = A0DC(:,1)*g3yi(:,1)**2
     &                    + two*g3yi(:,1)*A0DC(:,2)*g3yi(:,5)
     &                    + A0DC(:,3)*g3yi(:,2)**2
     &                    + A0DC(:,3)*g3yi(:,3)**2
     &                    + A0DC(:,3)*g3yi(:,4)**2
     &                    + A0DC(:,4)*g3yi(:,5)**2
c
            yiA0DCyj(:,4) = g1yi(:,1)*A0DC(:,1)*g2yi(:,1)
     &                    + g1yi(:,1)*A0DC(:,2)*g2yi(:,5)
     &                    + g1yi(:,2)*A0DC(:,3)*g2yi(:,2)
     &                    + g1yi(:,3)*A0DC(:,3)*g2yi(:,3)
     &                    + g1yi(:,4)*A0DC(:,3)*g2yi(:,4)
     &                    + g1yi(:,5)*A0DC(:,2)*g2yi(:,1)
     &                    + g1yi(:,5)*A0DC(:,4)*g2yi(:,5)
c
            yiA0DCyj(:,5) = g1yi(:,1)*A0DC(:,1)*g3yi(:,1)
     &                    + g1yi(:,1)*A0DC(:,2)*g3yi(:,5)
     &                    + g1yi(:,2)*A0DC(:,3)*g3yi(:,2)
     &                    + g1yi(:,3)*A0DC(:,3)*g3yi(:,3)
     &                    + g1yi(:,4)*A0DC(:,3)*g3yi(:,4)
     &                    + g1yi(:,5)*A0DC(:,2)*g3yi(:,1)
     &                    + g1yi(:,5)*A0DC(:,4)*g3yi(:,5)
c
            yiA0DCyj(:,6) = g2yi(:,1)*A0DC(:,1)*g3yi(:,1)
     &                    + g2yi(:,1)*A0DC(:,2)*g3yi(:,5)
     &                    + g2yi(:,2)*A0DC(:,3)*g3yi(:,2)
     &                    + g2yi(:,3)*A0DC(:,3)*g3yi(:,3)
     &                    + g2yi(:,4)*A0DC(:,3)*g3yi(:,4)
     &                    + g2yi(:,5)*A0DC(:,2)*g3yi(:,1)
     &                    + g2yi(:,5)*A0DC(:,4)*g3yi(:,5)
c
c.... ------------------------->  DC factor  <--------------------------
c
c	if ((ires .ne. 2) .or. (Jactyp .eq. 1)) then
c
c.... calculate 2-norm of Grad-local-V with respect to A0
c
c.... DC-mallet
c
            if (iDC .eq. 1) then
c
              fact = one
              if (ipord .eq. 2)  fact = 0.9
              if (ipord .eq. 3)  fact = 0.75
c
              gnorm = one / (
     &                giju(:,1)*yiA0DCyj(:,1)
     &              + two*giju(:,4)*yiA0DCyj(:,4)
     &              + two*giju(:,5)*yiA0DCyj(:,5)
     &              + giju(:,2)*yiA0DCyj(:,2) 
     &              + two*giju(:,6)*yiA0DCyj(:,6)
     &              + giju(:,3)*yiA0DCyj(:,3) 
     &              + epsM  )
c
c	            DC(:,intp)=dim((fact*sqrt(raLS*gnorm)),(rtLS*gnorm))
              DC(:,intp)=max(zero,(fact*sqrt(raLS*gnorm))-(rtLS*gnorm))
c
c.... flop count
c
!	    flops = flops + 46*npro
c
            endif
c
c.... DC-quadratic
c
            if (iDC .eq. 2) then
c
              gnorm = one / (
     &                giju(:,1)*yiA0DCyj(:,1)
     &              + two*giju(:,4)*yiA0DCyj(:,4)
     &              + two*giju(:,5)*yiA0DCyj(:,5)
     &              + giju(:,2)*yiA0DCyj(:,2) 
     &              + two*giju(:,6)*yiA0DCyj(:,6)
     &              + giju(:,3)*yiA0DCyj(:,3) 
     &              + epsM  )
         
c
              DC(:,intp) = two * rtLS * gnorm
c
c.... flop count
c
!	    flops = flops + 36*npro
c
            endif
c
c.... DC-min
c
            if (iDC .eq. 3) then
c
              fact = one
              if (ipord .eq. 2)  fact = pt5
c
              gnorm = one / (
     &                giju(:,1)*yiA0DCyj(:,1)
     &              + two*giju(:,4)*yiA0DCyj(:,4)
     &              + two*giju(:,5)*yiA0DCyj(:,5)
     &              + giju(:,2)*yiA0DCyj(:,2) 
     &              + two*giju(:,6)*yiA0DCyj(:,6)
     &              + giju(:,3)*yiA0DCyj(:,3) 
     &              + epsM  )
c
c	    DC(:,intp) = min( dim(fact * sqrt(raLS * gnorm),
              DC(:,intp) = min( max(zero,fact * sqrt(raLS * gnorm)-
     &                     rtLS * gnorm), two * rtLS * gnorm )
c
c.... flop count
c
!	    flops = flops + 48*npro
c
            endif
c
c	endif
c                  
        end subroutine calc_e3_dc_factor
c        
      end module e3_dc_func_m
