c       module bc_on_vi_m
c
c       	implicit none
c
c	contains
c
        subroutine itrBCvi (actual_vi, iBC, BC, ilwork)
c--apply flow BCs on sum_vi_area
c
c----------------------------------------------------------------------
c
c This program satisfies the boundary conditions on the Y-variables.
c
c input:
c  y      (nshg,nflow)   : y variables 
c  iBC    (nshg)        : Boundary Condition Code
c  BC     (nshg,ndofBC) : boundary condition constraint parameters
c  ylimit (3,nflow)     : (1,:) limiting flag
c                         (2,:) lower bound
c                         (3,:) upper bound
c output:
c  y      (nshg,nflow)   : Adjusted V value(s) corresponding to a 
c                           constraint d.o.f.
c  umesh  (numnp,nsd)    : mesh velocity. FOR ALE 
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension actual_vi(nshg,3),             iBC(nshg),
     &            ac(nshg,nflow),            BC(nshg,ndofBC)

        dimension ilwork(nlwork)           
     
        dimension umesh(numnp,nsd)    !FOR ALE   
        integer   istp

        real*8 tmp(nshg), y1(nshg),q1(nshg)
        dimension rn1(nshg), rmagn1(nshg)
        real*8 limitcount(nflow)
        integer locmax(1),locmin(1)
c
c.... ------------------------->  Velocity  <--------------------------
c.... 3D
c
c.... x1-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 1)
            actual_vi(:,1) =  BC(:,3)  - BC(:,4) * actual_vi(:,2)
     &                         - BC(:,5) * actual_vi(:,3)
          endwhere
c
c.... x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 2)
            actual_vi(:,2) = BC(:,3)  - BC(:,4) * actual_vi(:,1)
     &                        - BC(:,5) * actual_vi(:,3)
          endwhere
c
c.... x1-velocity and x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 3)
            actual_vi(:,1) =  BC(:,3)  - BC(:,4) * actual_vi(:,3)
            actual_vi(:,2) =  BC(:,5)  - BC(:,6) * actual_vi(:,3)
          endwhere
c
c.... x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 4)
            actual_vi(:,3) = BC(:,3) - BC(:,4) * actual_vi(:,1)
     &                       - BC(:,5) * actual_vi(:,2)
          endwhere
c
c.... x1-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 5)
            actual_vi(:,1) = BC(:,3) - BC(:,4) * actual_vi(:,2)
            actual_vi(:,3) = BC(:,5) - BC(:,6) * actual_vi(:,2)
          endwhere
c
c.... x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 6)
            actual_vi(:,2) = BC(:,3)  - BC(:,4) * actual_vi(:,1)
            actual_vi(:,3) = BC(:,5)  - BC(:,6) * actual_vi(:,1)
          endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 7)
            actual_vi(:,1) =  BC(:,3)
            actual_vi(:,2) =  BC(:,4)
            actual_vi(:,3) =  BC(:,5) 
          endwhere
c
c       endif
c
c.... end of velocity
c

c.... communications
c 
c        if (numpe > 1) then
c           call commu (actual_vi, ilwork, nsd, 'out')
c           if(ires.ne.2) call commu (ac, ilwork, nflow, 'out')
c        endif
c
c       slave has masters value, for abc we need to rotate it
c
c        if(iabc==1) then        !are there any axisym bc's
c           call rotabc(y, iBC, 'out')
c           if(ires.ne.2) call rotabc(ac, iBC, 'out')
c        endif
c     
c
c.... return
c
        return
        end subroutine itrBCvi
c       end module bc_on_vi_m