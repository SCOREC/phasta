        subroutine calc_gc_normal ( x,     shpb,
     &                              ienb,  iBCB,  normal)
c
        use BLparameters
        include "common.h"
c
        dimension x(numnp,nsd),
     &            shpb(nshl,ngaussb),
     &            ienb(npro,nshl),
     &            iBCB(npro,ndiBCB)
c
        real*8, dimension(nshg, nsd) :: normal
c
        integer   calc_factor(npro)
        integer   ielm, counter
c
c.... collect only all the first nenbl vertices has BLflag
c
       calc_factor(:) = 0
       do ielm = 1, npro
         counter = 0
c
         do i = 1, nenbl
           if ( BLflag(ienb(ielm, i)) .eq. 1 ) then
             counter = counter + 1
           endif
         enddo ! end loop bounadry vertice in an element
c
         if (counter .eq. nenbl) then
           calc_factor(ielm) = 1
         endif
       enddo
c
c.... assemble the normal vector
c
        call calc_normal(x, shpb, calc_factor, ienb, iBCB, normal)
c
c.... end
c
        return
        end

c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine setBLbc( disp, iBC, BC )
c
        include "common.h"
c
c.... please only pass mesh elas BC (i, ndof+2:ndof+5) into this subroutine
c
        dimension disp(nsd), BC(4)
        integer   iBC
c
            select case (ibits(iBC,14,3))
            case (1) ! x1 direction
              BC(1) = disp(1)
            case (2) ! x2 direction
              BC(1) = disp(2)
            case (3) ! x1 & x2 direction
              BC(1) = disp(1)
              BC(3) = disp(2)
            case (4) ! x3 direction
              BC(1) = disp(3)
            case (5) ! x1 & x3 direction
              BC(1) = disp(1)
              BC(3) = disp(3)
            case (6) ! x2 & x3 direction
              BC(1) = disp(2)
              BC(3) = disp(3)
            case (7) ! x1 & x2 & x3 direction
              BC(1) = disp(1)
              BC(2) = disp(2)
              BC(3) = disp(3)
            end select
c.... end
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine iBCelas_to_dbl( iBCelas, flag )
c
        real*8  flag
        integer iBCelas
c
        flag = REAL(ibits(iBCelas,14,3))
c
c.... end
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine dbl_to_iBCelas( flag, iBCelas )
c
        real*8  flag
        integer iBCcase, iBCelas
c
        iBCcase = INT(flag + 0.5)
c
          select case (iBCcase)
          case (1) ! x1 direction
            iBCelas = ibset(iBCelas, 14)
            iBCelas = ibclr(iBCelas, 15)
            iBCelas = ibclr(iBCelas, 16)
          case (2) ! x2 direction
            iBCelas = ibclr(iBCelas, 14)
            iBCelas = ibset(iBCelas, 15)
            iBCelas = ibclr(iBCelas, 16)
          case (3) ! x1 & x2 direction
            iBCelas = ibset(iBCelas, 14)
            iBCelas = ibset(iBCelas, 15)
            iBCelas = ibclr(iBCelas, 16)
          case (4) ! x3 direction
            iBCelas = ibclr(iBCelas, 14)
            iBCelas = ibclr(iBCelas, 15)
            iBCelas = ibset(iBCelas, 16)
          case (5) ! x1 & x3 direction
            iBCelas = ibset(iBCelas, 14)
            iBCelas = ibclr(iBCelas, 15)
            iBCelas = ibset(iBCelas, 16)
          case (6) ! x2 & x3 direction
            iBCelas = ibclr(iBCelas, 14)
            iBCelas = ibset(iBCelas, 15)
            iBCelas = ibset(iBCelas, 16)
          case (7) ! x1 & x2 & x3 direction
            iBCelas = ibset(iBCelas, 14)
            iBCelas = ibset(iBCelas, 15)
            iBCelas = ibset(iBCelas, 16)
          end select
c
c.... end
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
       subroutine layered_mesh_correction( x,  xold,  iBC,  BC,
     &                     colm,     rowp,    meshq,
     &                     HBrg,     eBrg,    yBrg,
     &                     Rcos,     Rsin,    iper,   ilwork,
     &                     shp,      shgl,    shpb,   shglb,
     &                     shpif)
c
        include "common.h"
c
        dimension x(numnp,nsd),     xold(numnp,nsd),
     &            iBC(nshg),        BC(nshg,ndofBC),
     &            ilwork(nlwork),   iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)
c
        real*8, dimension(maxtop, maxsh, maxqpt) :: shpif
c
        dimension HBrg(Kspace+1,Kspace),   eBrg(Kspace+1),
     &            yBrg(Kspace),            Rcos(Kspace),
     &            Rsin(Kspace)
c
        integer  colm(nshg+1),        rowp(nnz*nshg)
c
        real*8   umesh(numnp,nsd),    meshq(numel),
     &           disp(numnp, nsd),    elasDy(nshg,nelas)
c
          lhs = 1
          iprec=lhs
c
          call itrBCElas(umesh,  disp,  iBC,
     &                   BC(:,ndof+2:ndof+5),
     &                   iper,   ilwork)
c

          call SolGMRElas (x,        disp,    iBC,    BC,
     &                     colm,     rowp,    meshq,
     &                     HBrg,     eBrg,    yBrg,
     &                     Rcos,     Rsin,    iper,   ilwork,
     &                     shp,      shgl,    shpb,   shglb,
     &                     shpif,    elasDy)
c
          call itrCorrectElas(disp, elasDy)
c
          call itrBCElas(umesh,  disp,  iBC,
     &                   BC(:,ndof+2:ndof+5),
     &                   iper,   ilwork)
c
          x = x + disp
          xold = x
c
          if (myrank .eq. 0) write(*,*) "layered mesh corrected"
c
c.... end
c
        return
        end

