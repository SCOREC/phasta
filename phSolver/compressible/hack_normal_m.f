      module hack_normal_m
c
c------------------------------------------------------------------------------
c  calculating the normal at quadrature points based on coordinates for sphere
c  problem
c------------------------------------------------------------------------------
        implicit none
c
c
      contains
        subroutine hack_normal(nv, xl, shp, nshl, nenl)
c-------------------------------------------------------------------------------
          use e3if_func_m, only:sum_qpt
          use propar_m, only: npro
          use global_const_m, only: nsd
          implicit none
c
          real*8, dimension(npro,nsd),intent(inout) :: nv
          real*8, dimension(npro,nenl,nsd),intent(in) :: xl
          real*8, dimension(npro,nshl), intent(in) :: shp
          integer, intent(in) :: nshl, nenl
c          
          real*8, dimension(npro,nsd) :: xl_qpt
          real*8, dimension(npro) :: temp_len
          integer :: iel
          
c          
c... get the global coords at qudrature points
          do iel = 1,npro
            xl_qpt(iel,1) = sum_qpt(nshl,xl(iel,:,1),shp(iel,:))
            xl_qpt(iel,2) = sum_qpt(nshl,xl(iel,:,2),shp(iel,:))
            xl_qpt(iel,3) = sum_qpt(nshl,xl(iel,:,3),shp(iel,:))
          
c... get the distrance from orgin to qpt
            temp_len(iel) = sqrt(xl_qpt(iel,1)*xl_qpt(iel,1)
     &                    + xl_qpt(iel,2)*xl_qpt(iel,2)
     &                    + xl_qpt(iel,3)*xl_qpt(iel,3))
c... get the normal since center of the sphere is origin
            nv(iel,1)  = xl_qpt(iel,1) / temp_len(iel)
            nv(iel,2)  = xl_qpt(iel,2) / temp_len(iel)
            nv(iel,3)  = xl_qpt(iel,3) / temp_len(iel)
          enddo            
c        
        end subroutine hack_normal     
      end module hack_normal_m     
