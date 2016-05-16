      module e3if_vi_m
c
        use e3if_defs_m
c
        implicit none
c
      contains
c
      subroutine calc_vi(p,n)
c
        real*8, pointer, intent(in) :: p(:), n(:,:)
c
        real*8, parameter :: alpha = 1.d0
        real*8, parameter :: beta  = 1.d-6
c
        integer :: isd
c
        do isd = 1,nsd
c
c          vi(:,isd) = beta * p**alpha * n(:,isd)
C          vi(:,isd) = n(:,isd)
      vi(:,1) =  0.0d0
      vi(:,2) =  0.0d0
      vi(:,3) =  0.0d0
c
        enddo
c
      end subroutine calc_vi
c
      subroutine calc_vi_area_node(sum_vi_area_l,shp,nshl)
c
        real*8, dimension(:,:,:), intent(inout) :: sum_vi_area_l
        real*8, dimension(:,:),   intent(in)    :: shp
        integer, intent(in) :: nshl
c
        integer :: n,isd
c
        do n = 1,nshl
          do isd = 1,nsd
            sum_vi_area_l(:,n,isd) = sum_vi_area_l(:,n,isd) 
     &                             + shp(:,n)*vi(:,isd)*area(:)
          enddo
          sum_vi_area_l(:,n,nsd+1) = sum_vi_area_l(:,n,nsd+1) + shp(:,n)*area(:)
        enddo
c
      end subroutine calc_vi_area_node
c
      end module e3if_vi_m
