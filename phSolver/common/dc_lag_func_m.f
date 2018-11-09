      module dc_lag_func_m
c-------------------------------------------------------------------------------
c  functions and subroutines for the lagging DC
c-------------------------------------------------------------------------------
        implicit none
c
      contains
c
        subroutine alloc_init_dc_lag
c........................................................................
c  allocation and initialization of the global data structure for DC lag
c.......................................................................
          use conpar_m, only: nshg
          use dc_lag_data_m, only: sum_dc_lag_vol, sum_vol,
     &                             dc_lag_g, dc_lag_itr
          use number_def_m
          implicit none
c
          allocate(sum_dc_lag_vol(nshg))
          allocate(sum_vol(nshg))
          allocate(dc_lag_itr(nshg))
c
          dc_lag_itr = zero
c        
        end subroutine alloc_init_dc_lag
c
c
        subroutine dealloc_dc_lag
c.......................................................................
c deallocation
c......................................................................
          use dc_lag_data_m, only: sum_dc_lag_vol, sum_vol,
     &                             dc_lag_g, dc_lag_itr
          implicit none
c
          deallocate(sum_dc_lag_vol)
          deallocate(sum_vol) 
          deallocate(dc_lag_g)
          deallocate(dc_lag_itr)
c                   
        end subroutine dealloc_dc_lag
c
c
        subroutine get_dc_lag_qt(dc_lag_qt, dc_lag_l, shp, nshl)
c.......................................................................
c get the lagging DC at each quadrature point
c.......................................................................
          use propar_m, only: npro
          use number_def_m
          implicit none
c
          real*8, dimension(npro),intent(out) :: dc_lag_qt
          real*8, dimension(npro,nshl),intent(in) :: dc_lag_l
          real*8, dimension(npro,nshl),intent(in) :: shp
          integer, intent(in)::nshl
c
          integer :: iel, n
c
          dc_lag_qt = zero
c        
          do iel = 1,npro
            do n = 1,nshl
              dc_lag_qt(iel) = dc_lag_qt(iel)+shp(iel,n)*dc_lag_l(iel,n)
            enddo
          enddo
c
        end subroutine get_dc_lag_qt
        
c
      end module dc_lag_func_m
