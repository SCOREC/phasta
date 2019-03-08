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
          use conpar_m, only: numel
          use dc_lag_data_m, only: dc_lag_g, dc_calc_flag
          use number_def_m
          implicit none
c
          allocate(dc_lag_g(numel))
c
          dc_lag_g = zero
          dc_calc_flag = 0
c        
        end subroutine alloc_init_dc_lag
c
c
        subroutine dealloc_dc_lag
c.......................................................................
c deallocation
c......................................................................
          use dc_lag_data_m, only: dc_lag_g
          implicit none
c
          deallocate(dc_lag_g)
c                   
        end subroutine dealloc_dc_lag
c
c
      end module dc_lag_func_m
