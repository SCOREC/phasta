      module dc_lag_data_m
c-------------------------------------------------------------------------------
c  data structures for the discontinuious capturing lagging
c-------------------------------------------------------------------------------
        implicit none
c
        real*8, dimension(:), allocatable :: dc_lag_g !global dc_lag for each element
        real*8, dimension(:), allocatable :: dc_lag_pre ! block level dc_lag for each element, used in preprocessing only
        real*8, dimension(:), allocatable :: dc_lag_blk ! block level dc_lag
        integer :: dc_calc_flag !flag of the preprocessing for dc lag                   
      end module dc_lag_data_m
