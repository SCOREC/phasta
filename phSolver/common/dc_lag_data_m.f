      module dc_lag_data_m
c-------------------------------------------------------------------------------
c  data structures for the discontinuious capturing lagging
c-------------------------------------------------------------------------------
        implicit none
        real*8, dimension(:), allocatable :: sum_dc_lag_vol
        real*8, dimension(:), allocatable :: sum_vol
        real*8, dimension(:,:), allocatable :: sum_dc_lag_l
        real*8, dimension(:,:), allocatable :: sum_vol_l
c
        real*8, dimension(:), allocatable :: dc_lag_g !global dc_lag for each node
        real*8, dimension(:), allocatable :: dc_lag_itr !global dc_lag for each node calculated at each Newton iteration
        real*8, dimension(:,:), allocatable :: dc_lag_l ! blk level dc_lag
        real*8, dimension(:), allocatable :: dc_lag_qt  ! dc_lag at one qudraduture point at element level
c
        real*8, dimension(:), allocatable :: vol_elm ! volumn of each interior element
        real*8 :: vol_factor ! volume of master element for each interior element type
c                        
      end module dc_lag_data_m
