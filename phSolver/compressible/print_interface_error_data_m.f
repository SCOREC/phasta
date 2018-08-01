      module print_interface_error_data_m
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
        implicit none
c
        real*8, dimension(5) :: error_if_flux    ! surface averaged error for flux
        real*8, dimension(5) :: int_err_if_flux  ! surface intgeral of error for flux
        real*8, dimension(5) :: error_if_tan     ! surface averaged error for tangential quantities
        real*8, dimension(5) :: int_err_if_tan   ! surface intgeral of error for tangential quantities
c
        real*8, dimension(:,:), allocatable :: int_err_if_flux_blk  ! surface intgeral of error for flux at blk level
        real*8, dimension(:,:), allocatable :: int_err_if_tan_blk   ! surface intgeral of error for tangential quantities at blk level
c      
      end module print_interface_error_data_m
