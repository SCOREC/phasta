      module print_interface_error_data_m
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
        implicit none
c
        real*8, dimension(5) :: error_if_flux    ! surface averaged error for flux
        real*8, dimension(5) :: int_err_if_flux  ! surface intgeral of error for flux
        real*8, dimension(5) :: error_if_tan_1   ! surface averaged error for tangential quantities along x direction
        real*8, dimension(5) :: error_if_tan_2   ! surface averaged error for tangential quantities along y direction
        real*8, dimension(5) :: error_if_tan_3   ! surface averaged error for tangential quantities along z direction
        real*8, dimension(5) :: int_err_if_tan_1  ! surface intgeral of error for tangential quantities
        real*8, dimension(5) :: int_err_if_tan_2  ! surface intgeral of error for tangential quantities
        real*8, dimension(5) :: int_err_if_tan_3  ! surface intgeral of error for tangential quantities
        real*8 :: int_area                       ! area of interface
        real*8, dimension(5) :: int_err_if_flux_rank   !sum up all ranks
        real*8, dimension(5) :: int_err_if_tan_1_rank  !sum up all ranks
        real*8, dimension(5) :: int_err_if_tan_2_rank  !sum up all ranks
        real*8, dimension(5) :: int_err_if_tan_3_rank  !sum up all ranks
        real*8 :: int_area_rank  !sum up all ranks
c
        real*8, dimension(:,:), allocatable :: int_err_if_flux_blk  ! surface intgeral of error for flux at blk level
        real*8, dimension(:,:), allocatable :: int_err_if_tan_1_blk   ! surface intgeral of error for tangential quantities at blk level
        real*8, dimension(:,:), allocatable :: int_err_if_tan_2_blk   ! surface intgeral of error for tangential quantities at blk level
        real*8, dimension(:,:), allocatable :: int_err_if_tan_3_blk   ! surface intgeral of error for tangential quantities at blk level
        real*8, dimension(:), allocatable :: int_area_blk           ! area of interface at blk level
c
        integer :: err_flag                                                 
c      
      end module print_interface_error_data_m
