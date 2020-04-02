      module print_interface_error_data_m
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
        implicit none
c
        real*8, dimension(5) :: error_if_flux    ! surface averaged error for flux
        real*8, dimension(5) :: int_err_if_flux  ! surface intgeral of error for flux
        real*8, dimension(3) :: error_if_tan   ! surface averaged error for tangential quantities, sum up x,y,z direction in the L2 norm
        real*8, dimension(3) :: int_err_if_tan  ! surface intgeral of error for tangential quantities, sum up x,y,z direction
        real*8 :: int_area                       ! area of interface
c
        real*8, dimension(5) :: error_flux_nomalized
        real*8, dimension(5) :: int_flux
        real*8, dimension(3) :: error_tan_nomalized
        real*8, dimension(3) :: int_y
c        
        real*8, dimension(5) :: int_err_if_flux_rank   !sum up all ranks
        real*8, dimension(3) :: int_err_if_tan_rank  !sum up all ranks
        real*8 :: int_area_rank  !sum up all ranks
c
        real*8, dimension(5) :: int_flux_rank !sum up all ranks
        real*8, dimension(3) :: int_y_rank !sum up all ranks
        
c
        real*8, dimension(:,:), allocatable :: int_err_if_flux_blk  ! surface intgeral of error for flux at blk level
        real*8, dimension(:,:), allocatable :: int_err_if_tan_blk   ! surface intgeral of error for tangential quantities at blk level, sum up x, y, z direction
        real*8, dimension(:), allocatable :: int_area_blk           ! area of interface at blk level
c
        real*8, dimension(:,:), allocatable :: int_flux_blk  ! surface intgeral of flux at blk level
        real*8, dimension(:,:), allocatable :: int_y_blk  ! surface intgeral of flux at blk level        
c
        integer :: err_flag                                                 
c      
      end module print_interface_error_data_m
