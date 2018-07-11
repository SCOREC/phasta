      module print_kappa_dc_data_m
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
        implicit none
c
        real*8, dimension(:,:,:), allocatable :: kappa_dc ! a global 3 by 3 matrix
        real*8, dimension(:), allocatable :: kappa_dc_f ! normal of the 3 by 3 matrix
        real*8, dimension(:), allocatable :: pe_t_dc ! peclet number for DC in energy equation
c
        real*8, dimension(:,:,:), allocatable :: kappa_dc_blk ! a blk level 3 by 3 matrix
        real*8, dimension(:), allocatable :: kappa_dc_f_blk
        real*8, dimension(:), allocatable :: pe_t_dc_blk 
c      
      end module print_kappa_dc_data_m
