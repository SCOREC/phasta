      module interface_continuity_data_m
c-------------------------------------------------------------------------------      
c... storing the data used for enforcing strong interface continuity 
c... of some variables
c--------------------------------------------------------------------------------
        implicit none
c
        integer, dimension(:), allocatable :: i_if_pair, i_if_con
        integer, parameter :: i_con_field = 5 ! temporaryly using this for the first stage of 
                                              ! development to determine which solution field is continuous
c      
      end module interface_continuity_data_m
