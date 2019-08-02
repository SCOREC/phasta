      module timing_m
c-------------------------------------------------------------------
c... storing the data structure for the timing of GMRES solver
c-------------------------------------------------------------------
        implicit none
c
        real*8 :: tot_time, tot_start, tot_stop
        real*8 :: GMRES_time, GMRES_start, GMRES_stop
        real*8 :: commu_time, commu_start, commu_stop
        real*8 :: tmp_res
c
      end module timing_m
