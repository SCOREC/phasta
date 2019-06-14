      module interface_velocity_m
c----------------------------------------------------------------------
c    aims to calculate the interface
c    velocity at the global level
c----------------------------------------------------------------------      
        implicit none
c
        real*8, dimension(:,:), allocatable :: vi_normal_global
c
      contains
c
        subroutine calc_interface_vi(vi, u, p, nv)
c.......................................................................
c... calculate interface velocity based on the input thermodynamic variables
c.......................................................................
        use dgifinp_m
        use global_const_m, only: nsd
        use e3if_param_m, only: time
        use number_def_m
        implicit none
c
        real*8, dimension(nsd), intent(inout) :: vi
        real*8, dimension(nsd), intent(in) :: u, nv
        real*8, intent(in) :: p
c        
        real*8, dimension(nsd) :: v1
        real*8 :: c1, t1
c... determine ramping the phase change rate or not
        select case (vi_ramping)
        case (no_ramp)
          c1 = one
        case (linear_ramp)
          t1 = ramp_time
          v1 = vi_mag
c
          if (time <= t1) then
            c1 = time/t1
          else
            c1 = one
          endif
c          
        case default
          call error ('ERROR in e3if_vi:',' vi_ramping is not supported,',vi_ramping)
        end select        
c... determine the phase change model, only have the zero, constant and vielle's law here
        select case (phase_change_model)
        case (no_vi)
          vi = zero
        case (const_vi)
          vi(1) = c1 * vi_mag * nv(1)
          vi(2) = c1 * vi_mag * nv(2)
          vi(3) = c1 * vi_mag * nv(3)
        case (vieilles_burning)
          vi(1) = c1*burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * nv(1)
          vi(2) = c1*burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * nv(2)
          vi(3) = c1*burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * nv(3)
        case default
          call error ('ERROR in e3if_vi:',' phase_change_model is not supported.',phase_change_model)
        end select
c... adding flow velocity to the interface velocity
        vi(1) = vi(1) + u(1)
        vi(2) = vi(2) + u(2)
        vi(3) = vi(3) + u(3)
c                
        end subroutine calc_interface_vi
c
        subroutine set_interface_velocity(umesh,  y,    BC, iBC,
     &                                    ilwork, nlwork)
c.......................................................................        
c... calculate the interface velocity at the global level for every node
c.......................................................................
          use global_const_m, only: nsd
          use conpar_m, only: nshg, ndof
          use genpar_m, only: ndofBC
          use weighted_normal_data_m, only:w_normal_global
          use interface_pair_data_m
          use interfaceflag
          use workfc_m, only:numpe
          use number_def_m
          include "mpif.h"
c
          real*8,  dimension(nshg,nsd), intent(inout) :: umesh
          real*8,  dimension(nshg,ndof), intent(in) :: y
          real*8,  intent(in) ::  BC(nshg,ndofBC)
          integer, intent(in) :: iBC(nshg)
          integer, intent(in) :: nlwork !! this could be improved by further modulize the common
          integer, intent(in) :: ilwork(nlwork)
c
          real*8,  dimension(nshg,nsd) :: actual_vi
          integer :: inode, jnode, ierr
c
          actual_vi = zero
c... calculate the interface velocity for each node
c... assuming phase 0 is the master node of the interface pair and always be
c... the lighter phase (air)          
          do inode = 1, nshg
            if ((ifFlag(inode) .eq. 1).and.(i_if_pair(inode) .ne. inode)) then ! this is a interface slave (liquid)
              jnode = i_if_pair(inode) !master node number
              call calc_interface_vi(actual_vi(inode,:), y(inode,1:3),
     &                             y(jnode,4), w_normal_global(jnode,:)) ! check the slots of the y,
                                                                         ! make sure it is consistant with the other part of the code
              actual_vi(jnode,:) = actual_vi(inode,:) !cp to master(gas) side
            endif
          enddo
c... satisfying the BC
        call itrBCvi ( actual_vi ,iBC ,BC )
c... communication from master rank to all slave rank (like in itrBC, dealing y array) !! needs double check
        if (numpe > 1) then
           call commu (actual_vi, ilwork, nsd, 'out')
           call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        endif        
c... setting the mesh velocity on the interface to be the interface velocity
        do inode = 1, nshg
          if ( ifFlag(inode) .eq. 1 ) then
            umesh(inode,:) = actual_vi(inode,:)
c
c.... the following line is moved to solve mesh part
c.... since in restart case, we should update interface mesh BC
c.... before we solve the mesh
c            BC(inode,ndof+2:ndof+4) = umesh(inode,:) * dt
          endif
        enddo
c                  
        end subroutine set_interface_velocity

        
      end module interface_velocity_m
