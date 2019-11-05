      module if_velocity_m
c
c----------------------------------------
c    aims to calculate the interface
c    velocity at the global level
c----------------------------------------
        use workfc_m
        use number_def_m
        use pointer_data
        use blkdat_m
        use interfaceflag
        use interface_pair_data_m
        use hack_vp_m
c	use bc_on_vi_m
c
        implicit none
c
        real*8, pointer :: sum_vi_area(:,:)    ! interface velocity weighted by interface area
c
      contains
c
      subroutine init_sum_vi_area(nshg,nsd)
        integer, intent(in) :: nshg,nsd
        if (associated(sum_vi_area)) 
     &    deallocate (sum_vi_area)
        allocate (sum_vi_area(nshg,nsd+1))
        sum_vi_area(:,1:nsd) = zero
        sum_vi_area(:,1+nsd) = one
      end subroutine init_sum_vi_area
c
      subroutine destruct_sum_vi_area
        if (associated(sum_vi_area))
     &    deallocate(sum_vi_area)
      end subroutine destruct_sum_vi_area
c
      subroutine set_if_velocity 
     & (
     &  BC, iBC, umesh, disp, x, dt, ilwork,
     &  nshg, ndofBC, nsd, nelblif, nlwork, ndof,
     &  v)
c
        use dgifinp_m
        include "mpif.h"
c
        real*8,  intent(inout) ::  BC(nshg,ndofBC)
        integer, intent(inout) :: iBC(nshg)
        real*8,  dimension(nshg,nsd), intent(inout)    :: umesh, disp
        real*8,  dimension(nshg,nsd), intent(in)    :: x
        real*8, intent(in) :: dt
        integer, intent(in)    :: ilwork(nlwork)
        integer, intent(in) :: nshg, ndofBC, nsd, nelblif, nlwork, ndof
c
        integer :: iblk, iel, npro,inode, i0, i1, n, ierr
        integer, pointer :: ienif0(:,:), ienif1(:,:)
        real*8, dimension(nshg,3) :: actual_vi
        real*8, dimension(nshg,3),intent(in) :: v
c
        integer:: tot_front_edge_rank 
        real*8 :: sum_edge_x, sum_edge_x_rank !sum over all rank
c
c        if (numpe > 1) then
c          call commu (sum_vi_area(:,1:3), ilwork, nsd, 'in ')
c          call commu (sum_vi_area(:,4), ilwork, 1, 'in ')
c          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c        endif
c
c        if (numpe > 1) then
c          call commu (sum_vi_area(:,1:3), ilwork, nsd, 'out')
c          call commu (sum_vi_area(:,4), ilwork, 1, 'out')
c          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c        endif
c
c... calculate the avg x of the front edge
        call sum_front_edge_location(sum_edge_x,x)
c... communication
        if (numpe > 1) then
          call MPI_ALLREDUCE(sum_edge_x, sum_edge_x_rank, 1,
     &                       MPI_DOUBLE_PRECISION, MPI_SUM,
     &                       MPI_COMM_WORLD, ierr)
c
          call MPI_ALLREDUCE(tot_front_edge, tot_front_edge_rank, 1,
     &                       MPI_INTEGER, MPI_SUM,
     &                       MPI_COMM_WORLD, ierr)
        else
          sum_edge_x_rank = sum_edge_x
          tot_front_edge_rank = tot_front_edge
        endif
c
         if(tot_front_edge_rank .ne. 0) then
           burn_edge_avg_x = sum_edge_x_rank / tot_front_edge_rank
         else
           call error ('set_if_velocity', 'zero nodes on edge', tot_front_edge_rank) 
         endif 
c...
        actual_vi = zero
        do inode = 1, nshg
          if ( ifFlag(inode) .eq. 1 ) then
c            write(*,*) "rank",myrank,"i",inode,"x=",x(inode,:)
c            actual_vi(inode,:) = sum_vi_area(inode,:) / sum_vi_area(inode,nsd+1)
              if(burn_info(inode) .eq.2) then
                actual_vi(inode,1) = (-one*vi_mag  + v(i_if_pair(inode),1))
                actual_vi(inode,2) = (v(i_if_pair(inode),2))
                actual_vi(inode,3) = (v(i_if_pair(inode),3))
              else
                actual_vi(inode,1) = -one*(x(inode,1)/burn_edge_avg_x)
     &                               * vi_mag
     &                               + v(i_if_pair(inode),1) ! linearly distribute the phase change
                                                             ! rate along x axis
                actual_vi(inode,2) = (v(i_if_pair(inode),2))
                actual_vi(inode,3) = (v(i_if_pair(inode),3))
              endif
c
          endif
        enddo
c
        call itrBCvi ( actual_vi ,iBC ,BC )
c
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
c
100   format(a,'[',i2,'] ',i6,3f7.3,x,7e14.6)
200   format(a,'[',i2,'] ',i6,3e14.6)
      end subroutine set_if_velocity
c
      end module if_velocity_m
