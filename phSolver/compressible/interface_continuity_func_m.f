      module interface_continuity_func_m
c-------------------------------------------------------------------------------      
c... storing the subroutines used for enforcing strong interface continuity 
c... of some variables
c--------------------------------------------------------------------------------
        implicit none
c
      contains
c
        subroutine alloc_init_interface_continuity
c................................................................................
c... allocation and initialization of the pairing info
c................................................................................
          use interface_continuity_data_m
          use conpar_m, only:nshg
          use elmpar_m, only:nelblif
          use blkdat_m, only:lcblkif,iblkif_topology 
          use pointer_data, only:mienif0, mienif1
          implicit none
c
          integer :: inode, iblk, iel, npro, ipord, nenbl_if,
     &               itpid, nshlb_if, i, j, j_pair  
c... allocation
          allocate(i_if_pair(nshg))
c... initialization
c          do iblk = 1, nelblk
c                  !nenl   = lcblk(5,iblk)   ! no. of vertices per element
c                  iel    = lcblk(1,iblk)
c                  !lelCat = lcblk(2,iblk)
c                  !lcsyst = lcblk(3,iblk)
c                  !iorder = lcblk(4,iblk)
c                  nshl   = lcblk(10,iblk)
c                  !mater  = lcblk(7,iblk)
c                  !ndofl  = lcblk(8,iblk)
c                  !nsymdl = lcblk(9,iblk)
c                  npro   = lcblk(1,iblk+1) - iel
c                  !ngauss = nint(lcsyst)
c            ien => mien(iblk)%p
c            do iel = 1, npro
c              do ishl = 1, nshl
c                i_if_pair(ien(iel, ishl)) = ien(iel, ishl)
c                    enddo
c            enddo
c          enddo
c  
          do inode = 1,nshg
            i_if_pair(inode) = inode
          enddo
c          
          do iblk = 1, nelblif
            iel     = lcblkif(1, iblk)
            npro    = lcblkif(1,iblk+1) - iel
            !lcsyst0 = lcblkif(3, iblk)    ! element0 type
            !lcsyst1 = lcblkif(4, iblk)    ! element1 type
            ipord   = lcblkif(5, iblk)    ! polynomial order
            !nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
            !nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
            nenbl_if= lcblkif(8, iblk)    ! number of vertices on the interface
            !mater0  = lcblkif(9, iblk)
            !mater1  = lcblkif(10,iblk)
            !nshl0   = lcblkif(iblkif_nshl0,iblk)
            !nshl1   = lcblkif(iblkif_nshl1,iblk)
            itpid   = lcblkif(iblkif_topology,iblk)
            !ngaussif = nintif(itpid)
c... print out error msg if it is not linear element
            if(ipord .eq. 1) then
              nshlb_if   = nenbl_if ! only work with linear element
            else if(ipord .gt. 1) then
              write(*,*) "need to implement for higher order"
              call error('alloc_init_interface_continuity','higher order', ipord)
            endif
c... only support tet and wedge with triangle on the interface, otherwise, print out
c... error msg
            if(itpid .le. 4) then
              do i = 1, npro
                do j = 1, nshlb_if
                   j_pair = mod(2*j+1,3)+1 !flip the order for the local 2nd and 3rd node
                   i_if_pair(mienif1(iblk)%p(i,j)) = mienif0(iblk)%p(i,j_pair) ! making node in phase 0 the master node and
                                                                               ! flip the order for the local 2nd and 3rd node
                enddo
              enddo              
            else
              write(*,*) "interface topology is not supported"
              call error('alloc_init_interface_continuity','topology', itpid)
            endif
c            
          enddo          
c	
        end subroutine alloc_init_interface_continuity
c
        subroutine dealloc_interface_continuity
          use interface_continuity_data_m 
          implicit none
c
          deallocate(i_if_pair)
c                    
        end subroutine dealloc_interface_continuity
c
        subroutine itr_interface_continuity(y, ac)
c...............................................................................
c... ensuring the continuous interface field during the correct and update stage for
c... each non-linear iteration
          use interface_continuity_data_m
          use interfaceflag
          use conpar_m, only:ndof,nshg
          use genpar_m, only:ires
          implicit none
c
          real*8, dimension(nshg,ndof) :: y, ac
          integer :: i,j
c... handle the continuous field across interface(no communications)
c... Notice the order of the solution field is changed at itrdrv level:
c... y(:,1:3) - velocity, y(:,4) - pressure, y(:,5) - temperature 
          do j = 1,nshg
            if ( (ifFlag(j) .eq. 1) .and. 
     &           (i_if_pair(j) .ne. j) ) then !if j is interface pair slave
              i = i_if_pair(j)
              y(j,i_con_field) = y(i,i_con_field)
c              
              if(ires.ne.2) then
                ac(j,i_con_field) = ac(i,i_con_field)
              endif         
            endif
          enddo 
c          
        endsubroutine itr_interface_continuity
c     
      end module interface_continuity_func_m
