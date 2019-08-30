      module interface_pair_func_m
c-------------------------------------------------------------------------------      
c... storing the subroutines used for interface pair
c... of some variables
c--------------------------------------------------------------------------------
        implicit none
c
      contains
c
        subroutine alloc_init_interface_pair
c................................................................................
c... allocation and initialization of the pairing info
c................................................................................
          use interface_pair_data_m
          use interfaceflag
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
          do inode = 1,nshg
            i_if_pair(inode) = inode
            !if ((ifFlag(inode) .eq. 1)) then
            !  i_if_con(inode) = ibset(i_if_con(inode), 1) ! hacking the bit map to make T continuous
            !endif
          enddo
c          
          do iblk = 1, nelblif
            iel     = lcblkif(1, iblk)
            npro    = lcblkif(1,iblk+1) - iel
            !lcsyst0 = lcblkif(3, iblk)    ! element0 type
            !lcsyst1 = lcblkif(4, iblk)    ! element1 type
            ipord   = lcblkif(5, iblk)    ! polynomial order
            !nenl0   = lcblkif(6, iblk)    ! number of vertices per
            !element0
            !nenl1   = lcblkif(7, iblk)    ! number of vertices per
            !element1
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
              call error('alloc_init_interface_pair','higher order', ipord)
            endif
c... only support tet and wedge with triangle on the interface, otherwise, print out
c... error msg
            if(itpid .le. 4) then  !see L40 in global_param.f for details
              do i = 1, npro
                do j = 1, nshlb_if
                   j_pair = mod(2*j+1,3)+1 !flip the order for the local 2nd and 3rd node
                   i_if_pair(mienif0(iblk)%p(i,j)) = mienif1(iblk)%p(i,j_pair) ! making node in phase 1 the master node and
                                                                               ! flip the order for the local 2nd and 3rd node
                enddo
              enddo              
            else
              write(*,*) "interface element type is not supported"
              call error('alloc_init_interface_pair','topology',itpid)
            endif
c            
          enddo          
c       
        end subroutine alloc_init_interface_pair
c
        subroutine dealloc_interface_pair
          use interface_pair_data_m 
          implicit none
c
          deallocate(i_if_pair)
c                    
        end subroutine dealloc_interface_pair
c     
      end module interface_pair_func_m
