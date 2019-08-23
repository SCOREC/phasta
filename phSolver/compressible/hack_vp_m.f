      module hack_vp_m
c---------------------------------------------------
c... controling the phase change rate on the surface
c---------------------------------------------------
        implicit none
c
        integer, allocatable, dimension(:) :: burn_info !determing if one
                                               !interface is burning nor not:
                                               ! 0 : not interface
                                               ! 1 : no burn interfac
                                               ! 2 : burn interface
        integer, allocatable, dimension(:) :: burn_element ! determine if one
                                                ! marco element if
                                                ! burning or not
        integer :: tot_burn_node
        integer, parameter :: burn_face = 6
        integer, parameter :: burn_edge = 2
        integer, parameter :: iso_face = 26
        integer, parameter :: iso_edge = 21
        real*8 :: burn_face_avg_x
      contains
c
        subroutine find_burn_face
c...............................................
c  Find the node on the burning interface and
c  counting the tot number of it
c...............................................
          use interfaceflag
          use conpar_m, only: nshg
          use m2gfields ! read m2g fields
          implicit none
c
          integer inode
c
          tot_burn_node = 0
          do inode = 1, nshg
            if (ifFlag(inode).eq.1) then
              if ( (m2gClsfcn(inode,1) .eq. 2) .and.
     &             (m2gClsfcn(inode,2) .eq. burn_face) ) then !this is
                                                              !burn
                                                              !interface
                tot_burn_node = tot_burn_node +1
                burn_info(inode) = 2
              else if ( (m2gClsfcn(inode,1) .eq. 1) .and.
     &                 (m2gClsfcn(inode,2) .eq. burn_edge) ) then !burn
                                                                  !edge
                tot_burn_node = tot_burn_node +1
                burn_info(inode) = 2
              else
                burn_info(inode) = 1
              endif
            endif
          enddo
        end subroutine find_burn_face
c
        subroutine calc_burn_face_location(x)
c...............................................
c  Find the avg x coord of the burning interface                                                         
c...............................................
          use conpar_m, only: nshg,numnp
          use number_def_m
          use global_const_m, only: nsd
          implicit none
c
          real*8, dimension(numnp,nsd),intent(in) :: x
c
          integer inode         
c
          burn_face_avg_x = zero
c
          do inode = 1, nshg
            if(burn_info(inode).eq.2) then
              burn_face_avg_x = burn_face_avg_x + (1/tot_burn_node)
     &                                          * x(inode,1)
            endif
          enddo
c
        end subroutine calc_burn_face_location   
      end module hack_vp_m

