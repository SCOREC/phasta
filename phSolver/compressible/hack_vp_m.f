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
        integer, allocatable, dimension(:) :: front_edge_info ! the edge between burning and
                                                              ! non-burn surface.
                                                              ! 0: not on the edge
                                                              ! 1: on the edge
c
        integer :: tot_burn_node, tot_front_edge
        integer, parameter :: burn_face_num = 2
        integer, parameter :: burn_edge_num = 2
        integer, parameter :: burn_face(burn_face_num) = (/6,73/)
        integer, parameter :: burn_edge(burn_edge_num) = (/55,61/)
        integer, parameter :: front_edge = 61
c        integer, parameter :: iso_face = 26
c        integer, parameter :: iso_edge = 21
        real*8 :: burn_edge_avg_x
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
          integer inode,j
c
          tot_burn_node = 0
          tot_front_edge = 0
c
          do inode = 1, nshg
            if (ifFlag(inode).eq.1) then
              if (m2gClsfcn(inode,1) .eq. 2) then
                do j =1, burn_face_num 
                  if (m2gClsfcn(inode,2) .eq. burn_face(j))  then !this is
                                                                  !burn
                                                                  !interface
                    tot_burn_node = tot_burn_node +1
                    burn_info(inode) = 2
                  endif
                enddo
              else if (m2gClsfcn(inode,1) .eq. 1) then
c... adding the info to find the front edge
                if(m2gClsfcn(inode,2) .eq. front_edge) then !front edge
                  tot_front_edge = tot_front_edge + 1
                  front_edge_info(inode) = 1
                endif
c...
                do j =1, burn_edge_num 
                  if (m2gClsfcn(inode,2) .eq. burn_edge(j)) then !burn
                                                                 !edge
                    tot_burn_node = tot_burn_node +1
                    burn_info(inode) = 2
                  endif
                enddo
              else
                burn_info(inode) = 1
              endif
            endif
          enddo
        end subroutine find_burn_face
c
        subroutine sum_front_edge_location(sum_edge_x,x)
c...............................................
c  Find the avg x coord of the front edge                                                         
c...............................................
          use conpar_m, only: nshg,numnp
          use number_def_m
          use global_const_m, only: nsd
          implicit none
c
          real*8, dimension(numnp,nsd),intent(in) :: x
          real*8, intent(out) :: sum_edge_x
c
          integer inode         
c
          sum_edge_x = zero
c
          do inode = 1, nshg
            if(front_edge_info(inode) .eq. 1) then !this is front edge
                sum_edge_x = sum_edge_x + x(inode,1)
            endif      
          enddo
c
        end subroutine sum_front_edge_location   
      end module hack_vp_m

