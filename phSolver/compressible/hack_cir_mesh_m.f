      module hack_cir_mesh_m
c--------------------------------------------------------
c... hacking so that the mesh displacment (velocity) in 
c... the tangential direction on the circular surface are zero
c--------------------------------------------------------------
        implicit none
c
      contains
c
        subroutine enfore_cir_mesh(disp, umesh)
c.............................................................
c... enforing the mesh displacement in the tangentail direction
c... of the desired faces and edges to be zero
c ..............................................................
          use m2gfields ! read m2g fields
          use hackcir_m
          use global_const_m
          use conpar_m
          use number_def_m
          implicit none
c
          real*8, dimension(numnp,nsd), intent(inout) :: disp
          real*8, dimension(numnp,nsd), intent(inout) :: umesh
c
          integer :: i,j
c         
          if( (cir_num_face_tag .gt. 0) .or.
     &         (cir_num_edge_tag .gt. 0) ) then ! double check
            do i = 1, nshg
              if ( (m2gClsfcn(i,1) .eq. 3 ).or.
     &             (m2gClsfcn(i,1) .eq. 0)  ) then ! region or vertex
                cycle
              else if (m2gClsfcn(i,1) .eq. 2) then ! face
                do j = 1, cir_num_face_tag
                  if (m2gClsfcn(i,2) .eq. cir_face_tag(j)) then
                    select case(cir_axis_flag)
                    case(1) ! long axis aligned with x
                      disp(i,2) = zero
                      disp(i,3) = zero
                      umesh(i,2) = zero
                      umesh(i,3) = zero
                    case(2) ! long axis aligned with y
                      disp(i,1) = zero
                      disp(i,3) = zero
                      umesh(i,1) = zero
                      umesh(i,3) = zero
                    case(3) ! long axis aligned with z
                      disp(i,1) = zero
                      disp(i,2) = zero
                      umesh(i,1) = zero
                      umesh(i,2) = zero
                    case default
                      call error ('enfore_cir_tan  ',
     &                     'wrong axis', cir_axis_flag)
                    end select                     
                  endif
                enddo ! end loop all desired faces
              else if (m2gClsfcn(i,1) .eq. 1) then ! edge
                do j = 1, cir_num_edge_tag
                  if (m2gClsfcn(i,2) .eq. cir_edge_tag(j)) then
                    select case(cir_axis_flag)
                    case(1) ! long axis aligned with x
                      disp(i,2) = zero
                      disp(i,3) = zero
                      umesh(i,2) = zero
                      umesh(i,3) = zero
                    case(2) ! long axis aligned with y
                      disp(i,1) = zero
                      disp(i,3) = zero
                      umesh(i,1) = zero
                      umesh(i,3) = zero
                    case(3) ! long axis aligned with z
                      disp(i,1) = zero
                      disp(i,2) = zero
                      umesh(i,1) = zero
                      umesh(i,2) = zero
                    case default
                      call error ('enfore_cir_tan  ',
     &                     'wrong axis', cir_axis_flag)
                    end select                        
                  endif
               enddo ! end loop over all desired edges
              endif ! end checking entity dimension
ccc.. for testing
c              select case(cir_axis_flag)
c                    case(1) ! long axis aligned with x
c                      disp(i,2) = zero
c                      disp(i,3) = zero
c                      umesh(i,2) = zero
c                      umesh(i,3) = zero
c                    case(2) ! long axis aligned with y
c                      disp(i,1) = zero
c                      disp(i,3) = zero
c                      umesh(i,1) = zero
c                      umesh(i,3) = zero
c                    case(3) ! long axis aligned with z
c                      disp(i,1) = zero
c                      disp(i,2) = zero
c                      umesh(i,1) = zero
c                      umesh(i,2) = zero
c                    case default
c                      call error ('enfore_cir_tan  ',
c     &                     'wrong axis', cir_axis_flag)
c                    end select
            enddo ! end loop over nshg
          else
            call error ('enfore_cir_tan  ',
     &                  'No face or edge specified', 
     &                  cir_num_face_tag)
          endif ! end of double checking                     
c        
        end subroutine enfore_cir_mesh
c 
      end module hack_cir_mesh_m
