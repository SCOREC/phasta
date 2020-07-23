      module hack_cir_mesh_m
c--------------------------------------------------------
c... hacking so that the mesh displacment (velocity) in 
c... the tangential direction on the circular surface are zero
c--------------------------------------------------------------
        implicit none
c
      contains
c
        subroutine enfore_cir_mesh(iBC, BC)
c.............................................................
c... enforing the mesh displacement in the tangentail direction
c... of the desired faces and edges to be zero
c ..............................................................
          use m2gfields ! read m2g fields
          use hackcir_m
          use global_const_m
          use conpar_m
          use genpar_m
          use number_def_m
          implicit none
c
          integer*4, dimension(nshg), intent(inout) :: iBC
          real*8, dimension(nshg,ndofBC), intent(inout) :: BC
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
c... change the ibc and bc array value for the mesh
                      if ( (ibits(iBC(i),14,3) .eq. 2) .or.
     &                     (ibits(iBC(i),14,3) .eq. 4) ) then
                         iBC(i) = ibset(iBC(i),15) !turn on the y,z components
                         iBC(i) = ibset(iBC(i),16) ! of the mesh
                         BC(i,7:10) = zero !change the bc value to make
                                           !y,z component to be zero
                      endif
c... change the ibc and bc array value for the flow
                      if ( (ibits(iBC(i),3,3) .eq. 2) .or.
     &                     (ibits(iBC(i),3,3) .eq. 4) ) then
                         iBC(i) = ibset(iBC(i),4) !turn on the y,z components
                         iBC(i) = ibset(iBC(i),5) ! of the mesh
                         BC(i,3:6) = zero !change the bc value to make
                                           !y,z component to be zero
                      endif
                    case(2) ! long axis aligned with y
c... change the ibc and bc array value for the mesh
                      if ( (ibits(iBC(i),14,3) .eq. 1) .or.
     &                     (ibits(iBC(i),14,3) .eq. 4) ) then
                         iBC(i) = ibset(iBC(i),14) !turn on the x,z components
                         iBC(i) = ibset(iBC(i),16) ! of the mesh
                         BC(i,7:10) = zero !change the bc value to make
                                           !x,z component to be zero
                      endif
c... change the ibc and bc array value for the flow
                      if ( (ibits(iBC(i),3,3) .eq. 1) .or.
     &                     (ibits(iBC(i),3,3) .eq. 4) ) then
                         iBC(i) = ibset(iBC(i),3) !turn on the x,z components
                         iBC(i) = ibset(iBC(i),5) ! of the mesh
                         BC(i,3:6) = zero !change the bc value to make
                                           !x,z component to be zero
                      endif
                    case(3) ! long axis aligned with z
c... change the ibc and bc array value for the mesh
                      if ( (ibits(iBC(i),14,3) .eq. 1) .or.
     &                     (ibits(iBC(i),14,3) .eq. 2) ) then
                         iBC(i) = ibset(iBC(i),14) !turn on the x,y components
                         iBC(i) = ibset(iBC(i),15) ! of the mesh
                         BC(i,7:10) = zero !change the bc value to make
                                           !x,y component to be zero
                      endif
c... change the ibc and bc array value for the flow
                      if ( (ibits(iBC(i),3,3) .eq. 1) .or.
     &                     (ibits(iBC(i),3,3) .eq. 2) ) then
                         iBC(i) = ibset(iBC(i),3) !turn on the x,y components
                         iBC(i) = ibset(iBC(i),4) ! of the mesh
                         BC(i,3:6) = zero !change the bc value to make
                                           !x,y component to be zero
                      endif
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
c... change the ibc and bc array value for the flow
                      if ( (ibits(iBC(i),3,3) .eq. 2) .or.
     &                     (ibits(iBC(i),3,3) .eq. 4) ) then
                         iBC(i) = ibset(iBC(i),4) !turn on the y,z components
                         iBC(i) = ibset(iBC(i),5) ! of the mesh
                         BC(i,3:6) = zero !change the bc value to make
                                           !y,z component to be zero
                      endif
                   case(2) ! long axis aligned with y
c... change the ibc and bc array value for the flow
                      if ( (ibits(iBC(i),3,3) .eq. 1) .or.
     &                     (ibits(iBC(i),3,3) .eq. 4) ) then
                         iBC(i) = ibset(iBC(i),3) !turn on the x,z components
                         iBC(i) = ibset(iBC(i),5) ! of the mesh
                         BC(i,3:6) = zero !change the bc value to make
                                           !x,z component to be zero
                      endif
                    case(3) ! long axis aligned with z
c... change the ibc and bc array value for the flow
                      if ( (ibits(iBC(i),3,3) .eq. 1) .or.
     &                     (ibits(iBC(i),3,3) .eq. 2) ) then
                         iBC(i) = ibset(iBC(i),3) !turn on the x,y components
                         iBC(i) = ibset(iBC(i),4) ! of the mesh
                         BC(i,3:6) = zero !change the bc value to make
                                           !x,y component to be zero
                      endif
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
