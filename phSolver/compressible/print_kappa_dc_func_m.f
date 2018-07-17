      module print_kappa_dc_func_m
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
        implicit none
c
      contains
c
        subroutine alloc_init_print_kappa_dc
c...
c...
          use print_kappa_dc_data_m
          use conpar_m, only:numel
          use number_def_m
          implicit none
c
          allocate(kappa_dc(numel, 3, 3))
          allocate(kappa_dc_f(numel))
          allocate(pe_t_dc(numel))
c
          kappa_dc = zero
          kappa_dc_f = zero
          pe_t_dc = zero
c
        end subroutine alloc_init_print_kappa_dc
c
        subroutine dealloc_print_kappa_dc
c...
c...
          use print_kappa_dc_data_m
          implicit none
c
          deallocate(kappa_dc)
          deallocate(kappa_dc_f)
          deallocate(pe_t_dc)                  
c
        end subroutine  dealloc_print_kappa_dc
c
        subroutine alloc_int_print_kappa_dc_blk
c...
c...
          use print_kappa_dc_data_m
          use propar_m, only: npro
          use number_def_m
          implicit none
c
          allocate(kappa_dc_blk(npro, 3, 3))
          allocate(kappa_dc_f_blk(npro))
          allocate(pe_t_dc_blk(npro))
c
          kappa_dc_blk = zero
          kappa_dc_f_blk = zero
          pe_t_dc_blk = zero          
c
        end subroutine alloc_int_print_kappa_dc_blk
c
        subroutine dealloc_print_kappa_dc_blk
c...
c...
          use print_kappa_dc_data_m
          implicit none
c
          deallocate(kappa_dc_blk)
          deallocate(kappa_dc_f_blk)
          deallocate(pe_t_dc_blk)                  
c
        end subroutine  dealloc_print_kappa_dc_blk 
c
        subroutine calc_kappa_dc(A0, giju,  dc,  rho,
     &                           cp, u1,   u2,  u3,
     &                           xl, shape,intp,nshl,
     &                           ngauss )
c......................
c......................
          use propar_m, only: npro
          use number_def_m
          use global_const_m, only: nsd
          use elmpar_m, only: nenl
          use print_kappa_dc_data_m
          implicit none
c
          real*8, dimension(npro),intent(in) :: A0, rho, cp, u1, u2,
     &                                          u3
          real*8, dimension(npro,ngauss),intent(in) :: dc
          real*8, dimension(npro,6), intent(in) :: giju
          real*8, dimension(npro,nshl, nsd), intent(in) :: xl
          real*8, dimension(npro,nshl), intent(in) :: shape
          integer :: intp, nshl, ngauss
c
          real*8, dimension(npro,nsd) :: x_qt, norm
          real*8, dimension(npro) :: length, u_n
          real*8, dimension(npro,3,3) :: giju_f
          real*8, dimension(npro,nsd) :: temp1
          real*8, dimension(npro) :: temp2
          integer :: iel, i, j
c... get the coords of the quadrature points
          x_qt = zero
          do iel = 1, npro
            do i = 1, nshl
              do j = 1, nsd
                x_qt(iel,j) = x_qt(iel,j) + shape(iel,i)* xl(iel,i,j)
              enddo
            enddo
          enddo
c... get the out normal of the quadrature points
          do iel = 1, npro
            length(iel) = sqrt(dot_product( x_qt(iel,:), x_qt(iel,:) ) )
            do i = 1, nsd
              norm(iel,i) = x_qt(iel,i)/length(iel) 
            enddo
          enddo
c... get the velocity at the normal direction
          u_n(:) = u1(:)*norm(:,1) + u2(:)*norm(:,2) + u3(:)*norm(:,3)
c... convert the gij_u into 3 by 3
          giju_f(:,1,1) = giju(:,1)
          giju_f(:,1,2) = giju(:,4)
          giju_f(:,1,3) = giju(:,5)
          giju_f(:,2,1) = giju(:,4)
          giju_f(:,2,2) = giju(:,2)
          giju_f(:,2,3) = giju(:,6)
          giju_f(:,3,1) = giju(:,5)
          giju_f(:,3,2) = giju(:,6)
          giju_f(:,3,3) = giju(:,3)
c... local kappa_dc
          do iel = 1, npro
            kappa_dc_blk(iel,:,:) = kappa_dc_blk(iel,:,:) + ( dc(iel,intp)
     &                             * A0(iel)*giju_f(iel,:,:) )/ngauss
          enddo
c          
          do iel  = 1, npro
            temp1(iel,:) = matmul(giju_f(iel,:,:), norm(iel,:))
            temp2(iel) = dc(iel,intp) * A0(iel)
     &                 * dot_product(norm(iel,:), temp1(iel,:))
            kappa_dc_f_blk(iel) = kappa_dc_f_blk(iel) + temp2(iel)/ngauss
          enddo
c... local pe number
          do iel  = 1, npro
            pe_t_dc_blk(iel) = pe_t_dc_blk(iel) + rho(iel)*abs(u_n(iel))
     &                                          * cp(iel)
     &                                          *length(iel)/ temp2(iel)/ngauss
          enddo                                                                           
c
        end subroutine calc_kappa_dc        
      
c
      end module print_kappa_dc_func_m   
