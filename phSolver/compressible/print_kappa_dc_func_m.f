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
          implicit none
c
          allocate(kappa_dc_blk(npro, 3, 3))
          allocate(kappa_dc_f_blk(npro))
          allocate(pe_t_dc_blk(npro))          
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
      end module print_kappa_dc_func_m   
