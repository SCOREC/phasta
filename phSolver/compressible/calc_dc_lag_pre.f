      subroutine calc_dc_lag_pre(yold,  acold,  xold,  umeshold,
     &                           shp,   shgl)
c.......................................................................
c  Calculating the numerical viscousity introduced by discontinous 
c  capturing (DC) term using the desired fields(converged solution of last
c  time step) as a pre-processing.
c  The final output is a global numerical viscousity field for every element
c.......................................................................
        use pointer_data
        use e3_param_m
        use solid_data_m
        use e3_solid_func_m
        use timedataC
        use eqn_state_m
        use e3_solid_m
        use ifbc_def_m
        use dc_lag_data_m
c     
        include "common.h"
        include "mpif.h"
c
        real*8, dimension(nshg,ndof), intent(in) :: yold, acold
        real*8, dimension(numnp, nsd), intent(in) :: xold,  umeshold
        real*8, dimension(MAXTOP,maxsh,MAXQPT), intent(in) :: shp
        real*8, dimension(MAXTOP,nsd,maxsh,MAXQPT), intent(in) :: shgl
c
        real*8, dimension(nshg,nflow) :: tmp_res, tmp_rmes
        real*8, dimension(nshg,idflx) :: tmp_qres
        real*8, dimension(nshg,10) :: tmp_rerr
        real*8, dimension(nshg,nflow,nflow) :: tmp_BDiag 
        real*8, dimension(:,:), allocatable :: tmpshp
        real*8, dimension(:,:,:), allocatable :: tmpshgl
        real*8, dimension(:,:,:), allocatable :: tmp_egmass
c
c.... set the parameters for flux and surface tension calculations(cp from solgmr)
c
        idflx = 0 
        if(idiff >= 1)  idflx= idflx + (nflow-1) * nsd
        if (isurf == 1) idflx=idflx + nsd
c.... set up parameters(cp from elmgmr)
        ires   = 1        
c
c..... NEED TO ADD the global reconstruction of the diffusive flux 
c......vector, q, and lumped mass matrix, rmass in future
c
c... initialize the arrays(cp from elmgmr)
c
        tmp_res    = zero ! dummy array
        tmp_rmes   = zero ! to avoid trap_uninitialized, dummy array
        if (iprec .ne. 0) then 
          tmp_BDiag = zero ! dummy array
        endif
c          
        tmp_qres = zero !dummy array
        tmp_rerr = zero !dummy array
c          
        flxID = zero
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk          ! used in timeseries
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mater  = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
          if(lhs.eq.1) then
            nedof = nflow*nshl ! ensuring the value of nedof(total dof per element)
                               ! is consistent
                               ! with the element toplogy, not only accounting
                               ! for the max number of shp per element
                               ! (which reads in from geobc file) times nflow for
                               ! better memory efficiency
            allocate (tmp_egmass(npro,nedof,nedof)) !dummy array
            tmp_egmass = zero
          else
            allocate (tmp_egmass(1,1,1))
          endif          
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
          e3_malloc_ptr => e3_malloc
          e3_mfree_ptr => e3_mfree
c... for dc lagging in the pre-processing
          allocate(dc_lag_pre(npro))
          dc_lag_pre = zero
c... assigning material properties pointers
          select case (mat_eos(mater,1))
          case (ieos_ideal_gas,ieos_ideal_gas_2)
            getthm6_ptr => getthm6_ideal_gas
            getthm7_ptr => getthm7_ideal_gas
          case (ieos_ideal_gas_mixture)
            getthm6_ptr => getthm7_ideal_gas_mixture
            getthm7_ptr => getthm7_ideal_gas_mixture
          case (ieos_liquid_1)
            getthm6_ptr => getthm6_liquid_1
            getthm7_ptr => getthm7_liquid_1
          case (ieos_solid_1)
            getthm6_ptr => getthm6_solid_1
            getthm7_ptr => getthm7_solid_1
            iblk_solid = iblk 
            e3_malloc_ptr => e3_malloc_solid
            e3_mfree_ptr => e3_mfree_solid
          case (ieos_noble_abel)
            getthm6_ptr => getthm6_noble_abel
            getthm7_ptr => getthm7_noble_abel
         case default
            call error ('getthm  ', 'wrong material', mater)
          end select
c.. allocate arrays at element level
          if (associated(e3_malloc_ptr)) call e3_malloc_ptr          
c
          call AsIGMR (yold,                acold,
     &                 xold,                mxmudmi(iblk)%p,
     &                 tmpshp,
     &                 tmpshgl,             mien(iblk)%p,
     &                 mater,               tmp_res,
     &                 tmp_rmes,            tmp_BDiag,
     &                 tmp_qres,            tmp_egmass,
     &                 tmp_rerr,            umeshold )
c
c... mapping the numerical viscousity of dc from local to global
          do i = 1,npro
            dc_lag_g(mieMap(iblk)%p(i)) = dc_lag_pre(i)
          enddo
c... deallocation
          if (associated(e3_mfree_ptr)) call e3_mfree_ptr
c... for dc lagging
          deallocate(dc_lag_pre)
c
          deallocate ( tmp_egmass )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c                    
c.... end of interior element loop
c
        enddo
c             
      end subroutine calc_dc_lag_pre
