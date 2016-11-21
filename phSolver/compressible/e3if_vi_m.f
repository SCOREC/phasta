      module e3if_vi_m
c
        use dgifinp_m
        use e3if_param_m
c
        implicit none
c
      contains
c
      subroutine calc_vi(p,n,u)
c
        real*8, pointer, intent(in) :: p(:), n(:,:), u(:,:)
        real*8 :: vi_max
c
        integer :: isd,i
        real*8 :: vimag,v1,t1,c1
        real*8 :: pvap,cprod,cdest
c
        select case (vi_ramping)
        case (no_ramp)
          c1 = one
        case (linear_ramp)
          t1 = ramp_time
          v1 = vimag
c
          if (time <= t1) then
            c1 = time/t1
          else
            c1 = one
          endif
        case default
          call error ('ERROR in e3if_vi:',' vi_ramping is not supported,',vi_ramping)
        end select
c
        select case (phase_change_model)
        case (no_vi)
          vi = zero
        case (const_vi)
          vi = vi_mag * n
        case (vieilles_burning)
c
          vi(:,1) = burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * n(:,1)
          vi(:,2) = burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * n(:,2)
          vi(:,3) = burn_rate_coeff*(p/burn_rate_pref)**burn_rate_exp * n(:,3)
c
        case (cavitation)
c
c >>> NOTE: HARD CODED NUMBERS HERE
c
          pvap = 110.0d3
          cdest = 1.d-4
          cprod = 1.d-4
c
c          vi(:,1) = -nv0(:,1)*(cdest*min(zero,pres0(i)-pvap) + cprod*max(zero,pres0-pvap))
c          vi(:,2) = -nv0(:,2)*(cdest*min(zero,pres0(i)-pvap) + cprod*max(zero,pres0-pvap))
c          vi(:,3) = -nv0(:,3)*(cdest*min(zero,pres0(i)-pvap) + cprod*max(zero,pres0-pvap))
          vi(:,1) = nv1(:,1)*(cprod*max(zero,pres0-pvap))
          vi(:,2) = nv1(:,2)*(cprod*max(zero,pres0-pvap))
          vi(:,3) = nv1(:,3)*(cprod*max(zero,pres0-pvap))
c      vi_max = 10.0d0
c      vi(:,1) = nv1(:,1)*(min(vi_max,cprod*max(zero,pres0-pvap)))
c      vi(:,2) = nv1(:,2)*(min(vi_max,cprod*max(zero,pres0-pvap)))
c      vi(:,3) = nv1(:,3)*(min(vi_max,cprod*max(zero,pres0-pvap)))
c
        case default
          call error ('ERROR in e3if_vi:',' phase_change_model is not supported.',phase_change_model)
        end select
c
c... add flow velocity
c
        vi(:,1) = c1*(vi(:,1) + u(:,1))
        vi(:,2) = c1*(vi(:,2) + u(:,2))
        vi(:,3) = c1*(vi(:,3) + u(:,3))
c        vi(:,1) = c1 * (vi(:,1) + pt50 * (u0(:,1)+u1(:,1)))
c        vi(:,2) = c1 * (vi(:,2) + pt50 * (u0(:,2)+u1(:,2)))
c        vi(:,3) = c1 * (vi(:,3) + pt50 * (u0(:,3)+u1(:,3)))
c      write(*,*) (norm2(vi(i,:)),i=1,npro)
c
      end subroutine calc_vi
c
      subroutine calc_vi_area_node(sum_vi_area_l,shp,WdetJif,nshl)
c
        real*8, dimension(:,:,:), intent(inout) :: sum_vi_area_l
        real*8, dimension(:,:),   intent(in)    :: shp
        real*8, dimension(:),     intent(in)    :: WdetJif
        integer, intent(in) :: nshl
c
        integer :: n,isd
c
        do n = 1,nshl
          do isd = 1,nsd
            sum_vi_area_l(:,n,isd) = sum_vi_area_l(:,n,isd) 
     &                             + shp(:,n)*vi(:,isd)*area(:)
c     &                             + shp(:,n)*vi(:,isd)*WdetJif(:)
          enddo
          sum_vi_area_l(:,n,nsd+1) = sum_vi_area_l(:,n,nsd+1) + shp(:,n)*area(:)
c          sum_vi_area_l(:,n,nsd+1) = sum_vi_area_l(:,n,nsd+1) + shp(:,n)*WdetJif(:)
        enddo
c
      end subroutine calc_vi_area_node
c
      end module e3if_vi_m
