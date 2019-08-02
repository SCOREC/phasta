      subroutine SparseAp(iper, ilwork, iBC,  nvar, 
     &                    col,  row,    lhsK, p)   

C============================================================================
C
C "SparseAp": This routine performs the product of a sparsley stored matrix
C             lhsK(nvar*nvar, nnz_tot) and a vector p(nshg, nvar). The
C             results of the product is returned in p(nshg, nvar).
C
C Nahid Razmara, Spring 2000. (Sparse Matrix)
C============================================================================

c
      use timing_m
      include "common.h"
c
c
c.... Data declaration
c
	integer	col(nshg+1),	row(nnz*nshg), 
     &          iper(nshg),     ilwork(nlwork)
	real*8	lhsK(nvar*nvar,nnz_tot)
        real*8 	p(nshg,nvar),	q(nshg,nvar)
c
	real*8	tmp1,	tmp2,	tmp3,	tmp4,	tmp5
c     
c.... communicate:: copy the master's portion of uBrg to each slave
c
c      commu_start = MPI_Wtime()
      if (numpe > 1) then
         call commu (p, ilwork, nvar, 'out')
      endif
c      commu_stop = MPI_Wtime()
c      commu_time  = commu_time + commu_stop-commu_start
c
c.... local periodic boundary conditions (no communications)
c
        GMRES_start = MPI_Wtime()
        do j=1,nvar
           p(:,j)=p(iper(:),j)
        enddo
c
c       slave has masters value, for abc we need to rotate it
c        (if this is a vector only no SCALARS)
        if((iabc==1)) !are there any axisym bc's
     &     call rotabc(p(1,2), iBC,  'out')
c
c.... clear the vector
c
	q=zero
c
c.... Perform Matrix-Vector (AP) product
c
 	do i = 1, nshg 
	    do k = col(i), col(i+1)-1
c
              j = row(k) 
c
               do ii = 1, nvar
                  do jj = 1, nvar
                     id = (jj - 1) * nvar + ii
                     q(i,ii) = q(i,ii) + lhsK(id,k)*p(j,jj) 
                  enddo
               enddo
c
	    enddo
	enddo

        

c
	p =  q
c
c
c.... -------------------->   communications <-------------------------
c
c
        if((iabc==1))           !are there any axisym bc's
     &       call rotabc(p(1,2), iBC, 'in ')
c
        GMRES_stop = MPI_Wtime()
        GMRES_time = GMRES_time + GMRES_stop - GMRES_start
c
        if (numpe > 1) then
c
c.... send slave's copy of uBrg to the master
c
c           commu_start = MPI_Wtime()   
           call commu (p  , ilwork, nvar, 'in ')
c           commu_stop = MPI_Wtime()
c           commu_time = commu_time + commu_stop-commu_start
c     
c.... nodes treated on another processor are eliminated
c
           GMRES_start = MPI_Wtime()     
           numtask = ilwork(1)
           itkbeg = 1

           do itask = 1, numtask

              iacc   = ilwork (itkbeg + 2)
              numseg = ilwork (itkbeg + 4)

              if (iacc .eq. 0) then
                 do is = 1,numseg
                    isgbeg = ilwork (itkbeg + 3 + 2*is)
                    lenseg = ilwork (itkbeg + 4 + 2*is)
                    isgend = isgbeg + lenseg - 1
                    p(isgbeg:isgend,:) = zero
                 enddo
              endif
            
              itkbeg = itkbeg + 4 + 2*numseg

           enddo
           GMRES_stop = MPI_Wtime()
           GMRES_time = GMRES_time + GMRES_stop - GMRES_start
        endif
c
	return
	end

