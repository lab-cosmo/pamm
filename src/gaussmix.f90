! Gaussian Mixtur Model library
!
! Copyright (C) 2014, Piero Gasparotto and Michele Ceriotti
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
! Functions:
!    estep: Expectation step 
!    mstep: Maximizastion step

      MODULE gaussmix
         USE gaussian
      IMPLICIT NONE
      
      CONTAINS
      
         SUBROUTINE estep(nk,nsamples,vwda,clusters,lpks,pnk,loglike)
            ! Expectation step: get a new responsibility matrix and a new
            ! likelihood
            !
            ! Args:
            !    nk: number of gaussians 
            !    nsamples: number of points
            !    vwda: array containing the data
            !    clusters: array containing gaussians parameters
            !    lpks: array containing the logarithm of the Pk for each gaussian
            !    pnk: responsibility matrix
            !    loglike: logarithm of the likelihood
            
            INTEGER, INTENT(IN) :: nk
            INTEGER, INTENT(IN) :: nsamples
            DOUBLE PRECISION, DIMENSION(3,nsamples), INTENT(IN) :: vwda
            TYPE(gauss_type), DIMENSION(nk), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(nk,nsamples), INTENT(INOUT) :: pnk
            DOUBLE PRECISION, DIMENSION(nk), INTENT(INOUT) :: lpks
            DOUBLE PRECISION, INTENT(OUT) :: loglike
           
            ! variables used during the log-sum-exp trick
            DOUBLE PRECISION zmax,som,lognormfactor
            
            INTEGER i,j,k
    
            DO k=1,nk
               ! probably we can remove this, because the gaussians are already prepared
               call gauss_prepare(clusters(k))
               IF(lpks(k)==0.0d0)THEN
                  lpks(k)=-10000 ! we don't want to get NaN
               ENDIF
               DO i=1,nsamples
                  ! the gaussian density will often be so small as to underlflow to zero, so 
                  ! we will work with the logharitms of the densities rather then the densities themselves
                  pnk(k,i)=gauss_logeval(clusters(k),vwda(:,i))+lpks(k)
               ENDDO
            ENDDO
            
            ! We now have the unnormalized logs of the nk probalities for eanch point (pnk).
            ! To normalize we have to calculate the 'mixture weight' and 
            ! we will use the log-sum-exp formula
            
            loglike=0.0d0 ! initialize the log-likelihood to zero
            
            DO i=1,nsamples               
               zmax=MAXVAL(pnk(:,i),1)
               som=0.0d0
               DO k=1,nk
                  som=som+dexp(pnk(k,i)-zmax)
               ENDDO
               ! mixture weight
               lognormfactor=zmax+dlog(som)
               DO k=1,nk
                  pnk(k,i)=dexp(pnk(k,i)-lognormfactor) ! normalization
               ENDDO
               loglike=loglike+lognormfactor
            ENDDO      
	     END SUBROUTINE
	     
	     SUBROUTINE mstep(nk,nsamples,vwda,clusters,lpks,pnk,smooth)
            ! Maximization step: get the new gaussians parameters after the
            ! estep.
            ! 
            ! Args:
            !    nk: number of gaussians 
            !    nsamples: number of points
            !    vwda: array containing the data
            !    clusters: array containing gaussians parameters
            !    lpks: array containing the logarithm of the Pk for each gaussian
            !    pnk: responsibility matrix
            !    smooth: smoothing factor to stabilize the covariance matrix
            
            INTEGER, INTENT(IN) :: nk
            INTEGER, INTENT(IN) :: nsamples
            DOUBLE PRECISION, DIMENSION(3,nsamples), INTENT(IN) :: vwda
            TYPE(gauss_type), DIMENSION(nk), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(nk,nsamples), INTENT(IN) :: pnk
            DOUBLE PRECISION, DIMENSION(nk), INTENT(INOUT) :: lpks            
            DOUBLE PRECISION, INTENT(IN) :: smooth

            DOUBLE PRECISION, DIMENSION(3,nsamples) :: vtmp
            DOUBLE PRECISION som,som_n,sompnk,ctr
            INTEGER d,m,k
            
            DO k=1,nk
               sompnk=0.0d0
               sompnk=sum(pnk(k,:))
               ! estimate the fraction of all data points in the component k
               lpks(k)=dlog(sompnk/nsamples)
               
			      DO d=1,3
			         som_n=0.0d0
			         vtmp(d,:)=pnk(k,:)*vwda(d,:)
			         som_n=sum(vtmp(d,:))
			         ! estimate the new mean for the gaussian k
			         clusters(k)%mean(d)=som_n/sompnk
			         som_n=0.0d0
			         DO m=1,3
			            som_n=sum(pnk(k,:)*(vwda(d,:)-clusters(k)%mean(d))* &
			                    (vwda(m,:)-clusters(k)%mean(m)))
			            ! estimate the new covariance matrix for the gaussian k
			            clusters(k)%cov(m,d)=som_n/sompnk
			         ENDDO
			      ENDDO 
			      
			      ! to stabilize the covariance matrix
			      IF(smooth.ne.0.0d0)THEN
			         ctr=(clusters(k)%cov(1,1) + clusters(k)%cov(2,2) + clusters(k)%cov(3,3)) &
			             /3*smooth;
			         clusters(k)%cov=clusters(k)%cov*(1.-smooth)
			         DO d=1,3
			            clusters(k)%cov(d,d)=clusters(k)%cov(d,d)+ctr
			         ENDDO
			      ENDIF
            ENDDO
	     END SUBROUTINE 

      END MODULE gaussmix
