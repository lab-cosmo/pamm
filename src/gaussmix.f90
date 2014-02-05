! This performs ...
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
!    X_functions: Description

      MODULE gaussmix
         USE gaussian
      IMPLICIT NONE
      
      CONTAINS
      
         SUBROUTINE estep(Nk,NSamples,vwda,clusters,lpk,pnk,loglike)
            ! Initialize the gaussian calculating...
            ! ...
            ! Args:
            !    param: descript 
            INTEGER, INTENT(IN) :: Nk
            INTEGER, INTENT(IN) :: NSamples
            DOUBLE PRECISION, DIMENSION(3,NSamples), INTENT(IN) :: vwda
            TYPE(GAUSS_TYPE), DIMENSION(Nk), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(Nk,NSamples), INTENT(INOUT) :: pnk
            DOUBLE PRECISION, DIMENSION(Nk), INTENT(INOUT) :: lpk
            DOUBLE PRECISION, INTENT(OUT) :: loglike
           
            DOUBLE PRECISION zmax,som,lognormfactor
            INTEGER i,j,k
    
            DO k=1,Nk
               ! get the inverse of the cov matrix and the normalization factor
               call gauss_prepare(clusters(k))
               IF(lpk(k)==0.0d0)THEN
                  lpk(k)=-1000
               ENDIF
               write(*,*) k, Nk, clusters(k)%mean
               DO i=1,NSamples
                  pnk(k,i)=gauss_logeval(clusters(k),vwda(:,i))+lpk(k)
               ENDDO
            ENDDO
            
            loglike=0.0d0
            ! normalization
            DO i=1,NSamples               
               zmax=MAXVAL(pnk(:,i),1)
               som=0.0d0
               DO k=1,Nk
                  som=som+dexp(pnk(k,i)-zmax)
               ENDDO
               lognormfactor=zmax+dlog(som)
               DO k=1,Nk
                  pnk(k,i)=dexp(pnk(k,i)-lognormfactor)
               ENDDO
               loglike=loglike+lognormfactor
            ENDDO
            
	     END SUBROUTINE
	     
	     SUBROUTINE mstep(Nk,NSamples,vwda,clusters,lpk,pnk,smooth)
            ! Initialize the gaussian calculating...
            ! ...
            ! Args:
            !    param: descript 
            INTEGER, INTENT(IN) :: Nk
            INTEGER, INTENT(IN) :: NSamples
            DOUBLE PRECISION, DIMENSION(3,NSamples), INTENT(IN) :: vwda
            TYPE(GAUSS_TYPE), DIMENSION(Nk), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(Nk,NSamples), INTENT(IN) :: pnk
            DOUBLE PRECISION, DIMENSION(Nk), INTENT(INOUT) :: lpk            
            DOUBLE PRECISION, INTENT(IN) :: smooth
            
            DOUBLE PRECISION, DIMENSION(3,3) :: mID
            DOUBLE PRECISION, DIMENSION(3,NSamples) :: vtmp
            DOUBLE PRECISION som,som_n,sompnk,ctr
            INTEGER d,m,k
            
            mID = 1.0d0
            
            DO k=1,Nk
               sompnk=0.0d0
               sompnk=sum(pnk(k,:))
               lpk(k)=dlog(sompnk/NSamples)
               
			      DO d=1,3
			         som_n=0.0d0
			         vtmp(d,:)=pnk(k,:)*vwda(d,:)
			         som_n=sum(vtmp(d,:))
			         clusters(k)%mean(d)=som_n/sompnk
			         som_n=0.0d0
			         DO m=1,3
			            som_n=sum(pnk(k,:)*(vwda(d,:)-clusters(k)%mean(d))* &
			                    (vwda(m,:)-clusters(k)%mean(m)))
			            
			            clusters(k)%cov(m,d)=som_n/sompnk
			         ENDDO
			      ENDDO 
			      
			      ! to stabilize
			      ctr=(clusters(k)%cov(1,1) + clusters(k)%cov(2,2) + clusters(k)%cov(3,3))/3*smooth;
			      clusters(k)%cov=clusters(k)%cov*(1.-smooth)
			      DO d=1,3
			         clusters(k)%cov(d,d)=clusters(k)%cov(d,d)+ctr
			      ENDDO
            ENDDO
            
	     END SUBROUTINE
	    

      END MODULE gaussmix
