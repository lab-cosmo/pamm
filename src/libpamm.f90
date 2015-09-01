! Probabilistic Analysis of Molecular Motifs - Utility Library 
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


      MODULE libpamm
      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: twopi = (2.0d0*3.14159265358979d0)
      
      ! Structure that contains the parameters needed to define and
      ! estimate a gaussian
      TYPE gauss_type
         INTEGER D                                             ! dimensionality of the Gaussian
         DOUBLE PRECISION weight                               ! weight associated with the Gaussian cluster (not included in the normalization!)
         DOUBLE PRECISION lnorm                                ! logarithm of the normalization factor
         DOUBLE PRECISION det                                  ! determinant of the covariance matrix
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: mean   
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cov  ! convariance matrix
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: icov ! inverse convariance matrix
      END TYPE
      
      CONTAINS
		
      SUBROUTINE gauss_prepare(gpars)
         ! Initialize all the parameters of the gaussian
         !
         ! Args:
         !    gpars: gauss_type variable to initialize

         TYPE(gauss_type), INTENT(INOUT) :: gpars

         ! calculate the determinant of the covariance matrix
         gpars%det = detmatrix(gpars%D,gpars%cov)
         ! calculate the inverse of the convariance matrix
         CALL invmatrix(gpars%D,gpars%cov,gpars%icov)

         ! calculate the  logarithm of the normalization factor
         gpars%lnorm = dlog(1.0d0/dsqrt((twopi**gpars%D)*gpars%det))
      END SUBROUTINE gauss_prepare

      ! probably this is no more needed
      DOUBLE PRECISION FUNCTION gauss_logeval(gpars, x)
         ! Return the logarithm of the multivariate gaussian density
         !
         ! Args:
         !    gpars: gaussian parameters
         !    x: point in wich estimate the log of the gaussian
      
         TYPE(gauss_type), INTENT(IN) :: gpars
         DOUBLE PRECISION, INTENT(IN) :: x(gpars%D)
         DOUBLE PRECISION dv(gpars%D),tmpv(gpars%D)
         DOUBLE PRECISION xcx

         dv=x-gpars%mean
         tmpv = matmul(dv,gpars%icov)
         xcx = dot_product(dv,tmpv)

         gauss_logeval = gpars%lnorm - 0.5d0*xcx
      END FUNCTION gauss_logeval

      DOUBLE PRECISION FUNCTION gauss_eval(gpars, x)
         ! Return the multivariate gaussian density
         ! Args:
         !    gpars: gaussian parameters
         !    x: point in wich estimate the value of the gaussian
         
         TYPE(gauss_type), INTENT(IN) :: gpars
         DOUBLE PRECISION, INTENT(IN) :: x(gpars%D)

         gauss_eval = dexp(gauss_logeval(gpars,x))
      END FUNCTION gauss_eval
      
      SUBROUTINE pamm_p(x, pnks, nk, clusters, sig, alpha) 
         ! Computes for a configuration x the posterior probabilities 
         ! for it to belong to each of the PAMM clusters. 
         !
         ! Args:
         !    x: The point in wich calculate the probabilities
         !    alpha: The smoothing factor (defaults to one)
         !    nk: The number of gaussians in the mixture
         !    clusters: The array containing the structures with the gaussians parameters
         !    pks: The array containing the gaussians Pk
         !    pnks: The conditional probability of the point p given k
         !    bgsig: Background to be added to the mixture probability
         
         INTEGER, INTENT(IN) :: nk
         TYPE(gauss_type), INTENT(IN) :: clusters(nk)
         DOUBLE PRECISION, INTENT(IN) :: x(clusters(1)%d)
         DOUBLE PRECISION, INTENT(IN), OPTIONAL :: alpha
         DOUBLE PRECISION, INTENT(IN), OPTIONAL :: sig
         DOUBLE PRECISION, INTENT(OUT) :: pnks(nk)

         DOUBLE PRECISION pnormpk, palpha, mxpk, bgsig        
         INTEGER k
         
         palpha=1.0d0
         bgsig=0.00000001d0
         IF (PRESENT(alpha)) palpha = alpha
         IF (PRESENT(sig)) bgsig = sig
         
         pnks=0.0d0
         ! normalization factor (mixture weight)
         pnormpk=0.0d0 + bgsig ! add a background
         
         mxpk=-1d100
         DO k=1,nk
            ! optionally apply a smoothing based on alpha
            pnks(k) = gauss_logeval(clusters(k),x)            
            if (pnks(k).gt.mxpk) mxpk=pnks(k)
         ENDDO
         
         DO k=1,nk
            ! optionally apply a smoothing based on alpha
            pnks(k) = (dexp(pnks(k)-mxpk)*clusters(k)%weight)**palpha
            ! calculate the mixture weight
            pnormpk = pnormpk+pnks(k)
         ENDDO
         ! skip cases in which the probability is tooooo tiny
         ! IF (pnormpk.NE.0.0d0) 
         pnks = pnks/pnormpk ! normalize
      END SUBROUTINE pamm_p

      
      SUBROUTINE readclusters(fileid,Nk,clusters)
         ! Load the gaussian clusters from the stream fileid
         !
         ! Args:
         !    fileid: the file containing the gaussians parameters
         !    Nk: numeber of gaussians
         !    clusters: array of type_gaussian in wich we store the gaussians parameters

         INTEGER, INTENT(IN) :: fileid
         INTEGER, INTENT(OUT) :: Nk
         TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:), INTENT(OUT)  :: clusters         

         CHARACTER(LEN=1024) :: dummybuffer         
         INTEGER k,D

         ! skip the first two comment lines , for now...         
         READ(fileid,'(A)') dummybuffer
         DO WHILE ( dummybuffer(1:1) == '#' ) 
            READ(fileid,'(A)') dummybuffer
         ENDDO
         READ(dummybuffer, *) D, Nk
         
         IF (ALLOCATED(clusters)) DEALLOCATE(clusters)
         ALLOCATE(clusters(nk))
         DO k=1,Nk
           clusters(k)%D=D
           ALLOCATE(clusters(k)%mean(D))
           ALLOCATE(clusters(k)%cov(D,D))
           ALLOCATE(clusters(k)%icov(D,D))
           ! read first the mean
           READ(fileid,*) clusters(k)%weight, clusters(k)%mean, clusters(k)%cov
           ! call the routine to prepare the gaussians
           CALL gauss_prepare(clusters(k))
         ENDDO

      END SUBROUTINE readclusters

      SUBROUTINE writeclusters(outf,comments, nk, clusters) 
         ! Write out the gaussian model informations to a file. The file format
         ! is as follows:
         !
         ! # One or more comment lines
         ! dimensionality nclusters 
         ! weight1 mean11 mean12 mean13 ... mean1d cov111 cov112 ... cov11d ... cov1dd
         ! weight2 ....
         !
         ! Args:
         !    outf: the file id to which the gaussians parameters will be written
         !    comments: a string containing the comments section (including already the #)
         !    nk: number of clusters 
         !    clusters: type_gaussian container in wich we store the gaussian parameters

         INTEGER, INTENT(IN) :: outf
         CHARACTER(LEN=1024), INTENT(IN) :: comments
         INTEGER, INTENT(IN) :: nk
         TYPE(GAUSS_TYPE), DIMENSION(Nk), INTENT(IN)  :: clusters

         INTEGER k,i,j
         
         WRITE(outf, '(A)') trim(adjustl(comments))
         WRITE(outf, '(I12, I12)') clusters(1)%D, nk
         DO k=1,nk
           ! now read the pk and go to the next line
           write(outf,'(A1,ES21.8E4)',ADVANCE='NO') " ", clusters(k)%weight
           ! write first the mean
           DO i=1,clusters(1)%D
              write(outf,'(A1,ES21.8E4)',ADVANCE='NO') " ", clusters(k)%mean(i)
           ENDDO
           ! write the covariance matrix
           DO i=1,clusters(1)%D
              DO j=1,clusters(1)%D
                write(outf,'(A1,ES21.8E4)',ADVANCE='NO') " ", clusters(k)%cov(i,j)
              ENDDO
           ENDDO
           write(outf,*) ""
        ENDDO
      END SUBROUTINE writeclusters
      
      SUBROUTINE invmatrix(D,M,IM)
         ! inversion of a square matrix using lapack
         INTEGER, INTENT(IN) :: D
         DOUBLE PRECISION, DIMENSION(D,D), INTENT(IN) :: M
         DOUBLE PRECISION, DIMENSION(D,D), INTENT(OUT) :: IM

         DOUBLE PRECISION, DIMENSION(D) :: WORK
         INTEGER, DIMENSION(D) :: IPIV
         INTEGER info,error

         IM=M
         info = 0
         error = 0
         IPIV(:) = 0
         ! call lapack
         call DGETRF(D,D,IM,D,IPIV,info)
         WORK(:) = 0.0d0
         ! call lapack
         call DGETRI(D,IM,D,IPIV,WORK,D,info)
      END SUBROUTINE

      DOUBLE PRECISION FUNCTION detmatrix(D,M)
        ! determinant of a square matrix
         INTEGER, INTENT(IN) :: D
         DOUBLE PRECISION, DIMENSION(D,D) , INTENT(IN) :: M

         DOUBLE PRECISION, DIMENSION(D,D) :: matrix
         DOUBLE PRECISION :: mt, temp
         INTEGER :: i, j, k, l
         LOGICAL :: DetExists = .TRUE.

         matrix=M
         l = 1
         !Convert to upper triangular form
         DO k = 1, D-1
            IF (matrix(k,k) == 0) THEN
               DetExists = .FALSE.
               DO i = k+1,D
                  IF (matrix(i,k) /= 0) THEN
                     DO j=1,D
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                     END DO
                     DetExists = .TRUE.
                     l=-l
                     EXIT
                  ENDIF
               END DO
               IF (DetExists .EQV. .FALSE.) THEN
                   detmatrix=0.0d0
                   return
               END IF
            ENDIF
            DO j=k+1,D
                mt = matrix(j,k)/matrix(k,k)
                DO i=k+1,D
                    matrix(j,i) = matrix(j,i) - mt*matrix(k,i)
                END DO
            END DO
         ENDDO

         !Calculate determinant by finding product of diagonal elements
         detmatrix=l
         DO i=1,D
            detmatrix=detmatrix*matrix(i,i)
         ENDDO
      END FUNCTION detmatrix

      END MODULE libpamm
