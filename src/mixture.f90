! HB-Mixture library
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
!    separation: Calculates the distance between two position vectors (with PBC).
!    inv3x3: Invert a 3x3 matrix
!    gauss_prepare: Initialize all the parameters of the gaussian
!    gauss_logeval: Return the logarithm of the multivariate gaussian density
!    gauss_eval: Return the multivariate gaussian density
!    mixture_GetP: Return for each atom sh,sd and sa

      MODULE mixture
      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: dpigreco = (2.0d0*3.14159265358979d0)

      ! Structure that contains the parameters needed to define and
      ! estimate a gaussian
      TYPE gauss_type
         DOUBLE PRECISION lnorm ! logarithm of the normalization factor
         DOUBLE PRECISION det ! determinant of the covariance matrix
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: mean
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cov ! convariance matrix
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: icov ! inverse convariance matrix
      END TYPE

      CONTAINS

         SUBROUTINE collapsend(D,n1,n2,v1,v2)
            ! ndimensional version
            ! collapse v1 and v2 into v1
            !
            ! Args:
            !    D: dimensionality
            !    n1: size of the vector v1.
            !        n1 will be changed whend the routine will finish.
            !    n2: size of the vector v2
            !    v1: vector to wich I want to append v2
            !    v2: vector to append to v1

            INTEGER, INTENT(IN) :: D
            INTEGER, INTENT(INOUT) :: n1
            INTEGER, INTENT(IN) :: n2
            DOUBLE PRECISION, allocatable, dimension(:,:), INTENT(INOUT) :: v1
            DOUBLE PRECISION, dimension(:,:), INTENT(IN) :: v2

            DOUBLE PRECISION, allocatable, dimension(:,:) :: tmp_v

            IF(n1.EQ.0)THEN
               DEALLOCATE(v1)
               ALLOCATE(v1(D,n2))
               v1=v2
               n1=n2
            ELSE
               ALLOCATE(tmp_v(D,n1+n2))
               tmp_v(:,1:n1)=v1
               tmp_v(:,n1+1:n1+n2)=v2
               n1=n1+n2
               DEALLOCATE(v1)
               ALLOCATE(v1(D,n1))
               v1=tmp_v
               DEALLOCATE(tmp_v)
            ENDIF

         END SUBROUTINE collapsend

         SUBROUTINE collapse1d(n1,n2,v1,v2)
            ! 1D version
            ! collapse v1 and v2 into v1
            ! same as above, but 1D
            ! N.B. here I won't change n1, because I will change it calling
            ! N-Dim version
            INTEGER, INTENT(IN) :: n1
            INTEGER, INTENT(IN) :: n2
            DOUBLE PRECISION, allocatable, dimension(:), INTENT(INOUT) :: v1
            DOUBLE PRECISION, dimension(:), INTENT(IN) :: v2

            DOUBLE PRECISION, allocatable, dimension(:) :: tmp_v

            IF(n1.EQ.0)THEN
               DEALLOCATE(v1)
               ALLOCATE(v1(n2))
               v1=v2
            ELSE
               ALLOCATE(tmp_v(n1+n2))
               tmp_v(1:n1)=v1
               tmp_v(n1+1:n1+n2)=v2
               DEALLOCATE(v1)
               ALLOCATE(v1(n1+n2))
               v1=tmp_v
               DEALLOCATE(tmp_v)
            ENDIF

         END SUBROUTINE collapse1d

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
            IPIV(1) = 0
            IPIV(2) = 0
            IPIV(3) = 0
            ! call lapack
            call DGETRF(D,D,IM,D,IPIV,info)
            WORK(1) = 0.0d0
            WORK(2) = 0.0d0
            WORK(3) = 0.0d0
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

         SUBROUTINE gauss_prepare(D,gpars)
            ! Initialize all the parameters of the gaussian
            !
            ! Args:
            !    D: dimensionality
            !    gpars: gauss_type variable to initialize

            INTEGER, INTENT(IN) :: D
            TYPE(gauss_type), INTENT(INOUT) :: gpars

            ! calculate the determinant of the covariance matrix
            gpars%det = detmatrix(D,gpars%cov)
            ! calculate the inverse of the convariance matrix
            CALL invmatrix(D,gpars%cov,gpars%icov)

            ! calculate the  logarithm of the normalization factor
            gpars%lnorm = dlog(1.0d0/dsqrt((dpigreco**D)*gpars%det))
         END SUBROUTINE gauss_prepare

         ! probably this is no more needed
         DOUBLE PRECISION FUNCTION gauss_logeval(D,gpars, x)
            ! Return the logarithm of the multivariate gaussian density
            !
            ! Args:
            !    D: dimensionality
            !    gpars: gaussian parameters
            !    x: point in wich estimate the log of the gaussian

            INTEGER, INTENT(IN) :: D
            TYPE(gauss_type), INTENT(IN) :: gpars
            DOUBLE PRECISION, INTENT(IN) :: x(D)
            DOUBLE PRECISION dv(D),tmpv(D)
            DOUBLE PRECISION xcx

            dv=x-gpars%mean
            tmpv = matmul(dv,gpars%icov)
            xcx = dot_product(dv,tmpv)

            gauss_logeval = gpars%lnorm - 0.5d0*xcx
         END FUNCTION gauss_logeval

         DOUBLE PRECISION FUNCTION gauss_eval(D,gpars, x)
            ! Return the multivariate gaussian density
            ! Args:
            !    D: dimensionality
            !    gpars: gaussian parameters
            !    x: point in wich estimate the value of the gaussian

            INTEGER, INTENT(IN) :: D
            TYPE(gauss_type), INTENT(IN) :: gpars
            DOUBLE PRECISION, INTENT(IN) :: x(D)

            gauss_eval = dexp(gauss_logeval(D,gpars,x))
         END FUNCTION gauss_eval

         SUBROUTINE readgaussians(fileid,D,Nk,clusters,pks)
            ! Load the gaussians from the stream fileid
            !
            ! Args:
            !    fileid: the file containing the gaussians parameters
            !    D: dimensionality
            !    Nk: numeber of gaussians
            !    clusters: array of type_gaussian in wich we store the gaussians parameters
            !    lpks: logarithm of the Pks associated to the gaussians

            INTEGER, INTENT(IN) :: fileid
            INTEGER, INTENT(IN) :: D
            INTEGER, INTENT(INOUT) :: Nk
            TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: clusters
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: pks

            CHARACTER*1024 ::dummybuffer
            CHARACTER dummychar
            INTEGER k,i,j

            ! skip the first two comment lines , for now...
            READ(fileid,*) dummybuffer
            READ(fileid,*) dummybuffer
            ! read the number of gaussians form the file
            READ(fileid,*) Nk
            ! as always this the way that I know to allocate an dynamical array in a subroutine in fortran..
            DEALLOCATE(clusters)
            DEALLOCATE(pks)
            ALLOCATE(clusters(Nk),pks(Nk))
            DO k=1,Nk
              IF (.not.(ALLOCATED(clusters(k)%mean))) ALLOCATE(clusters(k)%mean(D))
              IF (.not.(ALLOCATED(clusters(k)%cov)))  ALLOCATE(clusters(k)%cov(D,D))
              IF (.not.(ALLOCATED(clusters(k)%icov))) ALLOCATE(clusters(k)%icov(D,D))
              ! read first the mean
              DO i=1,D
                 READ(fileid,'(A1,ES21.8E4)',ADVANCE='NO') dummychar,clusters(k)%mean(i)
              ENDDO
              ! read the covariance matrix
              DO i=1,D
                 DO j=1,D
                    READ(fileid,'(A1,ES21.8E4)',ADVANCE='NO') dummychar,clusters(k)%cov(i,j)
                 ENDDO
              ENDDO
              ! now read the pk and go to the next line
              READ(fileid,'(A1,ES21.8E4)') dummychar,pks(k)
              ! call the routine to prepare the gaussians
              CALL gauss_prepare(D,clusters(k))
            ENDDO

         END SUBROUTINE readgaussians

         SUBROUTINE writegaussianstofile(outputfile,D,nsamples,nminmax,tau,Nk,clusters,lpks)
            ! Write out the gaussian model informations to a file
            !
            ! Args:
            !    fileid: the file containing the gaussians parameters
            !    gaussp: type_gaussian container in wich we store the gaussian parameters
            !    lpk: logarithm of the Pk associated to the gaussian

            CHARACTER*1024, INTENT(IN) :: outputfile
            INTEGER, INTENT(IN) :: D
            INTEGER, INTENT(IN) :: nsamples
            INTEGER, INTENT(IN) :: nminmax
            INTEGER, INTENT(IN) :: Nk
            DOUBLE PRECISION, INTENT(IN) :: tau
            TYPE(GAUSS_TYPE), DIMENSION(Nk), INTENT(IN)  :: clusters
            DOUBLE PRECISION, DIMENSION(Nk), INTENT(IN) :: lpks

            INTEGER k,i,j

            OPEN(UNIT=12,FILE=trim(outputfile)//".gauss",STATUS='REPLACE',ACTION='WRITE')
            ! write a 2-lines header
	        WRITE(12,"(A31)",ADVANCE="NO") "# Quick Shift GM output. Ntot: "
	        WRITE(12,"(I12,A12,I11)",ADVANCE="NO") nsamples," , NVoroni: ",nminmax
	        WRITE(12,"(A8,ES21.8E4)") " , Tau: ", tau
	        WRITE(12,*) "# mean cov pk"
            ! read the number of gaussians form the file
            WRITE(12,*) Nk
            DO k=1,Nk
              ! write first the mean
              DO i=1,D
                 write(12,'(A1,ES21.8E4)',ADVANCE='NO') " ", clusters(k)%mean(i)
              ENDDO
              ! write the covariance matrix
              DO i=1,D
                 DO j=1,D
                   write(12,'(A1,ES21.8E4)',ADVANCE='NO') " ", clusters(k)%cov(i,j)
                 ENDDO
              ENDDO
              ! now read the pk and go to the next line
              write(12,'(A1,ES21.8E4)') " ", dexp(lpks(k))
           ENDDO
           CLOSE(UNIT=12)
         END SUBROUTINE writegaussianstofile

         SUBROUTINE ordergaussians(D,nk,clusters,pks,prif)
            ! Order the gaussians from the closest to prif
            ! Bubble-sort ordering is implemented here
            !
            ! Args:
            !    D: Dimensionality
            !    nk: number of gaussians to generate
            !    clusters: array containing gaussians parameters
            !    pks: The array containing the gaussians Pk
            !    prif: reference point

            INTEGER, INTENT(IN) :: D
            INTEGER, INTENT(IN) :: nk
            TYPE(gauss_type), DIMENSION(nk), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(nk), INTENT(INOUT) :: pks
            DOUBLE PRECISION, DIMENSION(D), INTENT(IN) :: prif

            TYPE(gauss_type) tmpgauss
            DOUBLE PRECISION distances(nk),tmpdistance,tmppk
			INTEGER j,i
			LOGICAL :: swapped = .TRUE.

			! calculate the distances
			DO i=1,nk
			   distances(i)=dot_product(clusters(i)%mean, prif)
            ENDDO
            ! now we can sort using the distances
            ! will use bubble sort
            DO j=nk-1,1,-1
               swapped = .FALSE.
               DO i = 1, j
                  IF (distances(i) > distances(i+1)) THEN
                     tmpdistance=distances(i)
                     distances(i)=distances(i+1)
                     distances(i+1)=tmpdistance
                     ! upgrade pks
                     tmppk=pks(i)
                     pks(i)=pks(i+1)
                     pks(i+1)=tmppk
                     ! upgrade also the clusters
                     tmpgauss=clusters(i)
                     clusters(i)=clusters(i+1)
                     clusters(i+1)=tmpgauss
                     swapped = .TRUE.
                  END IF
               END DO
               IF (.NOT. swapped) EXIT
            ENDDO
         END SUBROUTINE ordergaussians

         SUBROUTINE GetP(D,x,weight,alpha,nk,clusters,pks,pnks)
            ! Return for each atoms the sum of the ..
            !
            ! Args:
            !    D: Dimensionality
            !    x: The point in wich calculate the probabilities
            !    alpha: The smoothing factor
            !    nk: The number of gaussians in the mixture
            !    clusters: The array containing the structures with the gaussians parameters
            !    pks: The array containing the gaussians Pk
            !    pnks: The conditional probability of th point p given k

            INTEGER, INTENT(IN) :: D
            DOUBLE PRECISION, DIMENSION(D), INTENT(IN) :: x
            DOUBLE PRECISION, INTENT(IN) :: weight
            DOUBLE PRECISION, INTENT(IN) :: alpha
            INTEGER, INTENT(IN) :: nk
            TYPE(gauss_type), DIMENSION(nk), INTENT(IN) :: clusters
            DOUBLE PRECISION, DIMENSION(nk), INTENT(IN) :: pks
            DOUBLE PRECISION, DIMENSION(nk), INTENT(OUT) :: pnks

            DOUBLE PRECISION pnormpk !normalization factor
            INTEGER k

            pnks=0.0d0
            pnormpk=0.0d0 ! normalization factor (mixture weight)

            DO k=1,nk
               ! calculate the k probability (weighted) for the point (v,w,rad)
               ! and apply a smoothing elevating to alpha
               pnks(k) = (gauss_eval(D,clusters(k),x)*pks(k)*weight)**alpha
               ! calculate the mixture weight
               pnormpk = pnormpk+pnks(k)
            ENDDO
            ! skip cases in which the probability is tooooo tiny
            IF (pnormpk.NE.0.0d0) pnks = pnks/pnormpk ! normalization
         END SUBROUTINE GetP


      END MODULE mixture
