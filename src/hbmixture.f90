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
!    hbmixture_GetGMMP: Return for each atom sh,sd and sa

      MODULE hbmixture
      IMPLICIT NONE
      
      DOUBLE PRECISION, PARAMETER :: dpigreco = (2.0d0*3.14159265358979d0)

      ! Types used by bitwise operators to control the atom type
      ! they must be power of 2
      INTEGER, PARAMETER :: TYPE_NONE=0
      INTEGER, PARAMETER :: TYPE_H=1
      INTEGER, PARAMETER :: TYPE_DONOR=2
      INTEGER, PARAMETER :: TYPE_ACCEPTOR=4
      
      ! Structure that contains the parameters needed to define and
      ! estimate a gaussian
      TYPE gauss_type
         DOUBLE PRECISION lnorm ! logarithm of the normalization factor
         DOUBLE PRECISION det ! determinant of the covariance matrix
         DOUBLE PRECISION, DIMENSION(3) :: mean
         DOUBLE PRECISION, DIMENSION(3,3) :: cov ! convariance matrix
         DOUBLE PRECISION, DIMENSION(3,3) :: icov ! inverse convariance matrix
      END TYPE

      CONTAINS

         ! A few utility functions

         SUBROUTINE separation(cell_h, cell_ih, ri, rj, r)
            ! Calculates the distance between two position vectors (with PBC).
            !
            ! Note that minimum image convention is used, so only the image of
            ! atom j that is the shortest distance from atom i is considered.
            !
            ! Also note that while this may not work if the simulation
            ! box is highly skewed from orthorhombic, as
            ! in this case it is possible to return a distance less than the
            ! nearest neighbour distance. However, this will not be of
            ! importance unless the cut-off radius is more than half the
            ! width of the shortest face-face distance of the simulation box,
            ! which should never be the case.
            !
            ! Args:
            !    cell_h: The simulation box cell vector matrix.
            !    cell_ih: The inverse of the simulation box cell vector matrix.
            !    ri: The position vector of atom i.
            !    rj: The position vector of atom j
            !    r: The distance between the atoms i and j.

            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ri
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rj
            DOUBLE PRECISION, INTENT(OUT) :: r

            INTEGER k
            ! The separation in a basis where the simulation box
            ! is a unit cube.
            DOUBLE PRECISION, DIMENSION(3) :: sij
            DOUBLE PRECISION, DIMENSION(3) :: rij

            sij = matmul(cell_ih, ri-rj)
            DO k = 1, 3
               ! Finds the smallest separation of all the images of atom i and j
               sij(k) = sij(k) - dnint(sij(k)) ! Minimum Image Convention
            ENDDO
            ! TO IMPROVE WRITING THE FORMULAS EXPLICITLY
            rij = matmul(cell_h, sij)
            r = dsqrt(dot_product(rij, rij))

         END SUBROUTINE
         
         ! routine needed to invert a matrix.   
         SUBROUTINE inv3x3(M,IM)
            ! Invert a 3x3 matrix
            ! 
            ! Args:
            !    M: The 3x3 symmetric matrix to invert
            !    IM: The inverse matrix of M
            
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: M
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: IM
            
            ! cofactor maatrix
            DOUBLE PRECISION, DIMENSION(3,3) :: cof
            DOUBLE PRECISION d
            INTEGER i,j
            
            cof=0.0d0
            IM=0.0d0
            
            cof(1,1)=(M(2,2)*M(3,3)-M(2,3)*M(3,2))
            cof(1,2)=(M(1,3)*M(3,2)-M(1,2)*M(3,3))
            cof(1,3)=(M(1,2)*M(2,3)-M(2,2)*M(1,3))
            cof(2,1)=(M(2,3)*M(3,1)-M(2,1)*M(3,3))
            cof(2,2)=(M(1,1)*M(3,3)-M(1,3)*M(3,1))
            cof(2,3)=(M(1,3)*M(2,1)-M(1,1)*M(2,3))
            cof(3,1)=(M(2,1)*M(3,2)-M(2,2)*M(3,1))
            cof(3,2)=(M(1,2)*M(3,1)-M(1,1)*M(3,2))
            cof(3,3)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))
            
            ! claculate the determinant
            d=M(1,1)*cof(1,1)+M(1,2)*cof(2,1)+M(1,3)*cof(3,1)
            
            DO i=1,3
               DO j=1,3
                  IM(i,j)=cof(i,j)/d
               ENDDO
            ENDDO

         END SUBROUTINE inv3x3

         SUBROUTINE gauss_prepare(gpars)
            ! Initialize all the parameters of the gaussian
            ! 
            ! Args:
            !    gpars: gauss_type variable to initialize
             
            TYPE(gauss_type), INTENT(INOUT) :: gpars
            
            ! calculate the determinant of the covariance matrix
            gpars%det = gpars%cov(1,1)*(gpars%cov(2,2)*gpars%cov(3,3)-gpars%cov(3,2)*gpars%cov(2,3)) - &
                        gpars%cov(1,2)*(gpars%cov(2,1)*gpars%cov(3,3)-gpars%cov(2,3)*gpars%cov(3,1)) + &
                        gpars%cov(1,3)*(gpars%cov(2,1)*gpars%cov(3,2)-gpars%cov(2,2)*gpars%cov(3,1))
            
            ! calculate the inverse of the convariance matrix      
            CALL inv3x3(gpars%cov,gpars%icov)
            
            ! calculate the  logarithm of the normalization factor
            gpars%lnorm = dlog(1.0d0/dsqrt((dpigreco**3)*gpars%det))
         END SUBROUTINE gauss_prepare

         DOUBLE PRECISION FUNCTION gauss_logeval(gpars, x)
            ! Return the logarithm of the multivariate gaussian density
            ! 
            ! Args:
            !    gpars: gaussian parameters
            !    x: point in wich estimate the log of the gaussian 
            
            TYPE(gauss_type), INTENT(IN) :: gpars
            DOUBLE PRECISION, INTENT(IN) :: x(3)
            DOUBLE PRECISION dv(3),tmpv(3)
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
            DOUBLE PRECISION, INTENT(IN) :: x(3)
            
            gauss_eval = dexp(gauss_logeval(gpars,x))
         END FUNCTION gauss_eval


         ! the main routine for evaluating cluster weights
         SUBROUTINE hbmixture_GetGMMP(natoms,cell,icell,alpha,wcutoff,positions, &
                                      masktypes,nk,clusters,pks,sph,spd,spa)
            ! Return for each atoms the sum of the ..
            !
            ! Args:
            !    natoms: The number of atoms in the system.
            !    cell_h: The simulation box cell vector matrix.
            !    cell_ih: The inverse of the simulation box cell vector matrix.
            !    alpha: The smoothing factor
            !    wcutoff: The cutoff in w
            !    positions: The array containing the atoms coordiantes
            !    masktypes: The containing the atoms type
            !    nk: The number of gaussian from gaussian mixture model
            !    clusters: The array containing the structures with the gaussians parameters
            !    pks: The array containing the gaussians Pk
            !    sph:
            !    spd:
            !    spa:

            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: icell
            DOUBLE PRECISION, INTENT(IN) :: alpha
            DOUBLE PRECISION, INTENT(IN) :: wcutoff
            DOUBLE PRECISION, DIMENSION(3,natoms), INTENT(IN) :: positions
            INTEGER, DIMENSION(natoms), INTENT(IN) :: masktypes
            INTEGER, INTENT(IN) :: nk
            TYPE(gauss_type), DIMENSION(nk), INTENT(IN) :: clusters
            DOUBLE PRECISION, DIMENSION(nk), INTENT(IN) :: pks
            DOUBLE PRECISION, DIMENSION(nk,natoms), INTENT(OUT) :: sph, spa, spd

            DOUBLE PRECISION, DIMENSION(3) :: vwd
            DOUBLE PRECISION, DIMENSION(nk) :: pnk
            DOUBLE PRECISION pnormpk
            INTEGER ih,ia,id,k
            DOUBLE PRECISION rah, rdh

            ! initialize to zero the result vectors
            spa=0.0d0
            spd=0.0d0
            sph=0.0d0

            DO ih=1,natoms ! loop over H
               IF (IAND(masktypes(ih),TYPE_H).EQ.0) CYCLE ! test if it is an hydrogen
               DO id=1,natoms ! loop over D
                  ! test if it is a donor
                  IF (IAND(masktypes(id),TYPE_DONOR).EQ.0 .OR. ih.EQ.id) CYCLE
                  ! calculate the D-H distance
                  CALL separation(cell,icell,positions(:,ih),positions(:,id),rdh)
                  IF(rdh .gt. wcutoff) CYCLE  ! if the D-H distance is greater than the cutoff,
                                              ! we can already discard the D-H pair
                  DO ia=1,natoms ! loop over A
                     ! test if it is an acceptor
                     IF (IAND(masktypes(ia),TYPE_ACCEPTOR).EQ.0 &
                         .OR. (ia.EQ.id).OR.(ia.EQ.ih)) CYCLE
                     ! calculate the A-H distance
                     CALL separation(cell,icell,positions(:,ih),positions(:,ia),rah)
                     ! calculate w
                     vwd(2)=rah+rdh
                     IF(vwd(2).GT.wcutoff) CYCLE
                     ! calculate the PTC, v
                     vwd(1)=rdh-rah
                     ! calculate the A-D distance
                     CALL separation(cell,icell,positions(:,id),positions(:,ia),vwd(3))

                     pnk=0.0d0
                     pnormpk=0.0d0 ! normalization factor (mixture weight)

                     DO k=1,Nk
                        ! calculate the k probability for the point (v,w,rad)
                        ! and apply a smoothing elvating to alpha
                        ! and take care about the 1/J weight!
                        pnk(k) = ( gauss_eval(clusters(k), vwd)*pks(k)/ &
                                   ((vwd(2)*vwd(2)-vwd(1)*vwd(1))*vwd(3)) )**alpha
                        ! calculate the mixture weight
                        pnormpk = pnormpk+pnk(k)
                     ENDDO
                     IF (pnormpk.eq.0.0d0) CYCLE   ! skip cases in which the probability is tooooo tiny
                     pnk = pnk/pnormpk ! normalization

                     sph(:,ih) = sph(:,ih) + pnk(:)
                     spa(:,ia) = spa(:,ia) + pnk(:)
                     spd(:,id) = spd(:,id) + pnk(:)
                  ENDDO
               ENDDO
            ENDDO


         END SUBROUTINE hbmixture_GetGMMP

      END MODULE hbmixture
