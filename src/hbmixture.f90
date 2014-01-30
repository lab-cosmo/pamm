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

      MODULE hbmixture
         USE distance
         USE gaussian
      IMPLICIT NONE

      ! Types used by bitwise operators
      ! the must be power of 2
      INTEGER, PARAMETER :: TYPE_NONE=0
      INTEGER, PARAMETER :: TYPE_H=1
      INTEGER, PARAMETER :: TYPE_DONOR=2
      INTEGER, PARAMETER :: TYPE_ACCEPTOR=4

      CONTAINS

         ! TODO : general cell
         ! SUBROUTINE hbmixture_Getvw(cell_h, cell_ih, rH, rD, rA, vw)
         SUBROUTINE hbmixture_Getvw(L, rH, rD, rA, vw)
            ! Calculates the proton transfer coordinates
            !DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            !DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            ! BOX lenght
            DOUBLE PRECISION, INTENT(IN) :: L
            ! position vectors
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rH ! Hydrogen
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rD ! Donor (D)
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rA ! Acceptor (A)
            DOUBLE PRECISION, DIMENSION(2), INTENT(OUT) :: vw

            ! separation D-H , A-H
            DOUBLE PRECISION rDH,rAH

            CALL separation_cubic(L,rD,rH,rDH)
            CALL separation_cubic(L,rA,rH,rAH)

            ! PTC
            vw(1) = rDH - rAH
            vw(2) = rDH + rAH

         END SUBROUTINE hbmixture_Getvw

         SUBROUTINE hbmixture_GetGMMP(natoms,lbox,alpha,wcutoff,positions, &
                                      masktypes,nk, clusters, sph, spd, spa)
            ! Calculate the probabilities
            ! ...
            ! Args:
            !    param: descript
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, INTENT(IN) :: lbox
            DOUBLE PRECISION, INTENT(IN) :: alpha
            DOUBLE PRECISION, INTENT(IN) :: wcutoff
            DOUBLE PRECISION, DIMENSION(natoms,3), INTENT(IN) :: positions
            INTEGER, DIMENSION(natoms), INTENT(IN) :: masktypes
            INTEGER, INTENT(IN) :: nk
            TYPE(GAUSS_TYPE), DIMENSION(nk), INTENT(IN) :: clusters
            DOUBLE PRECISION, DIMENSION(natoms,nk), INTENT(OUT) :: sph, spa, spd

            DOUBLE PRECISION, DIMENSION(2) :: vw
            DOUBLE PRECISION, DIMENSION(nk) :: pnk
            DOUBLE PRECISION pnormpk
            INTEGER ih,ia,id,k
            DOUBLE PRECISION rah, rdh

            ! initialize
            spa=0.0d0
            spd=0.0d0
            sph=0.0d0
            DO ih=1,natoms ! loop over H
               IF (IAND(masktypes(ih),TYPE_H).EQ.0) CYCLE
               DO id=1,natoms
                  IF (IAND(masktypes(id),TYPE_DONOR).EQ.0 .OR. ih.EQ.id) CYCLE

                  CALL separation_cubic(lbox,positions(ih,:),positions(id,:),rdh)
                  IF(rdh .gt. wcutoff) CYCLE  ! if one of the distances is greater than the cutoff, we can already discard the D-H pair

                  DO ia=1,natoms
                     IF (IAND(masktypes(ia),TYPE_ACCEPTOR).EQ.0 &
                         .OR. (ia.EQ.id).OR.(ia.EQ.ih)) CYCLE

                     CALL separation_cubic(lbox,positions(ih,:),positions(ia,:),rah)
                     vw(2)=rah+rdh
                     IF(vw(2).GT.wcutoff) CYCLE
                     vw(1)=rdh-rah
                     !write(*,*)vw(1)," ",vw(2)
                     
                     pnk=pnk*0.0d0
                     pnormpk=0.0d0 ! Normalization factor
                     DO k=1,Nk
                        pnk(k) = gauss_eval(clusters(k), vw)**alpha
                        pnormpk = pnormpk+pnk(k)
                     ENDDO
                     ! Normalize
                     pnk = pnk/pnormpk
                     sph(ih,:) = sph(ih,:) + pnk(:)
                     spa(ia,:) = spa(ia,:) + pnk(:)
                     spd(id,:) = spd(id,:) + pnk(:)
                     !write(*,*) vw(1),vw(2),"//",ih,"",ia,"",id,"//",sph(ih,1),"",spa(ia,1),"",spd(id,1)
                  ENDDO
               ENDDO
            ENDDO

         END SUBROUTINE hbmixture_GetGMMP

      END MODULE hbmixture
