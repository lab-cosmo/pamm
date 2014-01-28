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
                                      masktypes,nk,clusters,probabilities)
            ! Calculate the probabilities 
            ! ...
            ! Args:
            !    param: descript 
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, INTENT(IN) :: lbox
            DOUBLE PRECISION, INTENT(IN) :: alpha
            DOUBLE PRECISION, INTENT(IN) :: wcutoff
            DOUBLE PRECISION, DIMENSION(natoms), INTENT(IN) :: positions
            INTEGER, DIMENSION(natoms), INTENT(IN) :: masktypes
            INTEGER, INTENT(IN) :: nk
            TYPE(GAUSS_TYPE), DIMENSION(nk), INTENT(IN) :: clusters
            DOUBLE PRECISION, DIMENSION(natoms,nk), INTENT(OUT) :: probabilities
            
            DOUBLE PRECISION, DIMENSION(2) :: vw
            DOUBLE PRECISION, DIMENSION(nk) :: pnk
            DOUBLE PRECISION pnormpk
            INTEGER i,j,m,k
            
            ! initialize
            probabilities=probabilities*0.0d0
            DO i=1,natoms ! loop over H
               IF (IAND(masktypes(i),TYPE_H).EQ.0) CYCLE
               DO j=1,natoms
                  IF (IAND(masktypes(j),TYPE_DONOR).EQ.0) CYCLE
                  IF (i.EQ.j) CYCLE
                  DO m=1,natoms
                     IF (IAND(masktypes(m),TYPE_ACCEPTOR).EQ.0) CYCLE
                     IF ((m.EQ.i).OR.(m.EQ.j)) CYCLE
                     ! calculate v and w
                     CALL hbmixture_Getvw(lbox, positions(i), &
                                          positions(j), positions(m), vw)
                     IF(vw(2).GT.wcutoff) CYCLE
                     pnk=pnk*0.0d0
                     pnormpk=0.0d0 ! Normalization factor
                     DO k=1,Nk
                        pnk(k) = gauss_eval(clusters(k), vw)**alpha
                        pnormpk = pnormpk+pnk(k)
                     ENDDO
                     ! Normalize
                     pnk = pnk/pnormpk
                     DO k=1,Nk
                        probabilities(i,k)=probabilities(i,k)+pnk(k)! P hydrongen
                        probabilities(j,k)=probabilities(j,k)+pnk(k)! donor
                        probabilities(m,k)=probabilities(m,k)+pnk(k)! acceptor
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

         END SUBROUTINE hbmixture_GetGMMP

      END MODULE hbmixture
