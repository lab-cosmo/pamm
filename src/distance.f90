! This contains the algorithms needed to calculate the distance between atoms.
!
! Copyright (C) 2014, Piero Gasparotto and Michele Ceriotti
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!
! Functions:
!    vector_separation: Calculates the vector separating two atoms.
!    separation: Calculates the square distance between two vectors.

      MODULE distance
      IMPLICIT NONE

      CONTAINS

         SUBROUTINE vector_separation(cell_h, cell_ih, ri, rj, rij, r2)
            ! Calculates the vector separating two atoms.
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
            !    rij: The vector separating atoms i and j.
            !    r2: The square of the distance between atoms i and j.

            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ri
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rj
            DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: rij
            DOUBLE PRECISION, INTENT(OUT) :: r2

            INTEGER k
            DOUBLE PRECISION, DIMENSION(3) :: sij
            ! The separation in a basis where the simulation box
            ! is a unit cube.

            sij = matmul(cell_ih, ri - rj)
            DO k = 1, 3
               ! Finds the smallest separation of all the images of atom i and j
               sij(k) = sij(k) - dnint(sij(k))
            ENDDO
            rij = matmul(cell_h, sij)
            r2 = dot_product(rij,rij)

         END SUBROUTINE

         SUBROUTINE separation(cell_h, cell_ih, ri, rj, r2)
            ! Calculates the squared distance between two position vectors.
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
            !    r2: The square of the distance between atoms i and j.

            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_h
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell_ih
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ri
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rj
            DOUBLE PRECISION, INTENT(OUT) :: r2

            INTEGER k
            ! The separation in a basis where the simulation box
            ! is a unit cube.
            DOUBLE PRECISION, DIMENSION(3) :: sij
            DOUBLE PRECISION, DIMENSION(3) :: rij

            sij = matmul(cell_ih, ri - rj)
            DO k = 1, 3
               ! Finds the smallest separation of all the images of atom i and j
               sij(k) = sij(k) - dnint(sij(k))
            ENDDO
            rij = matmul(cell_h, sij)
            r2 = dot_product(rij, rij)

         END SUBROUTINE

         SUBROUTINE separation_cubic(L,ri,rj,r)
            ! Cubic box version.
            ! Calculates the squared distance between two position vectors.
            !
            ! Note that minimum image convention is used, so only the image of
            ! atom j that is the shortest distance from atom i is considered.
            !
            ! Args:
            !    L: The box lenght
            !    ri: The position vector of atom i.
            !    rj: The position vector of atom j
            !    rij: The vector separating atoms i and j.
            !    r: The distance between atoms i and j.

            DOUBLE PRECISION, INTENT(IN) :: L
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ri
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rj
            DOUBLE PRECISION, INTENT(OUT) :: r

            DOUBLE PRECISION, DIMENSION(3) :: rij
            INTEGER i
            rij = ri-rj
            DO i=1,3
               ! MIC
               rij(i) = rij(i)-dnint(rij(i)/L)*L
            END DO
            r = dsqrt(dot_product(rij,rij))
         END SUBROUTINE

         SUBROUTINE vector_separation_cubic(L,ri,rj,rij,r)
            ! Cubic box version.
            ! Calculates the squared distance between two position vectors.
            !
            ! Note that minimum image convention is used, so only the image of
            ! atom j that is the shortest distance from atom i is considered.
            !
            ! Args:
            !    L: The box lenght
            !    ri: The position vector of atom i.
            !    rj: The position vector of atom j
            !    rij: The vector separating atoms i and j.
            !    r: The distance between atoms i and j.

            DOUBLE PRECISION, INTENT(IN) :: L
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ri
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rj
            DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: rij
            DOUBLE PRECISION, INTENT(OUT) :: r

            INTEGER i
            rij = ri-rj
            DO i=1,3
               ! MIC
               rij(i) = rij(i)-dnint(rij(i)/L)*L
            END DO
            r = dsqrt(dot_product(rij,rij))
         END SUBROUTINE

      END MODULE distance
