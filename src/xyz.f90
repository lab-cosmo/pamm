! This file contains the routines used to work with XYZ files
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
!    xyz_read: Get one frame off a xyz file

      MODULE xyz
        IMPLICIT NONE
! Constant to convert from bohrradius to angstrom
      DOUBLE PRECISION, PARAMETER :: bohr=0.529177219217
      CONTAINS
      SUBROUTINE xyz_read(ufile,natoms,header,labels,positions,endf)
         INTEGER, INTENT(OUT) :: natoms
         DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)  :: positions
         CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: labels
         CHARACTER(LEN=1024), INTENT(OUT) :: header
         INTEGER, INTENT(IN) :: ufile
         INTEGER, INTENT(OUT) :: endf
         INTEGER i
         
         READ(ufile,*,IOSTAT=endf) natoms
         IF(endf>0) STOP "*** Error occurred while reading file. ***"
         IF(endf<0 .or. natoms==0) return
         
         ! Get the snapshot header
         READ(ufile,'(A)',IOSTAT=endf) header
         IF(endf>0) STOP "*** Error occurred while reading file. ***"
         IF(endf<0) return
         IF ( ALLOCATED(labels) .and. SIZE(labels)/=natoms ) DEALLOCATE(labels) ! atom number changed
         IF (.not.(ALLOCATED(labels))) ALLOCATE(labels(natoms))
         IF ( ALLOCATED(positions) .and. SIZE(positions)/=natoms*3 ) DEALLOCATE(positions) ! atom number changed
         IF (.not.(ALLOCATED(positions))) ALLOCATE(positions(3,natoms))
         DO i=1,natoms ! label, x, y, z
            READ(ufile,*,IOSTAT=endf) labels(i),positions(1,i),positions(2,i),positions(3,i)
            IF(endf>0) STOP "*** Error occurred while reading file. ***"
            IF(endf<0) return
         END DO
         endf=0
      END SUBROUTINE xyz_read
		
      SUBROUTINE xyz_write(ufile,natoms,header,labels,positions)
         INTEGER, INTENT(IN) :: ufile, natoms
         CHARACTER(LEN=1024), INTENT(IN) :: header
         CHARACTER(LEN=4), DIMENSION(natoms), INTENT(IN) :: labels
         DOUBLE PRECISION, DIMENSION(3,natoms), INTENT(IN) :: positions
         INTEGER i
         ! header
         WRITE(ufile,"(I4)") natoms
         WRITE(ufile,"(A)") trim(header)
         ! body
         DO i=1,natoms
            WRITE(ufile,"(A2,A1,F11.7,A1,F11.7,A1,F11.7)") &
                  trim(labels(i)), " ", positions(1,i), " ", &
                  positions(2,i), " ", positions(3,i)
         ENDDO
      END SUBROUTINE xyz_write
		
      ! A few utility functions
      SUBROUTINE pbcdist(cell_h, cell_ih, ri, rj, r)
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
         rij = matmul(cell_h, sij)
         r = dsqrt(dot_product(rij, rij))
      END SUBROUTINE
      
      END MODULE xyz
		
		
