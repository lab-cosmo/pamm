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
         USE mixture
      IMPLICIT NONE

      ! Constant to convert from bohrradius to angstrom
      DOUBLE PRECISION, PARAMETER :: bohr=0.529177219217

      CONTAINS

         SUBROUTINE xyz_read(ufile,natoms,header,labels,positions,endf)
            INTEGER, INTENT(OUT) :: natoms
            DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)  :: positions
            CHARACTER*4, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: labels
            CHARACTER*1024, INTENT(OUT) :: header
            INTEGER, INTENT(IN) :: ufile
            INTEGER, INTENT(OUT) :: endf
            INTEGER i

            READ(ufile,*,IOSTAT=endf) natoms
            IF(endf>0) error STOP "*** Error occurred while reading file. ***"
            IF(endf<0) return

            ! Get the snapshot header
            READ(ufile,'(A)',IOSTAT=endf) header
		    IF(endf>0) error STOP "*** Error occurred while reading file. ***"
		    IF(endf<0) return

            IF ( ALLOCATED(labels) .and. SIZE(labels)/=natoms ) DEALLOCATE(labels) ! atom number changed
            IF (.not.(ALLOCATED(labels))) ALLOCATE(labels(natoms))

            IF ( ALLOCATED(positions) .and. SIZE(positions)/=natoms ) DEALLOCATE(positions) ! atom number changed
            IF (.not.(ALLOCATED(positions))) ALLOCATE(positions(3,natoms))

            DO i=1,natoms ! label, x, y, z
                READ(ufile,*,IOSTAT=endf) labels(i),positions(1,i),positions(2,i),positions(3,i)
                IF(endf>0) error STOP "*** Error occurred while reading file. ***"
                IF(endf<0) return
            END DO
            endf=0
         END SUBROUTINE xyz_read

         SUBROUTINE xyz_write(ufile,natoms,header,labels,positions)
            INTEGER, INTENT(IN) :: ufile, natoms

            CHARACTER*1024, INTENT(IN) :: header
            CHARACTER*4, DIMENSION(natoms), INTENT(IN) :: labels
            DOUBLE PRECISION, DIMENSION(3,natoms), INTENT(IN) :: positions
            INTEGER i

            ! header
            WRITE(ufile,"(I4)") natoms
            WRITE(ufile,"(A)") header

            ! body
            DO i=1,natoms
               WRITE(ufile,"(A2,A1,F11.7,A1,F11.7,A1,F11.7)") &
                    trim(labels(i)), " ", positions(1,i), " ", &
                    positions(2,i), " ", positions(3,i)
            ENDDO
         END SUBROUTINE xyz_write
      END MODULE xyz
