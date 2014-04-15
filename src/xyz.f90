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
      DOUBLE PRECISION, PARAMETER :: bohr=0.5291772192171717171717171717171717
      
      CONTAINS

         SUBROUTINE xyz_read(skipna,mode,nptmode,convert,ufile,natoms,positions,labels,cell,icell,endf)
            ! Get the coordinates of all the atoms
            !
            ! Args:
            !    mode: The flag that tell if we have to skip or not the snapshot
            !    ufile: ID of the input file
            !    natoms: The number of atoms in the system.
            !    positions: The array containing the atoms coordinates.
            
            LOGICAL, INTENT(IN) :: skipna ! to skip the reading of the line
                                          ! containing the number of atoms.
                                          ! I need this just for the first snapshot
            INTEGER, INTENT(IN) :: mode 
            LOGICAL, INTENT(IN) :: nptmode
            LOGICAL, INTENT(IN) :: convert
            INTEGER, INTENT(IN) :: ufile
            INTEGER, INTENT(OUT) :: natoms
            DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)  :: positions
            CHARACTER*4, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: labels
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(INOUT) :: cell
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(INOUT) :: icell
            INTEGER, INTENT(OUT) :: endf

            CHARACTER*30 dummy1,dummy2
            INTEGER i,nc

            IF (mode.eq.0) THEN
               ! Discard this timesnapshot
               nc=natoms+2
               IF(skipna) nc=nc-1
               DO i=1,nc
                  READ(ufile,*,IOSTAT=endf) dummy1
                  IF(endf>0) error STOP "*** Error occurred while reading file. ***"
                  IF(endf<0) return
               END DO
            ELSE
               ! Get the snapshot
            
               ! Get the atom number
               IF(.NOT.skipna)READ(ufile,*,IOSTAT=endf) natoms 
               IF(endf>0) error STOP "*** Error occurred while reading file. ***"
			   IF(endf<0) return
			   
		       ! Get the cell
		       IF ((cell(1,1)==(0.0d0) .OR. cell(2,2)==(0.0d0) .OR. cell(3,3)==(0.0d0)).or.(nptmode)) THEN
		          ! read the box size
		          READ(ufile,*,IOSTAT=endf) dummy1,dummy2,cell(1,1),cell(2,2),cell(3,3)
		          IF(endf>0) error STOP "*** Error occurred while reading file. ***"
		          IF(endf<0) return
		          ! convert the box if necessary
		          IF(convert) cell=cell*bohr
		          CALL invmatrix(3,cell,icell) 
		       ELSE
		          READ(ufile, '(A)') dummy1
               ENDIF

               IF (.not.(ALLOCATED(labels))) ALLOCATE(labels(natoms))
               IF (.not.(ALLOCATED(positions))) ALLOCATE(positions(3,natoms))
               
               DO i=1,natoms ! label, x, y, z
                  READ(ufile,*,IOSTAT=endf) labels(i),positions(1,i),positions(2,i),positions(3,i)
                  IF(endf>0) error STOP "*** Error occurred while reading file. ***"
                  IF(endf<0) return
               END DO
               ! convert if necessary
		       IF(convert) positions=positions*bohr
            ENDIF
            endf=0
         END SUBROUTINE xyz_read
         
         SUBROUTINE xyz_write(ufile,natoms,cell,ts,labels,positions,sh,sd,sa)
            INTEGER, INTENT(IN) :: ufile
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell
            INTEGER, INTENT(IN) :: ts
            CHARACTER*4, DIMENSION(natoms), INTENT(INOUT) :: labels
            DOUBLE PRECISION, DIMENSION(3,natoms), INTENT(IN) :: positions
            DOUBLE PRECISION, DIMENSION(natoms), INTENT(IN) :: sh, sd, sa
            
            INTEGER i

            ! header
            WRITE(ufile,"(I4)") natoms
            WRITE(ufile,"(a,F11.7,a,F11.7,a,F11.7,a,I10)") &
                 "# CELL(abc): ",cell(1,1)," ",cell(2,2)," ",cell(3,3)," Step: ",ts
            ! body
            DO i=1,natoms
               WRITE(ufile,"(A2,A1,F11.7,A1,F11.7,A1,F11.7)", ADVANCE='NO') &
                    trim(labels(i)), " ", positions(1,i), " ", &
                    positions(2,i), " ", positions(3,i)
               WRITE(ufile,"(3(A1,ES21.8E4))", ADVANCE='NO') " ", sh(i), &
                       " ", sd(i), " ", sa(i)
               WRITE(ufile,*)
            ENDDO
         END SUBROUTINE xyz_write

      END MODULE xyz
