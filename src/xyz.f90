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
!    xyz_GetInfo: Get the general info from the file (natoms,cell,nsteps)
!    xyz_AdjustCELL: Get the cell and the inverse cell
!    xyz_GetLabels: Get the atoms label
!    xyz_GetSnap: Get the atoms coordinates

      MODULE xyz
         USE hbmixture
      IMPLICIT NONE

      CONTAINS
      
         SUBROUTINE xyz_GetInfo(filename,natoms,cell)
            ! Get the number of atoms, the box lenght and the number of steps.
            !
            ! Args:
            !    filename: The filename.
            !    natoms: The number of atoms in the system.
            !    Lbox: The box lenght
            !    NSteps: The number of steps.

            CHARACTER*70, INTENT(IN) :: filename
            INTEGER, INTENT(OUT) :: natoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: cell

            CHARACTER*1024 :: cmdbuffer
            CHARACTER*30 dummy1,dummy2


            INTEGER nlines

            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) natoms
            ! we assume an orthorombic box
            READ(11,*) dummy1,dummy2,cell(1,1),cell(2,2),cell(3,3)
            CLOSE(UNIT=11)
         END SUBROUTINE xyz_GetInfo
         
         SUBROUTINE xyz_AdjustCELL(ufile,cell,icell)
            ! Get the number of atoms, the box lenght and the number of steps.
            !
            ! Args:
            !    ufile: ID of the input file
            !    cell: The simulation box cell vector matrix.
            !    icell: The inverse of the simulation box cell vector matrix.
            
            INTEGER, INTENT(IN) :: ufile
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: cell
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: icell

            CHARACTER*30 dummy1,dummy2
            INTEGER pos,ierr
			
			pos=FTELL(ufile)
			
			cell=0.0d0
            ! discard the first line (number of atoms)
            READ(ufile,*) dummy1
            ! we assume an orthorombic box
            READ(ufile,*) dummy1,dummy2,cell(1,1),cell(2,2),cell(3,3)
            ! get the inverse of the cell
            CALL inv3x3(cell,icell)
            CALL FSEEK(ufile, pos, 0, ierr)
         END SUBROUTINE xyz_AdjustCELL

         SUBROUTINE xyz_GetLabels(filename,labels)
            ! Get the atoms label
            !
            ! Args:
            !    filename: The filename.
            !    labels: The vector containing the atoms label

            CHARACTER*70, INTENT(IN) :: filename
            CHARACTER*4,ALLOCATABLE, DIMENSION(:),INTENT(OUT) :: labels

            CHARACTER*1024 :: cmdbuffer
            INTEGER natoms,i

            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) natoms
            ALLOCATE(labels(natoms))
            READ(11,*) cmdbuffer
            DO i=1,natoms
               READ(11,*) cmdbuffer
               labels(i)=trim(cmdbuffer)
            ENDDO
            CLOSE(UNIT=11)

         END SUBROUTINE

         SUBROUTINE xyz_GetSnap(mode,ufile,natoms,positions,endf)
            ! Get the coordinates of all the atoms
            !
            ! Args:
            !    mode: The flag that tell if we have to skip or not the snapshot
            !    ufile: ID of the input file
            !    natoms: The number of atoms in the system.
            !    positions: The array containing the atoms coordinates.
            
            INTEGER, INTENT(IN) :: mode 
            INTEGER, INTENT(IN) :: ufile
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(3,NAtoms), INTENT(OUT)  :: positions
            LOGICAL, INTENT(OUT) :: endf

            CHARACTER*1024 dummy
            INTEGER i,ierr
            
            endf = .false.
            
            IF (mode.eq.0) THEN
               ! Discard this timesnapshot
               DO i=1,NAtoms+2
                  READ(ufile,*,IOSTAT=ierr) dummy
                  IF(ierr<0) endf=.true.
                  IF(endf) cycle
               END DO
            ELSE
               ! Get the snapshot
               ! skip the header
               READ(ufile,*,IOSTAT=ierr) dummy ! NAtoms
               IF(ierr<0) endf=.true.
               READ(ufile,*,IOSTAT=ierr) dummy ! Cell info
               IF(ierr<0) endf=.true.
               DO i=1,natoms ! label, x, y, z
                  READ(ufile,*,IOSTAT=ierr) dummy,positions(1,i),positions(2,i),positions(3,i)
                  IF(ierr<0) endf=.true.
                  IF(endf) EXIT
               END DO
            ENDIF

         END SUBROUTINE
         
      END MODULE xyz
