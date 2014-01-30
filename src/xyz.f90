! This file contains some routines to work with XYZ files
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
!
! This contains the functions that calculate the potential, forces and
! virial tensor of a single-component LJ system.
! Includes functions which calculate the long-range correction terms for a
! simulation with a sharp nearest-neighbour cut-off.
!
! Functions:
!    XYZ_functions: Description ...

      MODULE xyz
      IMPLICIT NONE

      CONTAINS
      
         SUBROUTINE xyz_GetInfo(filename,NAtoms,cell,NSteps) !! TO TEST (because of wc -l)!
            ! Get the number of atoms, the box lenght and the number of steps.
            !
            ! Args:
            !    filename: The filename.
            !    NAtoms: The number of atoms.
            !    Lbox: The box lenght
            !    NSteps: The number of steps.

            CHARACTER*70, INTENT(IN) :: filename
            INTEGER, INTENT(OUT) :: NAtoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: cell
            INTEGER, INTENT(OUT) :: NSteps

            CHARACTER*1024 :: cmdbuffer
            CHARACTER*30 dummy1,dummy2


            INTEGER Nlines

            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) NAtoms
            READ(11,*) dummy1,dummy2,cell(1,1),cell(2,2),cell(3,3)
            CLOSE(UNIT=11)
            ! we assume to be in a unix system and instead of reading all the file
            ! we use wc -l (to test if it is faster!!!!!!!!)
            ! I used wc because normally we ar using really big file..
            ! get the line number and save to a temp file
            cmdbuffer='wc -l '//filename//'>> tmp.tmp'
            CALL system(cmdbuffer)
            OPEN(UNIT=11,FILE='tmp.tmp',STATUS='OLD',ACTION='READ')
            ! read the line numbers
            READ(11,*) Nlines
            CLOSE(UNIT=11)
            !remove the temp file
            cmdbuffer="rm tmp.tmp"
            CALL system(cmdbuffer)
            NSteps = Nlines/(NAtoms+2)

         END SUBROUTINE xyz_GetInfo

         SUBROUTINE xyz_GetLabels(filename,labels)
            ! Get the line numbers of the file and return the number
            ! of atoms and the number of steps.
            !
            ! Args:
            !    filename: The filename.
            !    labels: The vector containing the atoms' label

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

         SUBROUTINE xyz_GetSnap(mode,ufile,natoms,positions)
            ! Get the positions for the chosen step
            !
            ! Args:
            INTEGER, INTENT(IN) :: mode
            INTEGER, INTENT(IN) :: ufile
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(3,NAtoms), INTENT(OUT)  :: positions

            CHARACTER*1024 dummy
            INTEGER i, ierr
            IF (mode.eq.0) THEN
               ! Discard snap
               DO i=1,NAtoms+2
                  READ(ufile,*) dummy
               END DO
            ELSE
               ! Get the snap
               ! skip the header
               READ(ufile,*) dummy ! NAtoms
               READ(ufile,*) dummy ! Cell info
               DO i=1,natoms
                  READ(ufile,*) dummy,positions(1,i),positions(2,i),positions(3,i)
               END DO
            ENDIF

         END SUBROUTINE
         
      END MODULE xyz
