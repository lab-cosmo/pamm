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
!    xyz_GetLabels: Get the atoms label
!    xyz_GetSnap: Get the atoms coordinates

      MODULE xyz
      IMPLICIT NONE

      CONTAINS
      
         SUBROUTINE xyz_GetInfo(filename,natoms,cell,nsteps)
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
            INTEGER, INTENT(OUT) :: nsteps

            CHARACTER*1024 :: cmdbuffer
            CHARACTER*30 dummy1,dummy2


            INTEGER nlines

            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) natoms
            ! we assume an orthorombic box
            READ(11,*) dummy1,dummy2,cell(1,1),cell(2,2),cell(3,3)
            CLOSE(UNIT=11)
            ! we assume to be in a unix system and instead of reading all the file
            ! we use wc -l (to test if it is faster!!!!!!!!)
            ! I used wc because normally we ar using really big file..
            
            ! get the line number and save to a temp file
            cmdbuffer='wc -l '//filename//'>> tmp.tmp'
            ! execute the bash command
            CALL system(cmdbuffer)
            OPEN(UNIT=11,FILE='tmp.tmp',STATUS='OLD',ACTION='READ')
            ! read the line numbers
            READ(11,*) nlines
            CLOSE(UNIT=11)
            ! remove the temp file
            cmdbuffer="rm tmp.tmp"
            CALL system(cmdbuffer)
            nsteps = nlines/(natoms+2)
         END SUBROUTINE xyz_GetInfo

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

         SUBROUTINE xyz_GetSnap(mode,ufile,natoms,positions)
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

            CHARACTER*1024 dummy
            INTEGER i
            IF (mode.eq.0) THEN
               ! Discard this timesnapshot
               DO i=1,NAtoms+2
                  READ(ufile,*) dummy
               END DO
            ELSE
               ! Get the snapshot
               ! skip the header
               READ(ufile,*) dummy ! NAtoms
               READ(ufile,*) dummy ! Cell info
               DO i=1,natoms ! label, x, y, z
                  READ(ufile,*) dummy,positions(1,i),positions(2,i),positions(3,i)
               END DO
            ENDIF

         END SUBROUTINE
         
      END MODULE xyz
