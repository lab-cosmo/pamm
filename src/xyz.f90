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
      
         SUBROUTINE xyz_GetInfo(filename,NAtoms,Lbox,NSteps) !! TO TEST!
            ! Get the number of atoms, the box lenght and the number of steps.
            !
            ! Args:
            !    filename: The filename.
            !    NAtoms: The number of atoms.
            !    Lbox: The box lenght
            !    NSteps: The number of steps.

            CHARACTER*70, INTENT(IN) :: filename
            INTEGER, INTENT(OUT) :: NAtoms
            DOUBLE PRECISION, INTENT(OUT) :: Lbox
            INTEGER, INTENT(OUT) :: NSteps

            CHARACTER*1024 :: cmdbuffer
            CHARACTER*30 dummy1,dummy2,dummy3
            
            
            INTEGER Nlines
            
            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) NAtoms	
            READ(11,*) dummy1,dummy2,dummy3,Lbox
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

         END SUBROUTINE

         SUBROUTINE xyz_GetNts(filename,NSteps) !! TO TEST!
            ! Get the line numbers of the file and return the number of steps.
            !
            ! Args:
            !    filename: The filename.
            !    NSteps: The number of steps.

            CHARACTER*70, INTENT(IN) :: filename
            INTEGER, INTENT(OUT) :: NSteps

            CHARACTER*1024 :: cmdbuffer
            INTEGER Nlines,NAtoms
            
            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) NAtoms	
            CLOSE(UNIT=11)
            ! we assume to be in a unix system and instead of reading all the file
            ! we use wc -l (to test if it is faster!!!!!!!!)
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

         END SUBROUTINE

         SUBROUTINE xyz_GetNAtoms(filename,NAtoms)
            ! Get the the number of atoms.
            !
            ! Args:
            !    filename: The filename.
            !    NAtoms: The number of atoms.

            CHARACTER*70, INTENT(IN) :: filename
            INTEGER, INTENT(OUT) :: NAtoms

            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) NAtoms	
            CLOSE(UNIT=11)      

         END SUBROUTINE
         
         SUBROUTINE xyz_GetDimBox(filename,lbox)
            ! Get the box lenght
            ! Simplest case : CUBIC BOX
            ! To improve for the general case
            !
            ! Args:
            !    filename: The filename.
            !    lbox: The box lenght

            CHARACTER*70, INTENT(IN) :: filename
            DOUBLE PRECISION, INTENT(OUT) :: lbox

            CHARACTER*30 dummy1,dummy2,dummy3

            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) dummy1
            READ(11,*) dummy1,dummy2,dummy3,lbox
            CLOSE(UNIT=11)

         END SUBROUTINE

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
            INTEGER NAtoms,i
            
            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename,STATUS='OLD',ACTION='READ')
            READ(11,*) NAtoms
            ALLOCATE(labels(NAtoms))
            READ(11,*) cmdbuffer
            DO i=1,NAtoms
               READ(11,*) cmdbuffer
               labels(i)=trim(cmdbuffer)
            ENDDO
            CLOSE(UNIT=11)

         END SUBROUTINE

         ! We won't never use this routine.
         !! TO CANCEL
         SUBROUTINE xyz_CountAtom(NAtoms,labels,id,num)
            ! Get the line numbers of id atoms in labels.
            !
            ! Args:
            !    NAtoms: The number of atoms.
            !    labels: The vector containing the atoms' label
            !    id: the atom type to serch for
            !    num: the numbers of id atoms in labels

            INTEGER, INTENT(IN) :: NAtoms
            CHARACTER*4, DIMENSION(NAtoms), INTENT(IN) :: labels
            CHARACTER*4, INTENT(IN) :: id
            INTEGER, INTENT(OUT) :: num

            INTEGER i

            num=0
            DO i=1,NAtoms
               IF (trim(labels(i)).eq.trim(id)) THEN
                  num=num+1
               ENDIF
            ENDDO
            
         END SUBROUTINE

         SUBROUTINE xyz_GetSnap(mode,filename,NAtoms,pos,newpos,positions)
            ! Get the positions for the chosen step
            !
            ! Args:
            INTEGER, INTENT(IN) :: mode
            CHARACTER*70, INTENT(IN) :: filename
            INTEGER, INTENT(IN) :: NAtoms
            INTEGER, INTENT(IN) :: pos
            INTEGER, INTENT(OUT) :: newpos
            DOUBLE PRECISION, DIMENSION(NAtoms,3), INTENT(OUT)  :: positions

            CHARACTER*1024 dummy
            INTEGER i, ierr
            ! open the file and read the atom number
            OPEN(UNIT=11,FILE=filename)
            ! go to the desired position
            CALL FSEEK(11, pos, 0, ierr)
            IF (mode.eq.0) THEN
               ! Discard snap
               DO i=1,NAtoms+2
                  READ(11,*) dummy
               END DO
            ELSE
               ! Get the snap
               ! skip the header
               READ(11,*) dummy ! NAtoms
               READ(11,*) dummy ! Cell info
               DO i=1,NAtoms
                  READ(11,*) dummy,positions(i,1),positions(i,2),positions(i,3)
                  !test what we read
                  !WRITE(*,*) trim(dummy),positions(i,1),positions(i,2),positions(i,3)
               END DO  
            ENDIF
            newpos = FTELL(11)
            CLOSE(UNIT=11)

         END SUBROUTINE

         SUBROUTINE xyz_ExtractAtoms(NAtoms,positions,labels,N,idatom,ext_coord)
            ! Extrat idatom type atoms
            !
            ! Args:
            INTEGER, INTENT(IN) :: NAtoms
            DOUBLE PRECISION, DIMENSION(NAtoms,3), INTENT(IN)  :: positions
            CHARACTER*4, DIMENSION(NAtoms) :: labels
            INTEGER, INTENT(IN) :: N
            CHARACTER, INTENT(IN) :: idatom
            DOUBLE PRECISION, DIMENSION(N,4), INTENT(OUT)  :: ext_coord

            INTEGER i,counter

            counter=0
            DO i=1,NAtoms
               IF(labels(i).eq.idatom)THEN
                  counter=counter+1
                  ext_coord(counter,1) = i
                  ext_coord(counter,2) = positions(i,1)
                  ext_coord(counter,3) = positions(i,2)
                  ext_coord(counter,4) = positions(i,3)
               ENDIF
            END DO 
         END SUBROUTINE

      END MODULE xyz
