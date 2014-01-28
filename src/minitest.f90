! The main program which runs our driver test case potentials
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

      PROGRAM minitest
         USE distance
         USE gaussian
         USE xyz
         USE hbmixture
      IMPLICIT NONE
      
      ! Conversion constant from bohrradius to angstrom 
      DOUBLE PRECISION, PARAMETER :: bohr=0.5291772192171717171717171717171717
      
      CHARACTER*70 :: filename, gaussianfile, outputfile
      CHARACTER*1024 :: cmdbuffer
      ! system parameters
      INTEGER natoms
      DOUBLE PRECISION lbox, alpha, wcutoff
      INTEGER nsteps, startstep, delta
      ! gaussians
      INTEGER Nk
      ! array of gaussians
      TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:) :: clusters
      ! for the parser
      INTEGER ccmd
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string
      ! for a faster reading
      ! counters
      INTEGER i,ts
      ! parmater for seeking in the input file
      INTEGER pos,newpos
      
      LOGICAL verbose,convert
      
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: positions
      ! mask to define what is what
      INTEGER, ALLOCATABLE, DIMENSION(:) :: masktypes
      CHARACTER*4, ALLOCATABLE, DIMENSION(:) :: labels
      CHARACTER*4, DIMENSION(4) :: vtacc,vtdon,vtH
      
      DOUBLE PRECISION dummyd1,dummyd2
      INTEGER dummyi1,dummyi2
      
      !! default values
      DO i=1,4
         vtacc(i)="NULL"
         vtdon(i)="NULL"
         vtH(i)="NULL"
      ENDDO
      filename="NULL"
      gaussianfile="NULL"
      outputfile="out-Pad.dat"
      Nk=-1
      ccmd=0
      lbox=-1.0d0
      nsteps=-1
      natoms=-1
      delta=1
      startstep=1
      convert = .false.
      verbose = .false.
      !!
      
      !!!!! PARSER
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-i") THEN ! input xyz
            ccmd = 1
         ELSEIF (cmdbuffer == "-gn") THEN ! number of gaussians
            ccmd = 2
         ELSEIF (cmdbuffer == "-gf") THEN ! file containing gaussian parmeters
            ccmd = 3
         ELSEIF (cmdbuffer == "-o") THEN ! output file
            ccmd = 4
         ELSEIF (cmdbuffer == "-l") THEN ! box lenght
            ccmd = 5
         ELSEIF (cmdbuffer == "-na") THEN ! number of atoms
            ccmd = 6
         ELSEIF (cmdbuffer == "-ns") THEN ! number of stpes
            ccmd = 7
         ELSEIF (cmdbuffer == "-ss") THEN ! starting step
            ccmd = 8
         ELSEIF (cmdbuffer == "-gf") THEN ! delta
            ccmd = 9
         ELSEIF (cmdbuffer == "-a") THEN ! smoothing factor, alpha
            ccmd = 10
         ELSEIF (cmdbuffer == "-ct") THEN ! cutoff for w
            ccmd = 11
         ELSEIF (cmdbuffer == "-tacc") THEN ! acceptor types
            ccmd = 12
         ELSEIF (cmdbuffer == "-tdon") THEN ! donor types
            ccmd = 13
         ELSEIF (cmdbuffer == "-tH") THEN ! hydrogen types
            ccmd = 14
         ELSEIF (cmdbuffer == "-h") THEN ! help
            WRITE(*,*) ""
            WRITE(*,*) " HB-mixture test program."
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (cmdbuffer == "-c") THEN ! convert from bohrradius to angstrom
            convert = .true.
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = .true.
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " Wrong usage."
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 1) THEN ! input file
               filename=trim(cmdbuffer) 
            ELSEIF (ccmd == 2) THEN ! number of gaussian
               READ(cmdbuffer,*) Nk 
            ELSEIF (ccmd == 3) THEN ! gaussian file
               gaussianfile=trim(cmdbuffer) 
            ELSEIF (ccmd == 4) THEN ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 5) THEN ! box lenght
               READ(cmdbuffer,*) lbox 
            ELSEIF (ccmd == 6) THEN ! numbers of atoms
               READ(cmdbuffer,*) natoms
            ELSEIF (ccmd == 7) THEN ! numbers of steps
               READ(cmdbuffer,*) nsteps
            ELSEIF (ccmd == 8) THEN ! starting step
               READ(cmdbuffer,*) startstep
            ELSEIF (ccmd == 9) THEN ! delta
               READ(cmdbuffer,*) delta
            ELSEIF (ccmd == 10) THEN ! smoothing factor, alpha
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 11) THEN ! cutoff for w
               READ(cmdbuffer,*) wcutoff
            ELSEIF (ccmd == 12) THEN ! accettor types
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtacc(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtacc(par_count)
            ELSEIF (ccmd == 13) THEN ! donor types
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtdon(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtdon(par_count)
            ELSEIF (ccmd == 14) THEN ! hydrogen types
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtH(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtH(par_count)
            ENDIF           
         ENDIF
      ENDDO
      !!!!! END PARSER
      
      ! Mandatory parameters
      IF ((ccmd == 0).OR.(filename.EQ."NULL").OR.(gaussianfile.EQ."NULL").OR.(Nk.EQ.-1)) THEN
         WRITE(*,*) ""
         WRITE(*,*) " Wrong usage. Insert the right parameters! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF
      
      ! Check the steps and the box parameters 
      CALL xyz_GetInfo(filename,dummyi1,dummyd1,dummyi2)
      ! control nstpes
      IF (nsteps == -1) THEN
         nsteps=dummyi2
      ENDIF
      ! control natoms
      IF (natoms == -1) THEN
         natoms=dummyi1
      ENDIF
      ! control lbox
      IF (lbox.LT.(0.0d0)) THEN
         lbox=dummyd1
      ENDIF

      ! we can allocate the vectors now
      ALLOCATE(positions(natoms,3))
      ALLOCATE(labels(natoms))
      
      ! get the labels of the atoms
      CALL xyz_GetLabels(filename,labels)
      ! define what is acceptor,donor and hydrogen
      ALLOCATE(masktypes(natoms))
      ! set to TYPE_NONE
      masktypes=masktypes*TYPE_NONE
      DO i=1,natoms
         ! set the mask using BITWISE OR OPERATOR
         IF(testtype(labels(i),vtH)) masktypes(i)=IOR(masktypes(i),TYPE_H)
         IF(testtype(labels(i),vtdon)) masktypes(i)=IOR(masktypes(i),TYPE_DONOR)
         IF(testtype(labels(i),vtacc)) masktypes(i)=IOR(masktypes(i),TYPE_ACCEPTOR)
      ENDDO
      
      ! Read gaussian parameters from the gaussian file
      OPEN(UNIT=12,FILE=gaussianfile,STATUS='OLD',ACTION='READ')
      !decomment this if the first line is a comment string
      !READ(12,*) cmdbuffer
      ALLOCATE(clusters(Nk))
      DO i=1,Nk     
         READ(12,*) clusters(i)%mean(1), clusters(i)%mean(2), &
                    clusters(i)%cov(1,1), clusters(i)%cov(1,2), &
                    clusters(i)%cov(2,1), clusters(i)%cov(2,2), &
                    clusters(i)%pk
            ! calculate once the Icovs matrix and the norm_const
            CALL gauss_prepare(clusters(i))
         ENDDO
      CLOSE(UNIT=12)
      
      ! Loop over the trajectory
      pos=0
      DO ts=1,nsteps
         IF ((MODULO(ts,delta)==0) .AND. (ts>=startstep)) THEN
            ! read this snapshot
            IF(verbose)THEN
               WRITE(*,*) "Step: ",ts
            ENDIF
            CALL XYZ_GetSnap(1,filename,NAtoms,pos,newpos,positions)
            !!!!!!! HBMIXTURE HERE! !!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ELSE
            ! discard this snapshot
            CALL XYZ_GetSnap(0,filename,NAtoms,pos,newpos,positions)
         ENDIF
         ! pointer for the position in the coordinates file
         ! update for the next pass
         pos=newpos 
      ENDDO
      ! end the loop over the trajectory
      
      DEALLOCATE(positions)
      DEALLOCATE(labels)
      DEALLOCATE(masktypes)
      DEALLOCATE(clusters)
      
      CONTAINS
      
         SUBROUTINE helpmessage
            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: minitest [-h] -i filename -gn Ngaussians -gf gaussianfile "
            WRITE(*,*) "                  [-o outputfile] [-l box_lenght] [-na Natoms] "
            WRITE(*,*) "                  [-ns total_steps] [-ss starting_step] [-ev delta] "
            WRITE(*,*) "                  [-a smoothing_factor] [-ct cutoff] "
            WRITE(*,*) "                  [-tacc Accettor_type1,Accettor_type2,...] "
            WRITE(*,*) "                  [-tdon Donor_type1,Donor_type2,...] "
            WRITE(*,*) "                  [-tH   Hydrogen_type1,Hydrogen_type2,...] "
            WRITE(*,*) "                  [-c] [-v] "
            WRITE(*,*) ""
         END SUBROUTINE helpmessage
      
         LOGICAL FUNCTION testtype(id,vtype)
            CHARACTER*4, INTENT(IN) :: id
            CHARACTER*4, DIMENSION(4), INTENT(IN) :: vtype
            INTEGER i
            testtype=.false.
            DO i=1,4
               IF(trim(id).EQ.trim(vtype(i)))THEN
                  testtype=.true.
                  EXIT
               ENDIF
            ENDDO
         END FUNCTION
      
      END PROGRAM minitest
