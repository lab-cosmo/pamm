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
         USE matrixinverse
         USE gaussian
         USE xyz
         USE hbmixture
      IMPLICIT NONE

      ! Conversion constant from bohrradius to angstrom
      DOUBLE PRECISION, PARAMETER :: bohr=0.5291772192171717171717171717171717

      CHARACTER*70 :: filename, gaussianfile, outputfile
      CHARACTER*1024 :: cmdbuffer,tmp
      ! system parameters
      INTEGER natoms
      !DOUBLE PRECISION, DIMENSION(3) :: box
      DOUBLE PRECISION cell(3,3), icell(3,3), dummycell(3,3)
      DOUBLE PRECISION box(3)
      DOUBLE PRECISION alpha, wcutoff
      INTEGER nsteps, startstep, delta
      ! gaussians
      INTEGER Nk
      ! array of gaussians
      TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:) :: clusters
      ! vector that will contain the probabilities calculated using hb-mixture library
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: spa, spd, sph
      ! for the parser
      INTEGER ccmd
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string
      ! for a faster reading
      ! counters
      INTEGER i,ts,k
      ! parmater for seeking in the input file
      INTEGER pos,newpos


      LOGICAL verbose,convert,ptcm
      INTEGER errdef

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
      cell=0.0d0
      nsteps=-1
      natoms=-1
      delta=1
      startstep=1
      alpha=1.0d0
      wcutoff=5.0d0
      convert = .false.
      verbose = .false.
      ptcm = .false.
      errdef=0 
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
         ELSEIF (cmdbuffer == "-L") THEN ! box lenght
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
         ELSEIF (cmdbuffer == "-ev") THEN ! starting step
            ccmd = 15
         ELSEIF (cmdbuffer == "-h") THEN ! help
            WRITE(*,*) ""
            WRITE(*,*) " HB-mixture test program."
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (cmdbuffer == "-c") THEN ! convert from bohrradius to angstrom
            convert = .true.
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = .true.
         ELSEIF (cmdbuffer == "-P") THEN ! PTC mode
            ptcm = .true.
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " Wrong usage. Insert the right parameters!"
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
            ELSEIF (ccmd == 5) THEN ! box dimensions
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) box(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) box(par_count)
               cell(1,1)=box(1)
               cell(2,2)=box(2)
               cell(3,3)=box(3)
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
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtacc(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtacc(par_count)
            ELSEIF (ccmd == 13) THEN ! donor types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtdon(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtdon(par_count)
            ELSEIF (ccmd == 14) THEN ! hydrogen types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtH(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtH(par_count)
            ELSEIF (ccmd == 15) THEN ! cutoff for w
               READ(cmdbuffer,*) delta
            ENDIF
         ENDIF
      ENDDO
      !!!!! END PARSER

      ! Mandatory parameters
      IF (filename.EQ."NULL") THEN
         WRITE(*,*) ""
         WRITE(*,*) " Error: insert an input file! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF
      
      IF (.NOT.(errdef.EQ.3)) THEN
         WRITE(*,*) ""
         WRITE(*,*) " Error: insert hydrongen, donor and acceptor species! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF
      
      IF (.NOT.ptcm) THEN 
         ! Mandatory parameters
         IF ((gaussianfile.EQ."NULL").OR.(Nk.EQ.-1)) THEN
            WRITE(*,*) ""
            WRITE(*,*) " Error: insert the gaussians parameters! "
            CALL helpmessage
            CALL EXIT(-1)
         ENDIF
      ENDIF
      
      ! Check the steps and the box parameters
      
      CALL xyz_GetInfo(filename,dummyi1,dummycell,dummyi2) 
      ! control nstpes
      IF (nsteps == -1) THEN
         nsteps=dummyi2
      ENDIF
      ! control natoms
      IF (natoms == -1) THEN
         natoms=dummyi1
      ENDIF
      ! control the cell
      IF (cell(1,1)==(0.0d0) .OR. cell(2,2)==(0.0d0) .OR. cell(3,3)==(0.0d0)) THEN
         cell=dummycell
         IF(convert) cell=cell*bohr
      ENDIF
      
      ! get the inverse of the cell
      CALL inv(cell,icell)
      
      ! we can allocate the vectors now
      ALLOCATE(positions(3,natoms))
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
      
      IF (.NOT.ptcm) THEN 
         ! Read gaussian parameters from the gaussian file
         OPEN(UNIT=12,FILE=gaussianfile,STATUS='OLD',ACTION='READ')
         !decomment this if the first line is a comment string
         !READ(12,*) cmdbuffer
         ALLOCATE(clusters(Nk))
         DO i=1,Nk
            READ(12,*) clusters(i)%mean(1), clusters(i)%mean(2), clusters(i)%mean(3), &
                       clusters(i)%cov(1,1), clusters(i)%cov(2,1), clusters(i)%cov(3,1), &
                       clusters(i)%cov(1,2), clusters(i)%cov(2,2), clusters(i)%cov(3,2), &
                       clusters(i)%cov(1,3), clusters(i)%cov(2,3), clusters(i)%cov(3,3), &
                       dummyd2
            !! the follow because I'm using the gaussianmix of Michele (that save log(pk))
            clusters(i)%pk=exp(dummyd2) 
            ! calculate once the Icovs matrix and the norm_const
            CALL gauss_prepare(clusters(i))     
            
         ENDDO
         CLOSE(UNIT=12)

         ! we can now define and inizialize the probabilities vector
         ALLOCATE(spa(nk,natoms), spd(nk,natoms), sph(nk,natoms))
         spa=0.0d0
         spd=0.0d0
         sph=0.0d0
         ! outputfile
         OPEN(UNIT=7,FILE=outputfile,STATUS='REPLACE',ACTION='WRITE')
      ENDIF
      
      OPEN(UNIT=11,FILE=filename)
      ! Loop over the trajectory
      pos=0
      DO ts=1,nsteps
         IF ((MODULO(ts,delta)==0) .AND. (ts>=startstep)) THEN
            ! read this snapshot
            IF(verbose) WRITE(*,*) "Step: ",ts
            CALL xyz_GetSnap(1,11,natoms,positions)
            IF(convert) positions=positions*bohr
            
            IF(ptcm)THEN
               CALL write_vwd(natoms,cell,icell,wcutoff,masktypes,positions)
            ELSE
               !!!!!!! HBMIXTURE HERE! !!!!!!!
               CALL hbmixture_GetGMMP(natoms,cell,icell,alpha,wcutoff,positions,masktypes, &
                                      nk,clusters,sph,spd,spa)
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               ! write results to a formatted output
               CALL write_output(7,natoms,nk,cell,ts,positions,sph(1,:),spd(1,:),spa(1,:))
            ENDIF

         ELSE
            ! discard this snapshot
            CALL XYZ_GetSnap(0,11,natoms,positions)
         ENDIF
         ! pointer for the position in the coordinates file
         ! update for the next pass
         pos=newpos
      ENDDO
      ! end the loop over the trajectory

      DEALLOCATE(positions)
      DEALLOCATE(labels)
      DEALLOCATE(masktypes)
      CLOSE(UNIT=11)
      IF (.NOT.ptcm) THEN 
         DEALLOCATE(clusters)
         DEALLOCATE(spa, sph, spd)
         CLOSE(UNIT=7)
      ENDIF

      CONTAINS

         SUBROUTINE helpmessage
            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: minitest [-h] [-P] -i filename  "
            WRITE(*,*) "                   -tacc Accettor_type1,Accettor_type2,... "
            WRITE(*,*) "                   -tdon Donor_type1,Donor_type2,... "
            WRITE(*,*) "                   -tH   Hydrogen_type1,Hydrogen_type2,... "
            WRITE(*,*) "                  [-gn Ngaussians] [-gf gaussianfile] "
            WRITE(*,*) "                  [-o outputfile] [-L lx,ly,lz] [-na Natoms] "
            WRITE(*,*) "                  [-ns total_steps] [-ss starting_step] [-ev delta] "
            WRITE(*,*) "                  [-a smoothing_factor] [-ct cutoff] "
            WRITE(*,*) "                  [-c] [-v] "
            WRITE(*,*) ""
         END SUBROUTINE helpmessage

         SUBROUTINE write_output(filen,natoms,nk,cell,ts,positions,sh,sd,sa)
            INTEGER, INTENT(IN) :: filen
            INTEGER, INTENT(IN) :: natoms
            INTEGER, INTENT(IN) :: nk
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell
            INTEGER, INTENT(IN) :: ts
            DOUBLE PRECISION, DIMENSION(3,natoms), INTENT(IN) :: positions
            DOUBLE PRECISION, DIMENSION(natoms), INTENT(IN) :: sh, sd, sa

            ! header
            WRITE(filen,"(I4)") natoms
            WRITE(filen,"(a,F11.7,a,F11.7,a,F11.7,a,I10)") &
                 "# CELL(abc): ",cell(1,1)," ",cell(2,2)," ",cell(3,3)," Step: ",ts
            ! body
            DO i=1,natoms
               WRITE(filen,"(A2,A1,F11.7,A1,F11.7,A1,F11.7)", advance='no') &
                    trim(labels(i)), " ", positions(1,i), " ", &
                    positions(2,i), " ", positions(3,i)
               WRITE(filen,"(3(A1,ES18.9))", advance='no') " ", sh(i), &
                       " ", sd(i), " ", sa(i)
               WRITE(filen,*)
            ENDDO
         END SUBROUTINE write_output
         
         SUBROUTINE write_vwd(natoms,cell,icell,wcutoff,masktypes,positions)
            ! Calculate the probabilities
            ! ...
            ! Args:
            !    param: descript
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: icell
            DOUBLE PRECISION, INTENT(IN) :: wcutoff
            INTEGER, DIMENSION(natoms), INTENT(IN) :: masktypes
            DOUBLE PRECISION, DIMENSION(3,natoms), INTENT(IN) :: positions
            
            INTEGER ih,id,ia
            DOUBLE PRECISION,DIMENSION(3) :: vwd
            DOUBLE PRECISION rah,rdh
            
            DO ih=1,natoms ! loop over H
               IF (IAND(masktypes(ih),TYPE_H).EQ.0) CYCLE
               DO id=1,natoms
                  IF (IAND(masktypes(id),TYPE_DONOR).EQ.0 .OR. ih.EQ.id) CYCLE
                  CALL separation(cell,icell,positions(:,ih),positions(:,id),rdh)
                  IF(rdh.GT.wcutoff) CYCLE  ! if one of the distances is greater 
                                   !than the cutoff, we can already discard the D-H pair
                  DO ia=1,natoms  
                     IF (IAND(masktypes(ia),TYPE_ACCEPTOR).EQ.0 &
                        .OR. (ia.EQ.id).OR.(ia.EQ.ih)) CYCLE

                     CALL separation(cell,icell,positions(:,ih),positions(:,ia),rah)
                     vwd(2)=rdh+rah
                     IF(vwd(2).GT.wcutoff) CYCLE
                     vwd(1)=rdh-rah
                     ! Calculate the distance donor-acceptor
                     CALL separation(cell,icell,positions(:,id),positions(:,ia),vwd(3))
                     WRITE(*,*) " ",vwd(1)," ",vwd(2)," ",vwd(3)
                     !write(*,*) ia,id,ih," ",rah,rdh,rad
                  ENDDO
               ENDDO
            ENDDO
         END SUBROUTINE write_vwd

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
