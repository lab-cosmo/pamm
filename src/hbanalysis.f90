! The main program which use the mixture library to HB.
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

      PROGRAM hbanalysis
         USE xyz
         USE mixture
      IMPLICIT NONE

      ! Types used by bitwise operators to control the atom type
      ! they must be power of 2
      INTEGER, PARAMETER :: TYPE_NONE=0
      INTEGER, PARAMETER :: TYPE_H=1
      INTEGER, PARAMETER :: TYPE_DONOR=2
      INTEGER, PARAMETER :: TYPE_ACCEPTOR=4

      CHARACTER*1024 :: filename, gaussianfile !, outputfile
      CHARACTER*1024 :: cmdbuffer,tmp
      ! system parameters
      INTEGER natoms
      DOUBLE PRECISION cell(3,3), icell(3,3), dummycell(3,3)
      DOUBLE PRECISION alpha, wcutoff
      INTEGER nsteps, startstep, delta
      INTEGER D ! dimensionality of the dataset
      ! gaussians
      INTEGER Nk
      ! array of gaussians
      TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:) :: clusters
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pks
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)  :: pnks
      ! vector that will contain the probabilities calculated using hb-mixture library
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: spa, spd, sph
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sa, sd, sh
      DOUBLE PRECISION, DIMENSION(3) :: vwR
      INTEGER ia,id,ih
      DOUBLE PRECISION rah,rdh,weight
      ! for the parser
      INTEGER ccmd
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string
      ! for a faster reading
      ! counters
      INTEGER i,ts,k

      LOGICAL convert,ptcm,nptm,weighted!,verbose
      INTEGER errdef

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: positions
      ! mask to define what is what
      INTEGER, ALLOCATABLE, DIMENSION(:) :: masktypes
      INTEGER vtghb(12),nghb
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vghb
      CHARACTER*4, ALLOCATABLE, DIMENSION(:) :: labels
      CHARACTER*1024 :: header, dummyc


      CHARACTER*4, DIMENSION(4) :: vtacc,vtdon,vtH

      DOUBLE PRECISION dummyd1,dummyd2
      INTEGER dummyi1,dummyi2
      INTEGER endf



      !! default values
      DO i=1,4
         vtacc(i)="NULL"
         vtdon(i)="NULL"
         vtH(i)="NULL"
      ENDDO
      filename="NULL"
      gaussianfile="NULL"
      !outputfile="out-Pad.dat"
      nk=-1
      ccmd=0
      cell=0.0d0
      nsteps=-1
      delta=1
      startstep=1
      alpha=1.0d0
      wcutoff=5.0d0
      convert = .false.
      !verbose = .false.
      ptcm = .false.
      nptm = .false.
      weighted= .false.   ! don't use the weights
      endf = 0
      errdef=0
      vtghb=-1
      D=3
      !!

      !!!!! PARSER
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-i") THEN ! input xyz
            ccmd = 1
         ELSEIF (cmdbuffer == "-gf") THEN ! file containing gaussian parmeters
            ccmd = 2
         !ELSEIF (cmdbuffer == "-o") THEN ! output file
         !   ccmd = 4
         ELSEIF (cmdbuffer == "-l") THEN ! box lenght
            ccmd = 3
         ELSEIF (cmdbuffer == "-ns") THEN ! number of stpes
            ccmd = 4
         ELSEIF (cmdbuffer == "-ss") THEN ! starting step
            ccmd = 5
         ELSEIF (cmdbuffer == "-ev") THEN ! delta
            ccmd = 6
         ELSEIF (cmdbuffer == "-a") THEN ! smoothing factor, alpha
            ccmd = 7
         ELSEIF (cmdbuffer == "-ct") THEN ! cutoff for w
            ccmd = 8
         ELSEIF (cmdbuffer == "-ta") THEN ! acceptor types
            ccmd = 9
         ELSEIF (cmdbuffer == "-td") THEN ! donor types
            ccmd = 10
         ELSEIF (cmdbuffer == "-th") THEN ! hydrogen types
            ccmd = 11
         ELSEIF (cmdbuffer == "-ghb") THEN ! gaussians used to describe the HB
            ccmd = 12
         ELSEIF (cmdbuffer == "-npt") THEN ! npt mode
            nptm = .true.
         ELSEIF (cmdbuffer == "-w") THEN ! weighted  mode
            weighted = .true.
         ELSEIF (cmdbuffer == "-h") THEN ! help
            WRITE(*,*) ""
            WRITE(*,*) " HB-mixture test program."
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (cmdbuffer == "-c") THEN ! convert from bohrradius to angstrom
            convert = .true.
         !ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            !verbose = .true.
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
            ELSEIF (ccmd == 2) THEN ! gaussian file
               gaussianfile=trim(cmdbuffer)
            !ELSEIF (ccmd == 4) THEN ! output file
            !   outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN ! box dimensions
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) cell(par_count,par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) cell(par_count,par_count)
            ELSEIF (ccmd == 4) THEN ! numbers of steps
               READ(cmdbuffer,*) nsteps
            ELSEIF (ccmd == 5) THEN ! starting step
               READ(cmdbuffer,*) startstep
            ELSEIF (ccmd == 6) THEN ! delta
               READ(cmdbuffer,*) delta
            ELSEIF (ccmd == 7) THEN ! smoothing factor, alpha
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 8) THEN ! cutoff for w
               READ(cmdbuffer,*) wcutoff
            ELSEIF (ccmd == 9) THEN ! accettor types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtacc(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtacc(par_count)
            ELSEIF (ccmd == 10) THEN ! donor types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtdon(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtdon(par_count)
            ELSEIF (ccmd == 11) THEN ! hydrogen types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtH(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtH(par_count)
            ELSEIF (ccmd == 12) THEN ! gaussians describing the HB
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtghb(par_count)
                  par_count = par_count + 1
                  IF(par_count>15)THEN
                     WRITE(*,*) "To much HBgaussians for me... Change the source code!"
                     CALL EXIT(-1)
                  ENDIF
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtghb(par_count)
               ! parameters readed. Now I know the number of the parameters passed
               ! and so I can allocate the vector I will use later
               nghb=par_count
               ALLOCATE(vghb(nghb))
               vghb(1:nghb)=vtghb(1:nghb)
            ENDIF
         ENDIF
      ENDDO
      !!!!! END PARSER

      ! Mandatory parameters
!      IF (filename.EQ."NULL") THEN
!         WRITE(*,*) ""
!         WRITE(*,*) " Error: insert an input file! "
!         CALL helpmessage
!         CALL EXIT(-1)
!      ENDIF

      IF (.NOT.(errdef.EQ.3)) THEN
         WRITE(*,*) ""
         WRITE(*,*) " Error: insert hydrongen, donor and acceptor species! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      IF (.NOT.ptcm) THEN
         ! Mandatory parameters
         IF (gaussianfile.EQ."NULL") THEN
            WRITE(*,*) ""
            WRITE(*,*) " Error: insert the file containing the gaussians parameters! "
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (vtghb(1).EQ.-1) THEN
            ! the user didn't passed the numbers of the gaussians to use
            ! we will use the first
            !IF(verbose) WRITE(*,*) " Will use the first gaussian to describe the HB. "
            ALLOCATE(vghb(1))
            vghb(1)=1
            nghb=1
         ENDIF
      ENDIF

      CALL invmatrix(3,cell,icell)

      ! Loop over the trajectory
      ts=0
      DO
         ts=ts+1
         IF (nsteps.NE.-1 .AND. nsteps.LT.ts) EXIT  ! time to go!

         IF (ts-startstep .lt. 0  .or. MODULO(ts-startstep,delta)/=0) THEN
            ! skip the frame
            READ(5,*,IOSTAT=i) natoms
            READ(5,'(A)') header
            DO i=1,natoms
               READ(5,'(A)',IOSTAT=endf) header
            ENDDO
            IF(endf>0) error STOP "*** Error occurred while reading file. ***"
            IF(endf<0) EXIT
            CYCLE
         ELSE
            CALL xyz_read(5,natoms,header,labels,positions,endf)
            IF(endf<0) EXIT
            IF (ts-startstep == 0) THEN ! first step, must allocate stuff
                ! we can now define and inizialize the masks and Gaussians
                IF (.NOT.ptcm) THEN
                    ! Read gaussian parameters from the gaussian file
                    OPEN(UNIT=12,FILE=gaussianfile,STATUS='OLD',ACTION='READ')
                    ! this allocation is just to avoid the rising of errors
                    ALLOCATE(clusters(1),pks(1))
                    CALL readgaussians(12,D,nk,clusters,pks)
                    CLOSE(UNIT=12)

                    ALLOCATE(spa(nk,natoms), spd(nk,natoms), sph(nk,natoms))
                    ALLOCATE(sa(natoms), sd(natoms), sh(natoms), pnks(nk))
                    spa=0.0d0
                    spd=0.0d0
                    sph=0.0d0
                ENDIF
                ALLOCATE(masktypes(natoms))
                ! set to TYPE_NONE
                masktypes=TYPE_NONE
                DO i=1,natoms
                   ! set the mask using BITWISE OR OPERATOR
                   IF(testtype(labels(i),vtH)) masktypes(i)=IOR(masktypes(i),TYPE_H)
                   IF(testtype(labels(i),vtdon)) masktypes(i)=IOR(masktypes(i),TYPE_DONOR)
                   IF(testtype(labels(i),vtacc)) masktypes(i)=IOR(masktypes(i),TYPE_ACCEPTOR)
                ENDDO
            ENDIF

            IF (nptm .or. cell(1,1) == 0.0d0) THEN ! Variable cell! Try to read the cell parameters from the header
               READ(header, *) dummyc,dummyc,cell(1,1),cell(2,2),cell(3,3)
               CALL invmatrix(3,cell,icell)
            END IF

            ! here the core of the program
            DO ih=1,natoms ! loop over H
               IF (IAND(masktypes(ih),TYPE_H).EQ.0) CYCLE
               DO id=1,natoms
                  IF (IAND(masktypes(id),TYPE_DONOR).EQ.0 .OR. ih.EQ.id) CYCLE
                  ! Distance donor-hydrogen
                  CALL separation(cell,icell,positions(:,ih),positions(:,id),rdh)
                  IF(rdh.GT.wcutoff) CYCLE  ! if one of the distances is greater
                                            ! than the cutoff, we can already discard the D-H pair
                  DO ia=1,natoms

                     IF (IAND(masktypes(ia),TYPE_ACCEPTOR).EQ.0 &
                        .OR. (ia.EQ.id).OR.(ia.EQ.ih)) CYCLE
                     ! Distance acceptor-hydrogen
                     CALL separation(cell,icell,positions(:,ih),positions(:,ia),rah)
                     vwR(2)=rdh+rah
                     IF(vwR(2).GT.wcutoff) CYCLE
                     vwR(1)=rdh-rah
                     ! Calculate the distance donor-acceptor (R)
                     CALL separation(cell,icell,positions(:,id),positions(:,ia),vwR(3))

                     weight=1.0d0
                     ! weight = 1/J = 1/((w*w-v*v)*R)
                     IF(weighted) weight=1.0d0/((vwR(2)*vwR(2)-vwR(1)*vwR(1))*vwR(3))

                     IF(ptcm)THEN ! -P mode
                        ! write out : v,w,R  and weight
                        WRITE(*,"(3(A1,ES21.8E4))",ADVANCE = "NO")  " ",vwR(1)," ",vwR(2)," ",vwR(3)
                        WRITE(*,"(A1,ES21.8E4)") " ", weight
                     ELSE ! hbmix
                        ! call the library
                        CALL GetP(3,vwR,weight,alpha,Nk,clusters,pks,pnks)
                        ! get the pnks
                        ! cumulate the probabilities of all the triplets in wich an the atoms ar involved
                        sph(:,ih) = sph(:,ih) + pnks(:)
                        spa(:,ia) = spa(:,ia) + pnks(:)
                        spd(:,id) = spd(:,id) + pnks(:)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

            ! write the results

            IF(.not.ptcm)THEN
               ! sum over the gaussians of interest
               sa=0.0d0
               sd=0.0d0
               sh=0.0d0
               DO i=1,nghb
                  sh(:)=sh(:)+sph(vghb(i),:)
                  sd(:)=sd(:)+spd(vghb(i),:)
                  sa(:)=sa(:)+spa(vghb(i),:)
               ENDDO
               ! write results to a formatted output
               ! use positions as a dummy vector to store sh, sd, sa
               positions(1,:) = sh
               positions(2,:) = sd
               positions(3,:) = sa
               CALL xyz_write(6,natoms,header,labels,positions)
            ENDIF
         ENDIF
      ENDDO
      ! end the loop over the trajectory

      DEALLOCATE(positions)
      DEALLOCATE(labels)
      DEALLOCATE(masktypes)
      IF (.NOT.ptcm) THEN
         DEALLOCATE(clusters,pks,pnks)
         DEALLOCATE(spa, sph, spd, sa, sh, sd, vghb)
         !CLOSE(UNIT=7)
      ENDIF

      CONTAINS

         SUBROUTINE helpmessage
            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: hbanalysis [-h] [-P] -i filename [-l lx,ly,lz] " ![-o outputfile]
            WRITE(*,*) "                     -ta A1,A2,... -td D1,D2,... -th H1,H2,... "
            WRITE(*,*) "                    [-gf gaussianfile] [-ghb 1,2,..] [-a smoothing_factor] "
            WRITE(*,*) "                    [-ct cutoff] [-ev delta] "
            WRITE(*,*) "                    [-ns total_steps] [-ss starting_step] [-npt] [-c] " ![-v] "
            WRITE(*,*) ""
            WRITE(*,*) " Description ... Two modalities :  "
            WRITE(*,*) ""
            WRITE(*,*) " 1. Extract v,w,R and eventually the weight from the data."
            WRITE(*,*) "    This is the pre-processing mode: call it with the -P flag."
            WRITE(*,*) " 2. Analyze the hydrogen bond patterns using a probabilistic model "
            WRITE(*,*) "    based on a mixture of gaussians "
            WRITE(*,*) ""
            WRITE(*,*) "   -h                   : Print this message "
            WRITE(*,*) "   -P                   : Pre-processing mode "
            WRITE(*,*) "   -i  Input_file       : File containin the input data (XYZ format)"
            WRITE(*,*) "   -o  Output_file      : Output file. "
            WRITE(*,*) "   -l  lx,ly,lz         : Box lenghts (orthorombic box) "
            WRITE(*,*) "   -ta A1,A2,...        : Namelist of acceptors "
            WRITE(*,*) "   -ta D1,D2,...        : Namelist of donors "
            WRITE(*,*) "   -ta H1,H2,...        : Namelist of hydrogens "
            WRITE(*,*) "   -gf Gaussians_file   : File containing gaussians data "
            WRITE(*,*) "   -gh 1,2,...          : Index of the gaussians that describe the HB "
            WRITE(*,*) "   -a  smoothing_factor : Apply a smoothing factor to the gaussians "
            WRITE(*,*) "   -ct cutoff           : Apply a distance cut-off during the calculation "
            WRITE(*,*) "   -ev delta            : Stride while reading data from the file "
            WRITE(*,*) "   -ns total_steps      : Apply a distance cut-off during the calculation "
            WRITE(*,*) "   -ss starting_step    : Step from wich to start reading the data "
            WRITE(*,*) "   -npt                 : NPT mode "
            WRITE(*,*) "   -c                   : Convert from atomic units (Bohr) to Angstrom "
            WRITE(*,*) "   -w                   : Data points must be weighted "
            ! WRITE(*,*) "   -v                   : Verobose mode "
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

      END PROGRAM hbanalysis
