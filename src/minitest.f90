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
         USE xyz
         USE hbmixture
      IMPLICIT NONE

      CHARACTER*70 :: filename, gaussianfile, outputfile
      CHARACTER*1024 :: cmdbuffer,tmp
      ! system parameters
      INTEGER natoms
      DOUBLE PRECISION cell(3,3), icell(3,3), dummycell(3,3)
      DOUBLE PRECISION box(3)
      DOUBLE PRECISION alpha, wcutoff
      INTEGER nsteps, startstep, delta
      ! gaussians
      INTEGER Nk
      ! array of gaussians
      TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:) :: clusters
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pks
      ! vector that will contain the probabilities calculated using hb-mixture library
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: spa, spd, sph
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sa, sd, sh
      ! for the parser
      INTEGER ccmd
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string
      ! for a faster reading
      ! counters
      INTEGER i,ts,k

      LOGICAL verbose,convert,ptcm1,ptcm2,nptm
      INTEGER errdef

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: positions
      ! mask to define what is what
      INTEGER, ALLOCATABLE, DIMENSION(:) :: masktypes
      INTEGER vtghb(12),nghb
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vghb
      CHARACTER*4, ALLOCATABLE, DIMENSION(:) :: labels
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
      ptcm1 = .false.
      ptcm2 = .false.
      nptm = .false.
      endf = 0
      errdef=0
      vtghb=-1
      !!

      !!!!! PARSER
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-i") THEN ! input xyz
            ccmd = 1
         ELSEIF (cmdbuffer == "-gn") THEN ! total number of gaussians
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
         ELSEIF (cmdbuffer == "-ev") THEN ! delta
            ccmd = 9
         ELSEIF (cmdbuffer == "-a") THEN ! smoothing factor, alpha
            ccmd = 10
         ELSEIF (cmdbuffer == "-ct") THEN ! cutoff for w
            ccmd = 11
         ELSEIF (cmdbuffer == "-ta") THEN ! acceptor types
            ccmd = 12
         ELSEIF (cmdbuffer == "-td") THEN ! donor types
            ccmd = 13
         ELSEIF (cmdbuffer == "-th") THEN ! hydrogen types
            ccmd = 14
         ELSEIF (cmdbuffer == "-ghb") THEN ! gaussians used to describe the HB
            ccmd = 15
         ELSEIF (cmdbuffer == "-npt") THEN ! npt mode
            nptm = .true.
         ELSEIF (cmdbuffer == "-h") THEN ! help
            WRITE(*,*) ""
            WRITE(*,*) " HB-mixture test program."
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (cmdbuffer == "-c") THEN ! convert from bohrradius to angstrom
            convert = .true.
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = .true.
         ELSEIF (cmdbuffer == "-P1") THEN ! PTC mode
            ptcm1 = .true.
         ELSEIF (cmdbuffer == "-P2") THEN ! Perimetric coordinates mode
            ptcm2 = .true.
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
            ELSEIF (ccmd == 15) THEN ! gaussians describing the HB
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

      IF (ptcm1 .and. ptcm2) THEN
         WRITE(*,*) ""
         WRITE(*,*) " Error: you have to chose! -P1 or -P2, not both!! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

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

      IF (.NOT.ptcm1 .and. .NOT.ptcm2) THEN
         ! Mandatory parameters
         IF ((gaussianfile.EQ."NULL").OR.(Nk.EQ.-1)) THEN
            WRITE(*,*) ""
            WRITE(*,*) " Error: insert the gaussians parameters! "
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (vtghb(1).EQ.-1) THEN
            ! the user didn't passed the numbers of the gaussians to use
            ! we will use the first
            IF(verbose) WRITE(*,*) " Will use the first gaussian to describe the HB. "
            ALLOCATE(vghb(1))
            vghb(1)=1
            nghb=1
         ENDIF
      ENDIF
      
      OPEN(UNIT=11,FILE=filename)
      CALL xyz_read(1,nptm,convert,11,natoms,positions,labels,cell,icell,endf)
      CLOSE(UNIT=11)

      ! define what is acceptor,donor and hydrogen
      ALLOCATE(masktypes(natoms))
      ! set to TYPE_NONE
      masktypes=TYPE_NONE
      DO i=1,natoms
         ! set the mask using BITWISE OR OPERATOR
         IF(testtype(labels(i),vtH)) masktypes(i)=IOR(masktypes(i),TYPE_H)
         IF(testtype(labels(i),vtdon)) masktypes(i)=IOR(masktypes(i),TYPE_DONOR)
         IF(testtype(labels(i),vtacc)) masktypes(i)=IOR(masktypes(i),TYPE_ACCEPTOR)
      ENDDO

      IF (.NOT.ptcm1 .and. .NOT.ptcm2) THEN
         ! HB-mixture mode
         ! Read gaussian parameters from the gaussian file
         OPEN(UNIT=12,FILE=gaussianfile,STATUS='OLD',ACTION='READ')
         ! skip the first two comment lines
         READ(12,*) cmdbuffer
         READ(12,*) cmdbuffer
         ! read the number of gaussians form the file
         READ(12,*) dummyi1
         ! test the inserted value with the one from the file
         IF(Nk.ne.dummyi1)THEN
            WRITE(0,*) "Number of Gaussians on command line does not match init.", &
                       "Will read what I can."
         ENDIF
         ALLOCATE(clusters(Nk))
         ALLOCATE(pks(Nk))
         DO i=1,Nk
            CALL readgaussfromfile(12,clusters(i),pks(i))
         ENDDO
         CLOSE(UNIT=12)

         ! we can now define and inizialize the probabilities vector
         ALLOCATE(spa(nk,natoms), spd(nk,natoms), sph(nk,natoms))
         ALLOCATE(sa(natoms), sd(natoms), sh(natoms))
         spa=0.0d0
         spd=0.0d0
         sph=0.0d0
         ! outputfile
         OPEN(UNIT=7,FILE=outputfile,STATUS='REPLACE',ACTION='WRITE')
      ENDIF

      OPEN(UNIT=11,FILE=filename)
      ! Loop over the trajectory
      ts=0
      DO
         
         ts=ts+1
         IF(nsteps.NE.-1)THEN
            IF(nsteps.LT.ts) EXIT
         ENDIF
         IF ((MODULO(ts,delta)==0) .AND. (ts>=startstep)) THEN
            ! read this snapshot
            IF(verbose) WRITE(*,*) "Step: ",ts
            CALL xyz_read(1,nptm,convert,11,natoms,positions,labels,cell,icell,endf)
            IF(endf<0)EXIT
            IF(ptcm1)THEN
               
               CALL write_vwd(natoms,cell,icell,wcutoff,masktypes,positions)
            ELSEIF(ptcm2)THEN
               CALL write_xyz(natoms,cell,icell,wcutoff,masktypes,positions)
            ELSE
               !!!!!!! HBMIXTURE HERE! !!!!!!!
               CALL hbmixture_GetGMMP(natoms,cell,icell,alpha,wcutoff,positions,masktypes, &
                                      nk,clusters,pks,sph,spd,spa)
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! sum all the gaussians describing the HB
               sa=0.0d0
               sd=0.0d0
               sh=0.0d0
               DO i=1,nghb
                  sh(:)=sh(:)+sph(vghb(i),:)
                  sd(:)=sd(:)+spd(vghb(i),:)
                  sa(:)=sa(:)+spa(vghb(i),:)
               ENDDO
               ! write results to a formatted output
               CALL xyz_write(7,natoms,cell,ts,labels,positions,sh,sd,sa)
            ENDIF

         ELSE
            ! discard this snapshot
            CALL xyz_read(0,nptm,convert,11,natoms,positions,labels,cell,icell,endf)
            IF(endf<0)EXIT
         ENDIF
      ENDDO
      ! end the loop over the trajectory

      DEALLOCATE(positions)
      DEALLOCATE(labels)
      DEALLOCATE(masktypes)
      CLOSE(UNIT=11)
      IF (.NOT.ptcm1) THEN
         DEALLOCATE(clusters,pks)
         DEALLOCATE(spa, sph, spd, sa, sh, sd, vghb)
         CLOSE(UNIT=7)
      ENDIF

      CONTAINS

         SUBROUTINE helpmessage
            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: minitest [-P1/-P2] -i filename [-o outputfile] [-l lx,ly,lz] "
            WRITE(*,*) "                   -ta acc1,acc2,... -td don1,don2,... -th hyd1,hyd2,... "
            WRITE(*,*) "                  [-gn Ngaussians] [-gf gaussianfile] [-ghb ghb1,ghb2,..]"
            WRITE(*,*) "                  [-a smoothing_factor] [-ct cutoff] [-ev delta] "
            WRITE(*,*) "                  [-na Natoms] [-ns total_steps] [-ss starting_step]"
            WRITE(*,*) "                  [-npt] [-c] [-v] [-h] "
            WRITE(*,*) ""
         END SUBROUTINE helpmessage

         SUBROUTINE readgaussfromfile(fileid,gaussp,pk)
            ! Read a line from the file and get the paramters for the related gaussian
            !
            ! Args:
            !    fileid: the file containing the gaussians parameters
            !    gaussp: type_gaussian container in wich we store the gaussian parameters
            !    lpk: logarithm of the Pk associated to the gaussian

            INTEGER, INTENT(IN) :: fileid
            TYPE(gauss_type) , INTENT(INOUT) :: gaussp
            DOUBLE PRECISION, INTENT(INOUT) :: pk

            READ(fileid,*) gaussp%mean(1), gaussp%mean(2), gaussp%mean(3), &
                           gaussp%cov(1,1), gaussp%cov(2,1), gaussp%cov(3,1), &
                           gaussp%cov(1,2), gaussp%cov(2,2), gaussp%cov(3,2), &
                           gaussp%cov(1,3), gaussp%cov(2,3), gaussp%cov(3,3), &
                           pk

            CALL gauss_prepare(gaussp)

         END SUBROUTINE readgaussfromfile

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

         SUBROUTINE write_xyz(natoms,cell,icell,cutoff,masktypes,positions)
            ! Calculate the probabilities
            ! ...
            ! Args:
            !    param: descript
            INTEGER, INTENT(IN) :: natoms
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: cell
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: icell
            DOUBLE PRECISION, INTENT(IN) :: cutoff
            INTEGER, DIMENSION(natoms), INTENT(IN) :: masktypes
            DOUBLE PRECISION, DIMENSION(3,natoms), INTENT(IN) :: positions

            INTEGER ih,id,ia
            DOUBLE PRECISION,DIMENSION(3) :: vwd
            DOUBLE PRECISION rah,rdh,rad

            DO ih=1,natoms ! loop over H
               IF (IAND(masktypes(ih),TYPE_H).EQ.0) CYCLE
               DO id=1,natoms
                  IF (IAND(masktypes(id),TYPE_DONOR).EQ.0 .OR. ih.EQ.id) CYCLE
                  CALL separation(cell,icell,positions(:,ih),positions(:,id),rdh)
                  IF(rdh.GT.(cutoff)) CYCLE
                  DO ia=1,natoms
                     IF (IAND(masktypes(ia),TYPE_ACCEPTOR).EQ.0 &
                        .OR. (ia.EQ.id).OR.(ia.EQ.ih)) CYCLE

                     CALL separation(cell,icell,positions(:,ih),positions(:,ia),rah)
                     IF(rah.GT.(cutoff)) CYCLE
                     ! Calculate the distance donor-acceptor
                     CALL separation(cell,icell,positions(:,id),positions(:,ia),rad)
                     vwd(1)=rdh+rah-rad ! x
                     vwd(2)=rdh-rah+rad ! y
                     vwd(3)=-rdh+rah+rad ! z
                     !IF(rdh.GT.cutoff .or. rah.GT.cutoff .or. rad.GT.cutoff) CYCLE
                     WRITE(*,*) " ",vwd(1)," ",vwd(2)," ",vwd(3)
                     !write(*,*) ia,id,ih," ",rah,rdh,rad
                  ENDDO
               ENDDO
            ENDDO
         END SUBROUTINE write_xyz

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
