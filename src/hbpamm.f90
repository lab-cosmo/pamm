! An application of PAMM to the hydrogen bond. 
! Uses libpamm for the evaluation of sDHA
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

      PROGRAM hbpamm
         USE xyz
         USE libpamm
      IMPLICIT NONE
   
      ! Types used by bitwise operators to control the atom type
      ! they must be power of 2
      INTEGER, PARAMETER :: TYPE_NONE=0
      INTEGER, PARAMETER :: TYPE_H=1
      INTEGER, PARAMETER :: TYPE_DONOR=2
      INTEGER, PARAMETER :: TYPE_ACCEPTOR=4
      INTEGER, PARAMETER :: MAXPARS = 4
      CHARACTER(LEN=1024) :: filename, clusterfile 
      CHARACTER(LEN=1024) :: cmdbuffer

      !system parameters
      INTEGER natoms
      DOUBLE PRECISION cell(3,3), icell(3,3)
      DOUBLE PRECISION alpha, mucutoff

      !to remove
      ! gaussians
      INTEGER Nk
      ! array of gaussians
      TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:) :: clusters
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pks
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)  :: pnks
      ! vector that will contain the probabilities calculated using hb-mixture library
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: spa, spd, sph
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sa, sd, sh
      DOUBLE PRECISION, DIMENSION(3) :: x  ! [nu,mu,r]
      INTEGER ia,id,ih
      DOUBLE PRECISION rah,rdh,wfactor
      ! for the parser
      INTEGER ccmd
      INTEGER commas(MAXPARS), par_count  ! stores the index of commas in the parameter string
      ! for a faster reading
      ! counters
      INTEGER i,ts
      LOGICAL convert,dopamm,nptm,weighted!,verbose
      INTEGER errdef
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: positions
      ! mask to define what is what
      INTEGER, ALLOCATABLE, DIMENSION(:) :: masktypes
      INTEGER vtghb(12),nghb
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vghb
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: labels
      CHARACTER(LEN=1024) :: header, dummyc
      CHARACTER(LEN=4), DIMENSION(4) :: vtacc,vtdon,vtH
      INTEGER endf
      !default values
      DO i=1,4
         vtacc(i)="NULL"
         vtdon(i)="NULL"
         vtH(i)="NULL"
      ENDDO
      filename="NULL"
      clusterfile="NULL"
      nk=-1
      ccmd=0
      cell=0.0d0
      alpha=1.0d0
      mucutoff=5.0d0
      convert = .false.
      dopamm = .false.
      nptm = .false.
      weighted = .false.   ! don't use wfactor
      endf = 0
      errdef=0
      vtghb=-1
      !!!!

      !!!!! PARSER
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-h") THEN ! help
            ccmd = 1
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (cmdbuffer == "-gf") THEN ! file containing gaussian parmeters
            ccmd = 2
         ELSEIF (cmdbuffer == "-l") THEN ! box length
            ccmd = 3
         ELSEIF (cmdbuffer == "-a") THEN ! smoothing factor, alpha
            ccmd = 4
         ELSEIF (cmdbuffer == "-ct") THEN ! cutoff for mu
            ccmd = 5
         ELSEIF (cmdbuffer == "-ta") THEN ! acceptor types
            ccmd = 6
         ELSEIF (cmdbuffer == "-td") THEN ! donor types
            ccmd = 7
         ELSEIF (cmdbuffer == "-th") THEN ! hydrogen types
            ccmd = 8
         ELSEIF (cmdbuffer == "-ghb") THEN ! gaussians used to describe the HB
            ccmd = 9
         ELSEIF (cmdbuffer == "-npt") THEN ! npt mode
            nptm = .true.
         ELSEIF (cmdbuffer == "-w") THEN ! weighted  mode
            weighted = .true.
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " Wrong usage. Insert the right parameters!"
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 2) THEN ! gaussian file
               dopamm=.true.
               clusterfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN ! box dimensions
               ! read the box dimensions
               ! if not specified here, they will be read from the 
               ! the input stream, and in that case the format has to be:
               ! # CELL(abc): XX YY ZZ
               ! in here just an orthorombic cell can be used, so we read
               ! just the diagonal of the matrix.
               ! It is trivial to extend the code to the general case. 
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) cell(par_count,par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) cell(par_count,par_count)
            ELSEIF (ccmd == 4) THEN ! smoothing factor, alpha
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 5) THEN ! cutoff in the mu=d(AH)+d(dH) distance
               READ(cmdbuffer,*) mucutoff
            ELSEIF (ccmd == 6) THEN ! acceptor types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtacc(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtacc(par_count)
            ELSEIF (ccmd == 7) THEN ! donor types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtdon(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtdon(par_count)
            ELSEIF (ccmd == 8) THEN ! hydrogen types
               errdef=errdef+1
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtH(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vtH(par_count)
            ELSEIF (ccmd == 9) THEN ! gaussians describing the HB
               ! more then 1 gaussian can be used to describe the HB,
               ! but the standard case is to use just one gaussian
               ! and the default is to use the first gaussian in the list
               ! that will be passed through the clusterfile
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  IF (par_count .ge. MAXPARS) STOP "*** Too many HB clusters specified on command line. ***"
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vtghb(par_count)
                  par_count = par_count + 1
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

      !!!! CHECK for mandatory parameters !!!!!!
      ! The user, at least, must specifies what is what
      IF (.NOT.(errdef.EQ.3)) THEN
         WRITE(*,*) "", errdef
         WRITE(*,*) " Error: insert the hydrongen, donor and acceptor species! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! dopamm will be true just specifying the -gf flug
      IF (dopamm) THEN
         IF (clusterfile.EQ."NULL") THEN
            ! the user did something wrong in the GM specifications
            WRITE(*,*) ""
            WRITE(*,*) " Error: insert the file containing the gaussians parameters! "
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (vtghb(1).EQ.-1) THEN
            ! the user didn't specify indices of clusters that
            ! describe the HB into the PAMM framework.
            ! The first gaussian in the list will be used.
            ALLOCATE(vghb(1))
            vghb(1)=1
            nghb=1
         ENDIF
      ENDIF
      !!!! END CHECK .. we needed just few things !!!

      ! Invert the cell. Needed for PBC.
      CALL invmatrix(3,cell,icell)

      ! Loop over the input trajectory      
      ts=0 ! timestep
      DO
         ts=ts+1
         ! unit=5  --> read frome the standard input
         CALL xyz_read(5,natoms,header,labels,positions,endf)
         IF(endf<0) EXIT ! time to go!
         IF (ts==1) THEN ! first step, must allocate stuff
            ! inizialize the mask types and Gaussian structures
            IF (dopamm) THEN
               ! PAMM mode
               ! Read gaussian parameters from the gaussian file
               OPEN(UNIT=12,FILE=clusterfile,STATUS='OLD',ACTION='READ')
               ! read the gaussian model informations from a file.
               CALL readclusters(12,nk,clusters)
               CLOSE(UNIT=12)
               ALLOCATE(spa(nk,natoms), spd(nk,natoms), sph(nk,natoms))
               ALLOCATE(sa(natoms), sd(natoms), sh(natoms), pnks(nk))
               spa=0.0d0
               spd=0.0d0
               sph=0.0d0
            ENDIF
            ALLOCATE(masktypes(natoms))
            masktypes=TYPE_NONE
!!!!!!! should I put this part out of the if and reset the every step???
!!!!!!! this will make sense only if the order atoms in the input stream changes
!!!!!!! between different timestep
            DO i=1,natoms
               ! set the mask using BITWISE OR OPERATOR
               IF(testtype(labels(i),vtH)) masktypes(i)=IOR(masktypes(i),TYPE_H)
               IF(testtype(labels(i),vtdon)) masktypes(i)=IOR(masktypes(i),TYPE_DONOR)
               IF(testtype(labels(i),vtacc)) masktypes(i)=IOR(masktypes(i),TYPE_ACCEPTOR)
            ENDDO
         ENDIF
         IF (nptm .or. cell(1,1) == 0.0d0) THEN 
            ! NPT mode: this means ariable cell!
            ! Try to read the cell parameters the header in input stream
            READ(header, *) dummyc,dummyc,cell(1,1),cell(2,2),cell(3,3)
            CALL invmatrix(3,cell,icell)
         END IF
         
         IF (dopamm) THEN
            ! every step the collective variables values has to be reinitialized
            sph=0.0d0
            spa=0.0d0
            spd=0.0d0
         ENDIF
         
         ! here is the core of hbpamm
         ! we have to generate all the possible [nu,mu,r]
         DO ih=1,natoms ! loop over hydrogen types
            IF (IAND(masktypes(ih),TYPE_H).EQ.0) CYCLE
            DO id=1,natoms ! loop over donor types
               IF (IAND(masktypes(id),TYPE_DONOR).EQ.0 .OR. ih.EQ.id) CYCLE
               ! Get the distance D-H
               CALL pbcdist(cell,icell,positions(:,ih),positions(:,id),rdh)
               IF(rdh.GT.mucutoff) CYCLE  ! if the D-H distances is greater than the pre-defined
                                          ! cutoff, we do not need to check A-H
               DO ia=1,natoms ! loop over acceptor types
                  IF (IAND(masktypes(ia),TYPE_ACCEPTOR).EQ.0 &
                     .OR. (ia.EQ.id).OR.(ia.EQ.ih)) CYCLE
                  ! Get the distance acceptor-hydrogen
                  CALL pbcdist(cell,icell,positions(:,ih),positions(:,ia),rah)
                  x(2)=rdh+rah ! calculation of mu
                  IF(x(2).GT.mucutoff) CYCLE ! cutoff check in mu
                  x(1)=rdh-rah ! calculation of nu
                  ! Calculate the distance donor-acceptor (R)
                  CALL pbcdist(cell,icell,positions(:,id),positions(:,ia),x(3))

                  wfactor=1.0d0
                  IF(weighted) wfactor=1.0d0/((x(2)-x(1))*(x(2)+x(1))*x(3))

                  IF (dopamm) THEN
                     ! PAMM mode. We apply here the gaussian mixture model.
                     ! call lipbamm and compute the PAMM probability
                     CALL pamm_p(x, pnks, nk, clusters, alpha)
!!!!!! Where do I put the wfactor?? 
!!!!!! It should be passed inside pamm_p as in the previus GetP
                      
                     ! CALL GetP(3,x,wfactor,alpha,Nk,clusters,pks,pnks)
                     ! cumulate the probabilities conditional probabilieties
                     sph(:,ih) = sph(:,ih) + pnks(:)
                     spa(:,ia) = spa(:,ia) + pnks(:)
                     spd(:,id) = spd(:,id) + pnks(:)
                  ELSE
                     ! Pre-PAMM mode : write out x=[nu,mu,r] and the associated wfactor
                     WRITE(*,"(3(A1,ES21.8E4))",ADVANCE = "NO")  " ",x(1)," ",x(2)," ",x(3)
                     IF (weighted) WRITE(*,"(A1,ES21.8E4)",ADVANCE = "NO") " ", wfactor
                     write(*,*) "" ! go to a new line
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         IF(dopamm)THEN
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
            positions(1,:) = sh(:)
            positions(2,:) = sd(:)
            positions(3,:) = sa(:)
            CALL xyz_write(6,natoms,header,labels,positions)
         ENDIF
      ENDDO
      ! end the loop over the trajectory

      DEALLOCATE(positions)
      DEALLOCATE(labels)
      DEALLOCATE(masktypes)
      IF (dopamm) THEN
         DEALLOCATE(clusters,pks,pnks)
         DEALLOCATE(spa, sph, spd, sa, sh, sd, vghb)
      ENDIF

      CONTAINS

         SUBROUTINE helpmessage
            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: hbpamm -ta A1,A2,... -td D1,D2,... -th H1,H2,... "
            WRITE(*,*) "                [-h] [-l lx,ly,lz] [-w] [-npt] [-ct mucutoff] "
            WRITE(*,*) "                [-gf clusterfile] [-ghb 1,2,..] [-a smoothing_factor] "
            WRITE(*,*) "                 < input.xyz > output "
            WRITE(*,*) ""
            WRITE(*,*) " Description:  "
            WRITE(*,*) " hbpamm analyzes hydrogen-bonding patterns from simulation data"
            WRITE(*,*) " given in XYZ format. It can 1. pre-process a trajectory to compute the "
            WRITE(*,*) " local descriptors to be fed to the mixture model code, or "
            WRITE(*,*) " 2. use the parameters from a previously-fitted mixture model to "
            WRITE(*,*) " perform the actual analysis. "
            WRITE(*,*) ""
            WRITE(*,*) " Common options: "
            WRITE(*,*) "   -h                   : Print this message "
            WRITE(*,*) "   -l  lx,ly,lz         : Box lenghts (orthorombic box) "
            WRITE(*,*) "   -ta A1,A2,...        : Namelist of acceptors "
            WRITE(*,*) "   -td D1,D2,...        : Namelist of donors "
            WRITE(*,*) "   -th H1,H2,...        : Namelist of hydrogens "
            WRITE(*,*) "   -ct mucutoff         : Ignore DHA triplets with d(DH)+d(AH)>cutoff "
            WRITE(*,*) "   -npt                 : NPT mode. read cell date from the XYZ header "
            WRITE(*,*) "                          Format: # CELL: axx ayy azz           "
            WRITE(*,*) "   -w                   : Computes a weight for each DHA triplet "
            WRITE(*,*) "                          to account for the uniform-density phase space volume "
            WRITE(*,*) ""
            WRITE(*,*) " Pre-processing options: "
            WRITE(*,*) " Without further options, hbanalysis will just print  "
            WRITE(*,*) " nu=d(DH)-d(AH),mu=d(DH)+d(AH),r=d(DA) for each DHA triplet found "
            WRITE(*,*) " in the XYZ input. "
            WRITE(*,*) ""
            WRITE(*,*) " Mixture model analysis options: "
            WRITE(*,*) "   -gf Gaussians_file   : Activates the analysis and specifies the file "
            WRITE(*,*) "                          containing Gaussian clusters data "
            WRITE(*,*) "   -gh 1,2,...          : Indices of the gaussian(s) that describe the HB [default:1]"
            WRITE(*,*) "   -a  smoothing_factor : Apply a smoothing factor to the Gaussian model [default:1]"
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

      END PROGRAM hbpamm
