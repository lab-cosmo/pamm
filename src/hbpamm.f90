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
      INTEGER, ALLOCATABLE, DIMENSION(:) :: masktypes ! mask to define what is what
      CHARACTER(LEN=1024) :: filename, clusterfile 
      CHARACTER(LEN=1024) :: cmdbuffer

      !system parameters
      INTEGER natoms
      DOUBLE PRECISION cell(3,3), icell(3,3)
      DOUBLE PRECISION alpha, mucutoff, zeta
      
      ! Gaussian-mixture model
      INTEGER Nk ! number of gaussians in the mixture
      TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:) :: clusters
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pnks
      ! collective variables describing the HB
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: spa, spd, sph, sad
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: sa, sd, sh
      ! point in the high-dimensional description
      DOUBLE PRECISION, DIMENSION(3) :: x  ! [nu,mu,r]
      DOUBLE PRECISION rah,rdh,wfactor  ! distances and weight 
      INTEGER ia,id,ih ! simples indexes
      ! for a faster reading
      ! counters
      INTEGER i,ts
      LOGICAL convert,dopamm,dosad,nptm,weighted
      INTEGER delta
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: positions 
      INTEGER nghb ! number of gaussians describing the HB
      INTEGER, DIMENSION(MAXPARS) :: vghb ! indexes of the gaussians describing the HB
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: labels
      CHARACTER(LEN=1024) :: header, dummyc
      ! PARSER
      INTEGER ccmd     ! parser index
      INTEGER isep1, isep2, par_count  ! temporary indices for parsing command line arguments
      CHARACTER(LEN=4), DIMENSION(MAXPARS) :: vta,vtd,vth
      INTEGER endf
      
!!!!!!! Default value of the parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      vta(:)        = "NULL"
      vtd(:)        = "NULL"
      vth(:)        = "NULL"
      filename      = "NULL"
      clusterfile   = "NULL"
      vghb          = 0        
      nk            = -1        ! number of gaussians in the mixture
      mucutoff      = 5.0d0     ! cutoff in the mu coordinates to save cputime
      ccmd          = 0         ! no parameters specified
      delta         = 1         ! read every step 
      cell          = 0.0d0  
      alpha         = 1.0d0
      zeta          = 0.0d0
      convert       = .false.
      dopamm        = .false.
      dosad         = .false.
      nptm          = .false.
      weighted      = .false.   ! don't use wfactor by default
      endf          = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!! Command line parser !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
         ELSEIF (cmdbuffer == "-ev") THEN ! delta
            ccmd = 10
         ELSEIF (cmdbuffer == "-z") THEN ! smoothing factor, alpha
            ccmd = 11
         ELSEIF (cmdbuffer == "-sad") THEN ! hb lifetime statistics
            dosad = .true.
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
               isep1 = 0                  
               DO WHILE (index(cmdbuffer(isep1+1:), ',') > 0)
                  isep2 = index(cmdbuffer(isep1+1:), ',') + isep1
                  READ(cmdbuffer(isep1+1:isep2-1),*) cell(par_count,par_count)
                  par_count = par_count + 1
                  isep1=isep2
               ENDDO
               READ(cmdbuffer(isep1+1:),*) cell(par_count,par_count)
            ELSEIF (ccmd == 4) THEN ! smoothing factor, alpha
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 11) THEN ! zero factor, zeta
               READ(cmdbuffer,*) zeta
            ELSEIF (ccmd == 5) THEN ! cutoff in the mu=d(AH)+d(dH) distance
               READ(cmdbuffer,*) mucutoff
            ELSEIF (ccmd == 6) THEN ! acceptor types
               par_count = 1
               isep1 = 0                  
               DO WHILE (index(cmdbuffer(isep1+1:), ',') > 0)
                  IF (par_count .ge. MAXPARS) &
                     STOP &
         "*** Too many acceptor types specified on command line. ***"
                  isep2 = index(cmdbuffer(isep1+1:), ',') + isep1
                  READ(cmdbuffer(isep1+1:isep2-1),*) vta(par_count)
                  par_count = par_count + 1
                  isep1=isep2
               ENDDO
               READ(cmdbuffer(isep1+1:),*) vta(par_count)
            ELSEIF (ccmd == 7) THEN ! donor types
               par_count = 1
               isep1 = 0                  
               DO WHILE (index(cmdbuffer(isep1+1:), ',') > 0)
                  IF (par_count .ge. MAXPARS) &
                     STOP &
       "*** Too many donor types specified on command line. ***"
                  isep2 = index(cmdbuffer(isep1+1:), ',') + isep1
                  READ(cmdbuffer(isep1+1:isep2-1),*) vtd(par_count)
                  par_count = par_count + 1
                  isep1=isep2
               ENDDO
               READ(cmdbuffer(isep1+1:),*) vtd(par_count)
            ELSEIF (ccmd == 8) THEN ! hydrogen types
               par_count = 1
               isep1 = 0                  
               DO WHILE (index(cmdbuffer(isep1+1:), ',') > 0)
                  IF (par_count .ge. MAXPARS) &
                     STOP &
          "*** Too many hydrogen types specified on command line. ***"
                  isep2 = index(cmdbuffer(isep1+1:), ',') + isep1
                  READ(cmdbuffer(isep1+1:isep2-1),*) vth(par_count)
                  par_count = par_count + 1
                  isep1=isep2
               ENDDO
               READ(cmdbuffer(isep1+1:),*) vth(par_count)
            ELSEIF (ccmd == 9) THEN ! gaussians describing the HB
               ! more then 1 gaussian can be used to describe the HB,
               ! but the standard case is to use just one gaussian
               ! and the default is to use the first gaussian in the list
               ! that will be passed through the clusterfile               
               par_count = 1
               isep1 = 0                  
               DO WHILE (index(cmdbuffer(isep1+1:), ',') > 0)
                  IF (par_count .ge. MAXPARS) &
                     STOP &
            "*** Too many HB clusters specified on command line. ***"
                  isep2 = index(cmdbuffer(isep1+1:), ',') + isep1
                  READ(cmdbuffer(isep1+1:isep2-1),*) vghb(par_count)
                  par_count = par_count + 1
                  isep1=isep2
               ENDDO
               READ(cmdbuffer(isep1+1:),*) vghb(par_count)
            ELSEIF (ccmd == 10) THEN ! delta
               READ(cmdbuffer,*) delta
            ENDIF
         ENDIF
      ENDDO
!!!!!!!!!!! END PARSER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!! CHECK for mandatory parameters !!!!!!
      ! The user, at least, must specifies what is what
      IF ((vta(1).EQ."NULL").AND.(vtd(1).EQ."NULL").AND.(vth(1).EQ."NULL")) THEN
         WRITE(*,*) " Error: insert the hydrongen, donor and acceptor species! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! dopamm will be true just specifying the -gf flag
      IF (dopamm) THEN
         IF (clusterfile.EQ."NULL") THEN
            ! the user did something wrong in the GM specifications
            WRITE(*,*) " Error: insert the file containing the gaussians parameters! "
            CALL helpmessage
            CALL EXIT(-1)
         ENDIF
      ENDIF
      IF (dosad .and. .not. dopamm) THEN
         WRITE(*,*) " Error: cannot compute lifetime statistics without cluster data! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF
      ! End of the checks.

      ! Invert the cell. Needed for PBC.
      CALL invmatrix(3,cell,icell)

      ! Loop over the input trajectory      
      ts=0 ! timestep
      DO
         ts=ts+1
         IF (MODULO(ts,delta)/=0) THEN
            ! skip the frame
            READ(5,*,IOSTAT=endf) natoms
            IF(endf>0) STOP "*** Error occurred while reading file. ***"
            IF(endf<0) EXIT
            READ(5,'(A)',IOSTAT=endf) header
            IF(endf>0) STOP "*** Error occurred while reading file. ***"
            IF(endf<0) EXIT
            DO i=1,natoms
               READ(5,'(A)',IOSTAT=endf) header
            ENDDO
            IF(endf>0) STOP "*** Error occurred while reading file. ***"
            IF(endf<0) EXIT
            CYCLE
         ELSE
            ! unit=5  --> read frome the standard input
            CALL xyz_read(5,natoms,header,labels,positions,endf)
            IF(endf<0) EXIT ! time to go!
           
            IF (ALLOCATED(masktypes).and.SIZE(masktypes)/=natoms) DEALLOCATE(masktypes)
            IF (.not. ALLOCATED(masktypes)) THEN ! first step, must allocate stuff
               ! inizialize the mask types
               ALLOCATE(masktypes(natoms))
               masktypes=TYPE_NONE
               DO i=1,natoms
                  ! set the mask using BITWISE OR OPERATOR
                  IF(testtype(labels(i),vth)) masktypes(i)=IOR(masktypes(i),TYPE_H)
                  IF(testtype(labels(i),vtd)) masktypes(i)=IOR(masktypes(i),TYPE_DONOR)
                  IF(testtype(labels(i),vta)) masktypes(i)=IOR(masktypes(i),TYPE_ACCEPTOR)
               ENDDO
            ENDIF

            IF (dopamm .and. .not. ALLOCATED(spa)) THEN
               ! PAMM mode
               ! Read gaussian parameters from the gaussian file
               ! Non-periodic PAMM version
               OPEN(UNIT=12,FILE=clusterfile,STATUS='OLD',ACTION='READ')
               ! read the gaussian model informations from a file.
               CALL readclusters(12,nk,clusters)
               CLOSE(UNIT=12)
               IF(clusters(1)%D/=3) STOP "*** This analysis requires a 3D description of the HB! ***"
               ! The clusters should be sorted out if no specifications on the
               ! Gaussians describing the HB are given by the user
               IF (vghb(1).EQ.0) THEN
                  ! the user didn't specify aniything regarding the indices 
                  ! of the clusters that should describe the HB.
                  ! We will assume here there is just one cluster describing the HB
                  ! and that specific cluster will be the one closest to an hard coded
                  ! reference point
                  CALL sortclusters(nk, clusters)
                  vghb(1)=1
                  nghb=1
               ENDIF         
            ENDIF
            IF (ALLOCATED(spa) .and. SIZE(sa)/=natoms) THEN ! reallocate everything upon natoms change
               DEALLOCATE(spa,spd,sph,sa,sd,sh,pnks)
               IF (dosad) DEALLOCATE(sad)
            ENDIF
            IF (dopamm  .and. .not. ALLOCATED(spa)) THEN
               ALLOCATE(spa(nk,natoms), spd(nk,natoms), sph(nk,natoms))
               ALLOCATE(sa(natoms), sd(natoms), sh(natoms), pnks(nk))
               spa=0.0d0
               spd=0.0d0
               sph=0.0d0
               IF (dosad) THEN
                  ALLOCATE(sad(natoms,natoms))
                  sad=0.0d0
               ENDIF
            ENDIF
            
            IF (nptm .or. cell(1,1) == 0.0d0) THEN 
               ! NPT mode: this means variable cell!
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
            IF (dosad) sad = 0.0d0
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
                        CALL pamm_p(x, pnks, nk, clusters, alpha, zeta)
                        ! sum the probabilities from the gaussians describing the HB
                        sph(:,ih) = sph(:,ih) + pnks(:)
                        spa(:,ia) = spa(:,ia) + pnks(:)
                        spd(:,id) = spd(:,id) + pnks(:)
                        IF (dosad) THEN
                           DO i=1,nghb
                             sad(ia,id) = sad(ia,id) + pnks(i)
                           ENDDO
                        ENDIF
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
               IF (dosad) THEN
                  DO id=1,natoms ! loop over donor types
                     IF (IAND(masktypes(id),TYPE_DONOR).EQ.0) CYCLE
                     DO ia=1,natoms ! loop over acceptor types
                        IF (IAND(masktypes(ia),TYPE_ACCEPTOR).EQ.0 &
                            .OR. (ia.EQ.id)) CYCLE
                        WRITE(*,'(A1,ES13.5E4)', ADVANCE="NO") " ", sad(ia,id) 
                     ENDDO
                  ENDDO
                  WRITE(*,*) ""
               ELSE
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
            ENDIF
         ENDIF
      ENDDO
      ! end of the loop over the trajectory

      DEALLOCATE(positions)
      DEALLOCATE(labels)
      DEALLOCATE(masktypes)
      IF (dopamm) THEN
         DEALLOCATE(clusters,pnks)
         DEALLOCATE(spa, sph, spd, sa, sh, sd)
      ENDIF

      CONTAINS

         SUBROUTINE helpmessage
            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: hbpamm -ta A1,A2,... -td D1,D2,... -th H1,H2,... "
            WRITE(*,*) "                [-h] [-l lx,ly,lz] [-w] [-npt] [-ct mucutoff] [-ev delta] "
            WRITE(*,*) "                [-gf clusterfile] [-ghb 1,2,..] [-a smoothing_factor] "
            WRITE(*,*) "                [-z zeta_factor] "
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
            WRITE(*,*) "   -ev delta            : Stride while reading data from the XYZ input "
            WRITE(*,*) "   -npt                 : NPT mode. read cell data from the XYZ header "
            WRITE(*,*) "                          Header format: # CELL: axx ayy azz           "
            WRITE(*,*) "   -w                   : Computes a weight for each DHA triplet "
            WRITE(*,*) "                          to account for the uniform-density phase space volume "
            WRITE(*,*) ""
            WRITE(*,*) " Pre-processing options: "
            WRITE(*,*) " Without further options, hbanalysis will just print  "
            WRITE(*,*) " nu=d(DH)-d(AH),mu=d(DH)+d(AH),r=d(DA) for each DHA triplet found "
            WRITE(*,*) " in the XYZ input. "
            WRITE(*,*) ""
            WRITE(*,*) " Mixture model analysis options: "
            WRITE(*,*) "   -gf Gaussians_file    : Activates the analysis and specifies the file "
            WRITE(*,*) "                           containing Gaussian clusters data "
            WRITE(*,*) "   -ghb 1,2,...          : Indices of the gaussian(s) that describe the HB [default:1]"
            WRITE(*,*) "   -a smoothing_factor : Apply a smoothing factor to the Gaussian model [default:1]"
            WRITE(*,*) "   -z zeta_factor : Probabilities below this threshold are counted as 'no cluster' [default:0]"
            WRITE(*,*) "   -sad                 : Computes and print HB statistics for each donor/acceptor  "
            WRITE(*,*) "                          pair. Will generate a HUGE output file! "
            WRITE(*,*) ""
         END SUBROUTINE helpmessage

         LOGICAL FUNCTION testtype(id,vtype)
            CHARACTER*4, INTENT(IN) :: id
            CHARACTER*4, DIMENSION(MAXPARS), INTENT(IN) :: vtype
            INTEGER i
            testtype=.false.
            DO i=1,MAXPARS
               IF(trim(id).EQ.trim(vtype(i)))THEN
                  testtype=.true.
                  EXIT
               ENDIF
            ENDDO
         END FUNCTION

      SUBROUTINE sortclusters(nk,clusters)
         ! Sort the gaussians according to their distance from
         ! an hard-coded reference point
         !
         ! Args:
         !    nk: number of gaussian clusters
         !    clusters: array containing gaussians parameters

         INTEGER, INTENT(IN) :: nk
         TYPE(gauss_type), DIMENSION(nk), INTENT(INOUT) :: clusters

         TYPE(gauss_type) tmpgauss
         DOUBLE PRECISION distances(nk),prif(3),tmpdistance
         INTEGER j,i
         LOGICAL :: swapped = .TRUE.
         
         ! Set the reference to the HB cluster found in a b3lyp+vdW calculation
         prif(1)=-0.82d0
         prif(2)=2.82d0
         prif(3)=2.74d0
         ! calculate the distances of the means from the reference
         DO i=1,nk
            distances(i)=dot_product(clusters(i)%mean-prif,clusters(i)%mean-prif)
         ENDDO
         ! now we can sort using the distances
         ! will use bubble sort
         DO j=nk-1,1,-1
            swapped = .FALSE.
            DO i = 1, j
               IF (distances(i) > distances(i+1)) THEN
                  tmpdistance=distances(i)
                  distances(i)=distances(i+1)
                  distances(i+1)=tmpdistance
                  ! swap the clusters
                  tmpgauss=clusters(i)
                  clusters(i)=clusters(i+1)
                  clusters(i+1)=tmpgauss
                  swapped = .TRUE.
               END IF
            END DO
            IF (.NOT. swapped) EXIT
         ENDDO
      END SUBROUTINE sortclusters

      END PROGRAM hbpamm
