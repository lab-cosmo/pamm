! This file contain the main program for the Gaussian Mixture Model
! clustering. Starting from a set of 3D data points it will return the
! Nk gaussians better describing the clusters.
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
! EXPRESS OR IMPLIED, INCLUDIng BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRIngEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISIng FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALIngS IN THE SOFTWARE.
!
! Functions:
!    GetNlines: Get the number of lines from a file.
!    generatefromscratch: Iniatialize from scratch. Guess starting values for
!                         the gaussians parameters.
!    readgaussfromfile: Get the gaussian paramters from file
!    helpmessage: Banner containing the help message
!    ordergaussians: Order gaussians

      PROGRAM getmodel
         USE mixture
      IMPLICIT NONE

      CHARACTER*1024 :: outputfile                            ! The output file prefix
      DOUBLE PRECISION, ALLOCATABLE :: prif(:)                             ! Reference point needed to find out the important cluster
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: distmm ! similarity matrix
      DOUBLE PRECISION, DIMENSION(3) :: diff                  ! temp vector used to store distances

      INTEGER D
      INTEGER npc                                             ! number of point in a cluster. Used to define the gaussians covariances
      INTEGER Nk                                              ! Number of gaussians in the mixture
      INTEGER delta                                           ! Number of data to skeep (for a faster calculation)
      INTEGER Nlines                                          ! Number of lines of the input data file
      INTEGER nsamples,nsamplesin                             ! Total number points
      INTEGER nminmax                                         ! Number of samples extracted using minmax
      DOUBLE PRECISION :: dij,dmin
      INTEGER jmax
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma2, wj, probnmm
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npvoronoi,iminij,pnlist,nlist
      INTEGER seed     ! seed for the random number generator

      ! Array of Gaussians containing the gaussians parameters
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters
      ! Array containing the logarithm of the fractions (Pk) for each gaussian
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: lpks
      ! Array containing the input data pints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vwad , Ymm
      ! quick shift, roots and path to reach the root (used to speedup the calculation)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: idxroot,qspath
      ! Array containing the nk probalities for each of nsamples points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pnk ! responsibility matrix

      ! PARSER
      CHARACTER*1024 :: cmdbuffer   ! String used for reading text lines from files
      INTEGER ccmd                  ! Index used to control the input parameters
      LOGICAL verbose ! flag for verbosity
      LOGICAL weighted ! flag for using weigheted data
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string
	  DOUBLE PRECISION vpar(5)
	  DOUBLE PRECISION dummyr1,tau,meannew(3), tau2

      INTEGER i,j,k,counter,dummyi1 ! Counters and dummy variable


      !!!!!!! Iniatialze the parameters !!!!!!!
      outputfile="out"
      ccmd=0              ! no parameters specified
      delta=1             ! read every point
      Nk=0                ! number of gaussians
      nsamplesin=-1       ! total number of points
      nminmax=-1          ! number of samples extracted with minmax
      seed=12345          ! seed for the random number generator
      tau=-1              ! quick shift cut-off
      verbose = .false.   ! no verbosity
      weighted= .false.   ! don't use the weights
      D=-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-o") THEN      ! output file
            ccmd = 2
         ELSEIF (cmdbuffer == "-d") THEN      ! dimensionality
            ccmd = 9
         ELSEIF (cmdbuffer == "-ev") THEN     ! delta
            ccmd = 3
         ELSEIF (cmdbuffer == "-seed") THEN   ! seed for the random number genarator
            ccmd = 4
         ELSEIF (cmdbuffer == "-tau") THEN    ! threshold to differentiate different clusters
            ccmd = 5
         ELSEIF (cmdbuffer == "-maxdata") THEN ! N samples (total)
            ccmd = 6
         ELSEIF (cmdbuffer == "-nkde") THEN ! N samples extracted with minmax
            ccmd = 7
         ELSEIF (cmdbuffer == "-ref") THEN    ! point from wich calculate the distances to order
            ccmd = 8                          ! the gaussians in the output
         ELSEIF (cmdbuffer == "-w") THEN      ! use weights
            weighted = .true.
         ELSEIF (cmdbuffer == "-v") THEN      ! verbosity flag
            verbose = .true.
         ELSEIF (cmdbuffer == "-h") THEN      ! help flag
            CALL helpmessage
            CALL EXIT(-1)
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " Wrong usage. Insert the right parameters!"
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 2) THEN ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 9) THEN
               READ(cmdbuffer,*) D
               ALLOCATE(prif(D))
               prif=0.0d0
               prif(1)=-1.0d0
            ELSEIF (ccmd == 3) THEN ! delta
               READ(cmdbuffer,*) delta
            ELSEIF (ccmd == 4) THEN ! read the seed
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 5) THEN ! cut-off for quick-shift
               READ(cmdbuffer,*) tau
            ELSEIF (ccmd == 6) THEN
               READ(cmdbuffer,*) nsamplesin
            ELSEIF (ccmd == 7) THEN
               READ(cmdbuffer,*) nminmax
            ELSEIF (ccmd == 8) THEN ! vrif,wrif,dADrif
               IF (D<0) error STOP "Dimensionality must be given before the reference point. "
               par_count = 1
               commas(1) = 0   !TODO FIX EVERYWHERE! THERE IS NO NEED TO HAVE AN ARRAY FOR THE COMMAS!
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) prif(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) prif(par_count)
            ENDIF
         ENDIF
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL SRAND(seed) ! initialize the random number generator

      ! Check the input parameters
      IF (ccmd.EQ.0) THEN
         ! ccmd==0 : no parameters inserted
         WRITE(*,*) ""
         WRITE(*,*) " Wrong usage. Insert the right parameters!"
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! The dimensionalty can't be hard coded by default
      IF (D.EQ.-1) THEN
         WRITE(*,*) ""
         WRITE(*,*) " Wrong usage. Insert the dimensionality!"
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! Nminmax should be smaller than samples
      IF ((nminmax.NE.-1).AND.(nsamplesin.NE.-1).AND.(nminmax.GT.nsamplesin)) THEN
         ! ccmd==0 : no parameters inserted
         WRITE(*,*) ""
         WRITE(*,*) " Wrong usage. You probably misunderstood the meaning of nminmax! "
         WRITE(*,*) " Nminmax should be smaller than samples. "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! get the data from the standard input
      CALL readinput(D,nsamples,nsamplesin,delta,vwad,weighted,wj)

      ! If not specified, the number voronoi polyhedras
      ! are set to the square of the total number of points
      IF (nminmax.EQ.-1) nminmax=sqrt(float(nsamples))

	  ALLOCATE(iminij(nsamples))
	  ALLOCATE(pnlist(nminmax+1),nlist(nsamples))
	  ALLOCATE(Ymm(D,nminmax),npvoronoi(nminmax),probnmm(nminmax),sigma2(nminmax))
	  ALLOCATE(idxroot(nminmax),qspath(nminmax),distmm(nminmax,nminmax))

      ! Extract nminmax points on which the kernel density estimation is to be
      ! evaluated. Also partitions the nsamples points into the Voronoi polyhedra
      ! of the sampling points.
      IF(verbose) THEN
         WRITE(*,*) "NSamples: " , nsamples
         WRITE(*,*) "Selecting ", nminmax, " points using MINMAX"
      ENDIF

      CALL getvoronoi(nsamples,nminmax,vwad,Ymm,npvoronoi,iminij)

      ! Generate the neighbour list
      IF(verbose) write(*,*) "Generating neighbour list"
      CALL getnlist(nsamples,nminmax,npvoronoi,iminij, pnlist,nlist)

      ! Definition of the similarity matrix between Voronoi centers
      distmm=0.0d0
      sigma2=1e10 ! adaptive sigma of the kernel density estimator
      IF(verbose) write(*,*) "Computing similarity matrix"
      DO i=1,nminmax
         DO j=1,i-1
            ! distance between two voronoi centers
            distmm(i,j) = dot_product( Ymm(:,i) - Ymm(:,j) , Ymm(:,i) - Ymm(:,j) )
            if (distmm(i,j) < sigma2(i)) sigma2(i) = distmm(i,j)
            ! the symmetrical one
            distmm(j,i) = distmm(i,j)
            if (distmm(i,j) < sigma2(j)) sigma2(j) = distmm(i,j)
         ENDDO
      ENDDO

      IF(tau.EQ.-1)THEN
         ! tau set to 3*<sig>
         tau2=9.0d0*SUM(sigma2)/nminmax
         tau=dsqrt(tau2)
      ENDIF
      tau2=tau*tau ! we always work with squared distances....

      IF(verbose) write(*,*) "Computing kernel density on reference points."
      IF(verbose) write(*,*) "Tau : ", tau
      ! computes the KDE on the Voronoi centers using the neighbour list
      probnmm = 0.0d0
      DO i=1,nminmax
         DO j=1,nminmax

            ! check if the polyhedra are too far away
            IF (distmm(i,j)/sigma2(j)>25.0d0) CYCLE

            ! cycle just inside the polyhedra thanx to the neighbour list
            DO k=pnlist(j)+1,pnlist(j+1)
               probnmm(i)=probnmm(i)+ wj(nlist(k))* &
                          fkernel(sigma2(j),Ymm(:,i),vwad(:,nlist(k)))
            ENDDO
         ENDDO
      ENDDO

	  IF(verbose) write(*,*) "Running quick shift"

	  idxroot=0
	  ! Start quick shift
	  DO i=1,nminmax
	     IF(idxroot(i).NE.0) CYCLE
	     !idxroot(i)=i
	     qspath=0
	     qspath(1)=i
	     counter=1
	     DO WHILE(qspath(counter).NE.idxroot(qspath(counter)))
	        idxroot(qspath(counter))= &
	            GethigherNN(nminmax,qspath(counter),tau2,probnmm,distmm)
            IF(idxroot(idxroot(qspath(counter))).NE.0) EXIT
            counter=counter+1
            qspath(counter)=idxroot(qspath(counter-1))
	     ENDDO
	     DO j=1,counter
	        idxroot(qspath(j))=idxroot(idxroot(qspath(counter)))
	     ENDDO
	  ENDDO

	  IF(verbose) write(*,*) "Writing out"
	  qspath=0
	  qspath(1)=idxroot(1)
	  Nk=1
	  OPEN(UNIT=11,FILE=trim(outputfile)//".clusters",STATUS='REPLACE',ACTION='WRITE')
	  DO i=1,nminmax
	     ! write out the clusters
	     dummyi1=0
	     DO k=1,Nk
	        IF(idxroot(i).EQ.qspath(k))THEN
	           dummyi1=k
	           EXIT
	        ENDIF
	     ENDDO
	     IF(dummyi1.EQ.0)THEN
	        Nk=Nk+1
	        qspath(Nk)=idxroot(i)
	        dummyi1=Nk
	     ENDIF
	     WRITE(11,"(3(A1,ES15.4E4))",ADVANCE = "NO") " ", Ymm(1,i)," ", Ymm(2,i)," ", Ymm(3,i)
	     WRITE(11,"(A1,I4,A1,ES15.4E4)") " ", dummyi1 , " ", probnmm(i)
	  ENDDO
	  CLOSE(UNIT=11)

	  ! now qspath contains the indexes of Nk gaussians

	  !! covariance
	  ALLOCATE(clusters(Nk),lpks(Nk))
	  ! calculate the gaussians covariance from the data in the clusters
	  DO k=1,Nk
	     IF (.not.(ALLOCATED(clusters(k)%mean))) ALLOCATE(clusters(k)%mean(D))
        IF (.not.(ALLOCATED(clusters(k)%cov)))  ALLOCATE(clusters(k)%cov(D,D))
        IF (.not.(ALLOCATED(clusters(k)%icov))) ALLOCATE(clusters(k)%icov(D,D))
	     clusters(k)%mean=Ymm(:,qspath(k))
	     meannew = 0.0d0
	     npc=0
	     clusters(k)%cov = 0.0d0
	     DO i=1,nminmax
	        IF(idxroot(i).NE.qspath(k)) CYCLE
	        ! to improve using the also the information from the KDE
	        meannew = meannew + Ymm(:,i)
	        clusters(k)%cov(1,1) = clusters(k)%cov(1,1) + Ymm(1,i) * Ymm(1,i)
	        clusters(k)%cov(1,2) = clusters(k)%cov(1,2) + Ymm(1,i) * Ymm(2,i)
	        clusters(k)%cov(1,3) = clusters(k)%cov(1,3) + Ymm(1,i) * Ymm(3,i)
	        clusters(k)%cov(2,2) = clusters(k)%cov(2,2) + Ymm(2,i) * Ymm(2,i)
	        clusters(k)%cov(2,3) = clusters(k)%cov(2,3) + Ymm(2,i) * Ymm(3,i)
	        clusters(k)%cov(3,3) = clusters(k)%cov(3,3) + Ymm(3,i) * Ymm(3,i)
	        npc=npc+1
	     ENDDO
	     meannew = meannew/npc
	     clusters(k)%cov = clusters(k)%cov /npc
	     clusters(k)%cov(1,1) = clusters(k)%cov(1,1) - meannew(1)*meannew(1)
	     clusters(k)%cov(1,2) = clusters(k)%cov(1,2) - meannew(1)*meannew(2)
	     clusters(k)%cov(1,3) = clusters(k)%cov(1,3) - meannew(1)*meannew(3)
	     clusters(k)%cov(2,2) = clusters(k)%cov(2,2) - meannew(2)*meannew(2)
	     clusters(k)%cov(2,3) = clusters(k)%cov(2,3) - meannew(2)*meannew(3)
	     clusters(k)%cov(3,3) = clusters(k)%cov(3,3) - meannew(3)*meannew(3)
	     clusters(k)%cov(2,1) = clusters(k)%cov(1,2)
	     clusters(k)%cov(3,1) = clusters(k)%cov(1,3)
	     clusters(k)%cov(3,2) = clusters(k)%cov(2,3)
	     lpks(k)=LOG(FLOAT(npc)/nminmax)
	  ENDDO
	  ! write gaussians

	  ! oreder gaussians from the closest to the reference point
	  CALL ordergaussians(D,Nk,clusters,lpks,prif)

!	  OPEN(UNIT=11,FILE=trim(outputfile)//".gauss",STATUS='REPLACE',ACTION='WRITE')
!	  !WRITE(11,*) "# Mean-Shift output (Sig,Err Conv, Err Clusters): " , dsqrt(twosig2/2.0d0), errc, errclusters
!	  WRITE(11,"(A31)",ADVANCE="NO") "# Quick Shift GM output. Ntot: "
!	  WRITE(11,"(I12,A12,I11)",ADVANCE="NO") nsamples," , NVoroni: ",nminmax
!	  WRITE(11,"(A8,ES15.4E4)") " , Tau: ", tau
!	  WRITE(11,*) "# mean cov pk"
!	  WRITE(11,*) Nk
!	  DO k=1,Nk
!	     CALL writegausstofile(11,clusters(k),lpks(k))
!	  ENDDO

!	  CLOSE(UNIT=11)

	  CALL writegaussianstofile(outputfile,D,nsamples,nminmax,tau,Nk,clusters,lpks)

	  DEALLOCATE(clusters,lpks)
	  DEALLOCATE(vwad,wj)
	  DEALLOCATE(idxroot,qspath,distmm)
	  DEALLOCATE(pnlist,nlist,iminij)
	  DEALLOCATE(Ymm,npvoronoi,probnmm,sigma2)

	  CALL EXIT(0)
	  ! end of the main programs



      ! functions and subroutines !!!!!!!!!!!!!!!!!!!!

      CONTAINS

         SUBROUTINE helpmessage
            ! Banner to print out for helping purpose
            !

            WRITE(*,*) ""
            WRITE(*,*) " USAGE: getmodel [-h] -d D [-w] [-o output] [-seed seedrandom] "
            WRITE(*,*) "                 [-tau tau] [-ev delta] [-nsamples NTot] [-nminmax Nminmax] "
            WRITE(*,*) "                 [-rif v,w,R] [-v] "
            WRITE(*,*) ""
            WRITE(*,*) " Clusterize the data and define a mixture of gaussians describing them. "
            WRITE(*,*) " It is mandatory to specify the dimensionality of the data and the data "
            WRITE(*,*) " file must be passed through the standard input in the format: "
            WRITE(*,*) " x11 x12 x13 ... x1D [w1] "
            WRITE(*,*) " x21 x22 x23 ... x2D [w2] "
            WRITE(*,*) ""
            WRITE(*,*) " For the other options a default is defined.  "
            WRITE(*,*) ""
            WRITE(*,*) "   -h                : Print this message "
            WRITE(*,*) "   -d D              : Dimensionality "
            WRITE(*,*) "   -w                : Reads weights for the sample points "
            WRITE(*,*) "   -o output         : Prefix for output files [out]. This will produce : "
            WRITE(*,*) ""
            WRITE(*,*) "                            output.clusters (clusterized data) "
            WRITE(*,*) "                            output.gauss (gaussians) "
            WRITE(*,*) ""
            WRITE(*,*) "   -seed seed        : Seed to initialize the random number generator. [12345]"
            WRITE(*,*) "   -tau tau          : Quick shift cutoff [automatic] "
            WRITE(*,*) "   -ev delta         : Stride reading data frome the file [1] "
            WRITE(*,*) "   -maxsamples nmax  : Maximum number of samples read from the input [read all]"
            WRITE(*,*) "   -nkde nkde        : Number of points to evaluate KDE [sqrt(nsamples)]"
            WRITE(*,*) "   -ref X            : Reference point for ordering the clusters [ (-1,0,0,...) ]"
            WRITE(*,*) "   -v                : Verbose output "
            WRITE(*,*) ""
         END SUBROUTINE helpmessage

         SUBROUTINE readinput(D,nsamples,naspamplesin,delta,vout,weighted,wj)
            INTEGER, INTENT(IN) :: D
            INTEGER, INTENT(OUT) :: nsamples
            INTEGER, INTENT(IN) :: naspamplesin
            INTEGER, INTENT(IN) :: delta
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: vout
            LOGICAL, INTENT(IN) :: weighted
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: wj


            INTEGER, PARAMETER :: nbuff = 100000
            DOUBLE PRECISION, DIMENSION(D,nbuff) :: vbuff
            DOUBLE PRECISION, DIMENSION(nbuff) :: wbuff
            INTEGER io_status,counter,i

            nsamples=0
            counter=0
            i=0
            ! if I don't allocate the vectors I get a segmentation fault
            ALLOCATE(vout(D,10))
            IF(weighted) ALLOCATE(wj(10))
            vout=0.0d0
            DO
               i=i+1
               IF(weighted) THEN
                  READ(5,*, IOSTAT=io_status) vbuff(:,counter+1),wbuff(counter+1)
               ELSE
                  READ(5,*, IOSTAT=io_status) vbuff(:,counter+1)
               ENDIF
               IF(io_status<0) EXIT
               IF(io_status>0) error STOP "*** Error occurred while reading file. ***"
               IF(MODULO(i,delta).NE.0) CYCLE
               counter=counter+1
               IF(nsamples.GT.-1)THEN
                  IF(nsamples+counter.EQ.naspamplesin) EXIT
               ENDIF
               IF(counter.EQ.nbuff) THEN
                  IF(weighted) CALL collapse1d(nsamples,counter,wj,wbuff(1:counter))
                  CALL collapsend(D,nsamples,counter,vout,vbuff(:,1:counter))
                  counter=0
                  vpar=0
               ENDIF
            END DO

            IF(counter>0) THEN
               IF(weighted) CALL collapse1d(nsamples,counter,wj,wbuff(1:counter))
               CALL collapsend(D,nsamples,counter,vout,vbuff(:,1:counter))
            ENDIF

            IF(.NOT.weighted)THEN
               ALLOCATE(wj(nsamples))
               wj=1.0d0
            ENDIF

         END SUBROUTINE readinput

         SUBROUTINE getvoronoi(nsamples,nminmax,vwad,Ymm,npvoronoi,iminij)
            ! Select nminmax points from nsamples using minmax and
            ! the voronoi polyhedra around them.
            !
            ! Args:
            !    nsamples: total points number
            !    nminmax: number of voronoi polyhedra
            !    vwda: array containing the data
            !    Ymm: array that will contain the voronoi centers
            !    npvoronoi: array cotaing the number of points inside each voroni polyhedra
            !    iminij: array containg to wich polyhedra every point belong to.

            INTEGER, INTENT(IN) :: nsamples
            INTEGER, INTENT(IN) :: nminmax
            DOUBLE PRECISION, DIMENSION(3,nsamples), INTENT(IN) :: vwad
            DOUBLE PRECISION, DIMENSION(3,nminmax), INTENT(OUT) :: Ymm
            INTEGER, DIMENSION(nminmax), INTENT(OUT) :: npvoronoi
            INTEGER, DIMENSION(nsamples), INTENT(OUT) :: iminij

            INTEGER i,j
            DOUBLE PRECISION :: diff(3)
            DOUBLE PRECISION :: dminij(nsamples), dij, dmax

            iminij=0
            Ymm=0.0d0
            npvoronoi=0
            ! choose randomly the first point
            Ymm(:,1)=vwad(:,int(RAND()*nsamples))
            dminij = 1.0d99
            iminij = 1
            DO i=2,nminmax
               dmax = 0.0d0
               DO j=1,nsamples
                  dij = dot_product( Ymm(:,i-1) - vwad(:,j) , &
                                     Ymm(:,i-1) - vwad(:,j) )
                  IF (dminij(j)>dij) THEN
                     dminij(j) = dij
                     iminij(j) = i-1 ! also keeps track of the Voronoi attribution
                  ENDIF
                  IF (dminij(j) > dmax) THEN
                     dmax = dminij(j)
                     jmax = j
                  ENDIF
               ENDDO
               Ymm(:,i) = vwad(:, jmax)
               IF(verbose .AND. (modulo(i,1000).EQ.0)) &
                  write(*,*) i,"/",nminmax
            ENDDO

            ! finishes Voronoi attribution
            DO j=1,nsamples
               dij = dot_product( Ymm(:,nminmax) - vwad(:,j) , &
                                  Ymm(:,nminmax) - vwad(:,j) )
               IF (dminij(j)>dij) THEN
                  dminij(j) = dij
                  iminij(j) = nminmax
               ENDIF
            ENDDO

            ! Number of points in each voronoi polyhedra
            npvoronoi=0.0d0
            DO j=1,nsamples
               npvoronoi(iminij(j))=npvoronoi(iminij(j))+1
            ENDDO
         END SUBROUTINE getvoronoi

         SUBROUTINE getnlist(nsamples,nminmax,npvoronoi,iminij, pnlist,nlist)
            ! Build a neighbours list: for every voronoi center keep track of his
            ! neighboroud that correspond to all the points inside the voronoi
            ! polyhedra.
            !
            ! Args:
            !    nsamples: total points number
            !    nminmax: number of voronoi polyhedra
            !    weights: array cotaing the number of points inside each voroni polyhedra
            !    iminij: array containg to wich polyhedra every point belong to
            !    pnlist: pointer to neighbours list
            !    nlist: neighbours list

            INTEGER, INTENT(IN) :: nsamples
            INTEGER, INTENT(IN) :: nminmax
            INTEGER, DIMENSION(nminmax), INTENT(IN) :: npvoronoi
            INTEGER, DIMENSION(nsamples), INTENT(IN) :: iminij
            INTEGER, DIMENSION(nminmax+1), INTENT(OUT) :: pnlist
            INTEGER, DIMENSION(nsamples), INTENT(OUT) :: nlist

            INTEGER i,j
            INTEGER :: tmpnidx(nminmax)
            DOUBLE PRECISION :: dij, dmax

            pnlist=0
            nlist=0
            tmpnidx=0

            ! pointer to the neighbourlist
            pnlist(1)=0
            DO i=1,nminmax
               pnlist(i+1)=pnlist(i)+npvoronoi(i)
               tmpnidx(i)=pnlist(i)+1  ! temporary array to use while filling up the neighbour list
            ENDDO

            DO j=1,nsamples
               i=iminij(j) ! this is the Voronoi center the sample j belongs to
               nlist(tmpnidx(i))=j ! adds j to the neighbour list
               tmpnidx(i)=tmpnidx(i)+1 ! advances the pointer
            ENDDO
         END SUBROUTINE getnlist

         INTEGER FUNCTION GethigherNN(nminmax,idx,tau,probnmm,distmm)
            ! Return the index of the closest point higher in P
            !
            ! Args:
            !    nminmax: number of voronoi polyhedra
            !    idx: current point
            !    tau: cut-off in the jump
            !    probnmm: density estimations
            !    distmm: distances matrix

            INTEGER, INTENT(IN) :: nminmax
            INTEGER, INTENT(IN) :: idx
            DOUBLE PRECISION, INTENT(IN) :: tau
            DOUBLE PRECISION, DIMENSION(nminmax), INTENT(IN) :: probnmm
            DOUBLE PRECISION, DIMENSION(nminmax,nminmax), INTENT(IN) :: distmm

            INTEGER j
            DOUBLE PRECISION dmin

            dmin=1.0d10
            GethigherNN=idx
	        DO j=1,nminmax
	           IF(probnmm(j)>probnmm(idx))THEN
	              IF((distmm(idx,j).LT.dmin) .AND. (distmm(idx,j).LT.tau))THEN
	                 dmin=distmm(idx,j)
	                 GethigherNN=j
	              ENDIF
	           ENDIF
	        ENDDO
	        !write(*,*) "jump distance", dmin
         END FUNCTION GethigherNN

		 DOUBLE PRECISION FUNCTION fkernel(sig2,vc,vp)
            ! Calculate the (non-normalized) gaussian kernel
            !
            ! Args:
            !    sig2: sig**2
            !    vc: voronoi center's vector
            !    vp: point's vector

            DOUBLE PRECISION, INTENT(IN) :: sig2
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: vc
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: vp

            fkernel=(1/dsqrt(dpigreco*sig2))*dexp(-dot_product(vc-vp,vc-vp)*0.5/sig2)
         END FUNCTION fkernel

      END PROGRAM getmodel
