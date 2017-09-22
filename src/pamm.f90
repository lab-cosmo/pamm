! This file contain the main program for the PAMM clustering in
! both PERIODIC and NON PERIODIC space.
! Starting from a set of data points in high dimension it will first perform
! a non-parametric partitioning of the probability density and return the
! Nk multivariate Gaussian/Von Mises distributions better describing the clusters.
! Can also be run in post-processing mode, where it will read tuples and
! classify them based the model file specified in input.
!
! Copyright (C) 2016, Piero Gasparotto, Robert Meissner and Michele Ceriotti
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
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
! TODO: (i) we need only to store Hiinv and normkernel for each grid point
!           -> remove Di(:), Hi(:), ... to save memory
!           -> if periodic data is used we probably need to save Hi(:) too

      PROGRAM pamm
      USE libpamm
      USE random
      IMPLICIT NONE

      CHARACTER(LEN=1024) :: outputfile, clusterfile            ! The output file prefix
      CHARACTER(LEN=1024) :: gridfile, neighfile                ! The output file prefix
      CHARACTER(LEN=1024) :: cmdbuffer, comment                 ! String used for reading text lines from files

      LOGICAL periodic                                          ! flag for using periodic data
      LOGICAL verbose                                           ! flag for verbosity
      LOGICAL fpost                                             ! flag for postprocessing
      LOGICAL weighted                                          ! flag for using weigheted data
      LOGICAL savevor, saveadj, saveidxs, readgrid              ! additional IN/OUT logical flags
      LOGICAL saveneigh, readneigh                              ! gabriel graph IN/OUT logical flags
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: gabriel           ! gabriel graph matrix
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: mergeornot

      INTEGER ccmd                                              ! Index used to control the PARSER input parameters
      INTEGER endf                                              ! end file state for reading in data from file
      INTEGER D                                                 ! Dimensionality of problem
      INTEGER Nk                                                ! Number of gaussians in the mixture
      INTEGER nsamples                                          ! Total number points
      INTEGER ngrid                                             ! Number of samples extracted using minmax
      INTEGER seed                                              ! seed for the random number generator
      INTEGER nmsopt                                            ! number of mean-shift optimizations of the cluster centers
      INTEGER nbootstrap                                        ! number of bootstrap cycles
      INTEGER gs                                                ! flag for gabriel clustering (and also number of neigh shells)
      INTEGER rndidx                                            ! random sample point index
      INTEGER nbssample                                         ! number of sample points used for voronoi in bootstrap
      INTEGER nbstot                                            ! accumulator for nbssample
      INTEGER isep1, isep2, par_count                           ! temporary indices for parsing command line arguments
      INTEGER i,j,jmax,k,nn,counter                             ! counters
      INTEGER dummyi1                                           ! dummy variables

      ! neighbor list number of points in voronoi, voronoi
      ! association, pointer, ..., sample point index of grid point
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ni, iminij, pnlist, nlist, idxgrid
      ! quick shift, roots and path to reach the root (used to speedup the calculation)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: idxroot, qspath
      ! macrocluster
      INTEGER, ALLOCATABLE, DIMENSION(:) :: macrocl,sortmacrocl,clustercenters
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ineigh,idmindist

      DOUBLE PRECISION normwj                                   ! accumulator for wj
      DOUBLE PRECISION lim, delta                               ! lim and delta for bisectioning in adaptive bandwidth estimation
      DOUBLE PRECISION tmppks,normpks                           ! variables to set GM covariances
      DOUBLE PRECISION thrpcl                                   ! parmeter controlling the merging of the outlier clusters
      DOUBLE PRECISION fpoints                                  ! use either a fraction of sample points
      DOUBLE PRECISION fspread                                  ! or a fraction of the global avg. variance
      DOUBLE PRECISION nlocal                                   ! local number of points
      DOUBLE PRECISION tune                                     ! tuning used in bisectioning to find flocal
      DOUBLE PRECISION qs                                    ! scaling factor for QS
      DOUBLE PRECISION kdecut2                                  ! squared cutoff for KDE
      DOUBLE PRECISION msw
      DOUBLE PRECISION alpha                                    ! cluster smearing
      DOUBLE PRECISION zeta                                     ! background for clustering
      DOUBLE PRECISION thrmerg                                  ! threshold for adjacency cluster merging
      DOUBLE PRECISION dummd1,dummd2                            ! dummy variables
      DOUBLE PRECISION lnK                                      ! logarithm of kernel

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma2     ! adaptive localizations
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: qscut2     ! adaptive quick shift cutoff
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wj, lwj    ! weight of each sample point
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wi, lwi    ! accumulator for wj in each voronoi
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wlocal     ! local weights around grid point for grid points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wlocal2    ! local weights around grid point for sample points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: flocal     ! local number of points around localized grid
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: mindist    ! distance to closest voronoi
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Di         ! local dimensionality
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: period     ! Periodic lenght in each dimension
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: prelerr    ! relative error of probability
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pabserr    ! absolute error of probability
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dij        ! distance vector between two points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: normkernel ! normalization for ln(K(Y))
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: logdetHi   ! logarithm of bandwidth matrix determinant
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: prob       ! probabilities at grid points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rgrid, msmu, tmpmsmu, pcluster, px, tmps2

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x, y     ! Array containing the input data and grid points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: distmm   ! similarity matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: probboot ! bootstrap probabilities
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Q        ! global covariance matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Qi       ! local covariance matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Qiinv    ! inversed local covariance matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Hi       ! bandwidth matrix

      ! heavy bandwidth matrix for kernel density estimation
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Hiinv ! inversed bandwidth matrices

      ! Array of Von Mises distributions
      TYPE(vm_type), ALLOCATABLE, DIMENSION(:) :: vmclusters
      ! Array of Gaussians containing the gaussians parameters
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters

!!!!!!! Default value of the parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      outputfile = "out"
      clusterfile = "NULL"
      gridfile = "NULL"
      fpost = .false.
      alpha = 1.0d0
      zeta = 0.0d0
      thrpcl = 0.0d0
      fpoints = 0.15d0       ! fraction of points to be used as standard
      fspread = -1.0d0       ! fraction of avg. variance to be used (default: off)
      ccmd = 0               ! no parameters specified
      Nk = 0                 ! number of gaussians
      nmsopt = 0             ! number of mean-shift refinements
      ngrid = -1             ! number of samples extracted with minmax
      seed = 12345           ! seed for the random number generator
      thrmerg = 0.8d0        ! merge different clusters
      qs = 1.0d0             ! quick shift cut-off
      gs = -1                ! don't use clustering based on gabriel graph
      verbose = .FALSE.      ! no verbosity
      weighted = .FALSE.     ! don't use the weights
      nbootstrap = 0         ! do not use bootstrap
      savevor  = .FALSE.     ! don't print out the Voronoi
      saveidxs = .FALSE.     ! don't save the indexes of the grid points
      saveadj = .FALSE.      ! save adjacency
      readgrid = .FALSE.     ! don't read the grid from the standard input
      saveneigh = .FALSE.    ! don't save gabriel graph
      readneigh = .FALSE.    ! don't read gabriel graph


      D=-1
      periodic=.FALSE.
      CALL random_init(seed) ! initialize random number generator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-a") THEN                ! cluster smearing
            ccmd = 1
         ELSEIF (cmdbuffer == "-o") THEN            ! output file
            ccmd = 2
         ELSEIF (cmdbuffer == "-gf") THEN           ! file containing Vn parmeters
            ccmd = 3
         ELSEIF (cmdbuffer == "-seed") THEN         ! seed for the random number genarator
            ccmd = 4
         ELSEIF (cmdbuffer == "-qs") THEN           ! scale to differentiate cluster in quickshift
            ccmd = 5
         ELSEIF (cmdbuffer == "-nms") THEN          ! mean-shift steps
            ccmd = 6
         ELSEIF (cmdbuffer == "-ngrid") THEN        ! N of grid points
            ccmd = 7
         ELSEIF (cmdbuffer == "-bootstrap") THEN    ! estimate error of kde using bootstrap
            ccmd = 8
         ELSEIF (cmdbuffer == "-d") THEN            ! dimensionality
            ccmd = 9
         ELSEIF (cmdbuffer == "-fspread") THEN         ! fraction of global variance used for bandwidth estimation
            ccmd = 10
         ELSEIF (cmdbuffer == "-fpoints") THEN      ! fraction of points used for bandwidth estimation
            ccmd = 11
         ELSEIF (cmdbuffer == "-p") THEN            ! use periodicity
            ccmd = 12
         ELSEIF (cmdbuffer == "-z") THEN            ! add a background to the probability mixture
            ccmd = 13
         ELSEIF (cmdbuffer == "-savegrid") THEN     ! save the indices of grid points
            saveidxs= .TRUE.
         ELSEIF (cmdbuffer == "-readgrid") THEN     ! read the grid points from the standard input
            readgrid= .TRUE.
            ccmd = 14
         ELSEIF (cmdbuffer == "-savevoronois") THEN ! save the Voronoi associations
            savevor= .TRUE.
         ELSEIF (cmdbuffer == "-adj") THEN          ! do cluster merging using adjacency criterion
            saveadj= .TRUE.
            ccmd = 15
         ELSEIF (cmdbuffer == "-merger") THEN       ! cluster with a pk loewr than this are merged with the NN
            ccmd = 16
         ELSEIF (cmdbuffer == "-saveneigh") THEN    ! save gabriel graph
            saveneigh= .TRUE.
         ELSEIF (cmdbuffer == "-readneigh") THEN    ! read gabriel graph
            readneigh= .TRUE.
            ccmd = 17
         ELSEIF (cmdbuffer == "-gs") THEN           ! read gabriel graph neighbor shells
            ccmd = 18
         ELSEIF (cmdbuffer == "-w") THEN            ! use weights
            weighted = .TRUE.
         ELSEIF (cmdbuffer == "-v") THEN            ! verbosity flag
            verbose = .TRUE.
         ELSEIF (cmdbuffer == "-h") THEN            ! help flag
            CALL helpmessage
            CALL EXIT(-1)
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " No parameters specified!"
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 1) THEN                 ! read the cluster smearing
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 2) THEN                 ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN                 ! model file
               fpost=.true.
               clusterfile=trim(cmdbuffer)
            ELSEIF (ccmd == 4) THEN                 ! read the seed for the rng
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 5) THEN                 ! read cutoff for quickshift
               READ(cmdbuffer,*) qs
               IF (qs<0) STOP &
                 "The QS scaling should be positive!"
            ELSEIF (ccmd == 6) THEN                 ! read the number of mean-shift refinement steps
               READ(cmdbuffer,*) nmsopt
            ELSEIF (ccmd == 7) THEN                 ! number of grid points
               READ(cmdbuffer,*) ngrid
            ELSEIF (ccmd == 8) THEN                 ! read the N of bootstrap iterations
               READ(cmdbuffer,*) nbootstrap
               IF (nbootstrap<0) STOP &
                 "The number of iterations should be positive!"
            ELSEIF (ccmd == 9) THEN                 ! read the dimensionality
               READ(cmdbuffer,*) D
               ALLOCATE(period(D))
               period=-1.0d0
            ELSEIF (ccmd == 10) THEN                ! read fractional variance for bandwidth estimation
               READ(cmdbuffer,*) fspread
            ELSEIF (ccmd == 11) THEN                ! read fraction of points for bandwidth estimation
               READ(cmdbuffer,*) fpoints
            ELSEIF (ccmd == 12) THEN                ! read the periodicity in each dimension
               IF (D<0) STOP &
                 "Dimensionality (-d) must precede the periodic lenghts (-p). "
               par_count = 1
               isep1 = 0
               DO WHILE (index(cmdbuffer(isep1+1:), ',') > 0)
                  isep2 = index(cmdbuffer(isep1+1:), ',') + isep1
                  READ(cmdbuffer(isep1+1:isep2-1),*) period(par_count)
                  ! really brute, I know.
                  ! In the case the user will insert 6.28 or 3.14 as periodicity
                  ! the programm will automatically use a better accurancy for pi
                  IF (period(par_count) == 6.28d0) period(par_count) = twopi
                  IF (period(par_count) == 3.14d0) period(par_count) = twopi/2.0d0
                  par_count = par_count + 1
                  isep1=isep2
               ENDDO
               READ(cmdbuffer(isep1+1:),*) period(par_count)
               IF (period(par_count) == 6.28d0) period(par_count) = twopi
               IF (period(par_count) == 3.14d0) period(par_count) = twopi/2.0d0
               periodic=.true.
               IF (par_count/=D) STOP "Check the number of periodic dimensions (-p)!"
            ELSEIF (ccmd == 13) THEN ! read zeta
               READ(cmdbuffer,*) zeta
            ELSEIF (ccmd == 14) THEN                ! read the file containing the grid indexes
               gridfile=trim(cmdbuffer)
            ELSEIF (ccmd == 15) THEN                ! read the threashold for cluster adjancency merging
               READ(cmdbuffer,*) thrmerg
            ELSEIF (ccmd == 16) THEN
               READ(cmdbuffer,*) thrpcl
            ELSEIF (ccmd == 17) THEN                ! read the file containing gabriel graph
               neighfile=trim(cmdbuffer)
            ELSEIF (ccmd == 18) THEN
               READ(cmdbuffer,*) gs                 ! set the neighbor shell for gabriel shift
            ENDIF
         ENDIF
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL SRAND(seed) ! initialize the random number generator

      ! dimensionalty can't be hard coded by default
      IF (D.EQ.-1) THEN
         WRITE(*,*) ""
         WRITE(*,*) " Wrong usage. Insert the dimensionality!"
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! #####################  POST-PROCESSING MODE  ###################
      ! This modality will run just specifying the -gf flag.
      ! The program will just compute the pamm probalities
      ! for each given point
      IF (fpost) THEN
         IF (clusterfile.EQ."NULL") THEN
            ! the user did something wrong in the GM specifications
            WRITE(*,*) &
          " Error: insert the file containing the cluster parameters! "
            CALL helpmessage
            CALL EXIT(-1)
         ENDIF
         OPEN(UNIT=12,FILE=clusterfile,STATUS='OLD',ACTION='READ')
         ! read the model informations from a file.
         IF(periodic)THEN
            ! PERIODIC version
            CALL readvmclusters(12,nk,vmclusters)
            CLOSE(12)
            ALLOCATE(pcluster(nk), px(vmclusters(1)%D))

            DO WHILE (.true.) ! read from the stdin
              READ(*,*,IOSTAT=endf) px
              IF(endf>0) STOP "*** Error occurred while reading file. ***"
              IF(endf<0) EXIT
              ! compute the pamm probability for the point px
              CALL pamm_p_vm(px, pcluster, nk, vmclusters, alpha, zeta)
              !!! decomment if you want to print out
              !!! just the number of the cluster with
              !!! the higher probability
              dummyi1=1
              DO i=1,nk
                 IF (pcluster(i)>pcluster(dummyi1)) dummyi1=i
              ENDDO
              ! write out the number of the
              ! cluster with the highest probability
              WRITE(*,*) px,DLOG(pcluster),dummyi1
            ENDDO
            DEALLOCATE(vmclusters)
         ELSE
            ! NON-PERIODIC version
            CALL readclusters(12,nk,clusters)
            CLOSE(12)
            ALLOCATE(pcluster(nk), px(clusters(1)%D))

            DO WHILE (.true.) ! read from the stdin
              READ(*,*,IOSTAT=endf) px
              IF(endf>0) STOP "*** Error occurred while reading file. ***"
              IF(endf<0) EXIT
              ! compute the pamm probability for the point px
              CALL pamm_p(px, pcluster, nk, clusters, alpha, zeta)
              dummyi1=1
              DO i=1,nk
                 IF (pcluster(i)>pcluster(dummyi1)) dummyi1=i
              ENDDO
              ! write out the number of the
              ! cluster with the highest probability
              WRITE(*,*) px,DLOG(pcluster(dummyi1)),dummyi1
            ENDDO
            DEALLOCATE(clusters)
         ENDIF

         DEALLOCATE(pcluster)
         ! done, go home
         CALL EXIT(-1)
      ENDIF

      ! #####################  CLUSTERING MODE  ###################
      ! get the data from standard input
      CALL readinput(D, weighted, nsamples, x, normwj, wj)

      ! normalize weights -- wj contains the weights of data points
      wj = wj/normwj
      normwj = 1.0d0

      ! If not specified, the number of voronoi polyhedra
      ! are set to the square root of the total number of points
      IF (ngrid.EQ.-1) ngrid = int(sqrt(float(nsamples)))

      ! Initialize the arrays, since now I know the number of
      ! points and the dimensionality
      CALL allocatevectors(D,nsamples,nbootstrap,ngrid,iminij,pnlist,nlist, &
                           y,ni,idmindist,mindist,prob,probboot,idxroot,idxgrid,qspath, &
                           distmm,msmu,tmpmsmu,pabserr,prelerr,normkernel, &
                           wi,lwi,lwj,Q,Qi,logdetHi,Hi,Hiinv,Qiinv,dij, &
                           wlocal,wlocal2,flocal,ineigh,rgrid,sigma2,qscut2, &
                           gabriel,gs,tmps2,Di)

      ! Extract ngrid points on which the kernel density estimation is to be
      ! evaluated. Also partitions the nsamples points into the Voronoi polyhedra
      ! of the sampling points.
      IF(verbose) THEN
         WRITE(*,*) " NSamples: ", nsamples
         WRITE(*,*) " Selecting ", ngrid, " points using MINMAX"
      ENDIF


      ! #############################################################
      ! #                 Selects the grid points                   #
      ! #############################################################
      IF(readgrid)THEN
         ! Reads the grid from file

         IF (gridfile.EQ."NULL") THEN
            WRITE(*,*) &
          " Error: insert the file containing the grid! "
            CALL helpmessage
            CALL EXIT(-1)
         ENDIF

         OPEN(UNIT=12,FILE=gridfile,STATUS='OLD',ACTION='READ')
         ! read the grid from a file
         DO i=1,ngrid
            READ(12,*) idxgrid(i)
         ENDDO
         CLOSE(UNIT=12)

         WRITE(*,*) " Building the Voronoi associations"

         ! does the voronoi associations
         CALL getvoro(D,period,nsamples,ngrid,x,wj,y,ni,iminij,ineigh,wi,idxgrid)
      ELSE
         ! creates the grid with FPS, and at the same time gets the Voronoi association
         CALL mkgrid(D,period,nsamples,ngrid,x,wj,y,ni,iminij,ineigh,wi, &
                  saveidxs,idxgrid,outputfile)
      ENDIF

      ! error check of voronoi association
      ! -- wi contains the weights of the grid points (wi = sum_{j\in Vi} wj) --
      DO i=1,ngrid
        IF (wi(i).EQ.0.0d0) STOP &
          " Error: voronoi has no points associated with - probably two points are perfectly overlapping"
      ENDDO

      ! print out voronoi associations
      IF(savevor) CALL savevoronois(nsamples,iminij,outputfile)

      ! Generate neighbour list between voronoi sets
      IF(verbose) write(*,*) " Generating neighbour list"
        CALL getnlist(nsamples,ngrid,ni,iminij,pnlist,nlist)

      ! precomputes log weights
      lwi = DLOG(wi)
      lwj = DLOG(wj)

      ! #############################################################
      ! #     Generates distance matrix between grid points         #
      ! #############################################################
      IF(verbose) WRITE(*,*) &
        " Precalculate distance matrix between grid points"
      ! distance to closest voronoi
      mindist=HUGE(0.0d0) !  of the kernel density estimator
      DO i=1,ngrid
        IF(verbose .AND. (modulo(i,1000).EQ.0)) &
          WRITE(*,*) i,"/",ngrid
        DO j=1,i-1
          ! distance between two voronoi centers
          distmm(i,j) = pammr2(D,period,y(:,i),y(:,j))
          distmm(j,i) = distmm(i,j)
          IF (distmm(i,j).LT.mindist(i)) THEN
            mindist(i) = distmm(i,j)
            idmindist(i)=j
          ENDIF
          IF (distmm(i,j).LT.mindist(j)) THEN
            mindist(j) = distmm(i,j)
            idmindist(j)=i
          ENDIF
        ENDDO
        ! set distance to myself super far away
        if (gs>0) distmm(i,i) = HUGE(0.0d0)
      ENDDO

      ! #############################################################
      ! #            Computes Gabriel graphs between points         #
      ! #############################################################
      IF (gs > 0) THEN ! PROBABLY THIS SHOULD GO
        IF(readneigh) THEN
          IF(verbose) WRITE(*,*) &
            " Reading gabriel neighbors graph from file"
          ! read gabriel graph
          IF (neighfile.EQ."NULL") THEN
            WRITE(*,*) &
            " Error: insert the file containing the gabriel graph! "
            CALL helpmessage
            CALL EXIT(-1)
          ENDIF
          OPEN(UNIT=12,FILE=neighfile,STATUS='OLD',ACTION='READ')
          ! read graph from a file
          DO i=1,ngrid
             READ(12,*) (gabriel(i,j),j=1,ngrid)
             IF(MODULO(i,100).EQ.0) WRITE(*,*) i,'/',ngrid
          ENDDO
          CLOSE(UNIT=12)
        ELSE
          IF(verbose) WRITE(*,*) &
            " Finding gabriel neighbors between grid points"
          gabriel = .TRUE.
          nn = 0
          !$omp parallel do private(i,j,k) shared(distmm,gabriel,nn)
          DO i=1,ngrid
            gabriel(i,i) = .FALSE.
            DO j=1,ngrid
              IF (.NOT.gabriel(i,j)) CYCLE
              DO k=1,ngrid
                IF (distmm(i,j).GE.(distmm(i,k) + distmm(j,k))) THEN
                  gabriel(i,j) = .FALSE.
                  gabriel(j,i) = .FALSE.
                  EXIT
                ENDIF
              ENDDO
            ENDDO
            nn = nn+1
            IF(MODULO(nn,100).EQ.0) WRITE(*,*) nn,'/',ngrid
          ENDDO
        ENDIF

        IF(saveneigh) THEN
          IF(verbose) WRITE(*,*) &
            " Storing gabriel neighbors graph to file"
          OPEN(UNIT=12,FILE=trim(outputfile)//".neigh",STATUS='REPLACE',ACTION='WRITE')
          DO i=1,ngrid
            WRITE(12,*) gabriel(i,:)
          ENDDO
          CLOSE(UNIT=12)
        ENDIF
      ENDIF ! ENDS GABRIEL STUFF -- PROBABLY THIS SHOULD GO

      ! #############################################################
      ! #            Computes localization weights                  #
      ! #############################################################

      ! Now set localizations. It can be set either in terms of a fraction of the
      ! number samples, or directly as a fraction of the variance of the data
      ! delta to stop bisectioning
      delta = normwj/DBLE(nsamples)

      ! only one of the methods can be used at a time
      IF(fspread.GT.0.0d0) fpoints = -1.0d0

      ! estimate the global covariance from the grid
      CALL covariance(D,period,ngrid,normwj,wi,y,Q)
      WRITE(*,*) "Global eff. dim. ", effdim(D,Q)

      IF(periodic) THEN
        tune = SUM(period**2)
        sigma2 = tune
      ELSE
        tune = 0.0d0
        DO i=1,D
          tune = tune+Q(i,i)
        ENDDO
        sigma2 = tune
      ENDIF

      ! initialize the localization based on fraction of data spread
      IF (fspread.GT.0) sigma2 = sigma2*fspread**2

      IF(verbose) WRITE(*,*) &
        " Estimating kernel density bandwidths"
      DO i=1,ngrid
        IF(verbose .AND. (modulo(i,100).EQ.0)) &
          WRITE(*,*) i,"/",ngrid

        ! get initial flocal using sigma2 equal to tune
        CALL localization(D,period,ngrid,sigma2(i),y,wi,y(:,i),wlocal,flocal(i))

        ! estimate localization based either on fspread or fpoints
        IF(fpoints.GT.0) THEN
        ! ************************************************
        ! *** localization based on fraction of points ***
        ! ************************************************

          ! check if delta is smaller than weigths in voronoi
          lim = fpoints
          IF (fpoints.LE.wi(i)) THEN
            lim = wi(i)+delta
            WRITE(*,*) " Warning: localization smaller than voronoi, increase grid size (meanwhile adjusted localization)!"
          ENDIF

          ! quick approach to ntarget if necessary
          IF (flocal(i).LT.lim) THEN
            DO WHILE(flocal(i).LT.lim)
              ! approach the desired value
              sigma2(i)=sigma2(i)+tune
              CALL localization(D,period,ngrid,sigma2(i),y,wi,y(:,i),wlocal,flocal(i))
            ENDDO
          ENDIF


          ! fine tuning of localization using bisectioning
          j = 1
          DO WHILE(.TRUE.)
            ! fine tuning
            IF(flocal(i).GT.lim) THEN
              sigma2(i) = sigma2(i)-tune/2.0d0**j
            ELSE
              sigma2(i) = sigma2(i)+tune/2.0d0**j
            ENDIF
            CALL localization(D,period,ngrid,sigma2(i),y,wi,y(:,i),wlocal,flocal(i))
            ! exit loop if sigma gives flocal between lim +/- delta
            IF ((flocal(i).LE.lim+delta) .AND. (flocal(i).GE.lim-delta)) EXIT
            ! adjust scaling factor for new sigma
            j = j+1
          ENDDO
        ELSE
        ! ************************************************
        ! *** localization based on fraction of spread ***
        ! ************************************************

          ! consistency check if localization is too small
          IF (sigma2(i).LT.mindist(i)) THEN
            sigma2(i) = mindist(i)
            CALL localization(D,period,ngrid,sigma2(i),y,wi,y(:,i),wlocal,flocal(i))
            WRITE(*,*) " Warning: localization smaller than Voronoi diameter, increase grid size (meanwhile adjusted localization)!"
          ENDIF
        ENDIF

        ! ************************************************
        ! ***  bandwidth estimation from localization  ***
        ! ************************************************

        ! estimate local covariance using grid approximation
        CALL covariance(D,period,ngrid,flocal(i),wlocal,y,Qi)

        ! get local number of points
        nlocal = flocal(i)*DBLE(nsamples)

        ! estimate local dimensionality
        Di(i) = effdim(D,Qi)
        ! oracle shrinkage of covariance matrix
        CALL oracle(D,nlocal,Qi)
        ! inverse local covariance matrix and store it
        CALL invmatrix(D,Qi,Qiinv)

        ! estimate bandwidth from normal (Scott's) reference rule
        Hi = (4.0d0 / ( Di(i)+2.0d0) )**( 2.0d0 / (Di(i)+4.0d0) ) &
           * nlocal**( -2.0d0 / (Di(i)+4.0d0) ) * Qi
        
!        Hi = (4.0d0 / ( DBLE(D)+2.0d0) )**( 2.0d0 / (DBLE(D)+4.0d0) ) &
!           * nlocal**( -2.0d0 / (DBLE(D)+4.0d0) ) * Qi

        ! inverse of the bandwidth matrix
        CALL invmatrix(D,Hi,Hiinv(:,:,i))

        ! estimate logarithmic determinant of local BW H's
        logdetHi(i) = logdet(D,Hi)
        ! estimate the logarithmic normalization constants
        normkernel(i) = DBLE(D)*DLOG(twopi) + logdetHi(i)

        ! adaptive QS cutoff from local covariance
        qscut2(i) = 0.d0
        DO j=1,D
          qscut2(i) = qscut2(i) + Qi(j,j)
        ENDDO
      ENDDO

      ! scale adaptive QS cutoff
      qscut2 = qscut2 * qs**2

      ! #############################################################
      ! #            Computes Kernel Density Estimation             #
      ! #############################################################
      IF(verbose) WRITE(*,*) &
        " Computing kernel density on reference points"
      ! TODO: (1) if we have a mixture of non-periodic and periodic data one could split
      !       this procedure for each dimension...
      !       (2) using gaussians for periodic data
      !           a Gaussian distribution is approximately a van Mises distribution
      !           if van Mises kernel is sufficiently small ...

      ! global cutoff for kde
      kdecut2 = 9.0d0 * (DSQRT(DBLE(D))+1.0d0)**2
      ! setting initial probability to the smallest possible value
      prob = -HUGE(0.0d0)
      DO i=1,ngrid
        IF(verbose .AND. (modulo(i,100).EQ.0)) &
          WRITE(*,*) i,"/",ngrid
        DO j=1,ngrid
          ! renormalize the distance taking into accout the anisotropy of the multidimensional data
          dummd1 = mahalanobis(D,period,y(:,i),y(:,j),Hiinv(:,:,j))
          IF (dummd1.GT.kdecut2) THEN
            ! assume distribution in far away grid point is narrow
            ! and store sum of all contributions in grid point
            ! exponent of the gaussian
            ! natural logarithm of kernel
            lnK = -0.5d0 * (normkernel(j) + dummd1) + lwi(j)
            IF(prob(i).GT.lnK) THEN
              prob(i) = prob(i) + DLOG(1.0d0+DEXP(lnK-prob(i)))
            ELSE
              prob(i) = lnK + DLOG(1.0d0+DEXP(prob(i)-lnK))
            ENDIF
          ELSE
            ! cycle just inside the polyhedra using the neighbour list
            DO k=pnlist(j)+1,pnlist(j+1)
              ! this is the self correction
              IF(nlist(k).EQ.idxgrid(i)) CYCLE
              ! exponent of the gaussian
              dummd1 = mahalanobis(D,period,y(:,i),x(:,nlist(k)),Hiinv(:,:,j))
              ! weighted natural logarithm of kernel
              lnK = -0.5d0 * (normkernel(j) + dummd1) + lwj(nlist(k))
              IF(prob(i).GT.lnK) THEN
                prob(i) = prob(i) + DLOG(1.0d0+DEXP(lnK-prob(i)))
              ELSE
                prob(i) = lnK + DLOG(1.0d0+DEXP(prob(i)-lnK))
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      prob=prob-DLOG(normwj)


      ! #############################################################
      ! #            Computes the Statistical Error on KDE          #
      ! #############################################################
      IF(nbootstrap > 0) THEN
        ! uses bootstrapping to compute the error
        probboot = -HUGE(0.0d0)
        ! open output file for bootstrap containing cluster assignation for each bootstrap run
        OPEN(UNIT=11,FILE=trim(outputfile)//".bs",STATUS='REPLACE',ACTION='WRITE')
        DO nn=1,nbootstrap
          IF(verbose) WRITE(*,*) &
                " Bootstrapping, run ", nn
          ! rather than selecting nsel random points, we select a random
          ! number of points from each voronoi. this makes it possible
          ! to apply some simplifications and avoid computing distances
          ! from far-away voronoi
          nbstot = 0
          DO j=1,ngrid
            ! here we select points and assign them to grid points (i.e. this is an "inside out" version of the KDE code)
            ! we want to loop over grid points and know how many points we should pick from a bootstrapping sample.
            ! this is given by a binomial distribution -- the total number of samples will not be *exactly* nsamples, but will be close enough
            nbssample = random_binomial(nsamples, DBLE(ni(j))/DBLE(nsamples))
            IF (nbssample .eq. 0) CYCLE

            ! calculate "scaled" weight for contribution from far away voronoi
            ! we take into account the fact that we might have selected a number of samples different from ni(j)
            dummd2 = DLOG(DBLE(nbssample)/ni(j)) * lwi(j)

            nbstot = nbstot+nbssample
            DO i=1,ngrid
              ! this is the distance between the grid point from which we're sampling (j) and the one on which we're accumulating the KDE (i)
              dummd1 = mahalanobis(D,period,y(:,i),y(:,j),Hiinv(:,:,j))
              IF (dummd1.GT.kdecut2) THEN
                ! if the two cells are far apart, we just compute an "average contribution" from the far away Voronoi
                lnK = -0.5d0 * (normkernel(j) + dummd1) + dummd2
                IF(probboot(i,nn).GT.lnK) THEN
                  probboot(i,nn) = probboot(i,nn) + DLOG(1.0d0+DEXP(lnK-probboot(i,nn)))
                ELSE
                  probboot(i,nn) = lnK + DLOG(1.0d0+DEXP(probboot(i,nn)-lnK))
                ENDIF
              ELSE
                ! do the actual bootstrapping selection for this Voronoi
                DO k=1,nbssample
                  rndidx = int(ni(j)*random_uniform())+1
                  rndidx = nlist(pnlist(j)+rndidx)
                  IF ( rndidx.EQ.idxgrid(i) ) CYCLE
                  dummd1 = mahalanobis(D,period,y(:,i),x(:,rndidx),Hiinv(:,:,j))
                  lnK = -0.5d0 * (normkernel(j) + dummd1) + lwj(rndidx)
                  IF(probboot(i,nn).GT.lnK) THEN
                    probboot(i,nn) = probboot(i,nn) + DLOG(1.0d0+DEXP(lnK-probboot(i,nn)))
                  ELSE
                    probboot(i,nn) = lnK + DLOG(1.0d0+DEXP(probboot(i,nn)-lnK))
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          ! normalizes the probability estimate, keeping into account that we might have used a different number of sample points than nsamples
          probboot(:,nn) = probboot(:,nn)-(DLOG(normwj)+DLOG(DBLE(nbstot)/DBLE(nsamples)))

          !@@@@@@@@@@@@@@@

          ! #############################################################
          ! #                       Runs quick-shift                    #
          ! #############################################################
          IF(verbose) WRITE(*,*) " Starting Quick-Shift"
          idxroot=0
          DO i=1,ngrid
             IF(idxroot(i).NE.0) CYCLE
             IF(verbose .AND. (modulo(i,1000).EQ.0)) &
                   WRITE(*,*) i,"/",ngrid
             qspath=0
             qspath(1)=i
             counter=1
             DO WHILE(qspath(counter).NE.idxroot(qspath(counter)))
                IF (gs > 0) THEN
                  idxroot(qspath(counter)) = gs_next(ngrid,qspath(counter),probboot(:,nn),distmm,gabriel,gs)
                ELSE
                  idxroot(qspath(counter)) = qs_next(ngrid,qspath(counter),idmindist(qspath(counter)), &
                                             probboot(:,nn),distmm,qscut2(qspath(counter)))
                ENDIF
                IF(idxroot(idxroot(qspath(counter))).NE.0) EXIT
                counter=counter+1
                qspath(counter)=idxroot(qspath(counter-1))
             ENDDO
             DO j=1,counter
                ! we found a new root, and we now set this point as the root
                ! for all the point that are in this qspath
                idxroot(qspath(j))=idxroot(idxroot(qspath(counter)))
             ENDDO
          ENDDO

          ! get the cluster centers
          CALL unique(ngrid,idxroot,clustercenters)
          ! get the number of the clusters
          Nk=SIZE(clustercenters)

!          WRITE(comment,'(I0.5)') nn
!          OPEN(UNIT=11,FILE=trim(outputfile)//"-bs"//trim(comment)//".grid",STATUS='REPLACE',ACTION='WRITE')
          DO i=1,ngrid
!             DO j=1,D
!               WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", y(j,i)
!             ENDDO

             !print out grid file with additional information on probability, errors, localization, weights in voronoi, dim
!             WRITE(11,"(A1,I5,A1,ES18.7E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4)") &
!                 " " , MINLOC(ABS(clustercenters-idxroot(i)),1) ,      &
!                                                  " " , prob(i)

             ! for the merging we just need the cluster assignation of each bootstrap run
             WRITE(11,"(A1,I5)",ADVANCE="NO") " ",MINLOC(ABS(clustercenters-idxroot(i)),1)
          ENDDO
          WRITE(11,*) " "
!          CLOSE(UNIT=11)
        ENDDO ! ends loop on bootstrapping iterations
        ! computes the bootstrap error from the statistics of the nbootstrap KDE runs
        pabserr=0.0d0
        DO i=1,ngrid
          DO j=1,nbootstrap
            IF (probboot(i,j).GT.prob(i)) THEN
              pabserr(i) = pabserr(i) + DEXP(2.0d0*(probboot(i,j)+DLOG(1.0d0-DEXP(prob(i)-probboot(i,j)))))
            ELSE
              pabserr(i) = pabserr(i) + DEXP(2.0d0*(prob(i)+DLOG(1.0d0-DEXP(probboot(i,j)-prob(i)))))
            ENDIF
          ENDDO
          pabserr(i) = DLOG(DSQRT(pabserr(i) / (nbootstrap-1.0d0)))
          prelerr(i) = pabserr(i) - prob(i)
        ENDDO
        CLOSE(UNIT=11)
      ELSE
        ! uses a binomial-distribution ansatz to estimate the error
        DO i=1,ngrid
          prelerr(i) =DLOG(DSQRT(((((mindist(i)*twopi)**(-Di(i)))/DEXP(prob(i)))-1.0d0)/nsamples))


          !prelerr(i)= DSQRT(( ( (sigma2(i)**(-Di(i))) * &
          !                      (twopi**(-Di(i)/2.0d0))/ &
          !                       DEXP(prob(i)) )-1.0d0)/normwj)
          !IF(prelerr(i) .NE. prelerr(i))THEN
            ! if we get a NaN then we compute directly the log of the relativer error
            ! using a first order exapansion
          !  prelerr(i)=DLOG((sigma2(i)**(-Di(i)))*(twopi**(-Di(i)/2.0d0))/normwj) &
          !          -0.5d0*prob(i)
          !ELSE
          !  prelerr(i)=DLOG(prelerr(i))
          !ENDIF
          pabserr(i)=prelerr(i)+prob(i)
        ENDDO
      ENDIF

      ! #############################################################
      ! #                       Runs quick-shift                    #
      ! #############################################################
      IF(verbose) WRITE(*,*) " Starting Quick-Shift"
      idxroot=0
      DO i=1,ngrid
         IF(idxroot(i).NE.0) CYCLE
         IF(verbose .AND. (modulo(i,1000).EQ.0)) &
               WRITE(*,*) i,"/",ngrid
         qspath=0
         qspath(1)=i
         counter=1
         DO WHILE(qspath(counter).NE.idxroot(qspath(counter)))
            IF (gs > 0) THEN
              idxroot(qspath(counter)) = gs_next(ngrid,qspath(counter),prob,distmm,gabriel,gs)
            ELSE
              idxroot(qspath(counter)) = qs_next(ngrid,qspath(counter),idmindist(qspath(counter)), &
                                         prob,distmm,qscut2(qspath(counter)))
            ENDIF
            IF(idxroot(idxroot(qspath(counter))).NE.0) EXIT
            counter=counter+1
            qspath(counter)=idxroot(qspath(counter-1))
         ENDDO
         DO j=1,counter
            ! we found a new root, and we now set this point as the root
            ! for all the point that are in this qspath
            idxroot(qspath(j))=idxroot(idxroot(qspath(counter)))
         ENDDO
      ENDDO

      ! get the cluster centers
      CALL unique(ngrid,idxroot,clustercenters)
      ! get the number of the clusters
      Nk=SIZE(clustercenters)
      ! get sum of the probs, the normalization factor
      normpks=logsumexp(ngrid,idxroot/idxroot,prob,1)
      ! check if there are outliers that should be merged to the others
      ALLOCATE(mergeornot(Nk))
      DO k=1,Nk
         ! compute the relative weight of the cluster
         dummd1=DEXP(logsumexp(ngrid,idxroot,prob,clustercenters(k))-normpks)
         mergeornot(k)=dummd1.LT.thrpcl
      ENDDO
      ! merge the outliers
      DO i=1,Nk
        IF(.NOT.mergeornot(i)) CYCLE
        dummyi1=clustercenters(i)
        dummd1=HUGE(0.0d0)
        DO j=1,Nk
          IF(mergeornot(j)) CYCLE
          ! the mahalanobis does strange things. Let's use the L2norm.
          dummd2=pammr2(D,period,y(:,idxroot(dummyi1)),&
                          y(:,idxroot(clustercenters(j))))
          IF(dummd2.LT.dummd1) THEN
             dummd1=dummd2
             clustercenters(i)=clustercenters(j)
          ENDIF
        ENDDO
        WHERE(idxroot.EQ.dummyi1) idxroot=clustercenters(i)
      ENDDO
      IF(COUNT(mergeornot).GT.0)THEN
         DEALLOCATE(clustercenters)
         ! get the new cluster centers
         CALL unique(ngrid,idxroot,clustercenters)
         IF(verbose) WRITE(*,*) Nk-size(clustercenters), &
                " clusters where merged into other clusters"
         Nk=size(clustercenters)
         ! get the real maxima in the cluster, considering the errorbar
         DO i=1,Nk
            dummyi1=clustercenters(i)
            clustercenters(i) = getidmax(ngrid,idxroot,prob,pabserr,clustercenters(i))
            ! reassign the proper cluster root to each cluster points
            WHERE(idxroot.EQ.dummyi1) idxroot=clustercenters(i)
         ENDDO
      ENDIF

      IF(verbose) write(*,*) "Writing out"
      OPEN(UNIT=11,FILE=trim(outputfile)//".grid",STATUS='REPLACE',ACTION='WRITE')
      DO i=1,ngrid
         DO j=1,D
           WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", y(j,i)
         ENDDO

         CALL invmatrix(D,Hiinv(:,:,i),Hi)

         !print out grid file with additional information on probability, errors, localization, weights in voronoi, dim
         WRITE(11,"(A1,I5,A1,ES18.7E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4)") &
             " " , MINLOC(ABS(clustercenters-idxroot(i)),1) ,      &
                                              " " , prob(i) ,      &
                                              " " , pabserr(i),    &
                                              " " , prelerr(i),    &
                                              " " , sigma2(i),     &
                                              " " , flocal(i),     &
                                              " " , wi(i),         &
                                              " " , Di(i),         &
                                              " " , trmatrix(D,Hi)/DBLE(D)

      ENDDO

      CLOSE(UNIT=11)

      ! now we can procede and complete the definition of probability model
      ! now qspath contains the indexes of Nk gaussians
      IF(periodic) THEN
         ALLOCATE(vmclusters(Nk))
      ELSE
         ALLOCATE(clusters(Nk))
      ENDIF

      DO k=1,Nk
         IF(verbose) WRITE(*,*) k, "/", Nk
         IF(periodic)THEN
            ALLOCATE(vmclusters(k)%mean(D))
            ALLOCATE(vmclusters(k)%cov(D,D))
            ALLOCATE(vmclusters(k)%icov(D,D))
            ALLOCATE(vmclusters(k)%period(D))
            vmclusters(k)%period=period
            vmclusters(k)%mean=y(:,clustercenters(k))
         ELSE
            ALLOCATE(clusters(k)%mean(D))
            ALLOCATE(clusters(k)%cov(D,D))
            ALLOCATE(clusters(k)%icov(D,D))
            clusters(k)%mean=y(:,clustercenters(k))
         ENDIF

         ! optionally do a few mean-shift steps to find a better estimate
         ! of the cluster mode

         DO j=1,nmsopt
            msmu=0.0d0
            tmppks=-HUGE(0.0d0)

            DO i=1,ngrid
               dummd1 = mahalanobis(D,period,y(:,i),y(:,clustercenters(k)),&
                                    Hiinv(:,:,clustercenters(k)))
               msw = -0.5d0 * (normkernel(clustercenters(k)) + dummd1) + prob(i)
               CALL pammrij(D,period,y(:,i),y(:,clustercenters(k)),tmpmsmu)

               msmu = msmu + DEXP(msw)*tmpmsmu

               ! log-sum-exp
               IF(tmppks.GT.msw) THEN
                tmppks = tmppks + DLOG(1.0d0+DEXP(msw-tmppks))
               ELSE
                tmppks = msw + DLOG(1.0d0+DEXP(tmppks-msw))
               ENDIF
            ENDDO

            IF(periodic)THEN
               vmclusters(k)%mean = vmclusters(k)%mean + msmu / DEXP(tmppks)
            ELSE
               clusters(k)%mean = clusters(k)%mean + msmu / DEXP(tmppks)
            ENDIF
         ENDDO

         ! compute the gaussians covariance from the data in the clusters

         IF(periodic)THEN
            CALL getlcovclusterp(D,period,ngrid,nsamples,prob,y,idxroot,clustercenters(k),vmclusters(k)%cov)

  !write(*,*) vmclusters(k)%cov
            ! If we have a cluster with one point we compute the weighted covariance with
            ! the points in the Voronoi
            IF(COUNT(idxroot.EQ.clustercenters(k)).EQ.1) THEN
              !CALL getcovclusterp(D,period,nsamples,wj,x,iminij,clustercenters(k),vmclusters(k)%cov)
              WRITE(*,*) " Warning: single point cluster!!! "
            ENDIF
            vmclusters(k)%weight= &
             DEXP(logsumexp(ngrid,idxroot,prob,clustercenters(k))-normpks)
            vmclusters(k)%D=D
            !CALL oracle(D,logsumexp(ngrid,idxroot,prob,clustercenters(k))*DBLE(nsamples),vmclusters(k)%cov)
         ELSE
            CALL getlcovcluster(D,period,ngrid,prob,y,idxroot,clustercenters(k),clusters(k)%cov)
  !write(*,*) clusters(k)%cov
            ! If we have a cluster with one point we compute the weighted covariance with
            ! the points in the Voronoi
            IF(COUNT(idxroot.EQ.clustercenters(k)).EQ.1) THEN
              !CALL getcovcluster(D,period,nsamples,wj,x,iminij,clustercenters(k),clusters(k)%cov)
              WRITE(*,*) " Warning: single point cluster!!! "
            ENDIF
            clusters(k)%weight=&
             DEXP(logsumexp(ngrid,idxroot,prob,clustercenters(k))-normpks)
            clusters(k)%D=D
            CALL oracle(D,logsumexp(ngrid,idxroot,prob,clustercenters(k))*DBLE(nsamples),clusters(k)%cov)
         ENDIF
      ENDDO

      IF(periodic)THEN
         ! write the VM distributions
         ! write a 2-lines header containig a bit of information
         WRITE(comment,*) "# PAMMv2 clusters analysis. NSamples: ", nsamples, " NGrid: ", &
                   ngrid, " QSLambda: ", qs, ACHAR(10), &
                   "# Dimensionality/NClusters//Pk/Mean/Covariance/Period"

         OPEN(UNIT=12,FILE=trim(outputfile)//".pamm",STATUS='REPLACE',ACTION='WRITE')
         CALL writevmclusters(12, comment, nk, vmclusters)
         CLOSE(UNIT=12)
         DEALLOCATE(vmclusters)
      ELSE
         ! write the Gaussians
         ! write a 2-lines header
         WRITE(comment,*) "# PAMMv2 clusters analysis. NSamples: ", nsamples, " NGrid: ", &
                   ngrid, " QSLambda: ", qs, ACHAR(10), "# Dimensionality/NClusters//Pk/Mean/Covariance"

         OPEN(UNIT=12,FILE=trim(outputfile)//".pamm",STATUS='REPLACE',ACTION='WRITE')

         CALL writeclusters(12, comment, nk, clusters)
         CLOSE(UNIT=12)
         ! maybe I should deallocate better..
         DEALLOCATE(clusters)
      ENDIF

      DEALLOCATE(x,wj,Di)
      DEALLOCATE(period)
      DEALLOCATE(idxroot,qspath,distmm,idxgrid)
      DEALLOCATE(pnlist,nlist,iminij)
      DEALLOCATE(y,ni,prob,sigma2,rgrid,wi)
      DEALLOCATE(msmu,tmpmsmu)
      DEALLOCATE(Q,Qi,Hi,Hiinv,Qiinv,normkernel)
      DEALLOCATE(dij,tmps2)
      DEALLOCATE(wlocal,flocal,ineigh,mindist,idmindist)
      DEALLOCATE(prelerr,pabserr,clustercenters,mergeornot)
      IF(saveadj) DEALLOCATE(macrocl,sortmacrocl)
      IF(nbootstrap>0) DEALLOCATE(probboot)

      CALL EXIT(0)
      ! end of the main programs



!!!!! FUCTIONS and SUBROUTINES !!!!!!!!!!!!!!!!!!!!

      CONTAINS

      SUBROUTINE helpmessage
         ! Banner to print out for helping purpose
         !

         WRITE(*,*) ""
         WRITE(*,*) " USAGE: pamm [-h] -d D [-p 6.28,6.28,...] [-w] [-o output] [-ngrid ngrid] "
         WRITE(*,*) "             [-qs scale] [-kde err] [-z zeta_factor] [-a smoothing_factor] "
         WRITE(*,*) "             [-seed seedrandom] [-rif -1,0,0,...] [-v] "
         WRITE(*,*) ""
         WRITE(*,*) " Applies the PAMM clustering to a high-dimensional data set. "
         WRITE(*,*) " It is mandatory to specify the dimensionality of the data, which "
         WRITE(*,*) " must be passed through the standard input in the format: "
         WRITE(*,*) " x11 x12 x13 ... x1D [w1] "
         WRITE(*,*) " x21 x22 x23 ... x2D [w2] "
         WRITE(*,*) ""
         WRITE(*,*) " For other options a default is defined.  "
         WRITE(*,*) ""
         WRITE(*,*) "   -h                : Print this message "
         WRITE(*,*) "   -d  D             : Dimensionality "
         WRITE(*,*) "   -w                : Reads weights for the sample points [default: no weight] "
         WRITE(*,*) "   -o  output        : Prefix for output files [out]. This will produce : "
         WRITE(*,*) "                         output.grid (clusterized grid points) "
         WRITE(*,*) "                         output.pamm (cluster parameters) "
         WRITE(*,*) "   -savevoronois     : Save Voronoi associations. This will produce:"
         WRITE(*,*) "                         output.voronoislinks (points + associated Voronoi) "
         WRITE(*,*) "   -savegrid         : Save the indexes of the point selected as grid points "
         WRITE(*,*) "                       This will produce the file output.idxs "
         WRITE(*,*) "   -readgrid g.idxs  : Read from file the indexes of the points to be selected "
         WRITE(*,*) "                       as grid points "
         WRITE(*,*) "   -qs scaling       : Scaling factor used during the QS clustering "
         WRITE(*,*) "   -ngrid ngrid      : Number of grid points to evaluate KDE [sqrt(nsamples)]"
         WRITE(*,*) "   -bootstrap N      : Number of iteretions to do when using bootstrapping "
         WRITE(*,*) "                       to refine the KDE on the grid points"
         WRITE(*,*) "   -nms nms          : Do nms mean-shift steps with a Gaussian width qscut/5 to "
         WRITE(*,*) "                       optimize cluster centers [0] "
         WRITE(*,*) "   -kderr target     : Target relative error in the KDE [0.1] "
         WRITE(*,*) "   -seed seed        : Seed to initialize the random number generator. [12345]"
         WRITE(*,*) "   -p P1,...,PD      : Periodicity in each dimension [ (6.28,6.28,6.28,...) ]"
         WRITE(*,*) "   -merger threshold : Probability threshold below which a cluster should be merged "
         WRITE(*,*) "                       to the closest neighbour [0] "
         WRITE(*,*) "   -fpoints fraction : Fraction of the total number of samples to be included in "
         WRITE(*,*) "                       each local baloon during the localization step. This is activated"
         WRITE(*,*) "                       by default and can't be used together with -fspread [0.14] "
         WRITE(*,*) "   -fspread fraction : Fraction of the covariance of the total dataset to be used a the size for "
         WRITE(*,*) "                       each local baloon during the localization step. This is deactivated"
         WRITE(*,*) "                       by default and can't be used together with -fpoints [0.14] "
         WRITE(*,*) "   -v                : Verbose output "
         WRITE(*,*) ""
         WRITE(*,*) " Post-processing mode (-gf): this reads high-dim data and computes the "
         WRITE(*,*) " cluster probabilities associated with them, given the output of a "
         WRITE(*,*) " previous PAMM analysis. "
         WRITE(*,*) ""
         WRITE(*,*) "   -gf               : File to read reference Gaussian clusters from"
         WRITE(*,*) "   -a                : Additional smearing of clusters "
         WRITE(*,*) "   -z zeta_factor    : Probabilities below this threshold are counted as 'no cluster' [default:0]"
         WRITE(*,*) ""

      END SUBROUTINE helpmessage

      DOUBLE PRECISION FUNCTION logsumexp(ngrid,v1,prob,clusterid)
         INTEGER, INTENT(IN) :: ngrid,v1(ngrid),clusterid
         DOUBLE PRECISION, INTENT(IN) :: prob(ngrid)

         INTEGER ii
         DOUBLE PRECISION :: tmpv(ngrid)

         ! select just the probabilities of the element
         ! that belongs to the cluster tgt
         tmpv = prob
         WHERE (v1 .EQ. clusterid)
            tmpv = prob
         ELSEWHERE
            tmpv = -HUGE(0.0d0)
         END WHERE
         ! do log-sum-exp trick
         logsumexp = -HUGE(0.0d0)
         DO ii=1,ngrid
            IF(prob(ii).EQ.-HUGE(0.0d0)) CYCLE
            IF(logsumexp.GT.tmpv(ii)) THEN
               logsumexp = logsumexp + DLOG(1.0d0+DEXP(tmpv(ii)-logsumexp))
            ELSE
               logsumexp = tmpv(ii) + DLOG(1.0d0+DEXP(logsumexp-tmpv(ii)))
            ENDIF
         ENDDO
      END FUNCTION logsumexp

      INTEGER FUNCTION getidmax(ngrid,v1,prob,abser,clusterid)
         INTEGER, INTENT(IN) :: ngrid,v1(ngrid),clusterid
         DOUBLE PRECISION, INTENT(IN) :: prob(ngrid),abser(ngrid)

         DOUBLE PRECISION :: tmpv(ngrid)

         ! select just the probabilities of the element
         ! that belongs to the cluster tgt
         tmpv = prob
         WHERE (v1 .EQ. clusterid)
            tmpv = DEXP(prob)+DEXP(abser)
         ELSEWHERE
            tmpv = -HUGE(0.0d0)
         END WHERE
         getidmax = MAXLOC(tmpv,1)
      END FUNCTION getidmax

      DOUBLE PRECISION FUNCTION median(ngrid,a)
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, INTENT(IN) :: a(ngrid)

         INTEGER :: l
         DOUBLE PRECISION, DIMENSION(SIZE(a,1)) :: ac

         IF ( SIZE(a,1) < 1 ) THEN
         ELSE
           ac = a
           ! this is not an intrinsic: peek a sort algo from
           ! Category:Sorting, fixing it to work with real if
           ! it uses integer instead.
           CALL sort(ac,ngrid)

           l = SIZE(a,1)
           IF ( mod(l, 2) == 0 ) THEN
               median = (ac(l/2+1) + ac(l/2))/2.0
           ELSE
               median = ac(l/2+1)
           END IF
         END IF
      END FUNCTION median

      SUBROUTINE allocatevectors(D,nsamples,nbootstrap,ngrid,iminij,pnlist,nlist, &
                                 y,ni,idmindist,mindist,prob,probboot,idxroot,idxgrid,qspath, &
                                 distmm,msmu,tmpmsmu,pabserr,prelerr,normkernel, &
                                 wi,lwi,lwj,Q,Qi,logdetHi,Hi,Hiinv,Qiinv,dij, &
                                 wlocal,wlocal2,flocal,ineigh,rgrid,sigma2,qscut2, &
                                 gabriel,gs,tmps2,Di)

         INTEGER, INTENT(IN) :: D,nsamples,nbootstrap,ngrid,gs
         INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT):: iminij,pnlist,nlist,idxroot,idxgrid,qspath
         INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: ni,ineigh,idmindist
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: prob,msmu,tmpmsmu,wi,lwi,lwj,logdetHi
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: pabserr,prelerr,normkernel,wlocal,wlocal2,flocal
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: dij,sigma2,qscut2,rgrid,tmps2,Di,mindist
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: Q,Qi,Qiinv,Hi,probboot
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: y,distmm
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT) :: Hiinv
         LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: gabriel


         IF (ALLOCATED(iminij))     DEALLOCATE(iminij)
         IF (ALLOCATED(pnlist))     DEALLOCATE(pnlist)
         IF (ALLOCATED(nlist))      DEALLOCATE(nlist)
         IF (ALLOCATED(y))          DEALLOCATE(y)
         IF (ALLOCATED(ni))         DEALLOCATE(ni)
         IF (ALLOCATED(mindist))    DEALLOCATE(mindist)
         IF (ALLOCATED(idmindist))    DEALLOCATE(idmindist)
         IF (ALLOCATED(prob))       DEALLOCATE(prob)
         IF (ALLOCATED(probboot))   DEALLOCATE(probboot)
         IF (ALLOCATED(idxroot))    DEALLOCATE(idxroot)
         IF (ALLOCATED(idxgrid))    DEALLOCATE(idxgrid)
         IF (ALLOCATED(distmm))     DEALLOCATE(distmm)
         IF (ALLOCATED(msmu))       DEALLOCATE(msmu)
         IF (ALLOCATED(tmpmsmu))    DEALLOCATE(tmpmsmu)
         IF (ALLOCATED(pabserr))    DEALLOCATE(pabserr)
         IF (ALLOCATED(prelerr))    DEALLOCATE(prelerr)
         IF (ALLOCATED(wi))         DEALLOCATE(wi)
         IF (ALLOCATED(lwi))        DEALLOCATE(lwi)
         IF (ALLOCATED(lwj))        DEALLOCATE(lwj)
         IF (ALLOCATED(normkernel)) DEALLOCATE(normkernel)
         IF (ALLOCATED(Q))          DEALLOCATE(Q)
         IF (ALLOCATED(Qi))         DEALLOCATE(Qi)
         IF (ALLOCATED(logdetHi))   DEALLOCATE(logdetHi)
         IF (ALLOCATED(Hi))         DEALLOCATE(Hi)
         IF (ALLOCATED(Hiinv))      DEALLOCATE(Hiinv)
         IF (ALLOCATED(Qiinv))      DEALLOCATE(Qiinv)
         IF (ALLOCATED(dij))        DEALLOCATE(dij)
         IF (ALLOCATED(wlocal))     DEALLOCATE(wlocal)
         IF (ALLOCATED(wlocal2))    DEALLOCATE(wlocal2)
         IF (ALLOCATED(flocal))     DEALLOCATE(flocal)
         IF (ALLOCATED(sigma2))     DEALLOCATE(sigma2)
         IF (ALLOCATED(qscut2))     DEALLOCATE(qscut2)
         IF (ALLOCATED(tmps2))      DEALLOCATE(tmps2)
         IF (ALLOCATED(rgrid))      DEALLOCATE(rgrid)
         IF (ALLOCATED(Di))         DEALLOCATE(Di)
         IF (ALLOCATED(gabriel))    DEALLOCATE(gabriel)



         ! Initialize the arrays, since now I know the number of
         ! points and the dimensionality
         ALLOCATE(iminij(nsamples))
         ALLOCATE(pnlist(ngrid+1), nlist(nsamples), lwj(nsamples))
         ALLOCATE(y(D,ngrid), ni(ngrid), mindist(ngrid), idmindist(ngrid))
         ALLOCATE(prob(ngrid), sigma2(ngrid), qscut2(ngrid), rgrid(ngrid))
         ALLOCATE(idxroot(ngrid), qspath(ngrid), distmm(ngrid,ngrid))
         ALLOCATE(msmu(D), tmpmsmu(D),logdetHi(ngrid))
         ALLOCATE(pabserr(ngrid),prelerr(ngrid),normkernel(ngrid),wi(ngrid),lwi(ngrid))
         ! bootstrap probability density array will be allocated if necessary
         IF(nbootstrap > 0) ALLOCATE(probboot(ngrid,nbootstrap))
         ! Allocate variables for local bandwidth estimate
         ALLOCATE(Q(D,D),Qi(D,D))
         ALLOCATE(Hi(D,D),Hiinv(D,D,ngrid),Qiinv(D,D))
         ALLOCATE(dij(D),Di(ngrid))
         ALLOCATE(idxgrid(ngrid),tmps2(ngrid))
         ALLOCATE(wlocal(ngrid))
         ALLOCATE(wlocal2(nsamples))
         ALLOCATE(flocal(ngrid))
         ALLOCATE(ineigh(ngrid))
         ! allocate gabriel graph if desired
         IF (gs > 0) ALLOCATE(gabriel(ngrid,ngrid))
      END SUBROUTINE allocatevectors

      SUBROUTINE localization(D,period,N,s2,x,w,y,wl,num)
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: s2
         DOUBLE PRECISION, INTENT(IN) :: x(D,N)
         DOUBLE PRECISION, INTENT(IN) :: y(D)
         DOUBLE PRECISION, INTENT(IN) :: w(N)
         DOUBLE PRECISION, INTENT(OUT) :: wl(N)
         DOUBLE PRECISION, INTENT(OUT) :: num

         INTEGER ii
         DOUBLE PRECISION xy(D,N)

         DO ii=1,D
           xy(ii,:) = x(ii,:)-y(ii)
           IF (period(ii) > 0.0d0) THEN
             ! scaled lenght
             xy(ii,:) = xy(ii,:)/period(ii)
             ! Finds the smallest separation between the images of the vector elements
             xy(ii,:) = xy(ii,:) - DNINT(xy(ii,:)) ! Minimum Image Convention
             ! Rescale back the length
             xy(ii,:) = xy(ii,:)*period(ii)
           ENDIF
         ENDDO
         ! estimate weights for localization as product from
         ! spherical gaussian weights and weights in voronoi
         wl = DEXP(-0.5d0/s2*SUM(xy*xy,1))*w
         ! estimate local number of sample points
         num = SUM(wl)
      END SUBROUTINE localization

      INTEGER FUNCTION  findMinimum(x, startidx, endidx )
         INTEGER, INTENT(IN) :: startidx, endidx
         DOUBLE PRECISION, DIMENSION(1:), INTENT(IN) :: x

         DOUBLE PRECISION :: minimum
         INTEGER :: location
         INTEGER :: i

         minimum  = x(startidx)   ! assume the first is the min
         location = startidx      ! record its position
         DO i = startidx+1, endidx    ! start with next elements
            IF (x(i) < minimum) THEN  !   if x(i) less than the min?
               minimum  = x(i)        !      Yes, a new minimum found
               location = i           !      record its position
            END IF
         END DO
         findMinimum = Location       ! return the position
      END FUNCTION  findMinimum

      SUBROUTINE swap(aa, bb)
      !  This subroutine swaps the values of its two formal arguments.
         DOUBLE PRECISION, INTENT(INOUT) :: aa, bb
         DOUBLE PRECISION                :: temp
         temp = aa
         aa = bb
         bb = temp
      END SUBROUTINE swap

      SUBROUTINE swapi(aa, bb)
      !  This subroutine swaps the values of its two formal arguments.
         INTEGER, INTENT(INOUT) :: aa, bb
         INTEGER                :: temp
         temp = aa
         aa = bb
         bb = temp
      END SUBROUTINE swapi

      SUBROUTINE unique(ns,vin,vout)
        INTEGER, INTENT(IN) :: ns
        INTEGER, DIMENSION(:), INTENT(IN) :: vin
        INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: vout
        INTEGER :: kk,ii,jj,zz(ns)
        kk = 1
        zz(1) = vin(1)
        outer: do ii=2,ns
          do jj=1,kk
            if (zz(jj) == vin(ii)) then
              ! Found a match so start looking again
              cycle outer
            end if
          end do
          ! No match found so add it to the output
          kk = kk + 1
          zz(kk) = vin(ii)
        end do outer

        IF (ALLOCATED(vout)) DEALLOCATE(vout)
        ALLOCATE(vout(kk))
        vout(1:kk)=zz(1:kk)
      END SUBROUTINE unique

      SUBROUTINE sort(x, nn)
         ! This subroutine receives an array x() and sorts it into ascending order.
         INTEGER, INTENT(IN) :: nn
         DOUBLE PRECISION, INTENT(INOUT) :: x(nn)

         INTEGER  :: i
         INTEGER  :: location

         DO i = 1, nn-1 ! except for the last
            location = findMinimum(x, i, nn) ! find min from this to last
            CALL  swap(x(i), x(location)) ! swap this and the minimum
         END DO
      END SUBROUTINE sort

      SUBROUTINE argsort(x, sidx, nn)
         ! This subroutine receives an array x() and sorts it into ascending order.
         INTEGER, INTENT(IN) :: nn
         DOUBLE PRECISION, INTENT(IN) :: x(nn)
         INTEGER, INTENT(OUT) :: sidx(nn)

         INTEGER  i,location

         DO i = 1, nn
            sidx(i) = i
         ENDDO

         DO i = 1, nn-1 ! except for the last
            location = findMinimum(x, i, nn) ! find min from this to last
            CALL  swapi(sidx(i), sidx(location)) ! swap this and the minimum
         END DO
      END SUBROUTINE argsort

      ! computes a weighted covariance
      SUBROUTINE covariance(D,period,N,wnorm,w,x,Q)
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: wnorm
         DOUBLE PRECISION, INTENT(IN) :: w(N)
         DOUBLE PRECISION, INTENT(IN) :: x(D,N)
         DOUBLE PRECISION, INTENT(OUT) :: Q(D,D)

         DOUBLE PRECISION xm(D)         ! mean of each dimension
         DOUBLE PRECISION xxm(D,N)      ! difference of x and xm
         DOUBLE PRECISION xxmw(D,N)     ! weighted difference of x and xm

         DOUBLE PRECISION sumsin,sumcos     ! sum of cos and sin

!         DOUBLE PRECISION sumcos,sumsin

         INTEGER ii

         DO ii=1,D
           ! find the mean for periodic or non periodic data
           IF (period(ii) > 0.0d0) THEN
             sumsin = SUM(w*SIN(x(ii,:)*twopi/period(ii)))/wnorm
             sumcos = SUM(w*COS(x(ii,:)*twopi/period(ii)))/wnorm
             xm(ii) = ATAN2(sumsin,sumcos)
           ELSE
             xm(ii) = SUM(x(ii,:)*w)/wnorm
           ENDIF

           xxm(ii,:) = x(ii,:) - xm(ii)
           IF (period(ii) > 0.0d0) THEN
             ! this is the correct way
             xxm(ii,:) = xxm(ii,:) - DNINT(xxm(ii,:)/period(ii)) * period(ii)

!             ! scaled length
!              = xxm(ii,:)/period(ii)
!             ! Finds the smallest separation between the images of the vector elements
!             xxm(ii,:) = xxm(ii,:) - DNINT(xxm(ii,:)) ! Minimum Image Convention
!             ! Rescale back the length
!             xxm(ii,:) = xxm(ii,:)*period(ii)
           ENDIF
           xxmw(ii,:) = xxm(ii,:) * w/wnorm
         ENDDO
         CALL DGEMM("N", "T", D, D, N, 1.0d0, xxm, D, xxmw, D, 0.0d0, Q, D)
         Q = Q / (1.0d0-SUM((w/wnorm)**2.0d0))
      END SUBROUTINE covariance

      SUBROUTINE getlcovcluster(D,period,N,prob,x,clroots,idcl,Q)
      ! log version
         INTEGER, INTENT(IN) :: D,N,idcl
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: prob(N)
         DOUBLE PRECISION, INTENT(IN) :: x(D,N)
         INTEGER, INTENT(IN) :: clroots(N)
         DOUBLE PRECISION, INTENT(OUT) :: Q(D,D)

         DOUBLE PRECISION :: ww(N),normww

         ! get the norm of the weights
         normww = logsumexp(ngrid,clroots,prob,idcl)

         ! select just the probabilities of the element
         ! that belongs to the cluster tgt
         WHERE (clroots .EQ. idcl)
            ww = DEXP(prob - normww)
         ELSEWHERE
            ww = 0.0d0
         END WHERE

         CALL covariance(D,period,N,1.0d0,ww,x,Q)
      END SUBROUTINE getlcovcluster

      SUBROUTINE getlcovclusterp(D,period,N,Ntot,prob,x,clroots,idcl,Q)
      ! log version
         INTEGER, INTENT(IN) :: D,N,Ntot,idcl
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: prob(N)
         DOUBLE PRECISION, INTENT(IN) :: x(D,N)
         INTEGER, INTENT(IN) :: clroots(N)
         DOUBLE PRECISION, INTENT(OUT) :: Q(D,D)

         DOUBLE PRECISION :: ww(N),totnormp,nlk,R2,Re2,xx(D,N)
         INTEGER :: i
         ! get the total sum of the p(i)
         totnormp = logsumexp(ngrid,clroots/clroots,prob,1)
!WRITE(*,*) "clust , ", idcl
!WRITE(*,*) "totp: ", totnormp

         ! select just the probabilities of the element
         ! that belongs to the cluster target and normalize them by totnormp
         WHERE (clroots .EQ. idcl)
            ww = DEXP(prob - totnormp)
         ELSEWHERE
            ww = 0.0d0
         END WHERE
         ! Rescale the weights to the size of the dataset,
         ! in order to obtain the fraction of samples corresponding to the gridpoint i
         ww = Ntot * ww
         ! get the amount of points that belong to this cluster
         nlk = SUM(ww)
!WRITE(*,*) x(1,1),x(2,1)
!WRITE(*,*) "nk: ", nlk
!DO i=1,D
!WRITE(*,*) x(1,i),DCOS(x(1,i)),ww(i),ww(i)*DCOS(x(1,i))
!ENDDO
!STOP
         Q = 0.0d0
         DO i=1,D
           ! let's apply the M.I.C. to the angle
           xx(i,:) = x(i,:) - DNINT(x(i,:)/period(i)) * period(i)
           ! see the wiki:
           ! https://en.wikipedia.org/wiki/Von_Mises_distribution
           ! Section: Estimation of parameters
           R2 = (SUM(ww*DCOS(xx(i,:)))/nlk)**2 + (SUM(ww*DSIN(xx(i,:)))/nlk)**2
!WRITE(*,*) "dim, rr2: ",i, R2
           ! get the unbiased estimator
           Re2 = (nlk/(nlk-1.0d0))*(R2-(1.0d0/nlk))
!WRITE(*,*) "dim, urr2: ",i, R2
           ! one could iterate this, but this first approximation should already be enough
           ! See: https://en.wikipedia.org/wiki/Von_Mises%E2%80%93Fisher_distribution
           ! Section: Estimation of parameters
           ! and get the inverse of the concetration paramater
           ! that correspond to a variance
           Q(i,i) = 1.0d0/(DSQRT(Re2)*(2.0d0-Re2) / (1.0d0 - Re2))
!WRITE(*,*) "dim, k: ",i, DSQRT(Re2)*(2.0d0-Re2) / (1.0d0 - Re2)
         ENDDO
      END SUBROUTINE getlcovclusterp

      SUBROUTINE getcovcluster(D,period,N,prob,x,clroots,idcl,Q)
      ! non log version
         INTEGER, INTENT(IN) :: D,N,idcl
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: prob(N)
         DOUBLE PRECISION, INTENT(IN) :: x(D,N)
         INTEGER, INTENT(IN) :: clroots(N)
         DOUBLE PRECISION, INTENT(OUT) :: Q(D,D)

         DOUBLE PRECISION :: ww(N)

         ! select just the probabilities of the element
         ! that belongs to the cluster tgt
         ww=0.0d0
         WHERE (clroots .EQ. idcl)
            ww = prob
         END WHERE

         CALL covariance(D,period,N,SUM(ww),ww,x,Q)
      END SUBROUTINE getcovcluster

      SUBROUTINE readinput(D, fweight, nsamples, xj, totw, wj)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         LOGICAL, INTENT(IN) :: fweight
         INTEGER, INTENT(OUT) :: nsamples
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: xj
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: wj
         DOUBLE PRECISION, INTENT(OUT) :: totw

         ! uses a buffer to read the input reallocating the arrays when needed
         INTEGER, PARAMETER :: nbuff = 100000
         DOUBLE PRECISION :: vbuff(D,nbuff), wbuff(nbuff)
         DOUBLE PRECISION, ALLOCATABLE :: vtmp(:,:), wtmp(:)

         INTEGER io_status, counter

         nsamples = 0
         totw = 0.0d0
         counter = 0

         ! initial dummy allocation
         ALLOCATE(xj(D,1),wj(1),vtmp(D,1),wtmp(1))
         xj=0.0d0
         DO
            IF(fweight) THEN
               READ(5,*, IOSTAT=io_status) vbuff(:,counter+1), wbuff(counter+1)
            ELSE
               READ(5,*, IOSTAT=io_status) vbuff(:,counter+1)
            ENDIF

            IF(io_status<0 .or. io_status==5008) EXIT    ! also intercepts a weird error given by some compilers when reading past of EOF
            IF(io_status>0) STOP "*** Error occurred while reading file. ***"

            IF(fweight) THEN
               totw=totw+wbuff(counter+1)
            ELSE
               wbuff(counter+1)=1.0d0
               totw=totw+wbuff(counter+1)
            ENDIF

            counter=counter+1

            ! grow the arrays and dump the buffers
            IF(counter.EQ.nbuff) THEN
               DEALLOCATE(wtmp,vtmp)
               ALLOCATE(wtmp(nsamples+counter), vtmp(D,nsamples+counter))
               wtmp(1:nsamples) = wj
               vtmp(:,1:nsamples) = xj
               wtmp(nsamples+1:nsamples+counter) = wbuff
               vtmp(:,nsamples+1:nsamples+counter) = vbuff

               DEALLOCATE(wj, xj)
               ALLOCATE(wj(nsamples+counter), xj(D,nsamples+counter))
               wj=wtmp
               xj=vtmp

               nsamples=nsamples+counter
               counter=0
            ENDIF
         END DO

         IF(counter>0) THEN
            DEALLOCATE(wtmp,vtmp)
            ALLOCATE(wtmp(nsamples+counter), vtmp(D,nsamples+counter))
            wtmp(1:nsamples) = wj
            vtmp(:,1:nsamples) = xj
            wtmp(nsamples+1:nsamples+counter) = wbuff(1:counter)
            vtmp(:,nsamples+1:nsamples+counter) = vbuff(:,1:counter)

            DEALLOCATE(wj, xj)
            ALLOCATE(wj(nsamples+counter), xj(D,nsamples+counter))
            wj=wtmp
            xj=vtmp

            nsamples=nsamples+counter
            counter=0
         ENDIF
      END SUBROUTINE readinput

      SUBROUTINE readinputprobs(D, ngrid, yy, prb, ae, re, rgr)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(OUT) :: ngrid
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: yy
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: prb
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: ae
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: re
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: rgr
         ! uses a buffer to read the input reallocating the arrays when needed
         INTEGER, PARAMETER :: nbuff = 30000
         DOUBLE PRECISION :: vbuff(D,nbuff), prbuff(nbuff), aebuff(nbuff)
         DOUBLE PRECISION :: rebuff(nbuff), rgrbuff(nbuff)
         DOUBLE PRECISION, ALLOCATABLE :: vtmp(:,:), prtmp(:)
         DOUBLE PRECISION, ALLOCATABLE :: aetmp(:), retmp(:), rgrtmp(:)
         DOUBLE PRECISION tmparray(4)

         INTEGER io_status, counter

         ngrid = 0
         counter = 0
         tmparray = 0.0d0
         ! initial dummy allocation
         IF (ALLOCATED(yy)) DEALLOCATE(yy)
         IF (ALLOCATED(prb)) DEALLOCATE(prb)
         IF (ALLOCATED(ae)) DEALLOCATE(ae)
         IF (ALLOCATED(re)) DEALLOCATE(re)
         IF (ALLOCATED(rgr)) DEALLOCATE(rgr)
         ALLOCATE(yy(D,1),prb(1),ae(1),re(1),rgr(1),vtmp(D,1))
         yy=0.0d0
         DO
            READ(5,*, IOSTAT=io_status) vbuff(:,counter+1), tmparray(:)
            prbuff(counter+1)  = tmparray(1)
            aebuff(counter+1)  = tmparray(2)
            rebuff(counter+1)  = tmparray(3)
            rgrbuff(counter+1) = tmparray(4)*tmparray(4)

            IF(io_status<0 .or. io_status==5008) EXIT    ! also intercepts a weird error given by some compilers when reading past of EOF
            IF(io_status>0) STOP "*** Error occurred while reading file. ***"

            counter=counter+1

            ! grow the arrays and dump the buffers
            IF(counter.EQ.nbuff) THEN
               IF (ALLOCATED(vtmp)) DEALLOCATE(vtmp)
               IF (ALLOCATED(prtmp)) DEALLOCATE(prtmp)
               IF (ALLOCATED(aetmp)) DEALLOCATE(aetmp)
               IF (ALLOCATED(retmp)) DEALLOCATE(retmp)
               IF (ALLOCATED(rgrtmp)) DEALLOCATE(rgrtmp)
               ALLOCATE(vtmp(D,ngrid+counter),prtmp(ngrid+counter))
               ALLOCATE(aetmp(ngrid+counter),retmp(ngrid+counter))
               ALLOCATE(rgrtmp(ngrid+counter))
               vtmp(:,1:ngrid) = yy
               prtmp(1:ngrid)  = prb
               aetmp(1:ngrid)  = ae
               retmp(1:ngrid)  = re
               rgrtmp(1:ngrid) = rgr
               vtmp(:,ngrid+1:ngrid+counter)   = vbuff
               prtmp(ngrid+1:ngrid+counter)    = prbuff
               aetmp(ngrid+1:ngrid+counter)    = aebuff
               retmp(ngrid+1:ngrid+counter)    = rebuff
               rgrtmp(ngrid+1:ngrid+counter)   = rgrbuff

               DEALLOCATE(yy, prb, ae, re, rgr)
               ALLOCATE(prb(ngrid+counter), yy(D,ngrid+counter))
               yy  = vtmp
               prb = prtmp
               ae  = aetmp
               re  = retmp
               rgr = rgrtmp

               ngrid=ngrid+counter
               counter=0
            ENDIF
         END DO

         IF(counter>0) THEN
            IF (ALLOCATED(vtmp)) DEALLOCATE(vtmp)
            IF (ALLOCATED(prtmp)) DEALLOCATE(prtmp)
            IF (ALLOCATED(aetmp)) DEALLOCATE(aetmp)
            IF (ALLOCATED(retmp)) DEALLOCATE(retmp)
            IF (ALLOCATED(rgrtmp)) DEALLOCATE(rgrtmp)
            ALLOCATE(vtmp(D,ngrid+counter),prtmp(ngrid+counter))
            ALLOCATE(aetmp(ngrid+counter),retmp(ngrid+counter))
            ALLOCATE(rgrtmp(ngrid+counter))
            vtmp(:,1:ngrid) = yy
            prtmp(1:ngrid)  = prb
            aetmp(1:ngrid)  = ae
            retmp(1:ngrid)  = re
            rgrtmp(1:ngrid) = rgr
            vtmp(:,ngrid+1:ngrid+counter)   = vbuff(:,1:counter)
            prtmp(ngrid+1:ngrid+counter)    = prbuff(1:counter)
            aetmp(ngrid+1:ngrid+counter)    = aebuff(1:counter)
            retmp(ngrid+1:ngrid+counter)    = rebuff(1:counter)
            rgrtmp(ngrid+1:ngrid+counter)   = rgrbuff(1:counter)

            DEALLOCATE(yy, prb, ae, re, rgr)
            ALLOCATE(prb(ngrid+counter), yy(D,ngrid+counter))
            ALLOCATE(ae(ngrid+counter), re(ngrid+counter))
            ALLOCATE(rgr(ngrid+counter))
            yy  = vtmp
            prb = prtmp
            ae  = aetmp
            re  = retmp
            rgr = rgrtmp

            ngrid=ngrid+counter
            counter=0
         ENDIF
      END SUBROUTINE readinputprobs

      SUBROUTINE mkgrid(D,period,nsamples,ngrid,x,wj,y,ni,iminij, &
                        ineigh,wi,saveidx,idxgrid,ofile)
         ! Select ngrid grid points from nsamples using minmax and
         ! the voronoi polyhedra around them.
         !
         ! Args:
         !    nsamples: total points number
         !    ngrid: number of grid points
         !    x: array containing the data samples
         !    y: array that will contain the grid points
         !    ni: array cotaing the number of samples inside the Voronoj polyhedron of each grid point
         !    iminij: array containg the neighbor list for data samples

         INTEGER, INTENT(IN) :: D
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         INTEGER, INTENT(IN) :: nsamples
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, DIMENSION(D,nsamples), INTENT(IN) :: x
         DOUBLE PRECISION, DIMENSION(nsamples), INTENT(IN) :: wj

         DOUBLE PRECISION, DIMENSION(D,ngrid), INTENT(OUT) :: y
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: ni
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: ineigh
         INTEGER, DIMENSION(nsamples), INTENT(OUT) :: iminij
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(OUT) :: wi
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: idxgrid
         CHARACTER(LEN=1024), INTENT(IN) :: ofile
         LOGICAL, INTENT(IN) :: saveidx

         INTEGER i,j,irandom
         DOUBLE PRECISION :: dminij(nsamples), dij, dmax, dneigh

         iminij=0
         y=0.0d0
         ni=0
         ! choose randomly the first point
         irandom=int(RAND()*nsamples)
         idxgrid(1)=irandom
         IF(saveidx) THEN
            OPEN(UNIT=12,FILE=trim(ofile)//".idxs",STATUS='REPLACE',ACTION='WRITE')
            WRITE(12,"((I9))") irandom
         ENDIF
         y(:,1)=x(:,irandom)
         dminij = HUGE(0.0d0)
         iminij = 1
         ineigh = 0
         DO i=2,ngrid
            dmax = 0.0d0
            dneigh = HUGE(0.0d0)
            DO j=1,nsamples
               dij = pammr2(D,period,y(:,i-1),x(:,j))
               IF (dminij(j)>dij) THEN
                  dminij(j) = dij
                  iminij(j) = i-1 ! also keeps track of the Voronoi attribution
               ENDIF
               IF (dminij(j) > dmax) THEN
                  dmax = dminij(j)
                  jmax = j
               ENDIF
               IF ((dneigh > dij) .and. (dij .ne. 0.0d0)) THEN
                  dneigh = dij
                  ineigh(i-1) = j ! store index of closest sample neighbor to grid point
               ENDIF
            ENDDO
            y(:,i) = x(:, jmax)
            IF(saveidx) THEN
               WRITE(12,"((I9))") jmax
            ENDIF
            idxgrid(i)=jmax
            IF(verbose .AND. (modulo(i,1000).EQ.0)) &
               write(*,*) i,"/",ngrid
         ENDDO

         ! finishes Voronoi attribution
         dneigh = HUGE(0.0d0)
         DO j=1,nsamples
            dij = pammr2(D,period,y(:,ngrid),x(:,j))
            IF (dminij(j)>dij) THEN
               dminij(j) = dij
               iminij(j) = ngrid
            ENDIF
            IF ((dneigh > dij) .and. (dij .ne. 0.0d0)) THEN
               dneigh = dij
               ineigh(ngrid) = j ! store index of closest sample neighbor to grid point
            ENDIF
         ENDDO

         ! Assign neighbor list pointer of voronois
         ! Number of points in each voronoi polyhedra
         ni = 0
         wi = 0.0d0
         DO j=1,nsamples
            ni(iminij(j))=ni(iminij(j))+1
            wi(iminij(j))=wi(iminij(j))+wj(j)
         ENDDO

         IF (saveidxs) THEN
            CLOSE(12)
         ENDIF
      END SUBROUTINE mkgrid

      SUBROUTINE getvoro(D,period,nsamples,ngrid,x,wj,y,ni,iminij, &
                         ineigh,wi,idxgrid)
         IMPLICIT NONE

         ! Select ngrid grid points from nsamples using minmax and
         ! the voronoi polyhedra around them.
         !
         ! Args:
         !    nsamples: total points number
         !    ngrid: number of grid points
         !    x: array containing the data samples
         !    y: array that will contain the grid points
         !    ni: array cotaing the number of samples inside the Voronoj polyhedron of each grid point
         !    iminij: array containg the neighbor list for data samples

         INTEGER, INTENT(IN) :: D
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         INTEGER, INTENT(IN) :: nsamples
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, DIMENSION(D,nsamples), INTENT(IN) :: x
         DOUBLE PRECISION, DIMENSION(nsamples), INTENT(IN) :: wj
         INTEGER, DIMENSION(ngrid), INTENT(IN) :: idxgrid

         DOUBLE PRECISION, DIMENSION(D,ngrid), INTENT(OUT) :: y
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: ni
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: ineigh
         INTEGER, DIMENSION(nsamples), INTENT(OUT) :: iminij
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(OUT) :: wi

         INTEGER i,j,irandom
         DOUBLE PRECISION :: dminij(nsamples), dij, dmax, dneigh

         iminij=0
         y=0.0d0
         ni=0
         ! start from first point of user provided grid
         irandom=idxgrid(1)
         y(:,1)=x(:,irandom)
         dminij = HUGE(0.0d0)
         iminij = 1
         ineigh = 0
         DO i=2,ngrid
            dmax = 0.0d0
            dneigh = HUGE(0.0d0)
            DO j=1,nsamples
               dij = pammr2(D,period,y(:,i-1),x(:,j))
               IF (dminij(j)>dij) THEN
                  dminij(j) = dij
                  iminij(j) = i-1 ! also keeps track of the Voronoi attribution
               ENDIF
               IF (dminij(j) > dmax) THEN
                  dmax = dminij(j)
               ENDIF
               IF ((dneigh > dij) .and. (dij .ne. 0.0d0)) THEN
                  dneigh = dij
                  ineigh(i-1) = j ! store index of closest sample neighbor to grid point
               ENDIF
            ENDDO
            y(:,i) = x(:, idxgrid(i))
            IF(verbose .AND. (modulo(i,1000).EQ.0)) &
               write(*,*) i,"/",ngrid
         ENDDO

         ! finishes Voronoi attribution
         dneigh = HUGE(0.0d0)
         DO j=1,nsamples
            dij = pammr2(D,period,y(:,ngrid),x(:,j))
            IF (dminij(j)>dij) THEN
               dminij(j) = dij
               iminij(j) = ngrid
            ENDIF
            IF ((dneigh > dij) .and. (dij .ne. 0.0d0)) THEN
               dneigh = dij
               ineigh(ngrid) = j ! store index of closest sample neighbor to grid point
            ENDIF
         ENDDO

         ! Assign neighbor list pointer of voronois
         ! Number of points in each voronoi polyhedra
         ni = 0
         wi = 0.0d0
         DO j=1,nsamples
            ni(iminij(j))=ni(iminij(j))+1
            wi(iminij(j))=wi(iminij(j))+wj(j)
         ENDDO
      END SUBROUTINE getvoro

      SUBROUTINE getnlist(nsamples,ngrid,ni,iminij, pnlist,nlist)
         ! Build a neighbours list: for every voronoi center keep track of his
         ! neighboroud that correspond to all the points inside the voronoi
         ! polyhedra.
         !
         ! Args:
         !    nsamples: total points number
         !    ngrid: number of voronoi polyhedra
         !    weights: array cotaing the number of points inside each voroni polyhedra
         !    iminij: array containg to wich polyhedra every point belong to
         !    pnlist: pointer to neighbours list
         !    nlist: neighbours list

         INTEGER, INTENT(IN) :: nsamples
         INTEGER, INTENT(IN) :: ngrid
         INTEGER, DIMENSION(ngrid), INTENT(IN) :: ni
         INTEGER, DIMENSION(nsamples), INTENT(IN) :: iminij
         INTEGER, DIMENSION(ngrid+1), INTENT(OUT) :: pnlist
         INTEGER, DIMENSION(nsamples), INTENT(OUT) :: nlist

         INTEGER i,j
         INTEGER :: tmpnidx(ngrid)

         pnlist=0
         nlist=0
         tmpnidx=0

         ! pointer to the neighbourlist
         pnlist(1)=0
         DO i=1,ngrid
            pnlist(i+1)=pnlist(i)+ni(i)
            tmpnidx(i)=pnlist(i)+1  ! temporary array to use while filling up the neighbour list
         ENDDO

         DO j=1,nsamples
            i=iminij(j) ! this is the Voronoi center the sample j belongs to
            nlist(tmpnidx(i))=j ! adds j to the neighbour list
            tmpnidx(i)=tmpnidx(i)+1 ! advances the pointer
         ENDDO
      END SUBROUTINE getnlist

      DOUBLE PRECISION FUNCTION cls_link(ngrid, idcls, distmm, prob, rgrid, ia, ib, &
                                         errors, linknoerr)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngrid, idcls(ngrid), ia, ib
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid, ngrid), prob(ngrid), rgrid(ngrid)
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: errors
         DOUBLE PRECISION, INTENT(OUT) :: linknoerr

         INTEGER i, j
         DOUBLE PRECISION mxa, mxb, mxab, pab, emxa, emxb, emxab, g1, g2
         mxa   = 0.0d0
         mxb   = 0.0d0
         mxab  = 0.0d0
         emxa  = 0.0d0
         emxb  = 0.0d0
         emxab = 0.0d0

         linknoerr = 0.0d0
         DO i=1,ngrid
            IF (idcls(i)/=ia) CYCLE
            IF (prob(i).gt.mxa) THEN
               mxa  = prob(i)    ! also gets the probability density at the mode of cluster a
               emxa = errors(i) ! and the absolute error associated
            ENDIF
            DO j=1,ngrid
               IF (idcls(j)/=ib) CYCLE
               IF (prob(j).gt.mxb) THEN
                  mxb  = prob(j)
                  emxb = errors(j)
               ENDIF
               ! Ok, we've got a matching pair
               IF (dsqrt(distmm(i,j))<dsqrt(rgrid(i))+dsqrt(rgrid(j))) THEN
                  ! And they are close together!
                  pab = (prob(i)+prob(j))/2
                  IF (pab .gt. mxab) THEN
                     mxab  = pab
                     emxab = DSQRT((errors(i))**2+(errors(j))**2)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO


         IF(mxab.EQ.0)THEN
            cls_link = 0.0d0
         ELSE
            g1 = (mxab+emxab)/min(max(mxa-emxa,0.0d0),max(mxb-emxb,0.0d0))
            g2 = (mxab-emxab)/min(max(mxa+emxa,0.0d0),max(mxb+emxb,0.0d0))
            cls_link = min(1.0d0,max(g1,g2))
            linknoerr = mxab/min(mxa,mxb)
         ENDIF
      END FUNCTION


      INTEGER FUNCTION qs_next(ngrid,idx,idxn,probnmm,distmm,lambda)
         ! Return the index of the closest point higher in P
         !
         ! Args:
         !    ngrid: number of grid points
         !    idx: current point
         !    qscut: cut-off in the jump
         !    probnmm: density estimations
         !    distmm: distances matrix

         INTEGER, INTENT(IN) :: ngrid
         INTEGER, INTENT(IN) :: idx,idxn
         DOUBLE PRECISION, INTENT(IN) :: lambda
         DOUBLE PRECISION, INTENT(IN) :: probnmm(ngrid)
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid,ngrid)

         INTEGER j
         DOUBLE PRECISION dmin

         dmin = HUGE(0.0d0)

         qs_next = idx
         If (probnmm(idxn).GT.probnmm(idx) ) THEN
           qs_next=idxn
         ENDIF
         DO j=1,ngrid
            IF ( probnmm(j).GT.probnmm(idx) ) THEN
               IF ((distmm(idx,j).LT.dmin) .AND. (distmm(idx,j).LT.lambda)) THEN
                 dmin = distmm(idx,j)
                 qs_next = j
               ENDIF
            ENDIF
         ENDDO
      END FUNCTION qs_next

      INTEGER FUNCTION gs_next(ngrid,idx,probnmm,distmm,gabriel,nn)
         ! Return the index of the closest point higher in P
         !
         ! Args:
         !    ngrid: number of grid points
         !    idx: current point
         !    probnmm: density estimations
         !    distmm: distance matrix (squared)
         !    gabriel: gabriel graph
         !    nn: cut-off in the jumps

         INTEGER, INTENT(IN) :: ngrid
         INTEGER, INTENT(IN) :: idx
         INTEGER, INTENT(IN) :: nn                      ! number of neighbor shells
         DOUBLE PRECISION, INTENT(IN) :: probnmm(ngrid)
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid,ngrid)
         LOGICAL, INTENT(IN) :: gabriel(ngrid,ngrid)

         INTEGER i,j
         DOUBLE PRECISION dmin
         LOGICAL neighs(ngrid), nneighs(ngrid)

         ! neighbors and neighbors of neighbors and ...
         neighs = gabriel(idx,:)
         ! loop over the nn neighbor shells
         DO i=2,nn
            nneighs = .FALSE.
            DO j=1,ngrid
              IF (neighs(j)) nneighs = nneighs .OR. gabriel(j,:)
            ENDDO
            ! add new neighbors to neighbor array
            neighs = neighs .OR. nneighs
         ENDDO

         gs_next = idx
         dmin = HUGE(0.0d0)
         DO j=1,ngrid
            IF ( probnmm(j).GT.probnmm(idx) ) THEN
               IF ( (distmm(idx,j).LT.dmin) .AND. neighs(j) ) THEN
                  gs_next = j
                  dmin = distmm(idx,j)
               ENDIF
            ENDIF
         ENDDO
      END FUNCTION gs_next


!      INTEGER FUNCTION qs_next(D,period,N,i,cutoff,prob,M,y,multi)
!         ! Return the index of the closest point higher in P
!         !
!         ! Args:
!         !    N:     number of grid points
!         !    i:     index of current point
!         !    cut:   spherical cut-off in the jump
!         !    prob:  densities on grid
!         !    M:     distances matrix
!         !           upper triangular is mahalanobis distance using spherical covariance
!         !           lower triangular is mahalanobis distance using local covariance

!         INTEGER, INTENT(IN) :: D
!         DOUBLE PRECISION, INTENT(IN) :: period(D)
!         INTEGER, INTENT(IN) :: N
!         INTEGER, INTENT(IN) :: i
!         DOUBLE PRECISION, INTENT(IN) :: cutoff
!         DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: prob
!         DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: M
!         DOUBLE PRECISION, DIMENSION(D,N), INTENT(IN) :: y
!         DOUBLE PRECISION, INTENT(IN) :: multi
!
!         INTEGER j
!         DOUBLE PRECISION dmin
!
!         ! set dmin to highest possible 64-bit double
!         dmin = 1.0d308
!         ! inverse of the spherical cutoff

!         qs_next = i
!         DO j=1,N
!            IF ( prob(j).GT.prob(i) ) THEN
!               IF ((pammr2(D,period,y(i,:),y(j,:))/cutoff).LT.multi) THEN
!                 IF (M(j,i).LT.dmin) THEN
!                   WRITE(*,*) "QS: ", i,j, pammr2(D,period,y(i,:),y(j,:))/cutoff,M(j,i)
!                   dmin = M(j,i)
!                   qs_next = j
!                 ENDIF
!               ENDIF
!            ENDIF
!         ENDDO
!
!!         qs_next = i
!!         DO j=1,N
!!            IF ( prob(j).GT.prob(i) ) THEN
!!               IF ( (M(j,i).LT.dmin) .AND. (M(j,i).LT.cutoff) ) THEN
!!                  dmin = M(j,i)
!!                  qs_next = j
!!               ENDIF
!!            ENDIF
!!         ENDDO
!
!      END FUNCTION qs_next

      DOUBLE PRECISION FUNCTION fmultiVM(D,dlocal,period,x,y,icov,cov)
         ! Return the multivariate gaussian density
         ! Args:
         !    gpars: gaussian parameters
         !    x: point in wich estimate the value of the gaussian

         INTEGER , INTENT(IN) :: D
         DOUBLE PRECISION, INTENT(IN) :: dlocal
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: x(D)
         DOUBLE PRECISION, INTENT(IN) :: y(D)
         DOUBLE PRECISION, INTENT(IN) :: icov(D,D),cov(D,D)
         DOUBLE PRECISION dv(D)

         DOUBLE PRECISION dumm,dumm1,ev(D)
         INTEGER jj,effD

         !! Here if the concetration parameter is big enaugh, then
         !! the Vm distrib can be seen as a gaussian..
         !! Let's exploit this

         ! check the diagonal of Hi: if the biggest element
         ! is smaller than 0.6, we can safely use a multivariate gaussian


         dumm=0.0d0
         DO jj=1,D
            dumm1=DSQRT(cov(jj,jj))
            IF(dumm.GT.dumm1) dumm=dumm1
            ev(jj)=cov(jj,jj)
         ENDDO

         ! sort the diagonal
         CALL sort(ev, D)
         ! get the local dimensionality
         effD=NINT(REAL(dlocal))

         IF(dumm.LT.0.6d0)THEN
            ! TODO: removed periodicity from multivariate kernel, need to replace this function
            fmultiVM=fmultikernel(D,x,y,icov,1.0d0/DSQRT((twopi**DBLE(D))*detmatrix(D,cov)))

         ELSE
            ! productkernels
            fmultiVM=1.0d0
            CALL pammrij(D, period, x, y, dv)
            DO i = D, (D-effD), -1
               IF(DSQRT(ev(jj)).LT.0.6d0)THEN
                 fmultiVM = fmultiVM * (1.0d0/((twopi*ev(jj))**0.5d0))* &
                       DEXP(-0.5d0*(dv(jj)**2.0d0)/ev(jj))
               ELSE
                 fmultiVM = fmultiVM * fvmkernel(1.0d0/ev(jj),dv(jj))
               ENDIF
            END DO
         ENDIF

      END FUNCTION fmultiVM

      DOUBLE PRECISION FUNCTION fmultikernel(D,x,y,icov,norm)
         ! Return the multivariate gaussian density
         ! Args:
         !    gpars: gaussian parameters
         !    x: point in wich estimate the value of the gaussian

         INTEGER , INTENT(IN) :: D
         DOUBLE PRECISION, INTENT(IN) :: x(D)
         DOUBLE PRECISION, INTENT(IN) :: y(D)
         DOUBLE PRECISION, INTENT(IN) :: icov(D,D)
         DOUBLE PRECISION, INTENT(IN) :: norm
         DOUBLE PRECISION dv(D),tmpv(D),xcx

         dv = x - y
         tmpv = MATMUL(dv,icov)
         xcx = -0.5d0 * DOT_PRODUCT(dv,tmpv)

         fmultikernel = DEXP(xcx) * norm

      END FUNCTION fmultikernel

!      DOUBLE PRECISION FUNCTION fkernel(D,period,sig2,vc,vp)
!            ! Calculate the gaussian kernel
!            ! The normalization has to be done outside
!            !
!            ! Args:
!            !    sig2: sig**2
!            !    vc: voronoi center's vector
!            !    vp: point's vector

!            INTEGER, INTENT(IN) :: D
!            DOUBLE PRECISION, INTENT(IN) :: period(D)
!            DOUBLE PRECISION, INTENT(IN) :: sig2
!            DOUBLE PRECISION, INTENT(IN) :: vc(D)
!            DOUBLE PRECISION, INTENT(IN) :: vp(D)


!            fkernel=(1/( (twopi*sig2)**(dble(D)/2) ))* &
!                    dexp(-pammr2(D,period,vc,vp)*0.5/sig2)
!
!      END FUNCTION fkernel

      DOUBLE PRECISION FUNCTION fvmkernel(kkk,dist)
            ! Calculate the univariate von Mises kernel
            !
            ! Args:
            !    sig2: sig**2
            !    dist: distance between the two points

            DOUBLE PRECISION, INTENT(IN) :: kkk
            DOUBLE PRECISION, INTENT(IN) :: dist

            fvmkernel=DEXP(DCOS(dist)*kkk) / &
                      (BESSI0(kkk)*twopi)

      END FUNCTION fvmkernel

      SUBROUTINE savevoronois(nsamples,iminij,prvor)
         ! Store Voronoi data in a file
         !
         ! Args:
         !    nsamples   : total points number
         !    iminij     : voronoi link
         !    prvor      : prefix for the outputfile

         INTEGER, INTENT(IN) :: nsamples
         INTEGER, DIMENSION(nsamples), INTENT(IN) :: iminij
         CHARACTER(LEN=1024), INTENT(IN) :: prvor

         INTEGER j

         ! write out the voronoi links
         OPEN(UNIT=12,FILE=trim(prvor)//".voronoislinks",STATUS='REPLACE',ACTION='WRITE')
         ! header
         WRITE(12,*) "# sample point , voronoi association"

         DO j=1,nsamples
            ! write the
            WRITE(12,"((A1,I9))",ADVANCE="NO") " ", j
            ! write the Voronoi associated
            WRITE(12,"((A1,I9))") " ", iminij(j)
         ENDDO

         CLOSE(UNIT=12)

      END SUBROUTINE savevoronois

      SUBROUTINE savegrid(D,ngrid,y,prb,aer,rer,rgr,ofile)
         ! Store Voronoi data in a file
         !
         ! Args:
         !    D          : Dimensionality of a point
         !    ngrid      : number of grid points
         !    y          : grid points
         !    prb        : KDE
         !    aer        : absolute errors on the KDE
         !    rer        : relative errors on the KDE
         !    rgr        : rgrid distances (square distances)
         !    ofile      : prefix for the outputfile

         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, INTENT(IN) :: y(D,ngrid)
         DOUBLE PRECISION, INTENT(IN) :: prb(ngrid)
         DOUBLE PRECISION, INTENT(IN) :: aer(ngrid)
         DOUBLE PRECISION, INTENT(IN) :: rer(ngrid)
         DOUBLE PRECISION, INTENT(IN) :: rgr(ngrid)
         CHARACTER(LEN=1024), INTENT(IN) :: ofile

         INTEGER i,j

         ! write out the voronoi links
         OPEN(UNIT=12,FILE=trim(ofile)//".probs",STATUS='REPLACE',ACTION='WRITE')

         DO j=1,ngrid
            ! write first the point
            DO i=1,D
               WRITE(12,"((A1,ES15.4E4))",ADVANCE="NO") " ", y(i,j)
            ENDDO
            WRITE(12,"((A1,ES20.8E4))",ADVANCE="NO") " ", prb(j)
            WRITE(12,"((A1,ES20.8E4))",ADVANCE="NO") " ", aer(j)
            WRITE(12,"((A1,ES20.8E4))",ADVANCE="NO") " ", rer(j)
            WRITE(12,"((A1,ES20.8E4))") " ", DSQRT(rgr(j))
         ENDDO

         CLOSE(UNIT=12)
      END SUBROUTINE savegrid

   END PROGRAM pamm
