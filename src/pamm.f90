! This file contain the main program for the PAMM clustering in 
! both PERIODIC and NON PERIODIC space.
! Starting from a set of data points in high dimension it will first perform
! a non-parametric partitioning of the probability density and return the
! Nk multivariate Gaussian/Von Mises distributions better describing the clusters.
! Can also be run in post-processing mode, where it will read tuples and 
! classify them based the model file specified in input.
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
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

      PROGRAM pamm
      USE libpamm
      USE random
      IMPLICIT NONE
      
      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

      CHARACTER(LEN=1024) :: outputfile, clusterfile          ! The output file prefix
      DOUBLE PRECISION, ALLOCATABLE :: period(:)              ! Periodic lenght in each dimension
      LOGICAL periodic                                        ! flag for using periodic data
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: distmm ! similarity matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: diff     ! temp vector used to store distances

      INTEGER D                                               ! Dimensionality of problem
      INTEGER Nk                                              ! Number of gaussians in the mixture
      INTEGER nsamples                                        ! Total number points
      INTEGER ngrid                                           ! Number of samples extracted using minmax

      INTEGER jmax,ii,jj
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma2, rgrid, wj, prob, bigp, &
                                                     msmu, tmpmsmu, pcluster, px, &
                                                     oldsigma2
      DOUBLE PRECISION :: normwj                              ! accumulator for wj
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npvoronoi, iminij, pnlist, nlist
      INTEGER seed                                            ! seed for the random number generator
      INTEGER mode                                            ! KDE mode selection

      ! variable to set the covariance matrix
      DOUBLE PRECISION tmppks,normpks,tmpkernel
      DOUBLE PRECISION tmpad, tmpadb
    
      ! Array of Von Mises distributions
      TYPE(vm_type), ALLOCATABLE, DIMENSION(:) :: vmclusters

      ! Array of Gaussians containing the gaussians parameters
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters
      ! Array containing the input data pints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x , y
      ! quick shift, roots and path to reach the root (used to speedup the calculation)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: idxroot, idcls, qspath
      ! cluster connectivity matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: clsadj
      INTEGER, ALLOCATABLE, DIMENSION(:) :: macrocl,sortmacrocl
      LOGICAL isthere
      
      ! BOOTSTRAP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: probboot
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: proberr
      INTEGER, ALLOCATABLE, DIMENSION (:) :: nbsisel
      INTEGER nbootstrap, rndidx, rngidx, nn, nbssample, nbstot
      DOUBLE PRECISION tmperr, tmpcheck, qserr, concentrationK
      ! Variables for local bandwidth estimation
      DOUBLE PRECISION nlocal, lfac, lfac2, prefac, R2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xij, normkernel, wQ, wQgrid
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Q, Qlocal, Sw, Sinv, xtmp, xtmpw, ytmp,ytmpw
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Hi, Hiinv
      
      
      ! IN/OUT probs
      LOGICAL saveprobs, savevor, skipvoronois, pilotboot, saveadj
      
      ! PARSER
      CHARACTER(LEN=1024) :: cmdbuffer, comment   ! String used for reading text lines from files
      INTEGER ccmd                                ! Index used to control the input parameters
      INTEGER nmsopt                              ! number of mean-shift optimizations of the cluster centers
      LOGICAL verbose, fpost                      ! flag for verbosity
      LOGICAL weighted                            ! flag for using weigheted data
      INTEGER isep1, isep2, par_count             ! temporary indices for parsing command line arguments
      INTEGER neblike                             ! iterations for neblike path search
      DOUBLE PRECISION lambda, msw, alpha, zeta, thrmerg, dummd1, dummd2, convchk, lambda2

      ! DOUBLE PRECISION maxrgrid, minrgrid
      DOUBLE PRECISION mixbeta

      INTEGER i,j,k,m,n,ikde,counter,dummyi1,dummyi2,endf ! Counters and dummy variable

      DOUBLE PRECISION  bigperr, tmpbigp, gnorm
!!!!!!! Default value of the parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      outputfile="out"
      clusterfile="NULL"
      fpost=.false.
      alpha=1.0d0
      zeta=0.0d0
      pilotboot= .false.
      ccmd=0                 ! no parameters specified
      Nk=0                   ! number of gaussians
      nmsopt=0               ! number of mean-shift refinements
      ngrid=-1               ! number of samples extracted with minmax
      seed=12345             ! seed for the random number generator
      thrmerg=0.3d0          ! merge different clusters
      lambda=-1              ! quick shift cut-off
      verbose = .false.      ! no verbosity
      weighted= .false.      ! don't use the weights
      lfac = -1.0d0          ! localization factor
      prefac = 0.0d0         ! prefector of gaussian for localization
      neblike=0              ! don't use neb paths
      nbootstrap=0           ! do not use bootstrap
      qserr=8                ! threshold to accept a move in qs
      saveprobs= .false.     ! don't print out the probs
      savevor  = .false.     ! don't print out the Voronoi
      skipvoronois = .false. ! don't read the Voronoi associations 
      saveadj = .false.      ! save adjacency
            
      
      mixbeta = 0.1
      
      D=-1
      periodic=.false.
      CALL random_init(seed) ! initialize random number generator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-a") THEN       ! cluster smearing
            ccmd = 1
         ELSEIF (cmdbuffer == "-o") THEN           ! output file
            ccmd = 2
         ELSEIF (cmdbuffer == "-gf") THEN      ! file containing Vn parmeters
            ccmd = 3
         ELSEIF (cmdbuffer == "-seed") THEN    ! seed for the random number genarator
            ccmd = 4
         ELSEIF (cmdbuffer == "-l") THEN       ! threshold to differentiate different clusters
                ccmd = 5
         ELSEIF (cmdbuffer == "-nms") THEN     ! N of mean-shift steps
            ccmd = 6
         ELSEIF (cmdbuffer == "-ngrid") THEN   ! N of grid points
            ccmd = 7
         ELSEIF (cmdbuffer == "-bootstrap") THEN   ! refine the kde using bootstrapping
            ccmd = 8
         ELSEIF (cmdbuffer == "-d") THEN       ! dimensionality
            ccmd = 9
         ELSEIF (cmdbuffer == "-p") THEN       ! use periodicity
            ccmd = 11
         ELSEIF (cmdbuffer == "-z") THEN       ! add a background to the probability mixture
            ccmd = 12
         ELSEIF (cmdbuffer == "-neblike") THEN  ! use neb paths between the qs points
            ccmd = 13
         ELSEIF (cmdbuffer == "-loc") THEN  ! refine adptively sigma
             ccmd = 14
         ELSEIF (cmdbuffer == "-readprobs") THEN  ! read a file containing the grid points + info
            ccmd = 15
         ELSEIF (cmdbuffer == "-readvoronis") THEN  ! read a file containing the Voronoi associations
            ccmd = 16
         ELSEIF (cmdbuffer == "-skipvoronois") THEN  ! save the KDE estimates in a file
            skipvoronois= .true.
         ELSEIF (cmdbuffer == "-saveprobs") THEN  ! save the KDE estimates in a file
            saveprobs= .true.
         ELSEIF (cmdbuffer == "-savevoronois") THEN  ! save the Voronoi associations
            savevor= .true.
         ELSEIF (cmdbuffer == "-adj") THEN  ! save the Voronoi associations
            saveadj= .true.
            ccmd = 18
         ELSEIF (cmdbuffer == "-w") THEN       ! use weights
            weighted = .true.
         ELSEIF (cmdbuffer == "-v") THEN       ! verbosity flag
            verbose = .true.
         ELSEIF (cmdbuffer == "-h") THEN       ! help flag
            CALL helpmessage
            CALL EXIT(-1)
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " No parameters specified!"
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 1) THEN            ! read the cluster smearing
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 2) THEN            ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN            ! model file
               fpost=.true.
               clusterfile=trim(cmdbuffer)
            ELSEIF (ccmd == 1) THEN            ! read the cluster smearing
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 18) THEN           ! read the cluster smearing
               READ(cmdbuffer,*) thrmerg
               IF (thrmerg<0) STOP "Put a positive number!"
            ELSEIF (ccmd ==5) THEN ! read the seed
               READ(cmdbuffer,*) lambda
            ELSEIF (ccmd == 4) THEN ! localization parameter
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 14) THEN ! localization parameter
               READ(cmdbuffer,*) lfac
            ELSEIF (ccmd == 6) THEN ! read the number of mean-shift steps
               READ(cmdbuffer,*) nmsopt
            ELSEIF (ccmd == 7) THEN ! number of grid points
               READ(cmdbuffer,*) ngrid
            ELSEIF (ccmd == 8) THEN            ! read the num of bootstrap iterations
               READ(cmdbuffer,*) nbootstrap
               IF (nbootstrap<0) STOP "The number of iterations should be positive!"
            ELSEIF (ccmd == 9) THEN            ! read the dimensionality
               READ(cmdbuffer,*) D
               ALLOCATE(period(D))
               period=-1.0d0
            ELSEIF (ccmd == 11) THEN ! read the periodicity in each dimension
               IF (D<0) STOP "Dimensionality (-d) must precede the periodic lenghts (-p). "
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
            ELSEIF (ccmd == 12) THEN ! read zeta 
               READ(cmdbuffer,*) zeta
            ELSEIF (ccmd == 13) THEN ! activate neb like behaviour 
               READ(cmdbuffer,*) neblike
            ELSEIF (ccmd == 14) THEN ! read threshold for the qs assignation 
               READ(cmdbuffer,*) qserr
            ENDIF
         ENDIF
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL SRAND(seed) ! initialize the random number generator

      ! Check the input parameters
      IF (ccmd.EQ.0) THEN
         ! ccmd==0 : no parameters inserted
         WRITE(*,*) ""
         WRITE(*,*) &
  " Could not interpret command line arguments. Please check for errors!"
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
       
      ! POST-PROCESSING MODE
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
              !!dummyi1=1
              !!DO i=1,nk
              !!   IF (pcluster(i)>pcluster(dummyi1)) dummyi1=i
              !!ENDDO
              !!WRITE(*,*) px,dummyi1
              WRITE(*,*) px,pcluster(:)
            ENDDO
            DEALLOCATE(vmclusters)
         ELSE
            ! NON-PERIODIC version
            CALL readclusters(12,nk,clusters)
            CLOSE(12)
            ALLOCATE(pcluster(nk), px(clusters(1)%D))
            
            DO WHILE (.true.) 
              READ(*,*,IOSTAT=endf) px
              IF(endf>0) STOP "*** Error occurred while reading file. ***"
              IF(endf<0) EXIT
              CALL pamm_p(px, pcluster, nk, clusters, alpha, zeta)
              dummyi1=1
              DO i=1,nk
                 IF (pcluster(i)>pcluster(dummyi1)) dummyi1=i
              ENDDO
              ! write out the number of the 
              ! cluster with the highest probability
              WRITE(*,*) px,dummyi1 ! ,pcluster(dummyi1)
            ENDDO
            DEALLOCATE(clusters)
         ENDIF
         
         DEALLOCATE(pcluster)
         ! done, go home
         CALL EXIT(-1)
      ENDIF
      
      ! CLUSTERING MODE
      ! get the data from standard input
      CALL readinput(D, weighted, nsamples, x, normwj, wj)
      ! "renormalizes" the weight so we can consider them sort of sample counts
      IF (weighted) THEN  
         wj = wj * nsamples/sum(wj)
      ENDIF

      ! If not specified, the number voronoi polyhedra
      ! are set to the square of the total number of points
      IF (ngrid.EQ.-1) ngrid=int(sqrt(float(nsamples)))

      ! Initialize the arrays, since now I know the number of
      ! points and the dimensionality
      ALLOCATE(iminij(nsamples), nbsisel(nsamples))
      ALLOCATE(pnlist(ngrid+1), nlist(nsamples))
      ALLOCATE(y(D,ngrid), npvoronoi(ngrid), prob(ngrid), sigma2(ngrid), rgrid(ngrid))
      ALLOCATE(idxroot(ngrid), idcls(ngrid), qspath(ngrid), distmm(ngrid,ngrid))
      ALLOCATE(diff(D), msmu(D), tmpmsmu(D), oldsigma2(ngrid))
      ALLOCATE(proberr(ngrid),normkernel(ngrid),bigp(ngrid))
      ! bootstrap probability density array will be allocated if necessary
      IF(nbootstrap > 0) ALLOCATE(probboot(ngrid,nbootstrap))
      ! Allocate variables for local bandwidth estimate
      ALLOCATE(Q(D,D),Qlocal(D,D),Hi(D,D,ngrid),Hiinv(D,D,ngrid))
      ALLOCATE(xtmp(D,nsamples),xtmpw(D,nsamples),xij(D),ytmp(D,ngrid),ytmpw(D,ngrid))
      ALLOCATE(wQ(nsamples),Sw(D,D),Sinv(D,D),wQgrid(ngrid))

      ! Extract ngrid points on which the kernel density estimation is to be
      ! evaluated. Also partitions the nsamples points into the Voronoi polyhedra
      ! of the sampling points.
      IF(verbose) THEN
         WRITE(*,*) "NSamples: ", nsamples
         WRITE(*,*) "Selecting ", ngrid, " points using MINMAX"
      ENDIF

      CALL mkgrid(D,period,nsamples,ngrid,x,y,npvoronoi,iminij,wj)
  
      ! Generate the neighbour list
      IF(verbose) write(*,*) "Generating neighbour list"
      CALL getnlist(nsamples,ngrid,npvoronoi,iminij, pnlist,nlist)
      ! Definition of the distance matrix between grid points
      distmm=0.0d0
      rgrid=1d100 ! "voronoi radius" of grid points (squared)
      
      IF(verbose) WRITE(*,*) & 
        "Computing similarity matrix"
      !$omp parallel &
      !$omp default (none) &
      !$omp shared (ngrid,distmm,D,period,y,x,rgrid,verbose) &
      !$omp private (i,j)
      !$omp DO
      DO i=1,ngrid
#ifdef _OPENMP
         IF(verbose .AND. (modulo(i,100).EQ.0)) &
               WRITE(*,*) i,"/",ngrid," thread n. : ",omp_get_thread_num() 
#else
         IF(verbose .AND. (modulo(i,100).EQ.0)) &
               WRITE(*,*) i,"/",ngrid
#endif
         DO j=1,i-1
            ! distance between two voronoi centers
            ! also computes the nearest neighbor distance (squared) for each grid point
            distmm(i,j) = pammr2(D,period,y(:,i),y(:,j))
            if (distmm(i,j) < rgrid(i)) rgrid(i) = distmm(i,j)
            ! the symmetrical one
            distmm(j,i) = distmm(i,j)
            if (distmm(i,j) < rgrid(j)) rgrid(j) = distmm(i,j)   
         ENDDO
      ENDDO
      !$omp ENDDO
      !$omp END PARALLEL


      IF(verbose) WRITE(*,*) &
        "Local estimate of kernel widths"   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                    !!
!!!            Calculate local covariance matrix on grid               !!
!!!                 utilizing a two pass algorithm                     !!
!!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! If not specified, all localization weights are set to one
      IF (lfac.GT.0.0d0) prefac = -0.5d0/(lfac*lfac)

      DO i=1,ngrid
        ! localization
        DO ii=1,D
          ytmp(ii,:) = y(ii,:)-y(ii,i)
          IF (period(ii)<=0.0d0) CYCLE    
          ! minimum image convention  
          ytmp(ii,:) = ytmp(ii,:) / period(ii)
          ytmp(ii,:) = ytmp(ii,:) - DNINT(ytmp(ii,:))
          ytmp(ii,:) = ytmp(ii,:) * period(ii)
        ENDDO
        ! estimate weights for localization as product from 
        ! spherical gaussian weights and points in voronoi
        wQgrid = EXP(prefac*SUM(ytmp*ytmp,1))  
        ! estimate local number of sample points
        nlocal = SUM(npvoronoi*wQgrid)
        ! normalize weights
        wQgrid = wQgrid / SUM(wQgrid)
        ! add the voronoi weights
        wQgrid = wQgrid * npvoronoi/nsamples
        ! renormalize again
        wQgrid = wQgrid / SUM(wQgrid)
        
        ! estimate the mean
        DO ii=1,D
          IF (period(ii)>0.0d0) THEN
            dummd1 = SUM(SIN(y(ii,:))*wQgrid)
            dummd2 = SUM(COS(y(ii,:))*wQgrid)
            ytmp(ii,:) = y(ii,:) - ATAN2(dummd1,dummd2)
            ! minimum image convention
            ytmp(ii,:) = ytmp(ii,:) / period(ii)
            ytmp(ii,:) = ytmp(ii,:) - DNINT(ytmp(ii,:))
            ytmp(ii,:) = ytmp(ii,:) * period(ii)
          ELSE
            ytmp(ii,:) = y(ii,:) - SUM(y(ii,:)*wQgrid)
          ENDIF
          ytmpw(ii,:) = ytmp(ii,:)*wQgrid
        ENDDO
        
        ! estimate covariance matrix
        CALL DGEMM("N", "T", D, D, ngrid, 1.0d0, ytmp, D, ytmpw, D, 0.0d0, Q, D)
        Q = Q / (1.0d0-SUM(wQgrid**2.0d0))
        
        ! assign bandwidth using scotts rule
        Hi(:,:,i) = Q * nlocal**(-1.0d0/(D+4.0d0))
        ! inverse of the bandwidth matrix
        CALL invmatrix(D,Hi(:,:,i),Hiinv(:,:,i))
        
        ! estimate the normalization constants
        IF(periodic) THEN
          dummd1 = 1.0d0
          DO ii=1,D
             dummd1=dummd1*BESSI0(Hiinv(ii,ii,i))*twopi
          ENDDO
          normkernel(i) = 1.0d0/dummd1
        ELSE
          normkernel(i) = 1.0d0/DSQRT((twopi**DBLE(D))*detmatrix(D,Hi(:,:,i)))
        ENDIF
        
      ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                    !! 
!!!          Full covariance matrix estimation using DGEMM             !!
!!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!      ! If not specified, all localization weights are set to one
!      IF (lfac.GT.0.0d0) prefac = -0.5d0/(lfac*lfac) 

!      !$omp parallel &
!      !$omp default (none) &
!      !$omp shared (ngrid,nsamples,D,period,y,x,verbose,Hi,Hiinv,periodic,prefac,twopi) &
!      !$omp private (i,ii,Qlocal,nlocal,wQ,dummd1,dummd2,xtmp,xtmpw,normkernel)
!      !$omp DO
!      DO i=1,ngrid
!#ifdef _OPENMP
!         IF(verbose .AND. (modulo(i,100).EQ.0)) &
!               WRITE(*,*) i,"/",ngrid," thread n. : ",omp_get_thread_num() 
!#else
!         IF(verbose .AND. (modulo(i,100).EQ.0)) &
!               WRITE(*,*) i,"/",ngrid
!#endif  
!          
!        ! localization
!        DO ii=1,D
!          xtmp(ii,:) = x(ii,:)-y(ii,i)
!          IF (period(ii)<=0.0d0) CYCLE    
!          ! minimum image convention     
!          xtmp(ii,:) = xtmp(ii,:) / period(ii)
!          xtmp(ii,:) = xtmp(ii,:) - DNINT(xtmp(ii,:))
!          xtmp(ii,:) = xtmp(ii,:) * period(ii)
!        ENDDO
!        wQ = EXP(prefac*SUM(xtmp*xtmp,1))
!        nlocal = SUM(wQ) 
!        
!        ! estimate the mean
!        DO ii=1,D
!          IF (period(ii)>0.0d0) THEN
!            dummd1 = SUM(SIN(x(ii,:))*wQ)/nlocal
!            dummd2 = SUM(COS(x(ii,:))*wQ)/nlocal
!            xtmp(ii,:) = x(ii,:) - ATAN2(dummd1,dummd2)
!            ! minimum image convention
!            xtmp(ii,:) = xtmp(ii,:) / period(ii)
!            xtmp(ii,:) = xtmp(ii,:) - DNINT(xtmp(ii,:))
!            xtmp(ii,:) = xtmp(ii,:) * period(ii)
!          ELSE
!            xtmp(ii,:) = x(ii,:) - SUM(x(ii,:)*wQ)/nlocal
!          ENDIF
!          xtmpw(ii,:) = xtmp(ii,:)*wQ
!        ENDDO
!        !WRITE(*,*) 'local mean using sample points: ', SUM(x(ii,:)*wQ)/nlocal
!        
!        ! estimate covariance/variances
!        IF (periodic) THEN
!          Qlocal = 0.0d0
!          DO ii=1,D
!            CALL DGEMM("N", "T", 1, 1, nsamples, 1.0d0, xtmp(ii,:), 1, xtmpw(ii,:), 1, 0.0d0, Qlocal(ii,ii), 1) 
!          ENDDO
!        ELSE 
!          CALL DGEMM("N", "T", D, D, nsamples, 1.0d0, xtmp, D, xtmpw, D, 0.0d0, Qlocal, D)
!        ENDIF
!        Qlocal = Qlocal / (nlocal-1.0d0)
!        !WRITE(*,*) 'local cov from samples: ', Q

!        ! estimate local bandwidth using Scotts rule of thumb
!        Hi(:,:,i) = Qlocal * nlocal**(-1.0d0/(D+4.0d0))
!        
!        ! estimate the normalization constants
!        IF(periodic) THEN
!          dummd1 = 1.0d0
!          DO ii=1,D
!             dummd1=dummd1*BESSI0(Hiinv(ii,ii,i))*twopi
!          ENDDO
!          normkernel(i) = 1.0d0/dummd1
!        ELSE
!          normkernel(i) = 1.0d0/DSQRT((twopi**DBLE(D))*detmatrix(D,Hi(:,:,i)))
!        ENDIF
!      ENDDO
!      !$omp ENDDO
!      !$omp END PARALLEL 
!      
      ! set the lambda to be used in QS
      IF(lambda.LT.0)THEN
        ! compute the median of the NN distances
        lambda=4.0d0*median(ngrid,DSQRT(rgrid(:)))
        lambda2=lambda*lambda
      ENDIF
      
      IF(verbose) WRITE(*,*) &
        "Computing kernel density on reference points"
      
      prob = 0.0d0
      bigp = 0.0d0
      
      DO i=1,ngrid
        DO j=1,ngrid
          
          dummd1 = 0.0d0
          IF(periodic .and. (nbootstrap.eq.0)) tmpcheck=0.0d0
          
          IF (mahalanobis(D,period,y(:,i),y(:,j),Hiinv(:,:,j))>36.0d0) THEN
            ! assume distribution in far away grid point is narrow
            IF(periodic)THEN
              dummd2=1.0d0
              CALL pammrij(D, period, y(:,i), y(:,j), xij)
              DO jj=1,D
                 dummd2 = dummd2 * fkernelvm(Hiinv(jj,jj,j),xij(jj))
              ENDDO
              ! TODO: change npvoronoi to the sum of weights
              tmpkernel = npvoronoi(j) * dummd2
              IF(nbootstrap.eq.0) &
                tmpcheck = fmultikernel(D,period,y(:,i),y(:,j),Hiinv(:,:,j)) &
                         * npvoronoi(j) &
                         + tmpcheck 
            ELSE
              tmpkernel = fmultikernel(D,period,y(:,i),y(:,j),Hiinv(:,:,j))*npvoronoi(j)
            ENDIF     
            dummd1 = tmpkernel
          ELSE
            ! cycle just inside the polyhedra using the neighbour list
            DO k=pnlist(j)+1,pnlist(j+1)
              IF(periodic)THEN
                dummd2=1.0d0
                CALL pammrij(D, period, y(:,i),x(:,nlist(k)), xij)
                DO jj=1,D
                   dummd2 = dummd2 * fkernelvm(Hiinv(jj,jj,j),xij(jj))
                ENDDO
                tmpkernel = wj(nlist(k)) * dummd2
                IF(nbootstrap.eq.0) &
                  tmpcheck = fmultikernel(D,period,y(:,i),x(:,nlist(k)),Hiinv(:,:,j)) &
                           * wj(nlist(k)) &
                           + tmpcheck 
              ELSE
                tmpkernel = wj(nlist(k))*fmultikernel(D,period,y(:,i),x(:,nlist(k)),Hiinv(:,:,j))
              ENDIF     
              dummd1 = dummd1 + tmpkernel    
            ENDDO
          ENDIF
            
          ! get the non-normalized prob, used to compute the estimate of the error
          IF(periodic .and. (nbootstrap.eq.0)) THEN
            ! If we are using VM kernels and we are not estimating the error using BS
            ! then estimate the error using a gaussian kernel
            bigp(i) = bigp(i) + tmpcheck
          ELSE
            bigp(i) = bigp(i) + dummd1
          ENDIF
          ! We now have to normalize the kernel
          prob(i) = prob(i) + dummd1*normkernel(j)
          
        ENDDO          
      ENDDO
      prob=prob/normwj
      bigp=bigp/normwj
   
      ! use bootstrapping to estimate the error of our KDE

      IF(nbootstrap > 0) probboot = 0.0d0
  
      ! START of bootstrapping run
      DO nn=1,nbootstrap
        IF(verbose) WRITE(*,*) &
              "Bootstrapping, run ", nn

        ! rather than selecting nsel random points, we select a random 
        ! number of points from each voronoi. this makes it possible 
        ! to apply some simplifications and avoid computing distances 
        ! from far-away voronoi
        nbstot=0
        DO j=1,ngrid
          nbssample=random_binomial(nsamples, DBLE(npvoronoi(j))/DBLE(nsamples))
          nbstot = nbstot+nbssample
                 
          DO i=1,ngrid
            ! localization: do not compute KDEs for points that belong to far away Voronoi
            ! one could make this a bit more sophisticated, but this should be enough
            IF (mahalanobis(D,period,y(:,i),y(:,j),Hiinv(:,:,j))>36.0d0) CYCLE  
            
            tmpkernel=0.0d0
            DO k=1,nbssample
              rndidx = int(npvoronoi(j)*random_uniform())+1
              rndidx = nlist(pnlist(j)+rndidx)
              ! the vector that contains the Voronoi assignations is iminij   
              ! do not compute KDEs for points that belong to far away Voronoi
              IF(periodic)THEN
                dummd2=1.0d0
                CALL pammrij(D, period, y(:,i),x(:,rndidx), xij)
                DO jj=1,D
                   dummd2 = dummd2 & 
                          * fkernelvm(Hiinv(jj,jj,iminij(rndidx)),xij(jj))
                ENDDO
                tmpkernel = tmpkernel + dummd2 
              ELSE
                tmpkernel = tmpkernel &
                          + fmultikernel(D,period,y(:,i),x(:,rndidx), & 
                                         Hiinv(:,:,iminij(rndidx)))
              ENDIF
            ENDDO
            ! we have to normalize it
            probboot(i,nn) = probboot(i,nn) + tmpkernel*normkernel(j)
          ENDDO

        ENDDO
        probboot(:,nn) = probboot(:,nn)/nbstot
        
      ENDDO 

      ! END of bootstrapping run
      
      ! ---
      ! Computing the error
      ! ---
      
      ! get the error from the bootstrap if bootstrap was set
      IF(nbootstrap > 0) THEN
        proberr = 0.0d0
        DO i=1,ngrid
          ! this is the SD of the probability density (relative error)
          proberr(i) = DSQRT( SUM( (probboot(i,:) - prob(i))**2.0d0 ) / (nbootstrap-1.0d0) ) / prob(i)
        ENDDO 
      ELSE
        DO i=1,ngrid  
          proberr(i) = ( ( (1.0d0/bigp(i)) - 1.0d0 ) * (1.0d0/normwj) )**(0.5d0) 
        ENDDO
      ENDIF

      ! set the lambda to be used in the old version of QS 
      DO i=1,ngrid
        sigma2(i) = trmatrix(D,Hi(:,:,i))/D
      ENDDO
      IF(lambda.LT.0)THEN
         lambda=5.0d0 * SUM( DSQRT(sigma2) )/ngrid
         lambda2=lambda*lambda
      ENDIF
      
      idxroot=0
      ! Start quick shift

      DO i=1,ngrid
         IF(idxroot(i).NE.0) CYCLE
         qspath=0
         qspath(1)=i
         counter=1         
         DO WHILE(qspath(counter).NE.idxroot(qspath(counter)))
            idxroot(qspath(counter))= &
           ! qs_next(ngrid,qspath(counter),prob,distmm,rgrid,kderr,qserr,lambda2, & 
           !         verbose,neblike,nbootstrap,proberr)
            qs_next(ngrid,qspath(counter),lambda2,prob,distmm)
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
      
      IF(verbose) write(*,*) "Writing out"
      qspath=0
      qspath(1)=idxroot(1)
      Nk=1
      normpks=0.0d0
      OPEN(UNIT=11,FILE=trim(outputfile)//".grid",STATUS='REPLACE',ACTION='WRITE')
      DO i=1,ngrid
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
         idcls(i)=dummyi1 ! stores the cluster index
         DO j=1,D
           WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", y(j,i)
         ENDDO
         !print out the squared absolute error
         WRITE(11,"(A1,I4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4)") " ", dummyi1 , &
                                                        " " , prob(i), " " , &
                                                   proberr(i), " ",  trmatrix(D,Hi(:,:,i))/DBLE(D) 
         ! accumulate the normalization factor for the pks
         normpks=normpks+prob(i)
      ENDDO
      CLOSE(UNIT=11)
      
      ! builds the cluster adjacency matrix
      IF (verbose) WRITE(*,*) "Building cluster adjacency matrix"
      ALLOCATE(clsadj(Nk, Nk))
      ALLOCATE(macrocl(Nk))
      clsadj = 0.0d0
      DO i=1, Nk
         ! initialize each cluster to itself in the macrocluster assignation 
         macrocl(i)=i
         DO j=1,i-1
             clsadj(i,j) = cls_link(ngrid, idcls, distmm, prob, rgrid, i, j, proberr)
             clsadj(j,i) = clsadj(i,j)
         ENDDO
      ENDDO
      IF(saveadj)THEN
         OPEN(UNIT=11,FILE=trim(outputfile)//".adj",STATUS='REPLACE',ACTION='WRITE')
         DO i=1, Nk
             DO j=1, Nk
                 WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", clsadj(i,j)
             ENDDO
             WRITE(11,*) ""
         ENDDO
         CLOSE(11)
      ENDIF
      
      ! Let's print out the macroclusters
      DO i=1, Nk
         DO j=1, Nk
            ! Put a threshold under which there is no link between the clusters
            ! now it is just a default
            IF(i.EQ.j) CYCLE ! discard yourself
            IF(clsadj(i,j) .GT. thrmerg) THEN 
               IF(macrocl(j).EQ.j) THEN
                  ! the point is still initialized to himself 
                  macrocl(j)=macrocl(i)
               ELSE
                  ! it was already assigned
                  ! lets change also all the values that I may have changed before
                  DO k=1,j-1
                    IF(k.EQ.i) CYCLE ! I'll fix it later 
                    IF(macrocl(k).EQ.macrocl(i)) macrocl(k)=macrocl(j)
                  ENDDO
                  macrocl(i)=macrocl(j)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
      ! Count unique macroclusters and order them
      ALLOCATE(sortmacrocl(Nk))
      sortmacrocl=0
      dummyi1=0
      DO i=1, Nk
        isthere=.false.
        DO j=1, Nk
           IF( (.NOT.(sortmacrocl(j).EQ.0)) .AND. (macrocl(i).EQ.j)) THEN
              ! position j has already been set to something 
              ! and the value at the jth position corrispond to my cluster idx
              isthere=.true.
              macrocl(i)=sortmacrocl(macrocl(i))
              EXIT
           ENDIF 
        ENDDO
        IF(.NOT. isthere) THEN
           dummyi1=dummyi1+1
           ! increase the number of macroclusters found
           sortmacrocl(macrocl(i))=dummyi1
           ! rewrite the macrocluster assignation with a proper index
           macrocl(i)=sortmacrocl(macrocl(i))
        ENDIF
      ENDDO
      
      IF (verbose) WRITE(6,"((A6,I7,A15))") " Found ",dummyi1," macroclusters."
      IF (saveadj) THEN
        OPEN(UNIT=11,FILE=trim(outputfile)//".macrogrid",STATUS='REPLACE',ACTION='WRITE')
        DO i=1,ngrid
           DO j=1,D
             WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", y(j,i)
           ENDDO
           WRITE(11,"(A1,I4,A1,I4)") " ", idcls(i) , " ", macrocl(idcls(i))
        ENDDO
        CLOSE(UNIT=11)
      ENDIF
      
      ! now we can procede and complete the definition of probability model
      ! now qspath contains the indexes of Nk gaussians
      IF(periodic) THEN
         ALLOCATE(vmclusters(Nk))
      ELSE
         ALLOCATE(clusters(Nk))
      ENDIF
      
      DO k=1,Nk
         IF(periodic)THEN
            ALLOCATE(vmclusters(k)%mean(D))
            ALLOCATE(vmclusters(k)%cov(D,D))
            ALLOCATE(vmclusters(k)%icov(D,D))
            ALLOCATE(vmclusters(k)%period(D))
            vmclusters(k)%period=period
            vmclusters(k)%mean=y(:,qspath(k))
         ELSE
            ALLOCATE(clusters(k)%mean(D))
            ALLOCATE(clusters(k)%cov(D,D))
            ALLOCATE(clusters(k)%icov(D,D))
            clusters(k)%mean=y(:,qspath(k))
         ENDIF
         ! optionally do a few mean-shift steps to find a better estimate 
         ! of the cluster mode
         DO j=1,nmsopt
            msmu=0.0d0
            tmppks=0.0d0
            
            DO i=1,ngrid
               ! should correct the Gaussian evaluation with a Von Mises distrib in the case of periodic data
               IF(periodic)THEN
                  msw = prob(i)*exp(-0.5*pammr2(D,period,y(:,i),vmclusters(k)%mean)/(lambda2/16.0d0))
                  CALL pammrij(D,period,y(:,i),vmclusters(k)%mean,tmpmsmu)
               ELSE
                  msw = prob(i)*exp(-0.5*pammr2(D,period,y(:,i),clusters(k)%mean)/(lambda2/16.0d0))
                  CALL pammrij(D,period,y(:,i),clusters(k)%mean,tmpmsmu)
               ENDIF
               
               msmu = msmu + msw*tmpmsmu
               tmppks = tmppks + msw
            ENDDO
            
            IF(periodic)THEN
               vmclusters(k)%mean = vmclusters(k)%mean + msmu / tmppks
            ELSE
               clusters(k)%mean = clusters(k)%mean + msmu / tmppks
            ENDIF
         ENDDO
         
         ! compute the gaussians covariance from the data in the clusters
         IF(periodic)THEN
            vmclusters(k)%cov = 0.0d0
         ELSE
            clusters(k)%cov = 0.0d0
         ENDIF
         
         tmppks=0.0d0
         DO i=1,ngrid
            IF(idxroot(i).NE.qspath(k)) CYCLE
            tmppks=tmppks+prob(i)
            DO ii=1,D
               DO jj=1,D
                  IF(periodic)THEN
!                     ! Minimum Image Convention
                     dummd1=DSQRT(pammr2(1,period(ii),y(ii,i),vmclusters(k)%mean(ii)))
                     dummd2=DSQRT(pammr2(1,period(jj),y(jj,i),vmclusters(k)%mean(jj)))
                     
                     vmclusters(k)%cov(ii,jj)= vmclusters(k)%cov(ii,jj)+prob(i)* &
                                               dummd1 * dummd2
                                            
                  ELSE
                     dummd1=DSQRT(pammr2(1,period(ii),y(ii,i),clusters(k)%mean(ii)))
                     dummd2=DSQRT(pammr2(1,period(jj),y(jj,i),clusters(k)%mean(jj)))
                     
                     clusters(k)%cov(ii,jj)=clusters(k)%cov(ii,jj)+prob(i)* &
                                               dummd1 * dummd2
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         
         IF(periodic)THEN
            vmclusters(k)%cov=vmclusters(k)%cov/tmppks
            vmclusters(k)%weight=tmppks/normpks
            vmclusters(k)%D=D
         ELSE
            clusters(k)%cov=clusters(k)%cov/tmppks
            clusters(k)%weight=tmppks/normpks
            clusters(k)%D=D
         ENDIF
      ENDDO

      IF(periodic)THEN
         ! write the VM distributions
         ! write a 2-lines header containig a bit of information
         WRITE(comment,*) "# PAMMv2 clusters analysis. NSamples: ", nsamples, " NGrid: ", &
                   ngrid, " QSLambda: ", lambda, ACHAR(10), & 
                   "# Dimensionality/NClusters//Pk/Mean/Covariance/Period"

         OPEN(UNIT=12,FILE=trim(outputfile)//".pamm",STATUS='REPLACE',ACTION='WRITE')
         CALL writevmclusters(12, comment, nk, vmclusters)
         CLOSE(UNIT=12)
         DEALLOCATE(vmclusters)
      ELSE
         ! write the Gaussians       
         ! write a 2-lines header
         WRITE(comment,*) "# PAMMv2 clusters analysis. NSamples: ", nsamples, " NGrid: ", &
                   ngrid, " QSLambda: ", lambda, ACHAR(10), "# Dimensionality/NClusters//Pk/Mean/Covariance"
         
         OPEN(UNIT=12,FILE=trim(outputfile)//".pamm",STATUS='REPLACE',ACTION='WRITE')
         
         CALL writeclusters(12, comment, nk, clusters)
         CLOSE(UNIT=12)
         ! maybe I should deallocate better..
         DEALLOCATE(clusters)
      ENDIF
      
      DEALLOCATE(x,wj)
      DEALLOCATE(period)
      DEALLOCATE(idxroot,qspath,distmm)
      DEALLOCATE(pnlist,nlist,iminij,bigp)
      DEALLOCATE(y,npvoronoi,prob,sigma2,rgrid,oldsigma2)
      DEALLOCATE(diff,msmu,tmpmsmu)
      DEALLOCATE(Q,Qlocal,Hi,Hiinv,normkernel)
      DEALLOCATE(xtmp,xtmpw,xij,ytmp,ytmpw)
      DEALLOCATE(wQ,Sw,Sinv,wQgrid)
      DEALLOCATE(macrocl,sortmacrocl)
      IF(nbootstrap>0) DEALLOCATE(probboot,proberr)

      CALL EXIT(0)
      ! end of the main programs



!!!!! FUCTIONS and SUBROUTINES !!!!!!!!!!!!!!!!!!!!

      CONTAINS

      SUBROUTINE helpmessage
         ! Banner to print out for helping purpose
         !

         WRITE(*,*) ""
         WRITE(*,*) " USAGE: pamm [-h] -d D [-p 6.28,6.28,...] [-w] [-o output] [-ngrid ngrid] "
         WRITE(*,*) "             [-l lambda] [-kde err] [-z zeta_factor] [-a smoothing_factor] "
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
         WRITE(*,*) "   -d D              : Dimensionality "
         WRITE(*,*) "   -w                : Reads weights for the sample points [default: no weight] "
         WRITE(*,*) "   -o output         : Prefix for output files [out]. This will produce : "
         WRITE(*,*) "                         output.grid (clusterized grid points) "
         WRITE(*,*) "                         output.pamm (cluster parameters) "
         WRITE(*,*) "   -l lambda         : Quick shift cutoff [automatic] "
         WRITE(*,*) "                       (Not used with the -neblike flag active)"
         WRITE(*,*) "   -targeterror      : Relative target error use to adaptively refine the KDE [default: 0.1]"
         WRITE(*,*) "   -qserr err        : Relative threshold used during quickshift [8]"
         WRITE(*,*) "   -ngrid ngrid      : Number of grid points to evaluate KDE [sqrt(nsamples)]"
         WRITE(*,*) "   -bootstrap N      : Number of iteretions to do when using bootstrapping "
         WRITE(*,*) "                       to refine the KDE on the grid points"
         WRITE(*,*) "   -neblike N        : Try to improve the clustering using neblike path search algorithm [0] "
         WRITE(*,*) "                       N is the number of iterations "
         WRITE(*,*) "   -nms nms          : Do nms mean-shift steps with a Gaussian width lambda/5 to "
         WRITE(*,*) "                       optimize cluster centers [0] "
         WRITE(*,*) "   -seed seed        : Seed to initialize the random number generator. [12345]"
         WRITE(*,*) "   -p P1,...,PD      : Periodicity in each dimension [ (6.28,6.28,6.28,...) ]"
         WRITE(*,*) "   -savevoronois     : Save Voronoi associations. This will produce:"
         WRITE(*,*) "                         output.voronoislinks (points + associated Voronoi) "
         WRITE(*,*) "                         output.voronois (Voronoi centers + info) "
         WRITE(*,*) "   -saveprobs        : Save Voronoi associations. This will produce:"
         WRITE(*,*) "                         output.voronoislinks (points + associated Voronoi) "
         WRITE(*,*) "                         output.voronois (Voronoi centers + info) "
         WRITE(*,*) "   -l sigma          : Localization width for local bayesian run [default: off] "
         WRITE(*,*) "   -adj threshold    : Set the threshold to merge adjcent clusters and "
         WRITE(*,*) "                       write out the adjacency matrix "
         WRITE(*,*) "   -v                : Verbose output "
         WRITE(*,*) ""
         WRITE(*,*) " Post-processing mode (-gf): this reads high-dim data and computes the "
         WRITE(*,*) " cluster probabilities associated with them, given the output of a "
         WRITE(*,*) " previous PAMM analysis. "
         WRITE(*,*) ""
         WRITE(*,*) "   -gf               : File to read reference Gaussian clusters from"
         WRITE(*,*) "   -a                : Additional smearing of clusters "
         WRITE(*,*) "   -z zeta_factor : Probabilities below this threshold are counted as 'no cluster' [default:0]"
         WRITE(*,*) ""
      END SUBROUTINE helpmessage
      
      DOUBLE PRECISION FUNCTION median(ngrid,a)
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, intent(in) :: a(ngrid)
         
         INTEGER :: l
         DOUBLE PRECISION, dimension(size(a,1)) :: ac
         
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
      
      INTEGER FUNCTION  findMinimum(x, startidx, endidx )
         INTEGER, INTENT(IN) :: startidx, endidx
         DOUBLE PRECISION, DIMENSION(1:), INTENT(IN) :: x
         
         DOUBLE PRECISION :: minimum
         INTEGER :: location
         INTEGER :: i
   
         minimum  = x(startidx)		! assume the first is the min
         location = startidx			! record its position
         DO i = startidx+1, endidx		! start with next elements
            IF (x(i) < minimum) THEN	!   if x(i) less than the min?
               minimum  = x(i)		!      Yes, a new minimum found
               location = i                !      record its position
            END IF
         END DO
         findMinimum = Location        	! return the position
      END FUNCTION  findMinimum
      
      SUBROUTINE swap(aa, bb)
      !  This subroutine swaps the values of its two formal arguments.
         DOUBLE PRECISION, INTENT(INOUT) :: aa, bb
         DOUBLE PRECISION                :: temp
         temp = aa
         aa = bb
         bb = temp
      END SUBROUTINE swap
      
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
         totw=0.0d0
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

      SUBROUTINE mkgrid(D,period,nsamples,ngrid,x,y,npvoronoi,iminij,wj)
         ! Select ngrid grid points from nsamples using minmax and
         ! the voronoi polyhedra around them.
         !
         ! Args:
         !    nsamples: total points number
         !    ngrid: number of grid points
         !    x: array containing the data samples
         !    y: array that will contain the grid points
         !    npvoronoi: array cotaing the number of samples inside the Voronoj polyhedron of each grid point
         !    iminij: array containg the neighbor list for data samples

         INTEGER, INTENT(IN) :: D
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         INTEGER, INTENT(IN) :: nsamples
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, DIMENSION(D,nsamples), INTENT(IN) :: x
         DOUBLE PRECISION, DIMENSION(D,ngrid), INTENT(OUT) :: y
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: npvoronoi
         INTEGER, DIMENSION(nsamples), INTENT(OUT) :: iminij
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: wj      

         INTEGER i,j
         DOUBLE PRECISION :: dminij(nsamples), dij, dmax

         iminij=0
         y=0.0d0
         npvoronoi=0
         ! choose randomly the first point
         y(:,1)=x(:,int(RAND()*nsamples))
         dminij = 1.0d99
         iminij = 1         
         DO i=2,ngrid
            dmax = 0.0d0
            
            !$omp parallel &
            !$omp default (none) &
            !$omp shared (jmax,dmax,i,D,period,y,x,nsamples,dminij,iminij) &
            !$omp private (j,dij)
            
            !$omp DO
            
            DO j=1,nsamples
               dij = pammr2(D,period,y(:,i-1),x(:,j))
               IF (dminij(j)>dij) THEN
                  dminij(j) = dij
                  iminij(j) = i-1 ! also keeps track of the Voronoi attribution
               ENDIF
               IF (dminij(j) > dmax) THEN
                  !$omp CRITICAL
                  
                  dmax = dminij(j)
                  jmax = j
                  
                  !$omp END CRITICAL
               ENDIF
            ENDDO   
            
            !$omp ENDDO
            !$omp END PARALLEL         
            y(:,i) = x(:, jmax)
            IF(verbose .AND. (modulo(i,1000).EQ.0)) &
               write(*,*) i,"/",ngrid
         ENDDO

         ! finishes Voronoi attribution
         DO j=1,nsamples
            dij = pammr2(D,period,y(:,ngrid),x(:,j))
            IF (dminij(j)>dij) THEN
               dminij(j) = dij
               iminij(j) = ngrid
            ENDIF
         ENDDO

         ! Assign neighbor list pointer of voronois
         ! Number of points in each voronoi polyhedra
         npvoronoi=0
         DO j=1,nsamples
            npvoronoi(iminij(j))=npvoronoi(iminij(j))+1
         ENDDO
      END SUBROUTINE mkgrid

      SUBROUTINE getnlist(nsamples,ngrid,npvoronoi,iminij, pnlist,nlist)
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
         INTEGER, DIMENSION(ngrid), INTENT(IN) :: npvoronoi
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
            pnlist(i+1)=pnlist(i)+npvoronoi(i)
            tmpnidx(i)=pnlist(i)+1  ! temporary array to use while filling up the neighbour list
         ENDDO

         DO j=1,nsamples
            i=iminij(j) ! this is the Voronoi center the sample j belongs to
            nlist(tmpnidx(i))=j ! adds j to the neighbour list
            tmpnidx(i)=tmpnidx(i)+1 ! advances the pointer
         ENDDO
      END SUBROUTINE getnlist

      DOUBLE PRECISION FUNCTION cls_link(ngrid, idcls, distmm, prob, rgrid, ia, ib, errors)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngrid, idcls(ngrid), ia, ib
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid, ngrid), prob(ngrid), rgrid(ngrid)
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: errors
         
         INTEGER i, j
         DOUBLE PRECISION mxa, mxb, mxab, pab, emxa, emxb, emxab, tmpmin
         mxa = 0.0d0
         mxb = 0.0d0
         mxab = 0.0d0
         emxa = 0.0d0
         emxb = 0.0d0
         emxab = 0.0d0
         DO i=1, ngrid
            IF (idcls(i)/=ia) CYCLE
            IF (prob(i).gt.mxa) THEN
               mxa = prob(i)
               emxa = errors(i)*prob(i) ! also gets the probability density at the mode of cluster a
            ENDIF
            DO j=1,ngrid
               IF (idcls(j)/=ib) CYCLE
               IF (prob(j).gt.mxb) THEN
                  mxb = prob(j)
                  emxb = errors(j)*prob(j)
               ENDIF
               ! Ok, we've got a matching pair
               IF (dsqrt(distmm(i,j))<dsqrt(rgrid(i))+dsqrt(rgrid(j))) THEN
                  ! And they are close together!
                  pab = (prob(i)+prob(j))/2
                  IF (pab .gt. mxab) THEN
                     mxab = pab
                     emxab = DSQRT((errors(i)*prob(i))**2 &
                               +(errors(j)*prob(j))**2)
                  ENDIF
               ENDIF               
            ENDDO            
         ENDDO
         IF(mxa.EQ.0)THEN
            cls_link=0.0d0
         ELSE
            cls_link = (mxab+emxab)/min(max(mxa-emxa,0.0d0),max(mxb-emxb,0.0d0))
         ENDIF
      END FUNCTION
      
      SUBROUTINE build_path(ngrid, path, i0, i1, npath, distmm, rgrid) 
         ! builds a "straight" path across the grid, between points ia and ib
         ! by iteratively refining down to the resolution of the grid
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngrid, i0, i1
         INTEGER, INTENT(OUT) :: path(ngrid), npath
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid,ngrid), rgrid(ngrid)
         
         
         INTEGER l0, ia, ib, ic, i
         DOUBLE PRECISION :: dab, dthresh, dac, dbc, dmin, dcomb
         LOGICAL, DIMENSION(ngrid) :: pactive      
         
         pactive = .TRUE. ! these points can be added to the path         
         path(1) = i0
         path(2) = i1
         pactive(i0) = .FALSE.
         pactive(i1) = .FALSE.
         npath = 2
         
         ! path segments are marked as "inactive" by making the index negative
         DO WHILE (path(npath-1)>0) 
            l0 = 1            
            DO WHILE(path(l0)<0) ! select the first "active" segment
               l0 = l0+1
            ENDDO
            
            ! refines the segment. we look for a point that minimizes (d(a,c)+d(b,c))/d(a,b)-1 + abs(d(a,c)-d(b,c))/d(a,b)
            ! abs(d(a,c)-d(b,c))/d(a,b) is a measure of "centrality" of the point
            ! d(a,c)+d(b,c))/d(a,b)-1 is a measure of how aligned is each point with its neighbors
                           
            dmin=1d100
            ia = path(l0)
            ib = path(l0+1)
            dab = dsqrt(distmm(ia,ib))
            dthresh = dsqrt(rgrid(ia))+dsqrt(rgrid(ib)) ! do not refine below grid granularity
            ic=ia
                          
            DO i=1,ngrid               
               IF (pactive(i).eqv..FALSE.) CYCLE
               dac=dsqrt(distmm(ia,i))
               IF (dac<dthresh .or. dac>dab) CYCLE               ! discard points below the target granularity or too far away
               dbc=dsqrt(distmm(ib,i))
               IF (dbc<dthresh .or. dbc>dab) CYCLE              
               dcomb = (dac+dbc-dab)+dabs(dac-dbc)  ! find a point that is close to the midpoint of ia-ib
               IF (dcomb < dmin ) THEN
                  dmin=dcomb
                  ic = i
               ENDIF
            ENDDO
            
            ! found a new point to refine the path?
            IF (ic /= ia .and. ic /= ib) THEN
            !   IF ((prob(ic)-prob(idx)/prob(idx))<0) &
            !      write(*,*) "going down", prob(ic), prob(idx)                  
               ! add the new point to the path
               pactive(ic) = .FALSE.
               DO i=npath,l0+1,-1  ! shifts the path to make space for the new point
                  path(i+1)=path(i)
               ENDDO
               path(l0+1)=ic
               npath=npath+1            
            ELSE  ! mark the start point as inactive
               path(l0) = -path(l0)
            ENDIF
         ENDDO
         path = abs(path)
      END SUBROUTINE 
      
      SUBROUTINE neb_path(ngrid, path, npath, distmm, prob, w) 
         ! "refines" path to look a bit like a maximum probability path going through saddle points
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngrid, npath
         INTEGER, INTENT(INOUT) :: path(ngrid)
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid,ngrid), prob(ngrid), w
         
         
         INTEGER i, j, npi
         DOUBLE PRECISION :: dab, dac, dbc, dmin, dcomb, wi
         LOGICAL, DIMENSION(ngrid) :: pactive
         
         pactive = .TRUE.  ! new points can't be one of those already assigned to the path
         DO i=1,npath
            pactive(path(i)) = .FALSE.
         ENDDO
         DO i=2, npath-1           
           ! "climbing image" mode to try and reach the lowest probability point along the path
           wi = -w
           dmin = prob(1)
           npi = 1
           DO j=2, npath
               IF (prob(j)<dmin) THEN
                  dmin=prob(j)
                  npi = j
               ENDIF
           ENDDO
           IF (npi==i) wi = w  
            
           dmin = 1d100
           npi = path(i)
           dab = distmm(path(i-1), path(i+1))                      
           DO j=1,ngrid
              IF (.not. pactive(j)) CYCLE
              dac = distmm(j, path(i-1))
              IF (dac>dab) CYCLE
              dbc = distmm(j, path(i+1))
              IF (dbc>dab) CYCLE
              IF (wi<0) THEN
                 dcomb = (dac+dbc)/dab-1 + abs(dac-dbc)/dab + & 
                       wi * prob(j)/prob(path(i))
              ELSE
                  dcomb = (dac+dbc)/dab-1 + & 
                       wi * prob(j)/prob(path(i))
              ENDIF
              IF (dcomb<dmin) THEN ! found a better point, update the path
                 pactive(npi) = .TRUE.
                 pactive(j) = .FALSE.
                 dmin = dcomb
                 npi =j
              ENDIF
           ENDDO
           path(i) = npi
         ENDDO
      END SUBROUTINE
      
!      INTEGER FUNCTION qs_next(ngrid,idx,prob,distmm,rgrid,kderr,qserr,lambda2, &
!                               verbose,neblike,nbootstrap,errors)
!         ! Return the index of the closest point higher in P
!         ! 
!         ! Args:
!         !    ngrid:   number of grid points
!         !    idx:     current point
!         !    prob:    density estimations
!         !    distmm:  distances matrix
!         !    rgrid:  
!         !    kderr:
!         !    qserr:
!         !    lambda2: cut-off in the jump
!         !    verbose:
!         !    neblike
!         !    nbootstrap
!         !    errors
!
!
!         IMPLICIT NONE
!         INTEGER, INTENT(IN) :: ngrid
!         DOUBLE PRECISION, INTENT(IN), DIMENSION(ngrid) :: rgrid
!         INTEGER, INTENT(IN) :: idx,nbootstrap
!         DOUBLE PRECISION, INTENT(IN) :: kderr,qserr,lambda2
!         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: prob
!         DOUBLE PRECISION, DIMENSION(ngrid,ngrid), INTENT(IN) :: distmm
!         LOGICAL, INTENT(IN) :: verbose
!         INTEGER, INTENT(IN) :: neblike
!         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: errors
!
!         INTEGER i, j, ndx
!         
!         INTEGER :: npath
!         INTEGER :: path(ngrid)
!         LOGICAL :: fsaddle, fjump
!         DOUBLE PRECISION :: dmin
!
!
!         dmin=1.0d100
!         ! find the nearest point with higher probability         
!         ndx=idx
!         DO j=1,ngrid       
!            IF( (prob(j) .GT. prob(idx)) .AND. & 
!                ((prob(j)-errors(j)) .GT. (prob(idx)+errors(idx))) )THEN 
!               IF(distmm(idx,j) .LT. dmin)THEN
!                   ndx=j
!                   dmin=distmm(idx,j)
!               ENDIF
!            ENDIF
!         ENDDO
!           
!           
!           
!         qs_next=idx  
!         IF( (distmm(idx,ndx).LT.(lambda2)) ) THEN 
!            
!           qs_next = ndx   
!           IF (ndx/=idx) THEN         
!              npath = 0
!              ! builds a discrete near-linear path between the end points
!    
!    ! TO BE CHECKED. Somehow there could be a SIGSEV if the array is not initialized
!    ! to an initial value. There is an error somewhere that has to be debugged
!     path=idx
!     
!              CALL build_path(ngrid, path, idx, ndx, npath, distmm, rgrid) 
!              ! refine the path aiming for a "maximum probability path" 
!              IF (npath > 2) THEN  ! does a few iterations for path refinement
!                 DO i=1,3!neblike
!                    CALL neb_path(ngrid, path, npath, distmm, prob, 0.1d0)
!                 ENDDO
!              ENDIF            
!              
!              fsaddle = .false.
!              fjump = .false.
!              ! check if the path goes downhill to within accuracy
!              DO i=2,npath-1   
!                    IF ( ((prob(path(i))-prob(idx))/ &
!                       (prob(path(i))+prob(idx))<-0.03d0) )THEN 
!                    !    IF(errors(idx)/errors(path(i))<0.3d0) CYCLE 
!                        qs_next = idx
!                        fsaddle = .true.
!                        EXIT
!                    ENDIF
!              ENDDO
!              ! check if the path contains crazy jumps
!              DO i=1,npath-1
!                 IF (dsqrt( distmm(path(i),path(i+1)) ) >  &
!                    6*(dsqrt(rgrid(path(i)))+dsqrt(rgrid(path(i+1)))) ) THEN
!                    qs_next=idx ! abort jump!
!                    fjump = .true.
!                    EXIT
!                 ENDIF
!              ENDDO  
!              IF (verbose .and. (fsaddle .or. fjump)) THEN
!                 dmin = 0.0d0
!                 IF (fsaddle) write(*,*) "# SADDLE POINT DETECTED"
!                 IF (fjump) write(*,*) "# SADDLE POINT DETECTED"            
!                 write(*,*) 1, dmin, prob(path(1)), path(1), y(:,path(1))
!                 DO j=2,npath                        
!                    dmin = dmin + distmm(path(j),path(j-1))
!                    write(*,*) j, dmin, prob(path(j)), path(j), y(:,path(j))
!                 ENDDO
!              ENDIF   
!           ENDIF   
!         ENDIF
!         
!                     
!      END FUNCTION qs_next

      INTEGER FUNCTION qs_next(ngrid,idx,lambda2,probnmm,distmm)
         ! Return the index of the closest point higher in P
         !
         ! Args:
         !    ngrid: number of grid points
         !    idx: current point
         !    lambda: cut-off in the jump
         !    probnmm: density estimations
         !    distmm: distances matrix

         INTEGER, INTENT(IN) :: ngrid
         INTEGER, INTENT(IN) :: idx
         DOUBLE PRECISION, INTENT(IN) :: lambda2
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: probnmm
         DOUBLE PRECISION, DIMENSION(ngrid,ngrid), INTENT(IN) :: distmm

         INTEGER j
         DOUBLE PRECISION dmin

         dmin=1.0d10
         qs_next=idx
         DO j=1,ngrid
            IF(probnmm(j).GT.probnmm(idx))THEN
               IF((distmm(idx,j).LT.dmin) .AND. (distmm(idx,j).LT.lambda2))THEN
                  dmin=distmm(idx,j)
                  qs_next=j
               ENDIF
            ENDIF
         ENDDO
      END FUNCTION qs_next
      
      DOUBLE PRECISION FUNCTION fmultikernel(D,period,x,y,icov)
         ! Return the multivariate gaussian density
         ! Args:
         !    gpars: gaussian parameters
         !    x: point in wich estimate the value of the gaussian
         
         INTEGER , INTENT(IN) :: D
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: x(D)
         DOUBLE PRECISION, INTENT(IN) :: y(D)
         DOUBLE PRECISION, INTENT(IN) :: icov(D,D)
         DOUBLE PRECISION dv(D),tmpv(D),xcx
         
         CALL pammrij(D, period, x, y, dv)
         tmpv = matmul(dv,icov)
         xcx = dot_product(dv,tmpv)

         fmultikernel = dexp(-0.5d0*xcx)
      END FUNCTION fmultikernel

      DOUBLE PRECISION FUNCTION fkernel(D,period,sig2,vc,vp)
            ! Calculate the gaussian kernel
            ! The normalization has to be done outside
            !
            ! Args:
            !    sig2: sig**2
            !    vc: voronoi center's vector
            !    vp: point's vector

            INTEGER, INTENT(IN) :: D
            DOUBLE PRECISION, INTENT(IN) :: period(D)
            DOUBLE PRECISION, INTENT(IN) :: sig2
            DOUBLE PRECISION, INTENT(IN) :: vc(D)
            DOUBLE PRECISION, INTENT(IN) :: vp(D)


            fkernel=(1/( (twopi*sig2)**(dble(D)/2) ))* &
                    dexp(-pammr2(D,period,vc,vp)*0.5/sig2)
                    
      END FUNCTION fkernel
      
      DOUBLE PRECISION FUNCTION fkernelvm(kkk,dist)
            ! Calculate the univariate von Mises kernel
            !
            ! Args:
            !    sig2: sig**2
            !    dist: distance between the two points

            DOUBLE PRECISION, INTENT(IN) :: kkk
            DOUBLE PRECISION, INTENT(IN) :: dist
            !DOUBLE PRECISION vv
            
            ! this is just to avoid errors related to the numerical precision
            !vv=1.0d0/sig2
            !IF(DLOG(sig2).lt.(-2)) vv=100.0d0
           
            fkernelvm=DEXP(DCOS(dist)*kkk)!/ &
                      !(BESSI0(vv)*twopi)
                    
      END FUNCTION fkernelvm

      DOUBLE PRECISION FUNCTION genkernel(ngrid,D,period,sig2,prob,vp,vgrid)
            ! Calculate the (normalized) gaussian kernel
            ! in an arbitrary point
            !
            ! Args:
            !    ngrid: number of grid point
            !    D: dimensionality
            !    period: periodicity
            !    sig2: sig**2
            !    prob: probability of the grid points
            !    vp: point's vector
            !    vgrid: grid points

            INTEGER, INTENT(IN) :: D
            INTEGER, INTENT(IN) :: ngrid
            DOUBLE PRECISION, INTENT(IN) :: period(D)
            DOUBLE PRECISION, INTENT(IN) :: sig2
            DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: prob
            DOUBLE PRECISION, INTENT(IN) :: vp(D)
            DOUBLE PRECISION, INTENT(IN) :: vgrid(D,ngrid)
            
            INTEGER j
            DOUBLE PRECISION res!,norm
            
            res=0.0d0
            !norm=0.0d0    
            DO j=1,ngrid
               res=res+(prob(j)/( (twopi*sig2)**(dble(D)/2) ))* &
                   dexp(-0.5d0*pammr2(D,period,vgrid(:,j),vp)/sig2)
               !norm=norm+prob(j)
            ENDDO
            genkernel=res/ngrid 
      END FUNCTION genkernel
      
      SUBROUTINE getNpoint(D,period,np,r1,r2,listpoints)
         ! Get np points in a segment given the extremes
         !
         ! Args:
         !    D          : Dimensionality of a point
         !    period     : periodicity in each dimension
         !    np         : number of points to be generated
         !    r1         : first extreme
         !    r2         : second extreme
         !    listpoints : array of np points
         
         INTEGER, INTENT(IN) :: D,np
         DOUBLE PRECISION, DIMENSION(D), INTENT(IN) :: r1,r2,period
         DOUBLE PRECISION, DIMENSION(D,np), INTENT(OUT) :: listpoints

         DOUBLE PRECISION v12(D)
         INTEGER i

         v12=0.0d0
         CALL pammrij(D,period,r2,r1,v12)
         v12=v12/(np+1)
         
         listpoints=0.0d0
         DO i=1,np
            listpoints(:,i)= r1 + float(i)*v12
         ENDDO
      END SUBROUTINE getNpoint

      SUBROUTINE savevoronois(D,period,nsamples,ngrid,x,y,npvoronoi,iminij,wj,rgrid,prvor)     
         ! Store Voronoi data in a file
         ! 
         ! Args:
         !    D          : Dimensionality of a point
         !    period     : periodicity in each dimension
         !    nsamples   : total points number
         !    ngrid      : number of grid points
         !    x          : array containing the data samples
         !    y          : array that will contain the grid points
         !    npvoronoi  : array cotaing the number of samples inside the Voronoj polyhedron of each grid point
         !    iminij     : array containg the neighbor list for data samples
         !    prvor      : prefix for the outputfile
         
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: nsamples
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, INTENT(IN) :: period(D),rgrid(ngrid)
         DOUBLE PRECISION, DIMENSION(D,nsamples), INTENT(IN) :: x
         DOUBLE PRECISION, DIMENSION(D,ngrid), INTENT(IN) :: y
         INTEGER, DIMENSION(ngrid), INTENT(IN) :: npvoronoi
         INTEGER, DIMENSION(nsamples), INTENT(IN) :: iminij
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: wj
         CHARACTER(LEN=1024), INTENT(IN) :: prvor
         
         INTEGER i,j      

         ! write out the voronoi links
         OPEN(UNIT=12,FILE=trim(prvor)//".voronoislinks",STATUS='REPLACE',ACTION='WRITE')
         ! header
         WRITE(12,*) "# Dimensionality, NSamples // point, voronoiassociation, weight"
         WRITE(12,"((A2,I9,I9))") " #", D, nsamples

         DO j=1,nsamples
            ! write first the point
            DO i=1,D
               WRITE(12,"((A1,ES15.4E4))",ADVANCE="NO") " ", x(i,j)
            ENDDO
            ! write the Voronoi associated
            WRITE(12,"((A1,I9))",ADVANCE="NO") " ", iminij(j)
            ! write the weight
            WRITE(12,"((A1,ES15.4E4))") " ", wj(j)
         ENDDO         

         CLOSE(UNIT=12)

         ! write out the grid
         OPEN(UNIT=12,FILE=trim(prvor)//".voronois",STATUS='REPLACE',ACTION='WRITE')
         ! header
         WRITE(12,*) "# Dimensionality, NGrids // point, n. of points in the Voronoi, radius**2"
         WRITE(12,"((A2,I9,I9))") " #", D, ngrid
         DO i=1,ngrid
            ! write first the grid point
            DO j=1,D
               WRITE(12,"((A1,ES15.4E4))",ADVANCE="NO") " ", y(j,i)
            ENDDO
            ! write the number of points falling into this Voronoi
            WRITE(12,"((A1,I9))",ADVANCE="NO") " ", npvoronoi(i)
            WRITE(12,"((A1,ES15.4E4))") " ", rgrid(i)
         ENDDO
         CLOSE(UNIT=12)
      END SUBROUTINE savevoronois
      
      SUBROUTINE savegrid(outputfile,D,ngrid,y,prob,proberr,Hi)     
         ! Store results from adaptive runs in a file
         ! 
         ! Args:
         !    outputfile     : Dimensionality of a point
         ! ...
         
         CHARACTER(LEN=1024), INTENT(IN) :: outputfile
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, DIMENSION(D,ngrid), INTENT(IN) :: y
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: prob
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: proberr
         DOUBLE PRECISION, DIMENSION(D,D,ngrid), INTENT(IN) :: Hi
         
         INTEGER i,j    

         
         OPEN(UNIT=12,FILE=trim(outputfile)//".probs" &
           ,STATUS='REPLACE',ACTION='WRITE')

         ! header
         WRITE(12,*) "# Dimensionality, NGrids // point, probability, error, sigma"
         WRITE(12,"((A2,I9,I9))") " #", D, ngrid
         DO i=1,ngrid
           ! write first the grid point
           DO j=1,D
             !WRITE(12,"((A1,ES25.9E4))",ADVANCE="NO") " ", y(j,i)
             WRITE(12,"((A1,F10.5))",ADVANCE="NO") " ", y(j,i)
           ENDDO
           ! write the KDE
           WRITE(12,"((A1,ES25.9E4))",ADVANCE="NO") " ", prob(i)
           ! write the error
           WRITE(12,"((A1,ES25.9E4))",ADVANCE="NO") " ", proberr(i)
           WRITE(12,"((A1,ES25.9E4))") " ", DSQRT(trmatrix(D,Hi(:,:,i))/DBLE(D))   
         ENDDO
         CLOSE(UNIT=12)

      END SUBROUTINE savegrid
      
      SUBROUTINE savemat(outputfile,dmn,mat)
         ! Store results from adaptive runs in a file
         ! 
         ! Args:
         !    outputfile     : prefix for the outputfile
         !    matrixname     : name of the matrix to store
         !    dmn            : dimensions of the matrix 
         !                     (only squared matrices allowed)
         !    mat            : the matrix itself
         ! 
         
         CHARACTER(LEN=1024), INTENT(IN) :: outputfile          
         INTEGER, INTENT(IN) :: dmn                
         DOUBLE PRECISION, DIMENSION(dmn,dmn), INTENT(IN) :: mat
         
         INTEGER i,j    

        
         ! write out the grid
         OPEN(UNIT=12,FILE=trim(outputfile), &
              STATUS='REPLACE',ACTION='WRITE')

         DO i=1,dmn
           ! write first the grid point
           WRITE(12,*) " ", mat(i,:)
         ENDDO
         CLOSE(UNIT=12)

      END SUBROUTINE savemat
   END PROGRAM pamm
