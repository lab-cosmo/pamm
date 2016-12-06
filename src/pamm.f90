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

      PROGRAM pamm
      USE libpamm
      USE random
      IMPLICIT NONE
      
!      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

      CHARACTER(LEN=1024) :: outputfile, clusterfile, gridfile  ! The output file prefix
      DOUBLE PRECISION, ALLOCATABLE :: period(:)              ! Periodic lenght in each dimension
      LOGICAL periodic                                        ! flag for using periodic data
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: distmm ! similarity matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: diff     ! temp vector used to store distances

      INTEGER D                                               ! Dimensionality of problem
      INTEGER Nk                                              ! Number of gaussians in the mixture
      INTEGER nsamples                                        ! Total number points
      INTEGER ngrid                                           ! Number of samples extracted using minmax

      INTEGER jmax
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma2, rgrid, wj, prob, bigp, &
                                                     msmu, tmpmsmu, pcluster, px, tmps2, &
                                                     Di
      DOUBLE PRECISION :: normwj                              ! accumulator for wj
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: normvoro ! accumulator for wj in voronoi
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npvoronoi, iminij, pnlist, nlist
      INTEGER seed                                            ! seed for the random number generator
      
      ! variable to set the covariance matrix
      DOUBLE PRECISION tmppks,normpks,tmpkernel
    
      ! Array of Von Mises distributions
      TYPE(vm_type), ALLOCATABLE, DIMENSION(:) :: vmclusters

      ! Array of Gaussians containing the gaussians parameters
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters
      ! Array containing the input data pints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x , y
      ! quick shift, roots and path to reach the root (used to speedup the calculation)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: idxroot, idcls, qspath
      ! cluster connectivity matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: clsadj, clsadjel
      INTEGER, ALLOCATABLE, DIMENSION(:) :: macrocl,sortmacrocl
      INTEGER, ALLOCATABLE, DIMENSION(:) :: idxgrid
      DOUBLE PRECISION linkel,maxrer
      LOGICAL isthere
      
      ! BOOTSTRAP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: probboot
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: prelerr, pabserr
      INTEGER nbootstrap, rndidx, nn, nbssample, nbstot
      DOUBLE PRECISION refcov
      ! Variables for local bandwidth estimation
      INTEGER ntarget
      DOUBLE PRECISION nlocal, prefac , kderr
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xij, normkernel, wQ, pk
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Q, Qlocal, IM
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Hi, Hiinv
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xtmp, xtmpw
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wloc
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ineigh
      
      
      ! additional IN/OUT logical flags
      LOGICAL saveprobs, savevor, saveadj, saveidxs, readprobs, readjgrid, savecovs, accurate
      
      ! PARSER
      CHARACTER(LEN=1024) :: cmdbuffer, comment   ! String used for reading text lines from files
      INTEGER ccmd                                ! Index used to control the input parameters
      INTEGER nmsopt                              ! number of mean-shift optimizations of the cluster centers
      LOGICAL verbose, fpost                      ! flag for verbosity
      LOGICAL weighted                            ! flag for using weigheted data
      INTEGER isep1, isep2, par_count             ! temporary indices for parsing command line arguments
      DOUBLE PRECISION lambda, msw, alpha, zeta, thrmerg, lambda2

      ! Counters and dummy variable
      INTEGER i,j,k,ii,jj,counter,dummyi1,endf
      DOUBLE PRECISION dummd1,dummd2

!!!!!!! Default value of the parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      outputfile = "out"
      clusterfile = "NULL"
      gridfile = "NULL"
      fpost = .false.
      alpha = 1.0d0
      zeta = 0.0d0
      kderr = -1.0d0
      ccmd = 0               ! no parameters specified
      Nk = 0                 ! number of gaussians
      nmsopt = 0             ! number of mean-shift refinements
      ngrid = -1             ! number of samples extracted with minmax
      ntarget = -1           ! number of sample points for localization
      seed = 12345           ! seed for the random number generator
      thrmerg = 0.8d0        ! merge different clusters
      lambda = -1            ! quick shift cut-off
      verbose = .false.      ! no verbosity
      weighted = .false.     ! don't use the weights  
      prefac = 0.0d0         ! prefector of gaussian for localization
      nbootstrap = 0         ! do not use bootstrap
      saveprobs = .false.    ! don't print out the probs
      savevor  = .false.     ! don't print out the Voronoi
      saveidxs = .false.     ! don't save the indexes of the grid points
      saveadj = .false.      ! save adjacency
      readprobs = .FALSE.    ! don't read the probabilities from the standard input
      readjgrid = .FALSE.    ! don't read the grid from the standard input
      savecovs  = .FALSE.
      accurate  = .FALSE.    ! compute the covariance from the grid
      
      D=-1
      periodic=.false.
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
         ELSEIF (cmdbuffer == "-qslambda") THEN     ! cutoff to differentiate different clusters in QuickShift
                ccmd = 5
         ELSEIF (cmdbuffer == "-nms") THEN          ! mean-shift steps
            ccmd = 6
         ELSEIF (cmdbuffer == "-ngrid") THEN        ! N of grid points
            ccmd = 7
         ELSEIF (cmdbuffer == "-bootstrap") THEN    ! estimate the error in the kde using bootstrapping
            ccmd = 8
         ELSEIF (cmdbuffer == "-d") THEN            ! dimensionality
            ccmd = 9
         ELSEIF (cmdbuffer == "-kderr") THEN          ! targete rerror in the kde estimate
            ccmd = 10
         ELSEIF (cmdbuffer == "-p") THEN            ! use periodicity
            ccmd = 11
         ELSEIF (cmdbuffer == "-z") THEN            ! add a background to the probability mixture
            ccmd = 12
         ELSEIF (cmdbuffer == "-ntarget") THEN          ! refine adptively sigma
             ccmd = 14
         ELSEIF (cmdbuffer == "-readprobs") THEN    ! read the grid points from the standard input
            readprobs= .true.
         ELSEIF (cmdbuffer == "-accurate") THEN    ! read the grid points from the standard input
            accurate= .true.
         ELSEIF (cmdbuffer == "-readidxsgrid") THEN    ! read the grid points from the standard input
            readjgrid= .true.
            ccmd = 15
         ELSEIF (cmdbuffer == "-saveprobs") THEN    ! save the KDE estimates in a file
            saveprobs= .true.
         ELSEIF (cmdbuffer == "-savevoronois") THEN ! save the Voronoi associations
            savevor= .true.
         ELSEIF (cmdbuffer == "-saveidxs") THEN ! save the Voronoi associations
            saveidxs= .true.
         ELSEIF (cmdbuffer == "-savecovs") THEN ! save the Voronoi associations
            savecovs= .true.
         ELSEIF (cmdbuffer == "-adj") THEN          ! save the Voronoi associations
            saveadj= .true.
            ccmd = 18
         ELSEIF (cmdbuffer == "-w") THEN            ! use weights
            weighted = .true.
         ELSEIF (cmdbuffer == "-v") THEN            ! verbosity flag
            verbose = .true.
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
            ELSEIF (ccmd == 15) THEN                ! read the file containing the grid indexes                      
               gridfile=trim(cmdbuffer)
            ELSEIF (ccmd == 1) THEN                 ! read the cluster smearing
               READ(cmdbuffer,*) alpha             
            ELSEIF (ccmd == 18) THEN                ! read the threaslod for the cluster adjancency
               READ(cmdbuffer,*) thrmerg
            ELSEIF (ccmd ==5) THEN                  ! read the lambda
               READ(cmdbuffer,*) lambda
               IF (lambda<0) STOP "The QS cutoff should be positive!"
               lambda2=lambda*lambda
            ELSEIF (ccmd == 4) THEN                 ! read the seed for the rng
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 14) THEN                ! read the localization parameter
               READ(cmdbuffer,*) ntarget
            ELSEIF (ccmd == 6) THEN                 ! read the number of mean-shift refinement steps
               READ(cmdbuffer,*) nmsopt
            ELSEIF (ccmd == 7) THEN                 ! number of grid points
               READ(cmdbuffer,*) ngrid
            ELSEIF (ccmd == 8) THEN                 ! read the N of bootstrap iterations
               READ(cmdbuffer,*) nbootstrap
               IF (nbootstrap<0) STOP "The number of iterations should be positive!"
            ELSEIF (ccmd == 9) THEN                 ! read the dimensionality
               READ(cmdbuffer,*) D
               ALLOCATE(period(D))
               period=-1.0d0
            ELSEIF (ccmd == 10) THEN                 ! read the number of mean-shift refinement steps
               READ(cmdbuffer,*) kderr
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
      IF (readprobs) THEN
         CALL readinputprobs(D, ngrid, y, prob, pabserr, prelerr, rgrid)
         ! setup what is needed (mostly to random things)
         nsamples=ngrid
         normwj=ngrid
         CALL allocatevectors2(D,nsamples,nbootstrap,ngrid,accurate,iminij,pnlist,nlist, &
                               npvoronoi,probboot,idxroot,idcls,idxgrid,qspath, &
                               distmm, diff,msmu,tmpmsmu,normkernel, &
                               normvoro,bigp,Q,Qlocal,Hi,Hiinv,IM,xij,pk,wQ, &
                               xtmp,xtmpw,wloc,ineigh,wj,sigma2,tmps2,Di)
         wj=1.0d0
         sigma2=rgrid
         normkernel=1
         normvoro=1
         tmps2=0.0d0
         
         ! estimate Q from the grid
         CALL getweightedcov(D,ngrid,normwj,normvoro,y,Q)
         ! estimate localization from the global covariance
         ! using the Scott's Rule
         ! first we look for the highest value in the diagonal
         refcov=0.0d0
         CALL eigval(Q,D,pk) ! eigenvalues of the covariance matrix
         DO ii=1,D
             IF(pk(ii).GT.dummd1) refcov=pk(ii)
         ENDDO
         
         ! set the lambda to be used in QS
         IF(lambda.LT.0)THEN
         !lambda=5.0d0*median(ngrid,rgrid(:))
         lambda=refcov/5.0d0
         lambda2=lambda*lambda
         ENDIF  

         IF(verbose)THEN
           WRITE(*,*) " Ngrids: ", ngrid
           WRITE(*,*) " QS lambda: ", lambda
         ENDIF
         
         IF(verbose) WRITE(*,*) " Computing similarity matrix"
         DO i=1,ngrid
            distmm(i,i)=0.0d0
            IF(verbose .AND. (modulo(i,1000).EQ.0)) &
               WRITE(*,*) i,"/",ngrid
            DO j=1,i-1  
               ! distance between two voronoi centers
               ! also computes the nearest neighbor distance (squared) for each grid point
               distmm(i,j) = pammr2(D,period,y(:,i),y(:,j))
               distmm(j,i) = distmm(i,j)
            ENDDO
         ENDDO  
         
         GOTO 1111
      ELSE 
         CALL readinput(D, weighted, nsamples, x, normwj, wj)
         ! "renormalizes" the weight so we can consider them sort of sample counts
         IF (weighted) THEN  
             wj = wj * nsamples/sum(wj)
         ENDIF
      ENDIF

      ! If not specified, the number voronoi polyhedra
      ! are set to the square of the total number of points
      IF (ngrid.EQ.-1) ngrid=int(sqrt(float(nsamples)))
      
      ! If not specified, the target local number of sample points
      ! is set to the square of the total number of points
      IF (ntarget.EQ.-1) ntarget=int(sqrt(float(nsamples)))

      ! Initialize the arrays, since now I know the number of
      ! points and the dimensionality
      CALL allocatevectors(D,nsamples,nbootstrap,ngrid,accurate,iminij,pnlist,nlist, &
                           y,npvoronoi,prob,probboot,idxroot,idcls,idxgrid,qspath, &
                           distmm, diff,msmu,tmpmsmu,pabserr,prelerr,normkernel, &
                           normvoro,bigp,Q,Qlocal,Hi,Hiinv,IM,xij,pk,wQ, &
                           xtmp,xtmpw,wloc,ineigh,rgrid,sigma2,tmps2,Di)
      
      ! Extract ngrid points on which the kernel density estimation is to be
      ! evaluated. Also partitions the nsamples points into the Voronoi polyhedra
      ! of the sampling points.
      IF(verbose) THEN
         WRITE(*,*) " NSamples: ", nsamples
         WRITE(*,*) " Selecting ", ngrid, " points using MINMAX"
      ENDIF
      
      IF(readjgrid)THEN
         ! Read the grid 
         
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
         
         WRITE(*,*) " Building the Voronoi associations"
         
         ! do the voronoi associations
         CALL getvoro(D,period,nsamples,ngrid,x,wj,y,npvoronoi,iminij,ineigh,normvoro,idxgrid)
      ELSE
         CALL mkgrid(D,period,nsamples,ngrid,x,wj,y,npvoronoi,iminij,ineigh,normvoro, &
                  saveidxs,outputfile)
      ENDIF
      
      ! print out the voronois associations
      IF(savevor) CALL savevoronois(nsamples,iminij,outputfile)
      
      ! Generate the neighbour list
      IF(verbose) write(*,*) " Generating neighbour list"
      CALL getnlist(nsamples,ngrid,npvoronoi,iminij, pnlist,nlist)
      ! Definition of the distance matrix between grid points
      distmm=0.0d0
      rgrid=1d100 ! "voronoi radius" of grid points (squared)
      
      IF(verbose) WRITE(*,*) & 
        " Computing similarity matrix"
        
      DO i=1,ngrid
         IF(verbose .AND. (modulo(i,1000).EQ.0)) &
           WRITE(*,*) i,"/",ngrid
         DO j=1,i-1
           ! distance between two voronoi centers
           ! also computes the nearest neighbor distance (squared) for each grid point
           distmm(i,j) = pammr2(D,period,y(:,i),y(:,j))
    
           IF (distmm(i,j) < rgrid(i)) rgrid(i) = distmm(i,j)
           ! the symmetrical one
           distmm(j,i) = distmm(i,j)
           IF (distmm(i,j) < rgrid(j)) rgrid(j) = distmm(i,j)   
         ENDDO
      ENDDO 
      
      ! identity matrix
      IM = 0.0d0             
      DO ii=1,D
        IM(ii,ii) = 1.0d0    
      ENDDO
      
      IF(savecovs)THEN
        OPEN(UNIT=11,FILE=trim(outputfile)//".cov",STATUS='REPLACE',ACTION='WRITE')
        WRITE(11,*) "Local covariance (diagonal), Nlocal , Dlocal"
      ENDIF
      
      IF(verbose) WRITE(*,*) & 
        " Estimating localizations and bandwidths"
      
      ! estimate global covariance matrix
      IF(accurate)THEN
        ! estimate Q from the complete dataset
        CALL getweightedcov(D,nsamples,normwj,wj,x,Q)
      ELSE
        ! estimate Q from the grid
        CALL getweightedcov(D,ngrid,normwj,normvoro,y,Q)
      ENDIF
      
      ! get the biggest eigenvalue of Q
      CALL eigval(Q,D,pk) ! eigenvalues of the covariance matrix
      DO ii=1,D
        IF (pk(ii).GT.dummd1) dummd1 = pk(ii)
      ENDDO 
      
      ! use biggest eigenvalue of Q as initial guess
      sigma2 = dummd1
      
      ! estimate the localization for each grid 
      ! point based on the choice of nlocal
      DO i=1,ngrid
        IF(verbose .AND. (modulo(i,100).EQ.0)) &
          WRITE(*,*) i,"/",ngrid
          
        ! first estimate of nlocal
        IF (accurate) THEN          
          CALL getlocalweighted(D,nsamples,-0.5d0/sigma2(i),x,wj,y(:,i),wloc,nlocal)
        ELSE
          CALL getlocalweighted(D,ngrid,-0.5d0/sigma2(i),y,normvoro,y(:,i),wloc,nlocal)
        ENDIF  
       
        ! if nlocal is smaller than target value try to approach quickly to target value
        ! typically the initial sigma is big enough not to do this, however, nobody knows...
        IF (nlocal.LT.ntarget) THEN
          DO WHILE(nlocal.LT.ntarget)
            sigma2(i)=sigma2(i)+dummd1
            IF (accurate) THEN          
              CALL getlocalweighted(D,nsamples,-0.5d0/sigma2(i),x,wj,y(:,i),wloc,nlocal)
            ELSE
              CALL getlocalweighted(D,ngrid,-0.5d0/sigma2(i),y,normvoro,y(:,i),wloc,nlocal)
            ENDIF 
          ENDDO
        ENDIF
        
        j = 1
        ! fine tuning 
        ! TODO: we have to check if sigma gets to small if we want to
        !       do the fast evaluation on the grid maybe this is even
        !       not possible to evaluate on the grid...
        DO WHILE(.TRUE.)  
          IF(nlocal.GT.ntarget) THEN
            sigma2(i) = sigma2(i)-dummd1/2.0d0**j
          ELSE
            sigma2(i) = sigma2(i)+dummd1/2.0d0**j
          ENDIF
          
          IF (accurate) THEN          
            CALL getlocalweighted(D,nsamples,-0.5d0/sigma2(i),x,wj,y(:,i),wloc,nlocal)
          ELSE
            CALL getlocalweighted(D,ngrid,-0.5d0/sigma2(i),y,normvoro,y(:,i),wloc,nlocal)
          ENDIF  
          
          ! exit loop if sigma gives correct nlocal
          IF (ANINT(nlocal).EQ.ntarget) EXIT
          
          ! adjust scaling factor for new sigma
          j = j+1
        ENDDO

        prefac = -0.5d0/sigma2(i)
        ! estimate covariance matrix locally
        IF(accurate)THEN          
          ! estimate weights for localization
          ! using a baloon estimator centered on each grid
          CALL getlocalweighted(D,nsamples,prefac,x,wj,y(:,i),wloc,nlocal)
          ! estimate Q from the complete dataset
          CALL getweightedcov(D,nsamples,nlocal,wloc,x,Qlocal)
        ELSE
          ! estimate weights for localization
          ! using a baloon estimator centered on each grid
          CALL getlocalweighted(D,ngrid,prefac,y,normvoro,y(:,i),wloc,nlocal)
          ! estimate Q from the grid
          CALL getweightedcov(D,ngrid,nlocal,wloc,y,Qlocal)
        ENDIF
        
        ! estimate local dimensionality
        CALL eigval(Qlocal,D,pk) ! eigenvalues of the covariance matrix       
        pk = pk/SUM(pk)
        pk = pk*LOG(pk)
        ! we assume that 0*log(0) is zero
        ! thus we need to check for nan values 
        ! and set pk to zero for that value
        ! since log(x) for x <= 0 is nan
        WHERE( pk .ne. pk ) pk = 0.0d0
        Di(i) = EXP(-SUM(pk))
  
        ! apply oracle approximating shrinkage alogorithm on local Q
        dummd2 = ( (1.0d0-2.0d0/DBLE(D)) * trmatrix(D,Qlocal**2) & 
                  + trmatrix(D,Qlocal)**2.0d0 ) &
               / ( (nlocal + 1.0d0 - 2.0d0/DBLE(D)) &
                  * trmatrix(D,Qlocal**2) &
                  - (trmatrix(D,Qlocal)**2.0d0) / DBLE(D) )
      
        dummd1 = min(1.0d0,dummd2)
        
        ! store for each grid point the regularized covariance matrix
        Qlocal = (1.0d0-dummd1) * Qlocal + dummd1 * trmatrix(D,Qlocal)*IM / DBLE(D)
        
        ! apply scotts rule to oracle covariance matrix
        ! TODO: add dimensionality correction of Scotts rule
        Hi(:,:,i) = Qlocal * ( nlocal**(-1.0d0/(Di(i)+4.0d0)) )**2.0d0 
        
        IF(savecovs)THEN
          DO ii=1,D
            WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", Hi(ii,ii,i)
          ENDDO
          WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ",nlocal
          WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ",Di(i)
          WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ",sigma2(i)
          WRITE(11,*) " "
        ENDIF
        
        ! inverse of the bandwidth matrix
        CALL invmatrix(D,Hi(:,:,i),Hiinv(:,:,i))

        ! estimate the normalization constants
        normkernel(i) = 1.0d0/DSQRT((twopi**DBLE(D))*detmatrix(D,Hi(:,:,i)))
      ENDDO
      IF(savecovs) CLOSE(UNIT=11) 
      
      
      IF(verbose) WRITE(*,*) &
        " Computing kernel density on reference points"
      IF (periodic) THEN
        !!! Kernel Density estimation for periodic data
        prob = 0.0d0
        DO i=1,ngrid
          DO j=1,ngrid
            ! renormalize the distance taking into accout the anisotropy of the multidimensional data 
            IF (mahalanobis(D,period,y(:,i),y(:,j),Hiinv(:,:,j))>16.0d0) THEN
              ! assume distribution in far away grid point is narrow
              prob(i) = prob(i) + normvoro(j) & 
                      * fmultiVM(D,Di(j),period,y(:,i),y(:,j),Hiinv(:,:,j),Hi(:,:,j))
            ELSE
              ! cycle just inside the polyhedra using the neighbour list
              DO k=pnlist(j)+1,pnlist(j+1)
                prob(i) = prob(i) + wj(nlist(k)) & 
                        * fmultiVM(D,Di(j),period,y(:,i),x(:,nlist(k)),Hiinv(:,:,j),Hi(:,:,j))
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        prob=prob/normwj
      ELSE
        !!! Kernel Density estimation for non-periodic data
        prob = 0.0d0
        DO i=1,ngrid
          DO j=1,ngrid
            ! renormalize the distance taking into accout the anisotropy of the multidimensional data 
            IF (mahalanobis(D,period,y(:,i),y(:,j),Hiinv(:,:,j))>16.0d0) THEN
              ! assume distribution in far away grid point is narrow
              prob(i) = prob(i) + normvoro(j) &
                      * fmultigauss(D,period,y(:,i),y(:,j),Hiinv(:,:,j))
            ELSE
              ! cycle just inside the polyhedra using the neighbour list
              DO k=pnlist(j)+1,pnlist(j+1)
                prob(i) = prob(i) + wj(nlist(k)) &
                        * fmultigauss(D,period,y(:,i),x(:,nlist(k)),Hiinv(:,:,j))      
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        prob=prob/normwj
      ENDIF
      
      IF(saveprobs)THEN 
        comment=TRIM(outputfile)
        CALL savegrid(D,ngrid,y,prob,pabserr,prelerr,rgrid,comment)
      ENDIF
          
!      prob = 0.0d0
!      bigp = 0.0d0
!      
!      DO i=1,ngrid
!        DO j=1,ngrid
!          
!          dummd1 = 0.0d0
!          IF(periodic .and. (nbootstrap.eq.0)) tmpcheck=0.0d0
!          
!          ! renormalize the distance taking into accout the anisotropy of the multidimensional data 
!          IF (mahalanobis(D,period,y(:,i),y(:,j),  Hiinv(:,:,j))>16.0d0) THEN
!            ! assume distribution in far away grid point is narrow
!            IF(periodic)THEN
!              dummd2=1.0d0
!              CALL pammrij(D, period, y(:,i), y(:,j), xij)
!              DO jj=1,D
!                 dummd2 = dummd2 * fkernelvm(Hiinv(jj,jj,j),xij(jj))
!              ENDDO
!              tmpkernel = dummd2 * normvoro(j)
!              IF(nbootstrap.eq.0) &
!                tmpcheck = fmultikernel(D,period,y(:,i),y(:,j),Hiinv(:,:,j)) &
!                         * normvoro(j) &
!                         + tmpcheck 
!            ELSE
!              dummd1 = fmultikernel(D,period,y(:,i),y(:,j),Hiinv(:,:,j))*normvoro(j)
!            ENDIF     
!          ELSE
!            ! cycle just inside the polyhedra using the neighbour list
!            DO k=pnlist(j)+1,pnlist(j+1)
!              IF(periodic)THEN
!                dummd2=1.0d0
!                CALL pammrij(D, period, y(:,i),x(:,nlist(k)), xij)
!                DO jj=1,D
!                   dummd2 = dummd2 * fkernelvm(Hiinv(jj,jj,j),xij(jj))
!                ENDDO
!                tmpkernel = dummd2 * wj(nlist(k))
!                IF(nbootstrap.eq.0) &
!                  ! assume a gaussian kernel to estimate the error also with periodic data
!                  tmpcheck = fmultikernel(D,period,y(:,i),x(:,nlist(k)),Hiinv(:,:,j)) &
!                           * wj(nlist(k)) &
!                           + tmpcheck 
!              ELSE
!                tmpkernel = wj(nlist(k))*fmultikernel(D,period,y(:,i),x(:,nlist(k)),Hiinv(:,:,j))
!              ENDIF     
!              ! cumulate the non normalized value of the kernel
!              dummd1 = dummd1 + tmpkernel    
!            ENDDO
!          ENDIF
!            
!          ! get the non-normalized prob, used to compute the estimate of the error
!          IF(periodic .and. (nbootstrap.eq.0)) THEN
!            ! If we are using VM kernels and we are not estimating the error using BS
!            ! then estimate the error using a gaussian kernel
!            bigp(i) = bigp(i) + tmpcheck
!          ELSE
!            bigp(i) = bigp(i) + dummd1
!          ENDIF
!          ! We now have to normalize the kernel
!          prob(i) = prob(i) + dummd1*normkernel(j)
!          
!        ENDDO          
!      ENDDO
!      
!      prob=prob/normwj
!      bigp=bigp/normwj
!   
!      ! use bootstrapping to estimate the error of our KDE

      IF(nbootstrap > 0) probboot = 0.0d0
  
      ! START of bootstrapping run
      DO nn=1,nbootstrap
        IF(verbose) WRITE(*,*) &
              " Bootstrapping, run ", nn

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
            IF (mahalanobis(D,period,y(:,i),y(:,j),Hiinv(:,:,j))>16.0d0) THEN
              ! assume distribution in far away grid point is narrow
              IF(periodic)THEN
                dummd2=1.0d0
                CALL pammrij(D, period, y(:,i), y(:,j), xij)
                DO jj=1,D
                   dummd2 = dummd2 * fkernelvm(Hiinv(jj,jj,j),xij(jj))
                ENDDO
                tmpkernel = dummd2 * nbssample
              ELSE
                tmpkernel = fmultikernel(D,period,y(:,i),y(:,j),Hiinv(:,:,j)) &
                          * nbssample                    
              ENDIF     
            ELSE
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
                            * fkernelvm(Hiinv(jj,jj,j),xij(jj))
                  ENDDO
                  tmpkernel = tmpkernel + dummd2 
                ELSE
                  tmpkernel = tmpkernel &
                            + fmultikernel(D,period,y(:,i),x(:,rndidx), & 
                                           Hiinv(:,:,j))
                ENDIF
              ENDDO
            ENDIF
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
        prelerr = 0.0d0
        pabserr = 0.0d0
        DO i=1,ngrid
          pabserr(i) = DSQRT( SUM( (probboot(i,:) - prob(i))**2.0d0 ) / (nbootstrap-1.0d0) )
          prelerr(i) = pabserr(i) / prob(i)
        ENDDO 
      ELSE
        DO i=1,ngrid  
          !prelerr(i) = ( ( (1.0d0/bigp(i))-1.0d0)*(1.0d0/normwj) )**(0.5d0) 
          !pabserr(i) = ( ( 1.0d0 - bigp(i) ) * (bigp(i)/normwj) )**(0.5d0) 
          prelerr(i)=DSQRT((((((twopi/2.0d0)**(-Di(i)/2.0d0))*(DSQRT(rgrid(i))**(-Di(i))) &
                     *GAMMA(Di(i)/2.0d0))/prob(i))-1.0d0)/(normwj**2.0d0))
          ! I got here the relative error on Ni (point falling into the Voronoi i)
          ! that, propagating the error is equal to the relative error of prob(i).
          ! To get the absolute error we just need to do prelerr(i)*prob(i) 
          pabserr(i)=prelerr(i)*prob(i)
        ENDDO
      ENDIF

      IF(verbose) WRITE(*,*) &
        " Quick-Shift : localization based"  
      
      maxrer=MAXVAL(prelerr)
1111  idxroot=0
      ! Start quick shift
      ! lambda is based on localization
      IF(verbose) WRITE(*,*) " Starting Quick-Shift"
      DO i=1,ngrid
         IF(idxroot(i).NE.0) CYCLE
         IF(verbose .AND. (modulo(i,1000).EQ.0)) &
               WRITE(*,*) i,"/",ngrid
         qspath=0
         qspath(1)=i
         counter=1         
         DO WHILE(qspath(counter).NE.idxroot(qspath(counter)))
         
            ! We use an adpative scheme for the qslambda
            ! I renormalize all the rerr and I use that to amplify the qs lambda locally
!            idxroot(qspath(counter))= &
!                qs_next(ngrid,qspath(counter), & 
!                   (lambda2*(1.0d0+prelerr(i)/maxrer)),prob,distmm)
            idxroot(qspath(counter))= &
                qs_next(ngrid,qspath(counter), & 
                   sigma2(i),prob,distmm)

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
      
      IF(verbose) write(*,*) " Writing out"
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
         WRITE(11,"(A1,I4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4,A1,ES15.4E4)") & 
                                              " " , dummyi1 ,   &
                                              " " , prob(i) ,   &
                                              " " , pabserr(i), &
                                              " " , prelerr(i), &
                                              " ",  trmatrix(D,Hi(:,:,i))/DBLE(D) 
         ! accumulate the normalization factor for the pks
         normpks=normpks+prob(i)
      ENDDO
      CLOSE(UNIT=11)
      
      ! builds the cluster adjacency matrix
      IF(saveadj)THEN
         IF (verbose) WRITE(*,*) "Building cluster adjacency matrix"
         ALLOCATE(clsadj(Nk, Nk),clsadjel(Nk, Nk))
         ALLOCATE(macrocl(Nk))
         clsadj   = 0.0d0
         clsadjel = 0.0d0

         DO i=1, Nk
            IF(verbose .AND. (modulo(i,10).EQ.0)) WRITE(*,*) i,"/",Nk
            ! initialize each cluster to itself in the macrocluster assignation 
            macrocl(i)=i
            DO j=1,i-1
                clsadj(i,j) = cls_link(ngrid, idcls, distmm, prob, rgrid, i, j, &
                                       pabserr, linkel)
                clsadj(j,i) = clsadj(i,j)
                ! adjacency without considering the error
                clsadjel(i,j) = linkel
                clsadjel(j,i) = linkel
            ENDDO
         ENDDO
         
         OPEN(UNIT=11,FILE=trim(outputfile)//".adj",STATUS='REPLACE',ACTION='WRITE')
         OPEN(UNIT=12,FILE=trim(outputfile)//".adjel",STATUS='REPLACE',ACTION='WRITE')
         DO i=1, Nk
             DO j=1, Nk
                 WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", clsadj(i,j)
                 WRITE(12,"((A1,ES15.4E4))",ADVANCE = "NO") " ", clsadjel(i,j)
             ENDDO
             WRITE(11,*) ""
             WRITE(12,*) ""
         ENDDO
         CLOSE(11)
         CLOSE(12)

      
         ! Let's print out the macroclusters
         DO i=1, Nk
            DO j=1, Nk
               ! Put a threshold under which there is no link between the clusters
               ! now it is just a default
               IF(i.EQ.j) CYCLE ! discard yourself
               IF(clsadjel(i,j) .GT. thrmerg) THEN 
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
                  msw = prob(i)*exp(-0.5*pammr2(D,period,y(:,i),vmclusters(k)%mean)/(lambda2/25.0d0))
                  CALL pammrij(D,period,y(:,i),vmclusters(k)%mean,tmpmsmu)
               ELSE
                  msw = prob(i)*exp(-0.5*pammr2(D,period,y(:,i),clusters(k)%mean)/(lambda2/25.0d0))
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
            !! TODO : compute the covariance from the initial samples
            ! use the routine getweightedcov(D,nsamples,nlocal,wloc,x,Qlocal)
            tmppks=tmppks+prob(i)
            xij=0.0d0
            IF(periodic)THEN
               CALL pammrij(D,period,y(:,i),vmclusters(k)%mean,xij)
            ELSE
               CALL pammrij(D,period,y(:,i),clusters(k)%mean,xij)
            ENDIF
            
            DO ii=1,D
               DO jj=1,D
                  IF(periodic)THEN
                     vmclusters(k)%cov(ii,jj)= vmclusters(k)%cov(ii,jj)+prob(i)* &
                                               xij(ii)*xij(jj)               
                  ELSE
                     clusters(k)%cov(ii,jj)=clusters(k)%cov(ii,jj)+prob(i)* &
                                               xij(ii)*xij(jj)
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
      
      DEALLOCATE(x,wj,Di)
      DEALLOCATE(period)
      DEALLOCATE(idxroot,qspath,distmm,idxgrid)
      DEALLOCATE(pnlist,nlist,iminij,bigp)
      DEALLOCATE(y,npvoronoi,prob,sigma2,rgrid,normvoro)
      DEALLOCATE(diff,msmu,tmpmsmu)
      DEALLOCATE(Q,Qlocal,Hi,Hiinv,normkernel,IM)
      DEALLOCATE(xij,wQ,pk,tmps2)
      DEALLOCATE(xtmp,xtmpw,wloc,ineigh)
      IF(saveadj) DEALLOCATE(macrocl,sortmacrocl)
      IF(nbootstrap>0) DEALLOCATE(probboot,prelerr,pabserr)

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
         WRITE(*,*) "   -qslambda lambda  : Quick shift cutoff [automatic] "
         WRITE(*,*) "   -ngrid ngrid      : Number of grid points to evaluate KDE [sqrt(nsamples)]"
         WRITE(*,*) "   -bootstrap N      : Number of iteretions to do when using bootstrapping "
         WRITE(*,*) "                       to refine the KDE on the grid points"
         WRITE(*,*) "   -nms nms          : Do nms mean-shift steps with a Gaussian width lambda/5 to "
         WRITE(*,*) "                       optimize cluster centers [0] "
         WRITE(*,*) "   -kderr target     : Target relative error in the KDE [0.1] "
         WRITE(*,*) "   -seed seed        : Seed to initialize the random number generator. [12345]"
         WRITE(*,*) "   -p P1,...,PD      : Periodicity in each dimension [ (6.28,6.28,6.28,...) ]"
         WRITE(*,*) "   -savevoronois     : Save Voronoi associations. This will produce:"
         WRITE(*,*) "                         output.voronoislinks (points + associated Voronoi) "
         WRITE(*,*) "                         output.voronois (Voronoi centers + info) "
         WRITE(*,*) "   -saveprobs        : Save the KDE estimation on a file"
         WRITE(*,*) "   -readprobs        : Read the grid and the probabilities from the sdtin"
         WRITE(*,*) "   -readidxsgrid     : Read just the grid"
         WRITE(*,*) "   -loc sigma        : Localization width for local bayesian run [automatic] "
         WRITE(*,*) "   -adj threshold    : Set the threshold to merge adjcent clusters and "
         WRITE(*,*) "                       write out the adjacency matrix [default: off] "
         WRITE(*,*) "   -accurate         : Compute the covariances matrix from the entire dataset "
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
      
      SUBROUTINE allocatevectors(D,nsamples,nbootstrap,ngrid,accurate,iminij,pnlist,nlist, &
                                 y,npvoronoi,prob,probboot,idxroot,idcls,idxgrid,qspath, &
                                 distmm, diff,msmu,tmpmsmu,pabserr,prelerr,normkernel, &
                                 normvoro,bigp,Q,Qlocal,Hi,Hiinv,IM,xij,pk,wQ, &
                                 xtmp,xtmpw,wloc,ineigh,rgrid,sigma2,tmps2,Di)
                                 
         INTEGER, INTENT(IN) :: D,nsamples,nbootstrap,ngrid
         LOGICAL, INTENT(IN) :: accurate
         INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT):: iminij,pnlist,nlist,idxroot,idxgrid,qspath
         INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: npvoronoi,idcls,ineigh
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: prob,msmu,tmpmsmu,pk,bigp,normvoro
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: pabserr,prelerr,normkernel,wloc,wQ
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: diff,xij,sigma2,rgrid,tmps2,Di
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: Q,Qlocal,IM,probboot
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: xtmp,xtmpw,y,distmm
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT) :: Hi,Hiinv
         
         
         IF (ALLOCATED(iminij))     DEALLOCATE(iminij)
         IF (ALLOCATED(pnlist))     DEALLOCATE(pnlist)
         IF (ALLOCATED(nlist))      DEALLOCATE(nlist)
         IF (ALLOCATED(y))          DEALLOCATE(y)
         IF (ALLOCATED(npvoronoi))  DEALLOCATE(npvoronoi)
         IF (ALLOCATED(prob))       DEALLOCATE(prob)
         IF (ALLOCATED(probboot))     DEALLOCATE(probboot)
         IF (ALLOCATED(idxroot))    DEALLOCATE(idxroot)
         IF (ALLOCATED(idcls))      DEALLOCATE(idcls)
         IF (ALLOCATED(idxgrid))    DEALLOCATE(idxgrid)
         IF (ALLOCATED(distmm))     DEALLOCATE(distmm)
         IF (ALLOCATED(diff))       DEALLOCATE(diff)
         IF (ALLOCATED(msmu))       DEALLOCATE(msmu)
         IF (ALLOCATED(tmpmsmu))    DEALLOCATE(tmpmsmu)
         IF (ALLOCATED(pabserr))    DEALLOCATE(pabserr)
         IF (ALLOCATED(prelerr))    DEALLOCATE(prelerr)
         IF (ALLOCATED(normvoro))   DEALLOCATE(normvoro)
         IF (ALLOCATED(normkernel)) DEALLOCATE(normkernel)
         IF (ALLOCATED(bigp))       DEALLOCATE(bigp)
         IF (ALLOCATED(Q))          DEALLOCATE(Q)
         IF (ALLOCATED(Qlocal))     DEALLOCATE(Qlocal)
         IF (ALLOCATED(Hi))         DEALLOCATE(Hi)
         IF (ALLOCATED(Hiinv))      DEALLOCATE(Hiinv)
         IF (ALLOCATED(IM))         DEALLOCATE(IM)
         IF (ALLOCATED(xij))        DEALLOCATE(xij)
         IF (ALLOCATED(pk))         DEALLOCATE(pk)
         IF (ALLOCATED(wQ))         DEALLOCATE(wQ)
         IF (ALLOCATED(xtmp))       DEALLOCATE(xtmp)
         IF (ALLOCATED(xtmpw))      DEALLOCATE(xtmpw)
         IF (ALLOCATED(wloc))       DEALLOCATE(wloc)
         IF (ALLOCATED(sigma2))     DEALLOCATE(sigma2)
         IF (ALLOCATED(tmps2))      DEALLOCATE(tmps2)
         IF (ALLOCATED(rgrid))      DEALLOCATE(rgrid)
         IF (ALLOCATED(Di))    DEALLOCATE(Di)

         
         ! Initialize the arrays, since now I know the number of
         ! points and the dimensionality
         ALLOCATE(iminij(nsamples))
         ALLOCATE(pnlist(ngrid+1), nlist(nsamples))
         ALLOCATE(y(D,ngrid), npvoronoi(ngrid), prob(ngrid), sigma2(ngrid), rgrid(ngrid))
         ALLOCATE(idxroot(ngrid), idcls(ngrid), qspath(ngrid), distmm(ngrid,ngrid))
         ALLOCATE(diff(D), msmu(D), tmpmsmu(D))
         ALLOCATE(pabserr(ngrid),prelerr(ngrid),normkernel(ngrid),normvoro(ngrid),bigp(ngrid))
         ! bootstrap probability density array will be allocated if necessary
         IF(nbootstrap > 0) ALLOCATE(probboot(ngrid,nbootstrap))
         ! Allocate variables for local bandwidth estimate
         ALLOCATE(Q(D,D),Qlocal(D,D),Hi(D,D,ngrid),Hiinv(D,D,ngrid),IM(D,D))
         ALLOCATE(xij(D),pk(D),Di(ngrid))
         ALLOCATE(wQ(nsamples),idxgrid(ngrid),tmps2(ngrid))
         IF(accurate)THEN 
            ALLOCATE(xtmp(D,nsamples),xtmpw(D,nsamples),wloc(nsamples))
         ELSE
            ALLOCATE(xtmp(D,ngrid),xtmpw(D,ngrid),wloc(ngrid))
         ENDIF
         ALLOCATE(ineigh(ngrid))
      END SUBROUTINE allocatevectors
      
      SUBROUTINE allocatevectors2(D,nsamples,nbootstrap,ngrid,accurate,iminij,pnlist,nlist, &
                                 npvoronoi,probboot,idxroot,idcls,idxgrid,qspath, &
                                 distmm, diff,msmu,tmpmsmu,normkernel, &
                                 normvoro,bigp,Q,Qlocal,Hi,Hiinv,IM,xij,pk,wQ, &
                                 xtmp,xtmpw,wloc,ineigh,wj,sigma2,tmps2,Di)
         INTEGER, INTENT(IN) :: D,nsamples,nbootstrap,ngrid
         LOGICAL, INTENT(IN) :: accurate
         INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT):: iminij,pnlist,nlist,idxroot,idxgrid,qspath
         INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: npvoronoi,idcls,ineigh
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: msmu,tmpmsmu,pk,bigp,normvoro
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: normkernel,wloc,wQ,wj
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: diff,xij,sigma2,tmps2,Di
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: Q,Qlocal,IM,probboot
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: xtmp,xtmpw,distmm
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT) :: Hi,Hiinv
         
         
         IF (ALLOCATED(iminij))     DEALLOCATE(iminij)
         IF (ALLOCATED(pnlist))     DEALLOCATE(pnlist)
         IF (ALLOCATED(nlist))      DEALLOCATE(nlist)
         IF (ALLOCATED(npvoronoi))  DEALLOCATE(npvoronoi)
         IF (ALLOCATED(probboot))   DEALLOCATE(probboot)
         IF (ALLOCATED(idxroot))    DEALLOCATE(idxroot)
         IF (ALLOCATED(idcls))      DEALLOCATE(idcls)
         IF (ALLOCATED(idxgrid))    DEALLOCATE(idxgrid)
         IF (ALLOCATED(distmm))     DEALLOCATE(distmm)
         IF (ALLOCATED(diff))       DEALLOCATE(diff)
         IF (ALLOCATED(msmu))       DEALLOCATE(msmu)
         IF (ALLOCATED(tmpmsmu))    DEALLOCATE(tmpmsmu)
         IF (ALLOCATED(normvoro))   DEALLOCATE(normvoro)
         IF (ALLOCATED(normkernel)) DEALLOCATE(normkernel)
         IF (ALLOCATED(bigp))       DEALLOCATE(bigp)
         IF (ALLOCATED(Q))          DEALLOCATE(Q)
         IF (ALLOCATED(Qlocal))     DEALLOCATE(Qlocal)
         IF (ALLOCATED(Hi))         DEALLOCATE(Hi)
         IF (ALLOCATED(Hiinv))      DEALLOCATE(Hiinv)
         IF (ALLOCATED(IM))         DEALLOCATE(IM)
         IF (ALLOCATED(xij))        DEALLOCATE(xij)
         IF (ALLOCATED(pk))         DEALLOCATE(pk)
         IF (ALLOCATED(wQ))         DEALLOCATE(wQ)
         IF (ALLOCATED(xtmp))       DEALLOCATE(xtmp)
         IF (ALLOCATED(xtmpw))      DEALLOCATE(xtmpw)
         IF (ALLOCATED(wloc))       DEALLOCATE(wloc)
         IF (ALLOCATED(sigma2))     DEALLOCATE(sigma2)
         IF (ALLOCATED(tmps2))      DEALLOCATE(tmps2)
         IF (ALLOCATED(wj))         DEALLOCATE(wj)
         IF (ALLOCATED(Di))    DEALLOCATE(Di)

         
         ! Initialize the arrays, since now I know the number of
         ! points and the dimensionality
         ALLOCATE(iminij(nsamples))
         ALLOCATE(pnlist(ngrid+1), nlist(nsamples))
         ALLOCATE(npvoronoi(ngrid), sigma2(ngrid))
         ALLOCATE(idxroot(ngrid), idcls(ngrid), qspath(ngrid), distmm(ngrid,ngrid))
         ALLOCATE(diff(D), msmu(D), tmpmsmu(D))
         ALLOCATE(normkernel(ngrid),normvoro(ngrid),bigp(ngrid))
         ! bootstrap probability density array will be allocated if necessary
         IF(nbootstrap > 0) ALLOCATE(probboot(ngrid,nbootstrap))
         ! Allocate variables for local bandwidth estimate
         ALLOCATE(Q(D,D),Qlocal(D,D),Hi(D,D,ngrid),Hiinv(D,D,ngrid),IM(D,D))
         ALLOCATE(xij(D),pk(D),tmps2(ngrid),Di(ngrid))
         ALLOCATE(wQ(nsamples),wj(nsamples),idxgrid(ngrid))
         IF(accurate)THEN 
            ALLOCATE(xtmp(D,nsamples),xtmpw(D,nsamples),wloc(nsamples))
         ELSE
            ALLOCATE(xtmp(D,ngrid),xtmpw(D,ngrid),wloc(ngrid))
         ENDIF
         ALLOCATE(ineigh(ngrid))
      END SUBROUTINE allocatevectors2
      
      SUBROUTINE getlocalweighted(D,nps,prefac,vp,weights,pt,wloc,nloc)
         INTEGER, INTENT(IN) :: D,nps
         DOUBLE PRECISION, INTENT(IN) :: weights(nps),prefac
         DOUBLE PRECISION, INTENT(IN) :: vp(D,nps),pt(D)
         DOUBLE PRECISION, INTENT(OUT) :: wloc(nps),nloc
         
         INTEGER ii
         DOUBLE PRECISION tmpx(D,nps)
         
         DO ii=1,D
           tmpx(ii,:) = vp(ii,:)-pt(ii)
         ENDDO
         ! estimate weights for localization as product from 
         ! spherical gaussian weights and weights in voronoi
         wloc = EXP(prefac*SUM(tmpx*tmpx,1)) * weights
         ! estimate local number of sample points
         nloc = SUM(wloc)
      END SUBROUTINE getlocalweighted
      
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

      SUBROUTINE getweightedcov(D,nsamples,normweights,weights,samples,covm)
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: nsamples
         DOUBLE PRECISION, INTENT(IN) :: normweights
         DOUBLE PRECISION, INTENT(IN) :: weights(nsamples)
         DOUBLE PRECISION, INTENT(IN) :: samples(D,nsamples)
         DOUBLE PRECISION, INTENT(OUT) :: covm(D,D)
         
         DOUBLE PRECISION samplestmp(D,nsamples),wstmp(D,nsamples)
         INTEGER ii
         
         DO ii=1,D
           samplestmp(ii,:) = samples(ii,:) - SUM(samples(ii,:)*weights)/normweights
           wstmp(ii,:) = samplestmp(ii,:) * weights/normweights
         ENDDO
         CALL DGEMM("N", "T", D, D, nsamples, 1.0d0, samplestmp, D, wstmp, D, 0.0d0, covm, D)
         covm = covm / (1.0d0-SUM((weights/normweights)**2.0d0))   
      END SUBROUTINE getweightedcov
      
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

      SUBROUTINE mkgrid(D,period,nsamples,ngrid,x,wj,y,npvoronoi,iminij, &
                        ineigh,normvoro,saveidx,ofile)
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
         DOUBLE PRECISION, DIMENSION(nsamples), INTENT(IN) :: wj 
         
         DOUBLE PRECISION, DIMENSION(D,ngrid), INTENT(OUT) :: y
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: npvoronoi
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: ineigh
         INTEGER, DIMENSION(nsamples), INTENT(OUT) :: iminij
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(OUT) :: normvoro 
         CHARACTER(LEN=1024), INTENT(IN) :: ofile   
         LOGICAL, INTENT(IN) :: saveidx   

         INTEGER i,j,irandom
         DOUBLE PRECISION :: dminij(nsamples), dij, dmax, dneigh
         
         iminij=0
         y=0.0d0
         npvoronoi=0
         ! choose randomly the first point
         irandom=int(RAND()*nsamples)
         IF(saveidx) THEN
            OPEN(UNIT=12,FILE=trim(ofile)//".idxs",STATUS='REPLACE',ACTION='WRITE')
            WRITE(12,"((I9))") irandom
         ENDIF
         y(:,1)=x(:,irandom)
         dminij = 1.0d99
         iminij = 1     
         ineigh = 0    
         DO i=2,ngrid
            dmax = 0.0d0
            dneigh = 1.0d99
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
            IF(verbose .AND. (modulo(i,1000).EQ.0)) &
               write(*,*) i,"/",ngrid
         ENDDO

         ! finishes Voronoi attribution
         dneigh = 1.0d99
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
         npvoronoi = 0
         normvoro  = 0.0d0
         DO j=1,nsamples
            npvoronoi(iminij(j))=npvoronoi(iminij(j))+1
            normvoro(iminij(j))=normvoro(iminij(j))+wj(iminij(j))
         ENDDO
      END SUBROUTINE mkgrid
      
      SUBROUTINE getvoro(D,period,nsamples,ngrid,x,wj,y,npvoronoi,iminij, &
                         ineigh,normvoro,idxgrid)
         IMPLICIT NONE
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
         DOUBLE PRECISION, DIMENSION(nsamples), INTENT(IN) :: wj 
         
         DOUBLE PRECISION, DIMENSION(D,ngrid), INTENT(OUT) :: y
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: npvoronoi
         INTEGER, DIMENSION(ngrid), INTENT(OUT) :: ineigh
         INTEGER, DIMENSION(nsamples), INTENT(OUT) :: iminij
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(OUT) :: normvoro   
         INTEGER, DIMENSION(ngrid), INTENT(IN) :: idxgrid 

         INTEGER i,j
         DOUBLE PRECISION :: dminij(nsamples), dij, dmax, dneigh

         iminij=0
         y=0.0d0
         npvoronoi=0
         ! choose randomly the first point
         
         dminij = 1.0d99
         iminij = 1         
         DO i=1,ngrid
            IF(modulo(i,1000).EQ.0) WRITE(*,*) i,"/",ngrid
            dmax = 0.0d0
            dneigh = 1.0d99
            y(:,i)=x(:,idxgrid(i))
            DO j=1,nsamples
               dij = pammr2(D,period,y(:,i),x(:,j))
               IF (dminij(j)>dij) THEN
                  dminij(j) = dij
                  iminij(j) = i ! also keeps track of the Voronoi attribution
               ENDIF
               IF (dminij(j) > dmax) THEN
                  dmax = dminij(j)
               ENDIF
               IF ((dneigh > dij) .and. (dij .ne. 0.0d0)) THEN
                  dneigh = dij
                  ineigh(i) = j ! store index of closest sample neighbor to grid point
               ENDIF
            ENDDO           
         ENDDO

         ! Assign neighbor list pointer of voronois
         ! Number of points in each voronoi polyhedra
         npvoronoi = 0
         normvoro  = 0.0d0
         DO j=1,nsamples
            npvoronoi(iminij(j))=npvoronoi(iminij(j))+1
            normvoro(iminij(j))=normvoro(iminij(j))+wj(iminij(j))
         ENDDO
      END SUBROUTINE getvoro

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
         DO i=1, ngrid
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
      
      DOUBLE PRECISION FUNCTION fmultigauss(D,period,x,y,icov)
         ! Return the multivariate gaussian density
         ! Args:
         !    gpars: gaussian parameters
         !    x: point in wich estimate the value of the gaussian
         
         INTEGER , INTENT(IN) :: D
         DOUBLE PRECISION, INTENT(IN) :: period(D)
         DOUBLE PRECISION, INTENT(IN) :: x(D)
         DOUBLE PRECISION, INTENT(IN) :: y(D)
         DOUBLE PRECISION, INTENT(IN) :: icov(D,D)
         
         fmultigauss = fmultikernel(D,period,x,y,icov)/ &
                      DSQRT((twopi**DBLE(D))*detmatrix(D,Hi(:,:,i)))
         
      END FUNCTION fmultigauss
      
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
            fmultiVM=fmultikernel(D,period,x,y,icov)/ &
                      DSQRT((twopi**DBLE(D))*detmatrix(D,Hi(:,:,i)))
         ELSE
            ! productkernels
            fmultiVM=1.0d0
            CALL pammrij(D, period, x, y, dv)
            DO i = D, (D-effD), -1
               IF(DSQRT(ev(jj)).LT.0.6d0)THEN
                 fmultiVM = fmultiVM * (1.0d0/((twopi*ev(jj))**0.5d0))* &
                       dexp(-0.5d0*(dv(jj)**2.0d0)/ev(jj))
               ELSE
                 fmultiVM = fmultiVM * fkernelvm(1.0d0/ev(jj),dv(jj))
               ENDIF
            END DO
         ENDIF
         
      END FUNCTION fmultiVM
      
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
           
            fkernelvm=DEXP(DCOS(dist)*kkk) / &
                      (BESSI0(kkk)*twopi)
                    
      END FUNCTION fkernelvm

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
