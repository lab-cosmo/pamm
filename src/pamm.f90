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
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma2, rgrid, wj, probnmm, &
                                                     msmu, tmpmsmu, pcluster, px
      DOUBLE PRECISION :: normwj                              ! accumulator for wj
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npvoronoi, iminij, pnlist, nlist
      INTEGER seed                                            ! seed for the random number generator

      ! variable to set the covariance matrix
      DOUBLE PRECISION tmppks,normpks
    
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
      
      ! BOOTSTRAP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: probboot
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: errprobnmm
      INTEGER nbootstrap,rndidx,rngidx,nn
      DOUBLE PRECISION tmperr,normboot
      
      ! PARSER
      CHARACTER(LEN=1024) :: cmdbuffer, comment   ! String used for reading text lines from files
      INTEGER ccmd                                ! Index used to control the input parameters
      INTEGER nmsopt                              ! number of mean-shift optimizations of the cluster centers
      LOGICAL verbose, fpost                      ! flag for verbosity
      LOGICAL weighted                            ! flag for using weigheted data
!      LOGICAL adaptive                            ! flag for adaptively refine the sigmas
      LOGICAL neblike                             ! neblike path search
      INTEGER isep1, isep2, par_count             ! temporary indices for parsing command line arguments
      DOUBLE PRECISION lambda, lambda2, msw, alpha, zeta, kderr, dummd1,dummd2

      INTEGER i,j,k,ikde,counter,dummyi1,endf ! Counters and dummy variable

!!!!!!! Default value of the parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      outputfile="out"
      clusterfile="NULL"
      fpost=.false.
      alpha=1.0d0
      zeta=0.0d0
      ccmd=0              ! no parameters specified
      Nk=0                ! number of gaussians
      nmsopt=0            ! number of mean-shift refinements
      ngrid=-1            ! number of samples extracted with minmax
      seed=12345          ! seed for the random number generator
      kderr=0.1           ! target fractional error for KDE smoothing
      lambda=-1           ! quick shift cut-off
      verbose = .false.   ! no verbosity
      weighted= .false.   ! don't use the weights
!      adaptive= .false.   ! don't use the adaptive
      neblike= .false.    ! don't use neb paths
      nbootstrap=0        ! do not use bootstrap
      
      D=-1
      periodic=.false.
      CALL random_init(seed) ! initialize random number generator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-o") THEN           ! output file
            ccmd = 2
         ELSEIF (cmdbuffer == "-a") THEN       ! cluster smearing
            ccmd = 1
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
         ELSEIF (cmdbuffer == "-kde") THEN       ! dimensionality
            ccmd = 10
         ELSEIF (cmdbuffer == "-p") THEN       ! use periodicity
            ccmd = 11
         ELSEIF (cmdbuffer == "-z") THEN       ! add a background to the probability mixture
            ccmd = 12
         ELSEIF (cmdbuffer == "-w") THEN       ! use weights
            weighted = .true.
!         ELSEIF (cmdbuffer == "-adaptive") THEN  ! refine adptively sigma2
!            adaptive = .true.
         ELSEIF (cmdbuffer == "-neblike") THEN  ! refine adptively sigma2
            neblike = .true.
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
            ELSEIF (ccmd == 2) THEN            ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN            ! model file
               fpost=.true.
               clusterfile=trim(cmdbuffer)
            ELSEIF (ccmd == 1) THEN            ! read the cluster smearing
               READ(cmdbuffer,*) alpha
            ELSEIF (ccmd == 10) THEN           ! read the cluster smearing
               READ(cmdbuffer,*) kderr
               IF (kderr<0) STOP "The error should be a positive number!"
            ELSEIF (ccmd == 8) THEN            ! read the num of bootstrap iterations
               READ(cmdbuffer,*) nbootstrap
               IF (nbootstrap<0) STOP "The number of iterations should be positive!"
            ELSEIF (ccmd == 9) THEN            ! read the dimensionality
               READ(cmdbuffer,*) D
               ALLOCATE(period(D))
               period=-1.0d0
            ELSEIF (ccmd == 4) THEN ! read the seed
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 6) THEN ! read the number of mean-shift steps
               READ(cmdbuffer,*) nmsopt
            ELSEIF (ccmd == 5) THEN ! cut-off for quick-shift
               READ(cmdbuffer,*) lambda
            ELSEIF (ccmd == 7) THEN ! number of grid points
               READ(cmdbuffer,*) ngrid
            ELSEIF (ccmd == 12) THEN ! read zeta 
               READ(cmdbuffer,*) zeta
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
               
               IF (par_count/=D) STOP "Check the number of periods (-p)!"
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
      IF (weighted) THEN  ! "renormalizes" the weight so we can consider them sort of sample counts
         wj = wj * nsamples/sum(wj)
      ENDIF

      ! If not specified, the number voronoi polyhedras
      ! are set to the square of the total number of points
      IF (ngrid.EQ.-1) ngrid=int(sqrt(float(nsamples)))

      ! Initialize the arrays, since now I know the number of
      ! points and the dimensionality
      ALLOCATE(iminij(nsamples))
      ALLOCATE(pnlist(ngrid+1),nlist(nsamples))
      ALLOCATE(y(D,ngrid),npvoronoi(ngrid),probnmm(ngrid),sigma2(ngrid),rgrid(ngrid))
      ALLOCATE(idxroot(ngrid),idcls(ngrid),qspath(ngrid),distmm(ngrid,ngrid))
      ALLOCATE(diff(D), msmu(D), tmpmsmu(D))

      ! Extract ngrid points on which the kernel density estimation is to be
      ! evaluated. Also partitions the nsamples points into the Voronoi polyhedra
      ! of the sampling points.
      IF(verbose) THEN
         WRITE(*,*) "NSamples: ", nsamples
         WRITE(*,*) "Selecting ", ngrid, " points using MINMAX"
      ENDIF

      CALL mkgrid(D,period,nsamples,ngrid,x,y,npvoronoi,iminij)

      ! Generate the neighbour list
      IF(verbose) write(*,*) "Generating neighbour list"
      CALL getnlist(nsamples,ngrid,npvoronoi,iminij, pnlist,nlist)
      
      ! Definition of the distance matrix between grid points
      distmm=0.0d0
      rgrid=1d100 ! "voronoi radius" of grid points (squared)
      IF(verbose) write(*,*) "Computing similarity matrix"
      DO i=1,ngrid

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
      
      IF(lambda.EQ.-1)THEN
         ! set automatically the mean shift lambda set to 5*<sig>
         lambda2=SUM(rgrid)/ngrid
         lambda=5.0d0*dsqrt(lambda2)
      ENDIF
      lambda2=lambda*lambda ! we always work with squared distances....

      sigma2 = rgrid ! initially set KDE smearing to the nearest grid distance
      
      ikde = 0
!100   IF(verbose) WRITE(*,*) &
      IF(verbose) WRITE(*,*) &
          "Computing kernel density on reference points."

      IF(nbootstrap>0) THEN
          ALLOCATE(probboot(ngrid,nbootstrap))
          ALLOCATE(errprobnmm(ngrid))
          probboot = 0.0d0
          DO nn=1,nbootstrap
              IF(verbose) WRITE(*,*) &
                    "Bootstrapping, run ", nn
              DO i=1,ngrid
                  normboot=0.0d0
                  ! Instead of cycling on all the point we cycle on the Voronoi 
                  DO j=1,ngrid
                     ! do not compute KDEs for points that belong to far away Voronoi
                     IF (distmm(i,j)/sigma2(j)>36.0d0) CYCLE
                     ! cycle just inside the polyhedra using the neighbour list
                     !! BOOTSTRAPPING
                     rngidx=pnlist(j+1)-(pnlist(j)+1)
                     DO k=pnlist(j)+1,pnlist(j+1)
                         rndidx=int(rngidx*random_uniform())+pnlist(j)+1              
                         probboot(i,nn)=probboot(i,nn)+ wj(nlist(rndidx))* &
                           fkernel(D,period,sigma2(j),y(:,i),x(:,nlist(rndidx)))
                         normboot=normboot+wj(nlist(rngidx))
                     ENDDO
                  ENDDO         
                  probboot(i,nn)=probboot(i,nn)/normboot
              ENDDO
          ENDDO 
          ! Average the estimates and get an error bar
          errprobnmm = 0.0d0
          probnmm    = 0.0d0
          DO i=1,ngrid
              ! get the mean
              normboot=0.0d0
              tmperr=0.0d0
              DO nn=1,nbootstrap
                  normboot=normboot+probboot(i,nn)
                  tmperr=tmperr+probboot(i,nn)**2
              ENDDO
              probnmm(i)=normboot/nbootstrap
              ! get the s**2
              errprobnmm(i)=(tmperr-(normboot*normboot)/nbootstrap)/(nbootstrap-1.0d0)
              ! get the SME=s/sqrt(N)
              errprobnmm(i)=DSQRT(errprobnmm(i)/nbootstrap)
          ENDDO
      ELSE
          ! computes the KDE on the Voronoi centers using the neighbour list
          probnmm = 0.0d0
          DO i=1,ngrid
              DO j=1,ngrid
          
                 ! do not compute KDEs for points that belong to far away Voronoi
                 IF (distmm(i,j)/sigma2(j)>36.0d0) CYCLE
              
                 ! cycle just inside the polyhedra using the neighbour list
                 DO k=pnlist(j)+1,pnlist(j+1)
                     probnmm(i)=probnmm(i)+ wj(nlist(k))* &
                                fkernel(D,period,sigma2(j),y(:,i),x(:,nlist(k)))
                 ENDDO
              ENDDO
              probnmm(i)=probnmm(i)/normwj
          ENDDO 
      ENDIF

! I think this adaptive thing can be removed now
!!!#################################      
!!      IF (adaptive) THEN   
!!        ! compute sigma2 from the number of points in the 
!!        ! neighborhood of each grid point. the rationale is that 
!!        ! integrating over a Gaussian kernel with variance sigma2
!!        ! will cover a volume Vi=(2\pi sigma2)^(D/2).  
!!        ! If the estimate probability density on grid point i is pi
!!        ! then the probability to be within that (smooth) bin is ri=pi * Vi. 
!!        ! The number of points is then normwj*ri, and we can take this to 
!!        ! correspond to the mean of a binomial distribution. The (squared) relative error
!!        ! is then err2=(1-ri)/(ri normwj) so
!!        ! kderr^2 *normwj + 1 = 1/ri so
!!        ! Vi = 1/[pi*(1+kderr^2 normwj)] so 
!!        ! sigma2 = 1/(2pi) 1/rad[D/2](pi*(1+kderr^2 normwj))
!!        DO j=1,ngrid
!!            !IF(verbose) WRITE(*,*) "Update grid point ", j, sigma2(j)
!!            sigma2(j) = 1/twopi *1/( probnmm(j)*(1+normwj*kderr*kderr))**(D/2)
!!            ! kernel density estimation cannot become smaller than the distance with the nearest grid point
!!            IF (sigma2(j).lt.rgrid(j)) sigma2(j)=rgrid(j)
!!            !IF(verbose) WRITE(*,*) "Prob ", probnmm(j),  " new sigma ", sigma2(j), "rgrid", rgrid(j)      
!!        ENDDO
!!        ikde = ikde+1
!!        if (ikde<2) GOTO 100 ! seems one can actually iterate to self-consistency....
!!      ENDIF
      
      ! CLUSTERING, local maxima search
      IF(verbose) write(*,*) "Running quick shift"
!!!!!!      IF(verbose) write(*,*) "Lambda: ", lambda

      idxroot=0
      ! Start quick shift
      DO i=1,ngrid
         IF(idxroot(i).NE.0) CYCLE
         qspath=0
         qspath(1)=i
         counter=1         
         DO WHILE(qspath(counter).NE.idxroot(qspath(counter)))
            ! if we are using bootstrapping we use the errors to understand
            ! what's going on
            IF(nbootstrap>0) THEN
                idxroot(qspath(counter))= &
                 qs_next(ngrid,qspath(counter),probnmm,distmm,rgrid,kderr,lambda, & 
                         verbose,neblike,nbootstrap,errprobnmm)
            ELSE ! if not, we do it in the pammv1 way
                idxroot(qspath(counter))= &
                 qs_next(ngrid,qspath(counter),probnmm,distmm,rgrid,kderr,lambda, & 
                         verbose,neblike,nbootstrap)
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
         WRITE(11,"(A1,I4,A1,ES15.4E4,A1,ES15.4E4)") " ", dummyi1 , " ", probnmm(i), " ", errprobnmm(i) 
         ! accumulate the normalization factor for the pks
         normpks=normpks+probnmm(i)
      ENDDO
      CLOSE(UNIT=11)
      
      ! builds the cluster adjacency matrix
      IF (verbose) WRITE(*,*) "Building cluster adjacency matrix"
      ALLOCATE(clsadj(Nk, Nk))
      clsadj = 0.0d0
      DO i=1, Nk
          DO j=1,i-1
              clsadj(i,j) = cls_link(ngrid, idcls, distmm, probnmm, rgrid, i, j)
              clsadj(j,i) = clsadj(i,j)
          ENDDO
      ENDDO
      OPEN(UNIT=11,FILE=trim(outputfile)//".adj",STATUS='REPLACE',ACTION='WRITE')
      DO i=1, Nk
          DO j=1, Nk
              WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", clsadj(i,j)
          ENDDO
          WRITE(11,*) ""
      ENDDO
      CLOSE(11)
      
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
                  msw = probnmm(i)*exp(-0.5*pammr2(D,period,y(:,i),vmclusters(k)%mean)/(lambda2/25.0d0))
                  CALL pammrij(D,period,y(:,i),vmclusters(k)%mean,tmpmsmu)
               ELSE
                  msw = probnmm(i)*exp(-0.5*pammr2(D,period,y(:,i),clusters(k)%mean)/(lambda2/25.0d0))
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
            tmppks=tmppks+probnmm(i)
            DO ii=1,D
               DO jj=1,D
                  IF(periodic)THEN
                     ! Minimum Image Convention
                     dummd1 = (y(ii,i) - vmclusters(k)%mean(ii)) / period(ii)
                     dummd1 = dummd1 - dnint(dummd1)                          
                     dummd1 = dummd1 * period(ii)
                     
                     dummd2 = (y(jj,i)-vmclusters(k)%mean(jj)) / period(jj)
                     dummd2 = dummd2 - dnint(dummd2)                          
                     dummd2 = dummd2 * period(jj)
                     
                     vmclusters(k)%cov(ii,jj)= vmclusters(k)%cov(ii,jj)+probnmm(i)* &
                                               dummd1 * dummd2
                  ELSE
                     clusters(k)%cov(ii,jj)=clusters(k)%cov(ii,jj)+probnmm(i)* &
                       (y(ii,i)-clusters(k)%mean(ii))*(y(jj,i)-clusters(k)%mean(jj))
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
      DEALLOCATE(pnlist,nlist,iminij)
      DEALLOCATE(y,npvoronoi,probnmm,sigma2,rgrid)
      DEALLOCATE(diff,msmu,tmpmsmu)
      IF(nbootstrap>0) DEALLOCATE(probboot,errprobnmm)

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
         WRITE(*,*) "                            output.grid (clusterized grid points) "
         WRITE(*,*) "                            output.pamm (cluster parameters) "
         WRITE(*,*) "   -l lambda         : Quick shift cutoff [automatic] "
         WRITE(*,*) "                       (Not used with the -neblike flag active)"
         WRITE(*,*) "   -kde err          : Target fractional error for KDE smoothing [0.1]"
         WRITE(*,*) "   -ngrid ngrid      : Number of grid points to evaluate KDE [sqrt(nsamples)]"
         WRITE(*,*) "   -bootstrap N      : Number of iteretions to do when using bootstrapping "
         WRITE(*,*) "                       to refine the KDE on the grid points"
         WRITE(*,*) "   -neblike          : Try to improve the clustering using neblike path search algorithm "
         WRITE(*,*) "   -nms nms          : Do nms mean-shift steps with a Gaussian width lambda/5 to"
         WRITE(*,*) "                       optimize cluster centers [0] "
         WRITE(*,*) "   -seed seed        : Seed to initialize the random number generator. [12345]"
         WRITE(*,*) "   -p P1,...,PD      : Periodicity in each dimension [ (6.28,6.28,6.28,...) ]"
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

      SUBROUTINE mkgrid(D,period,nsamples,ngrid,x,y,npvoronoi,iminij)
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
            ENDDO            
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

      DOUBLE PRECISION FUNCTION cls_link(ngrid, idcls, distmm, probnmm, rgrid, ia, ib)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngrid, idcls(ngrid), ia, ib
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid, ngrid), probnmm(ngrid), rgrid(ngrid)
         
         INTEGER i, j
         DOUBLE PRECISION mxa, mxb, mxab, pab
         mxa = 0
         mxb = 0
         mxab = 0
         DO i=1, ngrid
            IF (idcls(i)/=ia) CYCLE
            IF (probnmm(i).gt.mxa) mxa = probnmm(i) ! also gets the probability density at the mode of cluster a
            DO j=1,ngrid
               IF (idcls(j)/=ib) CYCLE
               IF (probnmm(j).gt.mxb) mxb = probnmm(j)
               ! Ok, we've got a matching pair
               IF (dsqrt(distmm(i,j))<dsqrt(rgrid(i))+dsqrt(rgrid(j))) THEN
                  ! And they are close together!
                  pab = (probnmm(i)+probnmm(j))/2
                  IF (pab .gt. mxab) mxab = pab
               ENDIF               
            ENDDO            
         ENDDO
         cls_link = mxab/min(mxa,mxb) 
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
            !   IF ((probnmm(ic)-probnmm(idx)/probnmm(idx))<0) &
            !      write(*,*) "going down", probnmm(ic), probnmm(idx)                  
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
      
      SUBROUTINE neb_path(ngrid, path, npath, distmm, probnmm, w) 
         ! "refines" path to look a bit like a maximum probability path going through saddle points
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngrid, npath
         INTEGER, INTENT(INOUT) :: path(ngrid)
         DOUBLE PRECISION, INTENT(IN) :: distmm(ngrid,ngrid), probnmm(ngrid), w
         
         
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
           dmin = probnmm(1)
           npi = 1
           DO j=2, npath
               IF (probnmm(j)<dmin) THEN
                  dmin=probnmm(j)
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
                       wi * probnmm(j)/probnmm(path(i))
              ELSE
                  dcomb = (dac+dbc)/dab-1 + & 
                       wi * probnmm(j)/probnmm(path(i))
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
      
      INTEGER FUNCTION qs_next(ngrid,idx,probnmm,distmm,rgrid,kderr,lambda, &
                               verbose,neblike,nbootstrap,errors)
         ! Return the index of the closest point higher in P
         !
         ! Args:
         !    ngrid: number of grid points
         !    idx: current point
         !    lambda: cut-off in the jump
         !    probnmm: density estimations
         !    distmm: distances matrix
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ngrid
         DOUBLE PRECISION, INTENT(IN), DIMENSION(ngrid) :: rgrid
         INTEGER, INTENT(IN) :: idx,nbootstrap
         DOUBLE PRECISION, INTENT(IN) :: kderr,lambda
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: probnmm
         DOUBLE PRECISION, DIMENSION(ngrid,ngrid), INTENT(IN) :: distmm
         LOGICAL, INTENT(IN) :: verbose
         LOGICAL, INTENT(IN) :: neblike
         DOUBLE PRECISION, INTENT(IN), DIMENSION(ngrid), OPTIONAL :: errors

         INTEGER i, j, ndx
         
         INTEGER :: npath
         INTEGER :: path(ngrid)
         DOUBLE PRECISION :: dmin,relerr
            
         dmin=1.0d100
         
! find the nearest point with higher probability         
         IF(neblike)THEN
            ndx=idx
            DO j=1,ngrid       
               IF(probnmm(j)>probnmm(idx))THEN 
                  IF(nbootstrap>0)THEN     
                    ! ok, chek the error associated
                    relerr=DSQRT(errors(j)**2+errors(idx)**2)/(probnmm(j)-probnmm(idx))
                    IF(relerr>kderr) CONTINUE
                  ENDIF
                  IF(distmm(idx,j).lt.dmin)THEN
                      ndx=j
                      dmin=distmm(idx,j)
                  ENDIF
               ENDIF
            ENDDO
                     
            qs_next = ndx   
            IF (ndx/=idx) THEN         
               npath = 0
               ! builds a discrete near-linear path between the end points
               CALL build_path(ngrid, path, idx, ndx, npath, distmm, rgrid) 
               ! refine the path aiming for a "maximum probability path" 
               IF (npath > 2) THEN  ! does a few iterations for path refinement
                  CALL neb_path(ngrid, path, npath, distmm, probnmm, 0.1d0)
                  CALL neb_path(ngrid, path, npath, distmm, probnmm, 0.1d0)               
                  CALL neb_path(ngrid, path, npath, distmm, probnmm, 0.1d0)
                  CALL neb_path(ngrid, path, npath, distmm, probnmm, 0.1d0)
                  CALL neb_path(ngrid, path, npath, distmm, probnmm, 0.1d0)
                  CALL neb_path(ngrid, path, npath, distmm, probnmm, 0.1d0)
                  CALL neb_path(ngrid, path, npath, distmm, probnmm, 0.1d0)
               ENDIF            
               
               ! check if the path goes downhill to within accuracy
               DO i=2,npath-1
   !!!!!check using the error
   !!!!! p1-p2 -> associated error, err = SQRT(errp1**2+errp2**2) 
   !!!!! relative error -> err/(p1-p2)
   !!!!! that is the comparison I have to do!           
                  IF(nbootstrap>0)THEN
                     relerr=DSQRT(errors(path(i))**2+errors(idx)**2)/ &
                              (probnmm(path(i))-probnmm(idx))
                     IF(ABS(relerr)<kderr) THEN
                         qs_next = idx 
                         EXIT
                     ENDIF
                  ELSE
                     IF ((probnmm(path(i))-probnmm(idx))/ &
                        (probnmm(path(i))+probnmm(idx))<kderr) THEN 
                         qs_next = idx 
                         EXIT
                     ENDIF
                  ENDIF
                  !IF ((probnmm(path(i))-probnmm(idx))/ &
                  !    (probnmm(path(i))+probnmm(idx))<-3*kderr) THEN 
                  !   ! yey! we found a point lower in probability so we should not jump!
                  !   qs_next = idx       
                  !   !IF (verbose) THEN
                  !   !   dmin = 0.0d0
                  !   !   write(*,*) "# SADDLE POINT DETECTED"
                  !   !   write(*,*) 1, dmin, probnmm(path(1)), path(1)
                  !   !   DO j=2,npath                        
                  !   !      dmin = dmin + distmm(path(j),path(j-1))
                  !   !      write(*,*) j, dmin, probnmm(path(j)), path(j)                        
                  !   !   ENDDO
                  !   !ENDIF
                  !   EXIT
                  !ENDIF
               ENDDO
               ! check if the path contains crazy jumps
               DO i=1,npath-1
                  IF (dsqrt( distmm(path(i),path(i+1)) ) >  &
                     6*(dsqrt(rgrid(path(i)))+dsqrt(rgrid(path(i+1)))) ) THEN
                     qs_next=idx ! abort jump!
                     IF (verbose) THEN
                        dmin = 0.0d0
                        write(*,*) "# LONG JUMP DETECTED"
                        write(*,*) 1, dmin, probnmm(path(1)), path(1)
                        DO j=2,npath                        
                           dmin = dmin + distmm(path(j),path(j-1))
                           write(*,*) j, dmin, probnmm(path(j)), path(j)                        
                        ENDDO
                     ENDIF
                     EXIT
                  ENDIF
               ENDDO     
            ENDIF         
         ELSE
            qs_next=idx
            DO j=1,ngrid
               IF(probnmm(j)>probnmm(idx))THEN
                  IF(nbootstrap>0)THEN     
                    ! ok, chek the error associated
                    relerr=DSQRT(errors(j)**2+errors(idx)**2)/(probnmm(j)-probnmm(idx))
                    IF(relerr>kderr) CONTINUE
                  ENDIF
                  IF((distmm(idx,j).LT.dmin) .AND. (distmm(idx,j).LT.lambda))THEN
                     dmin=distmm(idx,j)
                     qs_next=j
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         
      END FUNCTION qs_next

      DOUBLE PRECISION FUNCTION fkernel(D,period,sig2,vc,vp)
            ! Calculate the (normalized) gaussian kernel
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

      DOUBLE PRECISION FUNCTION genkernel(ngrid,D,period,sig2,probnmm,vp,vgrid)
            ! Calculate the (normalized) gaussian kernel
            ! in an arbitrary point
            !
            ! Args:
            !    ngrid: number of grid point
            !    D: dimensionality
            !    period: periodicity
            !    sig2: sig**2
            !    probnmm: probability of the grid points
            !    vp: point's vector
            !    vgrid: grid points

            INTEGER, INTENT(IN) :: D
            INTEGER, INTENT(IN) :: ngrid
            DOUBLE PRECISION, INTENT(IN) :: period(D)
            DOUBLE PRECISION, INTENT(IN) :: sig2
            DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: probnmm
            DOUBLE PRECISION, INTENT(IN) :: vp(D)
            DOUBLE PRECISION, INTENT(IN) :: vgrid(D,ngrid)
            
            INTEGER j
            DOUBLE PRECISION res!,norm
            
            res=0.0d0
            !norm=0.0d0    
            DO j=1,ngrid
               res=res+(probnmm(j)/( (twopi*sig2)**(dble(D)/2) ))* &
                   dexp(-0.5d0*pammr2(D,period,vgrid(:,j),vp)/sig2)
               !norm=norm+probnmm(j)
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

   END PROGRAM pamm
