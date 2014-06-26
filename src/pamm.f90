! This file contain the main program for the PAMM clustering.
! Starting from a set of data points in high dimension it will first perform
! a non-parametric partitioning of the probability density and return the
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
!    readgaussfromfile: Get the gaussian paramters from file
!    helpmessage: Banner containing the help message
!    ordergaussians: Order gaussians

      PROGRAM pamm
      USE libpamm
      IMPLICIT NONE

      CHARACTER(LEN=1024) :: outputfile                       ! The output file prefix
      DOUBLE PRECISION, ALLOCATABLE :: xref(:)                ! Reference point needed to find out the important cluster
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: distmm ! similarity matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: diff     ! temp vector used to store distances

      INTEGER D                                               ! Dimensionality of problem
      INTEGER Nk                                              ! Number of gaussians in the mixture
      INTEGER nsamples                                        ! Total number points
      INTEGER ngrid                                         ! Number of samples extracted using minmax

      INTEGER jmax,ii,jj
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma2, wj, probnmm, msmu
      DOUBLE PRECISION :: normwj                              ! accumulator for wj
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npvoronoi, iminij, pnlist, nlist
      INTEGER seed                                            ! seed for the random number generator

      ! variable to set the covariance matrix
      DOUBLE PRECISION tmppks,normpks

      ! Array of Gaussians containing the gaussians parameters
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters
      ! Array containing the input data pints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x , y
      ! quick shift, roots and path to reach the root (used to speedup the calculation)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: idxroot, qspath

      ! PARSER
      CHARACTER(LEN=1024) :: cmdbuffer, comment   ! String used for reading text lines from files
      INTEGER ccmd                  ! Index used to control the input parameters
      INTEGER nmsopt   ! number of mean-shift optimizations of the cluster centers
      LOGICAL verbose  ! flag for verbosity
      LOGICAL weighted ! flag for using weigheted data
      INTEGER isep1, isep2, par_count  ! temporary indices for parsing command line arguments
      DOUBLE PRECISION lambda, lambda2, msw

      INTEGER i,j,k,counter,dummyi1 ! Counters and dummy variable


!!!!!!! Default value of the parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      outputfile="out"
      ccmd=0              ! no parameters specified
      Nk=0                ! number of gaussians
      nmsopt=0            ! number of mean-shift refinements
      ngrid=-1          ! number of samples extracted with minmax
      seed=12345          ! seed for the random number generator
      lambda=-1              ! quick shift cut-off
      verbose = .false.   ! no verbosity
      weighted= .false.   ! don't use the weights
      D=-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-o") THEN           ! output file
            ccmd = 2
         ELSEIF (cmdbuffer == "-d") THEN       ! dimensionality
            ccmd = 9
         ELSEIF (cmdbuffer == "-seed") THEN    ! seed for the random number genarator
            ccmd = 4
         ELSEIF (cmdbuffer == "-l") THEN       ! threshold to differentiate different clusters
            ccmd = 5
         ELSEIF (cmdbuffer == "-ngrid") THEN   ! N of grid points
            ccmd = 7
         ELSEIF (cmdbuffer == "-nms") THEN     ! N of mean-shift steps
            ccmd = 6
         ELSEIF (cmdbuffer == "-ref") THEN     ! point from wich calculate the distances to order
            ccmd = 8                           ! the gaussians in the output
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
            ELSEIF (ccmd == 2) THEN ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 9) THEN
               READ(cmdbuffer,*) D
               ALLOCATE(xref(D))
               xref=0.0d0
               xref(1)=-1.0d0
            ELSEIF (ccmd == 4) THEN ! read the seed
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 6) THEN ! read the seed
               READ(cmdbuffer,*) nmsopt
            ELSEIF (ccmd == 5) THEN ! cut-off for quick-shift
               READ(cmdbuffer,*) lambda
            ELSEIF (ccmd == 7) THEN ! number of grid points
               READ(cmdbuffer,*) ngrid
            ELSEIF (ccmd == 8) THEN
               IF (D<0) STOP "Dimensionality (-d) must be precede the reference point (-red). "
               par_count = 1
               isep1 = 0
               DO WHILE (index(cmdbuffer(isep1+1:), ',') > 0)
                  isep2 = index(cmdbuffer(isep1+1:), ',') + isep1
                  READ(cmdbuffer(isep1+1:isep2-1),*) xref(par_count)
                  par_count = par_count + 1
                  isep1=isep2
               ENDDO
               READ(cmdbuffer(isep1+1:),*) xref(par_count)
            ENDIF
         ENDIF
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL SRAND(seed) ! initialize the random number generator

      ! Check the input parameters
      IF (ccmd.EQ.0) THEN
         ! ccmd==0 : no parameters inserted
         WRITE(*,*) ""
         WRITE(*,*) " Could not interpret command line arguments. Please check for errors!"
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

      ! get the data from standard input
      CALL readinput(D, weighted, nsamples, x, normwj, wj)

      ! If not specified, the number voronoi polyhedras
      ! are set to the square of the total number of points
      IF (ngrid.EQ.-1) ngrid=int(sqrt(float(nsamples)))

      ALLOCATE(iminij(nsamples))
      ALLOCATE(pnlist(ngrid+1),nlist(nsamples))
      ALLOCATE(y(D,ngrid),npvoronoi(ngrid),probnmm(ngrid),sigma2(ngrid))
      ALLOCATE(idxroot(ngrid),qspath(ngrid),distmm(ngrid,ngrid))
      ALLOCATE(diff(D), msmu(D))

      ! Extract ngrid points on which the kernel density estimation is to be
      ! evaluated. Also partitions the nsamples points into the Voronoi polyhedra
      ! of the sampling points.
      IF(verbose) THEN
         WRITE(*,*) "NSamples: ", nsamples
         WRITE(*,*) "Selecting ", ngrid, " points using MINMAX"
      ENDIF

      CALL mkgrid(D,nsamples,ngrid,x,y,npvoronoi,iminij)

      ! Generate the neighbour list
      IF(verbose) write(*,*) "Generating neighbour list"
      CALL getnlist(nsamples,ngrid,npvoronoi,iminij, pnlist,nlist)

      ! Definition of the similarity matrix between grid points
      distmm=0.0d0
      sigma2=1e10 ! adaptive sigma of the kernel density estimator
      IF(verbose) write(*,*) "Computing similarity matrix"
      DO i=1,ngrid
         DO j=1,i-1
            ! distance between two voronoi centers
            distmm(i,j) = dot_product( y(:,i) - y(:,j) , y(:,i) - y(:,j) )
            if (distmm(i,j) < sigma2(i)) sigma2(i) = distmm(i,j)
            ! the symmetrical one
            distmm(j,i) = distmm(i,j)
            if (distmm(i,j) < sigma2(j)) sigma2(j) = distmm(i,j)
         ENDDO
      ENDDO

      IF(lambda.EQ.-1)THEN
         ! set automatically the mean shift lambda set to 5*<sig>
         lambda2=SUM(sigma2)/ngrid
         lambda=5.0d0*dsqrt(lambda2)
      ENDIF
      lambda2=lambda*lambda ! we always work with squared distances....

      IF(verbose) write(*,*) "Computing kernel density on reference points."
      IF(verbose) write(*,*) "Lambda : ", lambda

      ! computes the KDE on the Voronoi centers using the neighbour list
      probnmm = 0.0d0
      DO i=1,ngrid
         DO j=1,ngrid

            ! do not compute KDEs for points that belong to far away Voronoj
            IF (distmm(i,j)/sigma2(j)>36.0d0) CYCLE

            ! cycle just inside the polyhedra thanx to the neighbour list
            DO k=pnlist(j)+1,pnlist(j+1)
               probnmm(i)=probnmm(i)+ wj(nlist(k))* &
                          fkernel(D,sigma2(j),y(:,i),x(:,nlist(k)))
            ENDDO
         ENDDO
         probnmm(i)=probnmm(i)/normwj
      ENDDO

      IF(verbose) write(*,*) "Running quick shift"

      idxroot=0
      ! Start quick shift
      DO i=1,ngrid
         IF(idxroot(i).NE.0) CYCLE
         !idxroot(i)=i
         qspath=0
         qspath(1)=i
         counter=1
         DO WHILE(qspath(counter).NE.idxroot(qspath(counter)))
            idxroot(qspath(counter))= &
               qs_next(ngrid,qspath(counter),lambda2,probnmm,distmm)
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
          DO j=1,D
            WRITE(11,"((A1,ES15.4E4))",ADVANCE = "NO") " ", y(j,i)
          ENDDO
          WRITE(11,"(A1,I4,A1,ES15.4E4)") " ", dummyi1 , " ", probnmm(i)
          ! accumulate the normalization factor for the pks
          normpks=normpks+probnmm(i)
      ENDDO
      CLOSE(UNIT=11)

      ! now qspath contains the indexes of Nk gaussians

      !! covariance
      ALLOCATE(clusters(Nk))
      ! calculate the gaussians covariance from the data in the clusters
      DO k=1,Nk
         ALLOCATE(clusters(k)%mean(D))
         ALLOCATE(clusters(k)%cov(D,D))
         ALLOCATE(clusters(k)%icov(D,D))

         clusters(k)%mean=y(:,qspath(k))

         ! optionally do a few mean-shift steps to find a better estimate of the cluster mode
         DO j=1,nmsopt
            msmu=0.0d0
            tmppks=0.0d0
            DO i=1,ngrid
               msw =  probnmm(i)*exp(-0.5*sum((y(:,i)-clusters(k)%mean)*(y(:,i)-clusters(k)%mean))/(lambda2/25))
               msmu = msmu + msw*y(:,i)
               tmppks = tmppks + msw
            ENDDO
            clusters(k)%mean = msmu / tmppks
         ENDDO

         clusters(k)%cov = 0.0d0
         tmppks=0.0d0
         DO i=1,ngrid
            IF(idxroot(i).NE.qspath(k)) CYCLE
            tmppks=tmppks+probnmm(i)
            DO ii=1,D
               DO jj=1,D
                  clusters(k)%cov(ii,jj)=clusters(k)%cov(ii,jj)+probnmm(i)* &
                       (y(ii,i)-clusters(k)%mean(ii))*(y(jj,i)-clusters(k)%mean(jj))
               ENDDO
            ENDDO
         ENDDO
         clusters(k)%cov=clusters(k)%cov/tmppks
         clusters(k)%weight=tmppks/normpks
         clusters(k)%D=D
      ENDDO
      ! write gaussians

      ! oreder gaussians from the closest to the reference point
      CALL sortclusters(Nk, clusters, xref)


      ! write a 2-lines header
      WRITE(comment,*) "# PAMM clusters analysis. NSamples: ", nsamples, " NGrid: ", &
                ngrid, " QSLambda: ", lambda, ACHAR(10), "# Dimensionality/NClusters//Pk/Mean/Covariance "

      OPEN(UNIT=12,FILE=trim(outputfile)//".pamm",STATUS='REPLACE',ACTION='WRITE')

      CALL writeclusters(12, comment, nk, clusters)
      CLOSE(UNIT=12)
      DEALLOCATE(clusters)
      DEALLOCATE(x,wj)
      DEALLOCATE(idxroot,qspath,distmm)
      DEALLOCATE(pnlist,nlist,iminij)
      DEALLOCATE(y,npvoronoi,probnmm,sigma2)
      DEALLOCATE(diff)

      CALL EXIT(0)
      ! end of the main programs



         ! functions and subroutines !!!!!!!!!!!!!!!!!!!!

      CONTAINS

      SUBROUTINE helpmessage
         ! Banner to print out for helping purpose
         !

         WRITE(*,*) ""
         WRITE(*,*) " USAGE: pamm [-h] -d D [-w] [-o output] [-ngrid ngrid] [-l lambda] "
         WRITE(*,*) "              [-seed seedrandom] [-rif -1,0,0,...] [-v] "
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
         WRITE(*,*) "   -ngrid ngrid      : Number of grid points to evaluate KDE [sqrt(nsamples)]"
         WRITE(*,*) "   -nms nms          : Do nms mean-shift steps with a Gaussian width lambda/5 to"
         WRITE(*,*) "                       optimize cluster centers [0] "
         WRITE(*,*) "   -seed seed        : Seed to initialize the random number generator. [12345]"
         WRITE(*,*) "   -ref X            : Reference point for ordering the clusters [ (-1,0,0,...) ]"
         WRITE(*,*) "   -v                : Verbose output "
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
               wbuff(counter+1)=1.0d0
            ENDIF
            totw=totw+wbuff(counter+1)
            IF(io_status<0 .or. io_status==5008) EXIT    ! also intercepts a weird error given by some compilers when reading past of EOF
            IF(io_status>0) STOP "*** Error occurred while reading file. ***"

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

      SUBROUTINE mkgrid(D,nsamples,ngrid,x,y,npvoronoi,iminij)
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
               dij = dot_product( y(:,i-1) - x(:,j) , &
                                  y(:,i-1) - x(:,j) )
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
            dij = dot_product( y(:,ngrid) - x(:,j) , &
                               y(:,ngrid) - x(:,j) )
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

      INTEGER FUNCTION qs_next(ngrid,idx,lambda,probnmm,distmm)
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
         DOUBLE PRECISION, INTENT(IN) :: lambda
         DOUBLE PRECISION, DIMENSION(ngrid), INTENT(IN) :: probnmm
         DOUBLE PRECISION, DIMENSION(ngrid,ngrid), INTENT(IN) :: distmm

         INTEGER j
         DOUBLE PRECISION dmin

         dmin=1.0d10
         qs_next=idx
         DO j=1,ngrid
            IF(probnmm(j)>probnmm(idx))THEN
               IF((distmm(idx,j).LT.dmin) .AND. (distmm(idx,j).LT.lambda))THEN
                  dmin=distmm(idx,j)
                  qs_next=j
               ENDIF
            ENDIF
         ENDDO
      END FUNCTION qs_next

      DOUBLE PRECISION FUNCTION fkernel(D,sig2,vc,vp)
            ! Calculate the (normalized) gaussian kernel
            !
            ! Args:
            !    sig2: sig**2
            !    vc: voronoi center's vector
            !    vp: point's vector

            INTEGER, INTENT(IN) :: D
            DOUBLE PRECISION, INTENT(IN) :: sig2
            DOUBLE PRECISION, INTENT(IN) :: vc(D)
            DOUBLE PRECISION, INTENT(IN) :: vp(D)

                                       ! put here ** -D/2
            fkernel=(1/( (twopi*sig2)**(dble(D)/2) ))* &
                    dexp(-sum((vc-vp)*(vc-vp))*0.5/sig2)
      END FUNCTION fkernel


      SUBROUTINE sortclusters(nk,clusters,prif)
         ! Sort the gaussians from the closest to prif
         ! Bubble-sort ordering is implemented here
         !
         ! Args:
         !    nk: number of gaussian clusters
         !    clusters: array containing gaussians parameters
         !    prif: reference point

         INTEGER, INTENT(IN) :: nk
         TYPE(gauss_type), DIMENSION(nk), INTENT(INOUT) :: clusters
         DOUBLE PRECISION, DIMENSION(D), INTENT(IN) :: prif

         TYPE(gauss_type) tmpgauss
         DOUBLE PRECISION distances(nk),tmpdistance
         INTEGER j,i
         LOGICAL :: swapped = .TRUE.

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

   END PROGRAM pamm
