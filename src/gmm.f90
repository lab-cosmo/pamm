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

      PROGRAM gmm
         USE gaussmix
      IMPLICIT NONE

      CHARACTER*256 :: filename                               ! The input data file containing v,w and rad
      CHARACTER*256 :: outputfile                             ! The output file containing the calculated gaussians
      DOUBLE PRECISION :: prif(3)                             ! Reference point needed to find out the important gaussian
      DOUBLE PRECISION twosig2                                ! 2*sig^2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: distij ! similarity matrix
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: means               ! temporary vector to stor the centers of the gaussians
      INTEGER, ALLOCATABLE, DIMENSION(:) :: meansidx          ! vector to store the points corrispondig the gaussians centers
      DOUBLE PRECISION, DIMENSION(3) :: diff                  ! temp vector used to store distances
      
      LOGICAL goclust                                         ! variable used to chek if a cluster i new or not during quick shift
      INTEGER npc                                             ! number of point in a cluster. Used to define the gaussians covariances
      INTEGER Nk                                              ! Number of gaussians in the mixture
      INTEGER delta                                           ! Number of data to skeep (for a faster calculation)
      INTEGER Nlines                                          ! Number of lines of the input data file
      INTEGER nsamples                                        ! Total number points
      INTEGER nminmax                                         ! Number of samples extracted using minmax
      DOUBLE PRECISION :: dij, dmax, dmin
      INTEGER imin,jmax
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dminij, weights, probnmm
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iminij
      INTEGER seed     ! seed for the random number generator
      LOGICAL kernelmode

      ! Array of Gaussians containing the gaussians parameters
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters
      ! Array containing the logarithm of the fractions (Pk) for each gaussian
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: lpks
      ! Array containing the input data pints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vwad , Ymm
      INTEGER, ALLOCATABLE, DIMENSION(:) :: YmmIclust
      ! Array containing the nk probalities for each of nsamples points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pnk ! responsibility matrix

      ! PARSER
      CHARACTER*1024 :: cmdbuffer   ! String used for reading text lines from files
      INTEGER ccmd                  ! Index used to control the input parameters
      LOGICAL verbose ! flag for verbosity
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string
	  DOUBLE PRECISION vpar(5)
	  DOUBLE PRECISION dummyr1,tau,meannew(3)
	  LOGICAL testtree
	  INTEGER idxtemp,idx
	  
      INTEGER i,j,k,counter,dummyi1 ! Counters and dummy variable
      
      
      !!!!!!! Iniatialze the parameters !!!!!!!
      filename="NULL"
      outputfile="out"
      ccmd=0 ! no parameters specified
      delta=1 ! read every point
      Nk=0 ! number of gaussians
      nsamples=-1         ! total number of points
      nminmax=-1          ! number of samples extracted with minmax
      seed=1357           ! seed for the random number generator
      tau=0.7d0           ! quick shift cut-off
      prif(1)=-1.0d0      ! point from wich order the gaussians in the output
      prif(2)= 2.9d0
      prif(3)= 2.9d0
      kernelmode=.false.  ! no gaussian kernel during the density estimation
      twosig2=-1.1d0      ! 2*sig^2 (for the gaussian kernel)
      verbose = .false.   ! no verbosity
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-i") THEN          ! data file (v,w,rad)
            ccmd = 1
         ELSEIF (cmdbuffer == "-o") THEN      ! output file
            ccmd = 2
         ELSEIF (cmdbuffer == "-ev") THEN     ! delta
            ccmd = 3
         ELSEIF (cmdbuffer == "-seed") THEN   ! seed for the random number genarator
            ccmd = 4
         ELSEIF (cmdbuffer == "-tau") THEN    ! threshold to differentiate different clusters
            ccmd = 5
         ELSEIF (cmdbuffer == "-nsamples") THEN ! N samples (total)
            ccmd = 6
         ELSEIF (cmdbuffer == "-nminmax") THEN ! N samples extracted with minmax 
            ccmd = 7
         ELSEIF (cmdbuffer == "-rif") THEN    ! point from wich calculate the distances to order
            ccmd = 8                          ! the gaussians in the output
         ELSEIF (cmdbuffer == "-kernel") THEN ! Apply a kernel to defines the points' probality
            kernelmode = .true.
         ELSEIF (cmdbuffer == "-v") THEN      ! verbosity flag
            verbose = .true.
         ELSEIF (cmdbuffer == "-h") THEN      ! help flag
            WRITE(*,*) ""
            WRITE(*,*) " Quick Shift Gaussian-Mixture "
            CALL helpmessage
            CALL EXIT(-1)
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " Wrong usage. Insert the right parameters!"
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 1) THEN ! input data file
               filename=trim(cmdbuffer)
            ELSEIF (ccmd == 2) THEN ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN ! delta
               READ(cmdbuffer,*) delta
            ELSEIF (ccmd == 4) THEN ! read the seed
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 5) THEN ! cut-off for quick-shift
               READ(cmdbuffer,*) tau
            ELSEIF (ccmd == 6) THEN
               READ(cmdbuffer,*) nsamples
            ELSEIF (ccmd == 7) THEN
               READ(cmdbuffer,*) nminmax
            ELSEIF (ccmd == 8) THEN ! vrif,wrif,dADrif
               par_count = 1
               commas(1) = 0
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
      IF (ccmd.EQ.0 .OR. filename.EQ."NULL") THEN
         ! ccmd==0 : no parameters inserted
         ! filename=="NULL" no input datas specified
         WRITE(*,*) ""
         WRITE(*,*) " Wrong usage. Insert the right parameters!"
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! Get total number of data points in the file
      Nlines = GetNlines(filename)
      
      ! if the user don't want to use the total number of points
      ! in the file the -flag nsamples can be used
      IF(nsamples.NE.-1) delta=Nlines/nsamples    
      ! adjust nsamples
      nsamples = Nlines/delta
      IF(verbose) WRITE(*,*) "NSamples: " , nsamples
      
      ! Read the points from the input file    
      counter=0
      ALLOCATE(vwad(3,nsamples))
      OPEN(UNIT=12,FILE=filename,STATUS='OLD',ACTION='READ')
      DO i=1,Nlines
         IF(MODULO(i,delta)==0)THEN ! read a point every delta
            counter=counter+1       ! index for the real data set
            READ(12,*) vwad(1,counter),vwad(2,counter),vwad(3,counter)
         ELSE
            READ(12,*) dummyr1      ! discard the line
         ENDIF
      ENDDO
      CLOSE(UNIT=12)       
     
      ! If not specified, the number samples (voronoi polyhedras)
      ! are set to the square of the total number of points
      IF (nminmax.EQ.-1) nminmax=sqrt(float(nsamples))
      
      ALLOCATE(dminij(nsamples),iminij(nsamples))
      ALLOCATE(Ymm(3,nminmax),YmmIclust(nminmax),weights(nminmax))
      ALLOCATE(distij(nminmax,nminmax),probnmm(nminmax))
      ALLOCATE(means(3,nminmax),meansidx(nminmax))
      
      ! Extract the samples (centers of the voronoi polyhedras)
      ! Start MINMAX
      IF(verbose) write(*,*) "Selecting ", nminmax, " points using MINMAX"
      ! choose randomly the first point 
      Ymm(:,1)=vwad(:,int(RAND()*nsamples))
      dminij = 1.0d99
      iminij = 1
      DO i=2,nminmax
         dmax = 0.0d0
         DO j=1,nsamples
            dij = dot_product( Ymm(:,i-1) - vwad(:,j) , Ymm(:,i-1) - vwad(:,j) )
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
         IF(verbose .AND. (modulo(i,1000).EQ.0)) write(*,*) i,"/",nminmax
      ENDDO
      
      ! finishes Voronoi attribution
      DO j=1,nsamples
         dij = dot_product( Ymm(:,nminmax) - vwad(:,j) , Ymm(:,nminmax) - vwad(:,j) )
         IF  (dminij(j)>dij) THEN
            dminij(j) = dij
            iminij(j) = nminmax
         ENDIF
      ENDDO
      ! Voronoi weights and variances
      weights=0.0d0
      DO j=1,nsamples
         weights(iminij(j))=weights(iminij(j))+1 
      ENDDO
      
      ! Density estimation for quick shift
      probnmm=weights/nsamples
      
      ! Define tau from the data ... TODO
      
      dmax=0.0d0 
      twosig2=0.0d0
      DO i=1,nminmax-1
         dminij(i) = 1.0d10
         DO j=2,nminmax
            dij = dot_product( Ymm(:,i) - Ymm(:,j) , Ymm(:,i) - Ymm(:,j) )
            distij(i,j) = dij
            distij(j,i) = dij
            IF (dij<dminij(i)) dminij(i)=dij                            
         ENDDO
         IF (dmax < dminij(i)) dmax=dminij(i)
         twosig2=twosig2+dminij(i)
      ENDDO
      ! set sig=dNN/2
      twosig2=0.5d0*twosig2/nminmax
      
      write(*,*) "Avarage dNN: ", dsqrt(twosig2*2)
	
	  OPEN(UNIT=11,FILE="c-"//trim(outputfile)//".dat",STATUS='REPLACE',ACTION='WRITE')
	  
	  means=0.0d0
	  ! Start quick shift
	  DO i=1,nminmax
	     YmmIclust(i)=i
	     idx=0
	     DO WHILE(idx.NE.YmmIclust(i))
	     ! serch the closest neighbour inside tau with higher Pj
	        idx=YmmIclust(i)
	        dmin=1.0d10
	        DO j=1,nminmax
	           IF(idx.EQ.j)CYCLE
	           IF(probnmm(j)>probnmm(idx))THEN
	              IF((distij(idx,j).LT.dmin) .AND. (distij(idx,j).LT.tau))THEN
	                 dmin=distij(idx,j)
	                 YmmIclust(i)=j
	              ENDIF
	           ENDIF 
	        ENDDO
	     ENDDO
	     
	     ! check the gaussians
	     goclust=.false.
	     DO k=1,Nk
	        IF (meansidx(k)==YmmIclust(i)) THEN
	           goclust=.true.
	           YmmIclust(i)=k
	           EXIT
	        ENDIF
	     ENDDO
	     IF(.NOT.(goclust))THEN
	        ! A new gaussian has been found
	        Nk=Nk+1
	        meansidx(Nk)=YmmIclust(i)
	        means(:,Nk)=Ymm(:,YmmIclust(i))
	        YmmIclust(i)=Nk
	        WRITE(*,*) "Cluster ", Nk , means(:,Nk)
	     ENDIF
	     
	     ! write out the clusters
	     WRITE(11,"(3(A1,ES15.4E4))",ADVANCE = "NO") " ", Ymm(1,i)," ", Ymm(2,i)," ", Ymm(3,i)
	     WRITE(11,"(A1,I4,A1,ES15.4E4)") " ", YmmIclust(i) , " ", probnmm(i)
	  ENDDO
	  
	  CLOSE(UNIT=11)
	  
	  !! covariance
	  ALLOCATE(clusters(Nk),lpks(Nk))
	  ! calculate the gaussians covariance from the data in the clusters
	  DO k=1,Nk
	     meannew = 0.0d0
	     npc=0
	     clusters(k)%mean=means(:,k)
	     clusters(k)%cov = 0.0d0
	     DO i=1,nminmax
	        IF(YmmIclust(i).NE.k) CYCLE
	        meannew = meannew + weights(i)*Ymm(:,i)
	        clusters(k)%cov(1,1) = clusters(k)%cov(1,1) + weights(i)*Ymm(1,i) * Ymm(1,i)
	        clusters(k)%cov(1,2) = clusters(k)%cov(1,2) + weights(i)*Ymm(1,i) * Ymm(2,i)
	        clusters(k)%cov(1,3) = clusters(k)%cov(1,3) + weights(i)*Ymm(1,i) * Ymm(3,i)
	        clusters(k)%cov(2,2) = clusters(k)%cov(2,2) + weights(i)*Ymm(2,i) * Ymm(2,i)
	        clusters(k)%cov(2,3) = clusters(k)%cov(2,3) + weights(i)*Ymm(2,i) * Ymm(3,i)
	        clusters(k)%cov(3,3) = clusters(k)%cov(3,3) + weights(i)*Ymm(3,i) * Ymm(3,i)
	        npc=npc+weights(i)
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
	     lpks(k)=LOG(FLOAT(npc)/nsamples)
	  ENDDO
	  ! write gaussians
	  
	  ! oreder gaussians from the closest to the reference point
	  CALL ordergaussians(Nk,clusters,lpks,prif)
	  
	  OPEN(UNIT=11,FILE="g-"//trim(outputfile)//".dat",STATUS='REPLACE',ACTION='WRITE')
	  !WRITE(11,*) "# Mean-Shift output (Sig,Err Conv, Err Clusters): " , dsqrt(twosig2/2.0d0), errc, errclusters
	  WRITE(11,"(A31)",ADVANCE="NO") "# Quick Shift GM output. Ntot: "
	  WRITE(11,"(I12,A12,I11)") nsamples," , NVoroni: ",nminmax
	  WRITE(11,*) "# mean cov pk"
	  WRITE(11,*) Nk
	  DO k=1,Nk
	     CALL writegausstofile(11,clusters(k),lpks(k))
	  ENDDO
	  
	  CLOSE(UNIT=11)
	  DEALLOCATE(vwad,YmmIclust,clusters,lpks,Ymm,dminij,distij,probnmm)
	  DEALLOCATE(means,meansidx)
	  CALL EXIT(0)
	  ! end of the main programs
 
 
 
      ! functions and subroutines !!!!!!!!!!!!!!!!!!!!

      CONTAINS

         SUBROUTINE helpmessage
            ! Banner to print out for helping purpose
            !

            WRITE(*,*) ""
            WRITE(*,*) " USAGE: gmm [-h] -i filename [-o output] [-seed seedrandom] "
            WRITE(*,*) "            [-tau tau] [-ev delta] [-nsamples NTot] [-nminmax Nminmax] "
            WRITE(*,*) "            [-rif v,w,R] [-v] "
            WRITE(*,*) ""
            WRITE(*,*) " Clusterize the data and define the mixture of gaussians from that. "
            WRITE(*,*) " It is mandatory to specify the file containg the data. "
            WRITE(*,*) " For the other options a default is defined.  "
            WRITE(*,*) ""
            WRITE(*,*) "   -h                : Print this message "
            WRITE(*,*) "   -i inputfile      : File containin the input data "
            WRITE(*,*) "   -o output         : Name for the output. This will produce : "
            WRITE(*,*) ""
            WRITE(*,*) "                            c-output.dat (clusterized data) "
            WRITE(*,*) "                            g-output.dat (gaussians) "
            WRITE(*,*) ""
            WRITE(*,*) "   -seed seed        : Seed to initialize the random number generator "
            WRITE(*,*) "   -tau tau          : Quick shift cutoff "
            WRITE(*,*) "   -ev delta         : Stride reading data frome the file "
            WRITE(*,*) "   -nsamples ntot    : Number of points to read from the input file "
            WRITE(*,*) "   -nminmax nvoroni  : Number of Voronoi polyhedra "
            WRITE(*,*) "   -rif v,w,R        : Point from wich order the gaussians "
            WRITE(*,*) "   -v                : Verobose mode "
            WRITE(*,*) ""
         END SUBROUTINE helpmessage

         SUBROUTINE readgaussfromfile(fileid,gaussp,lpk)
            ! Read a line from the file and get the paramters for the related gaussian
            !
            ! Args:
            !    fileid: the file containing the gaussians parameters
            !    gaussp: type_gaussian container in wich we store the gaussian parameters
            !    lpk: logarithm of the Pk associated to the gaussian

            INTEGER, INTENT(IN) :: fileid
            TYPE(gauss_type) , INTENT(INOUT) :: gaussp
            DOUBLE PRECISION, INTENT(INOUT) :: lpk

            READ(fileid,*) gaussp%mean(1), gaussp%mean(2), gaussp%mean(3), &
                           gaussp%cov(1,1), gaussp%cov(2,1), gaussp%cov(3,1), &
                           gaussp%cov(1,2), gaussp%cov(2,2), gaussp%cov(3,2), &
                           gaussp%cov(1,3), gaussp%cov(2,3), gaussp%cov(3,3), &
                           lpk
            lpk=log(lpk)
            CALL gauss_prepare(gaussp)

         END SUBROUTINE readgaussfromfile
         
         SUBROUTINE writegausstofile(fileid,gaussp,lpk)
            ! Read a line from the file and get the paramters for the related gaussian
            !
            ! Args:
            !    fileid: the file containing the gaussians parameters
            !    gaussp: type_gaussian container in wich we store the gaussian parameters
            !    lpk: logarithm of the Pk associated to the gaussian

            INTEGER, INTENT(IN) :: fileid
            TYPE(gauss_type) , INTENT(INOUT) :: gaussp
            DOUBLE PRECISION, INTENT(INOUT) :: lpk
            
            WRITE(fileid,*) gaussp%mean(1)," ",gaussp%mean(2)," ",gaussp%mean(3)," ", &
                            ! write covariance matrix
                            gaussp%cov(1,1)," ",gaussp%cov(1,2)," ",gaussp%cov(1,3)," ", &
                            gaussp%cov(2,1)," ",gaussp%cov(2,2)," ",gaussp%cov(2,3)," ", &
                            gaussp%cov(3,1)," ",gaussp%cov(3,2)," ",gaussp%cov(3,3)," ", &
                            ! write Pk of the gaussian
                            dexp(lpk)

         END SUBROUTINE writegausstofile

         SUBROUTINE generatefromscratch(maxmin,nsamples,vwda,ng,clusters,lpks)
            ! Iniatialize from scratch. Guess starting values for
            ! the gaussians parameters.
            !
            ! Args:
            !    maxmin: logical to decide to if maxmin or not
            !    ng: number of gaussians to generate
            !    nsamples: number of points
            !    vwda: array containing the data
            !    rangevwd: array containing the extreme values in the input data set
            !    pnk: responsibility matrix
            !    loglike: logarithm of the likelihood
            !    clusters: array containing gaussians parameters
            !    lpks: array containing the logarithm of the Pk for each gaussian

            LOGICAL, INTENT(IN) :: maxmin
            INTEGER, INTENT(IN) :: ng
            INTEGER, INTENT(IN) :: nsamples
            DOUBLE PRECISION, DIMENSION(3,nsamples), INTENT(IN) :: vwda
            TYPE(gauss_type), DIMENSION(ng), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(ng), INTENT(INOUT) :: lpks

            INTEGER i,j,jmax
            DOUBLE PRECISION :: cov(3,3), m(3)
            ! maxmin
            DOUBLE PRECISION :: dminij(nsamples), dij, dmax

            ! gets initial covariance as the covariance of the whole data set
            cov = 0.0
            m = 0.0
            do i=1,nsamples
              m = m + vwda(:,i)
              cov(1,1) = cov(1,1) + vwda(1,i) * vwda(1,i)
              cov(1,2) = cov(1,2) + vwda(1,i) * vwda(2,i)
              cov(1,3) = cov(1,3) + vwda(1,i) * vwda(3,i)
              cov(2,2) = cov(2,2) + vwda(2,i) * vwda(2,i)
              cov(2,3) = cov(2,3) + vwda(2,i) * vwda(3,i)
              cov(3,3) = cov(3,3) + vwda(3,i) * vwda(3,i)
            enddo
            m = m/nsamples
            cov = cov/nsamples
            cov(1,1) = cov(1,1) - m(1)*m(1)
            cov(1,2) = cov(1,2) - m(1)*m(2)
            cov(1,3) = cov(1,3) - m(1)*m(3)
            cov(2,2) = cov(2,2) - m(2)*m(2)
            cov(2,3) = cov(2,3) - m(2)*m(3)
            cov(3,3) = cov(3,3) - m(3)*m(3)
            cov(2,1) = cov(1,2)
            cov(3,1) = cov(1,3)
            cov(3,2) = cov(2,3)

            ! initializes the means randomly
            DO i=1,ng
              clusters(i)%mean = vwda(:,int(RAND()*nsamples))
            ENDDO

            IF(maxmin)THEN
               ! initializes the means with minmax
               clusters(1)%mean = vwda(:,int(RAND()*nsamples))
               dminij = 1d99
               DO i=2,ng
                  dmax = 0.0
                  DO j=1,nsamples
                     dij = dot_product( clusters(i-1)%mean - vwda(:,j) , &
                                        clusters(i-1)%mean - vwda(:,j) )
                     dminij(j) = min(dminij(j), dij)
                     IF (dminij(j) > dmax) THEN
                        dmax = dminij(j)
                        jmax = j
                     ENDIF
                  ENDDO
                  clusters(i)%mean = vwda(:, jmax)
               ENDDO
            ENDIF

            ! initializes the covariance matrix of all the clusters
            DO i=1,ng
               clusters(i)%cov = cov
               lpks(i)=dlog(1.0d0/ng)
               CALL gauss_prepare(clusters(i))
            ENDDO
         END SUBROUTINE generatefromscratch

         INTEGER FUNCTION GetNlines(filein)
            ! Get the number of lines from the file.
            ! This is usefull to calculate the number of timesnapshots.
            !
            ! Args:
            !    filein: a text file

            CHARACTER*70, INTENT(IN) :: filein

            CHARACTER*1024 :: cmdbuffer
            ! call the shell command wc and save the output to a temporary file
			   cmdbuffer='wc -l '//filein//'>> tmp.tmp'
			   ! execute the bash command
            CALL system(cmdbuffer)
            ! read the result from the temporary file
            OPEN(UNIT=8,FILE='tmp.tmp',STATUS='OLD',ACTION='READ')
            READ(8,*) GetNlines
            CLOSE(UNIT=8)
            !remove the temporary files
            cmdbuffer="rm tmp.tmp"
            CALL system(cmdbuffer)
         END FUNCTION
         
         SUBROUTINE ordergaussians(ng,clusters,pks,prif)
            ! Order the gaussians from closest to prif
            !
            ! Args:
            !    ng: number of gaussians to generate
            !    clusters: array containing gaussians parameters
            !    prif: reference point
            
            INTEGER, INTENT(IN) :: ng
            TYPE(gauss_type), DIMENSION(ng), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(ng), INTENT(INOUT) :: pks
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: prif
            
            TYPE(gauss_type) tmpgauss
            DOUBLE PRECISION distances(ng),tmpdistance,tmppk
			INTEGER j,i
			LOGICAL :: swapped = .TRUE.
			
			! calculate the distances
			DO i=1,ng
			   distances(i)=dot_product(clusters(i)%mean, prif)		   
            ENDDO
            ! now we can sort using the distances
            ! will use bubble sort
            DO j=ng-1,1,-1
               swapped = .FALSE.
               DO i = 1, j
                  IF (distances(i) > distances(i+1)) THEN            
                     tmpdistance=distances(i)
                     distances(i)=distances(i+1)
                     distances(i+1)=tmpdistance
                     ! upgrade pks
                     tmppk=pks(i)
                     pks(i)=pks(i+1)
                     pks(i+1)=tmppk
                     ! upgrade also the clusters
                     tmpgauss=clusters(i)
                     clusters(i)=clusters(i+1)
                     clusters(i+1)=tmpgauss
                     swapped = .TRUE.
                  END IF
               END DO
               IF (.NOT. swapped) EXIT	   
            ENDDO         
         END SUBROUTINE ordergaussians

      END PROGRAM gmm
