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

      CHARACTER*1024 :: filename                               ! The input data file containing v,w and rad
      CHARACTER*1024 :: outputfile                             ! The output file containing the calculated gaussians
      DOUBLE PRECISION :: prif(3)                             ! Reference point needed to find out the important gaussian
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: distmm ! similarity matrix
      DOUBLE PRECISION, DIMENSION(3) :: diff                  ! temp vector used to store distances
      
      INTEGER npc                                             ! number of point in a cluster. Used to define the gaussians covariances
      INTEGER Nk                                              ! Number of gaussians in the mixture
      INTEGER delta                                           ! Number of data to skeep (for a faster calculation)
      INTEGER Nlines                                          ! Number of lines of the input data file
      INTEGER nsamples                                        ! Total number points
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
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string
	  DOUBLE PRECISION vpar(5)
	  DOUBLE PRECISION dummyr1,tau,meannew(3), tau2
	  
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
      tau=-1              ! quick shift cut-off
      prif(1)=-1.0d0      ! point from wich order the gaussians in the output
      prif(2)= 2.9d0
      prif(3)= 2.9d0
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
      ALLOCATE(vwad(3,nsamples),wj(nsamples))
      OPEN(UNIT=12,FILE=filename,STATUS='OLD',ACTION='READ')
      DO i=1,Nlines
         IF(MODULO(i,delta)==0)THEN ! read a point every delta
            counter=counter+1       ! index for the real data set
            READ(12,*) vwad(1,counter),vwad(2,counter),vwad(3,counter),wj(counter)
         ELSE
            READ(12,*) dummyr1      ! discard the line
         ENDIF
      ENDDO
      CLOSE(UNIT=12)       
     
      ! If not specified, the number samples (voronoi polyhedras)
      ! are set to the square of the total number of points
      IF (nminmax.EQ.-1) nminmax=sqrt(float(nsamples))
      
	  ALLOCATE(iminij(nsamples))
	  ALLOCATE(pnlist(nminmax+1),nlist(nsamples))
	  ALLOCATE(Ymm(3,nminmax),npvoronoi(nminmax),probnmm(nminmax),sigma2(nminmax))
	  ALLOCATE(idxroot(nminmax),qspath(nminmax),distmm(nminmax,nminmax))
      
      ! Extract nminmax points on which the kernel density estimation is to be
      ! evaluated. Also partitions the nsamples points into the Voronoi polyhedra
      ! of the sampling points.
      IF(verbose) write(*,*) "Selecting ", nminmax, " points using MINMAX"
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
	     meannew = 0.0d0
	     npc=0
	     clusters(k)%mean=Ymm(:,qspath(k))
	     clusters(k)%cov = 0.0d0
	     DO i=1,nminmax
	        IF(idxroot(i).NE.qspath(k)) CYCLE
	        meannew = meannew + probnmm(i)*Ymm(:,i)
	        clusters(k)%cov(1,1) = clusters(k)%cov(1,1) + probnmm(i)*Ymm(1,i) * Ymm(1,i)
	        clusters(k)%cov(1,2) = clusters(k)%cov(1,2) + probnmm(i)*Ymm(1,i) * Ymm(2,i)
	        clusters(k)%cov(1,3) = clusters(k)%cov(1,3) + probnmm(i)*Ymm(1,i) * Ymm(3,i)
	        clusters(k)%cov(2,2) = clusters(k)%cov(2,2) + probnmm(i)*Ymm(2,i) * Ymm(2,i)
	        clusters(k)%cov(2,3) = clusters(k)%cov(2,3) + probnmm(i)*Ymm(2,i) * Ymm(3,i)
	        clusters(k)%cov(3,3) = clusters(k)%cov(3,3) + probnmm(i)*Ymm(3,i) * Ymm(3,i)
	        npc=npc+1
	     ENDDO
	     meannew = meannew/npc
	     clusters(k)%cov = clusters(k)%cov /npc
	     clusters(k)%cov(1,1) = clusters(k)%cov(1,1) - probnmm(1)*probnmm(1)
	     clusters(k)%cov(1,2) = clusters(k)%cov(1,2) - probnmm(1)*probnmm(2)
	     clusters(k)%cov(1,3) = clusters(k)%cov(1,3) - probnmm(1)*probnmm(3)
	     clusters(k)%cov(2,2) = clusters(k)%cov(2,2) - probnmm(2)*probnmm(2)
	     clusters(k)%cov(2,3) = clusters(k)%cov(2,3) - probnmm(2)*probnmm(3)
	     clusters(k)%cov(3,3) = clusters(k)%cov(3,3) - probnmm(3)*probnmm(3)
	     clusters(k)%cov(2,1) = clusters(k)%cov(1,2)
	     clusters(k)%cov(3,1) = clusters(k)%cov(1,3)
	     clusters(k)%cov(3,2) = clusters(k)%cov(2,3)
	     lpks(k)=LOG(FLOAT(npc)/nminmax)
	  ENDDO
	  ! write gaussians
	  
	  ! oreder gaussians from the closest to the reference point
	  CALL ordergaussians(Nk,clusters,lpks,prif)
	  
	  OPEN(UNIT=11,FILE=trim(outputfile)//".gauss",STATUS='REPLACE',ACTION='WRITE')
	  !WRITE(11,*) "# Mean-Shift output (Sig,Err Conv, Err Clusters): " , dsqrt(twosig2/2.0d0), errc, errclusters
	  WRITE(11,"(A31)",ADVANCE="NO") "# Quick Shift GM output. Ntot: "
	  WRITE(11,"(I12,A12,I11)",ADVANCE="NO") nsamples," , NVoroni: ",nminmax
	  WRITE(11,"(A8,ES15.4E4)") " , Tau: ", tau
	  WRITE(11,*) "# mean cov pk"
	  WRITE(11,*) Nk
	  DO k=1,Nk
	     CALL writegausstofile(11,clusters(k),lpks(k))
	  ENDDO
	  
	  CLOSE(UNIT=11)
	  
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
            WRITE(*,*) " USAGE: getmodel [-h] -i filename [-o output] [-seed seedrandom] "
            WRITE(*,*) "                 [-tau tau] [-ev delta] [-nsamples NTot] [-nminmax Nminmax] "
            WRITE(*,*) "                 [-rif v,w,R] [-v] "
            WRITE(*,*) ""
            WRITE(*,*) " Clusterize the data and define a mixture of gaussians describing them. "
            WRITE(*,*) " It is mandatory to specify the file containg the data. "
            WRITE(*,*) " For the other options a default is defined.  "
            WRITE(*,*) ""
            WRITE(*,*) "   -h                : Print this message "
            WRITE(*,*) "   -i inputfile      : File containin the input data (XYZ format) "
            WRITE(*,*) "   -o output         : Name for the output. This will produce : "
            WRITE(*,*) ""
            WRITE(*,*) "                            output.clusters (clusterized data) "
            WRITE(*,*) "                            output.gauss (gaussians) "
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
            
            fkernel=dexp(-dot_product(vc-vp,vc-vp)*0.5/sig2)
         END FUNCTION fkernel

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

      END PROGRAM getmodel
