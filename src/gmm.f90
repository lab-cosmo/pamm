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
!    helpmessage: banner containing the help message

      PROGRAM gmm
         USE gaussmix
      IMPLICIT NONE

      CHARACTER*70 :: filename     ! The input data file containing v,w and rad
      CHARACTER*70 :: gaussianfile ! The input file containing the initial gaussians
      CHARACTER*70 :: outputfile   ! The output file containing the calculated gaussians
      DOUBLE PRECISION smooth      ! Smoothing Factor for the covariance matrix
      DOUBLE PRECISION errc        ! Convergence criterias
      DOUBLE PRECISION test        ! Test to check the convergence of the likelihood
      DOUBLE PRECISION loglike     ! Actual logLikelihood
      DOUBLE PRECISION oldloglike  ! Previous logLikelihood
      DOUBLE PRECISION bic         ! Bayesian Inference criteria
      DOUBLE PRECISION aic         ! Akaike Inference criteria
      DOUBLE PRECISION prif(3)     ! Reference point needed to find out the important gaussian
      INTEGER Nk       ! Number of gaussians in the mixture
      INTEGER delta    ! Number of data to skeep (for a faster calculation)
      INTEGER Nlines   ! Number of lines of the input data file
      INTEGER nsamples ! Number of points used during the calculation
      INTEGER seed     ! seed for the random number generator
      INTEGER best     ! index of the gaussian of interest
      LOGICAL maxmin ! flag for verbosity

      ! Array of Gaussians containing the gaussians parameters
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters
      ! Array containing the logarithm of the fractions (Pk) for each gaussian
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: lpks
      ! Array containing the input data pints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vwad
      ! Array containing the nk probalities for each of nsamples points
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pnk ! responsibility matrix

      ! PARSER
      CHARACTER*1024 :: cmdbuffer   ! String used for reading text lines from files
      INTEGER ccmd                  ! Index used to control the input parameters
      LOGICAL verbose ! flag for verbosity
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string

      INTEGER i,j,k,counter,dummyi1 ! Counters and dummy variable

      !!!!!!! Iniatialze the parameters !!!!!!!
      filename="NULL"
      gaussianfile="NULL"
      outputfile="gaussiansout.dat"
      ccmd=0 ! no parameters specified
      delta=1 ! read every point
      smooth=0.0d0 ! no smoothing by default
      verbose = .false. ! no verbosity
      maxmin = .false. ! no maxmin
      errc=1.0e-05
      Nk=-1 ! number of Gaussians in the mixture. Initialized at -1
      seed=1357
      test=111
      prif(1)=-1.0d0
      prif(2)= 2.9d0
      prif(3)= 2.9d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-i") THEN          ! data file (v,w,rad)
            ccmd = 1
         ELSEIF (cmdbuffer == "-gf") THEN     ! file containing intial gaussian parmeters
            ccmd = 2
         ELSEIF (cmdbuffer == "-o") THEN      ! output file
            ccmd = 3
         ELSEIF (cmdbuffer == "-ev") THEN     ! delta
            ccmd = 4
         ELSEIF (cmdbuffer == "-s") THEN      ! smoothing factor
            ccmd = 5
         ELSEIF (cmdbuffer == "-n") THEN      ! Nk
            ccmd = 6
         ELSEIF (cmdbuffer == "-err") THEN    ! convergence error in the likelihood
            ccmd = 7
         ELSEIF (cmdbuffer == "-seed") THEN   ! seed for the random number genarator
            ccmd = 8
         ELSEIF (cmdbuffer == "-rif") THEN    ! point from wich calculate the distances to order
            ccmd = 9                          ! the gaussians in the output
         ELSEIF (cmdbuffer == "-maxmin") THEN ! initialize using maxmin routine
            maxmin = .true.
         ELSEIF (cmdbuffer == "-h") THEN      ! help
            WRITE(*,*) ""
            WRITE(*,*) " GMM (Gaussian-Mixture Model)"
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (cmdbuffer == "-v") THEN      ! flag for verbose standard output
            verbose = .true.
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " Wrong usage. Insert the right parameters!"
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 1) THEN ! input data file
               filename=trim(cmdbuffer)
            ELSEIF (ccmd == 2) THEN ! initial gaussians parameters file
               gaussianfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 4) THEN ! delta
               READ(cmdbuffer,*) delta
            ELSEIF (ccmd == 5) THEN ! smoothing factor
               READ(cmdbuffer,*) smooth
            ELSEIF (ccmd == 6) THEN ! number of gaussians
               READ(cmdbuffer,*) Nk
            ELSEIF (ccmd == 7) THEN ! stop criteria
               READ(cmdbuffer,*) errc
            ELSEIF (ccmd == 8) THEN ! read the seed
               READ(cmdbuffer,*) seed
            ELSEIF (ccmd == 9) THEN ! vrif,wrif,dADrif
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
      IF (ccmd.EQ.0 .OR. Nk.LE.0 .OR. filename.EQ."NULL") THEN
         ! ccmd==0 : no parameters inserted
         ! Nk<0 : no Nk insterted
         ! filename=="NULL" no input datas specified
         WRITE(*,*) ""
         WRITE(*,*) " Wrong usage. Insert the right parameters!"
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF

      ! Get total number of data points in the file
      Nlines = GetNlines(filename)
      ! New number of points
      nsamples = Nlines/delta
      counter=0
      ALLOCATE(vwad(3,nsamples))
      OPEN(UNIT=12,FILE=filename,STATUS='OLD',ACTION='READ')
      DO i=1,Nlines
         IF(MODULO(i,delta)==0)THEN ! read a point every delta
            counter=counter+1 ! index for the real data set
            READ(12,*) vwad(1,counter),vwad(2,counter),vwad(3,counter)
         ENDIF
      ENDDO
      CLOSE(UNIT=12)

      ALLOCATE(pnk(Nk,nsamples))
      pnk=0.0d0 ! responsibility matrix initialized to zero

      ! Initializes Gaussians
      ALLOCATE(clusters(Nk), lpks(Nk)) ! Allocate the arrays with Nk specified by the user
      ! Generate from scractch the gaussians parameters
      CALL generatefromscratch(maxmin,nsamples,vwad,Nk,clusters,lpks)

      IF (gaussianfile.ne."NULL") THEN
         ! initializes Gaussians from file
         OPEN(UNIT=12,FILE=gaussianfile,STATUS='OLD',ACTION='READ')
         ! Skip the two header lines
         READ(12,*) cmdbuffer
         READ(12,*) cmdbuffer
         ! Read the number of gaussian from the file
         READ(12,*) dummyi1 ! use a dummy variables to check the equivalence with the input parameters
         ! Actually Nk=-1 or a value specified by the user with the command line
         IF(Nk.ne.dummyi1)THEN
            WRITE(0,*) "Number of Gaussians on command line does not match init.", &
                       "Will read what I can."
         ENDIF
         ! Nk=dummyi1  -> read all the gaussian from the file
         ! Nk<dummyi1  -> read only Nk from the file
         ! Nk>dummyi1  -> read all the gaussain from the file a the rest gaussians are
         !                generated from scratch
         DO i=1,min(Nk,dummyi1)
            CALL readgaussfromfile(12,clusters(i),lpks(i))
         ENDDO
         CLOSE(UNIT=12)
      ENDIF

      i=0
      ! Initialize loglike
      loglike=-1e-40
      DO WHILE(test>errc) ! convergence criteria
         i=i+1
         IF(verbose) WRITE(*,*) "Iteration N. ",i
         ! expectation step
         oldloglike=loglike
         CALL estep(Nk,nsamples,vwad,clusters,lpks,pnk,loglike)
         ! test the convergence
         test=abs((loglike-oldloglike)/loglike)
         ! maximization step
         CALL mstep(Nk,nsamples,vwad,clusters,lpks,pnk,smooth)
         ! Calculate BIC and AIC
         bic=-2*loglike+(10*Nk-1)*dlog(dfloat(nsamples))
         aic=-2*loglike+2*(10*Nk-1)
         IF(verbose)THEN
            WRITE(*,*) " LogLikelihood: ", loglike
            WRITE(*,*) " Error: ", test
            WRITE(*,*) " BIC: ", bic
            WRITE(*,*) " AIC: ", aic
         ENDIF
      ENDDO
      
      best=bestgaussian(Nk,clusters,prif)

      OPEN(UNIT=11,FILE=outputfile,STATUS='REPLACE',ACTION='WRITE')
      WRITE(11,"(A13,I11,A17,F16.7,A6,F16.7,A6,F16.7)") "# n.samples: ", nsamples, &
                                                         " loglikelihood: ", loglike/nsamples, " bic: ", &
                                                         bic/nsamples, " aic: ", aic/nsamples
      WRITE(11,*) "# mean cov pk"
      WRITE(11,*) Nk
      ! Write the first the most important gaussian
      CALL writegausstofile(11,clusters(best),lpks(best))
      ! and now all the others
      DO k=1,Nk      ! write mean
         IF(k.EQ.best) CYCLE
         CALL writegausstofile(11,clusters(k),lpks(k))
      ENDDO
      CLOSE(UNIT=11)

      DEALLOCATE(clusters,lpks,vwad,pnk)

      CONTAINS

         SUBROUTINE helpmessage
            ! Banner to print out for helping purpose
            !

            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: gmm [-h] [-v] -i filename -n gaussians number [-seed seedrandom] "
            WRITE(*,*) "             [-o outputfile] [-ev delta] [-err error] [-s smoothing_factor] "
            WRITE(*,*) "             [-gf gaussianfile] [-maxmin] [-rif vrif,wrif,dADrif]"
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
            cov(2,3) = cov(3,2)

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
         
         INTEGER FUNCTION bestgaussian(ng,clusters,prif)
            ! Return the index of the gaussian closest to prif
            !
            ! Args:
            !    ng: number of gaussians to generate
            !    clusters: array containing gaussians parameters
            !    prif: reference point
            
            INTEGER, INTENT(IN) :: ng
            TYPE(gauss_type), DIMENSION(ng), INTENT(INOUT) :: clusters
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: prif
            
            DOUBLE PRECISION old,new
			INTEGER i
			
			old=dsqrt(dot_product(clusters(1)%mean, prif))
			bestgaussian=1
			DO i=2,ng
			   new=dsqrt(dot_product(clusters(i)%mean, prif))
			   IF(new.LT.old)THEN
			      old=new
			      bestgaussian=i
			   ENDIF			   
            ENDDO

         END FUNCTION bestgaussian

      END PROGRAM gmm
