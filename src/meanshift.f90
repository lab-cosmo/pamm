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

      CHARACTER*70 :: filename     ! The input data file containing v,w and rad
      CHARACTER*70 :: outputfile   ! The output file containing the clusterized data
      CHARACTER*70 :: gaussianfile ! The output file containing the calculated gaussians
      DOUBLE PRECISION err1        ! Convergence criterias for every point
      DOUBLE PRECISION err2        ! Convergence criterias for the clusters
      DOUBLE PRECISION sig
      DOUBLE PRECISION prif(3)     ! Reference point needed to find out the important gaussian
      INTEGER Nk       ! Number of clusters
      INTEGER delta    ! Number of data to skeep (for a faster calculation)
      INTEGER Nlines   ! Number of lines of the input data file
      INTEGER nsamples ! Number of points used during the calculation
      INTEGER best     ! index of the gaussian of interest

      ! Array of Gaussians containing the gaussians parameters
      DOUBLE PRECISION :: twosig2
      DOUBLE PRECISION, DIMENSION(3,50) :: means
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pks
      TYPE(gauss_type), ALLOCATABLE, DIMENSION(:) :: clusters
      ! Array containing the input data pints
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vwad
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vwadIclust
      INTEGER npc
      
      ! PARSER
      CHARACTER*1024 :: cmdbuffer   ! String used for reading text lines from files
      INTEGER ccmd                  ! Index used to control the input parameters
      LOGICAL goclust ! flag for verbosity
      INTEGER commas(4), par_count  ! stores the index of commas in the parameter string

      INTEGER i,j,k,counter,indexg ! Counters
      DOUBLE PRECISION kernel,norm,test
      DOUBLE PRECISION, DIMENSION(3) :: xold,xnew,diff
      
      !!!!!!! Iniatialze the parameters !!!!!!!
      filename="NULL"
      outputfile="NULL"
      gaussianfile="clusters.dat"
      ccmd=0 ! no parameters specified
      nsamples=1 ! read every point
      err1=1.0e-05
      err2=1.0e-05
      Nk=0 ! number of clusters
      sig=0.2
      prif(1)=-1.0d0
      prif(2)= 2.9d0
      prif(3)= 2.9d0
      nsamples=-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!! Command line parser !!!!!!!!!!!!!
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-i") THEN            ! data file (v,w,rad)
            ccmd = 1
         ELSEIF (cmdbuffer == "-o") THEN        ! output file
            ccmd = 2
         ELSEIF (cmdbuffer == "-gf") THEN        ! output file
            ccmd = 8
         ELSEIF (cmdbuffer == "-nsamples") THEN ! N samples
            ccmd = 3
         ELSEIF (cmdbuffer == "-sig") THEN      ! sigma (gaussian)
            ccmd = 4
         ELSEIF (cmdbuffer == "-err") THEN      ! cenvergence for the single point
            ccmd = 5
         ELSEIF (cmdbuffer == "-errc") THEN     ! cenvergence bitwin the point
            ccmd = 6
         ELSEIF (cmdbuffer == "-rif") THEN    ! point from wich calculate the distances to order
            ccmd = 7                          ! the gaussians in the output
         ELSEIF (cmdbuffer == "-h") THEN      ! help
            WRITE(*,*) ""
            WRITE(*,*) " Mean-shift algorithm"
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
            ELSEIF (ccmd == 8) THEN ! gaussians file
               gaussianfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN ! nsamples
               READ(cmdbuffer,*) nsamples
            ELSEIF (ccmd == 4) THEN ! smoothing factor
               READ(cmdbuffer,*) sig
            ELSEIF (ccmd == 5) THEN ! number of gaussians
               READ(cmdbuffer,*) err1
            ELSEIF (ccmd == 6) THEN ! stop criteria
               READ(cmdbuffer,*) err2
            ELSEIF (ccmd == 7) THEN ! vrif,wrif,dADrif
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
      ! New number of points
      IF(nsamples.EQ.-1)THEN
         nsamples = Nlines
         delta=1
      ELSE
         delta=INT(Nlines/nsamples)
         !write(*,*) nsamples,delta
      ENDIF
      counter=0
      ALLOCATE(vwad(3,nsamples),vwadIclust(nsamples))
      OPEN(UNIT=12,FILE=filename,STATUS='OLD',ACTION='READ')
      DO i=1,Nlines
         IF(MODULO(i,delta)==0)THEN ! read a point every delta
            counter=counter+1 ! index for the real data set
            READ(12,*) vwad(1,counter),vwad(2,counter),vwad(3,counter)
         ELSE
            READ(12,*) kernel
         ENDIF
      ENDDO
      CLOSE(UNIT=12)
      
      twosig2=2.0d0*sig*sig
      kernel=0.0d0
      vwadIclust=0
      IF(outputfile.NE."NULL") OPEN(UNIT=11,FILE=outputfile,STATUS='REPLACE',ACTION='WRITE')
      DO i=1,nsamples ! run over all points 1
         xold=vwad(:,i)
         test=1e10
         indexg=1
         norm=0                  ! iteration
         DO WHILE(test.GT.err1)  ! convergence criteria
            xnew=0.0d0
            norm=0.0d0
            DO j=1,nsamples      ! sum over all atoms
               diff = vwad(:,j)-xold
               kernel=dot_product(diff,diff)/twosig2
               if (kernel .gt. 10) cycle  ! skips tiny contributions
               kernel=dexp( -kernel )
               xnew=xnew+(kernel*vwad(:,j))
               norm=norm+kernel
            ENDDO                
            xnew=xnew/norm
            diff = xnew-xold
            test=dsqrt(dot_product(diff,diff))
            xold=xnew
            !write(*,*) ">> Test: ",test,xold
         ENDDO  ! close the iteration cycle    
         goclust=.false.
         DO k=1,Nk
            diff = xnew-means(:,k)
            
            IF (dsqrt(dot_product(diff,diff)).LT.err2) THEN
               goclust=.true.
               indexg=k
               exit
            ENDIF
         ENDDO
         IF(.NOT.(goclust))THEN
            Nk=Nk+1
            means(:,Nk)=xnew
            indexg=Nk
            WRITE(*,*) xnew(1),xnew(2),xnew(3)
         ENDIF
         IF(outputfile.NE."NULL")THEN 
            WRITE(11,"(3(A1,ES15.4E4))",ADVANCE = "NO") " ", vwad(1,i)," ", vwad(2,i)," ", vwad(3,i)
            WRITE(11,"(A1,I3)") " ", indexg
         ENDIF
         vwadIclust(i)=indexg
      ENDDO ! close the run over all points 1
      IF(outputfile.NE."NULL") CLOSE(UNIT=11)
      
      ALLOCATE(clusters(Nk),pks(Nk))
      ! calculate the gaussians covariance from the data in the clusters
      DO k=1,Nk 
         xnew = 0.0d0
         npc=0
         clusters(k)%mean=means(:,k)
         clusters(k)%cov = 0.0d0
         DO i=1,nsamples
            IF(vwadIclust(i).NE.k) CYCLE
            xnew = xnew + vwad(:,i)
            clusters(k)%cov(1,1) = clusters(k)%cov(1,1) + vwad(1,i) * vwad(1,i)
            clusters(k)%cov(1,2) = clusters(k)%cov(1,2) + vwad(1,i) * vwad(2,i)
            clusters(k)%cov(1,3) = clusters(k)%cov(1,3) + vwad(1,i) * vwad(3,i)
            clusters(k)%cov(2,2) = clusters(k)%cov(2,2) + vwad(2,i) * vwad(2,i)
            clusters(k)%cov(2,3) = clusters(k)%cov(2,3) + vwad(2,i) * vwad(3,i)
            clusters(k)%cov(3,3) = clusters(k)%cov(3,3) + vwad(3,i) * vwad(3,i)
            npc=npc+1
         ENDDO
         xnew = xnew/npc
         clusters(k)%cov = clusters(k)%cov /npc
         clusters(k)%cov(1,1) = clusters(k)%cov(1,1) - xnew(1)*xnew(1)
         clusters(k)%cov(1,2) = clusters(k)%cov(1,2) - xnew(1)*xnew(2)
         clusters(k)%cov(1,3) = clusters(k)%cov(1,3) - xnew(1)*xnew(3)
         clusters(k)%cov(2,2) = clusters(k)%cov(2,2) - xnew(2)*xnew(2)
         clusters(k)%cov(2,3) = clusters(k)%cov(2,3) - xnew(2)*xnew(3)
         clusters(k)%cov(3,3) = clusters(k)%cov(3,3) - xnew(3)*xnew(3)
         clusters(k)%cov(2,1) = clusters(k)%cov(1,2)
         clusters(k)%cov(3,1) = clusters(k)%cov(1,3)
         clusters(k)%cov(3,2) = clusters(k)%cov(2,3)
         pks(k)=npc/nsamples
      ENDDO
      
      ! oreder gaussians from the closest to the reference point
      CALL ordergaussians(Nk,clusters,pks,prif)

      OPEN(UNIT=11,FILE=gaussianfile,STATUS='REPLACE',ACTION='WRITE')
      !WRITE(11,"(A13,I11,A17,F16.7,A6,F16.7,A6,F16.7)") "# n.samples: ", nsamples, &
      !                                                   " loglikelihood: ", loglike/nsamples, " bic: ", &
      !                                                   bic/nsamples, " aic: ", aic/nsamples
      WRITE(11,*) "#"
      WRITE(11,*) "# mean cov pk"
      WRITE(11,*) Nk
      ! Write out the gaussians
      DO k=1,Nk      ! write mean
         CALL writegausstofile(11,clusters(k),pks(k))
      ENDDO
      !CLOSE(UNIT=11)

      DEALLOCATE(vwad,vwadIclust,clusters,pks)

      CONTAINS

         SUBROUTINE helpmessage
            ! Banner to print out for helping purpose
            !

            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: meanshift [-h] -i filename [-o outputfile] -nsamples N "
            WRITE(*,*) "                   [-sig gausssig] [-err error] [-errc error] "
            WRITE(*,*) "                   [-rif vrif,wrif,dADrif]"
            WRITE(*,*) ""
         END SUBROUTINE helpmessage
         
         SUBROUTINE writegausstofile(fileid,gaussp,pk)
            ! Read a line from the file and get the paramters for the related gaussian
            !
            ! Args:
            !    fileid: the file containing the gaussians parameters
            !    gaussp: type_gaussian container in wich we store the gaussian parameters
            !    pk: Pk associated to the gaussian

            INTEGER, INTENT(IN) :: fileid
            TYPE(gauss_type) , INTENT(INOUT) :: gaussp
            DOUBLE PRECISION, INTENT(INOUT) :: pk
            
            WRITE(fileid,*) gaussp%mean(1)," ",gaussp%mean(2)," ",gaussp%mean(3)," ", &
                            ! write covariance matrix
                            gaussp%cov(1,1)," ",gaussp%cov(1,2)," ",gaussp%cov(1,3)," ", &
                            gaussp%cov(2,1)," ",gaussp%cov(2,2)," ",gaussp%cov(2,3)," ", &
                            gaussp%cov(3,1)," ",gaussp%cov(3,2)," ",gaussp%cov(3,3)," ", &
                            ! write Pk of the gaussian
                            pk

         END SUBROUTINE writegausstofile

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
