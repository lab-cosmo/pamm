! The main program which runs our driver test case potentials
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
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!

      PROGRAM gmm
         USE matrixinverse
         USE gaussian
         USE gaussmix
      IMPLICIT NONE
      
      CHARACTER*70 :: filename, gaussianfile, outputfile
      CHARACTER*1024 :: cmdbuffer
      DOUBLE PRECISION smooth,errc,dummyd1,loglike,oldloglike,test,bic,aic
      INTEGER delta,dummyi1,Nlines,NSamples
      ! gaussians
      INTEGER Nk
      ! array of gaussians
      TYPE(GAUSS_TYPE), ALLOCATABLE, DIMENSION(:) :: clusters
      DOUBLE PRECISION, DIMENSION(3,2) :: rangevwd
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vwad,pnk
      INTEGER ccmd
      INTEGER i,j,k,counter
    
      LOGICAL verbose

      !! DEFAULT
      filename="NULL"
      gaussianfile="NULL"
      outputfile="gaussiansout.dat"
      ccmd=0
      delta=1
      smooth=1.0e-5
      verbose = .false.
      Nk=-1
      ! Initialize loglike
      loglike=-1e-40
      rangevwd=0.0d0
      !!

      !!!!! PARSER
      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-i") THEN ! data file (v,w,rad)
            ccmd = 1
         ELSEIF (cmdbuffer == "-gf") THEN ! file containing gaussian parmeters
            ccmd = 2
         ELSEIF (cmdbuffer == "-o") THEN ! output file
            ccmd = 3
         ELSEIF (cmdbuffer == "-ev") THEN ! delta
            ccmd = 4
         ELSEIF (cmdbuffer == "-s") THEN ! smoothing factor
            ccmd = 5
         ELSEIF (cmdbuffer == "-n") THEN ! Nk
            ccmd = 6
         ELSEIF (cmdbuffer == "-err") THEN ! error
            ccmd = 7
         ELSEIF (cmdbuffer == "-h") THEN ! help
            WRITE(*,*) ""
            WRITE(*,*) " GMM (Gaussian-Mixture Model)"
            CALL helpmessage
            CALL EXIT(-1)
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = .true.
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) ""
               WRITE(*,*) " Wrong usage. Insert the right parameters!"
               CALL helpmessage
               CALL EXIT(-1)
            ELSEIF (ccmd == 1) THEN ! input file
               filename=trim(cmdbuffer)
            ELSEIF (ccmd == 2) THEN ! gaussian file
               gaussianfile=trim(cmdbuffer)
            ELSEIF (ccmd == 3) THEN ! output file
               outputfile=trim(cmdbuffer)
            ELSEIF (ccmd == 4) THEN ! delta
               READ(cmdbuffer,*) delta
            ELSEIF (ccmd == 5) THEN ! delta
               READ(cmdbuffer,*) smooth
            ELSEIF (ccmd == 6) THEN ! delta
               READ(cmdbuffer,*) Nk   
            ELSEIF (ccmd == 7) THEN ! delta
               READ(cmdbuffer,*) errc
            ENDIF
         ENDIF
      ENDDO
      !!!!! END PARSER
      
      Nlines = GetNlines(filename)
      NSamples = Nlines/delta
      counter=0
      !open (11, file=filename)
      !do while (.true.)
      !   read (11, *, end=999) x
      !enddo
      !999 continue
      ALLOCATE(vwad(3,NSamples))
      
      OPEN(UNIT=12,FILE=filename,STATUS='OLD',ACTION='READ')
      DO i=1,Nlines
         IF(MODULO(i,delta)==0)THEN
            counter=counter+1
            READ(12,*) vwad(1,counter),vwad(2,counter),vwad(3,counter)
            ! define the max,min values
            DO j=1,3
               ! min
               IF(vwad(j,counter).LT.rangevwd(j,1)) rangevwd(j,1)=vwad(j,counter)
               ! max
               IF(vwad(j,counter).GT.rangevwd(j,2)) rangevwd(j,2)=vwad(j,counter)
            ENDDO
         ENDIF
      ENDDO
      CLOSE(UNIT=12)
      
      ! Mandatory parameters
      IF (filename.EQ."NULL") THEN
         WRITE(*,*) ""
         WRITE(*,*) " Error: insert an input file! "
         CALL helpmessage
         CALL EXIT(-1)
      ENDIF
      
      IF (gaussianfile.EQ."NULL") THEN
         ! start from skratch
         IF(Nk==-1) THEN
            WRITE(*,*) "Insert a gaussian number or an initial file."
            CALL EXIT(-1)
         ENDIF
         ALLOCATE(clusters(Nk))
         CALL generatefromscratch(Nk,1,rangevwd,clusters)
      ELSE
         OPEN(UNIT=12,FILE=gaussianfile,STATUS='OLD',ACTION='READ')
         !decomment this if the first line is a comment string
         READ(12,*) cmdbuffer
         READ(12,*) cmdbuffer
         READ(12,*) dummyi1
         IF ((Nk.EQ.dummyi1).OR.(Nk.LT.dummyi1)) THEN
            ALLOCATE(clusters(dummyi1))
            DO i=1,Nk
               CALL readgaussfromfile(12,clusters(i))
            ENDDO
            IF (Nk.LT.dummyi1) THEN
               Nk=dummyi1
               IF(verbose) WRITE(*,*) "Your gaussians number is smaller than the one from the init file."
            ENDIF
         ELSE ! Nk is greater than that in the file
            IF(verbose)THEN
               WRITE(*,*) "Your gaussians number is bigger than the one from the init file."
               WRITE(*,*) "Number of gaussian started from scratch: ", (Nk-dummyi1)
            ENDIF
            ALLOCATE(clusters(Nk))
            DO i=1,dummyi1
               CALL readgaussfromfile(12,clusters(i))
            ENDDO
            ! generate from scratch the others
            CALL generatefromscratch(Nk,(dummyi1+1),rangevwd,clusters)
         ENDIF
         
         CLOSE(UNIT=12)
            
      ENDIF
      
      ALLOCATE(pnk(Nk,NSamples))
      
      i=0
      test=1e10
      DO WHILE(test<errc)
         i=i+1
         IF(verbose) WRITE(*,*) "Iteration N. ",i
         ! EXPECTATION STEP
         oldloglike=loglike
         CALL estep(Nk,NSamples,vwad,clusters,pnk,loglike)
         test=abs((loglike-oldloglike)/loglike)
         ! MAXIMIZATION STEP
         CALL mstep(Nk,NSamples,vwad,clusters,pnk,smooth)
         ! BIC e AIC
         bic=-2*loglike+(Nk*3+Nk*9+Nk)*dlog(dfloat(NSamples))
         aic=-2*loglike+2*(Nk*3+Nk*9+Nk)
         IF(verbose)THEN
            WRITE(*,*) " LogLikelihood: ", loglike
            WRITE(*,*) " Error: ", test
            WRITE(*,*) " BIC: ", bic
            WRITE(*,*) " AIC: ", aic
         ENDIF
      ENDDO
      
      OPEN(UNIT=11,FILE=outputfile,STATUS='REPLACE',ACTION='WRITE')
      WRITE(11,*) "# "
      WRITE(11,*) "# "
      WRITE(11,*) Nk
      DO k=1,Nk
         WRITE(11,*) clusters(k)%mean(1)," ",clusters(k)%mean(2)," ",clusters(k)%mean(3)," ", &
                     clusters(k)%cov(1,1)," ",clusters(k)%cov(1,2)," ",clusters(k)%cov(1,3)," ", &
                     clusters(k)%cov(2,1)," ",clusters(k)%cov(2,2)," ",clusters(k)%cov(2,3)," ", &
                     clusters(k)%cov(3,1)," ",clusters(k)%cov(3,2)," ",clusters(k)%cov(3,3)," ", &
                     dlog(clusters(k)%pk)
      ENDDO
      CLOSE(UNIT=11)
      
      DEALLOCATE(clusters)
      DEALLOCATE(vwad)
      DEALLOCATE(pnk)

      CONTAINS

         SUBROUTINE helpmessage
            WRITE(*,*) ""
            WRITE(*,*) " SYNTAX: gmm [-h] -i filename [-gf gaussianfile] [-o outputfile] "
            WRITE(*,*) "             [-n gaussians number] [-err error] [-ev delta] [-s smoothing_factor] [-v] "
            WRITE(*,*) ""
         END SUBROUTINE helpmessage
         
         SUBROUTINE readgaussfromfile(fileid,gaussp)
            ! Calculate the probabilities
            ! ...
            ! Args:
            !    param: descript
            INTEGER, INTENT(IN) :: fileid
            TYPE(GAUSS_TYPE) , INTENT(INOUT) :: gaussp
            
            DOUBLE PRECISION dummy
            
            READ(fileid,*) gaussp%mean(1), gaussp%mean(2), gaussp%mean(3), &
                           gaussp%cov(1,1), gaussp%cov(2,1), gaussp%cov(3,1), &
                           gaussp%cov(1,2), gaussp%cov(2,2), gaussp%cov(3,2), &
                           gaussp%cov(1,3), gaussp%cov(2,3), gaussp%cov(3,3), &
                           dummy
            gaussp%pk=exp(dummy) 
            ! calculate once the Icovs matrix and the norm_const
            CALL gauss_prepare(gaussp)
            
         END SUBROUTINE readgaussfromfile
         
         SUBROUTINE generatefromscratch(Ng,iniz,rangevwd,clusters) ! to improve
            ! Calculate the probabilities
            ! ...
            ! Args:
            !    param: descript
            INTEGER, INTENT(IN) :: Ng
            INTEGER, INTENT(IN) :: iniz
            DOUBLE PRECISION, DIMENSION(3,2), INTENT(IN) :: rangevwd
            TYPE(GAUSS_TYPE), DIMENSION(Nk), INTENT(INOUT) :: clusters
            INTEGER i,j
            
            DO i=iniz,Ng
               clusters(i)%cov=0.0d0
               DO j=1,3
                  ! random mean
                  clusters(i)%mean(j)=RAND()*(rangevwd(j,2)-rangevwd(j,1))+rangevwd(j,1) 
                  ! spherical gaussian
                  clusters(i)%cov(j,j)=1e-2
               ENDDO
               ! arbitrary Pk
               clusters(i)%pk=1e-1
            ENDDO
            
         END SUBROUTINE generatefromscratch
         
         INTEGER FUNCTION GetNlines(filen)
            CHARACTER*70, INTENT(IN) :: filen
            
            CHARACTER*1024 :: cmdbuffer
			
			   cmdbuffer='wc -l '//filen//'>> tmp.tmp'
            CALL system(cmdbuffer)
            OPEN(UNIT=8,FILE='tmp.tmp',STATUS='OLD',ACTION='READ')
            ! read the line numbers
            READ(8,*) Nlines
            CLOSE(UNIT=8)
            !remove the temp file
            cmdbuffer="rm tmp.tmp"
            CALL system(cmdbuffer)
         END FUNCTION

      END PROGRAM gmm
