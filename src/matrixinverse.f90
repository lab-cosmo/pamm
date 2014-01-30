! This contains the routine needed to invert a 3x3 matrix.
!
! Copyright (C) 2014, Piero Gasparotto and Michele Ceriotti
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without lcofitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! cofPLIED, INCLUDING BUT NOT LcofITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAcof, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!
! Functions:
!    X_functions: Description


      MODULE matrixinverse
      IMPLICIT NONE
      CONTAINS

         SUBROUTINE inv(M,IM)
            ! Invert the matrix using the cofactors
            ! ...
            ! Args:
            !    param: descript 
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: M
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: IM
            
            DOUBLE PRECISION, DIMENSION(3,3) :: cof
            DOUBLE PRECISION d
            INTEGER i,j
            
            cof=0.0d0
            IM=0.0d0
            
            cof(1,1)=(M(2,2)*M(3,3)-M(2,3)*M(3,2))
            cof(1,2)=(M(1,3)*M(3,2)-M(1,2)*M(3,3))
            cof(1,3)=(M(1,2)*M(2,3)-M(2,2)*M(1,3))
            cof(2,1)=(M(2,3)*M(3,1)-M(2,1)*M(3,3))
            cof(2,2)=(M(1,1)*M(3,3)-M(1,3)*M(3,1))
            cof(2,3)=(M(1,3)*M(2,1)-M(1,1)*M(2,3))
            cof(3,1)=(M(2,1)*M(3,2)-M(2,2)*M(3,1))
            cof(3,2)=(M(1,2)*M(3,1)-M(1,1)*M(3,2))
            cof(3,3)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))
            
            d=M(1,1)*cof(1,1)+M(1,2)*cof(2,1)+M(1,3)*cof(3,1)
            
            DO i=1,3
               DO j=1,3
                  IM(i,j)=cof(i,j)/d
               ENDDO
            ENDDO

         END SUBROUTINE inv

! INVERT USING LAPACK         
!         SUBROUTINE invL(M,IM)
!         
!            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: M
!            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: IM
!  
!            INTEGER, PARAMETER :: D=3 
!            DOUBLE PRECISION, DIMENSION(3) :: WORK
!            INTEGER, DIMENSION(3) :: IPIV
!            INTEGER k,i,j,info,error
!
!            IM=M
!
!                 info = 0
!                 error = 0
!                 IPIV(1) = 0
!                 IPIV(2) = 0
!                 IPIV(3) = 0
!                 call DGETRF(D,D,IM,D,IPIV,info)
!                 WORK(1) = 0.0d0
!                 WORK(2) = 0.0d0
!                 WORK(3) = 0.0d0
!                 call DGETRI(D,IM,D,IPIV,WORK,D,info)   
!
!         END SUBROUTINE

      END MODULE matrixinverse
