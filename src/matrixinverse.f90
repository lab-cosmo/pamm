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
         
         SUBROUTINE inv2x2(M,IM)
            ! Invert the matrix, simplest case
            ! ...
            ! Args:
            !    param: descript 
            DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN) :: M
            DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: IM
            
            DOUBLE PRECISION det
            
            IM=0.0d0
            det=M(1,1)*M(2,2) - M(1,2)*M(2,1)
            IM(1,1) = M(2,2)/det
            IM(1,2) = -M(1,2)/det
            IM(2,1) = -M(2,1)/det
            IM(2,2) = M(1,1)/det

         END SUBROUTINE inv2x2

         SUBROUTINE inv3x3(M,IM)
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

         END SUBROUTINE inv3x3

      END MODULE matrixinverse
