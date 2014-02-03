! This contains the algorithms needed to calculate 2D gaussian.
!
! Copyright (C) 2014, Piero Gasparotto and Michele Ceriotti
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!
! Functions:
!    X_functions: Description


      MODULE gaussian
         USE matrixinverse
      IMPLICIT NONE
      
      DOUBLE PRECISION, PARAMETER :: dpigreco = (2.0d0*3.14159265358979d0)**3

      TYPE gauss_type
         DOUBLE PRECISION pk
         DOUBLE PRECISION norm
         DOUBLE PRECISION, DIMENSION(3) :: mean
         DOUBLE PRECISION, DIMENSION(3,3) :: cov
         DOUBLE PRECISION, DIMENSION(3,3) :: icov
      END TYPE

      CONTAINS

         SUBROUTINE gauss_prepare(gpars)
            ! Initialize the gaussian calculating...
            ! ...
            ! Args:
            !    param: descript 
            TYPE(GAUSS_TYPE), INTENT(INOUT) :: gpars
            DOUBLE PRECISION det
            
            
            det = gpars%cov(1,1)*(gpars%cov(2,2)*gpars%cov(3,3)-gpars%cov(3,2)*gpars%cov(2,3)) - &
                  gpars%cov(1,2)*(gpars%cov(2,1)*gpars%cov(3,3)-gpars%cov(2,3)*gpars%cov(3,1)) + &
                  gpars%cov(1,3)*(gpars%cov(2,1)*gpars%cov(3,2)-gpars%cov(2,2)*gpars%cov(3,1))
                  
            CALL inv3x3(gpars%cov,gpars%icov)

            gpars%norm = gpars%pk/dsqrt(dpigreco*det)   ! includes the weight in
                                                        ! the normalization constant
         END SUBROUTINE gauss_prepare

         DOUBLE PRECISION FUNCTION gauss_eval(gpars, x)
            ! Get the value value of the gaussian in X
            TYPE(GAUSS_TYPE), INTENT(IN) :: gpars
            DOUBLE PRECISION, INTENT(IN) :: x(3)
            DOUBLE PRECISION dv(3),tmpv(3)
            DOUBLE PRECISION xcx

            dv=x-gpars%mean
            
            tmpv = matmul(dv,gpars%icov)
            
            xcx = dot_product(dv,tmpv)

            gauss_eval = gpars%norm * dexp(-0.5d0*xcx)
         END FUNCTION gauss_eval

      END MODULE
