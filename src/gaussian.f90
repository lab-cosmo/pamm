! This contains the routines needed to define and estimate a 3D gaussian.
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
! Types:
!    gauss_type: Structure defining a gaussian
!
! Functions:
!    gauss_prepare: Initialize all the parameters of the gaussian
!    gauss_logeval: Return the logarithm of the multivariate gaussian density
!    gauss_eval: Return the multivariate gaussian density


      MODULE gaussian
         USE matrixinverse
      IMPLICIT NONE
      
      DOUBLE PRECISION, PARAMETER :: dpigreco = (2.0d0*3.14159265358979d0)
      
      ! Structure that contains the parameters needed to define and
      ! estimate a gaussian
      TYPE gauss_type
         DOUBLE PRECISION lnorm ! logarithm of the normalization factor
         DOUBLE PRECISION det ! determinant of the covariance matrix
         DOUBLE PRECISION, DIMENSION(3) :: mean
         DOUBLE PRECISION, DIMENSION(3,3) :: cov ! convariance matrix
         DOUBLE PRECISION, DIMENSION(3,3) :: icov ! inverse convariance matrix
      END TYPE

      CONTAINS

         SUBROUTINE gauss_prepare(gpars)
            ! Initialize all the parameters of the gaussian
            ! 
            ! Args:
            !    gpars: gauss_type variable to initialize
             
            TYPE(gauss_type), INTENT(INOUT) :: gpars
            
            ! calculate the determinant of the covariance matrix
            gpars%det = gpars%cov(1,1)*(gpars%cov(2,2)*gpars%cov(3,3)-gpars%cov(3,2)*gpars%cov(2,3)) - &
                        gpars%cov(1,2)*(gpars%cov(2,1)*gpars%cov(3,3)-gpars%cov(2,3)*gpars%cov(3,1)) + &
                        gpars%cov(1,3)*(gpars%cov(2,1)*gpars%cov(3,2)-gpars%cov(2,2)*gpars%cov(3,1))
            
            ! calculate the inverse of the convariance matrix      
            CALL inv3x3(gpars%cov,gpars%icov)
            
            ! calculate the  logarithm of the normalization factor
            gpars%lnorm = dlog(1.0d0/dsqrt((dpigreco**3)*gpars%det))
         END SUBROUTINE gauss_prepare

         DOUBLE PRECISION FUNCTION gauss_logeval(gpars, x)
            ! Return the logarithm of the multivariate gaussian density
            ! 
            ! Args:
            !    gpars: gaussian parameters
            !    x: point in wich estimate the log of the gaussian 
            
            TYPE(gauss_type), INTENT(IN) :: gpars
            DOUBLE PRECISION, INTENT(IN) :: x(3)
            DOUBLE PRECISION dv(3),tmpv(3)
            DOUBLE PRECISION xcx

            dv=x-gpars%mean
            tmpv = matmul(dv,gpars%icov)
            xcx = dot_product(dv,tmpv)

            gauss_logeval = gpars%lnorm - 0.5d0*xcx
         END FUNCTION gauss_logeval

         DOUBLE PRECISION FUNCTION gauss_eval(gpars, x)
            ! Return the multivariate gaussian density
            ! Args:
            !    gpars: gaussian parameters
            !    x: point in wich estimate the value of the gaussian 
            
            TYPE(gauss_type), INTENT(IN) :: gpars
            DOUBLE PRECISION, INTENT(IN) :: x(3)
            
            gauss_eval = dexp(gauss_logeval(gpars,x))
         END FUNCTION gauss_eval

      END MODULE
