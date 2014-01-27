! This performs the calculations necessary to run a Lennard-Jones (LJ) 
! simulation.
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
! This contains the functions that calculate the potential, forces and
! virial tensor of a single-component LJ system.
! Includes functions which calculate the long-range correction terms for a
! simulation with a sharp nearest-neighbour cut-off.
!
! Functions:
!    X_functions: Description

      MODULE hbmixture
         USE distance
      IMPLICIT NONE

      INTEGER, PARAMETER :: TYPE_NONE=0 
      INTEGER, PARAMETER :: TYPE_H=1
      INTEGER, PARAMETER :: TYPE_ACCEPTOR=2
      INTEGER, PARAMETER :: TYPE_DONOR=4

      CONTAINS

         SUBROUTINE ninte
            ! Initialize the gaussian calculating...
            ! ...
            ! Args:
            !    param: descript 
            write(*,*) "ninte"

         END SUBROUTINE ninte

      END MODULE hbmixture
