module random
  implicit none
  contains
  
  subroutine random_init(seed)
    implicit none
    integer, intent(in) :: seed
    integer iseed,i
    integer, allocatable :: state(:)
      
    call random_seed(size=iseed) 
    allocate(state(iseed))
    do i=1,iseed
      state(i)=seed+i
    end do
    call random_seed(put=state)
    return
  end subroutine random_init
  
  double precision function random_uniform()
    call random_number(harvest=random_uniform)
  end function random_uniform
  
  double precision function random_gaussian()
    implicit none
    integer iset
    double precision gset,v1,v2,fac,rsq
    save iset,gset
    data iset/0/
    if (iset.eq.0) then
1      v1=2.d0*random_uniform()-1.d0
       v2=2.d0*random_uniform()-1.d0
       rsq=v1**2+v2**2  
       if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1 
       fac=sqrt(-2.d0*log(rsq)/rsq) 
       gset=v1 * fac
       random_gaussian=v2*fac
       iset=1
    else 
       random_gaussian=gset
       iset=0 
    endif
    return
  end function random_gaussian
  
  INTEGER FUNCTION random_binomial(n, p)

! FUNCTION GENERATES A RANDOM BINOMIAL VARIATE USING C.D.Kemp's method.
! This algorithm is suitable when many random variates are required
! with the SAME parameter values for n & p.

!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 <= REAL <= 1)
!    N = NUMBER OF BERNOULLI TRIALS
!           (1 <= INTEGER)

! Reference: Kemp, C.D. (1986). `A modal method for generating binomial
!            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN)    :: p
!     Local variables

INTEGER             :: r0, ru, rd
DOUBLE PRECISION                :: odds_ratio, u, p_r, pd, pu, zero = 0.0, one = 1.0

  r0 = floor((n+1)*p) 
  p_r = bin_prob(n, p, r0)
  odds_ratio = p / (one - p)

u = random_uniform()
u = u - p_r
IF (u < zero) THEN
  random_binomial = r0
  RETURN
END IF

pu = p_r
ru = r0
pd = p_r
rd = r0
DO
  rd = rd - 1
  IF (rd >= 0) THEN
    pd = pd * (rd+1) / (odds_ratio * FLOAT(n-rd))
    u = u - pd
    IF (u < zero) THEN
      random_binomial = rd
      RETURN
    END IF
  END IF

  ru = ru + 1
  IF (ru <= n) THEN
    pu = pu * (n-ru+1) * odds_ratio / FLOAT(ru)
    u = u - pu
    IF (u < zero) THEN
      random_binomial = ru
      RETURN
    END IF
  END IF
END DO

!     This point should not be reached, but just in case:

random_binomial = r0
RETURN

END FUNCTION random_binomial

DOUBLE PRECISION FUNCTION bin_prob(n, p, r)
!     Calculate a binomial probability

IMPLICIT NONE

INTEGER, INTENT(IN) :: n, r
DOUBLE PRECISION, INTENT(IN)    :: p

!     Local variable
REAL                :: one = 1.0

bin_prob = EXP( lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1)) &
                + r*LOG(p) + (n-r)*LOG(one - p) )
RETURN

END FUNCTION bin_prob



DOUBLE PRECISION FUNCTION lngamma(x)

! Logarithm to base e of the gamma function.
!
! Accurate to about 1.e-14.
! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 28 February 1988

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x

!       Local variables

DOUBLE PRECISION :: a1 = -4.166666666554424D-02, a2 = 2.430554511376954D-03,  &
                    a3 = -7.685928044064347D-04, a4 = 5.660478426014386D-04,  &
                    temp, arg, product, lnrt2pi = 9.189385332046727D-1,       &
                    pi = 3.141592653589793D0
LOGICAL          :: reflect

!       lngamma is not defined if x = 0 or a negative integer.

IF (x > 0.d0) GO TO 10
IF (x /= INT(x)) GO TO 10
lngamma = 0.d0
RETURN

!       If x < 0, use the reflection formula:
!               gamma(x) * gamma(1-x) = pi * cosec(pi.x)

10 reflect = (x < 0.d0)
IF (reflect) THEN
  arg = 1.d0 - x
ELSE
  arg = x
END IF

!       Increase the argument, if necessary, to make it > 10.

product = 1.d0
20 IF (arg <= 10.d0) THEN
  product = product * arg
  arg = arg + 1.d0
  GO TO 20
END IF

!  Use a polynomial approximation to Stirling's formula.
!  N.B. The real Stirling's formula is used here, not the simpler, but less
!       accurate formula given by De Moivre in a letter to Stirling, which
!       is the one usually quoted.

arg = arg - 0.5D0
temp = 1.d0/arg**2
lngamma = lnrt2pi + arg * (LOG(arg) - 1.d0 + &
                  (((a4*temp + a3)*temp + a2)*temp + a1)*temp) - LOG(product)
IF (reflect) THEN
  temp = SIN(pi * x)
  lngamma = LOG(pi/temp) - lngamma
END IF
RETURN
END FUNCTION lngamma

end module random
