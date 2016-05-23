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

end module random
