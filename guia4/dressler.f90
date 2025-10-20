program main
    implicit none
    integer, parameter :: utable = 11
    integer, parameter :: nbins = 10
    integer, parameter :: ngals = 5725
    integer :: i, jbin
    real, parameter :: den_min = 0.0, den_max = 2.0
    real :: dx, logden
    integer :: ty
    integer, dimension(nbins) :: ell=0, s0=0, sirr=0

    open(unit=utable, file='distances.dat', status='old')
    
    dx = (den_max - den_min)/real(nbins)

    do i=1, ngals
        read(utable, *) logden, ty
    
        jbin = int((logden-den_max)/dx) + 1
    
        if (ty==1) ell(jbin) = ell(jbin) + 1
        if (ty==2) s0(jbin) = s0(jbin) + 1
        if (ty==3) sirr(jbin) = sirr(jbin) + 1
    
    end do
    
    print*, ell
    print*, s0
    print*, sirr

    close(utable)

end program main