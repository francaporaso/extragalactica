program main
    implicit none
    
    integer, parameter :: ugals = 11, ulf = 12
    integer, parameter :: nrows = 472922
    integer, parameter :: nbins=20
    integer :: i
    integer :: jbin
    real :: M_r, r, z_gal, wi
    real :: dx    
    integer, dimension(nbins) :: bincount = 0
    real, dimension(nbins) :: phi = 0
    
    dx = (-16.0+23.0)/real(nbins)
    
    open(unit=ugals, file='gals_test.dat', status='old', action='read')
    do i=1, nrows
        read(ugals, *) M_r, r, z_gal, wi
        jbin = int((M_r+23.0)/dx) + 1
        
        bincount(jbin) = bincount(jbin) + 1
        phi(jbin) = phi(jbin) + wi
    end do
    close(ugals)
    
    open(unit=ulf, file='lumfunc_n20.dat', status='unknown')
    do i=1,nbins
        write(ulf, *) bincount(i), phi(i), -23.0+(i-1)*dx, -23.0+i*dx
    end do
    close(ulf)

end program main
