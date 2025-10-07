program main
    implicit none
    
    integer, parameter :: ugals = 11
    integer, parameter :: nrows = 472922
    integer, parameter :: nbins=50
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
    
    print *, bincount
    print *, phi
    
end program main
