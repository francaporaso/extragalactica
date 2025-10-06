program main
    implicit none
    
    integer, parameter :: ugals = 11
    integer, parameter :: nrows = 100
    integer :: i,j
    real, dimension(nrows) :: M_r, wi
    real, dimension(21) :: mbins
    real, dimension(20) :: phi

    open(unit=ugals, file='gals_test.dat', status='old', action='read')

    mbins = (/-23+i*0.35, i=0,21/)
    do i=1, nrows
        read(ugals, *) M_r(i), wi(i)
    end do

    do i=1,nrows
        do j=1,20
            if ((M_r(i)>mbins(j)).and.(M_r(i)<mbins(j+1))) then
                phi(j) = sum(wi)
            end if
        end do
    end do

    
    print *, phi

end program main
