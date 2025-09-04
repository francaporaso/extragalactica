module interpolate
    implicit none
    ! interpolar z de la tabla para obtener distancia
    integer, parameter :: unit=10, n=6000
    integer :: i, j
    real, dimension(n) :: z, chi
    
    open(unit=unit, file='../guia1/cosmo.dat', status='old', action='read')
    do i=1, n
        read(unit, *) z(i), chi(i)
    end do
    close(unit)
    
contains
    
    real function get_dL(z)
        real, intent(in) :: z
        integer :: j
        real :: z_obj, chi_itp
    
        call locate(z, n, z_obj, j)
        chi_itp = chi(j) + (z_obj - z(j))*(chi(j+1)-chi(j))/(z(j+1)-z(j))

        get_dL = (1.0+z_obj)*chi_itp
        
    end function get_chi

end module interpolate

program main
    use interpolate
    implicit none
    integer :: i
    real, dimension(10) :: dL

    z_obj = (/(i*0.1, i=1,10)/)
    do i = 1, 10
        dL(i) = get_dL(z_obj(i))
    end do

    print *, dL
end program main

include 'locate.f'
