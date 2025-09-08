module lumdist
    ! interpolar z de la tabla para obtener distancia

    include 'locate.f'
    
    implicit none
    
    integer, parameter :: unit=10, n=6000
    integer :: i, j
    real, dimension(n) :: z, chi
    real :: z_obj, chi_itp, dL

contains

    subroutine get_luminosity_distance(z_obj)

        open(unit=unit, file='../guia1/cosmo.dat', status='old', action='read')
        do i=1, n
            read(unit, *) z(i), chi(i)
        end do
        close(unit)
    
        z_obj = 0.54321

        call locate(z, n, z_obj, j)
        chi_itp = chi(j) + (z_obj - z(j))*(chi(j+1)-chi(j))/(z(j+1)-z(j))
    
        dL = (1.0+z_obj)*chi_itp
    
        print *, dL
    end subroutine

end module lumdist

