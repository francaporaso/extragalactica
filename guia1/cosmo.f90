program main
    implicit none
    real, parameter :: c = 299792.458 !km/s
    real, parameter :: H0 = 70.0 !km/s/Mpc
    real, parameter :: Om0=0.3, Ode0=0.7, Or0=0.0
    real, parameter :: Ok0 = 0.0
    real :: zi, Dz, integral, chi
    real :: dl, da
    integer :: i, n, u

    Dz = 0.001
    n = 6000
    u = 10

    open(unit=u, file='cosmo.dat', status='unknown')
    !write(u, '(A10, 1X, A10, 1X, A10, 1X, A10)')'z', 'chi', 'd_L', 'd_A'
    
    do i=1, n
        zi = 0.0 + Dz*(i-1)
        call qromb(E, 0.0, zi, integral)
        
        chi = (c/H0)*integral
        dl = (1.0+zi)*chi
        da = chi/(1.0+zi)
        
        write(u, '(F10.3, A, F10.3, A, F10.3, A, F10.3)') zi, ',', chi, ',', dl, ',', da   
    end do
    
    close(u)

contains

    real function E(z)
        real, intent(in) :: z
        E = 1.0/sqrt(Or0*(1.0+z)**4 + Om0*(1.0+z)**3 + Ok0*(1.0+z)**2 + Ode0)
    end function E
    
end program main

include '../qromb.f'
include '../trapzd.f'
include '../polint.f'

