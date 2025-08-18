program main
    implicit none

    real, parameter :: c = 3.0e8 !m/s
    real, parameter :: H0=70.0, Om0=0.3, Ode0=0.7, Or0=0.0, Ok0=0.0
    real :: z = 1.0, chi
    real, external :: E
    call qromb(E, 0.0, z, chi)
    print *, chi*c/H0
    ! call comoving_distance(z, chi)
    ! print *, chi
    ! call luminous_distance
    ! call angular_diameter_distance

end program main

real function E(z)
    real, intent(in) :: z
    E = sqrt(Or0*(1+z)**4 + Om0*(1+z)**3 + Ok0*(1+z)**2 + Ode0)    
end function E

! subroutine comoving_distance(z, chi)
!     real, intent(in) :: z
!     real, intent(out) :: chi
!     real, external :: E, qromb
    
!     chi = c/H0 * qromb(E, 0.0, z, 1.0)
! end subroutine comoving_distance

! subroutine luminous_distance(arg1,  arg2)
!     type1, intent(in) :: arg1
!     type2, intent(out) ::  arg2

    
! end subroutine luminous_distance

! subroutine angular_diameter_distance(arg1,  arg2)
!     type1, intent(in) :: arg1
!     type2, intent(out) ::  arg2
! end subroutine angular_diameter_distance