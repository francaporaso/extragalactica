program main
    implicit none

    type :: t_catalogue
        character(len=20) :: spec_id
        character(len=20) :: obj_id
        real :: ra, dec, redshift
        real :: pet_u, pet_g, pet_r, pet_i, pet_z
        real :: mod_u, mod_g, mod_r, mod_i, mod_z
        real :: ext_u, ext_g, ext_r, ext_i, ext_z
        real :: petR50_r, petR90_r
        real :: fracDev_r
        real :: velDisp
    end type t_catalogue
    
    integer, parameter :: ncols=24, nrows=20000
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: H0=70.0, c0=299792.458
    
    integer :: i, iostat, unit=10
    character(len=250) :: header
    type(t_catalogue), dimension(nrows) :: cat

    real :: mu50, M_r

    real, dimension(nrows) :: Mag_r

    open(unit=unit, file='SDSS_guia2_fmcaporaso.csv', status='old', action='read')    
    ! skip header
    read(unit, '(A)', iostat=iostat) header
    do i=1, nrows
        read(unit,*) cat(i)
    end do
    close(unit=unit)

    mu50 = surface_brightness(cat(1)%pet_r, cat(1)%petR50_r)
    M_r = absolute_magnitude(cat(1)%pet_r, cat(1)%redshift)
    do i=1, nrows
        Mag_r(i) = absolute_magnitude(cat(i)%pet_r, cat(i)%redshift)
    end do
    
    print *, cat(1)%pet_r, cat(1)%petR50_r, mu50, M_r
    ! print *, Mag_r
    
    ! open(unit=unit+1, file='../guia1/cosmo.dat', status='old', action='read')    
    ! do i=1, nrows
    !     read(unit+1,*) z(i), chi(i), dl(i), da(i)
    ! end do
    ! close(unit=unit+1)
    
contains

    real function surface_brightness(mag, rad)
        real, intent(in) :: mag, rad
        surface_brightness = mag + 2.5*log10(pi*rad**2)
    end function surface_brightness

    real function absolute_magnitude(mag, z)
        real, intent(in) :: mag, z
        ! debemos corregir mag aparente por extinción
        ! distancia luminosidad debe ser en pc
        ! M = m - ext + dist_mod

        absolute_magnitude = mag + 5.0 - 5.0*log10(c0*z/H0)   
    end function absolute_magnitude

end program

! ============= Definiciones
! tamaño propio: r = petroR50_r * dA 
! indice de concentracion: C = r90/r50
! distancia propia: d = c*z/H0 (z<<1)
! mag absoluta: M = m + 5 - 5*log10(d)
! AB syst?
! brillo superficial: mu50_r = m + 2.5*log10(pi**2*petroR50_r)
