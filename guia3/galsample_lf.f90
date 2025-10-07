module distancetab
    implicit none
    integer, parameter :: n=6000
    real, dimension(n) :: z, chi, lum_dist
end module distancetab

module cosmology
    implicit none
    real, parameter :: Om0 = 0.3
    real, parameter :: Ode0 = 0.7
    real, parameter :: H0 = 70.0 !km/s/Mpc
    real, parameter :: c = 299792.458 !km/s
end module cosmology

module catalogue
    implicit none
    real :: z_gal
    real :: pet_r, ext_r
    real :: r50
    real, dimension(5) :: rk_p
    real, dimension(5) :: rks_p
end module catalogue

program main
    use distancetab
    use cosmology
    use catalogue

    implicit none
    
    integer, parameter :: udata = 11, ucosmo = 12, utable = 13
    integer, parameter :: nrows = 685034
    real, parameter :: corr_ab_r = 0.01
    real, parameter :: rmax = 17.77, rmin = 14.5

    integer :: i,k
    real :: r
    real :: dL, dist_mod
    real :: M_abs_r
    real :: dL_max, zmax
    real :: Vmax, wi
    
    real, external :: luminosity_distance, integrand_vmax, get_z
    
    ! ============ open cosmo.dat table    
    open(unit=ucosmo, file='../guia1/cosmo.dat', status='old', action='read')
    do i=1, n
        read(ucosmo, *) z(i), chi(i), lum_dist(i)
    end do
    close(ucosmo)
    ! ============

    ! ============ get galaxies
    open(unit=utable, file='gals_test.dat', status='unknown')
    open(unit=udata, file='mgs.dat', status='old', action='read')
    do i=1,nrows
        read(udata, *) z_gal, pet_r, ext_r, r50, (rk_p(k), k=1,5), (rks_p(k), k=1,5)
        r = pet_r - ext_r
        if ((r<14.5).or.(r>17.77)) cycle
        if (r50<1.5) cycle
        if ((z_gal>0.15).or.(z_gal<=0.0)) cycle

        dL = luminosity_distance(z_gal) ! Mpc
        dist_mod = -(5.0*log10(dL)+25.0)

        M_abs_r = r + dist_mod + corr_ab_r + rks_p(3)
        if ((M_abs_r<-23).or.(M_abs_r>-16)) cycle

        ! === Vmax weights
        dL_max = 10.0**(-0.2*(M_abs_r - rmax + 25.0))
        zmax = get_z(dL_max)
        call qromb(integrand_vmax, 0.0, zmax, Vmax)
        wi = H0/(Vmax*c)

        write(utable,*) M_abs_r, r+corr_ab_r, z_gal, wi
    
    end do
    close(udata)
    close(utable)
    ! ============

end program main

real function get_z(dL_obj)
    use distancetab
    implicit none
    real, intent(in) :: dL_obj
    integer :: j

    call locate(lum_dist, n, dL_obj, j)
    get_z = z(j) + (dL_obj - lum_dist(j))*(z(j+1)-z(j))/(lum_dist(j+1)-lum_dist(j))
end function

! martinez, zandivarez 2002 (ec 10)
real function integrand_vmax(z_obj)
    use cosmology
    implicit none
    real, intent(in) :: z_obj
    real :: luminosity_distance

    integrand_vmax = (luminosity_distance(z_obj)/(1+z_obj))**2/sqrt(Om0*(1+z_obj)**3 + (1-Om0-Ode0)*(1+z_obj)**2 + Ode0)
end function

real function luminosity_distance(z_obj)
    use distancetab
    implicit none
    real, intent(in) :: z_obj
    integer :: j

    call locate(z, n, z_obj, j)
    luminosity_distance = lum_dist(j) + (z_obj - z(j))*(lum_dist(j+1)-lum_dist(j))/(z(j+1)-z(j))
end function

include '../locate.f'
include '../qromb.f'
include '../trapzd.f'
include '../polint.f'