module distancetab
    implicit none
    integer, parameter :: n=6000
    real, dimension(n) :: z, chi
end module distancetab

module cosmology
    implicit none
    real :: Om0 = 0.3
    real :: Ode0 = 0.7
    real :: H0 = 70.0 !km/s/Mpc
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

    integer :: i,k
    real :: r
    real :: dL, dist_mod
    real :: M_abs_r

    real :: luminosity_distance, angular_diam_distance

    ! ============ open cosmo.dat table    
    open(unit=ucosmo, file='../guia1/cosmo.dat', status='old', action='read')
    do i=1, n
        read(ucosmo, *) z(i), chi(i)
    end do
    close(ucosmo)
    ! ============

    ! ============ get galaxies
    open(unit=utable, file='gals_test.dat', status='unknown')
    write(utable, *) '# z_gal, M_abs_r, r, K00_{ugriz}, K01_{ugriz}'
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

        write(utable, *) z_gal, M_abs_r, r, (rk_p(k), k=1,5), (rks_p(k), k=1,5)
    end do
    close(udata)
    close(utable)
    ! ============

end program main

real function luminosity_distance(z_obj)
    use distancetab
    implicit none
    real, intent(in) :: z_obj
    integer :: j

    call locate(z, n, z_obj, j)
    luminosity_distance = (1.0+z_obj)*(chi(j) + (z_obj - z(j))*(chi(j+1)-chi(j))/(z(j+1)-z(j)))
end function

real function angular_diam_distance(z_obj)
    use distancetab
    implicit none
    real, intent(in) :: z_obj
    integer :: j

    call locate(z, n, z_obj, j)
    angular_diam_distance = (chi(j) + (z_obj - z(j))*(chi(j+1)-chi(j))/(z(j+1)-z(j)))/(1.0+z_obj)
end function

include '../locate.f'
