module cosmotab
    implicit none
    integer, parameter :: n=6000
    real, dimension(n) :: z, chi
end module cosmotab

module galsample
    implicit none
    real, dimension(5) :: M_abs_pet
    real, dimension(5) :: M_abs_mod
    real :: u_r, g_r
    real :: r50, r90
    real :: c9050 ! concentration index
    real :: mu50 ! surface brigtness
end module galsample

module catalogue
    implicit none
    character(len=20) :: spec_id
    character(len=20) :: obj_id
    real :: ra, dec, redshift
    real, dimension(5) :: petromag
    real, dimension(5) :: modelmag
    real, dimension(5) :: extintion
    real :: petR50_r, petR90_r
    real :: fracDev_r
    real :: velDisp 
end module catalogue

program main
    use cosmotab
    use catalogue
    use galsample

    implicit none
    
    integer :: i,k
    integer, parameter :: ucosmo = 10, usdss=11, utable=12
    integer, parameter :: nrows=20000
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: arcsec2rad = pi/(180.0*3600.0)
    real, dimension(5), parameter :: corr_ab = (/-0.036,0.012,0.010,0.028,0.040/)

    real :: luminosity_distance, angular_diam_distance
    real :: dL, dA, dist_mod, r

    ! ============ open cosmo.dat table    
    open(unit=ucosmo, file='../guia1/cosmo.dat', status='old', action='read')
    do i=1, n
        read(ucosmo, *) z(i), chi(i)
    end do
    close(ucosmo)
    
    ! ============ SDSS
    open(unit=usdss, file='SDSS_guia2_fmcaporaso.dat', status='old', action='read')    
    read(usdss, '(A)') !skip header
    
    ! ============ Nueva tabla
    open(unit=utable, file='gals.dat', status='unknown')
    write(utable, '(21(A12, 1X))') 'ra','dec','z','M_mod_u','M_mod_g','M_mod_r','M_mod_i','M_mod_z','M_&
        &pet_u','M_pet_g','M_pet_r','M_pet_i','M_pet_z','u_r','g_r','r50','r90','c9050','fracDeV_r','mu50','velDisp'

    do i=1, nrows
        read(usdss,*) spec_id,obj_id,ra,dec,redshift,(petromag(k),k=1,5),(model&
        &mag(k),k=1,5),(extintion(k),k=1,5),petR50_r,petR90_r,fracDev_r,velDisp

        r = petromag(3) - extintion(3)

        ! === gold sample
        if ((r<14.5).or.(r>17.77).or.(petR50_r<1.5)) cycle
        if (modelmag(1) < -100) cycle ! hay 1 gx con mod_u = -9999, excluido
 
        dL = luminosity_distance(redshift) ! Mpc
        dA = angular_diam_distance(redshift)*1000.0 ! kpc
        dist_mod = -(5*log10(dl) + 25)
        
        do k=1, 5
            petromag(k) = petromag(k)-extintion(k)
            modelmag(k) = modelmag(k)-extintion(k)
        end do

        do k=1,5
            M_abs_pet(k) = petromag(k)+dist_mod+corr_ab(k)
            M_abs_mod(k) = modelmag(k)+dist_mod+corr_ab(k)
        end do
        
        ! === Colores
        u_r = modelmag(1) - modelmag(3)
        g_r = modelmag(2) - modelmag(3)

        ! === Radios petro en kpc
        r50 = petR50_r*arcsec2rad*dA
        r90 = petR90_r*arcsec2rad*dA

        ! === Concentración
        c9050 = petR90_r / petR50_r

        ! === Brillo superficial
        mu50 = petromag(3) + 2.5*log10(2*pi*petR50_r**2) !

        ! === save to file
        write(utable, '(21(F12.7, 1X))') ra, dec, redshift, (M_abs_mod(k),k=1,5), (M_abs_&
        &pet(k),k=1,5), u_r, g_r, r50, r90, c9050, fracDeV_r, mu50, velDisp

    end do

    close(unit=usdss)
    close(unit=utable)

end program main

real function luminosity_distance(z_obj)
    use cosmotab
    implicit none
    real, intent(in) :: z_obj
    integer :: j 

    call locate(z, n, z_obj, j)
    luminosity_distance = (1.0+z_obj)*(chi(j) + (z_obj - z(j))*(chi(j+1)-chi(j))/(z(j+1)-z(j)))
end function

real function angular_diam_distance(z_obj)
    use cosmotab
    implicit none
    real, intent(in) :: z_obj
    integer :: j

    call locate(z, n, z_obj, j)
    angular_diam_distance = (chi(j) + (z_obj - z(j))*(chi(j+1)-chi(j))/(z(j+1)-z(j)))/(1.0+z_obj)
end function

include 'locate.f'
