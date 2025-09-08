module cosmotab
    implicit none

    integer, parameter :: n=6000
    real, dimension(n) :: z, chi
end module cosmotab

program main
    use cosmotab
    
    implicit none
    
    integer :: i, nskip=0
    integer, parameter :: ucosmo = 10, usdss=11, utable=12
    integer, parameter :: nrows=20000
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: arcsec2rad = pi/(180.0*3600.0)
    real, parameter :: corr_ab_u=-0.036, corr_ab_g=0.012, corr_ab_r=0.010, corr_ab_i=0.028, corr_ab_z=0.040

    character(len=20) :: spec_id
    character(len=20) :: obj_id
    real :: ra, dec, redshift
    real :: pet_u, pet_g, pet_r, pet_i, pet_z
    real :: mod_u, mod_g, mod_r, mod_i, mod_z
    real :: ext_u, ext_g, ext_r, ext_i, ext_z
    real :: petR50_r, petR90_r
    real :: fracDev_r
    real :: velDisp
    
    real :: M_abs_pet_u, M_abs_pet_g, M_abs_pet_r, M_abs_pet_i, M_abs_pet_z
    real :: M_abs_mod_u, M_abs_mod_g, M_abs_mod_r, M_abs_mod_i, M_abs_mod_z
    real :: u_r, g_r
    real :: r50, r90
    real :: conc
    real :: mu50

    real :: dl, dA, dist_mod

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
    write(utable, *) '# ra, dec, z, M_abs_mod_u, M_abs_mod_g, M_abs_mod_r, M_abs_mod_i, M_abs_mod_z, M_abs_&
        &pet_u, M_abs_pet_g, M_abs_pet_r, M_abs_pet_i, M_abs_pet_z, u_r, g_r, r50, r90, conc, fracDev_r, mu50, velDisp'

    do i=1, nrows
        read(usdss,*)spec_id,obj_id,ra,dec,redshift,pet_u,pet_g,pet_r,pet_i,pet_z,mod_u,mod_g,mod_r,mod_i,mod_z,ext_u,ext_g,ext_&
        &r,ext_i,ext_z,petR50_r,petR90_r,fracDev_r,velDisp
        
        if (.not.((mod_r>=14.5).and.(mod_r<=17.7).and.(petR50_r>1.5))) then
            nskip=nskip+1
            cycle
        end if
        
        ! === Magnitudes absolutas
        dl = luminosity_distance(redshift)
        dist_mod = -(5*log10(dl) + 25)

        M_abs_mod_u = mod_u - ext_u + dist_mod + corr_ab_u
        M_abs_mod_g = mod_g - ext_g + dist_mod + corr_ab_g
        M_abs_mod_r = mod_r - ext_r + dist_mod + corr_ab_r
        M_abs_mod_i = mod_i - ext_i + dist_mod + corr_ab_i
        M_abs_mod_z = mod_z - ext_z + dist_mod + corr_ab_z

        M_abs_pet_u = pet_u - ext_u + dist_mod + corr_ab_u
        M_abs_pet_g = pet_g - ext_g + dist_mod + corr_ab_g
        M_abs_pet_r = pet_r - ext_r + dist_mod + corr_ab_r
        M_abs_pet_i = pet_i - ext_i + dist_mod + corr_ab_i
        M_abs_pet_z = pet_z - ext_z + dist_mod + corr_ab_z

        ! === Colores
        u_r = mod_u - mod_r
        g_r = mod_g - mod_r

        ! === Radios petro en kpc
        dA = angular_diam_distance(redshift)*1000.0 ! kpc
        r50 = petR50_r*arcsec2rad*dA
        r90 = petR90_r*arcsec2rad*dA

        ! === Concentración
        conc = petR90_r / petR50_r

        ! === Brillo superficial
        mu50 = mod_r + 2.5*log10(pi*petR50_r**2)

        ! === save to file
        write(utable, '(21(F10.5, 1X))') ra, dec, redshift, M_abs_mod_u, M_abs_mod_g, M_abs_mod_r, M_abs_mod_i, M_abs_mod_z, M_abs_&
            &pet_u, M_abs_pet_g, M_abs_pet_r, M_abs_pet_i, M_abs_pet_z, u_r, g_r, r50, r90, conc, fracDeV_r, mu50, velDisp

    end do

    close(unit=usdss)
    close(unit=utable)

contains

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

end program main

include 'locate.f'
