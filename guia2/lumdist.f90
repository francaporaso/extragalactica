module cosmotab
    integer, parameter :: ucosmo=10, n=6000
    real, dimension(n) :: z, chi
end module cosmotab

program main
    use cosmotab
    
    implicit none
    
    integer :: i
    integer :: usdss = 11, nrows=20000
    real :: luminosity_distance

    character(len=20) :: spec_id
    character(len=20) :: obj_id
    real :: ra, dec, redshift
    real :: pet_u, pet_g, pet_r, pet_i, pet_z
    real :: mod_u, mod_g, mod_r, mod_i, mod_z
    real :: ext_u, ext_g, ext_r, ext_i, ext_z
    real :: petR50_r, petR90_r
    real :: fracDev_r
    real :: velDisp
    
    ! ============ open cosmo.dat table
    
    open(unit=ucosmo, file='../guia1/cosmo.dat', status='old', action='read')
    do i=1, n
        read(ucosmo, *) z(i), chi(i)
    end do
    close(ucosmo)
    
    ! ============ SDSS
    
    open(unit=usdss, file='SDSS_guia2_fmcaporaso.dat', status='old', action='read')    
    read(usdss, '(A)') !skip header
    
    do i=1, nrows
        read(usdss,*)spec_id,obj_id,ra,dec,redshift,pet_u,pet_g,pet_r,pet_i,pet_z,mod_u,mod_g,mod_r,mod_i,mod_z,ext_u,ext_g,ext_&
        &r,ext_i,ext_z,petR50_r,petR90_r,fracDev_r,velDisp
        
        if ((mod_r>=14.5).and.(mod_r<=17.7).and.(petR50_r>1.5)) then
            ! Calculito...          
        end if

    end do
    close(unit=usdss)

end program main
    
function luminosity_distance(z_obj)
    use cosmotab
    implicit none
    real, intent(in) :: z_obj
    integer :: j
    real :: luminosity_distance

    call locate(z, n, z_obj, j)
    luminosity_distance = (1.0+z_obj)*(chi(j) + (z_obj - z(j))*(chi(j+1)-chi(j))/(z(j+1)-z(j)))

end function

function absolute_magnitude(m_ap)
    real, intent(in) :: m_ap

end function


include 'locate.f'

! Mabs = mapa - ext + dm + CorrAB
! dm := modulo de distancia