program main
    implicit none
    
    integer, parameter :: ncols=24, nrows=20000
    integer :: i, iostat, unit=10
    !real :: z(nrows), chi(nrows), dl(nrows), da(nrows)
    character(len=200) :: header, line
    character(len=20) :: spec_id(nrows), photo_id(nrows)
    real :: ra(nrows), dec(nrows), redshift(nrows)
    real :: pet_u(nrows),pet_g(nrows),pet_r(nrows),pet_i(nrows),pet_z(nrows)
    real :: mod_u(nrows),mod_g(nrows),mod_r(nrows),mod_i(nrows),mod_z(nrows)
    real :: ext_u(nrows), ext_g(nrows), ext_r(nrows), ext_i(nrows), ext_z(nrows)
    real :: petroR50_r(nrows), petroR90_r(nrows)
    real :: fracDeV_r(nrows), velDisp(nrows)

    ! tamaño propio: l = theta * dA 
    
    open(unit=unit, file='SDSS_guia2_fmcaporaso.csv', status='old', action='read')    
    ! skip header
    read(unit, '(A)', iostat=iostat) header
    print *, header

    do i=1, nrows
        read(unit, '(A)', iostat=iostat) line
        if (iostat/=0) exit
        read(line, *) spec_id(i), photo_id(i), ra(i), dec(i), redshift(i),& 
        pet_u(i), pet_g(i), pet_r(i), pet_i(i), pet_z(i), mod_u(i), mod_g(i), mod_r(i), mod_i(i), mod_z(i),&
        ext_u(i), ext_g(i), ext_r(i), ext_i(i), ext_z(i), petroR50_r(i), petroR90_r(i), fracDeV_r(i), velDisp(i)
    end do
    close(unit=10)
    
    print *, ra(1), dec(1), redshift(1)
    
        ! open(unit=10, file='../guia1/cosmo.dat', status='old', action='read')    
        ! do i=1, nrows
        !     read(10,*) z(i), chi(i), dl(i), da(i)
        ! end do
        ! close(unit=10)
    
end program

! subroutine routine(arg1,  arg2)
!     type1, intent(in) :: arg1
!     type2, intent(out) ::  arg2

    
! end subroutine routine