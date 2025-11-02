module galaxias
    implicit none
    integer, parameter :: ngals=5725
    integer, allocatable :: cls_gal(:)
    integer, allocatable :: ig_gal(:) ! galaxy index
    integer, allocatable :: ty_gal(:) ! tipo morf-> 1:E 2:S0 3:S 4:Nosesabe
    real, allocatable :: ra_gal(:), dec_gal(:)
end module galaxias

module clusters
    implicit none
    integer, parameter :: ncluster=56
    integer :: cls_cl
    integer :: nmem
    real :: redshift_cl
    real :: ra_cl, dec_cl
end module clusters

module const
    implicit none
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: d2r = pi/180.0
    real, parameter :: H0 = 70.0 !km/s/Mpc
    real, parameter :: cvel = 299792.458 ! km/s
end module const

program main
    
    use const
    use galaxias
    use clusters
    use numerical

    implicit none
    
    integer, parameter :: ugals=11, ucluster=12, utable=13
    integer :: i, j, k
    integer, allocatable :: indx_gal(:)
    real :: sep, dist_cl, proj_density
    real, allocatable :: dist_gal(:)

    real :: a1, d1, a2, d2

    real, external :: angular_distance
    
    open(unit=ugals, file='galaxias.dat', status='old')
    read(ugals, *) ! read header
    
    open(unit=ucluster, file='centros.dat', status='old')
    read(ucluster, *) ! read header
    
    open(unit=utable, file='density_data.dat', status='unknown')
    
    do i=1, ncluster
        
        read(ucluster, *) cls_cl, nmem, redshift_cl, ra_cl, dec_cl
        dist_cl = redshift_cl*cvel/H0

        allocate(cls_gal(nmem), ig_gal(nmem), ra_gal(nmem), dec_gal(nmem), ty_gal(nmem))
        allocate(dist_gal(nmem), indx_gal(nmem))

        dist_gal = 9999.0
        
        do j=1, nmem
            read(ugals, *) cls_gal(j), ig_gal(j), ra_gal(j), dec_gal(j), ty_gal(j)
        end do
        
        do j=1, nmem
            a1 = ra_gal(j)
            d1 = dec_gal(j)
            do k=1, nmem
                a2 = ra_gal(k)
                d2 = dec_gal(k)
                if (j==k) cycle
                sep = angular_distance(a1, d1, a2, d2)
                dist_gal(k) = tan(sep)*dist_cl
            end do
            call indexx_sp(dist_gal, indx_gal)

            proj_density = 10.0/(pi*dist_gal(indx_gal(10))**2)
            write(utable, *) proj_density, ty_gal(indx_gal(10))
        end do
        
        deallocate(dist_gal, indx_gal)
        deallocate(cls_gal, ig_gal, ra_gal, dec_gal, ty_gal)
        
    end do
    
    close(ucluster)
    close(ugals)
    close(utable)

end program main

real function angular_distance(ra1, dec1, ra2, dec2)
    use const
    implicit none
    real, intent(in) :: ra1, dec1, ra2, dec2
    real :: dra, cos_sep

    ! dra = abs(ra2-ra1)
    dra = d2r*(ra2-ra1)
    cos_sep = sin(d2r*dec1)*sin(d2r*dec2) + cos(d2r*dec1)*cos(d2r*dec2)*cos(dra)
    angular_distance = acos(cos_sep)
    
end function angular_distance