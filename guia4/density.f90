module galaxias
    implicit none
    integer, parameter :: ngals=5725
    integer, allocatable :: cls_gal(:)
    integer, allocatable :: ig_gal(:) ! galaxy index
    integer, allocatable :: ty_gal(:) ! tipo morf-> 1:E 2:S 3:Irr 4:Nosesabe
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
    real :: sep, dist_cl
    real, allocatable :: dist_gal(:)

    real, external :: angular_distance
    
    open(unit=ugals, file='galaxias.dat', status='old')
    read(ugals, *) ! read header
    
    open(unit=ucluster, file='centros.dat', status='old')
    read(ucluster, *) ! read header
    
    !open(unit=utable, file='distances.dat', status='unknown')

    do i=1, ncluster
        
        read(ucluster, *) cls_cl, nmem, redshift_cl, ra_cl, dec_cl
        dist_cl = redshift_cl*cvel/H0

        allocate(cls_gal(nmem), ig_gal(nmem), ra_gal(i), dec_gal(nmem), ty_gal(nmem))
        allocate(dist_gal(nmem), indx_gal(nmem))
                
        do j=1, nmem        
            read(ugals, *) cls_gal(j), ig_gal(j), ra_gal(i), dec_gal(j), ty_gal(j)
        end do

        do j=1, nmem
            do k=1, nmem
                if (k==j) cycle
                sep = angular_distance(ra_gal(j), dec_gal(j), ra_gal(k), dec_gal(k))
                dist_gal(k) = tan(sep)*dist_cl
            end do

            call indexx_sp(dist_gal, indx_gal)
            print *, dist_gal(indx_gal(10))
        
        end do
        
        deallocate(cls_gal, ig_gal, ra_gal, dec_gal, ty_gal, dist_gal, indx_gal)
    
    end do
    
    close(ucluster)
    close(ugals)
    !close(utable)

end program main

real function angular_distance(ra1, dec1, ra2, dec2)
    use const
    implicit none
    real, intent(in) :: ra1, dec1, ra2, dec2
    real :: dra, cos_sep

    dra = abs(ra2-ra1)
    cos_sep = sin(d2r*dec1)*sin(d2r*dec2) + cos(d2r*dec1)*cos(d2r*dec2)*cos(d2r*dra)
    angular_distance = acos(cos_sep)
    
end function angular_distance