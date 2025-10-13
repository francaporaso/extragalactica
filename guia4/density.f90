module galaxias
    implicit none
    integer, parameter :: ngals=5725
    integer, dimension(ngals) :: cls_gal
    integer, dimension(ngals) :: ig_gal ! galaxy index
    integer, dimension(ngals) :: ty_gal ! tipo morf-> 1:E 2:S 3:Irr 4:Nosesabe
    real, dimension(ngals) :: ra_gal, dec_gal
end module galaxias

module clusters
    implicit none
    integer, parameter :: ncluster=56
    integer :: cls_cl
    integer :: nmember
    real :: redshift_cl
    real :: ra_cl, dec_cl
end module clusters

module const
    implicit none
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: d2r = pi/180.0
end module const

program main
    use const
    use galaxias
    use clusters
    implicit none
    integer, parameter :: ugals=11, ucluster=12
    integer :: i, j
    real, dimension(79) :: ra, dec

    open(unit=ugals, file='galaxias.dat', status='old')
    read(ugals, *) ! read header
    do i=1, ngals
        read(ugals, *) cls_gal(i), ig_gal(i), ra_gal(i), dec_gal(i), ty_gal(i)
    end do
    close(ugals)
   
    open(unit=ucluster, file='centros.dat', status='old')
    read(ucluster, *) ! read header
    read(ucluster, *) cls_cl, nmember, redshift_cl, ra_cl, dec_cl
    
    !probar para 1 cluster primero
    ! do i=1,79
    !     ra(i), dec(i)
    ! end do

    ! do i=1, 79
    !     a1 = ra(i)
    !     d1 = dec(i)
    !     do j=1,79
    !         if (j==i) cycle
    !         sep = angular_distance(a1, d1, ra(j), dec(j))
    !     end do
    ! end do
    
    close(ucluster)


end program main

subroutine angular_distance(ra1, dec1, ra2, dec2, sep)
    use const
    implicit none
    real, intent(in) :: ra1, dec1, ra2, dec2
    real, intent(out) :: sep
    real :: dra, cos_sep

    dra = abs(ra2-ra1)
    cos_sep = sin(d2r*dec1)*sin(d2r*dec2) + cos(d2r*dec1)*cos(d2r*dec2)*cos(d2r*dra)
    sep = acos(cos_sep)/d2r
    
end subroutine angular_distance