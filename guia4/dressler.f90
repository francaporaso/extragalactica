program main
    implicit none
    integer, parameter :: utable = 11
    integer, parameter :: nbins = 11
    integer, parameter :: ngals = 5725
    integer :: i, j, jbin
    real :: dx, proj_den, logproj_den
    integer :: ty
    real, parameter :: logden_min = -0.5, logden_max = 2.5
    integer, dimension(nbins) :: ell=0, s0=0, sirr=0, tot=0
    
    open(unit=utable, file='density_data.dat', status='old')
    open(unit=12, file='dressler_data.dat', status='unknown')

    dx = (logden_max - logden_min)/real(nbins)
    
    do i=1, ngals
        read(utable, *) proj_den, ty
        logproj_den = log10(proj_den)

        jbin = int((-logden_min+logproj_den)/dx) + 1
        
        tot(jbin) = tot(jbin) + 1
        if (ty==1) ell(jbin) = ell(jbin) + 1
        if (ty==2) s0(jbin) = s0(jbin) + 1
        if (ty==3) sirr(jbin) = sirr(jbin) + 1
        
    end do
    close(utable)
    
    do j=1, nbins
        write(12, *) logden_min+(j-1)*dx, logden_min+j*dx, logden_min+0.5*(2*j-1)*dx, tot(j), ell(j)/real(tot(j)), &
                     s0(j)/real(tot(j)), sirr(j)/real(tot(j))
    end do
    close(12)

end program main

!test anderson-darlin...