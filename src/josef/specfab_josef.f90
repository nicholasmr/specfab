subroutine mtframe_2d_josef(nlm, m,t, am,at, Emm,Emt, Exx,Exz, fabtype,  Eij_grain,alpha,n_grain)

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(out)   :: m(2),t(2), am,at, Emm,Emt, Exx,Exz, fabtype ! (m,t) = 2D vectors (x,z coords)
    real(kind=dp), intent(in)    :: Eij_grain(2), alpha
    integer, intent(in)          :: n_grain
    
    real(kind=dp), dimension(3)            :: m3,t3, e1,e2,e3, eigvals ! (m3,t3) = 3D (m,t)
    real(kind=dp), dimension(3), parameter :: ex=[1,0,0], ey=[0,1,0], ez=[0,0,1]
    
    call frame(nlm, 'e', e1,e2,e3, eigvals) ! eigen frame (a(e1) >= a(e2) >= a(e3))
    
    m3 = ez
    t3 = ex
    m = [0,1]
    t = [1,0]
    fabtype = 1
    
    ! Single max ?
    if ( (eigvals(1) >= 1./3) .and. (eigvals(2) <= 1./3)) then
        fabtype = 1
        m3 = e1
        m  = e1([1,3])
        am = eigvals(1)
        if (abs(e2(2)) < 1e-20) then ! Vanishing y-comp? Then this eigenvector is in the x--z plane (and is therefore "t") 
            t3 = e2
            t  = e2([1,3])
            at = eigvals(2)
        else
            t3 = e3
            t  = e3([1,3])
            at = eigvals(3)           
        end if
    end if
    
    ! Girdle ?
    if ( (eigvals(1) >= 1./3) .and. (eigvals(2) >= 1./3)) then
        fabtype = 0
        m3 = e3
        m  = e3([1,3])
        am = eigvals(3)
        if (abs(e2(2)) < 1e-20) then ! Vanishing y-comp? Then this eigenvector is in the x--z plane (and is therefore "t") 
            t3 = e2
            t  = e2([1,3])
            at = eigvals(2)
        else
            t3 = e1
            t  = e1([1,3])
            at = eigvals(1)
        end if
    end if
    
    if (m(2)<0) then
        m = -1.0*m
        t = -1.0*t
        m3 = -1.0*m3
        t3 = -1.0*t3
    end if
   
    if ( NORM2(m) .lt. 0.99999999999D0 ) then
        print *,'m=', m
        print *,'norm2(m)=',NORM2(m)
        print *,'e1=', e1
        print *,'e2=', e2
        print *,'e3=', e3
!        stop 'm not normalized!'
    end if
    
    Emm = Evw_tranisotropic(outerprod(m3,m3), tau_vv(m3),     nlm, Eij_grain,alpha,n_grain) ! Longitidinal
    Emt = Evw_tranisotropic(outerprod(m3,t3), tau_vw(m3,t3),  nlm, Eij_grain,alpha,n_grain) ! Shear

    Exx = Evw_tranisotropic(outerprod(ex,ex), tau_vv(ex),     nlm, Eij_grain,alpha,n_grain) ! Longitidinal
    Exz = Evw_tranisotropic(outerprod(ex,ez), tau_vw(ex,ez),  nlm, Eij_grain,alpha,n_grain) ! Shear

end


subroutine mtframe_3d_josef(nlm, m,t, am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype,  Eij_grain,alpha,n_grain)

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(out)   :: m(3),t(3), am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype
    real(kind=dp), intent(in)    :: Eij_grain(2), alpha
    integer, intent(in)          :: n_grain
    
    real(kind=dp), dimension(3)            :: e1,e2,e3, eigvals
    real(kind=dp), dimension(3), parameter :: ex=[1,0,0], ey=[0,1,0], ez=[0,0,1]
    
    call frame(nlm, 'e', e1,e2,e3, eigvals) ! eigen frame (a(e1) >= a(e2) >= a(e3))

    ! if isotropic, this will be the state
    m = ey
    t = ex
    fabtype = 1
    
    ! Single max ?
    if ((eigvals(1) >= 1./3) .and. (eigvals(2) <= 1./3)) then
        fabtype = 1
        m = e1
        t = e2 
        am  = eigvals(1)
        at1 = eigvals(2)
        at2 = eigvals(3)
    end if
    
    ! Girdle ?
    if ((eigvals(1) >= 1./3) .and. (eigvals(2) >= 1./3)) then
        fabtype = 0
        m = e3
        t = e2 
        am  = eigvals(3)
        at1 = eigvals(2)
        at2 = eigvals(1)
    end if
    
    if (m(1)<0.0) then
        m = -1.0*m
        t = -1.0*t
    end if
    
    Emm = Evw_tranisotropic(outerprod(m,m), tau_vv(m),   nlm, Eij_grain,alpha,n_grain) ! Longitidinal
    Emt = Evw_tranisotropic(outerprod(m,t), tau_vw(m,t), nlm, Eij_grain,alpha,n_grain) ! Shear

    Exx = Evw_tranisotropic(outerprod(ex,ex), tau_vv(ex),     nlm, Eij_grain,alpha,n_grain) ! Longitidinal
    Exy = Evw_tranisotropic(outerprod(ex,ey), tau_vw(ex,ey),  nlm, Eij_grain,alpha,n_grain) ! Shear
    Exz = Evw_tranisotropic(outerprod(ex,ez), tau_vw(ex,ez),  nlm, Eij_grain,alpha,n_grain) ! Shear

end

