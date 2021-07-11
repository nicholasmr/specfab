subroutine mtframe_josef_2d(nlm, m,t, am,at, Emm,Emt, Exx,Exz, fabtype,  Ecc,Eca,alpha,nprime)
    
    use specfabpy_const
    implicit none
    
    complex(kind=dp), intent(in) :: nlm(:)
    real(kind=dp), intent(out)   :: m(2),t(2), am,at, Emm,Emt, Exx,Exz, fabtype ! (m,t) = 2D vectors (x,z coords)
    real(kind=dp), intent(in)    :: Ecc, Eca, alpha
    integer, intent(in)          :: nprime
    
    call mtframe_2d_josef(nlm, m,t, am,at, Emm,Emt, Exx,Exz, fabtype,  Ecc,Eca,alpha,nprime)
end

subroutine mtframe_josef_3d(nlm, m,t, am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype,  Ecc,Eca,alpha,nprime)

    use specfabpy_const
    implicit none

    complex(kind=dp), intent(in) :: nlm(:)
    real(kind=dp), intent(out)   :: m(3),t(3), am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype 
    real(kind=dp), intent(in)    :: Ecc, Eca, alpha
    integer, intent(in)          :: nprime
    
    call mtframe_3d_josef(nlm, m,t, am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype,  Ecc,Eca,alpha,nprime)
end

