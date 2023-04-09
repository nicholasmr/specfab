subroutine mtframe_josef_2d(nlm, m,t, am,at, Emm,Emt, Exx,Exz, fabtype,  Eij_grain,alpha,n_grain)
    
    use specfabpy_const
    implicit none
    
    complex(kind=dp), intent(in) :: nlm(:)
    real(kind=dp), intent(out)   :: m(2),t(2), am,at, Emm,Emt, Exx,Exz, fabtype ! (m,t) = 2D vectors (x,z coords)
    real(kind=dp), intent(in)    :: Eij_grain(2), alpha
    integer, intent(in)          :: n_grain
    
    call mtframe_2d_josef(nlm, m,t, am,at, Emm,Emt, Exx,Exz, fabtype,  Eij_grain,alpha,n_grain)
end

subroutine mtframe_josef_3d(nlm, m,t, am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype,  Eij_grain,alpha,n_grain)

    use specfabpy_const
    implicit none

    complex(kind=dp), intent(in) :: nlm(:)
    real(kind=dp), intent(out)   :: m(3),t(3), am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype 
    real(kind=dp), intent(in)    :: Eij_grain(2), alpha
    integer, intent(in)          :: n_grain
    
    call mtframe_3d_josef(nlm, m,t, am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype,  Eij_grain,alpha,n_grain)
end

