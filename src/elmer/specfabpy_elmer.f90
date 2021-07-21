function ae2_to_a2(ae2)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: ae2(5)
    real(kind=dp)             :: ae2_to_a2(3,3)
    ae2_to_a2 = ae2_to_a2__sf(ae2)
end

function ae4_to_a4(ae2, ae4)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: ae2(5), ae4(9)
    real(kind=dp)             :: ae4_to_a4(3,3,3,3)
    ae4_to_a4 = ae4_to_a4__sf(ae2, ae4)
end

function a4_IBOF(a2)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: a2(3,3)
    real(kind=dp)             :: a4_IBOF(3,3,3,3)
    
    a4_IBOF = a4_IBOF__sf(a2)
end
