! D. A. Lilien and N. M. Rathmann, 2019-2025
! GOLF from D. A. Lilien, 2025, modified from Gillet-Chaulet and Gagliardini

module specfab_elmer
  
    use rheologies
    use golf 
  
    implicit none 
        
contains
  
  function vectorize9(M) result(vec)
      implicit none
      real(kind=dp) :: M(3,3)
      real(kind=dp) :: vec(9)
      vec = reshape(M, [9])
  end

  function outerprod_to_Mandel(A,B) result (M)

      implicit none
      real(kind=dp), intent(in) :: A(3,3), B(3,3)
      real(kind=dp)             :: M(6,6)
      real(kind=dp), parameter  :: s = sqrt(2.0d0) 

      M(1,:) = [  A(1,1)*B(1,1),   A(1,1)*B(2,2),   A(1,1)*B(3,3), s*A(1,1)*B(2,3), s*A(1,1)*B(1,3), s*A(1,1)*B(1,2)]
      M(2,:) = [  A(2,2)*B(1,1),   A(2,2)*B(2,2),   A(2,2)*B(3,3), s*A(2,2)*B(2,3), s*A(2,2)*B(1,3), s*A(2,2)*B(1,2)]
      M(3,:) = [  A(3,3)*B(1,1),   A(3,3)*B(2,2),   A(3,3)*B(3,3), s*A(3,3)*B(2,3), s*A(3,3)*B(1,3), s*A(3,3)*B(1,2)]
      M(4,:) = [s*A(2,3)*B(1,1), s*A(2,3)*B(2,2), s*A(2,3)*B(3,3), 2*A(2,3)*B(2,3), 2*A(2,3)*B(1,3), 2*A(2,3)*B(1,2)]
      M(5,:) = [s*A(1,3)*B(1,1), s*A(1,3)*B(2,2), s*A(1,3)*B(3,3), 2*A(1,3)*B(2,3), 2*A(1,3)*B(1,3), 2*A(1,3)*B(1,2)]
      M(6,:) = [s*A(1,2)*B(1,1), s*A(1,2)*B(2,2), s*A(1,2)*B(3,3), 2*A(1,2)*B(2,3), 2*A(1,2)*B(1,3), 2*A(1,2)*B(1,2)]
  end

  include "include/IBOF.f90"

  !-------------------
  ! VECTORIZED RHEOLOGIES
  !-------------------
      
  function Cmat_inverse_orthotropic(eps, A, n, m1,m2,m3, Eij) result(C)

      ! Vectorized inverse orthotropic rheology
      ! Returns 9x9 matrix "C" such that vec(tau) = matmul(C, vec(eps)), where vec(tau) and vec(eps) are 9x1 column vectors.

      implicit none
      real(kind=dp), intent(in) :: eps(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
      real(kind=dp)             :: C(9,9), lami(6), ci(6), cvi(6), gam, Mij(6,3,3), Fij(6,3,3), Ii(6), viscosity

      call rheo_params_orthotropic(Eij, n, lami, gam)
      call rheo_structs_orthotropic(eps,m1,m2,m3, 'R', Fij, Ii)
      Mij = Mij_orthotropic(m1,m2,m3)

      cvi(1:3) = 4.0d0/3 * lami(1:3)/gam
      cvi(4:6) =       2 * 1/lami(4:6)

      ci(1:3) = -3/4.0d0 * cvi(1:3)
      ci(4:6) = cvi(4:6)

      viscosity = A**(-1/n) * sum(cvi*Ii**2)**powlawexp_rev(n)
      C = viscosity * ( &
          + ci(1) * outerprod99(vectorize9(identity - 3*Mij(1,:,:)), vectorize9(Mij(1,:,:))) &
          + ci(2) * outerprod99(vectorize9(identity - 3*Mij(2,:,:)), vectorize9(Mij(2,:,:))) &
          + ci(3) * outerprod99(vectorize9(identity - 3*Mij(3,:,:)), vectorize9(Mij(3,:,:))) &
          + ci(4) * outerprod99(vectorize9(symmetricpart(Mij(4,:,:))), vectorize9(Mij(4,:,:))) &
          + ci(5) * outerprod99(vectorize9(symmetricpart(Mij(5,:,:))), vectorize9(Mij(5,:,:))) &
          + ci(6) * outerprod99(vectorize9(symmetricpart(Mij(6,:,:))), vectorize9(Mij(6,:,:))) &
      )
  end

  function Cmandel_inverse_orthotropic(eps, A, n, m1,m2,m3, Eij) result(C)

      ! Vectorized inverse orthotropic rheology
      ! Returns 6x6 Mandel matrix "C" such that mat_to_vec(tau) = matmul(Cmandel, mat_to_vec(eps)), where
      !       mat_to_vec(H) = [H_xx, H_yy, H_zz, sqrt(2)*H_yz, sqrt(2)*H_xz, sqrt(2)*H_xy] 
      ! defines the 6x1 Mandel vector (see mandel.f90).
      ! To get the usual 3x3 form of tau: tau = vec_to_mat( matmul(Cmandel, mat_to_vec(eps)) )

      implicit none
      real(kind=dp), intent(in) :: eps(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
      real(kind=dp)             :: C(6,6), lami(6), ci(6), cvi(6), gam, Mij(6,3,3), Fij(6,3,3), Ii(6), viscosity

      call rheo_params_orthotropic(Eij, n, lami, gam)
      call rheo_structs_orthotropic(eps,m1,m2,m3, 'R', Fij, Ii)
      Mij = Mij_orthotropic(m1,m2,m3)

      cvi(1:3) = 4.0d0/3 * lami(1:3)/gam
      cvi(4:6) =       2 * 1/lami(4:6)

      ci(1:3) = -3/4.0d0 * cvi(1:3)
      ci(4:6) = cvi(4:6)

      viscosity = A**(-1/n) * sum(cvi*Ii**2)**powlawexp_rev(n)

      C = viscosity * ( &
          + ci(1) * outerprod_to_Mandel(identity-3*Mij(1,:,:),Mij(1,:,:)) &
          + ci(2) * outerprod_to_Mandel(identity-3*Mij(2,:,:),Mij(2,:,:)) &
          + ci(3) * outerprod_to_Mandel(identity-3*Mij(3,:,:),Mij(3,:,:)) &
          + ci(4) * outerprod_to_Mandel(symmetricpart(Mij(4,:,:)), symmetricpart(Mij(4,:,:))) &
          + ci(5) * outerprod_to_Mandel(symmetricpart(Mij(5,:,:)), symmetricpart(Mij(5,:,:))) &
          + ci(6) * outerprod_to_Mandel(symmetricpart(Mij(6,:,:)), symmetricpart(Mij(6,:,:))) &
      )
  end

  subroutine Cmat_inverse_orthotropic_dimless(eps, n, m1,m2,m3, Eij, MinInVar, viscosity, C6)

      ! Same as Cmandel_inverse_orthotropic(), but returns the decomposed problem (tensorial and viscosity parts) for Lilien's Elmer/Ice interface.

      implicit none
      real(kind=dp), intent(in)  :: eps(3,3), n, m1(3),m2(3),m3(3), Eij(6), MinInVar
      real(kind=dp)              :: C(6,6), lami(6), ci(6), cvi(6), gam, Mij(6,3,3), Fij(6,3,3), Ii(6)
      real(kind=dp), intent(out) :: viscosity, C6(6,6) 

      call rheo_params_orthotropic(Eij, n, lami, gam)
      call rheo_structs_orthotropic(eps,m1,m2,m3, 'R', Fij, Ii)
      Mij = Mij_orthotropic(m1,m2,m3)

      cvi(1:3) = 4.0d0/3 * lami(1:3)/gam
      cvi(4:6) =       2 * 1/lami(4:6)

      ci(1:3) = -3/4.0d0 * cvi(1:3)
      ci(4:6) = cvi(4:6)

      viscosity = sum(cvi*Ii**2)
      if (viscosity .lt. MinInVar) viscosity = MinInVar
      viscosity = (2.0d0 * viscosity)**powlawexp_rev(n)

      C = + ci(1) * outerprod_to_Mandel(identity-3*Mij(1,:,:),Mij(1,:,:)) &
          + ci(2) * outerprod_to_Mandel(identity-3*Mij(2,:,:),Mij(2,:,:)) &
          + ci(3) * outerprod_to_Mandel(identity-3*Mij(3,:,:),Mij(3,:,:)) &
          + ci(4) * outerprod_to_Mandel(symmetricpart(Mij(4,:,:)), symmetricpart(Mij(4,:,:))) &
          + ci(5) * outerprod_to_Mandel(symmetricpart(Mij(5,:,:)), symmetricpart(Mij(5,:,:))) &
          + ci(6) * outerprod_to_Mandel(symmetricpart(Mij(6,:,:)), symmetricpart(Mij(6,:,:))) 
          
      C6 = C([1, 2, 3, 6, 4, 5], [1, 2, 3, 6, 4, 5])
      C6(1:3, 1:3) = 2.0_dp * C6(1:3, 1:3)
      C6(4:6, 1:3) = C6(4:6, 1:3) / sqrt(2.0_dp)
      C6(1:3, 4:6) = C6(1:3, 4:6) * sqrt(2.0_dp)
  end

  function rheo_rev_orthotropic_dimless(eps, n, m1,m2,m3, Eij) result(tau)

      implicit none
      real(kind=dp), intent(in) :: eps(3,3), n, m1(3),m2(3),m3(3), Eij(6)
      real(kind=dp)             :: tau(3,3), lami(6), gam, Gij(6,3,3), Ji(6)

      call rheo_params_orthotropic(Eij, n, lami, gam)
      call rheo_structs_orthotropic(eps,m1,m2,m3, 'R', Gij, Ji)
      
      lami(1:3) = 4.0d0/3 * lami(1:3)/gam
      lami(4:6) =       2 * 1/lami(4:6)

      tau = singleinner13(lami*Ji, Gij)
  end

  !-------------------
  ! STRUCTURE TENSORS
  !-------------------

  function a4_IBOF(a2) 

      implicit none
      
      real(kind=dp), intent(in) :: a2(3,3)
      real(kind=dp)             :: a4_IBOF(3,3,3,3)
      real(kind=dp)             :: ae2(6), ae4_IBOF(9)
      
      ! IBOF routine assumes dimension ae2(6), but ae2-methods assume ae2(5)
      ae2(1)=a2(1,1)
      ae2(2)=a2(2,2)
      ae2(3)=a2(3,3)
      ae2(4)=a2(1,2)
      ae2(5)=a2(2,3)
      ae2(6)=a2(1,3)
      
      call IBOF(ae2, ae4_IBOF)
      a4_IBOF = ae4_to_a4(a2_to_ae2(a2), ae4_IBOF)
  end

  function a2_to_ae2(a2) result(ae2)

      implicit none
      
      real(kind=dp), intent(in) :: a2(3,3)
      real(kind=dp)             :: ae2(5)
      
      ae2 = 0.0
      
      ae2(1) = a2(1, 1)
      ae2(2) = a2(2, 2)
      ae2(3) = a2(1, 2)
      ae2(4) = a2(2, 3)
      ae2(5) = a2(1, 3)
  end

  function a4_to_ae4(a4) result(ae4)

      !* 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
      
      implicit none
      
      real(kind=dp), intent(in) :: a4(3,3,3,3)
      real(kind=dp)             :: ae4(9)
      
      ae4 = 0.0 ! init
      
      ae4(1) = a4(1,1,1,1)
      ae4(2) = a4(2,2,2,2)
      ae4(3) = a4(1,1,2,2)
      ae4(4) = a4(1,1,2,3)
      ae4(5) = a4(1,2,2,3)
      ae4(6) = a4(1,1,1,3)
      ae4(7) = a4(1,1,1,2)
      ae4(8) = a4(2,2,2,3)
      ae4(9) = a4(1,2,2,2)
  end

  function ae2_to_a2(ae2) result(a2)

      implicit none

      real(kind=dp), intent(in) :: ae2(5)
      real(kind=dp)             :: a2(3,3)
      
      a2 = 0.0
      include "include/ae2_to_a2__body.f90"
  end

  function ae4_to_a4(ae2, ae4) result(a4)

      implicit none

      real(kind=dp), intent(in) :: ae2(5), ae4(9)
      Real(kind=dp)             :: a4(3,3,3,3)
  !    integer :: i,j,k,l

      a4 = 0.0 ! init
      include "include/ae4_to_a4__body.f90"
  end
  

  function Cmat_GOLF_dimless(a2) result(C)
  
      implicit None
      
      real(kind=dp), intent(in)  :: a2(3,3)
      real(kind=dp)              :: C(6,6), ae2(6)
      real(kind=dp)              :: fabricGrid(4878)!, invn, ss
      REAL(KIND=dp)              :: ai(3), Angle(3)
      integer, parameter         :: dim=3

      ae2(1) = a2(1, 1)
      ae2(2) = a2(2, 2)
      ae2(3) = a2(3, 3)
      ae2(4) = a2(1, 2)
      ae2(5) = a2(2, 3)
      ae2(6) = a2(1, 3)

      include "tabulated_vpsc/040010010.f90" ! beta = 0.04, gamma = 1
      CALL R2Ro(ae2,dim,ai,Angle)
      CALL OPILGGE_ai_nl(ai,Angle,FabricGrid,C)
  end

!!!  ORIGINAL VERSION FOR REFERENCE !!!

!  function Cmat_GOLF_dimless(a2, StrainRate, n, beta1, gamma1) result(C)
!  
!      implicit None
!      
!      real(kind=dp), intent(in)  :: a2(3,3), strainrate(3,3), n, beta1, gamma1
!!      real(kind=dp), intent(in), optional :: beta, gammav
!      real(kind=dp)              :: C(6,6), ae2(6), D(6), Stress(3,3)
!      real(kind=dp)              :: fabricGrid(4878), invn, ss
!      REAL(KIND=dp)              :: ai(3), Angle(3)
!      integer                    :: i, j, k, dim
!      character(len=64)          :: viscosityFile
!      LOGICAL                    :: FirstTime=.TRUE.
!       
!      INTEGER, PARAMETER :: INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /) ,&
!           INDj(1:6)=(/ 1, 2, 3, 2, 3, 1 /)

!      SAVE FirstTime, fabricGrid

!      dim=3

!      ! for some reason defined differently here
!      ae2(1) = a2(1, 1)
!      ae2(2) = a2(2, 2)
!      ae2(3) = a2(3, 3)
!      ae2(4) = a2(1, 2)
!      ae2(5) = a2(2, 3)
!      ae2(6) = a2(1, 3)
!!      if (present(beta)) then
!!              beta1 = beta * 1000
!!          else
!!              beta1 = 0.4 * 1000
!!      end if
!!      if (present(gammav)) then
!!              gamma1 = gammav * 100
!!          else
!!              gamma1 = 1.0 * 100
!!      end if

!      IF (FirstTime) THEN
!        !write(viscosityfile,'(a,I3,I3,a)') 'tabulated_vpsc/0', int(beta1*1000), int(gamma1*1000), '10.Va'
!        viscosityFile = 'tabulated_vpsc/040010010.Va'
!        OPEN( 1, File = viscosityFile)
!        DO i=1,813
!            READ( 1, '(6(e14.8))' ) fabricGrid( 6*(i-1)+1:6*(i-1)+6 )
!        END DO
!        CLOSE(1)
!        FirstTime=.False.
!      END IF

!      CALL R2Ro(ae2,dim,ai,Angle)
!      CALL OPILGGE_ai_nl(ai,Angle,FabricGrid,C)

!      ! inverse n
!      invn = (1.0 - n)/(2.0*n)

!      ! format the strain rate for the calculation of stress
!      D(1) = StrainRate(1,1)
!      D(2) = StrainRate(2,2)
!      D(3) = StrainRate(3,3)
!      D(4) = 2. * StrainRate(1,2)
!      D(5) = 2. * StrainRate(2,3)
!      D(6) = 2. * StrainRate(3,1)

!      Stress = 0.
!      DO k = 1, 2*dim
!        DO j = 1, 2*dim
!          Stress( INDi(k),INDj(k) ) = &
!          Stress( INDi(k),INDj(k) ) + C(k,j) * D(j)
!        END DO
!        IF (k > 3)  Stress( INDj(k),INDi(k) ) = Stress( INDi(k),INDj(k) )
!      END DO

!      ss = 0.0_dp
!      DO i = 1, 3
!        DO j = 1, 3
!          ss = ss + Stress(i,j)**2.
!        END DO
!      END DO
!      ss=ss/4.       ! pour avoir le meme resultat si Isotropic
!      ss = (2.*ss)**invn

!      C = C * ss
!  end

end module specfab_elmer
