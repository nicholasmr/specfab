!/*****************************************************************************/
! *
! *  Modified from Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini, Ga¨el Durand, David Lilien
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> tools for utilizing the GOLF anisotropic law
      !
      ! Definition of the interpolation grid 
      !
  
      Module GOLF
        use header

!     INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)  ! If not using
                                                         ! with Elmer   
      
      Real, Parameter :: kmin=0.002_dp ! valeur de ki mimum
      Integer, Parameter :: Ndiv=30    ! Ndiv+2 Number of points along ik1
      Integer, Parameter :: Ntot=813   ! Total number of points
      Integer, Parameter :: NetaI=4878 ! 6*4884 length of EtaI
      Integer, Parameter, Dimension(32) :: Nk2 = & 
                (/ -1,  46,  93, 139, 183, 226, 267, 307, 345, 382,&
                 417, 451, 483, 514, 543, 571, 597, 622, 645, 667,& 
                 687, 706, 723, 739, 753, 766, 777, 787, 795, 802,&
                 807, 811/) 

      contains
      

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!   Expression of the viscosity matrix in the reference frame    !!!!!!
!!!!!! 							              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  from  the 6 eta6 viscoisties of the matrice law  
!
!  S = eta6(k) tr(M_kd)M_k^D + eta6(k+3)(M_k d + d M_k)^D 
!     
!   get the Voigt eta36  expressed in the general reference frame
!      
!     Si = Aij dj 
!     Aij = eta(i,j) i=1,6; j=1,6 non-symetric matrix 
!      where S=(S11,S22,S33,S12,S23,S31)      
!      and   d=(d11,d22,d33,2 d12,2 d23,2 d31)      
!      Voigt notation as in ELMER
!

       Subroutine ViscGene(eta6,Angle,eta36)
       Implicit None
       Real(kind=dp), Intent(in), Dimension(3) :: Angle    ! Euler Angles        
       Real(kind=dp), Intent(in), Dimension(6) :: Eta6     ! 6 Viscosities of the
                                                  ! matrice law

       Real(kind=dp), Intent(out), Dimension(6,6) :: eta36  ! Viscosity matrix in 
                                                   ! the reference frame
        
       Real(kind=dp) :: p,t,o,st,ct
       Real(kind=dp), Dimension(3,3) :: Q
       Real(kind=dp) :: dQk 
       Real(kind=dp), Dimension(6), Parameter :: coef=(/1.,1.,1.,.0,.0,.0/)
       Integer, Dimension(6), Parameter :: ik=(/1,2,3,1,2,3/)
       Integer, Dimension(6), Parameter :: jk=(/1,2,3,2,3,1/)
       Integer :: k,m,n
       Integer :: i,j 
       

!  Angle = Phi, Theta, Omega
       
       p=Angle(1)
       t=Angle(2)
       o=Angle(3)

! terms of the rotation matrix from RG to RO

        ct = cos(t)
        Q(3,3) = ct
        Q(1,3) = sin(o)
        Q(2,3) = cos(o)
        Q(3,1) = sin(p)
        Q(3,2) = cos(p)

        Q(1,1) = Q(3,2)*Q(2,3) - Q(3,1)*Q(1,3)*ct 
        Q(1,2) = Q(3,1)*Q(2,3) + Q(3,2)*Q(1,3)*ct 
        Q(2,1) = - Q(3,2)*Q(1,3) - Q(3,1)*Q(2,3)*ct 
        Q(2,2) = - Q(3,1)*Q(1,3) + Q(3,2)*Q(2,3)*ct 

        st = sin(t)
        Q(1,3) = Q(1,3)*st
        Q(2,3) = Q(2,3)*st
        Q(3,1) = Q(3,1)*st
        Q(3,2) = -Q(3,2)*st

! 36  terms of the Voigt matrix in the reference frame 
        Do m=1,6
          Do n=1,6
            eta36(m,n) = 0.
          End Do
        End Do
        Do k=1,3
          dQk=Q(k,1)*Q(k,1)+Q(k,2)*Q(k,2)+Q(k,3)*Q(k,3)
          Do m=1,6
            Do n=1,6
              eta36(m,n)=eta36(m,n)+eta6(k)*Q(k,ik(n))*Q(k,jk(n))*&
               (Q(k,ik(m))*Q(k,jk(m))-1./3.*coef(m)*dQk)    
              eta36(m,n) = eta36(m,n) - 2./3.*eta6(k+3)*coef(m)* &
               Q(k,ik(n))*Q(k,jk(n))
            End Do
          End Do

          eta36(1,1)=eta36(1,1)+eta6(k+3)*Q(k,1)*Q(k,1)*2.
          eta36(1,4)=eta36(1,4)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(1,6)=eta36(1,6)+eta6(k+3)*Q(k,1)*Q(k,3)

          eta36(2,2)=eta36(2,2)+eta6(k+3)*Q(k,2)*Q(k,2)*2.
          eta36(2,4)=eta36(2,4)+eta6(k+3)*Q(k,2)*Q(k,1)
          eta36(2,5)=eta36(2,5)+eta6(k+3)*Q(k,2)*Q(k,3)
          
          eta36(3,3)=eta36(3,3)+eta6(k+3)*Q(k,3)*Q(k,3)*2.
          eta36(3,5)=eta36(3,5)+eta6(k+3)*Q(k,3)*Q(k,2)
          eta36(3,6)=eta36(3,6)+eta6(k+3)*Q(k,3)*Q(k,1)

          eta36(4,1)=eta36(4,1)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(4,2)=eta36(4,2)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(4,4)=eta36(4,4)+eta6(k+3)*(dQk-Q(k,3)*Q(k,3))*0.5
          eta36(4,5)=eta36(4,5)+eta6(k+3)*Q(k,1)*Q(k,3)*0.5
          eta36(4,6)=eta36(4,6)+eta6(k+3)*Q(k,2)*Q(k,3)*0.5

          eta36(5,2)=eta36(5,2)+eta6(k+3)*Q(k,2)*Q(k,3)
          eta36(5,3)=eta36(5,3)+eta6(k+3)*Q(k,2)*Q(k,3)
          eta36(5,4)=eta36(5,4)+eta6(k+3)*Q(k,1)*Q(k,3)*0.5
          eta36(5,5)=eta36(5,5)+eta6(k+3)*(dQk-Q(k,1)*Q(k,1))*0.5
          eta36(5,6)=eta36(5,6)+eta6(k+3)*Q(k,1)*Q(k,2)*0.5

          eta36(6,1)=eta36(6,1)+eta6(k+3)*Q(k,1)*Q(k,3)
          eta36(6,3)=eta36(6,3)+eta6(k+3)*Q(k,1)*Q(k,3)
          eta36(6,4)=eta36(6,4)+eta6(k+3)*Q(k,2)*Q(k,3)*0.5
          eta36(6,5)=eta36(6,5)+eta6(k+3)*Q(k,1)*Q(k,2)*0.5
          eta36(6,6)=eta36(6,6)+eta6(k+3)*(dQk-Q(k,2)*Q(k,2))*0.5
        End Do
         
         End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   sens = 1 : tri des ki pour avoir k1 < k2 < k3
!   sens =-1 : tri des visc apres calcul avec k1 < k2 < k3 
!
!
       Subroutine triki(ki0,ki,visc,ordre,sens)
       implicit none
       
       Real(kind=dp), Dimension(3) ::  ki0,ki
       real(kind=dp), Dimension(6) :: visc,b
       real(kind=dp) :: a
       Integer :: sens,i,j
       Integer, Dimension(3) :: ordre 
       
       

!
!    Passage pour trier les ki
!
       If (sens.EQ.1) Then
         Do i=1,3
           ki(i)=ki0(i)
           ordre(i)=i
         End Do
         Do j=2,3
          a=ki(j)
          Do i=j-1,1,-1
          If (ki(i).LE.a) Goto 20
          ki(i+1)=ki(i)
          ordre(i+1)=ordre(i)
          End Do
  20      Continue
         ki(i+1)=a
         ordre(i+1)=j
         End Do
!
!   Passage pour remettre les viscosite dans le bon ordre
!
       ElseIf (sens.EQ.-1) Then

         Do i=1,6
         b(i)=visc(i)
         End Do
         Do i=1,3
          visc(ordre(i))=b(i)
          visc(ordre(i)+3)=b(i+3)
        End Do

       Else
         Write(*,*)'triki.f : sens <> 1 ou -1' 
       Stop
       End If
       End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Interpolation quadratique d'une quantite Q definie en x1,x2,x3
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      Function InterP(t,x,Q)
!
      Implicit None
      Real(kind=dp), Dimension(3) :: x,Q
      Real(kind=dp) :: t,InterP,d12,d23
      Real(kind=dp) :: Ip

!
      d12=x(2)-x(1)
      d23=x(3)-x(2)
      Ip=Q(1)*(x(2)-t)*(x(3)-t)/((d12+d23)*d12)
      Ip=Ip+Q(2)*(t-x(1))*(x(3)-t)/(d12*d23)
      Ip=Ip+Q(3)*(t-x(1))*(t-x(2))/((d12+d23)*d23)

      InterP = Ip
      Return
      End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Interpolation quadratique d'une quantite Q definie en 9 points 
!
!         y3 ---     Q7 ---- Q8 ------ Q9
!                    |        |         | 
!                    |        |         |
!         y2 ---     Q4 ----- Q5 ----- Q6
!                    |        |         | 
!                    |        |         |
!         y1 ---     Q1 ----- Q2 ----- Q3
!
!                    |        |         |
!                   x1        x2       x3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      Function InterQ9(x,y,xi,yi,Q)
!
      Implicit None
      Real(KIND=dp), Dimension(3) ::  xi,yi,a
      Real(kind=dp), Dimension(9) ::  Q
      Real(kind=dp) ::  InterQ9
      Real(kind=dp) ::  Ip,x,y
      Integer  i

!
         a(1)=InterP(x,xi,Q(1))
         a(2)=InterP(x,xi,Q(4))
         a(3)=InterP(x,xi,Q(7))
         Ip=InterP(y,yi,a)

         InterQ9 = Ip

      Return
      End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!       Orthotropic Polycrystalline Ice Law from LGGE            !!!!!!
!!!!!! 							              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   give the Voigt viscosity matrix eta36 expressed 
!                                    in the general reference frame
!     Si = Aij dj 
!     Aij = eta36(i,j) i=1,6; j=1,6 non-symetric matrix 
!      where S=(S11,S22,S33,S12,S23,S31)      
!      and   d=(d11,d22,d33,2d12,2d23,2d31)      
!         as in Elmer
!
       
       Subroutine OPILGGE_ai_nl(ai,Angle,etaI,eta36)
       
       Implicit None
       Real(kind=dp), Intent(in), Dimension(3) :: ai       ! Texture parameters 
       Real(kind=dp), Dimension(3) :: ki       ! Texture parameters 
       Real(kind=dp), Intent(in), Dimension(3) :: Angle    ! Euler Angles               
       Real(kind=dp), Intent(in), Dimension(:) :: etaI ! Grid Viscosities 
       Real(kind=dp), Intent(out), Dimension(6,6) :: eta36
       Real(kind=dp), Dimension(6) :: eta6
       Real(kind=dp) :: aplusa
       Integer :: i
        
            Do i=1,3
              ki(i)=ai(i)
            End do

            Do i=1,3
             ki(i)=min(Max(ki(i),0._dp),1._dp)
            End do
            aplusa=ki(1)+ki(2)+ki(3)
            If (aplusa.GT.1._dp) then
                    do i=1,3
                     ki(i)=ki(i)/aplusa
                    end do
            endif
       
       ! Value of the 6 relative viscosities in the orthotropic frame

       Call ViscMat_ai(ki,eta6,etaI)

       ! Viscosities in the reference frame

       Call ViscGene(eta6,Angle,eta36)

       End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc									cccc
!ccccc       subroutine de calcul des viscosites Analytiques            cccc
!ccccc       Les viscosites ont ete calculees sur une grille            cccc
!ccccc       avec le sous-prog makeEtaI.f                               cccc
!ccccc         Elles passent dans etaI(6*814)=EtaI(4884)                cccc
!ccccc 									cccc      
!ccccc      !!! En entree a1,a2,a3 les valeurs propres du tenseur       cccc
!ccccc                                          d'orientation           cccc
!ccccc 									cccc      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       
       Subroutine ViscMat_ai(ai0,eta6,etaI)
       
       Implicit None
       Real(kind=dp), Intent(in), Dimension(3) :: ai0
       Real(kind=dp), Intent(out), Dimension(6) :: eta6
       Real(kind=dp), Intent(in), Dimension(NetaI) :: etaI
       Real(kind=dp), Dimension(3) :: a1i,a2i
       Real(kind=dp), Dimension(3) :: ai 
       Real(kind=dp), Dimension(9) :: etaN
       Real(kind=dp) :: Delta
       Real(kind=dp) ::  a1,a2
       Real(kind=dp), parameter ::  UnTier = 0.3333333333333333333333333333333333333333_dp  
       Integer, Dimension(3) :: ordre
       Integer :: i,j,n
       Integer :: ik1,ik2
       Integer :: N4,N5,N6
              
       
       Delta = UnTier / Ndiv
!
! tri des ki 
! pour avoir k1 < k2 < k3 
!
!
        Call triki(ai0,ai,eta6,ordre,1)
!
        a1 = ai(1)
        a2 = ai(2)
!
! Position de a1,a2 dans EtaI
! calcul des indices ik1,ik2      
! ik1,ik2 indice du noeud en haut a droite du pt (a1,a2) 
!
         ik1 = Int((a1 + Delta)/Delta) + 1
         ik2 = Int((a2 + Delta)/Delta) + 1

! Si ik1 + 2ik2 -3(Ndiv+1) >0 => on est sur la frontiere a2=a3
! ( a1+2a2=1)
!  => ik1=ik1-1 sauf si ik1=2 , ik2=ik2-1         
!         
         If ((ik1+2*ik2-3*(Ndiv+1)).GE.0) Then
          If ((ik1.NE.2).And.((ik1+2*ik2-3*(Ndiv+1)).NE.0).And.  &
          (abs((a1-Delta*(ik1-1))/a1).GT.1.0E-5))  ik1=ik1-1
          ik2=ik2-1 
         End If
         If (ik1.EQ.1) ik1=2
!
! Indice N4,N5 et N6 des points 4,5 et 6 dans EtaI
!  
         
         N4 = Nk2(ik1-1) + ik2 - ik1 + 3         
         N5 = Nk2(ik1) + ik2 - ik1 + 2         
         N6 = Nk2(ik1+1) + ik2 - ik1 + 1         

!
! Remplissage etaN(1 a 9)
!  7 - 8 - 9  
!  4 - 5 - 6 
!  1 - 2 - 3
!
       Do i=1,3
       a1i(i)=Delta*(ik1-3+i)
       a2i(i)=Delta*(ik2-3+i)
       End Do
!
       Do n=1,6
         Do i=1,3              
           etaN(1+(i-1)*3) = etaI(6*(N4-3+i)+n)
           etaN(2+(i-1)*3) = etaI(6*(N5-3+i)+n)
           etaN(3+(i-1)*3) = etaI(6*(N6-3+i)+n)
         End Do
!
! interpolation sur le Q9  
!
!        If ( (a1 < a1i(1)) .OR. &
!             (a1 > a1i(3)) .OR. &
!             (a2 < a2i(1)) .OR. &
!             (a2 > a2i(3)) ) Then
!           write(*,*)a1,a1i
!           write(*,*)a2,a2i
!           Stop
!         End If
            
         eta6(n)=InterQ9(a1,a2,a1i,a2i,etaN)

       End Do
!
! tri des eta
!
        Call triki(ai0,ai,eta6,ordre,-1)

       Return 
       End

!**************************************************************************************
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute Eigenvalues (ai) and euler angles (Euler)
!  of the second order orientation tensor (a2)
!   using lapack DGEEV lapack routine
!--------------------------------------------------------
      subroutine R2Ro(a2,dim,ai,Euler)

      use header    ! types d'Elmer

      implicit none

      Real(dp),dimension(6),intent(in) :: a2
      Real(dp),dimension(3),intent(out) :: ai,Euler
      Real(dp),dimension(3,3) :: A,EigenVec
      Real(dp) :: Dumy(1,3),EI(3),Work(24)
      Real(dp) :: norm
      integer :: dim               ! dimension  (2D-3D)
      integer :: i,infor 

      Do i=1,3
         A(i,i)=a2(i)
      End do
         A(1,2)=a2(4)
         A(2,1)=A(1,2)
         A(2,3)=a2(5)
         A(3,2)=A(2,3)
         A(1,3)=a2(6)
         A(3,1)=A(1,3)
      

      CALL DGEEV('N','V',3,A,3,ai,EI,Dumy,1,EigenVec,3,Work,24,infor)

     ! need a right handed orthonormal basis to compute euler angles;
     ! not guarantee by DGEEV.
      EigenVec(1,3)=EigenVec(2,1)*EigenVec(3,2)-EigenVec(3,1)*EigenVec(2,2)
      EigenVec(2,3)=EigenVec(3,1)*EigenVec(1,2)-EigenVec(1,1)*EigenVec(3,2)
      EigenVec(3,3)=EigenVec(1,1)*EigenVec(2,2)-EigenVec(2,1)*EigenVec(1,2)

     ! normalize
      norm=sqrt(EigenVec(1,3)**2+EigenVec(2,3)**2+EigenVec(3,3)**2)
      EigenVec(:,3)=EigenVec(:,3)/norm

      Euler(2)=Acos(EigenVec(3,3))
      if (abs(Euler(2)).gt.tiny(Euler(2))) then !3D euler angles 
        Euler(1)=ATAN2(EigenVec(1,3),-EigenVec(2,3))
        Euler(3)=ATAN2(EigenVec(3,1),EigenVec(3,2))
      else ! only one rotation of angle phi
        Euler(3)=0.0
        Euler(1)=ATAN2(EigenVec(2,1),EigenVec(1,1))
      end if

      RETURN
      END

      End Module Golf
