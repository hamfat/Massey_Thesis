!---------------------------------------------------------------------- 
! "dimless.f" 
!---------------------------------------------------------------------- 
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION V, N, v1, vca, vL, vK, v2, v4, gca, gL, gK, phin
      DOUBLE PRECISION  minf, v3, ninf, C, psi
      
      V = U(1)
      N = U(2)
! bifurcation parameters
      v1 = PAR(1)
	  v3 = PAR(2)
! model parameters
      vca = 80.d0;  
      vL = (-70.d0/vca);
      vK = (-90.d0/vca);
      v2 = (25.d0/vca);
	  !v3 = (-11.d0/vca);
      v4 = (14.5d0/vca);
	  gK = 3.1416E-13;
      gca = (1.57E-13/gk); 
      gL = ((7.854E-14)/gk);
	  C = 1.9635E-14;
      phin = 2.664d0;
	  psi=((C*phin)/gk)
!model equations and auxiliary functions
      minf = (0.5d0*(1+tanh((V-v1)/v2)));
      ninf = 0.5d0*(1+tanh((V-v3)/v4));
      F(1) = -(gL*(V-vL)+N*(V-vK)+gca*minf*(V-1.d0));
      F(2) = psi*cosh((V-v3)/(2.d0*v4))*(ninf-N);
      END SUBROUTINE FUNC
      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
    ! initial values for the bifurcation parameters 
      PAR(1)=-0.5d0
	  PAR(2)=-0.1375
    ! initial conditions
      U(1)=-0.18484d0
      U(2)=0.37227d0
      END SUBROUTINE STPNT
      SUBROUTINE BCND 
      END SUBROUTINE BCND
      SUBROUTINE ICND 
      END SUBROUTINE ICND
      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
      SUBROUTINE PVLS
      END SUBROUTINE PVLS
