! *************************************************************
!
! Purpose: Check convergence to steady state for the 2D NS Solver
! 
! *************************************************************

SUBROUTINE CalcConvergence(u,v,crit,critx,crity)
USE Parameters

  IMPLICIT NONE

  INTEGER ::i,j
  DOUBLE PRECISION,INTENT(IN), DIMENSION(2,Nx+1,Ny+1) :: u,v
  DOUBLE PRECISION, INTENT(INOUT) :: crit,critx,crity
  DOUBLE PRECISION :: critxsum,critysum

! *************************************************************
! Sum the error in u and v for n+1 and n comparison
! *************************************************************

  critxsum = 0.0D0
  critysum = 0.0D0

  DO i = 3,Nx

     DO j = 2,Ny

         critxsum = critxsum + ( u(2,i,j) - u(1,i,j) )**2.0D0

     END DO 
  END DO

  critx = sqrt(critxsum) ! L2 norm

  DO i = 2,Nx
     DO j = 3,Ny
         critysum = critysum + (v(2,i,j)-v(1,i,j))**2.0D0
     END DO
  END DO
  crity = sqrt(critysum)

! *************************************************************
! Assign Maximum Error crit
! *************************************************************

  IF (critx > crity) THEN
     crit = critx
  ELSE IF (critx < crity) THEN
     crit = crity
  END IF

END SUBROUTINE CalcConvergence
