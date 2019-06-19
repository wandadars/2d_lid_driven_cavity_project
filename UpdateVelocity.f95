! *************************************************************
!
! Purpose: Updates the velocity found in the poisson solver
!
! **************************************************************

SUBROUTINE UpdateVelocity(un,vn,u,v)
USE Parameters

IMPLICIT NONE

INTEGER :: i,j
  
DOUBLE PRECISION,INTENT(INOUT), DIMENSION(2,Nx+1,Ny+1) :: u,v
DOUBLE PRECISION,INTENT(IN), DIMENSION(Nx+1,Ny+1) :: un,vn

! *************************************************************
! Update Velocities
! *************************************************************

! STEP 1: velocity at only time step n-1 is stored in u(1), v(1). 
! Move the information at time step n to n-1
! STEP 2: put the velocity found in poisson solver un, vn, as u(2) and
! v(2).

  DO i= 3,Nx
     DO j = 2,Ny
         u(1,i,j) = u(2,i,j)
         u(2,i,j) = un(i,j)
     END DO
  END DO

  DO i = 2,Nx
     DO j = 3,Ny
         v(1,i,j) = v(2,i,j)
         v(2,i,j) = vn(i,j)
     END DO
  END DO

END SUBROUTINE UpdateVelocity
