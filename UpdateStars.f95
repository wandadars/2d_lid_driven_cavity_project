! *************************************************************
!
! Purpose: To Update the values of the intermediate velocity
!          field boundary values
!
! *************************************************************

SUBROUTINE UpdateStars(P,us,vs)
USE Parameters
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN), DIMENSION(Nx,Ny) :: P
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(Nx+1,Ny+1) :: us,vs

! Here we use the pressure field to update what the boundary values
! for the Star velocities are using extrapolation.



END SUBROUTINE UpdateStars
