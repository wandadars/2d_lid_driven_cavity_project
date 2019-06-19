! *************************************************************
! Purpose: Step through time for the 2D NS solver
!
! *************************************************************

SUBROUTINE TimeStep(u,v,crit,ts)
USE Parameters

IMPLICIT NONE

DOUBLE PRECISION, INTENT(INOUT),  DIMENSION(2,Nx+1,Ny+1) :: u,v
DOUBLE PRECISION, DIMENSION(Nx+1,Ny+1) :: us,vs,un,vn
DOUBLE PRECISION, INTENT(INOUT) :: ts,crit
DOUBLE PRECISION :: critx,crity

! *************************************************************
! Set the initial inputs if this is the first time step
! *************************************************************
IF (ts == 0.0D0) THEN
   CALL SetInitCond(u,v,us,vs,un,vn)
END IF

! *************************************************************
! Write out the convergence criteria with time step
! *************************************************************
OPEN(UNIT=10,FILE='convhist.dat',STATUS='REPLACE')
WRITE(10,*)'Solution Time,   Critx,   Crity,   Crit  '

! *************************************************************
! Loop here until steady state is reached
! *************************************************************
DO WHILE (crit > 1.0D-6)

! *************************************************************
! Calculate the Diffusive terms and get u* and v*
! (called us,vs)
! *************************************************************
CALL StarSolver(u,v,us,vs)

! *************************************************************
! Solve the pressure poisson equation to get u@(n+1) and v(n+1)
! *************************************************************
CALL PressPoisson(us,vs,un,vn)

! *************************************************************
! Update the velocity terms u and v
! *************************************************************
CALL UpdateVelocity(un,vn,u,v)

! *************************************************************
! Update the convergence criteria
! *************************************************************
CALL CalcConvergence(u,v,crit,critx,crity)
WRITE(*,*) 'time = ', ts
WRITE(*,*) 'convergence criteria = ', crit

! *************************************************************
! Write Convergence History
! *************************************************************
WRITE(10,'(1X,f8.3,3X,f8.3,3X,f8.3,3X,f8.3,/)') ts,critx,crity,crit
 
!**************************************************************
! Check for excessive iteration and stop if detected
! *************************************************************
IF (ts > 400.0D0*deltat) THEN
        crit = 1.0D-10
END IF

! *************************************************************
! Increment Time
! *************************************************************
ts = ts + deltat

END DO ! criteria is met to say steady state

CLOSE(10)
END SUBROUTINE TimeStep 
