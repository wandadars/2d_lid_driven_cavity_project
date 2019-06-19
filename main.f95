! *************************************************************
!
! Purpose:  main program to solve the 2D Navier Stokes
!           flow problem in a lid-driven cavity flow
!
! *************************************************************

PROGRAM main
USE Parameters

IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:) :: u,v
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:) :: us,vs
DOUBLE PRECISION :: ts,crit

! *************************************************************
! Read User Inputs
! *************************************************************
CALL ReadInput( )

! ************************************************************
! Allocate Arrays from User Input
! ************************************************************
ALLOCATE( u(2,Nx+1,Ny+1) )
ALLOCATE( v(2,Nx+1,Ny+1) )

ALLOCATE( us(Nx+1,Ny+1) )
ALLOCATE( vs(Nx+1,Ny+1) )

 ts = 0.0D0 ! sets the initial time as 0
 crit = 1.0D0 ! initial setting

! *************************************************************
! Call Time Stepping Loop
! This loop calculates the initial comditions, the convective terms
! the diffusive terms, the u*, v*, the pressure poisson equation
! and updated the time
! *************************************************************
CALL TimeStep(u,v,crit,ts)

! *************************************************************
! Write Solution
! *************************************************************
WRITE(*,*) 'SOLUTION FINISHED, Proceed to write step'

 CALL WriteSolution(Nx,Ny,deltax,deltay,Re,u,v)

DEALLOCATE( u )
DEALLOCATE( v )
DEALLOCATE( us )
DEALLOCATE( vs )

 END PROGRAM main
