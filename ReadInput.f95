SUBROUTINE ReadInput( )
USE Parameters

IMPLICIT NONE

INTEGER:: unitNum

!*******************************************************
! Read Input Data From Input File: 2D_Input.inp
!**********************************************************
unitNum = 1
open(unit=unitNum, file='2D_Input.inp', status = 'old', action='read')

READ(unitNum,*) xL, yL,Nx,Ny,Re,deltat

!**********************************************************
! ECHO DATA BACK TO USE
!**********************************************************

WRITE(*,*) 'Length of X-direction of 2D rectangular domain'
WRITE(*,*) xL

WRITE(*,*) 'Length of Y-direction of 2D rectangular domain'
WRITE(*,*) yL

WRITE(*,*) 'Number of points in the x direction, Nx'
WRITE(*,*) Nx

WRITE(*,*) 'Number of points in the y direction, Ny'
WRITE(*,*) Ny

WRITE(*,*) 'Flow Reynolds number'
WRITE(*,*) Re

WRITE(*,*) 'Time step'
WRITE(*,*) deltat

END SUBROUTINE

