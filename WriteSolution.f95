! *************************************************************
!
! Purpose: write u, v, at geometric centerline
!
! **************************************************************

SUBROUTINE WriteSolution(u,v)
USE Parameters

IMPLICIT NONE

INTEGER :: i,j,xcenter,ycenter
DOUBLE PRECISION, INTENT(IN), DIMENSION(2,Nx+1,Ny+1) :: u,v
DOUBLE PRECISION, DIMENSION(Nx+1,Ny+1) :: xp,yp,up,vp
DOUBLE PRECISION :: vcenter,ucenter,xc,yc

! **************************************************************
! Find u, v: requires interpolation
! **************************************************************
DO j = 1, Ny
  DO i = 1, Nx
     IF (j ==1) THEN
       yp(i,j) = 0.0D0
       up(i,j) = u(2,i,j)
     ELSE IF (j == Ny) THEN
       yp(i,j) = 1.0D0
       up(i,j) = u(2,i,j)
     ELSE
       yp(i,j) = yp(i,j-1) + deltay
       up(i,j) = (u(2,i,j+1) + u(2,i,j))/2.0D0 ! geometric avg
     END IF
  END DO
END DO



DO j = 1, Ny
  DO i = 1, Nx
    IF (i == 1) THEN
      xp(i,j) = 0.0D0
      vp(i,j) = v(2,i,j)
    ELSEIF (i == Ny) THEN
      xp(i,j) = 1.0D0
      vp(i,j) = v(2,i,j)
    ELSE
      xp(i,j) = xp(i-1,j) + deltax
      vp(i,j) = (v(2,i+1,j) + v(2,i,j))/2.0D0
    END IF
  END DO
END DO

! **************************************************************
! Find geometric center
! **************************************************************
xcenter = (Nx+1)/2.0
ycenter = (Ny+1)/2.0

! **************************************************************
! Vertical Centerline
! **************************************************************
OPEN(UNIT=2,FILE='vertcenter.dat',STATUS='REPLACE')
OPEN(UNIT=3,FILE='horizcenter.dat',STATUS='REPLACE')
OPEN(UNIT=4,FILE='velfield.dat',STATUS='REPLACE')

DO j=1,Ny
  WRITE(2,'(1x,2(1x,i5),5(1x,e23.12))') xcenter,j,xp(xcenter,j), &
                      yp(xcenter,j),up(xcenter,j),vp(xcenter,j)
END DO

DO i=1,Nx
  WRITE(3,'(1x,2(1x,i5),5(1x,e23.12))') i,ycenter,xp(i,ycenter), &
                      yp(i,ycenter),up(i,ycenter),vp(i,ycenter)
END DO

DO i = 1,Nx
  DO j=1,Ny
    WRITE(4,'(1x,2(1x,i5),5(1x,e23.12))') i,j,xp(i,j),yp(i,j),up(i,j),vp(i,j)
  END DO
END DO

CLOSE(2)
CLOSE(3)
CLOSE(4)

! **************************************************************
! WRITE FlOWFIELD DATA TO LEGACY FORMAT
! **************************************************************



END SUBROUTINE WriteSolution
