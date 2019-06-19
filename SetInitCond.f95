! *************************************************************
!
! Purpose: Set the initial conditions for the cavity driven flow problem
!
! *************************************************************

SUBROUTINE SetInitCond(u,v,us,vs,un,vn)
USE Parameters

  IMPLICIT NONE

  INTEGER :: i,j  
  DOUBLE PRECISION, INTENT(INOUT),  DIMENSION(2,Nx+1,Ny+1) :: u,v
  DOUBLE PRECISION, INTENT(INOUT), DIMENSION(Nx+1,Ny+1) :: us,vs,un,vn

! *************************************************************
! Define the bounds of the space
! *************************************************************
  deltax = xL/(Nx-1)
  deltay = yL/(Ny-1)

! *************************************************************
! Set all points 1 for n and n-1 time
! *************************************************************
! index for u and v are u(n,x,y) 
! where u(1,x,y) is the n-1 time and u(2,x,y) is the n time
DO i = 1, Nx+1

  DO j = 1,Ny+1
   
    u(1,i,j)=0.0D0
    us(i,j) = 0.0D0
    un(i,j) = 0.0D0

    u(2,i,j)=0.0D0

    v(1,i,1:j)=0.0D0
    vs(i,j)=0.0D0
    vn(i,j) = 0.0D0

    v(2,i,1:j)=0.0D0

  ENDDO

ENDDO

! *************************************************************
! Set the top row to 1
! *************************************************************
DO i = 1,Nx+1

  u(:,i,Ny+1) = 1.0D0
  us(i,Ny+1) = 1.0D0
  un(i,Ny+1) = 1.0D0 


ENDDO
! Print Initial Conditions

OPEN(UNIT = 1, FILE ='Initial_Conditions.dat',STATUS = 'REPLACE')

WRITE(1,'(1X,/, A, /)') 'Initial u velocity Elements:'
DO i =1,Nx+1

        WRITE(1,'(5000(1X,f10.3, 4X))'),(u(1,i,j),j=1,Nx+1)

        WRITE(1,*)

ENDDO

CLOSE(1)


 END SUBROUTINE SetInitCond
