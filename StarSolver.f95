! *************************************************************
!
! Purpose: Calculates coefficients for the diffusive term H, sets up
! the coefficients for the Thomas Algorithm
!
! *************************************************************

SUBROUTINE StarSolver(u,v,us,vs)
USE Parameters
IMPLICIT NONE

INTEGER :: i, j
DOUBLE PRECISION,INTENT(IN), DIMENSION(2,Nx+1,Ny+1) :: u, v
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:) :: A
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(Nx+1,Ny+1) :: us, vs
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: b, x
DOUBLE PRECISION, DIMENSION(2) :: vAvg, uAvg
DOUBLE PRECISION :: SecT, ThirdT, FourthT
DOUBLE PRECISION ::gamx, gamy, deltx, delty
INTEGER :: SysSize, pos, Error, vertjump


SysSize = ( (Nx-3)+1 )*( (Ny-2)+1) ! # of eqns to solve

! ALLOCATE VARIABLES
ALLOCATE( A(SysSize,SysSize) )
ALLOCATE( b(SysSize) )
ALLOCATE( x(SysSize) )
  
WRITE(*,*) 'Starting the Calculation for U Star'

! *************************************************************
! Establish constants
! *************************************************************
gamx = deltat/(Re*deltax**2.0D0)
gamy = deltat/(Re*deltay**2.0D0)
deltx = deltat/deltax
delty = deltat/deltay
vertjump = Nx-2
!*************************************************************

!************ Initialize Arrays to 0**************************

DO i = 1,SysSize

  DO j = 1, SysSize
    A(i,j) = 0
  ENDDO

  b(i) = 0
  x(i) = 0

ENDDO
!**************************************************************


!************************ Solve for U* *************************


! Construct Coefficient Matrix and RHS
pos = 1

! Step through every element of the domain
DO j = 2, Ny

  DO i = 3,Nx

    ! BOTTOM BOUNDARY
    IF( j == 2 ) THEN

      A(pos,pos) = 1 + gamx + 2*gamy 

      IF( i == 3 ) THEN ! Left Boundary
         A(pos,pos+1) = -0.5*gamx
         b(pos) = 0.5*gamx*us(i-1,j)
      ELSE IF( i == Nx ) THEN ! Right Boundary
         A(pos,pos-1) = -0.5*gamx
         b(pos) = 0.5*gamx*us(i+1,j)
      ELSE
         A(pos,pos+1) = -0.5*gamx
         A(pos,pos-1) = -0.5*gamx
      ENDIF 

      A(pos,pos+vertjump) = -(2.0/3.0)*gamy
   
      vAvg(2) = 0.25*( v(2,i,j)+v(2,i-1,j)+v(2,i,j+1)+v(2,i-1,j+1) )
      vAvg(1) = 0.25*( v(1,i,j)+v(1,i-1,j)+v(1,i,j+1)+v(1,i-1,j+1) )

      SecT = -0.75*deltx*u(2,i,j)*( u(2,i+1,j)-u(2,i-1,j))- &
              3*delty*vAvg(2)*( -u(2,i,j-1)+u(2,i,j) )

      ThirdT =  0.25*deltx*u(1,i,j)*(u(1,i+1,j)-u(1,i-1,j)) + &
                delty*vAvg(1)*( -u(1,i,j-1)+u(1,i,j) )

      FourthT = 0.5*gamx*( u(2,i+1,j)-2*u(2,i,j)+u(2,i-1,j) ) &
                +(2.0/3.0)*gamy*(2*u(2,i,j-1)-3*u(2,i,j)+u(2,i,j+1))

      b(pos) = b(pos)+ u(2,i,j)+SecT+ThirdT+FourthT-(4.0/3.0)*gamy*us(i,j-1)
      
      pos = pos+1
     
    ENDIF      

    ! INTERNAL NODES
    IF( j > 2 .AND. j < Ny ) THEN

      A(pos,pos) = 1 + gamx + gamy

      IF( i == 3 ) THEN  ! LEFT BOUNDARY

         A(pos,pos+1) = -0.5*gamx
         b(pos) = 0.5*gamx*us(i-1,j)

      ELSE IF( i == Nx ) THEN ! RIGHT BOUNDARY

         A(pos,pos-1) = -0.5*gamx
         b(pos) = 0.5*gamx*us(i+1,j)

      ELSE
         A(pos,pos+1) = -0.5*gamx
         A(pos,pos-1) = -0.5*gamx

      ENDIF

      A(pos,pos+vertjump) = -0.5*gamy
      A(pos,pos-vertjump) = -0.5*gamy

      vAvg(2) = 0.25*( v(2,i,j)+v(2,i-1,j)+v(2,i,j+1)+v(2,i-1,j+1) )
      vAvg(1) = 0.25*( v(1,i,j)+v(1,i-1,j)+v(1,i,j+1)+v(1,i-1,j+1) )

      SecT = 0.75*deltx*u(2,i,j)*(u(2,i+1,j)-u(2,i-1,j)) - &
             0.75*delty*vAvg(2)*(u(2,i,j+1)-u(2,i,j-1) )

      ThirdT =  0.25*deltx*u(1,i,j)*(u(1,i+1,j)-u(1,i-1,j))+ &
                0.25*delty*vAvg(1)*(u(1,i,j+1)-u(1,i,j-1) )  

      FourthT = 0.5*gamx*( u(2,i+1,j)-2*u(2,i,j)+u(2,i-1,j) )+ &
                0.5* gamy*(u(2,i,j+1)-2*u(2,i,j)+u(2,i,j-1) )

      b(pos) = b(pos) + u(2,i,j) +SecT+ThirdT+FourthT

      pos = pos+1


    ENDIF

    ! TOP BOUNDARY
    IF( j == Ny ) THEN

      A(pos,pos) = 1 + gamx + 2*gamy
      
      IF( i == 3 ) THEN ! LEFT Boundary

         A(pos,pos+1) = -0.5*gamx
         b(pos) = 0.5*gamx*us(i-1,j)

      ELSE IF( i == Nx ) THEN ! Right Boundary

         A(pos,pos-1) = -0.5*gamx
         b(pos) = 0.5*gamx*us(i+1,j)

      ELSE ! Internal
         A(pos,pos+1) = -0.5*gamx
         A(pos,pos-1) = -0.5*gamx

      ENDIF

      A(pos,pos-vertjump) = -(2.0/3.0)*gamy

      vAvg(2) = 0.25*( v(2,i,j)+v(2,i-1,j)+v(2,i,j+1)+v(2,i-1,j+1) )
      vAvg(1) = 0.25*( v(1,i,j)+v(1,i-1,j)+v(1,i,j+1)+v(1,i-1,j+1) )

      SecT = -0.75*deltx*u(2,i,j)*( u(2,i+1,j)-u(2,i-1,j)) &
             - 3.0*delty*vAvg(2)*(u(2,i,j+1)-u(2,i,j) )

      ThirdT =  0.25*deltx*u(1,i,j)*(u(1,i+1,j)-u(1,i-1,j))+ &
                delty*vAvg(1)*(u(1,i,j+1)-u(1,i,j) )

      FourthT = 0.5*gamx*( u(2,i+1,j)-2*u(2,i,j)+u(2,i-1,j)  )+ &
                (2.0/3.0)*gamy*(2*u(2,i,j+1)- 3*u(2,i,j)+u(2,i,j-1) )

      b(pos) = b(pos) + u(2,i,j) +SecT+ThirdT+FourthT+(4.0/3.0)*gamy*us(i,j+1)

      pos = pos+1


    ENDIF


  ENDDO


ENDDO


CALL cg_method(SysSize, A, b, x, Error )


! Now that we have the solution for U*, reconstruct 2D array of data

pos = 1

DO j = 2,Ny

  DO i = 3,Nx

        us(i,j) = x(pos)
        pos = pos + 1
  ENDDO

ENDDO



! Error Checking
OPEN( UNIT=1, FILE='U_and_V_Star_Data.dat', STATUS = 'REPLACE')

WRITE(1,'(1X,/, A, /)') 'U* Coefficient Matrix Elements:'
DO i =1, SysSize

        WRITE(1,'(5000(1X,f6.3,2X))'),(A(i,j),j=1,SysSize)
        WRITE(1,*)

ENDDO


WRITE(1,'(1X,/, A, /)') 'U* Elements:'
DO i =3, Nx

        WRITE(1,'(5000(1X,f10.3, 4X))'),(us(i,j),j=2,Ny)
        WRITE(1,*)

ENDDO



!********************  END SOLUTION FOR Us i.e. U*  ***************

WRITE(*,*) 'Starting the Calculation for V Star'


!************ Re-Initialize Arrays******************************

DO i = 1,SysSize

  DO j = 1, SysSize

    A(i,j) = 0

  ENDDO

  b(i) = 0
  x(i) = 0

ENDDO

!**************************************************************




!************************ Solve for V* **************************
vertjump = Nx-1

pos = 1

! Step through every element of the domain
DO j = 3, Ny

  DO i = 2,Nx

    ! LEFT BOUNDARY
    IF( i == 2 ) THEN

      A(pos,pos) = 1 + 2*gamx + gamy

      IF( j == 3 ) THEN ! Bottom Boundary

         A(pos,pos+vertjump) = -0.5*gamy
         b(pos) = 0.5*gamy*vs(i,j-1)

      ELSE IF( j == Ny ) THEN ! Top Boundary

         A(pos,pos-vertjump) = -0.5*gamy
         b(pos) = 0.5*gamy*vs(i,j+1)

      ELSE
         A(pos,pos+vertjump) = -0.5*gamy
         A(pos,pos-vertjump) = -0.5*gamy

      ENDIF

      A(pos,pos+1) = -(2.0/3.0)*gamx

      uAvg(2) = 0.25*( u(2,i,j)+u(2,i-1,j)+u(2,i,j+1)+u(2,i-1,j+1) )
      uAvg(1) = 0.25*( u(1,i,j)+u(1,i-1,j)+u(1,i,j+1)+u(1,i-1,j+1) )

      SecT = -3*deltx*uAvg(2)*(-v(2,i-1,j)+v(2,i,j) )- &
              0.75*delty*v(2,i,j)*(v(2,i,j+1)-v(2,i,j-1) )

      ThirdT = deltx*uAvg(1)*(-v(1,i-1,j)+v(1,i,j))+ &
               0.25*delty*v(1,i,j)*(v(1,i,j+1)-v(1,i,j-1))

      FourthT = (2.0/3.0)*gamx*(2*v(2,i-1,j)-3*v(2,i,j)+v(2,i+1,j) ) &
                + 0.5*gamy*(v(2,i,j+1)-2*v(2,i,j)+v(2,i,j-1) )

      b(pos) = b(pos)+v(2,i,j)+SecT+ThirdT+FourthT+(4.0/3.0)*gamx*vs(i-1,j)

      pos = pos+1

    ENDIF

    ! INTERNAL NODES
    IF( i > 2 .AND. i < Nx ) THEN

      A(pos,pos) = 1+ gamx + gamy

      IF( j == 3 ) THEN  ! BOTTOM BOUNDARY

         A(pos,pos+vertjump) = -0.5*gamy
         b(pos) = 0.5*gamy*vs(i,j-1)

      ELSE IF( j == Ny ) THEN ! TOP BOUNDARY

         A(pos,pos-vertjump) = -0.5*gamy
         b(pos) = 0.5*gamy*vs(i,j+1)

      ELSE
         A(pos,pos+vertjump) = -0.5*gamy
         A(pos,pos-vertjump) = -0.5*gamy

      ENDIF

      A(pos,pos+1) = -0.5*gamx
      A(pos,pos-1) = -0.5*gamx

      uAvg(2) = 0.25*( u(2,i,j)+u(2,i-1,j)+u(2,i,j+1)+u(2,i-1,j+1) )
      uAvg(1) = 0.25*( u(1,i,j)+u(1,i-1,j)+u(1,i,j+1)+u(1,i-1,j+1) )

      SecT = -0.75*deltx*uAvg(2)*(v(2,i+1,j)-v(2,i-1,j))- &
             0.75*delty*v(2,i,j)*( v(2,i,j+1)-v(2,i,j-1) )

      ThirdT = 0.25*deltx*uAvg(1)*(v(1,i+1,j)-v(1,i-1,j))+ &
               0.25*delty*v(1,i,j)*(v(1,i,j+1)-v(1,i,j-1))

      FourthT = 0.5*gamx*(v(2,i+1,j)-2*v(2,i,j)+v(2,i-1,j) ) &
                + 0.5*gamy*(v(2,i,j+1)-2*v(2,i,j)+v(2,i,j-1) )

      b(pos) = b(pos)+ v(2,i,j) +SecT+ThirdT+FourthT

      pos = pos+1


    ENDIF

    ! RIGHT BOUNDARY
    IF( i == Nx ) THEN

      A(pos,pos) = 1 +2*gamx + gamy

      IF( j == 3 ) THEN ! Bottom Boundary

         A(pos,pos+vertjump) = -0.5*gamy
         b(pos) = 0.5*gamy*vs(i,j-1)

      ELSE IF( j == Ny ) THEN ! Top Boundary

         A(pos,pos-vertjump) = -0.5*gamy
         b(pos) = 0.5*gamy*vs(i,j+1)

      ELSE
         A(pos,pos+vertjump) = -0.5*gamy
         A(pos,pos-vertjump) = -0.5*gamy

      ENDIF

      A(pos,pos-1) = -(2.0/3.0)*gamx

      uAvg(2) = 0.25*( u(2,i,j)+u(2,i-1,j)+u(2,i,j+1)+u(2,i-1,j+1) )
      uAvg(1) = 0.25*( u(1,i,j)+u(1,i-1,j)+u(1,i,j+1)+u(1,i-1,j+1) )

      SecT = -3*deltx*uAvg(2)*(v(2,i+1,j)-v(2,i,j) )- &
              0.75*delty*v(2,i,j)*(v(2,i,j+1)-v(2,i,j-1) )

      ThirdT = deltx*uAvg(1)*(v(1,i+1,j)-v(1,i,j))+ &
               0.25*delty*v(1,i,j)*(v(1,i,j+1)-v(1,i,j-1))

      FourthT = (2.0/3.0)*gamx*(2*v(2,i+1,j)-3*v(2,i,j)+v(2,i-1,j) ) &
                + 0.5*gamy*(v(2,i,j+1)-2*v(2,i,j)+v(2,i,j-1) )

      b(pos) = b(pos)+v(2,i,j)+SecT+ThirdT+FourthT+(4.0/3.0)*gamx*vs(i+1,j)

      pos = pos+1


    ENDIF


  ENDDO


ENDDO


CALL cg_method(SysSize, A, b, x, Error )

IF( Error /= 0 ) THEN
  WRITE(*,*) 'ERROR: Could not solve for V*. Inversion Method failed'
ENDIF

! Now that we have the solution for V*, reconstruct 2D array of data
pos = 1
DO j = 3,Ny

  DO i = 2,Nx

        vs(i,j) = x(pos)
        pos = pos + 1
  ENDDO

ENDDO



! Error Checking

WRITE(1,'(1X,/, A, /)') 'V* Coefficient Matrix Elements:'
DO i =1, SysSize

        WRITE(1,'(5000(1X,f6.3,2X))'),(A(i,j),j=1,SysSize)

        WRITE(1,*)


ENDDO




WRITE(1,'(1X,/, A, /)') 'V* Elements:'

DO i =2, Nx

        WRITE(1,'(5000(1X,f10.3, 4X))'),(vs(i,j),j=3,Ny)

        WRITE(1,*)


ENDDO

CLOSE(1)

!********************  END SOLUTION FOR Us i.e. V*  ***************


! DEALLOCATE ARRAYS
DEALLOCATE(A)
DEALLOCATE(b)
DEALLOCATE(x)



!**************************************************************

WRITE(*,*) 'Finished Solving for U and V Star'

  END SUBROUTINE StarSolver

