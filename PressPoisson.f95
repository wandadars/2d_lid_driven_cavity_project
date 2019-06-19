! *************************************************************
! 
! Purpose: Solve the Pressure Poisson Equation for the 2D NS equations
! 
! *************************************************************


SUBROUTINE PressPoisson(us,vs,un,vn)
USE Parameters

IMPLICIT NONE

INTEGER :: i,j

DOUBLE PRECISION,INTENT(INOUT), DIMENSION(Nx+1,Ny+1) :: us,vs,un,vn
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:) :: P  ! Pressure matrix

DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:) :: x    ! Solution vector ( used in solver, then converted to pressure array)
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:) :: A  ! Coefficient Matrix
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:) :: b    ! RHS Vector

DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:) :: AResized
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:) :: bResized

INTEGER :: counter      ! Used to navigate coefficient matrix
INTEGER :: vertJump     ! Used to move vertically in pressure array
INTEGER :: Error        ! Used to store error code from Solver

DOUBLE PRECISION ::  DeltaxFactor, DeltayFactor     ! Coefficient matrix factors
INTEGER :: SysSize   ! Size of system of equations being solved ( recall striking of rows and columns)

SysSize = (Nx-1)*(Ny-1)-1

ALLOCATE( x(SysSize+1) )
ALLOCATE( A( SysSize+1 ,SysSize+1) )
ALLOCATE( b(SysSize+1) )

ALLOCATE( AResized(SysSize,SysSize) )
ALLOCATE( bResized(SysSize) )
ALLOCATE( P(Nx,Ny) ) 

! *************************************************************
! Use the input u* and v* as us and vs
! Set the output as un and vn
! This function computes the pressure field, P,  necessary to make the input velocity
! field divergence free
! *************************************************************

! The active elements of the pressure array are within the following ranges: i =2, Nx-1, j =2, Ny-1

WRITE(*,*) 'Beginning Pressure Poisson Solver'

DeltaxFactor = 1.0/(deltax**2)        ! This is simply 1/(deltax^2)
DeltayFactor = 1.0/(deltay**2)        ! This is simply 1/(deltay^2)

! Write out results to screen for checking
!WRITE(*,*) Deltax, Deltay, DeltaxFactor, DeltayFactor


! Initialize the entire pressure array to zero
DO i = 1,Nx

        Do j = 1,Nx
                P(i,j) = 0       
         ENDDO
ENDDO

!--------Begin Coefficient Matrix Filling-----------

! NOTE: This algorithm fills a coefficient matrix assuming the 2D pressure array
! has been turned into a vector where the variables in the vector move down as
! the physical location of the pressure variable moves to the right. i.e. the
! array is stored in a vector using ROW-MAJOR order:
! (http://en.wikipedia.org/wiki/Row-major_order)


! Initialize arrays and vectors to zero

DO i = 1,(Nx-1)*(Ny-1)

        DO j = 1,(Nx-1)*(Ny-1)
                
                A(i,j) = 0

       ENDDO

        b(i) = 0
        x(i) = 0

ENDDO


vertjump = Nx-1 ! Number of spaces right to move to represent movement up


! The idea here is to step throught the entire pressure array and fill the
! coefficient array depending on the position within the array using IF
! statements.
counter = 1     ! As this variable increments it corresponds to moving right in the physical domain and then up when a wall is reached

DO j = 2,Ny

        DO i = 2,Nx
                
                ! Fill INTERNAL NODES
                IF(i >= 3 .AND. i <Nx .AND. j >= 3 .AND. j < Ny ) THEN

                        ! Derivative with respect to x
                        A(counter,counter+1) = 1*DeltaxFactor
                        A(counter,counter) = A(counter,counter) -2*DeltaxFactor
                        A(counter, counter-1) = 1*DeltaxFactor

                        ! Derivative with respect to y
                        A(counter,counter+vertjump) = 1*DeltayFactor
                        A(counter,counter) = A(counter,counter) -2*DeltayFactor
                        A(counter,counter-vertjump) = 1*DeltayFactor
                        
                        ! RHS vector
                         b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay  )
                        
                ENDIF


                ! Fill Array for BOTTOM INTERNAL BOUNDARY
                IF (i > 2 .AND. i < Nx .AND. j == 2) THEN
        
                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter) -2*DeltaxFactor
                        A(counter,counter+1) = 1*DeltaxFactor
                        A(counter,counter-1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) -1*DeltayFactor
                        A(counter,counter+vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay  ) -vs(i,j)/(deltat*Deltay)


                ENDIF

                ! Fill Coefficient Array for LEFT INTERNAL BOUNDARY

                 IF (i == 2 .AND. j > 2 .AND. j<Ny) THEN

                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter) -1*DeltaxFactor
                        A(counter,counter+1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) - 2*DeltayFactor
                        A(counter,counter+vertjump) = 1*DeltayFactor
                        A(counter,counter-vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay )-us(i,j)/(deltat*Deltax)

                ENDIF

                ! Fill Array for RIGHT INTERNAL BOUNDARY

                IF (i == Nx .AND. j > 2 .AND. j<Ny) THEN

                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter) -1*DeltaxFactor
                        A(counter,counter-1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) - 2*DeltayFactor
                        A(counter,counter+vertjump) = 1*DeltayFactor
                        A(counter,counter-vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay )+us(i,j)/(deltat*Deltax)

                ENDIF

                

                ! Fill Array for TOP INTERNAL BOUNDARY
                IF (i >2 .AND. i <Nx .AND. j == Ny) THEN

                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter) -2*DeltaxFactor
                        A(counter,counter+1) = 1*DeltaxFactor
                        A(counter,counter-1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) -1*DeltayFactor
                        A(counter,counter-vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay )+ vs(i,j)/(deltat*Deltay)

                ENDIF


                ! Fill BOTTOM LEFT CORNER NODE
                IF ( i == 2 .AND. j ==2 ) THEN

                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter)-1*DeltaxFactor
                        A(counter,counter+1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) -1*DeltayFactor
                        A(counter,counter+vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay )-us(i,j)/(deltat*Deltax) &
                        - vs(i,j)/(deltat*Deltay)


                ENDIF


                ! Fill Array for BOTTOM RIGHT CORNER

                 IF (i == Nx .AND. j== 2) THEN

                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter) -1*DeltaxFactor
                        A(counter,counter-1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) -1*DeltayFactor
                        A(counter,counter+vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay )+us(i,j)/(deltat*Deltax) &
                        - vs(i,j)/(deltat*Deltay)

                ENDIF



                ! Fill Array for TOP LEFT CORNER

                 IF (i == 2 .AND. j== Ny) THEN

                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter) -1*DeltaxFactor
                        A(counter,counter+1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) -1*DeltayFactor
                        A(counter,counter-vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay )-us(i,j)/(deltat*Deltax) &
                        + vs(i,j)/(deltat*Deltay)


                ENDIF



                ! Fill Array for TOP RIGHT CORNER

                IF (i == Nx .AND. j== Ny) THEN

                        ! Partial derivative with respect to x
                        A(counter,counter) = A(counter,counter) -1*DeltaxFactor
                        A(counter,counter-1) = 1*DeltaxFactor

                        ! Partial derivative with respect to y
                        A(counter,counter) = A(counter,counter) -1*DeltayFactor
                        A(counter,counter-vertjump) = 1*DeltayFactor

                        ! RHS Vector
                        b(counter)=1/deltat*((us(i+1,j)-us(i,j) )/Deltax+ &
                        (vs(i,j+1)-vs(i,j) )/Deltay )+us(i,j)/(deltat*Deltax) &
                        + vs(i,j)/(Deltay*Deltay)

                ENDIF


                ! Increment counter i.e. move to next node
                counter = counter + 1
        ENDDO

ENDDO

!-------------END MATRIX FILLING------------------------------

WRITE(*,*)'Pressure Coefficient Matrix Filling Complete'
! We have put together the coefficient matrix and the RHS vector. This problem
! is underdefined because we have not given any boundary conditions for
! pressure. Also recall that this is not really pressure; it is more of a
! correction term. With that in mind, we just fix a value for pressure somewhere
! in the domain and then the system of equations will be properly defined and
! solvable.  We pick P(2,2,) = 0. This choice means that we remove the first row
! and first column of the coefficient array and the first row of the RHS vector.

P(2,2) = 0

! Strike out first row and first column
DO j = 2, (Nx-1)*(Ny-1)

        DO i = 2, (Nx-1)*(Ny-1) 

                AResized(i-1,j-1) = A(i,j)

        ENDDO

        bResized(j-1)=b(j)

ENDDO


! At this point we are ready to solve for the pressure distribution that will
! make the input velocity field divergence free

WRITE(*,*) 'Calling Pressure Solver'

CALL cg_method(SysSize,AResized,bResized,x, Error )

WRITE(*,*) ' Solver Error Code: ', Error

! We now have the solution for the pressure field that will give a divergence
! free velocity field at the next timestep. The only issue is that the solution
! is in the form of a vector, and it would be much easier to use an array.

! Convert solution vector to array
counter = 1

DO j = 2,Ny

        DO i = 2,Nx
        
                IF( i == 2 .AND. j == 2) THEN
                        P(i,j) = 0
                ELSE
                        P(i,j) = x(counter)
                        counter = counter +1

                ENDIF

        ENDDO

ENDDO


!------------------------------- Solve for un-------------------------
DO j = 2,Ny

        DO i = 3,Nx

                un(i,j) = us(i,j) - (deltat/Deltax)*(P(i,j)-P(i-1,j))

               ! TroubleShoot Pressure Gradient Issues
               ! WRITE(*,*) P(i,j)-P(i-1,j)

        ENDDO

ENDDO
!--------------------------- END SOLVE FOR un-------------------------



!------------------------------- Solve for vn-------------------------
DO j = 3,Ny

        DO i = 2,Nx
                vn(i,j)=vs(i,j) - (deltat/Deltay)*(P(i,j)-P(i,j-1))
        ENDDO

ENDDO

!--------------------------- END SOLVE FOR vn-------------------------



! Error Checking
OPEN( UNIT=1, FILE='Coef_Matrix_And_RHS_And_Other_Data.dat', STATUS = 'REPLACE')


WRITE(1,'(1X,/, A, /)') 'Coefficient Matrix AResized Elements:'
WRITE(1,'(1X,/,A,I8, /)')'Array is Square with dimension', counter-1

DO i =1,counter-1

        WRITE(1,'(5000(1X,f10.3, 4X))'),(AResized(i,j),j=1,counter-1)
        WRITE(1,*)


ENDDO


WRITE(1,'(1X,/, A, /)') 'U Velocity Matrix Us Elements:'

DO i =2,Nx+1

        WRITE(1,'(5000(1X,f10.3, 4X))'), (us(i,j),j=2,Ny+1)
        WRITE(1,*)

ENDDO


WRITE(1,'(1X,/, A, /)') 'V Velocity Matrix Vs Elements:'

DO i =2,Nx+1

        WRITE(1,'(5000(1X,f10.3, 4X))'), (vs(i,j),j=2,Ny+1)
        WRITE(1,*)

ENDDO


WRITE(1,'(1X,/, A, /)') 'Resized RHS Vector Elements:'
DO j = 1,counter-1

        WRITE(1,'(1X,f10.3, 4X,/)') bResized(j)

ENDDO


WRITE(1,'(1X,/, A, /)') 'Pressure Matrix P Elements:'
WRITE(1,'(1X,/,A,I8, /)')'Array is Square with dimension', SysSize+1
DO i =2,Nx

        WRITE(1,'(5000(1X,f8.3, 4X))'),(P(i,j),j=2,Ny)

        WRITE(1,*)


ENDDO







WRITE(1,'(1X,/, A, /)') 'U Velocity Matrix U(n+1) Elements:'
DO i =2,Nx+1

        WRITE(1,'(5000(1X,f10.3, 4X))'), (un(i,j),j=2,Nx+1)

        WRITE(1,*)

ENDDO



WRITE(1,'(1X,/, A, /)') 'V Velocity Matrix V(n+1) Elements:'
DO i =2,Nx+1

        WRITE(1,'(5000(1X,f10.3, 4X))'), (vs(i,j),j=2,Ny+1)

        WRITE(1,*)

ENDDO



! Close Output File
CLOSE(1)


DEALLOCATE(x)
DEALLOCATE(A)
DEALLOCATE(b)
DEALLOCATE(P)
DEALLOCATE(AResized)
DEALLOCATE(bResized)

WRITE(*,*) 'Pressure Poisson Solver Completed'


END SUBROUTINE PressPoisson
