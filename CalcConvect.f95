! *************************************************************
!
! Purpose: Calculates the convection term
!
! *************************************************************

SUBROUTINE CalcConvect(u,v,Nx,Ny,deltax,deltay,H)

IMPLICIT NONE

INTEGER, INTENT(IN) :: Nx,Ny
INTEGER :: i,j

  DOUBLE PRECISION,INTENT(INOUT), DIMENSION(2,Nx+1,Ny+1) :: u,v,H
  DOUBLE PRECISION, INTENT(IN) :: deltax, deltay
  DOUBLE PRECISION, DIMENSION(2) :: vAvg, uAvg

  WRITE(*,*) 'Starting the Calculation for the Convective Terms'

! *************************************************************
! Build H values
! *************************************************************
! H(1,i,j) is the x term
! H(2,i,j) is the y term
! RECALL: u,v(1,i,j) is the n-1 term
!         u,v(2,i,j) is the n term
! Be careful not to confuse the first subscript in H with the first
! subscript in u,v

DO i = 3,Nx
 DO j = 2,Ny

   IF (j == Ny) THEN

    vAvg(2) = 0.25*(v(2,i,j)+v(2,i-1,j)+v(2,i,j+1)+v(2,i-1,j+1) )
    vAvg(1) = 0.25*(v(1,i,j)+v(1,i-1,j)+v(1,i,j+1)+v(1,i-1,j+1) )

    H(1,i,j) = -3.0D0/2.0D0*(u(2,i,j)*((u(2,i+1,j)-u(2,i,j))/deltax) + &
                           vAvg(2)*((u(2,i,j+1)-u(2,i,j))/deltay/2.0D0)) + &
              1.0D0/2.0D0*(u(1,i,j)*((u(1,i+1,j)-u(1,i,j))/deltax) + &
                           vAvg(1)*((u(1,i,j+1)-u(1,i,j))/deltay/2.0D0))
  ELSE

    vAvg(2) = 0.25*(v(2,i,j)+v(2,i-1,j)+v(2,i,j+1)+v(2,i-1,j+1) )
    vAvg(1) = 0.25*(v(1,i,j)+v(1,i-1,j)+v(1,i,j+1)+v(1,i-1,j+1) )
    WRITE(*,*) vAvg(1), vAvg(2)
    H(1,i,j) = -3.0D0/2.0D0*(u(2,i,j)*((u(2,i+1,j)-u(2,i,j))/deltax) + &
                           vAvg(2)*((u(2,i,j+1)-u(2,i,j))/deltay)) + &
              1.0D0/2.0D0*(u(1,i,j)*((u(1,i+1,j)-u(1,i,j))/deltax) + &
                           vAvg(1)*((u(1,i,j+1)-u(1,i,j))/deltay))
  END IF

 END DO ! j
END DO ! i


DO i = 2,Nx
 DO j = 3,Ny

    IF (i==Nx) THEN
   
     uAvg(2) = 0.25*(v(2,i,j)+v(2,i-1,j)+v(2,i,j+1)+v(2,i-1,j+1) )
     uAvg(1) = 0.25*(v(1,i,j)+v(1,i-1,j)+v(1,i,j+1)+v(1,i-1,j+1) )


  
     H(2,i,j) = -3.0D0/2.0D0*(uAvg(2)*((v(2,i+1,j)-v(2,i,j))/deltax/2.0D0) + &
                           v(2,i,j)*((v(2,i,j+1)-v(2,i,j))/deltay)) + &
              1.0D0/2.0D0*(uAvg(1)*((v(1,i+1,j)-v(2,i,j))/deltax/2.0D0) + &
                           v(1,i,j)*((v(1,i,j+1)-v(2,i,j))/deltay))
    ELSE
    uAvg(2) = 0.25*(v(2,i,j)+v(2,i-1,j)+v(2,i,j+1)+v(2,i-1,j+1) )
    uAvg(1) = 0.25*(v(1,i,j)+v(1,i-1,j)+v(1,i,j+1)+v(1,i-1,j+1) )


    H(2,i,j) = -3.0D0/2.0D0*(uAvg(2)*((v(2,i+1,j)-v(2,i,j))/deltax) + &
                           v(2,i,j)*((v(2,i,j+1)-v(2,i,j))/deltay)) + &
              1.0D0/2.0D0*(uAvg(1)*((v(1,i+1,j)-v(1,i,j))/deltax) + &
                           v(1,i,j)*((v(1,i,j+1)-v(1,i,j))/deltay))
    END IF

 END DO ! j
END DO ! i


! Error Checking
OPEN( UNIT=1, FILE='Convective_Matrix.dat', STATUS = 'REPLACE')


WRITE(1,'(1X,/, A, /)') 'Convective Matrix H Elements:'
WRITE(1,'(1X,/,A,f8.3,A,f8.3, /)')'Deltax and Deltay are:',deltax,'  ',deltay

DO i =3,Nx
        WRITE(1,'(5000(1X,f10.3, 4X))'),(H(1,i,j),j=2,Ny)
        WRITE(1,*)
ENDDO



WRITE(*,*) 'Finished Calculating the Convective Terms'

END SUBROUTINE CalcConvect
