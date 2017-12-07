#include "constants.h"

MODULE testSubroutines

CONTAINS

SUBROUTINE compPots()
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC) :: x, y, z, fx, fy, fz, totalU
	integer :: i
	
	z = -1.49_8
	x = 0.1_8
	y = 0.1_8
	DO i=0,20,1
		CALL totalForce(x, y, z + i*(0.01_8/20.0_8), fx, fy, fz, totalU)
		PRINT *, totalU, fx, fy, fz
	END DO
END SUBROUTINE compPots

SUBROUTINE trackEnergyGain(state, energy_start, energy_end, sympT, freq)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), intent(out) :: energy_start, energy_end
	real(kind=PREC), optional, intent(inout) :: sympT
	real(kind=PREC), optional, intent(in) :: freq
	
	integer :: i, numSteps

	numSteps = 1000e0/dt
!	totalU = 0.0_8
	
	CALL calcEnergy(state, energy_start)
!	energy = 0.0
	
	DO i=1,numSteps,1
		IF(present(sympT)) THEN
			CALL symplecticStep(state, dt, energy_end, sympT, freq)
		ELSE
			CALL symplecticStep(state, dt, energy_end)
			sympT = sympT+dt
		END IF
		IF (100.0_8*energy_end/(MASS_N*GRAV) > 38.0_8 + 5.0_8) THEN
!			PRINT *, "DEAD"
			EXIT
		END IF
	END DO
!	PRINT *, energy_start, energy_end
END SUBROUTINE trackEnergyGain

SUBROUTINE testEnergyGain(freq, height, sympT, eGain)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), intent(in) :: freq, height
	real(kind=PREC), intent(inout) :: sympT
	real(kind=PREC), intent(out) :: eGain

	real(kind=PREC), dimension(6) :: state
	real(kind=PREC) :: eStart, eEnd
	integer :: i, numSteps
	
	state = (/0.05_8, 0.0_8, height, 0.0_8, 0.0_8, 0.0_8/)
	CALL calcEnergy(state, eStart)

	numSteps = 10e0/dt
!	totalU = 0.0_8
	
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, eEnd, sympT, freq)
!		sympT = sympT+dt
!		PRINT *, sympT, eEnd
!		PRINT *, eEnd, state(3), height / 4.0_8, state(6)
		IF (state(3) >= -1.5_8 + ((height + 1.5) / 4.0_8) .AND. state(6) > 0) THEN
			eGain = eEnd - eStart
			EXIT
		END IF
	END DO
!	PRINT *, energy_start, energy_end
END SUBROUTINE testEnergyGain

SUBROUTINE trackAndPrint(state, sympT)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), optional, intent(in) :: sympT
	
	real(kind=PREC) :: pr, rdot, r, pphi, phidot, phi, ptheta, thetadot, theta,&
		ldot, l, phiOffset, prevPhi, totalU, totalKE, t, energy, &
		fx_dan, fy_dan, fz_dan, e_dan, fx_nate, fy_nate, fz_nate, e_nate
	integer :: i, numSteps, numPoints, modulus
		
!	numSteps = 250e0_8/dt
	numSteps = 1000e0/dt
	t = 0.0_8
	totalU = 0.0_8
	
!	CALL calcEnergy(state, energy)
	energy = 0.0
	
	DO i=1,numSteps,1
		IF(INT(dt*10_8*i)-INT(dt*10_8*(i-1)) .NE. 0) THEN
!		IF(1 .EQ. 1) THEN
!			energy = totalU + SUM(state(1,4:6)**2)/(2.0_8*MASS_N)

			PRINT *, dt*i, state(1), state(2), state(3),&
			state(4)/MASS_N, state(5)/MASS_N, state(6)/MASS_N, energy,&
			totalU!, fx, fy, fz

			!PRINT *, dt*i, energy
		END IF
		IF(present(sympT)) THEN
			CALL symplecticStep(state, dt, energy, t, 60.0_8)
		ELSE
			CALL symplecticStep(state, dt, energy)
!			CALL totalForce(state(1), state(2), state(3), fx_nate, fy_nate, fz_nate, e_nate)
!			CALL totalForceDan(state(1), state(2), state(3), fx_dan, fy_dan, fz_dan, e_dan)
!			PRINT *, t, (fx_nate-fx_dan)/fx_nate, (fy_nate-fy_dan)/fy_nate,&
!			(fz_nate-fz_dan)/fz_nate, (e_nate-e_dan)/e_nate
			t = t+dt
		END IF
	END DO
	!PRINT *, t, state(1), state(2), state(3), state(4), state(5), state(6)
END SUBROUTINE trackAndPrint

SUBROUTINE calcx0Mesh()
	USE forcesAndPotential
	USE constants
	!Calculate a field map on one plane for diagnostics
	IMPLICIT NONE
	real(kind=PREC) :: x, y, z, fx, fy, fz, totalUplus
	real(kind=PREC) :: x0, y0, z0
	integer :: xIt, yIt, zIt, aIt
	
	x0 = -2
	y0 = -2
	z0 = -2
	
	DO xIt=0,100,1
		DO yIt=0,100,1
			DO zIt=0,100,1
				x = x0 + 4 * xIt/100.0_8
				y = y0 + 4 * yIt/100.0_8
				z = z0 + 4 * zIt/100.0_8
				CALL totalPotential(x, y, z, totalUplus)
				PRINT *, x, y, z, totalUplus/MASS_N
			END DO
		END DO
	END DO
END SUBROUTINE calcx0Mesh

END MODULE
