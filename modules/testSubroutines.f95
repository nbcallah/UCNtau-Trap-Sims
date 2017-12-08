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

SUBROUTINE zOffDipCalc(t, z)
	real(kind=PREC), intent(in) :: t
	real(kind=PREC), intent(out) :: z
	
	real(kind=PREC) :: speed
	real(kind=PREC), dimension(10) :: dipHeights
	real(kind=PREC), dimension(10) :: dipEnds
	
	integer :: i
	
	IF (t > dipEnds(10)) THEN
		zOff = 0.01
	END IF
	
	dipHeights = (/0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010/)
	dipEnds = (/0.0, 20.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 280.0/)
	
	speed = 0.49_8/13.0_8
	
	DO i=1,10,1
		IF (dipEnds(i) > t) THEN
			EXIT
		END IF
	END DO
	
	z = dipHeights(i-1) - speed*(t-dipEnds(i-1))
	
	IF (z < dipHeights(i)) THEN
		z = dipHeights(i)
	END IF
END SUBROUTINE zOffDipCalc

SUBROUTINE trackDaggerHitTime(state)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), dimension(6) :: prevState

	real(kind=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta
	
	integer :: i, numSteps
	
	t = 0.0_8
	
	numSteps = 20e0/dt
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, energy)
		t = t + dt
	END DO
	
	DO
		prevState = state
		CALL symplecticStep(state, dt, energy)
		t = t + dt
		IF (SIGN(1.0_8, state(2)) .NE. SIGN(1.0_8, prevState(2))) THEN
			fracTravel = ABS(prevState(2))/(ABS(state(2)) + ABS(prevState(2)))
			predX = prevState(1) + fracTravel * (state(1) - prevState(1))
			predZ = prevState(3) + fracTravel * (state(3) - prevState(3))
			
			CALL zOffDipCalc(t - 20.0_8, zOff)
			IF (state(1) > 0.0) THEN
				zeta = 0.5_8 - SQRT(state(1)**2 &
				                  + (SQRT(state(2)**2 + (state(3) - zOff)**2) - 1.0_8)**2)
			ELSE
				zeta = 1.0_8 - SQRT(state(1)**2 &
				                  + (SQRT(state(2)**2 + (state(3) - zOff)**2) - 0.5_8)**2)
			END IF
			IF (ABS(predX) < .2 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
				PRINT *, t, predX, predZ - zOff
				EXIT
			END IF
			
			IF (t > 1000) THEN
				EXIT
			END IF
!			PRINT *, predX, predZ
!			PRINT *, state(1), state(2), state(3)
!			PRINT *, prevState(2), prevState(2) + fracTravel * (state(2) - prevState(2)), state(2)
		END IF
	END DO
END SUBROUTINE trackDaggerHitTime

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

SUBROUTINE testEnergyGain(freq, height, sympT, eStart, eEnd)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), intent(in) :: freq, height
	real(kind=PREC), intent(inout) :: sympT
	real(kind=PREC), intent(out) :: eStart, eEnd

	real(kind=PREC), dimension(6) :: state
	integer :: i, numSteps
	
	state = (/0.05_8, 0.0_8, height, 0.0_8, 0.0_8, 0.0_8/)
	CALL calcEnergy(state, eStart)

	numSteps = 10e0/dt
	
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, eEnd, sympT, freq)
		IF (state(3) >= -1.5_8 + ((height + 1.5) / 4.0_8) .AND. state(6) > 0) THEN
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
