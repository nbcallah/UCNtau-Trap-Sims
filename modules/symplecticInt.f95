#include "constants.h"

MODULE symplecticInt

CONTAINS

SUBROUTINE symplecticStep(state, stateRoundErr, deltaT, conservedEnergy, energy, t)
	!Take 1 step of length deltaT using symplectic integration scheme.
	!See Candy & Rozmus 1991 and McLachlan & Atela 1992
	USE constants
	USE forcesAndPotential
	USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
	IMPLICIT NONE
	real(kind=PREC), intent(inout), dimension(6) :: state
	real(kind=PREC), intent(inout), dimension(6) :: stateRoundErr
	real(kind=PREC), intent(in) :: deltaT, conservedEnergy
	real(kind=PREC), intent(out) :: energy
	real(kind=PREC), optional, intent(inout) :: t

	real(kind=PREC) :: fx, fy, fz, totalU
	real(kind=PREC), dimension(6) :: stateA
	integer :: n
	
	!Only gets compiled if VARIABLE isn't defined at the top
	!For each step of the integration,
		!Calc force
	!Unrolling the loop - first iteration gives energy but subsequent won't
	n = 1
	IF(PRESENT(t)) THEN
		CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU, t)
		energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N)! &
!			+ state(1)*GRAV*MASS_N
		t = t + a(n)*deltaT
	ELSE
		CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
		energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N)! &
!			+ state(1)*GRAV*MASS_N
	END IF
	!CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
	!energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N) &
	!	+ state(1)*GRAV*MASS_N - uZero
	!PRINT *, "Force: ", fx, fy, fx
	!Use integration constants (a's and b's)
	!Reference:
	state(4) = state(4) + b(n)*fx*deltaT
	state(5) = state(5) + b(n)*fy*deltaT
	state(6) = state(6) + b(n)*fz*deltaT
	state(1) = state(1) + a(n)*state(4)*deltaT/MASS_N
	state(2) = state(2) + a(n)*state(5)*deltaT/MASS_N
	state(3) = state(3) + a(n)*state(6)*deltaT/MASS_N

	DO n=2,4,1
		!Calc force
		IF(PRESENT(t)) THEN
			CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU, t)
			t = t + a(n)*deltaT
		ELSE
			CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
		END IF
		!CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
		!PRINT *, "Force: ", fx, fy, fx
		!Use integration constants (a's and b's)
		!Reference:
		!px=px+b(n)*fx*deltaT
		!py=py+b(n)*fy*deltaT
		!pz=pz+b(n)*fz*deltaT
		!x=x+a(n)*px*deltaT/MASS_N
		!y=y+a(n)*py*deltaT/MASS_N
		!z=z+a(n)*pz*deltaT/MASS_N
		state(4) = state(4) + b(n)*fx*deltaT
		state(5) = state(5) + b(n)*fy*deltaT
		state(6) = state(6) + b(n)*fz*deltaT
		state(1) = state(1) + a(n)*state(4)*deltaT/MASS_N
		state(2) = state(2) + a(n)*state(5)*deltaT/MASS_N
		state(3) = state(3) + a(n)*state(6)*deltaT/MASS_N
	END DO
END SUBROUTINE symplecticStep

!SUBROUTINE symplecticStep(x, y, z, px, py, pz, deltaT, conservedEnergy)
!	!Take 1 step of length deltaT using symplectic integration scheme.
!	!See Candy & Rozmus 1991 and McLachlan & Atela 1992
!	USE constants
!	USE forcesAndPotential
!	USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
!	IMPLICIT NONE
!	real(kind=PREC), intent(inout) :: x, y, z, px, py, pz
!	real(kind=PREC), intent(in) :: deltaT, conservedEnergy
!
!	real(kind=PREC) :: fx, fy, fz
!	integer :: n
!	
!	!If we want to use variable step size, define more parameters
!#ifdef VARIABLE
!	real(kind=PREC) :: currentDt, x0, y0, z0, px0, py0, pz0,&
!		prevX, prevY, prevZ, prevPx, prevPy, prevPz, distance, prevE, curE
!	integer :: numSteps, numStepsIt
!#endif
!	
!	!Only gets compiled if VARIABLE isn't defined at the top
!#ifndef VARIABLE
!	!For each step of the integration,
!	DO n=1,4,1
!		!Calc force
!		CALL totalForce(x, y, z, fx, fy, fz)
!		!Use integration constants (a's and b's)
!		px=px+b(n)*fx*deltaT
!		py=py+b(n)*fy*deltaT
!		pz=pz+b(n)*fz*deltaT
!		x=x+a(n)*px*deltaT/MASS_N
!		y=y+a(n)*py*deltaT/MASS_N
!		z=z+a(n)*pz*deltaT/MASS_N
!	END DO
!#endif
!
!	!If variable is defined at the top:
!#ifdef VARIABLE
!	!If we want variable step sizes, need to start at the beginning each time.
!	x0=x
!	y0=y
!	z0=z
!	px0=px
!	py0=py
!	pz0=pz
!	!double eStart
!	!double eEnd
!
!	!Start with nominal step size
!	numSteps = 1
!	currentDt = deltaT
!	
!	DO
!		! get energy to compare to next steps
!		x=x0
!		y=y0
!		z=z0
!		px=px0
!		py=py0
!		pz=pz0
!		!For each sub-step, integrate motion
!		DO numStepsIt=1,numSteps,1
!			DO n=1,4,1
!				CALL totalForce(x, y, z, fx, fy, fz)
!				px=px+b(n)*fx*currentDt
!				py=py+b(n)*fy*currentDt
!				pz=pz+b(n)*fz*currentDt
!				x=x+a(n)*px*currentDt/MASS_N
!				y=y+a(n)*py*currentDt/MASS_N
!				z=z+a(n)*pz*currentDt/MASS_N
!			END DO
!		END DO
!		!Calc energy to compare
!		prevE = curE
!		CALL calcEnergy((/x,y,z,px,py,pz/), curE)
!		IF (numSteps > 1) THEN
!			!distance = ((prevX-x)**2+(prevY-y)**2+(prevZ-z)**2)/0.1468995e0_8&
!					!+((prevPx-px)**2+(prevPy-py)**2+(prevPz-pz)**2)/(MASS_N*MASS_N*2.68842e0_8)
!			!If it's converged, exit
!			!IF(distance < 1.0d-8) THEN
!			!	EXIT
!			!END IF
!			!IF(ABS((curE-prevE)/(prevE)) < 1.0d-9) THEN
!			IF(numSteps == 4) THEN
!				!WRITE(ERROR_UNIT,*) numSteps
!				EXIT
!			END IF
!		END IF
!		IF (numSteps > 64) THEN
!			WRITE(ERROR_UNIT,*) "Integrator did not converge after 64 steps! Terminating"
!			EXIT
!		END IF
!		prevX=x
!		prevY=y
!		prevZ=z
!		prevPx=px
!		prevPy=py
!		prevPz=pz
!		!If not, try again with a smaller step size
!		currentDt = currentDt / 2.0e0_8
!		numSteps = numSteps * 2
!	END DO
!	
!#endif
!END SUBROUTINE symplecticStep

END MODULE
