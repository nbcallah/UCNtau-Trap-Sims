#include "constants.h"

MODULE trackGeometry

CONTAINS

SUBROUTINE randomPointSpring(x, y, z, px, py, pz, e)
	IMPLICIT NONE
	real(kind=PREC), intent(out) :: x, y, z, px, py, pz
	real(kind=PREC), intent(in) :: e
	real(kind=PREC) :: r

	!x = RAND(0)
	!y = RAND(0)
	CALL RANDOM_NUMBER(x)
	CALL RANDOM_NUMBER(y)
	z = 0.0e0_8
	px=MASS_N*x
	py=MASS_N*y
	pz=0.0e0_8
END SUBROUTINE randomPointSpring
	

SUBROUTINE randomPoint(x, y, z, px, py, pz, e)
	! Generate a random point inside the trap with an energy e
	USE constants
	USE forcesAndPotential
	USE wires
	IMPLICIT NONE
	real(kind=PREC), intent(out) :: x, y, z, px, py, pz
	real(kind=PREC), intent(in) :: e
	real(kind=PREC) :: l, totalU, pLen, pMag, pTrap, rej, u, theta
	
	DO
		DO
			DO
				CALL RANDOM_NUMBER(x)
				CALL RANDOM_NUMBER(y)
				CALL RANDOM_NUMBER(z)
!				x = (RAND(0)-0.5e0_8)*2.0e0_8*(torus_major_r+torus_minor_r)
!				y = (RAND(0)-0.5e0_8)*2.0e0_8*(torus_major_r+torus_minor_r)
!				z = (RAND(0)-0.5e0_8)*2.0e0_8*(torus_minor_r)
				x = (x-0.5e0_8)*2.0e0_8*(torus_major_r+torus_minor_r)
				y = (y-0.5e0_8)*2.0e0_8*(torus_major_r+torus_minor_r)
				z = (z-0.5e0_8)*2.0e0_8*(torus_minor_r)
				!l is distance from center-line of torus to point
				l = SQRT(z**2+(torus_major_r-SQRT(x**2+y**2))**2)
				IF (l < torus_minor_r) THEN
					EXIT
				END IF
			END DO
			CALL totalPotential(x, y, z, totalU)
			IF (totalU < e) THEN
				EXIT
			END IF
		END DO

		pMag = SQRT(2.0e0_8*MASS_N*(e-totalU))
		pTrap = SQRT(2.0e0_8*MASS_N*uTrapMax)

		CALL RANDOM_NUMBER(rej)
		!IF (RAND(0) < (pMag**2)/(pTrap**2)) THEN
		IF (rej < (pMag**2)/(pTrap**2)) THEN
!			px = RAND(0)
!			py = RAND(0)
!			pz = RAND(0)
			!CALL RANDOM_NUMBER(px)
			!CALL RANDOM_NUMBER(py)
			!CALL RANDOM_NUMBER(pz)
			CALL RANDOM_NUMBER(u)
			u = 2.0_8*u - 1.0_8
			CALL RANDOM_NUMBER(theta)
			theta = 2.0_8*PI*theta
			px = SQRT(1.0_8-u**2)*COS(theta)
			py = SQRT(1.0_8-u**2)*SIN(theta)
			pz = u
!			pLen = 1.0_8/SQRT(px**2+py**2+pz**2)
!			px = px * pLen
!			py = py * pLen
!			pz = pz * pLen

			!Scale by ~e-U to get a total energy of e
			px = SQRT(2.0e0_8*MASS_N*(e-totalU)) * px
			py = SQRT(2.0e0_8*MASS_N*(e-totalU)) * py
			pz = SQRT(2.0e0_8*MASS_N*(e-totalU)) * pz
			
			EXIT
		END IF
	END DO
END SUBROUTINE randomPoint

SUBROUTINE randomPointBox(x, y, z, px, py, pz, e)
	! Generate a random point inside the trap with an energy < e
	USE constants
	USE forcesAndPotential
	USE wires
	IMPLICIT NONE
	real(kind=PREC), intent(out) :: x, y, z, px, py, pz
	real(kind=PREC), intent(in) :: e
	real(kind=PREC) :: totalE, minU
	
	CALL calcEnergy((/-1.5_8*torus_major_r, 1.5_8*torus_major_r, 1.5_8 * torus_minor_r, 0.0_8, 0.0_8, 0.0_8/), minU, 0.0_8)
	
	!Make uniform distribution in phase space
	DO
		!Create a point inside a cube
		CALL RANDOM_NUMBER(x)
		CALL RANDOM_NUMBER(y)
		CALL RANDOM_NUMBER(z)
		x = x * 2.0_8 * 1.5_8 * torus_major_r - 1.5_8 * torus_major_r
		y = y * 2.0_8 * 1.5_8 * torus_major_r - 1.5_8 * torus_major_r
		z = z * 2.0_8 * 1.5_8 * torus_major_r - 1.5_8 * torus_major_r
		IF (z > 1.5_8 * torus_minor_r) THEN !Fairly distributes points in a cube, but rejects points outside of my trap's bounding box
			CYCLE
		END IF
		!Create a point inside a momentum-cube
		CALL RANDOM_NUMBER(px)
		CALL RANDOM_NUMBER(py)
		CALL RANDOM_NUMBER(pz)
		!Scale so that it has at most (e-minU) energy.
		px = px * 2.0_8 * SQRT(2.0_8*MASS_N*(e-minU)) - SQRT(2.0_8*MASS_N*(e-minU))
		py = py * 2.0_8 * SQRT(2.0_8*MASS_N*(e-minU)) - SQRT(2.0_8*MASS_N*(e-minU))
		pz = pz * 2.0_8 * SQRT(2.0_8*MASS_N*(e-minU)) - SQRT(2.0_8*MASS_N*(e-minU))
		CALL calcEnergy((/x, y, z, px, py, pz/), totalE)!, 0.0_8)
		IF (totalE < e) THEN
			EXIT
		END IF
	END DO
END SUBROUTINE randomPointBox

SUBROUTINE resetState(state)
	!Reset the 2 tracks back along their separation vectors.
	IMPLICIT NONE
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	real(kind=PREC), dimension(2,6) :: delta
	real(kind=PREC), dimension(1,6) :: rejection
	real(kind=PREC) :: lenA, lenB, lenrejAB
	
	!Calculate difference vectors
	delta(1,:) = (state(2,:)-state(1,:))&
	/METRIC_SCALE
	delta(2,:) = (state(3,:)-state(1,:))&
	/METRIC_SCALE
	
	!Calculate the rejection, or perp. part of A from B
	rejection(1,:) = delta(2,:)&
					-delta(1,:)&
						*(SUM(delta(1,:)*delta(2,:))/SUM(delta(1,:)**2))
	
	! Length is just the sum of the square of the difference vector components
	lenA = SQRT(SUM(delta(1,:)**2))
	lenrejAB = SQRT(SUM(rejection(1,:)**2))
	
	!Normalize and re-scale the difference to the initial delta
	delta(1,:) = delta(1,:) / lenA
	rejection(1,:) = rejection(1,:) / lenrejAB
	delta(1,:) = (1.0e-9_8&
		*METRIC_SCALE)&
		* delta(1,:)
	rejection(1,:) = (1.0e-9_8&
		*METRIC_SCALE)&
		* rejection(1,:)
	
	!Set the secondary trajectories to be ref+delta
	state(2,:) = state(1,:) + delta(1,:)
	state(3,:) = state(1,:) + rejection(1,:)
	
!	PRINT *, "Reset Separations:"
!	PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/METRIC_SCALE)**2))
!	PRINT *, SQRT(SUM(((state(1,:)-state(3,:))/METRIC_SCALE)**2))
END SUBROUTINE resetState

SUBROUTINE cross(a, b, res)
	IMPLICIT NONE
	real(kind=PREC), dimension(3), intent(in) :: a, b
	real(kind=PREC), dimension(3), intent(out) :: res
	res(1) = a(2)*b(3)-a(3)*b(2)
	res(2) = a(1)*b(3)-a(3)*b(1)
	res(3) = a(1)*b(2)-a(2)*b(1)
END SUBROUTINE cross

SUBROUTINE createInitialTrajectories(state)
	IMPLICIT NONE
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	real(kind=PREC), dimension(3) :: randomP, perp1, perp2, a, b
	real(kind=PREC) :: pNorm, pLen
	
	! Generate a random vector
!	randomP(1) = RAND(0)
!	randomP(2) = RAND(0)
!	randomP(3) = RAND(0)
	CALL RANDOM_NUMBER(randomP)
	
	!Normalize and scale it so it's similar to the reference trajectory
	pLen = SQRT(SUM(state(1,4:6)**2))
	pNorm = SQRT(SUM(randomP**2))

	randomP = randomP/pNorm
	randomP = pLen*randomP
	
	!Cross Product, generates 2 perpendicular vectors
	CALL cross(state(1,4:6), randomP, perp1)
	CALL cross(state(1,4:6), perp1, perp2)
	
	!Now that we have 2 perpindicular vectors, start with the initial p vector
	!For both trajectories. Then subtract delta from the vector and add
	!Delta in the perpindicular direction.
	
	!Now we have 2 perpindicular vectors to the trajectory in P space.
	!Want to subtract a from P and add b in perp direction
	!Constrains distance to be epsilon, size to still be P
	!N.B. Only approximate separation in metric space.
	a=EPSILON*EPSILON*4.0e-27_8*4.0e-27_8/(2.0e0_8*pNorm)
	b=EPSILON*4.0e-27_8*SQRT(1.0e0_8-EPSILON*EPSILON*4.0e-27_8*4.0e-27_8/(4.0e0_8*pNorm*pNorm))
	
	!Normalize and scale by b
	pNorm = SQRT(SUM(perp1**2))
	perp1 = perp1 / pNorm
	perp1 = perp1 * b
	
	pNorm = SQRT(SUM(perp2**2))
	perp2 = perp2 / pNorm
	perp2 = perp2 * b
	
	!Subtract off a and add in b in perp direction
	state(2,4:6) = (1.0e0_8-a/pLen)*state(1,4:6)+perp1
	state(3,4:6) = (1.0e0_8-a/pLen)*state(1,4:6)+perp2
	
	!keep real space coord.
	state(2,1:3) = state(1,1:3)
	state(3,1:3) = state(1,1:3)
END SUBROUTINE createInitialTrajectories

!SUBROUTINE createInitialTrajectories(state)
!	IMPLICIT NONE
!	real(kind=PREC), dimension(3,6), intent(inout) :: state
!	real(kind=PREC), dimension(3) :: randomP, perp1, perp2, a, b
!	real(kind=PREC) :: pNorm, pLen
!	
!	! Generate a random vector
!!	randomP(1) = RAND(0)
!!	randomP(2) = RAND(0)
!!	randomP(3) = RAND(0)
!	CALL RANDOM_NUMBER(randomP)
!	
!	!Normalize and scale it so it's similar to the reference trajectory
!	pLen = SQRT(SUM((state(1,4:6)**2)**2))
!	pNorm = SQRT(SUM(randomP**2))
!
!	randomP = randomP/pNorm
!	randomP = pLen*randomP
!	
!	!Cross Product, generates 2 perpendicular vectors
!	CALL cross(state(1,4:6), randomP, perp1)
!	CALL cross(state(1,4:6), perp1, perp2)
!	
!	!Now we have 2 perpindicular vectors to the trajectory in P space.
!	!Want to subtract a from P and add b in perp direction
!	!Constrains distance to be epsilon, size to still be P
!	a=EPSILON*EPSILON*4.0d-27*4.0d-27/(2.0e0_8*pNorm)
!	b=EPSILON*4.0d-27*SQRT(1.0e0_8-EPSILON*EPSILON*4.0d-27*4.0d-27/(4.0e0_8*pNorm*pNorm))
!	
!	!Normalize and scale by b
!	pNorm = SQRT(SUM(perp1**2))
!	perp1 = perp1 / pNorm
!	perp1 = perp1 * b
!	
!	pNorm = SQRT(SUM(perp2**2))
!	perp2 = perp2 / pNorm
!	perp2 = perp2 * b
!	
!	!Subtract off a and add in b in perp direction
!	state(2,4:6) = (1.0e0_8-a/pLen)*state(1,4:6)+perp1
!	state(3,4:6) = (1.0e0_8-a/pLen)*state(1,4:6)+perp2
!	
!	!keep real space coord.
!	state(2,1:3) = state(1,1:3)
!	state(3,1:3) = state(1,1:3)
!END SUBROUTINE createInitialTrajectories

END MODULE
