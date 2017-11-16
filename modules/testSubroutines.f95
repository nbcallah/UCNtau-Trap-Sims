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

SUBROUTINE trackAndPrint(state, sympT)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	real(kind=PREC), optional, intent(in) :: sympT
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	
	real(kind=PREC) :: pr, rdot, r, pphi, phidot, phi, ptheta, thetadot, theta,&
		ldot, l, phiOffset, prevPhi, totalU, totalKE, t, energy, &
		fx_dan, fy_dan, fz_dan, e_dan, fx_nate, fy_nate, fz_nate, e_nate
	integer :: i, numSteps, numPoints, modulus

!	DO i=1,11,1
!		!PRINT *, "track: ", i
!		CALL randomPoint(state(1,1), state(1,2), state(1,3),&
!		state(1,4), state(1,5), state(1,6), uTrapMax)
!	END DO
		
!	numSteps = 250e0_8/dt
	numSteps = 1000e0/dt
!	modulus = numSteps/50000
!	modulus = numSteps/100000
	t = 0.0_8
	stateRoundErr = 0.0_8
	totalU = 0.0_8
	
!	CALL calcEnergy(state(1,:), energy)
	energy = 0.0
	
	DO i=1,numSteps,1
		IF(INT(dt*1_8*i)-INT(dt*1_8*(i-1)) .NE. 0) THEN
!		IF(1 .EQ. 1) THEN
!			energy = totalU + SUM(state(1,4:6)**2)/(2.0_8*MASS_N)

			PRINT *, dt*i, state(1,1), state(1,2), state(1,3),&
			state(1,4)/MASS_N, state(1,5)/MASS_N, state(1,6)/MASS_N, energy,&
			totalU!, fx, fy, fz

			!PRINT *, dt*i, energy
		END IF
		IF(present(sympT)) THEN
			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, 0.0e0_8, energy, t)
		ELSE
			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, 0.0e0_8, energy)
!			CALL totalForce(state(1,1), state(1,2), state(1,3), fx_nate, fy_nate, fz_nate, e_nate)
!			CALL totalForceDan(state(1,1), state(1,2), state(1,3), fx_dan, fy_dan, fz_dan, e_dan)
!			PRINT *, t, (fx_nate-fx_dan)/fx_nate, (fy_nate-fy_dan)/fy_nate,&
!			(fz_nate-fz_dan)/fz_nate, (e_nate-e_dan)/e_nate
			t = t+dt
		END IF
	END DO
	!PRINT *, t, state(1,1), state(1,2), state(1,3), state(1,4), state(1,5), state(1,6)
END SUBROUTINE trackAndPrint
!
!SUBROUTINE kinetic(px, py, pz, ke)
!	!Calculate kinetic energy, given 3 momenta
!	IMPLICIT NONE
!	real(kind=PREC), intent(in) :: px, py, pz
!	real(kind=PREC), intent(out) :: ke
!	ke = (px*px+py*py+pz*pz)/(2.0e0_8*MASS_N)
!END SUBROUTINE kinetic
!
!SUBROUTINE calcFieldMesh()
!	USE forcesAndPotential
!	USE wires
!	!This subroutine calculates a value over the volume of the trap
!	IMPLICIT NONE
!	real(kind=PREC) :: x, y, z, totalU, totalN, dV, pSphere, eTrap, eOffset
!	integer :: xIt, yIt, zIt
!	
!	!Define trapping energy
!	eTrap = MASS_N * GRAV * 2.0e0_8*torus_major_r*TRAPPING_FRACTION
!	!Calculate number of trapped neutrons
!	pSphere = MASS_N*7.0e0_8
!
!	totalN=0.0_8
!	dV=(2.0e0_8*(torus_major_r+torus_minor_r)/200.0e0_8)&
!		*(2.0e0_8*(torus_major_r+torus_minor_r)/200.0e0_8)&
!		*(2.0e0_8*torus_minor_r/20.0e0_8)
!	
!	!Divide up volume
!	DO xIt=0,200,1
!		DO yIt=0,200,1
!			DO zIt=0,20,1
!				!Set coordinates to Iterator*Span/Iterations-Offset
!				x=(2.0e0_8*(torus_major_r+torus_minor_r)/200.0e0_8)*xIt-(torus_major_r+torus_minor_r)
!				y=(2.0e0_8*(torus_major_r+torus_minor_r)/200.0e0_8)*yIt-(torus_major_r+torus_minor_r)
!				z=(2.0e0_8*torus_minor_r/20.0e0_8)*zIt-torus_minor_r
!				!If it's inside the torus
!				IF (SQRT(z**2+(SQRT(x**2+y**2)-torus_major_r)**2) < torus_minor_r) THEN
!					CALL totalPotential(x, y, z, totalU)
!					totalU = totalU
!					WRITE(*,"(ES13.5, ES13.5, ES13.5, ES13.5)"), x, y, z, totalU
!					IF (totalU < eTrap) THEN
!						!WRITE(*,"(ES13.5, ES13.5, ES13.5, ES13.5)"), x, y, z, totalU
!						!Add a bit of neutrons assuming a constant phase space density
!						!Density based on 7 m/s sphere
!						totalN = totalN + dV*10.0e6_8*(2*MASS_N*(eTrap-totalU))**(3.0e0_8/2.0e0_8)/pSphere**3
!					END IF
!				END IF
!			END DO
!		END DO
!	END DO
!	!WRITE(*,"(ES13.5)"), totalN
!END SUBROUTINE calcFieldMesh
!
!SUBROUTINE lowFieldSearch()
!	USE forcesAndPotential
!	USE wires
!	USE constants
!	IMPLICIT NONE
!	
!	integer :: tIt, pIt, rIt, nEntries, numT, numP, numR
!	real(kind=PREC) :: theta, phi, r, x, y, z, bz, br, bt
!	
!	nEntries = 0
!	
!	numT = 10000
!	numP = 100
!	numR = 100
!	
!	DO tIt = 0,numT,1
!		DO pIt = 0,numP,1
!			DO rIt = 0,numR,1
!				theta = tIt*2.0_8*pi/numT
!				phi = pIt*2.0_8*pi/numP
!				r = rIt*torus_minor_r/numR
!				x = (torus_major_r+r*COS(phi))*COS(theta)
!				y = (torus_major_r+r*COS(phi))*SIN(theta)
!				z = r*SIN(phi)
!				CALL totalB(x, y, z, bz, br, bt)
!				IF(SQRT(bz*bz+br*br+bt*bt) < 1E-2) THEN
!					nEntries = nEntries + 1
!					PRINT *, x, y, z, bz, br, bt
!!					IF(nEntries > 20000) THEN
!!						CALL EXIT(0)
!!					END IF
!!				ELSE IF(MOD(tIt, 100) == 0 .AND. MOD(pIt, 10) == 0 .AND. MOD(rIt, 10) == 0) THEN
!!					PRINT *, x, y, z, bz, br, bt
!				END IF
!			END DO
!		END DO
!	END DO
!END SUBROUTINE lowFieldSearch
!
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
!
!SUBROUTINE testVecExp()
!	USE forcesAndPotential
!	USE wires
!	USE finiteWires
!	!Calculate a field map on one plane for diagnostics
!	IMPLICIT NONE
!	real(kind=PREC) :: x, y, z, x1, y1, z1, x2, y2, z2
!	real(kind=PREC) :: r, cTh, sTh
!	real(kind=PREC), dimension(3,3) :: dv, de
!	integer :: xIt, yIt, zIt
!	
!	x1 = -1.25_8
!	y1 = 0.25_8
!	z1 = 0.25_8
!	x2 = -0.75_8
!	y2 = 0.25_8
!	z2 = 0.25_8
!
!	!DO xIt=0,50,1
!		DO zIt=0,10,1
!			DO yIt = 0,10,1
!				x=-torus_major_r
!				y=(.2_8/10.0_8)*yIt-(0.25_8-(.2_8/10.0_8)/2.0_8)
!				z=(.2_8/10.0_8)*zIt-(0.25_8-(.2_8/10.0_8)/2.0_8)
!				r = SQRT(x*x + y*y)
!				cTh = x/r !Sin(theta)
!				sTh = y/r !Cos(theta)
!				CALL calcFinitePertDerivativesCondensed(r, cTh, sTh, z, dv, x1, y1, z1, x2, y2, z2, 0.5_8)
!				CALL calcFinitePertDerivatives(r, cTh, sTh, z, de, x1, y1, z1, x2, y2, z2, 0.5_8)
!				PRINT *, x, y, z, SUM(dv-de)!, fx/MASS_N, fy/MASS_N, fz/MASS_N
!				!PRINT *, x, y, z, fx+(totalUplus-totalUminus)/(2.0*delta)
!				!WRITE(*,"(I6, I6)"), xIt, zIt
!			END DO
!		END DO
!	!END DO
!END SUBROUTINE testVecExp
!
!SUBROUTINE getRMS()
!	USE forcesAndPotential
!	USE trackGeometry
!	USE symplecticInt
!	USE wires
!	USE constants
!	!Calculate RMS coordinates for a track
!	IMPLICIT NONE
!	real(kind=PREC), dimension(6) :: q, qrms
!	real(kind=PREC) :: u, ke
!	integer :: i, numSteps, n
!	
!	DO n=1,20,1
!		qrms = 0.0e0_8
!		! Get a random point
!		CALL randomPoint(q(1), q(2), q(3), q(4), q(5), q(6), uTrapMax)
!		!Track for 1000 seconds
!		numSteps = 1000.0_8 / dt
!		CALL totalPotential(q(1), q(2), q(3), u)
!		CALL kinetic(q(4), q(5), q(6), ke)
!		!PRINT *, u+ke
!		!WRITE(*,"(ES13.5, ES13.5, ES13.5, ES13.5, ES13.5, ES13.5)"),&
!		!q(1), q(2), q(3), q(4), q(5), q(6)
!		!For each step
!		DO i=0,numSteps,1
!!			IF (MODULO(i,10)==0) THEN
!!				WRITE(*,"(ES13.5, ES13.5, ES13.5, ES13.5, ES13.5, ES13.5)"),&
!!				q(1), q(2), q(3), q(4), q(5), q(6)
!!			END IF
!			!Integrate motion
!			!CALL symplecticStep(q(1), q(2), q(3), q(4), q(5), q(6), dt, 0.0e0_8)
!			!Add square values to average
!			qrms=(q**2+i*qrms)/(i+1)
!		END DO
!		!Sqrt for RMS
!		qrms = SQRT(qrms)
!		WRITE(*,"(ES13.5, ES13.5, ES13.5, ES13.5, ES13.5, ES13.5)"),&
!		qrms(1), qrms(2), qrms(3), qrms(4), qrms(5), qrms(6)
!		CALL totalPotential(q(1), q(2), q(3), u)
!		CALL kinetic(q(4), q(5), q(6), ke)
!		!PRINT *, u+ke
!	END DO
!END SUBROUTINE getRMS
!
!SUBROUTINE trackAndPrint(state, sympT)
!	USE symplecticInt
!	USE constants
!	USE wires
!	USE forcesAndPotential
!	IMPLICIT NONE
!	real(kind=PREC), dimension(3,6), intent(inout) :: state
!	real(kind=PREC), optional, intent(in) :: sympT
!	real(kind=PREC), dimension(3,6) :: stateRoundErr
!	
!	real(kind=PREC) pr, rdot, r, pphi, phidot, phi, ptheta, thetadot, theta,&
!		ldot, l, phiOffset, prevPhi, totalU, totalKE, t, energy,&
!		fx, fy, fz
!	integer :: i, numSteps, numPoints, modulus
!
!!	DO i=1,11,1
!!		!PRINT *, "track: ", i
!!		CALL randomPoint(state(1,1), state(1,2), state(1,3),&
!!		state(1,4), state(1,5), state(1,6), uTrapMax)
!!	END DO
!		
!	numSteps = 250e0_8/dt
!	modulus = numSteps/50000
!	t = 0.0_8
!	stateRoundErr = 0.0_8
!	totalU = 0.0_8
!	
!	CALL calcEnergy(state(1,:), energy)
!	
!	DO i=1,numSteps,1
!		IF(INT(dt*100.0_8*i)-INT(dt*100.0_8*(i-1)) .NE. 0) THEN
!		!IF (MODULO(i, modulus) == 0) THEN
!		!IF(t > 195.0 .AND. t < 205.0) THEN
!			!CALL calcEnergy(state(1,:), energy)
!			!CALL totalPotential(state(1,1), state(1,2), state(1,3), totalU)
!			!CALL totalForce(state(1,1), state(1,2), state(1,3), fx, fy, fz)
!			!CALL totalB(state(1,1), state(1,2), state(1,3), fx, fy, fz)
!			energy = totalU + SUM(state(1,4:6)**2)/(2.0_8*MASS_N)
!			PRINT *, dt*i, state(1,1), state(1,2), state(1,3),&
!			state(1,4)/MASS_N, state(1,5)/MASS_N, state(1,6)/MASS_N, energy,&
!			totalU!, fx, fy, fz
!			!PRINT *, dt*i, energy
!		END IF
!		IF(present(sympT)) THEN
!			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, 0.0e0_8, totalU, t)
!		ELSE
!			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, 0.0e0_8, totalU)
!			t = t+dt
!		END IF
!	END DO
!	!PRINT *, t, state(1,1), state(1,2), state(1,3), state(1,4), state(1,5), state(1,6)
!END SUBROUTINE trackAndPrint
!
!SUBROUTINE testEnergyExtract(state)
!	USE symplecticInt
!	USE constants
!	USE wires
!	USE forcesAndPotential
!	IMPLICIT NONE
!	real(kind=PREC), dimension(3,6), intent(inout) :: state
!	real(kind=PREC) :: totalU, energy, energySUB, fx, fy, fz
!	integer :: i
!	
!	DO i = 1,100,1
!		CALL totalForce(state(1,1), state(1,2), state(1,3), fx, fy, fz, totalU)
!		energy = totalU + SUM(state(1,4:6)*state(1,4:6))/(2.0_8*MASS_N) &
!				+ state(1,1)*GRAV*MASS_N - uZero
!		CALL calcEnergy(state(1,:), energySUB)
!		PRINT *, energySUB - energy
!		state(1,1) = state(1,1) + 0.0001_8
!		state(1,2) = state(1,2) + 0.0001_8
!		state(1,3) = state(1,3) + 0.0001_8
!	END DO
!	
!	PRINT *, energySUB, energy
!	!PRINT *, t, state(1,1), state(1,2), state(1,3), state(1,4), state(1,5), state(1,6)
!END SUBROUTINE testEnergyExtract
!
!SUBROUTINE testKandE()
!	USE ellipticInt
!	IMPLICIT NONE
!	real(kind=PREC) :: k, oldEllK, oldEllE, newEllK, newEllE, dummySum
!	integer(kind=8) :: i
!	
!	dummySum = 0.0_8
!	
!	DO i=0,200000000,1
!		k=.9e0_8*i/200000000.0e0_8
!		CALL ellipticKandE(k, newEllK, newEllE)
!		dummySum = dummySum + newEllK + newEllE
!		!CALL ellipticKandEOld(k, oldEllK, oldEllE)
!		!dummySum = dummySum + oldEllK + oldEllE
!		!PRINT *, k, newEllK, oldEllK, newEllK-oldEllK, newEllE, oldEllE, newEllE-oldEllE
!		!PRINT *, k, newEllE, oldEllE
!	END DO
!	PRINT *, dummySum
!END SUBROUTINE testKandE

END MODULE
