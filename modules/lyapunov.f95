#include "constants.h"

MODULE lyapunov

CONTAINS

SUBROUTINE calcLyapunov(state, results, t)
	USE trackGeometry
	USE symplecticInt
	USE wires
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	!state will be an array;
	!first index: 1-ref, 2-traj_a, 3-traj_b
	!second index: 1-x, 2-y, 3-z, 4-px, 5-py, 6-pz
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	!results: lyap1, lyap2, rmsx, y, z, px, py, pz, eStart, eEnd, maxE
	real(kind=PREC), dimension(11), intent(out) :: results
	real(kind=PREC), optional, intent(inout) :: t
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	real(kind=PREC), dimension(2,6) :: delta
	real(kind=PREC), dimension(1,6) :: rejection
	real(kind=PREC), dimension(2,nsep) :: separations
	real(kind=PREC), dimension(1,6) :: qrms
	
	real(kind=PREC) eStart, eMax, eCur
	real(kind=PREC) dumpVar, tStart
	integer :: i, j, numSteps, rmsNum
	
	!CALL randomPoint(state(1,1), state(1,2), state(1,3),&
		!state(1,4), state(1,5), state(1,6), uTrapMax)
	!CALL createInitialTrajectories(state)
	
	qrms = 0.0e0_8
	stateRoundErr = 0.0_8
	eMax = 0.0_8
	
!	PRINT *, "Starting Separations:"
!	PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/METRIC_SCALE)**2))
!	PRINT *, SQRT(SUM(((state(1,:)-state(3,:))/METRIC_SCALE)**2))
!	
	
	CALL calcEnergy(state(1,:), results(9))
	numSteps = liptime/dt
	rmsNum = 0
	DO i=1,nsep,1
!		PRINT *, "Starting Separations:"
!		PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/&
!			METRIC_SCALE)**2))&
!			,SQRT(SUM(((state(1,:)-state(3,:))/&
!			METRIC_SCALE)**2))
		DO j=1,numSteps,1
			IF(PRESENT(t)) THEN
				tStart = t
				CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, eStart, eCur, t)
				t = tStart
				CALL symplecticStep(state(2,:), stateRoundErr(2,:), dt, eStart, dumpVar, t)
				t = tStart
				CALL symplecticStep(state(3,:), stateRoundErr(3,:), dt, eStart, dumpVar, t)
			ELSE
				CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, eStart, eCur)
				CALL symplecticStep(state(2,:), stateRoundErr(2,:), dt, eStart, dumpVar)
				CALL symplecticStep(state(3,:), stateRoundErr(3,:), dt, eStart, dumpVar)
			END IF
			!CALL symplecticStep(state(3,1), state(3,2), state(3,3),&
				!state(3,4), state(3,5), state(3,6), dt, eStart)
			qrms(1,:)=(state(1,:)**2+rmsNum*qrms(1,:))/(rmsNum+1)
			IF (eCur > eMax) THEN
				eMax = eCur
			END IF
			rmsNum = rmsNum+1
		END DO
	
		delta(1,:) = (state(2,:)-state(1,:))&
		/METRIC_SCALE
		delta(2,:) = (state(3,:)-state(1,:))&
		/METRIC_SCALE

		separations(1,i) = SQRT(SUM(delta(1,:)**2))
		rejection(1,:) = delta(2,:)&
					-delta(1,:)&
						*(SUM(delta(1,:)*delta(2,:))/SUM(delta(1,:)**2))
		separations(2,i) = SQRT(SUM(rejection(1,:)**2))
		
		!WRITE(*,"(ES13.5, ES13.5)"), separations(1,i), separations(2,i)
		CALL resetState(state)
!		stateRoundErr(2,1) = 0.0_8
!		stateRoundErr(2,2) = 0.0_8
!		stateRoundErr(2,3) = 0.0_8
!		stateRoundErr(3,1) = 0.0_8
!		stateRoundErr(3,2) = 0.0_8
!		stateRoundErr(3,3) = 0.0_8
	END DO
	
	!CALL calcEnergy(state(1,:), results(10))
	results(10) = eCur
	results(11) = eMax
	
	qrms = SQRT(qrms)
	
	results(3:8) = qrms(1,:)
	results(1) = (1.0e0_8/(nsep*liptime))*SUM(LOG(separations(1,:)/EPSILON))
	results(2) = (1.0e0_8/(nsep*liptime))*SUM(LOG(separations(2,:)/EPSILON))
	
	!WRITE(*,&
	!	"(ES27.19E3,ES27.19E3,ES27.19E3,ES27.19E3,ES27.19E3,ES27.19E3,ES27.19E3,ES27.19E3)",&
	!	ADVANCE="NO"), &
	!	(1.0e0_8/(nsep*liptime))*SUM(LOG(separations(1,:)/EPSILON)), &
	!	(1.0e0_8/(nsep*liptime))*SUM(LOG(separations(2,:)/EPSILON)), &
	!	qrms(1,:)
	
	!PRINT *, (1.0e0_8/(nsep*liptime))*SUM(LOG(separations(1,:)/EPSILON)), &
	!	(1.0e0_8/(nsep*liptime))*SUM(LOG(separations(2,:)/EPSILON)), &
	!	qrms(1,:)

	!DO i=1,nsep,1
	!	WRITE(*,"(ES13.5, ES13.5)"), separations(1,i), separations(2,i)
	!END DO
END SUBROUTINE calcLyapunov

SUBROUTINE rampTrack(state, results, t)
	USE trackGeometry
	USE symplecticInt
	USE wires
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	!state will be an array;
	!first index: 1-ref, 2-traj_a, 3-traj_b
	!second index: 1-x, 2-y, 3-z, 4-px, 5-py, 6-pz
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	!results: 0, 0, startx, y, z, endx, y, z, eStart, eEnd, maxE
	real(kind=PREC), dimension(11), intent(out) :: results
	real(kind=PREC), optional, intent(inout) :: t
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	real(kind=PREC) :: eCur, eStart
	integer :: i, j, numSteps
	
	CALL calcEnergy(state(1,:), results(9), 0.0_8)
	results(1) = 0.0_8
	results(2) = 0.0_8
	results(3) = state(1,1)
	results(4) = state(1,2)
	results(5) = state(1,3)
	
	numSteps = 90.0_8/dt
	!DO i=1,nsep,1
	DO j=1,numSteps,1
		IF(PRESENT(t)) THEN
			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, eStart, eCur, t)
		ELSE
			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, eStart, eCur)
		END IF
!			IF(MOD(j,numSteps/50) == 0) THEN
!				PRINT *, t, state(1,1), state(1,2), state(1,3),&
!					state(1,4)/MASS_N, state(1,5)/MASS_N, state(1,6)/MASS_N, eCur,&
!					eCur
!			END IF
		CALL reflect(state(1,:))
	END DO
	!END DO
	results(11) = 0.0_8
	results(10) = eCur
	results(6) = state(1,1)
	results(7) = state(1,2)
	results(8) = state(1,3)
END SUBROUTINE rampTrack

SUBROUTINE reflect(state)
	USE wires
	USE constants
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC) :: u1, u2, theta, phi
	real(kind=PREC), dimension(3,3) :: rotM1, rotM2
	real(kind=PREC), dimension(3) :: norm
	integer :: refl
	
	refl = 0
!	PRINT *, state(1), state(2), state(3)
	IF(state(1) < 1.5_8 * (-torus_major_r) .and. state(4) < 0.0_8) THEN
		refl=1
		norm = (/1.0_8, 0.0_8, 0.0_8/)
		CALL RANDOM_NUMBER(u1)
		CALL RANDOM_NUMBER(u2)
		theta = ASIN(SQRT(u1))
		phi = 2.0_8 * PI * u2
		rotM1 = RESHAPE((/COS(theta), 0.0_8, SIN(theta),&
						0.0_8, 1.0_8, 0.0_8,&
						-SIN(theta), 0.0_8, COS(theta)/), SHAPE(rotM1))
		rotM2 = RESHAPE((/1.0_8, 0.0_8, 0.0_8,&
						0.0_8, COS(phi), -SIN(phi),&
						0.0_8, SIN(phi), COS(phi)/), SHAPE(rotM2))
	ELSE IF(state(1) > 1.5_8 * (torus_major_r) .and. state(4) > 0.0_8) THEN
		refl=1
		norm = (/-1.0_8, 0.0_8, 0.0_8/)
		CALL RANDOM_NUMBER(u1)
		CALL RANDOM_NUMBER(u2)
		theta = ASIN(SQRT(u1))
		phi = 2.0_8 * PI * u2
		rotM1 = RESHAPE((/COS(theta), 0.0_8, SIN(theta),&
						0.0_8, 1.0_8, 0.0_8,&
						-SIN(theta), 0.0_8, COS(theta)/), SHAPE(rotM1))
		rotM2 = RESHAPE((/1.0_8, 0.0_8, 0.0_8,&
						0.0_8, COS(phi), -SIN(phi),&
						0.0_8, SIN(phi), COS(phi)/), SHAPE(rotM2))
	ELSE IF(state(2) < 1.5_8 * (-torus_major_r) .and. state(5) < 0.0_8) THEN
		refl=1
		norm = (/0.0_8, 1.0_8, 0.0_8/)
		CALL RANDOM_NUMBER(u1)
		CALL RANDOM_NUMBER(u2)
		theta = ASIN(SQRT(u1))
		phi = 2.0_8 * PI * u2
		rotM1 = RESHAPE((/1.0_8, 0.0_8, 0.0_8,&
						0.0_8, COS(theta), -SIN(theta),&
						0.0_8, SIN(theta), COS(theta)/), SHAPE(rotM1))
		rotM2 = RESHAPE((/COS(phi), 0.0_8, SIN(phi),&
						0.0_8, 1.0_8, 0.0_8,&
						-SIN(phi), 0.0_8, COS(phi)/), SHAPE(rotM2))
	ELSE IF(state(2) > 1.5_8 * (torus_major_r) .and. state(5) > 0.0_8) THEN
		refl=1
		norm = (/0.0_8, -1.0_8, 0.0_8/)
		CALL RANDOM_NUMBER(u1)
		CALL RANDOM_NUMBER(u2)
		theta = ASIN(SQRT(u1))
		phi = 2.0_8 * PI * u2
		rotM1 = RESHAPE((/1.0_8, 0.0_8, 0.0_8,&
						0.0_8, COS(theta), -SIN(theta),&
						0.0_8, SIN(theta), COS(theta)/), SHAPE(rotM1))
		rotM2 = RESHAPE((/COS(phi), 0.0_8, SIN(phi),&
						0.0_8, 1.0_8, 0.0_8,&
						-SIN(phi), 0.0_8, COS(phi)/), SHAPE(rotM2))
	ELSE IF(state(3) < 1.5_8 * (-torus_minor_r) .and. state(6) < 0.0_8) THEN
		refl=1
		norm = (/0.0_8, 0.0_8, 1.0_8/)
		CALL RANDOM_NUMBER(u1)
		CALL RANDOM_NUMBER(u2)
		theta = ASIN(SQRT(u1))
		phi = 2.0_8 * PI * u2
		rotM1 = RESHAPE((/1.0_8, 0.0_8, 0.0_8,&
						0.0_8, COS(theta), -SIN(theta),&
						0.0_8, SIN(theta), COS(theta)/), SHAPE(rotM1))
		rotM2 = RESHAPE((/COS(phi), SIN(phi), 0.0_8,&
						-SIN(phi), COS(phi), 0.0_8,&
						0.0_8, 0.0_8, 1.0_8/), SHAPE(rotM2))
	ELSE IF(state(3) > 1.5_8 * (torus_minor_r) .and. state(6) > 0.0_8) THEN
		refl=1
		norm = (/0.0_8, 0.0_8, -1.0_8/)
		CALL RANDOM_NUMBER(u1)
		CALL RANDOM_NUMBER(u2)
		theta = ASIN(SQRT(u1))
		phi = 2.0_8 * PI * u2
		rotM1 = RESHAPE((/1.0_8, 0.0_8, 0.0_8,&
						0.0_8, COS(theta), -SIN(theta),&
						0.0_8, SIN(theta), COS(theta)/), SHAPE(rotM1))
		rotM2 = RESHAPE((/COS(phi), SIN(phi), 0.0_8,&
						-SIN(phi), COS(phi), 0.0_8,&
						0.0_8, 0.0_8, 1.0_8/), SHAPE(rotM2))
	END IF
	
	IF(refl == 1) THEN
!		PRINT *, "--------"
!		PRINT *, theta, phi
!		PRINT *, norm
		norm = MATMUL(rotM1, norm)
		norm = MATMUL(rotM2, norm)
!		PRINT *, norm
!		PRINT *, "--------"
!		PRINT *, SUM(norm**2)
		norm = norm * SQRT(state(4)**2 + state(5)**2 + state(6)**2)
		state(4) = norm(1)
		state(5) = norm(2)
		state(6) = norm(3)
	END IF
END SUBROUTINE reflect

SUBROUTINE escapeTime(state, results)
	USE trackGeometry
	USE symplecticInt
	USE wires
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	!state will be an array;
	!first index: 1-ref, 2-traj_a, 3-traj_b
	!second index: 1-x, 2-y, 3-z, 4-px, 5-py, 6-pz
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	!results: t, radius, rmsx, y, z, px, py, pz, eStart, eEnd, maxE
	real(kind=PREC), dimension(11), intent(out) :: results
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	real(kind=PREC), dimension(1,6) :: qrms
	
	real(kind=PREC) eStart, eMax, eCur
	real(kind=PREC) dumpVar
	integer :: i, j, numSteps, rmsNum
	
	qrms = 0.0e0_8
	stateRoundErr = 0.0_8
	eMax = 0.0_8
	
	CALL calcEnergy(state(1,:), results(9))
	numSteps = liptime/dt
	rmsNum = 0
	i = 0
	DO
		IF(state(1,1) > cleanerHeight .OR. state(1,1) < -1.5) THEN
			EXIT
		END IF
		IF(dt * i > 3000) THEN
		!IF(dt * i > 500) THEN
			!i = -1
			EXIT
		END IF
		CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, eStart, eCur)
		qrms(1,:)=(state(1,:)**2+rmsNum*qrms(1,:))/(rmsNum+1)
		IF (eCur > eMax) THEN
			eMax = eCur
		END IF
		rmsNum = rmsNum+1
		i = i+1
	END DO
	
	CALL calcEnergy(state(1,:), results(10))
	results(11) = eMax
	
	qrms = SQRT(qrms)
	
	results(3:8) = qrms(1,:)
	results(1) = dt * i
	!results(2) = SQRT(state(1,2)**2 + state(1,3)**2)
	results(2) = state(1,1)
END SUBROUTINE escapeTime

SUBROUTINE trackSumLyapunov(state)
	USE trackGeometry
	USE symplecticInt
	USE wires
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	!state will be an array;
	!first index: 1-ref, 2-traj_a, 3-traj_b
	!second index: 1-x, 2-y, 3-z, 4-px, 5-py, 6-pz
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	real(kind=PREC), dimension(2,6) :: delta
	real(kind=PREC), dimension(1,6) :: rejection
	real(kind=PREC), dimension(2,nsep) :: separations
	
	real(kind=PREC) eStart, U
	integer :: i, j, numSteps
	
	!CALL randomPoint(state(1,1), state(1,2), state(1,3),&
		!state(1,4), state(1,5), state(1,6), uTrapMax)
	!CALL createInitialTrajectories(state)
	
	separations = 0.0_8
	
	stateRoundErr = 0.0_8
	
!	PRINT *, "Starting Separations:"
!	PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/METRIC_SCALE)**2))
!	PRINT *, SQRT(SUM(((state(1,:)-state(3,:))/METRIC_SCALE)**2))
!	
	
	numSteps = liptime/dt
	DO i=1,nsep,1
!		PRINT *, "Starting Separations:"
!		PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/&
!			METRIC_SCALE)**2))&
!			,SQRT(SUM(((state(1,:)-state(3,:))/&
!			METRIC_SCALE)**2))
		DO j=1,numSteps,1
			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, 0.0_8, U)
			CALL symplecticStep(state(2,:), stateRoundErr(2,:), dt, 0.0_8, U)
			CALL symplecticStep(state(3,:), stateRoundErr(3,:), dt, 0.0_8, U)
		END DO
	
		delta(1,:) = (state(2,:)-state(1,:))&
		/METRIC_SCALE
		delta(2,:) = (state(3,:)-state(1,:))&
		/METRIC_SCALE

		separations(1,i) = SQRT(SUM(delta(1,:)**2))
		rejection(1,:) = delta(2,:)&
					-delta(1,:)&
						*(SUM(delta(1,:)*delta(2,:))/SUM(delta(1,:)**2))
		separations(2,i) = SQRT(SUM(rejection(1,:)**2))
		
		PRINT *, (1.0e0_8/(liptime))*SUM(LOG(separations(1,1:i)/EPSILON))
		
		!WRITE(*,"(ES13.5, ES13.5)"), separations(1,i), separations(2,i)
		CALL resetState(state)
		stateRoundErr(2,1) = 0.0_8
		stateRoundErr(2,2) = 0.0_8
		stateRoundErr(2,3) = 0.0_8
		stateRoundErr(3,1) = 0.0_8
		stateRoundErr(3,2) = 0.0_8
		stateRoundErr(3,3) = 0.0_8
	END DO
END SUBROUTINE trackSumLyapunov

SUBROUTINE trackLyapunov(state)
	USE trackGeometry
	USE symplecticInt
	USE wires
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	!state will be an array;
	!first index: 1-ref, 2-traj_a, 3-traj_b
	!second index: 1-x, 2-y, 3-z, 4-px, 5-py, 6-pz
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	real(kind=PREC), dimension(2,6) :: delta
	real(kind=PREC), dimension(1,6) :: rejection
	real(kind=PREC), dimension(2) :: separations
	real(kind=PREC) :: res1, res2, e, U
	
	integer :: i, j, numSteps
	
	stateRoundErr = 0.0_8
	
	numSteps = liptime/dt

	DO i=1,nsep,1
!		PRINT *, "Starting Separations:"
!		PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/&
!			METRIC_SCALE)**2))&
!			,SQRT(SUM(((state(1,:)-state(3,:))/&
!			METRIC_SCALE)**2))
		DO j=1,numSteps,1
			CALL symplecticStep(state(1,:), stateRoundErr(1,:), dt, 0.0_8, U)
			CALL symplecticStep(state(2,:), stateRoundErr(2,:), dt, 0.0_8, U)
			CALL symplecticStep(state(3,:), stateRoundErr(3,:), dt, 0.0_8, U)
			CALL calcEnergy(state(1,:), e)
			
			delta(1,:) = (state(2,:)-state(1,:))&
			/METRIC_SCALE
			delta(2,:) = (state(3,:)-state(1,:))&
			/METRIC_SCALE

			separations(1) = SQRT(SUM(delta(1,:)**2))
			rejection(1,:) = delta(2,:)&
						-delta(1,:)&
							*(SUM(delta(1,:)*delta(2,:))/SUM(delta(1,:)**2))
			separations(2) = SQRT(SUM(rejection(1,:)**2))
			
			res1 = (1.0e0_8/(liptime))*LOG(separations(1)/EPSILON)
			res2 = (1.0e0_8/(liptime))*LOG(separations(2)/EPSILON)
			!PRINT *, res1, res2, state(1,:), e
			PRINT *, res1, res2, e
		END DO
		!WRITE(*,"(ES13.5, ES13.5)"), separations(1,i), separations(2,i)
		CALL resetState(state)
		stateRoundErr(2,1) = 0.0_8
		stateRoundErr(2,2) = 0.0_8
		stateRoundErr(2,3) = 0.0_8
		stateRoundErr(3,1) = 0.0_8
		stateRoundErr(3,2) = 0.0_8
		stateRoundErr(3,3) = 0.0_8
	END DO
END SUBROUTINE trackLyapunov

SUBROUTINE henonHeilesSurf(state)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), dimension(6) :: stateRoundErr
	real(kind=PREC), dimension(6) :: prevState
	real(kind=PREC) :: totalU, U
	integer :: numCross
	
	numCross = 0
	stateRoundErr = 0.0_8
	
	DO
		prevState = state
		CALL symplecticStep(state, stateRoundErr, dt, 0.0_8, U)
		CALL henonHeilesU(state(1), state(2), 0.0_8, totalU)
		IF(state(1) == 0.0_8 .AND. MOD(numCross, 100) == 0) THEN
			numCross = numCross + 1
			IF(MOD(numCross, 10) == 0) THEN
				PRINT *, state(2), state(5)
			END IF
			!PRINT *, totalU + (state(4)**2+state(5)**2)/2.0_8
		END IF
		IF(prevState(1) < 0.0_8 .AND. state(1) > 0.0_8) THEN
			numCross = numCross + 1
!			PRINT *, state(2), state(5)
			IF(MOD(numCross, 10) == 0) THEN
				PRINT *, prevState(2)&
					+ABS(prevState(1))&
					*(state(2)-prevState(2))/(state(1)-prevState(1))&
					,prevState(5)&
					+ABS(prevState(1))&
					*(state(5)-prevState(5))/(state(1)-prevState(1))
			END IF
			!PRINT *, totalU + (state(4)**2+state(5)**2)/2.0_8
		END IF
		IF(prevState(1) > 0.0_8 .AND. state(1) < 0.0_8) THEN
			numCross = numCross + 1
!			PRINT *, state(2), state(5)
			IF(MOD(numCross, 10) == 0) THEN
				PRINT *, prevState(2)&
					+ABS(prevState(1))&
					*(state(2)-prevState(2))/(-state(1)+prevState(1))&
					,prevState(5)&
					+ABS(prevState(1))&
					*(state(5)-prevState(5))/(-state(1)+prevState(1))
			END IF
			!PRINT *, totalU + (state(4)**2+state(5)**2)/2.0_8
		END IF
		IF(numCross == 1000000) THEN
			EXIT
		END IF
	END DO
END SUBROUTINE henonHeilesSurf

END MODULE
