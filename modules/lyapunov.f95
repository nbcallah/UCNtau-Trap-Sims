#include "constants.h"

MODULE lyapunov

CONTAINS

SUBROUTINE calcLyapunov(state, results)
	USE trackGeometry
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	!state will be an array;
	!first index: 1-ref, 2-traj_a, 3-traj_b
	!second index: 1-x, 2-y, 3-z, 4-px, 5-py, 6-pz
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	!results: lyap1, lyap2, rmsx, y, z, px, py, pz, x0, y0, z0, px0, py0, pz0, eStart, eEnd, maxE
	real(kind=PREC), dimension(17), intent(out) :: results
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	real(kind=PREC), dimension(2,6) :: delta
	real(kind=PREC), dimension(1,6) :: rejection
	real(kind=PREC), dimension(2,nsep) :: separations
	real(kind=PREC), dimension(1,6) :: qrms
	
	real(kind=PREC) eStart, eMax, eCur
	real(kind=PREC) dumpVar, tStart
	integer :: i, j, numSteps, rmsNum
	
	eMax = 0.0_8
    !results: lyap1, lyap2, rmsx, y, z, px, py, pz, x0, y0, z0, px0, py0, pz0, eStart, eEnd, maxE
    results(9:14) = state(1,:)
	
	CALL calcEnergy(state(1,:), results(15))
	numSteps = liptime/dt
	rmsNum = 0
	DO i=1,nsep,1
		DO j=1,numSteps,1
            CALL symplecticStep(state(1,:), dt, eCur)
            CALL symplecticStep(state(2,:), dt, dumpVar)
            CALL symplecticStep(state(3,:), dt, dumpVar)
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
!        PRINT *, separations(1,i)
		rejection(1,:) = delta(2,:)&
					-delta(1,:)&
						*(SUM(delta(1,:)*delta(2,:))/SUM(delta(1,:)**2))
		separations(2,i) = SQRT(SUM(rejection(1,:)**2))
		CALL resetState(state)
	END DO
	
    !results: lyap1, lyap2, rmsx, y, z, px, py, pz, x0, y0, z0, px0, py0, pz0, eStart, eEnd, maxE
	results(16) = eCur
	results(17) = eMax
	
	qrms = SQRT(qrms)
	
	results(3:8) = qrms(1,:)
	results(1) = (1.0e0_8/(nsep*liptime))*SUM(LOG(separations(1,:)/EPSILON))
	results(2) = (1.0e0_8/(nsep*liptime))*SUM(LOG(separations(2,:)/EPSILON))
END SUBROUTINE calcLyapunov

SUBROUTINE calcCleanTime(state, results)
	USE trackGeometry
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	!state will be an array;
	!first index: 1-ref, 2-traj_a, 3-traj_b
	!second index: 1-x, 2-y, 3-z, 4-px, 5-py, 6-pz
	real(kind=PREC), dimension(3,6), intent(inout) :: state
    real(kind=PREC), dimension(10), intent(out) :: results
    !results: t, x0, y0, z0, px0, py0, pz0, eStart, eEnd, eMax
	real(kind=PREC), dimension(3,6) :: stateRoundErr
	real(kind=PREC), dimension(2,6) :: delta
	real(kind=PREC), dimension(1,6) :: rejection
	real(kind=PREC), dimension(2,nsep) :: separations
	real(kind=PREC), dimension(1,6) :: qrms
	
	real(kind=PREC) eStart, eMax, eCur, settlingTime, t
	real(kind=PREC) dumpVar, tStart
	integer :: i, j, numSteps, rmsNum
    
    t = 0.0_8
    
    settlingTime = 0.0_8
	
	eMax = 0.0_8
    !results: t, x0, y0, z0, px0, py0, pz0, eStart, eEnd, eMax
    results(2:7) = state(1,:)
	
	CALL calcEnergy(state(1,:), results(8))
    IF (results(8) < 0.3445*GRAV*MASS_N) THEN
        results(1) = -1
        results(9) = results(8)
        results(10) = results(8)
        RETURN
    END IF

    numSteps = settlingTime/dt
    DO i=1,numSteps,1
        CALL symplecticStep(state(1,:), dt, eCur)
        t = t + dt
    END DO
    
    DO
        IF (state(1,3) >= (-1.5_8 + 0.38_8) .AND. state(1,2) > 0.0_8) THEN
!            PRINT *, state(1,3), (-1.5_8 + 0.38_8), state(1,2)
            results(1) = t - settlingTime
            results(9) = eCur
            results(10) = eMax
            RETURN
        END IF
        IF (t-settlingTime > 400.0_8) THEN
            results(1) = 400.0_8
            results(9) = eCur
            results(10) = eMax
            RETURN
        END IF
        CALL symplecticStep(state(1,:), dt, eCur)
        t = t + dt
        IF (eCur > eMax) THEN
            eMax = eCur
        END IF
    END DO
END SUBROUTINE calcCleanTime

END MODULE
