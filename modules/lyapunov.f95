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

END MODULE
