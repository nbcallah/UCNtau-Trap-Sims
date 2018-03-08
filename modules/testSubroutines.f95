#include "constants.h"

MODULE testSubroutines

CONTAINS

SUBROUTINE zOffDipCalc(t, z)
    real(kind=PREC), intent(in) :: t
    real(kind=PREC), intent(out) :: z
    
    integer :: nDips = 5
    real(kind=PREC) :: speed
    real(kind=PREC), dimension(5) :: dipHeights
    real(kind=PREC), dimension(5) :: dipEnds
    real(kind=PREC) :: holdT
    
    integer :: i
    
    holdT = 1420_8

    dipHeights = (/0.49_8, 0.49_8, 0.380_8, 0.250_8, 0.010_8/) !3 dip
    dipEnds =     (/0.0_8, holdT, holdT+20.0_8,  holdT+40.0_8,  holdT+140.0_8/) !3 dip
    
    IF (t > dipEnds(nDips)) THEN
        z = 0.01_8
        RETURN
    END IF

    DO i=1,nDips,1
        IF (dipEnds(i) > t) THEN
            EXIT
        END IF
    END DO

    speed = SIGN(1.0_8, dipHeights(i-1) - dipHeights(i))*0.49_8/13.0_8
    
    z = dipHeights(i-1) - speed*(t-dipEnds(i-1))
    
    IF ((speed > 0 .AND. z < dipHeights(i)) .OR. (speed < 0 .AND. z > dipHeights(i))) THEN
        z = dipHeights(i)
    END IF
END SUBROUTINE zOffDipCalc

SUBROUTINE reflect(state, norm, tang)
    USE trackGeometry
    USE constants
    real(kind=PREC), dimension(6), intent(inout) :: state
    real(kind=PREC), dimension(3), intent(in) :: norm, tang
    real(kind=PREC) :: u1, u2, theta, phi, pN, pT, pTprime, pLen, pTarget
    real(kind=PREC), dimension(3) :: tangPrime, newPdir
    
    pTarget = SQRT(state(4)**2 + state(5)**2 + state(6)**2)
    
    CALL cross(norm, tang, tangPrime)
    
    CALL RANDOM_NUMBER(u1)
    CALL RANDOM_NUMBER(u2)
    theta = ASIN(SQRT(u1))
    phi = 2.0_8 * PI * u2

    pN = COS(theta)
    pT = SIN(theta)*COS(phi)
    pTprime = SIN(theta)*SIN(phi)
    
    newPdir = pN*norm + pT*tang + pTprime*tangPrime

    state(4) = newPdir(1)
    state(5) = newPdir(2)
    state(6) = newPdir(3)
    
    pLen = SQRT(state(4)**2 + state(5)**2 + state(6)**2)
    state(4) = state(4) * pTarget/pLen
    state(5) = state(5) * pTarget/pLen
    state(6) = state(6) * pTarget/pLen
END SUBROUTINE reflect

FUNCTION wavenum(ePerp, u) RESULT(ret)
    complex, intent(in) :: ePerp, u
    complex :: ret
    ret = sqrt((2*(MASS_N/HBAR)/HBAR)*(ePerp - u))
END FUNCTION wavenum

FUNCTION gamma(kn, knm1) RESULT(ret)
    complex, intent(in) :: kn, knm1
    complex :: ret
    ret = knm1/kn
END FUNCTION gamma

FUNCTION m(kn, knm1, z) RESULT(ret)
    complex, intent(in) :: kn, knm1, z
    complex, dimension(2,2) :: ret
    ret(1,1) = (1.0_8/2.0_8)*(1.0_8 + gamma(kn,knm1))*EXP((0.0_8,1.0_8)*(knm1-kn)*z);
    ret(1,2) = (1.0_8/2.0_8)*(1.0_8 - gamma(kn,knm1))*EXP(-(0.0_8,1.0_8)*(knm1+kn)*z);
    ret(2,1) = (1.0_8/2.0_8)*(1.0_8 - gamma(kn,knm1))*EXP((0.0_8,1.0_8)*(knm1+kn)*z);
    ret(2,2) = (1.0_8/2.0_8)*(1.0_8 + gamma(kn,knm1))*EXP(-(0.0_8,1.0_8)*(knm1-kn)*z);
END FUNCTION m

SUBROUTINE absorb(ePerp, prob)
    USE constants
    real(kind=PREC), intent(in) :: ePerp
    real(kind=PREC), intent(out) :: prob
    real(kind=8) :: voxide, woxide, vboron, wboron, vznd, wzns
    complex, dimension(4) :: pots
    complex, dimension(4) :: zs
    complex, dimension(2,2) :: mbar
    complex :: ePerp_c
    integer :: i
    
    ePerp_c = CMPLX(ePerp, 0.0_8)
    
    voxide = (2*PI*((HBAR/MASS_N)*HBAR))*ABORON*NBORONB2O3
    voxide = voxide + (2*PI*((HBAR/MASS_N)*HBAR))*AOXYGEN*NOXYGENB2O3
    woxide = (HBAR/2)*NBORONB2O3*2200*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN
    vboron = (2*PI*((HBAR/MASS_N)*HBAR))*ABORON*NBORON
    wboron = (HBAR/2)*NBORON*2200*SIGMABORON
    vzns = (2*PI*((HBAR/MASS_N)*HBAR))*AZINC*NZINC
    vzns = vzns + (2*PI*((HBAR/MASS_N)*HBAR))*ASULFUR*NSULFUR
    wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR
    
    pots(1) = (0.0_8,0.0_8)
    pots(2) = CMPLX(voxide, -woxide)
    pots(3) = CMPLX(vboron, -wboron)
    pots(4) = CMPLX(vzns, -wzns)
    zs(1) = (0.0_8, 0.0_8)
    !zs(2) = (3.5e-9, 0.0)
    zs(2) = (0.0_8, 0.0_8)
    !zs(3) = (3.5e-9 + 5.0e-9)
    zs(3) = (4.615385e-9_8)
    zs(4) = (10000e-9_8)
    mbar(1,1) = (1.0_8, 0.0_8)
    mbar(1,2) = (0.0_8, 0.0_8)
    mbar(2,1) = (0.0_8, 0.0_8)
    mbar(2,2) = (1.0_8, 0.0_8)
    
    DO i = 4,2,-1
        mbar = MATMUL(mbar, m(wavenum(ePerp_c, pots(i)), wavenum(ePerp_c, pots(i-1)), zs(i-1)))
    END DO
    
    prob = 1.0_8 - REALPART(CONJG(-mbar(2,1)/mbar(2,2))*(-mbar(2,1)/mbar(2,2)))
END SUBROUTINE absorb

SUBROUTINE fixedEffDaggerHitTime(state)
    USE symplecticInt
    USE constants
    USE forcesAndPotential
    IMPLICIT NONE
    real(kind=PREC), dimension(6), intent(inout) :: state
    real(kind=PREC), dimension(6) :: prevState

    real(kind=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta
    real(kind=PREC) :: theta, eStart
    real(kind=PREC) :: settlingTime
    
    real(kind=PREC) :: absProb, absU, deathTime
    
    integer :: i, numSteps, nHit, nHitHouseLow, nHitHouseHigh
    
    theta = ACOS(state(6)/SQRT(state(4)**2 + state(5)**2 + state(6)**2))
    CALL calcEnergy(state, eStart)
    
    t = 0.0_8
    nHit = 0
    nHitHouseLow = 0
    nHitHouseHigh = 0
    
    settlingTime = 50.0_8+20.0_8
    
    CALL RANDOM_NUMBER(deathTime)
    deathTime = -877.7_8*LOG(deathTime)
    IF (deathTime < 1400) THEN
        WRITE(1) energy, theta, 0.0_8, state(5)*state(5)/(2.0_8*MASS_N), &
            state(1), state(2), state(3), -1.0_8, nHit, nHitHouseLow, nHitHouseHigh, &
            (energy-eStart)/eStart, deathTime
        RETURN
    END IF
    
    numSteps = settlingTime/dt
    DO i=1,numSteps,1
        CALL symplecticStep(state, dt, energy)
        t = t + dt
    END DO
    
    DO
        prevState = state
        CALL symplecticStep(state, dt, energy)
        t = t + dt
        IF (t-settlingTime > deathTime) THEN
            WRITE(1) energy, theta, t-settlingTime, state(5)*state(5)/(2.0_8*MASS_N), &
                state(1), state(2), state(3), -1.0_8, nHit, nHitHouseLow, nHitHouseHigh, &
                (energy-eStart)/eStart, deathTime
            EXIT
        END IF
        IF (SIGN(1.0_8, state(2)) .NE. SIGN(1.0_8, prevState(2))) THEN
            fracTravel = ABS(prevState(2))/(ABS(state(2)) + ABS(prevState(2)))
            predX = prevState(1) + fracTravel * (state(1) - prevState(1))
            predZ = prevState(3) + fracTravel * (state(3) - prevState(3))
            
            CALL zOffDipCalc(t - settlingTime, zOff)
            IF (predX > 0.0_8) THEN
                zeta = 0.5_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 1.0_8)**2)
            ELSE
                zeta = 1.0_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 0.5_8)**2)
            END IF
            !TD offset from central axis: 6" ~0.1524m
            IF (predX > -0.3524_8 .AND. predX < 0.0476_8 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
                nHit = nHit + 1
                CALL absorb(state(5)*state(5)/(2.0_8*MASS_N), absProb)
                CALL RANDOM_NUMBER(absU)
                IF (absU < absProb) THEN
                    WRITE(1) energy, theta, t-settlingTime, state(5)*state(5)/(2.0_8*MASS_N), &
                        predX, 0.0_8, predZ, zOff, nHit, nHitHouseLow, nHitHouseHigh, &
                        (energy-eStart)/eStart, deathTime
                    EXIT
                END IF
                
                IF (prevState(2) > 0.0_8 .AND. prevState(5) < 0.0_8) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0.0_8 .AND. prevState(5) > 0.0_8) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE
                    PRINT *, "UHOH"
                END IF
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8) .AND. &
                    predZ < (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
                    ABS(predX + 0.1524_8) < (0.40_8 + 2.0179_8*(predZ + 1.5_8 - zOff - 0.2_8))/2.0_8) THEN
                IF (prevState(2) > 0.0_8 .AND. prevState(5) < 0.0_8) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0.0_8 .AND. prevState(5) > 0.0_8) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                END IF
                nHitHouseLow = nHitHouseLow + 1
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
                    predZ < (-1.5_8 + zOff + 0.2_8 + 0.2667_8) .AND. &
                    ABS(predX + 0.1524_8) < 0.69215_8/2.0_8) THEN
                IF (prevState(2) > 0.0_8 .AND. prevState(5) < 0.0_8) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0.0_8 .AND. prevState(5) > 0.0_8) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                END IF
                nHitHouseHigh = nHitHouseHigh + 1
            END IF
        END IF
    END DO
END SUBROUTINE fixedEffDaggerHitTime

END MODULE
