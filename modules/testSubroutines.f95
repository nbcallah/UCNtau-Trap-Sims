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
    
    integer :: nDips = 4
    real(kind=PREC) :: speed
    real(kind=PREC), dimension(4) :: dipHeights
    real(kind=PREC), dimension(4) :: dipEnds
    real(kind=PREC) :: holdT
    
    integer :: i
    
    holdT = 900
    
!    dipHeights = (/0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010/) !9 dip
!    dipHeights = (/0.49, 0.250, 0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010/) !9 dip PSE
    dipHeights = (/0.49, 0.380, 0.250, 0.010/) !3 dip
!    dipHeights = (/0.49_8, 0.380_8, 0.250_8, 0.01_8/)
!    dipEnds =     (/0.0,  40.0,  80.0,  100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 300.0/) !9 dip
!    dipEnds =     (/0.0_8,  200.0_8,  200.0+holdT, 200.0+holdT+20.0, 200.0+holdT+40.0, 200.0+holdT+50.0, &
!                    200.0+holdT+60.0, 200.0+holdT+70.0, 200.0+holdT+80.0, 200.0+holdT+90.0, &
!                    200.0+holdT+100.0, 200.0+holdT+120.0/) !9 dip PSE
    dipEnds =     (/0.0,  20.0,  40.0,  140.0/) !3 dip
!    dipEnds =     (/0.0_8,  40.0_8,  400.0_8, 500.0_8/)
    
    IF (t > dipEnds(nDips)) THEN
        z = 0.01
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
!std::complex<double> wavenum(double ePerp, std::complex<double> u) {
!    return std::sqrt((2*MASS_N/(HBAR*HBAR))*(ePerp - u));
!}

FUNCTION gamma(kn, knm1) RESULT(ret)
    complex, intent(in) :: kn, knm1
    complex :: ret
    ret = knm1/kn
END FUNCTION gamma
!std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1) {
!    return knm1/kn;
!}

FUNCTION m(kn, knm1, z) RESULT(ret)
    complex, intent(in) :: kn, knm1, z
    complex, dimension(2,2) :: ret
    ret(1,1) = (1.0/2.0)*(1.0 + gamma(kn,knm1))*EXP((0.0,1.0)*(knm1-kn)*z);
    ret(1,2) = (1.0/2.0)*(1.0 - gamma(kn,knm1))*EXP(-(0.0,1.0)*(knm1+kn)*z);
    ret(2,1) = (1.0/2.0)*(1.0 - gamma(kn,knm1))*EXP((0.0,1.0)*(knm1+kn)*z);
    ret(2,2) = (1.0/2.0)*(1.0 + gamma(kn,knm1))*EXP(-(0.0,1.0)*(knm1-kn)*z);
END FUNCTION m
!std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z) {
!    std::vector<std::complex<double>> res = {std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0)};
!    res[0] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1-kn)*z);
!    res[1] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1+kn)*z);
!    res[2] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1+kn)*z);
!    res[3] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1-kn)*z);
!    return res;
!}

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
    
    ePerp_c = CMPLX(ePerp, 0.0)
    
    voxide = (2*PI*((HBAR/MASS_N)*HBAR))*ABORON*NBORONB2O3
    voxide = voxide + (2*PI*((HBAR/MASS_N)*HBAR))*AOXYGEN*NOXYGENB2O3
    woxide = (HBAR/2)*NBORONB2O3*2200*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN
    vboron = (2*PI*((HBAR/MASS_N)*HBAR))*ABORON*NBORON
    wboron = (HBAR/2)*NBORON*2200*SIGMABORON
    vzns = (2*PI*((HBAR/MASS_N)*HBAR))*AZINC*NZINC
    vzns = vzns + (2*PI*((HBAR/MASS_N)*HBAR))*ASULFUR*NSULFUR
    wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR
    
    pots(1) = (0,0)
    pots(2) = CMPLX(voxide, -woxide)
    pots(3) = CMPLX(vboron, -wboron)
    pots(4) = CMPLX(vzns, -wzns)
    zs(1) = (0.0, 0.0)
    !zs(2) = (3.5e-9, 0.0)
    zs(2) = (0.0, 0.0)
    !zs(3) = (3.5e-9 + 5.0e-9)
    zs(3) = (20e-9)
    zs(4) = (10000e-9)
    mbar(1,1) = (1, 0)
    mbar(1,2) = (0, 0)
    mbar(2,1) = (0, 0)
    mbar(2,2) = (1, 0)
    
    DO i = 4,2,-1
        mbar = MATMUL(mbar, m(wavenum(ePerp_c, pots(i)), wavenum(ePerp_c, pots(i-1)), zs(i-1)))
    END DO
    
    prob = 1.0_8 - REALPART(CONJG(-mbar(2,1)/mbar(2,2))*(-mbar(2,1)/mbar(2,2)))
END SUBROUTINE absorb
!double absorbProbQuantOxide(double ePerp, double thickOxide, double thickBoron) {
!    const double voxide = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORONB2O3 + (2*M_PI*(HBAR*HBAR)/MASS_N)*AOXYGEN*NOXYGENB2O3;
!    const double woxide = (HBAR/2)*NBORONB2O3*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN;
!    const double vboron = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORON;
!    const double wboron = (HBAR/2)*NBORON*SIGMABORON;
!    const double vzns = (2*M_PI*(HBAR*HBAR)/MASS_N)*AZINC*NZINC + (2*M_PI*(HBAR*HBAR)/MASS_N)*ASULFUR*NSULFUR;
!    const double wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR;
!    
!    std::vector<std::complex<double>> pots = {std::complex<double>(0, 0),
!                                              std::complex<double>(voxide, -woxide),
!                                              std::complex<double>(vboron, -wboron),
!                                              std::complex<double>(vzns, -wzns)};
!    std::vector<std::complex<double>> mbar = {std::complex<double>(1,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(1,0)};
!    std::vector<double> zs = {0.0, thickOxide*1e-9, thickOxide*1e-9 + thickBoron*1e-9, 10000e-9};
!    
!    for(int i = pots.size()-1; i > 0; i--) {
!        mbar = matmul(mbar, m(k(ePerp, pots[i]), k(ePerp, pots[i-1]), zs[i-1]));
!    }
!    
!    return 1.0 - (std::conj(-mbar[2]/mbar[3])*-mbar[2]/mbar[3]).real();
!}

SUBROUTINE trackDaggerHitTime(state)
    USE symplecticInt
    USE constants
    USE forcesAndPotential
    IMPLICIT NONE
    real(kind=PREC), dimension(6), intent(inout) :: state
    real(kind=PREC), dimension(6) :: prevState

    real(kind=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta
    real(kind=PREC) :: settlingTime
    real(kind=4), dimension(50) :: hitT
    real(kind=4), dimension(50) :: hitE
    
    integer :: i, numSteps, nHit
    
    nHit = 0
        
    hitT = 0.0_8
    hitE = 0.0_8
    
    t = 0.0_8
    
    !settlingTime = 20.0_8 + 200.0_8
    settlingTime = 20.0_8 + 50.0_8
    
    numSteps = settlingTime/dt
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
            
            CALL zOffDipCalc(t - settlingTime, zOff)
            IF (predX > 0.0_8) THEN
                zeta = 0.5_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 1.0_8)**2)
            ELSE
                zeta = 1.0_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 0.5_8)**2)
            END IF
            !TD offset from central axis: 6" ~0.1524m
            IF (predX > -0.3524_8 .AND. predX < 0.0476_8 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
                nHit = nHit + 1
                hitT(nHit) = t - settlingTime
                hitE(nHit) = state(5)*state(5)/(2.0_8*MASS_N)
                IF (nHit .EQ. 50) THEN
                    EXIT
                END IF
                IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE
                    PRINT *, "UHOH"
                END IF
!                WRITE(1) t - (20.0_8 + 50.0_8), predX, predZ - zOff
!                EXIT
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8) .AND. &
                    predZ < (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
                    ABS(predX + 0.1524_8) < (0.40_8 + 2.0179_8*(predZ + 1.5_8 - zOff - 0.2_8))/2.0_8) THEN
                IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                END IF
!                PRINT *, "BOUNCE LOWER"
!                PRINT *, predX, predZ, zOff
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
                    predZ < (-1.5_8 + zOff + 0.2_8 + 0.2667_8) .AND. &
                    ABS(predX + 0.1524_8) < 0.69215_8/2.0_8) THEN
                IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                END IF
!                PRINT *, "BOUNCE UPPER"
!                PRINT *, predX, predZ, zOff
            END IF
            
            IF (t > 2000) THEN
                EXIT
            END IF
        END IF
    END DO
    WRITE(1) energy, hitT, hitE
END SUBROUTINE trackDaggerHitTime

SUBROUTINE fixedEffDaggerHitTime(state)
    USE symplecticInt
    USE constants
    USE forcesAndPotential
    IMPLICIT NONE
    real(kind=PREC), dimension(6), intent(inout) :: state
    real(kind=PREC), dimension(6) :: prevState

    real(kind=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta
    real(kind=PREC) :: settlingTime
    
    real(kind=PREC) :: absProb, absU, deathTime
    
    integer :: i, numSteps, nHit, nHitHouseLow, nHitHouseHigh
    
    t = 0.0_8
    nHit = 0
    nHitHouseLow = 0
    nHitHouseHigh = 0
    
    !settlingTime = 20.0_8 + 200.0_8
    settlingTime = 50.0_8
    
    CALL RANDOM_NUMBER(deathTime)
    deathTime = -877.7*LOG(deathTime)
    
    numSteps = settlingTime/dt
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
                    WRITE(1) energy, t-settlingTime, state(5)*state(5)/(2.0_8*MASS_N), &
                        predX, 0.0_8, predZ, zOff, nHit, nHitHouseLow, nHitHouseHigh
                    EXIT
                END IF
                
                IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE
                    PRINT *, "UHOH"
                END IF
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8) .AND. &
                    predZ < (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
                    ABS(predX + 0.1524_8) < (0.40_8 + 2.0179_8*(predZ + 1.5_8 - zOff - 0.2_8))/2.0_8) THEN
                IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                END IF
                nHitHouseLow = nHitHouseLow + 1
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
                    predZ < (-1.5_8 + zOff + 0.2_8 + 0.2667_8) .AND. &
                    ABS(predX + 0.1524_8) < 0.69215_8/2.0_8) THEN
                IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
                    CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
                    CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
                    state = prevState
                END IF
                nHitHouseHigh = nHitHouseHigh + 1
            END IF
            IF (t-settlingTime > deathTime) THEN
                EXIT
            END IF
        END IF
    END DO
END SUBROUTINE fixedEffDaggerHitTime

SUBROUTINE trackEnergyGain(state, energy_start, energy_end, sympT, freq)
    USE symplecticInt
    USE constants
    USE forcesAndPotential
    IMPLICIT NONE
    real(kind=PREC), dimension(6), intent(inout) :: state
    real(kind=PREC), intent(out) :: energy_start, energy_end
    real(kind=PREC), optional, intent(inout) :: sympT
    real(kind=PREC), optional, intent(in) :: freq
    real(kind=PREC) :: eta, prevEta, t, energy, t_end
    integer :: i, numSteps, triggered, rising
    
    t = 0.0_8
    
    triggered = 0
    rising = 0
    prevEta = 10.0_8
    eta = 10.0_8

    numSteps = 1000e0/dt
!    totalU = 0.0_8
    
!    CALL calcEnergy(state, energy_start)
    energy_start = 0.0_8
!    energy = 0.0
    
    DO i=1,numSteps,1
        IF(present(sympT)) THEN
            CALL symplecticStep(state, dt, energy, sympT, freq)
        ELSE
            CALL symplecticStep(state, dt, energy)
            t = t+dt
        END IF
!        IF (triggered .EQ. 0) THEN
        prevEta = eta
        IF (state(1) > 0) THEN
            eta = 1.0_8 - SQRT(state(1)**2 + (SQRT(state(2)**2 + state(3)**2) - 0.5_8)**2)
        END IF
        IF (state(1) <= 0) THEN
            eta = 0.5_8 - SQRT(state(1)**2 + (SQRT(state(2)**2 + state(3)**2) - 1.0_8)**2)
        END IF
        IF ((rising .EQ. 0) .AND. (prevEta < eta)) THEN
            rising = 1
        END IF
        IF ((rising .EQ. 1) .AND. (prevEta > eta) .AND. (triggered .EQ. 0)) THEN
!            PRINT *, t, energy
            triggered = 1
            rising = 0
            energy_start = energy
        END IF
        IF ((rising .EQ. 1) .AND. (prevEta > eta) .AND. (triggered .EQ. 1)) THEN
!            PRINT *, t, energy
            t_end = t
            rising = 0
            energy_end = energy
        END IF
!        END IF
        IF (100.0_8*energy/(MASS_N*GRAV) > 38.0_8 + 5.0_8) THEN
!            PRINT *, "DEAD"
            EXIT
        END IF
    END DO
!    PRINT *, energy_start, energy_end
!    PRINT *, t_end, energy_end
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
!    PRINT *, energy_start, energy_end
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
        
!    numSteps = 250e0_8/dt
    numSteps = 10e0/dt
    t = 0.0_8
    totalU = 0.0_8
    
!    CALL calcEnergy(state, energy)
    energy = 0.0
    
    DO i=1,numSteps,1
!        IF(t > 475.0_8 .AND. t < 480.0_8) THEN
!        IF(INT(dt*10_8*i)-INT(dt*10_8*(i-1)) .NE. 0) THEN
        IF(1 .EQ. 1) THEN
!            energy = totalU + SUM(state(1,4:6)**2)/(2.0_8*MASS_N)

            PRINT *, dt*i, state(1), state(2), state(3),&
            state(4)/MASS_N, state(5)/MASS_N, state(6)/MASS_N, energy,&
            totalU!, fx, fy, fz

            !PRINT *, dt*i, energy
        END IF
        IF(present(sympT)) THEN
            CALL symplecticStep(state, dt, energy, t, 60.0_8)
        ELSE
            CALL symplecticStep(state, dt, energy)
!            CALL totalForce(state(1), state(2), state(3), fx_nate, fy_nate, fz_nate, e_nate)
!            CALL totalForceDan(state(1), state(2), state(3), fx_dan, fy_dan, fz_dan, e_dan)
!            PRINT *, t, (fx_nate-fx_dan)/fx_nate, (fy_nate-fy_dan)/fy_nate,&
!            (fz_nate-fz_dan)/fz_nate, (e_nate-e_dan)/e_nate
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
