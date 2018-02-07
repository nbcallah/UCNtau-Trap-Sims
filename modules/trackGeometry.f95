#include "constants.h"

MODULE trackGeometry

CONTAINS

SUBROUTINE randomPointTrap(x,y,z,px,py,pz)
    USE forcesAndPotential
    USE constants
    IMPLICIT NONE
    real(kind=PREC), intent(out) :: x, y, z, px, py, pz
    real(kind=PREC) :: zeta, maj_r, min_r, totalU, energy, max_p, maxEnergy, &
                    target_p, p_reject_val, p_len, en_reject_val
    real(kind=PREC) :: u1, u2, phi, theta !Used for lambertian generation of track directions
!    maxEnergy = GRAV*MASS_N*0.38_8 !Assume cleaning height of 38cm
    maxEnergy = GRAV*MASS_N*0.345_8 !Assume cleaning height of 38cm - 3.5cm to the zero energy point
    !A UCN of energy 34.5cm can reach 38cm w.r.t. the bottom of the trap.
    max_p = SQRT(2.0_8*MASS_N*maxEnergy)
    
    DO
        CALL RANDOM_NUMBER(energy)
        energy = maxEnergy * energy
!        IF (energy < GRAV*MASS_N*0.021_8) THEN
!            CYCLE
!        END IF
!        IF (energy > GRAV*MASS_N*0.1_8) THEN
!            EXIT
!        END IF
        IF (energy > GRAV*MASS_N*0.015_8) THEN
            EXIT
        END IF
!        CALL RANDOM_NUMBER(en_reject_val)
!        IF (en_reject_val < energy/(GRAV*MASS_N)/0.1) THEN
!            EXIT
!        END IF
    END DO

    DO
        z = -1.464413669130002_8
        CALL RANDOM_NUMBER(x)
        x = x*0.15_8 - 0.075_8
        CALL RANDOM_NUMBER(y)
        y = y*0.15_8 - 0.075_8

        CALL totalPotential(x, y, z, totalU)
        IF (totalU < energy) THEN
            EXIT
        END IF
    END DO
        
    target_p = SQRT(2.0_8*MASS_N*(energy - totalU))

    CALL RANDOM_NUMBER(u1)
    CALL RANDOM_NUMBER(u2)
    theta = ASIN(SQRT(u1))
    phi = 2.0_8 * PI * u2

    px = SIN(theta)*COS(phi)
    py = SIN(theta)*SIN(phi)
    pz = COS(theta)

    p_len = SQRT(px*px + py*py + pz*pz)

    px = (target_p/p_len)*px
    py = (target_p/p_len)*py
    pz = (target_p/p_len)*pz
!    PRINT *, x, y, z, px, py, pz, (px*px + py*py + pz*pz)/(2.0_8*MASS_N), totalU, theta
END SUBROUTINE randomPointTrap
!
!SUBROUTINE randomPointTrap(x, y, z, px, py, pz)
!    USE forcesAndPotential
!    USE constants
!    IMPLICIT NONE
!    real(kind=PREC), intent(out) :: x, y, z, px, py, pz
!    real(kind=PREC) :: zeta, maj_r, min_r, totalU, energy, max_p, maxEnergy, &
!                    target_p, p_reject_val, p_len, en_reject_val
!    maxEnergy = GRAV*MASS_N*0.38_8 !Assume cleaning height of 38cm
!    max_p = SQRT(2.0_8*MASS_N*maxEnergy)
!    DO
!        CALL RANDOM_NUMBER(energy)
!        energy = maxEnergy * energy
!        IF (energy > GRAV*MASS_N*0.1_8) THEN
!            EXIT
!        END IF
!        CALL RANDOM_NUMBER(en_reject_val)
!        IF (en_reject_val < energy/(GRAV*MASS_N)/0.1) THEN
!            EXIT
!        END IF
!    END DO
!    max_p = SQRT(2.0_8*MASS_N*energy)
!    !Sample the energy distribution
!    DO !Reject on P
!        DO !Reject on E
!            DO !Reject on inside/outside trap
!                CALL RANDOM_NUMBER(x)
!                x = x*1.6_8 - 1.0_8
!                CALL RANDOM_NUMBER(y)
!                y = y*2.0_8 - 1.0_8
!                CALL RANDOM_NUMBER(z)
!                z = z*0.5 - 1.5_8
!                IF (x > 0) THEN
!                    maj_r = 1.0
!                    min_r = 0.5
!                ELSE
!                    maj_r = 0.5
!                    min_r = 1.0
!                END IF
!
!                zeta = min_r - SQRT(x*x + (SQRT(y*y + z*z) - maj_r)*(SQRT(y*y + z*z) - maj_r))
!                IF (zeta >= 0.00001 .AND. zeta < min_r) THEN
!                    EXIT
!                END IF
!            END DO
!            CALL totalPotential(x, y, z, totalU)
!            IF (totalU < energy) THEN
!                EXIT
!            END IF
!        END DO
!        target_p = SQRT(2.0_8*MASS_N*(energy - totalU))
!        CALL RANDOM_NUMBER(p_reject_val)
!        IF (p_reject_val < (target_p/max_p)**3) THEN
!            EXIT
!        END IF
!    END DO
!    
!    CALL RANDOM_NUMBER(px)
!    CALL RANDOM_NUMBER(py)
!    CALL RANDOM_NUMBER(pz)
!    
!    p_len = SQRT(px*px + py*py + pz*pz)
!    
!    px = (target_p/p_len)*px
!    py = (target_p/p_len)*py
!    pz = (target_p/p_len)*pz    
!END SUBROUTINE randomPointTrap

SUBROUTINE cross(a, b, res)
    IMPLICIT NONE
    real(kind=PREC), dimension(3), intent(in) :: a, b
    real(kind=PREC), dimension(3), intent(out) :: res
    res(1) = a(2)*b(3)-a(3)*b(2)
    res(2) = a(1)*b(3)-a(3)*b(1)
    res(3) = a(1)*b(2)-a(2)*b(1)
END SUBROUTINE cross

END MODULE
