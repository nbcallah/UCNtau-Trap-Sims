#include "constants.h"

MODULE forcesAndPotential
    real(kind=PREC) :: minU

CONTAINS

SUBROUTINE totalPotential(x, y, z, totalU, t)
    IMPLICIT NONE
    real(kind=PREC), intent(in) :: x, y, z
    real(kind=PREC), intent(out) :: totalU
    real(kind=PREC), optional, intent(in) :: t
    
    IF(PRESENT(t)) THEN
        CALL potential(x, y, z, totalU, t, 0.0)
    ELSE
        CALL potential(x, y, z, totalU, 0.0, 0.0)
    END IF
    totalU = totalU - minU
END SUBROUTINE totalPotential

SUBROUTINE totalForce(x, y, z, fx, fy, fz, totalU, t, freq)
    IMPLICIT NONE
    real(kind=PREC), intent(in) :: x, y, z
    real(kind=PREC), intent(out) :: fx, fy, fz
    real(kind=PREC), intent(out) :: totalU
    real(kind=PREC), optional, intent(in) :: t, freq
    
    IF(PRESENT(t)) THEN
        CALL force(x, y, z, fx, fy, fz, totalU, t, freq)
    ELSE
        CALL force(x, y, z, fx, fy, fz, totalU, 0, 0)
    END IF
    totalU = totalU - minU
END SUBROUTINE totalForce

!SUBROUTINE totalForceDan(x, y, z, fx, fy, fz, totalU, t)
!    IMPLICIT NONE
!    real(kind=PREC), intent(in) :: x, y, z
!    real(kind=PREC), intent(out) :: fx, fy, fz
!    real(kind=PREC), intent(out) :: totalU
!    real(kind=PREC), optional, intent(in) :: t
!    
!    IF(PRESENT(t)) THEN
!        CALL force_dan(x, y, z, fx, fy, fz, t)
!    ELSE
!        CALL force_dan(x, y, z, fx, fy, fz)
!    END IF
!
!END SUBROUTINE totalForceDan


SUBROUTINE calcEnergy(state, energy, t)
    IMPLICIT NONE
    real(kind=PREC), dimension(6), intent(in) :: state
    real(kind=PREC), intent(out) :: energy
    real(kind=PREC), optional, intent(in) :: t
    real(kind=PREC) :: totalU
    
    IF(PRESENT(t)) THEN
        CALL totalPotential(state(1), state(2), state(3), totalU, t)
    ELSE
        CALL totalPotential(state(1), state(2), state(3), totalU)
    END IF
    
    energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N)
END SUBROUTINE calcEnergy

END MODULE
