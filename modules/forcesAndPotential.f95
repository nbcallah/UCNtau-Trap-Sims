#include "constants.h"

MODULE forcesAndPotential

CONTAINS

SUBROUTINE totalPotential(x, y, z, totalU, t)
	IMPLICIT NONE
	real(kind=PREC), intent(in) :: x, y, z
	real(kind=PREC), intent(out) :: totalU
	real(kind=PREC), optional, intent(in) :: t
	
	IF(PRESENT(t)) THEN
		CALL potential(x, y, z, totalU, t)
	ELSE
		CALL potential(x, y, z, totalU)
	END IF
END SUBROUTINE totalPotential


SUBROUTINE totalForce(x, y, z, fx, fy, fz, totalU, t)
	IMPLICIT NONE
	real(kind=PREC), intent(in) :: x, y, z
	real(kind=PREC), intent(out) :: fx, fy, fz
	real(kind=PREC), intent(out) :: totalU
	real(kind=PREC), optional, intent(in) :: t
	
	IF(PRESENT(t)) THEN
		CALL force(x, y, z, fx, fy, fz, t)
	ELSE
		CALL force(x, y, z, fx, fy, fz)
	END IF

END SUBROUTINE totalForce


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
