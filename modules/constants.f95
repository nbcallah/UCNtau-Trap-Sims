#include "constants.h"

MODULE constants
    IMPLICIT NONE
    real(kind=PREC) :: PI
    real(kind=PREC), dimension(4) :: a
    real(kind=PREC), dimension(4) :: b
    real(kind=PREC) :: dt, liptime
    integer :: nsep
    SAVE
END MODULE