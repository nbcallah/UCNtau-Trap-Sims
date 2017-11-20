#include "modules/constants.h"

!SUBROUTINE totalForce(x, y, z, fx, fy, fz, loops)
!	!Calculate force from B-field and add in gravity
!	IMPLICIT NONE
!	real(kind=PREC), intent(in) :: x, y, z
!	real(kind=PREC), intent(out) :: fx, fy, fz
!	real(kind=PREC), dimension(3,NWIRES+1), intent(in) :: loops
!	real(kind=PREC) :: k
!	
!	k=MASS_N
!	fy=0
!	fz=0
!	fx=-k*x
!END SUBROUTINE totalForce

PROGRAM track
!	USE mpi
	USE constants
	USE testSubroutines
	USE forcesAndPotential
	USE trackGeometry

	IMPLICIT NONE
	real(kind=PREC) :: x, y, z, fx, fy, fz, totalU, energy
	real(kind=PREC), allocatable :: states(:,:,:)
	character(len=64) :: arg
	integer :: i
	integer :: seedLen
	integer, dimension(32) :: rngSeed

	
	IF (IARGC() .NE. 3) THEN
		PRINT *, "Error! Not enough or too many arguments!"
		PRINT *, "x y z"
		CALL EXIT(0)
	END IF
	
	CALL GETARG(1, arg)
	READ(arg,*) x
	CALL GETARG(2, arg)
	READ(arg,*) y
	CALL GETARG(3, arg)
	READ(arg,*) z
	
!	dt = 0.000025_8
	dt = 0.0005_8
	
	minU = -2.4283243003838247e-26_8
	
	allocate(states(1,3,6))
	
	PI=4.0e0_8*ATAN(1.0e0_8)
	a(1)=.5153528374311229364e0_8
	a(2)=-.085782019412973646e0_8
	a(3)=.4415830236164665242e0_8
	a(4)=.1288461583653841854e0_8
	b(1)=.1344961992774310892e0_8
	b(2)=-.2248198030794208058e0_8
	b(3)=.7563200005156682911e0_8
	b(4)=.3340036032863214255e0_8
!	PRINT *, SUM(a), SUM(b)
!	a(1)=(1.0_8/6.0_8)&
!		*(2.0_8+2.0_8**(1.0_8/3.0_8)&
!			+2.0_8**(-1.0_8/3.0_8))
!	a(2)=(1.0_8/6.0_8)&
!		*(1.0_8-2.0_8**(1.0_8/3.0_8)&
!			-2.0_8**(-1.0_8/3.0_8))
!	a(3)=a(2)
!	a(4)=a(1)
!	b(1)=0.0_8
!	b(2)=1.0_8/(2.0_8-2.0_8**(1.0_8/3.0_8))
!	b(3)=1.0_8/(1.0_8-2.0_8**(2.0_8/3.0_8))
!	b(4)=b(2)

	CALL RANDOM_SEED(size=seedLen)
	IF (seedLen > 32) THEN
		PRINT *, "Error! The requested length of seed is too long"
		CALL EXIT(0)
	END IF
	!I'm not going to care about proper types since it's just for seed values
	rngSeed(1) = 4434
	DO i=2,seedLen,1
		rngSeed(i) = MOD((48271*rngSeed(i-1)), 2147483647)
	END DO
	CALL RANDOM_SEED(put=rngSeed(1:seedLen))	
	
	DO i=1,10000,1
		CALL randomPointTrap(states(1,1,1), states(1,1,2), states(1,1,3), states(1,1,4), states(1,1,5), states(1,1,6))
!		CALL trackEnergyGain(states(1,:,:), 0.0_8)
		CALL calcEnergy(states(1,1,:), energy)
		PRINT *, states(1,1,1), states(1,1,2), states(1,1,3),&
			states(1,1,4), states(1,1,5), states(1,1,6), energy/(GRAV*MASS_N)
	END DO
	
	
	
!	states(1,1,:) = (/x,y,z,0.1_8*MASS_N,0.06542_8*MASS_N,0.0_8/)
	
!	CALL trackAndPrint(states(1,:,:), 0.0_8)
	
	
	
	
!	CALL compPots()
!	CALL calcx0Mesh()
	
!	PRINT *, states(1,1,:)
	
!	CALL force(x, y, z, fx, fy, fz, totalU)
!	CALL force(states(1,1,1), states(1,1,2), states(1,1,3), fx, fy, fz, totalU)
!	PRINT *, x, y, z, fx, fy, fz, totalU
!	CALL potential(x, y, z, totalU)
!	CALL potential(states(1,1,1), states(1,1,2), states(1,1,3), totalU)
!	PRINT *, x, y, z, 0.0, 0.0, 0.0, totalU
	
END PROGRAM track
