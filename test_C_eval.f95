#include "modules/constants.h"

PROGRAM track
    USE mpi
    USE constants
    USE testSubroutines
    USE forcesAndPotential
    USE trackGeometry
    USE lyapunov

    IMPLICIT NONE
    real(kind=PREC) :: x, y, z, fx, fy, fz, totalU, energy, sympT
    real(kind=PREC) :: energy_start, energy_end, maxEgain
    real(kind=PREC) :: freq, height
    real(kind=PREC), allocatable :: states(:,:,:)
    real(kind=PREC) :: res_lyap(17)
    real(kind=PREC) :: res_clean(10)
    character(len=256) :: arg
    character(len=256) :: fName
    character(len=256) :: rankString
    integer :: i, j, k
    integer :: seedLen, seedOff
    integer, dimension(32) :: rngSeed
    integer :: rank, size, tag, next, from, ierr, workerIt, trajPerWorker
    integer :: ntraj

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    
    IF (IARGC() .NE. 3) THEN
        PRINT *, "Error! Not enough or too many arguments!"
        PRINT *, "timestep n_traj OUTFILE"
        CALL MPI_FINALIZE(ierr)
        CALL EXIT(0)
    END IF
    
    CALL GETARG(1, arg)
    READ(arg,*) dt
    CALL GETARG(2, arg)
    READ(arg,*) ntraj
    CALL GETARG(3, fName)
!    READ(arg, FORMAT=A) fName
    
    write (rankString, "(I0)") rank
    
    fName = TRIM(fName) // TRIM(rankString)
    
    PRINT *, fName
    
    OPEN(UNIT=1,FILE=fName, FORM='UNFORMATTED')
    
    trajPerWorker = ntraj/size
    
!    minU = -2.4283243003838247e-26_8
    minU = -2.390245661413933e-26_8 !For lambda = 0.0508, brem=1.4
    minU = -2.390352484438862e-26_8 !For lambda = 0.05114, brem=1.35
    
    allocate(states(ntraj,3,6))
        
    PI=4.0e0_8*ATAN(1.0e0_8)
    a(1)=.5153528374311229364e0_8
    a(2)=-.085782019412973646e0_8
    a(3)=.4415830236164665242e0_8
    a(4)=.1288461583653841854e0_8
    b(1)=.1344961992774310892e0_8
    b(2)=-.2248198030794208058e0_8
    b(3)=.7563200005156682911e0_8
    b(4)=.3340036032863214255e0_8
!    PRINT *, SUM(a), SUM(b)
!    a(1)=(1.0_8/6.0_8)&
!        *(2.0_8+2.0_8**(1.0_8/3.0_8)&
!            +2.0_8**(-1.0_8/3.0_8))
!    a(2)=(1.0_8/6.0_8)&
!        *(1.0_8-2.0_8**(1.0_8/3.0_8)&
!            -2.0_8**(-1.0_8/3.0_8))
!    a(3)=a(2)
!    a(4)=a(1)
!    b(1)=0.0_8
!    b(2)=1.0_8/(2.0_8-2.0_8**(1.0_8/3.0_8))
!    b(3)=1.0_8/(1.0_8-2.0_8**(2.0_8/3.0_8))
!    b(4)=b(2)

    CALL RANDOM_SEED(size=seedLen)
    seedOff = 0
    IF (seedLen + seedOff > 32) THEN
        PRINT *, "Error! The requested length of seed is too long"
        CALL EXIT(0)
    END IF
    !I'm not going to care about proper types since it's just for seed values
!    rngSeed(1) = -2015418360
!    rngSeed(1) = -392767427
    rngSeed = (/-1945520552, 519016354, -404253796, 1561684179,&
                -1722369288, -492083488, -1625858952, 1054014135,&
                2043038395, 880511665, -981405493, -1547842263,&
                -1112416185, 1800698465, -116679498, -1701772841,&
                -754606736, -201409009, -736179651, 473055375,&
                -1274105917, 1798798780, 606230487, -106513495,&
                700856297, 1350392983, -618623404, 977451019,&
                526790775, 2063245412, 178787983, 1953263558/)
!    rngSeed(1) = 4434
!    DO i=2,seedLen,1
!        rngSeed(i) = MOD((48271*rngSeed(i-1)), 2147483647)
!    END DO
    CALL RANDOM_SEED(put=rngSeed(1+seedOff:seedLen+seedOff))
    
!    DO i=1,100,1
!        CALL zOffDipCalc(i*(250.0_8/100.0_8), z)
!        PRINT *, z, i*(250.0_8/100.0_8)
!    END DO
        
    DO i=1,ntraj,1
        CALL randomPointTrap(states(i,1,1), states(i,1,2), states(i,1,3),&
            states(i,1,4), states(i,1,5), states(i,1,6))
        CALL createInitialTrajectories(states(i,:,:))
    END DO
    
    liptime = 5.0
	nsep = 100
    
    DO i=trajPerWorker*rank+1,trajPerWorker*(rank+1),1
!        CALL calcLyapunov(states(i,:,:), res_lyap)
!        WRITE(1), res_lyap
        CALL calcCleanTime(states(i,:,:), res_clean)
        WRITE(1), res_clean
    END DO
    
    CALL MPI_FINALIZE(ierr)
    
END PROGRAM track
