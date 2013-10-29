Module parallel

!   ======================================================================
!   This module will export all of mpif.h or omp_lib as well as a parallel
!   global variable, p.
!   ======================================================================

    use constants
#ifdef _OPENMP
    use omp_lib
#endif

    Implicit None
    Private

!   ----------------
!   Public variables
!   ----------------

    Public IAmMom

!   --------------------------------
!   Public subroutines and functions
!   --------------------------------

    Public startParallel
    Public endParallel
    Public isParallel
    Public real_time
    Public threadNum
    Public maxProcs
    Public nodeDistribute
    Public WaitForMom
    Public reduce
    Public broadcast

!   ------------------------
!   Local private varialbles
!   ------------------------

#ifdef _MPI
    include "mpif.h"
    Integer, Save      :: nProcs      ! Total number of processors
    Integer, Save      :: myRank      ! The rank of this node
    Integer, Save      :: mpierr      ! Most recent MPI error
#endif

    Integer, Parameter :: momRank = 0 ! The rank of the mother node

#ifndef _MPI
#ifndef _OPENMP
!   The clock start time. Used to catch timer overflow. Serial only
    Real(KINDR), Save :: starttime
#endif
#endif

!   -------------------------------------
!   The variable to determine if I am Mom
!   -------------------------------------

    Logical, Save :: IAmMom

!   ---------------------
!   The public interfaces
!   ---------------------

    Interface reduce
        Module Procedure reduceInt
    End Interface reduce

    Interface broadcast
        Module Procedure broadcastReal
        Module Procedure broadcastReal2
        Module Procedure broadcastReal3
        Module Procedure broadcastInteger0
        Module Procedure broadcastInteger1
        Module Procedure broadcastInteger2
        Module Procedure broadcastCharacter
    End Interface broadcast

Contains

Subroutine startParallel ()

!   ==============================================
!   Starts MPI and sets some appropriate variables
!   ==============================================

#ifdef _MPI
!   Start MPI
    call MPI_INIT (mpierr)
    call mpicheck ('Error starting MPI')

!   Get the total number of tasks and the rank of this node
    call MPI_COMM_RANK (MPI_COMM_WORLD, myRank, mpierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, nProcs, mpierr)

!   Determine if this task is the mother or not
    if (myRank == momRank) then
        IAmMom = .true.
    else
        IAmMom = .false.
    end if
#else
    
!   OpenMP and serial
    IAmMom = .true.
#ifndef _OPENMP
!   Mark start time in case of counter overflow if serial
    starttime = 0
    starttime = real_time()
#endif
    
#endif

End Subroutine startParallel

Subroutine endParallel (error)

!   =================================================
!   Kill in the best way depending on parellel method
!   =================================================

    Logical, Intent(In) :: error

#ifdef _MPI
    Integer :: idummy
    if (error) then
        call MPI_ABORT (MPI_COMM_WORLD, mpierr, idummy)
    else
        call MPI_FINALIZE (mpierr)
        call mpicheck ('Error finalizing MPI')
    end if
#else
!   Serial and OpenMP
    if (error) stop
#endif

End Subroutine endParallel

Function isParallel() Result (inPar)

!   ===============================================
!   Returns if current region is parallel or serial
!   ===============================================

    Logical :: inPar

#ifdef _MPI
!   MPI
    call MPI_INITIALIZED(inPar, mpierr)
    call mpicheck ('Error checking if MPI init is called')
#else
#ifdef _OPENMP
!   OpenMP
    inPar = omp_in_parallel()
#else
!   Serial
    inPar = .false.
#endif
#endif

End Function isParallel

Function threadNum() Result (tn)

!   =========================================
!   Returns the thread number for this thread
!   =========================================

    Integer :: tn

#ifdef _MPI
!   MPI
    tn = myRank
#else
#ifdef _OPENMP
!   OpenMP
    tn = OMP_get_thread_num()
#else
!   Serial
    tn = momRank
#endif
#endif

End Function threadNum

Function maxProcs() Result (np)

!   ============================================
!   Returns the max number of processors allowed
!   ============================================

    Integer :: np

#ifdef _MPI
!   MPI
    np = nProcs
#else
#ifdef _OPENMP
!   OpenMP
    np = OMP_get_max_threads()
#else
!   Serial
    np = 1
#endif
#endif

End Function maxProcs

Subroutine nodeDistribute (total, low, high)

!   ========================================
!   Distributes a number over the processors
!   ========================================

    Integer, Intent(In)    :: total
    Integer, Intent(Out)   :: low
    Integer, Intent(Out)   :: high

#ifdef _MPI
    if (nProcs <= 1) then
        low  = 1
        high = total
    else
        low  = 1 + ( myRank * total )     / nProcs
        high = ( ( myRank + 1 ) * total ) / nProcs
    end if
#else
    low  = 1
    high = total
#endif

End Subroutine nodeDistribute

Subroutine WaitForMom

!   ========================================
!   Wait for the mother node using a barrier
!   ========================================

#ifdef _MPI
    call MPI_BARRIER (MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error in MPI barrier')
#else
!   For serial or OpenMP do nothing
    return
#endif

End Subroutine WaitForMom

Function real_time () Result(time)

!   ==================================================================
!   Returns the system time in best way for the parallel method chosen
!   ==================================================================

    Real(KINDR) :: time

#ifndef _OPENMP
#ifndef _MPI
!   Serial only
    Integer :: tcount, trate, tmax
#endif
#endif

#ifdef _OPENMP
!   When using OpenMP, use it's timer
    time = REAL(OMP_get_wtime(), KINDR)
#else
#ifdef _MPI
!   When using MPI, use it's timer
    time = REAL(MPI_WTIME(), KINDR)
#else
!   For serial calculations, use the intrinsic SYSTEM_CLOCK
    call SYSTEM_CLOCK (COUNT=tcount, COUNT_RATE=trate, COUNT_MAX=tmax)
    time = REAL(tcount, KINDR) / REAL(trate, KINDR)
!   If the timer has overflowed, add the max to the current time
    if (time < starttime) time = time + tmax
#endif
#endif

End Function real_time

!   ****************************************
!   REDUCING ROUTINES HIDDEN IN AN INTERFACE
!   ****************************************

Subroutine reduceInt (num, name)

!   =======================
!   Reduce an integer array
!   =======================

    Integer,      Intent(InOut) :: num
    Character(*), Intent(In)    :: name

#ifdef _MPI
!   Sum up the number on all nodes to return the total value
    call MPI_ALLREDUCE (MPI_IN_PLACE, num, 1, MPI_INTEGER, &
                        MPI_SUM, MPI_COMM_WORLD, mpierr)   
    call mpicheck ('Error reducing ' // name // ' from all nodes')
#else
!   Serial and OpenMP do nothing
    return
#endif

End Subroutine reduceInt

!   ********************************************
!   BROADCASTING ROUTINES HIDDEN IN AN INTERFACE
!   ********************************************

Subroutine broadcastReal (n, array, name)

!   ======================
!   Broadcast a real array
!   ======================

    Integer,            Intent(In)    :: n
    Real(KINDR),        Intent(InOut) :: array(n)
    Character(*),       Intent(In)    :: name

#ifdef _MPI
!   Broadcast for MPI
    call MPI_BCAST (array, n, MPI_DOUBLE_PRECISION, &
                    momRank, MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error broadcasting ' // name // ' to children')
#else
!   Do nothing for serial of OpenMP
    return
#endif

End Subroutine broadcastReal

Subroutine broadcastReal2 (n, n2, array, name)

!   ======================
!   Broadcast a real array
!   ======================

    Integer,            Intent(In)    :: n
    Integer,            Intent(In)    :: n2
    Real(KINDR),        Intent(InOut) :: array(n,n2)
    Character(*),       Intent(In)    :: name

#ifdef _MPI
!   Broadcast for MPI
    call MPI_BCAST (array, n*n2, MPI_DOUBLE_PRECISION, &
                    momRank, MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error broadcasting ' // name // ' to children')
#else
!   Do nothing for serial of OpenMP
    return
#endif

End Subroutine broadcastReal2

Subroutine broadcastReal3 (n, n2, n3, array, name)

!   ======================
!   Broadcast a real array
!   ======================

    Integer,            Intent(In)    :: n
    Integer,            Intent(In)    :: n2
    Integer,            Intent(In)    :: n3
    Real(KINDR),        Intent(InOut) :: array(n,n2,n3)
    Character(*),       Intent(In)    :: name

#ifdef _MPI
!   Broadcast for MPI
    call MPI_BCAST (array, n*n2*n3, MPI_DOUBLE_PRECISION, &
                    momRank, MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error broadcasting ' // name // ' to children')
#else
!   Do nothing for serial of OpenMP
    return
#endif

End Subroutine broadcastReal3

Subroutine broadcastInteger0 (val, name)

!   =========================
!   Broadcast a integer array
!   =========================

    Integer,            Intent(InOut) :: val
    Character(*),       Intent(In)    :: name

#ifdef _MPI
!   Broadcast for MPI
    call MPI_BCAST (val, 1, MPI_INTEGER, &
                    momRank, MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error broadcasting ' // name // ' to children')
#else
!   Do nothing for serial of OpenMP
    return
#endif

End Subroutine broadcastInteger0

Subroutine broadcastInteger1 (n, array, name)

!   =========================
!   Broadcast a integer array
!   =========================

    Integer,            Intent(In)    :: n
    Integer,            Intent(InOut) :: array(n)
    Character(*),       Intent(In)    :: name

#ifdef _MPI
!   Broadcast for MPI
    call MPI_BCAST (array, n, MPI_INTEGER, &
                    momRank, MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error broadcasting ' // name // ' to children')
#else
!   Do nothing for serial of OpenMP
    return
#endif

End Subroutine broadcastInteger1

Subroutine broadcastInteger2 (n, n2, array, name)

!   =========================
!   Broadcast a integer array
!   =========================

    Integer,            Intent(In)    :: n
    Integer,            Intent(In)    :: n2
    Integer,            Intent(InOut) :: array(n,n2)
    Character(*),       Intent(In)    :: name

#ifdef _MPI
!   Broadcast for MPI
    call MPI_BCAST (array, n*n2, MPI_INTEGER, &
                    momRank, MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error broadcasting ' // name // ' to children')
#else
!   Do nothing for serial of OpenMP
    return
#endif

End Subroutine broadcastInteger2

Subroutine broadcastCharacter (n, array, name)

!   ===========================
!   Broadcast a character array
!   ===========================

    Integer,            Intent(In)    :: n
    Character(n),       Intent(InOut) :: array
    Character(*),       Intent(In)    :: name

#ifdef _MPI
!   Broadcast for MPI
    call MPI_BCAST (array, n, MPI_CHARACTER, &
                    momRank, MPI_COMM_WORLD, mpierr)
    call mpicheck ('Error broadcasting ' // name // ' to children')
#else
!   Do nothing for serial of OpenMP
    return
#endif

End Subroutine broadcastCharacter

!   *********************
!   ONLY FOR MPI, PRIVATE
!   *********************

#ifdef _MPI
Subroutine mpicheck (msg)

!   =====================================================================
!   Checks the MPI return code for an error, and signals if there was one
!   =====================================================================

    Character(*), Intent(In) :: msg

    if (mpierr /= MPI_SUCCESS) call quit (msg)

End Subroutine mpicheck
#endif

End Module parallel
