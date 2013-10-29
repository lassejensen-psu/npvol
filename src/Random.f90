Module Random

!   ========================================================================
!   Provides a thread-safe random number generator.  A 4-element integer
!   "state" array is generated from the rng_seed, and this array is given to
!   the rng_uniform function to generate a new random number.
!   ========================================================================

    use constants

    Implicit None
    Public

Contains

    Function rng_seed(seed) Result(state)

!       ===============================================================
!       Seeds the RNG using a single integer and a default seed vector.
!       Generates thread-safe random numbers using the algorithm
!       of Marsaglia and Zaman.
!       ===============================================================

        Integer, Intent(In) :: seed
        Integer             :: state(4)

        state = (/ seed, 362436069, 16163801, 1131199299 /)

    End Function rng_seed

    Function rng_uniform3(state) Result(u)

!       =====================================
!       Draws a uniform real vector on [0,1].
!       =====================================

        Integer, Intent(InOut) :: state(4)
        Real(KINDR)            :: u(3)
        Integer                :: imz

!       Make a random number
        imz = state(1) - state(3)

!       Compensate for a negative number by overflowing to be positive
        if (imz < 0) imz = imz + HUGE(imz)

        state(1) = state(2)
        state(2) = state(3)
        state(3) = imz
        state(4) = 69069 * state(4) + 1013904243
        imz = imz + state(4)

!       Compensate for a negative number by overflowing to be positive
        if (imz < 0) imz = imz + HUGE(imz)

!       Normalize the random number between 0 and 1 by dividing by the
!       largest possible integer.
        u(1) = REAL(imz, KINDR) / REAL(HUGE(imz), KINDR)

!       Repeat for 2 and 3
        imz = state(1) - state(3)
        if (imz < 0) imz = imz + HUGE(imz)
        state(1) = state(2)
        state(2) = state(3)
        state(3) = imz
        state(4) = 69069 * state(4) + 1013904243
        imz = imz + state(4)
        if (imz < 0) imz = imz + HUGE(imz)
        u(2) = REAL(imz, KINDR) / REAL(HUGE(imz), KINDR)

        imz = state(1) - state(3)
        if (imz < 0) imz = imz + HUGE(imz)
        state(1) = state(2)
        state(2) = state(3)
        state(3) = imz
        state(4) = 69069 * state(4) + 1013904243
        imz = imz + state(4)
        if (imz < 0) imz = imz + HUGE(imz)
        u(3) = REAL(imz, KINDR) / REAL(HUGE(imz), KINDR)

    End Function rng_uniform3

End Module Random
