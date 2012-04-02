Program NPVol

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! Calulate the volume of a nanoparticle using Monte Carlo integration.|
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

    use constants
    use global
    use parallel
    use GetOptions
    use Random

    Implicit None

    Integer, Parameter :: NUM_PROBES_IN_LOOP = 100 ! Probes per sample
    Integer, Parameter :: DEFAULT_NSAMPLES = 10**6 ! Min samples by default
    Integer, Parameter :: MULTIPLIER = 50          ! For samples test

    Integer :: time     ! For seeding the random number generator
    Integer :: iAtom    ! This atom
    Integer :: iSamp    ! This sample
    Integer :: nHits    ! Number of hits total
    Integer :: nSamples ! Number of Monte Carlo sample points
    Integer :: hSum     ! Number of hits found for this sample
    Integer :: readstat ! Read status
    Integer :: lowSamp  ! Low sample for this processor
    Integer :: highSamp ! High sample for this processor

    Integer :: rng_state(4)               ! Random number generator state
    Integer :: iTouch(NUM_PROBES_IN_LOOP) ! # atoms this probe touches

    Real(KINDR) :: start_t      ! Start time
    Real(KINDR) :: tot_t        ! Total time
    Real(KINDR) :: maxbox(3)    ! Max box dimensions
    Real(KINDR) :: minbox(3)    ! Min box dimensions
    Real(KINDR) :: box(3)       ! The box size
    Real(KINDR) :: totvol       ! Total volume of the box
    Real(KINDR) :: probeRad     ! Radius of the probe
    Real(KINDR) :: atomRad      ! Current atom's radius
    Real(KINDR) :: atomProbeRad ! Radius of current atom + probe
    Real(KINDR) :: r(3)         ! Vector distance
    Real(KINDR) :: dist         ! Scalar distance
    Real(KINDR) :: probe(3,NUM_PROBES_IN_LOOP) ! The probe positions
    Real(KINDR) :: volume       ! Nanoparticle volume
    Real(KINDR) :: opt_vol      ! Volume during optimization
    Real(KINDR) :: radmult      ! The multiplier to optimize
    Real(KINDR) :: stepsize     ! Stepsize in optimization

    Logical :: loptimize                   ! Perform parameter optimization?
    Logical :: lfound(NUM_PROBES_IN_LOOP)  ! Is this probe in an atom?

    Character(1)   :: option   ! Command line options
    Character(80)  :: optarg   ! Command line arguments
    Character(80)  :: filename ! The XYZ filename

!   ----------------------------------------
!   Start the clock and parallelization info
!   ----------------------------------------

    call startParallel ()
    start_t = real_time()

!   ---------------------------
!   Interperet the command line
!   ---------------------------

    filename  = ''
    nSamples  = 0
    g%lbohr   = .false.
    loptimize = .false.

!   Loop over the command line arguments.
    do
        call GetOpt ('h?s:bo', option, optarg)
        select case (option)
            case ('>')
                Exit
            case ('!')
                call quit ("Unknown option '" // TRIM(optarg) // &
                           "'.  Use '-h' or '-?' for help.")
            case ('s')
!               The number of samples
                read(optarg, *, iostat=readstat) nSamples
                if (readstat /= 0) then
                    call quit ('Given number of samples must be an integer')
                end if
            case ('o')
                loptimize = .true.
            case ('b')
                g%lbohr = .true.
            case ('h')
                call Help ()
            case ('?')
                call Help ()
            case ('.')
                if (filename .eq. '') then
                    filename = optarg
                else
                    call quit ('Only one XYZ file allowed.' // &
                               "'.  Use '-h' or '-?' for help.")
                end if
            case default
                call quit ('Unknown error in options parsing.')
        end select
    end do

!   Make sure a filename was given
    if (filename .eq. '') then
        call quit ('You must include an XYZ file for the coordinates.')
    end if

!   ----------------------------------------------------
!   Read the XYZ file and assign the radii for each atom
!   ----------------------------------------------------

    call CollectCoordinates (filename)

!   --------------------------
!   Set up for the calculation
!   --------------------------

!   If a certain number of samples was not requested, base the number of samples
!   on the number of atoms in the system.
    if (nSamples == 0) nSamples = MAX(g%nAtoms * MULTIPLIER, DEFAULT_NSAMPLES)
!   The below is equivalent to rounding up to the nearest
!   NUM_PROBES_IN_LOOP and then dividing by NUM_PROBES_IN_LOOP.
!   This is done to account for the unrolling in the main loop.
    nSamples = CEILING(REAL(nSamples) / NUM_PROBES_IN_LOOP)

!   Write the number of samples that will be used
    if (IAmMom) then
        write(*,'(A14,1X,I10)') '# Samples    :', nSamples * NUM_PROBES_IN_LOOP
    end if

!   The probeRad is defined such that if it is in the center of a cube with a
!   sphere on each corner, it will be slightly intersecting with each sphere.
!   First, find the largest radius of the spheres we will encounter.
    probeRad = MAXVAL(g%rad)
!   Now, calculate the smallest circle that will not fit in the void
    probeRad = ( SQRT(TWO * ( TWO * probeRad )**2) - TWO * probeRad )
    probeRad = probeRad + probeRad * 0.2_KINDR

!   Grab the corner points of the box that completely contain the system
    maxbox = MAXVAL(g%xyz, DIM=2)
    minbox = MINVAL(g%xyz, DIM=2)

!   Define the absolute box size.
    box(:) = ABS(maxbox(:) - minbox(:))

!   Use this to determine the box volume
    totvol = PRODUCT(box(:))

!   Seed the random number generator based on the clock time
!   Set the samples on this processor
    call SYSTEM_CLOCK (COUNT=time)
!$OMP PARALLEL PRIVATE(rng_state, probe, dist, iAtom, r, iTouch, &
!$OMP                  lfound, hSum, atomRad, atomProbeRad)
    rng_state = rng_seed(932117 + time + threadNum() * 694257)

!   Distribute the samples over the nodes   
    call nodeDistribute (nSamples, lowSamp, highSamp) 

!   ------------------------------
!   Perform the volume integration
!   ------------------------------

!   Determine the volume using a monte carlo algorithm.
    nHits = 0
!$OMP DO REDUCTION(+:nHits)
    samples_: do iSamp = lowSamp, highSamp

!       Grab a random number for the x, y, and z value of the probe(s)
!       Place these points in the box.  This loop is unrolled for speed.
!       The routine is placed in another file to save space here & for clarity
        include 'probe_positions.fh'

!       Reset the found probe index and the touching index
        lfound = .false.
        iTouch = 0
        hSum   = 0

!       If these points are in any atom, count them as in the system
!       Use the square of the distance to test this, not the square root, since
!       the square root slows this loop by about a factor of two.
        atoms_: do iAtom = 1, g%nAtoms

!           Calculate this atom's radius squared and (probe+atom)**2
            atomRad = g%rad(iAtom)
            atomProbeRad = ( probeRad + atomRad )**2
            atomRad = atomRad * atomRad

!           This block checks if this probe is in an atom or if it
!           touches more than three atoms.  If so, this probe is marked as a
!           hit.  This code is not executed if the current probe has already
!           been flagged as in/touching the nanoparticle.  It has been
!           unrolled becase doing so is about 30% faster than a standard loop
!           from 1 to 100.
!
!           Because this code takes 1600 lines, it has been moved to
!           another file and is included with the preprocesser. Subroutines
!           are not used because they introduce overhead to call and would
!           defeat the purpose of unrolling the loop.
            include 'check_if_probe_touches_atom.fh'

!           If all probes are hit, skip to next sample
            if (hSum == NUM_PROBES_IN_LOOP) exit atoms_

        end do atoms_

    end do samples_
!$OMP END DO
!$OMP END PARALLEL
    call reduce (nHits, 'nHits')

!   --------------------
!   Print out the volume
!   --------------------

!   Finally, calculate the volume of the DIM system in nm**3
    volume = REAL(nHits, KINDR) / REAL(nSamples*NUM_PROBES_IN_LOOP, KINDR)
    volume = volume * totvol * BOHR2NM**3

!   Write the results to screen
    if (IAmMom) then
        write(*,'(A14,1X,ES10.4,1X,A)') 'Volume       :', volume, 'nm^3'
    end if

!   -------------------------------------------------------------------------
!   Shall we try to optimize the multiplier?
!   This will find a value to multiply the atomic radius by such that the sum
!   of atomic volumes given by spheres with that adjusted volumes.  This
!   multiplier can be used in other codes that use the summation method to
!   quickly approximate nanoparticle volume.
!   -------------------------------------------------------------------------

    if (loptimize) then

!       Set initial radmult and the optimization stepsize
        radmult = ONE
        stepsize = HALF

!       Calculate the initial volume with the summation method
        opt_vol = SUM(FOURTHIRD * PI * ( radmult * g%rad )**3)

!       Continue to optimize until we reach the correct volume to 5%
        do while (( ABS(volume - opt_vol) / volume ) > 0.05_KINDR)

!           Only change the step size if we went past the mark, as in
!           stepsize positive and higher than mark or
!           stepsize negative and lower than mark
!           When we change, we switch sign and halve value
            if ((ABS(stepsize) == stepsize .and. opt_vol > volume) .or. &
                (ABS(stepsize) /= stepsize .and. opt_vol < volume)) then
                stepsize = -stepsize * HALF
            end if

!           Calculate new stepsize
            radmult = radmult + stepsize

!           Recalculate the volume
            opt_vol = SUM(FOURTHIRD * PI * ( radmult * g%rad )**3)

        end do

!       Print out the optimized radmult and the volume
        if (IAmMom) then
            write(*,'(A14,1X,F10.4)')       'Multiplier   :', radmult
            write(*,'(A14,1X,ES10.4,1X,A)') 'Volume (Sum) :', opt_vol, 'nm^3'
        end if

    end if

!   End the clock
    tot_t = real_time() - start_t

!   Write the results to screen
    if (IAmMom) then
        write(*,'(A14,1X,F10.4,1X,A)')  'Elapsed Time :', tot_t, 'seconds'
    end if

    call quit ('NORMAL TERMINATION')

End Program NPVol
