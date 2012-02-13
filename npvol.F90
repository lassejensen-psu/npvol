Program NPVol

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! Calulate the volume of a nanoparticle using Monte Carlo integration.|
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

#ifdef _OPENMP
    use omp_lib ! Use the OpenMP module
#endif

! Increment preprocessor macro
#define INCREMENT(X) X = X + 1

    Implicit None

#ifdef _MPI
    include "mpif.h" ! Include the MPI headers
    Integer :: itmp
    Integer :: mpierr
    Integer :: myRank
    Integer :: nTasks
#endif

    Integer,     Parameter :: KINDR = KIND(0d0)        ! Double precision
    Integer,     Parameter :: NUM_PROBES_IN_LOOP = 100 ! Probes per sample
    Real(KINDR), Parameter :: ANGSTROM2NM = 1.0E-1_KINDR
    Real(KINDR), Parameter :: BOHR2ANGSTROM = 0.52917720859_KINDR
    Real(KINDR), Parameter :: ANGSTROM2BOHR = 1.0E0_KINDR / BOHR2ANGSTROM
    Real(KINDR), Parameter :: BOHR2NM       = BOHR2ANGSTROM * ANGSTROM2NM
    Real(KINDR), Parameter :: ZERO = 0.0E0_KINDR
    Real(KINDR), Parameter :: TWO  = 2.0E0_KINDR

    Integer :: time     ! The time to seed the random number generator
    Integer :: rate     ! Clock rate
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
    Real(KINDR) :: volume       ! Nanoparticle bolume
    Real(KINDR) :: probe(3,NUM_PROBES_IN_LOOP) ! The probe positions

    Logical :: isMom ! Is the mother processer
    Logical :: lfound(NUM_PROBES_IN_LOOP)  ! Is this probe in an atom?

    Character(1)   :: option   ! Command line options
    Character(80)  :: optarg   ! Command line arguments
    Character(80)  :: filename ! The XYZ filename
    Character(160) :: radline  ! The line containing the radii

!   Define a new type to hold the parameters of the nanoparticle
    Type GlobalType
        Logical                     :: lbohr
        Integer                     :: nAtoms
        Real(KINDR),    Allocatable :: rad(:)       ! Radius of each type
        Character(2),   Allocatable :: name(:)      ! The atom names
        Real(KINDR),    Allocatable :: xyz(:,:)     ! The coordinates
    End Type GlobalType

!   Create the global variable
    Type(GlobalType) :: g

!   ----------------------------------------
!   Start the clock and parallelization info
!   ----------------------------------------

#ifdef _MPI
    call MPI_INIT (mpierr)
    if (mpierr /= MPI_SUCCESS) call quit ('Error starting MPI')
    start_t = REAL(MPI_WTIME(), KINDR)

!   Get the total number of tasks and the rank of this node
    call MPI_COMM_RANK (MPI_COMM_WORLD, myRank, mpierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, nTasks, mpierr)

!   Determine if this task is the mother or not
    if (myRank == 0) then
        isMom = .true.
    else
        isMom = .false.
    end if
#else
#ifdef _OPENMP
    start_t = REAL(omp_get_wtime(), KINDR)
#else
    call SYSTEM_CLOCK (COUNT=time, COUNT_RATE=rate)
    start_t = REAL(time, KINDR) / REAL(rate, KINDR)
#endif
    isMom = .true. ! Always true for OpenMP and serial
#endif

!   ---------------------------
!   Interperet the command line
!   ---------------------------

    filename = ''
    nSamples = 0
    g%lbohr  = .false.

!   Loop over the command line arguments.
    do
        option = GetOpt('h?s:b')
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
    if (nSamples == 0) nSamples = MAX(g%nAtoms * 25, 10000)

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
#ifdef _MPI
    rng_state = rng_seed(932117 + time + myRank)

!   Determine how to distribute the samples over the nodes
    if (nTasks <= 1) then
        lowSamp  = 1
        highSamp = nSamples
    else
        lowSamp  = 1 + ( myRank * nSamples )      / nTasks
        highSamp = ( ( myRank + 1 ) *  nSamples ) / nTasks
    end if
#else
    rng_state = rng_seed(932117 + time)
    lowSamp  = 1
    highSamp = nSamples
#endif

!   ------------------------------
!   Perform the volume integration
!   ------------------------------

!   Determine the volume using a monte carlo algorithm.
    nHits = 0
!$OMP PARALLEL DO DEFAULT(SHARED)                              &
!$OMP             PRIVATE(probe, dist, iAtom, r, iTouch,       &
!$OMP                     lfound, hSum, atomRad, atomProbeRad) &
!$OMP             REDUCTION(+:nHits)
    samples_: do iSamp = lowSamp, highSamp

!       Grab a random number for the x, y, and z value of the probe(s)
!       Place these points in the box
        probe(:,1)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,2)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,3)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,4)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,5)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,6)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,7)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,8)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,9)   = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,10)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,11)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,12)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,13)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,14)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,15)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,16)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,17)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,18)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,19)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,20)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,21)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,22)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,23)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,24)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,25)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,26)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,27)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,28)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,29)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,30)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,31)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,32)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,33)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,34)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,35)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,36)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,37)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,38)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,39)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,40)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,41)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,42)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,43)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,44)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,45)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,46)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,47)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,48)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,49)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,50)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,51)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,52)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,53)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,54)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,55)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,56)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,57)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,58)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,59)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,60)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,61)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,62)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,63)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,64)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,65)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,66)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,67)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,68)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,69)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,70)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,70)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,71)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,72)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,73)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,74)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,75)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,76)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,77)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,78)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,79)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,80)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,80)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,81)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,82)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,83)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,84)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,85)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,86)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,87)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,88)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,89)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,90)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,90)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,91)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,92)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,93)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,94)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,95)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,96)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,97)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,98)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,99)  = box(:) * rng_uniform3(rng_state) + minbox(:)
        probe(:,100) = box(:) * rng_uniform3(rng_state) + minbox(:)

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

!           **********************
!           BEGIN MACRO DEFINITION
!           **********************

!           Define a macro that checks if this probe is in an atom or if it
!           touches more than three atoms.  If so, this probe is marked as a
!           hit.  This code is not executed if the current probe has already
!           been flagged as in/touching the nanoparticle.
!
!           Why the macro?
!           This is done as a way to unroll a static loop over some logic 
!           without cluttering this file with the same block over and over
!           and without incurring the overhead of an actual subroutine call
!           (which would defeat the purpose of unrolling the loop).  While
!           this is an incredibly ugly way to write code (kids, don't try
!           this at home!) it actually generates much faster code without
!           copy and pasting multiple times.
#define CHECK_IF_PROBE_TOUCHES_ATOM(A, P) \
            if (.not. lfound(P)) then;               \
                r(:) = g%xyz(:,A) - probe(:,P);      \
                dist = DOT_PRODUCT(r, r);            \
                if (dist <= atomRad) then;           \
                    INCREMENT(nHits);                \
                    INCREMENT(hSum);                 \
                    lfound(P) = .true.;              \
                else if (dist <= atomProbeRad) then; \
                    INCREMENT(iTouch(P));            \
                    if (iTouch(P) > 3) then;         \
                        INCREMENT(nHits);            \
                        INCREMENT(hSum);             \
                        lfound(P) = .true.;          \
                    end if;                          \
                end if;                              \
            end if

!           ********************
!           END MACRO DEFINITION
!           ********************

!           Unroll this loop 100 times.  Unrolling in this fashion is roughly
!           30% faster than using a standard loop from 1 to 100.
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 1)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 2)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 3)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 4)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 5)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 6)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 7)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 8)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 9)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 10)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 11)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 12)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 13)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 14)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 15)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 16)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 17)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 18)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 19)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 20)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 21)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 22)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 23)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 24)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 25)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 26)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 27)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 28)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 29)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 30)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 31)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 32)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 33)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 34)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 35)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 36)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 37)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 38)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 39)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 40)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 41)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 42)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 43)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 44)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 45)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 46)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 47)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 48)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 49)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 50)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 51)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 52)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 53)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 54)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 55)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 56)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 57)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 58)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 59)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 60)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 61)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 62)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 63)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 64)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 65)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 66)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 67)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 68)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 69)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 70)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 71)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 72)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 73)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 74)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 75)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 76)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 77)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 78)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 79)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 80)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 81)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 82)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 83)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 84)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 85)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 86)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 87)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 88)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 89)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 90)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 91)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 92)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 93)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 94)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 95)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 96)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 97)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 98)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 99)
            CHECK_IF_PROBE_TOUCHES_ATOM(iAtom, 100)

!           If all probes are hit, skip to next sample
            if (hSum == NUM_PROBES_IN_LOOP) exit atoms_

        end do atoms_

    end do samples_
!$OMP END PARALLEL DO
#ifdef _MPI
!   Collect the totals from all processes
    itmp = nHits
    call MPI_ALLREDUCE (itmp, nHits, 1, MPI_INTEGER, &
                        MPI_SUM, MPI_COMM_WORLD, mpierr)
    if (mpierr /= MPI_SUCCESS) call quit ('Could not reduce nHits')
#endif

!   --------
!   Clean up
!   --------

!   Finally, calculate the volume of the DIM system in nm**3
    volume = REAL(nHits, KINDR) / REAL(nSamples*NUM_PROBES_IN_LOOP, KINDR)
    volume = volume * totvol * BOHR2NM**3

!   End the clock
#ifdef _MPI
    tot_t = REAL(MPI_WTIME(), KINDR) - start_t
#else
#ifdef _OPENMP
    tot_t = REAL(omp_get_wtime(), KINDR) - start_t
#else
    call SYSTEM_CLOCK (COUNT=time, COUNT_RATE=rate)
    tot_t = ( REAL(time, KINDR) / REAL(rate, KINDR) ) - start_t
#endif
#endif

!   Write the results to screen
    if (isMom) then
        write(*,'(A14,1X,I10)')         '# Samples    :', nSamples
        write(*,'(A14,1X,ES10.4,1X,A)') 'Volume       :', volume, 'nm^3'
        write(*,'(A14,1X,F10.4,1X,A)')  'Elapsed Time :', tot_t, 'seconds'
    end if

    call quit ('NORMAL TERMINATION')

Contains

    Character Function GetOpt(optstr)

!   ========================================================================
!   Retrieves command line arguments.
!
!   This software is distributed under the following terms:
! 
!   Copyright (C) 2005 Dominik Epple
! 
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions
!   are met:
!   1. Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!   2. Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!   ========================================================================

        Character(*), Intent(In) :: optstr
        Integer                  :: argc
        Character(80)            :: arg
        Character(1)             :: okey
        Integer, Intrinsic       :: iargc
        Integer                  :: optstr_len
        Integer, Save            :: optstr_ind
        Integer, Save            :: opt_index = 1
        Logical                  :: lfound

!       Grab the number of arguments on the command line
!       This is a Fortran builtin
        argc = iargc()

!       If there are zero arguments left, then stop
        if (opt_index .gt. argc) then
           GetOpt = '>'
           return
        end if

!       Grab the next argument from the command line
!       This is a Fortran builtin
        call GetArg(opt_index, arg)

!       If the argument is an option specifier, look closer
        if (arg(1:1) .eq. '-') then

!           The option key is the Character after the dash
            okey = arg(2:2)
!           Initiallize the found switch and the index of the options string
            lfound = .false.
            optstr_ind = 1
!           Find options string length
            optstr_len = LEN_TRIM(optstr)

!           Search through the options string for the option
            do while(optstr_ind .le. optstr_len)
!               If the option is found, then change found switch and exit loop
                if (optstr(optstr_ind:optstr_ind) .eq. okey) then
                    lfound = .true.
!                   If the option carries an argument, grab the argument
                    if (optstr(optstr_ind + 1:optstr_ind + 1) .eq. ':') then
                        INCREMENT(optstr_ind)
                        INCREMENT(opt_index)
                        Call GetArg(opt_index, optarg)
                    end if
                    exit
                end if
                INCREMENT(optstr_ind)
            end do

!           If found, return the option.  If not found return warning.
            if (lfound) then
                getopt = okey
            else
                GetOpt = '!'
                optarg = arg
            end if

!       If the argument is not an option, it is a regular command line argument
        else
            GetOpt = '.'
            optarg = arg
        end if

        INCREMENT(opt_index)

        return

    End Function GetOpt

    Subroutine CollectCoordinates (filename)

!       ======================================================
!       Subroutine to collect the coordinates in the XYZ file.
!       ======================================================

        Character(80),  Intent(In) :: filename ! Name of the file
        Integer                    :: openstat ! File open status
        Integer                    :: readstat ! File read status
        Integer                    :: i        ! Counter
        Integer                    :: indx     ! Current index in a line
        Integer                    :: oldIndx  ! Current index in a line
        Integer                    :: nElem    ! Number of elements found
        Integer                    :: iElem    ! Element loop index
        Integer                    :: iAtom    ! Atom loop index
        Real(KINDR)                :: radii(100)    ! The radii detected
        Character(2)               :: elements(100) ! The elements detected
        Character(160)             :: line     ! A line of the file 
        Character(160)             :: radline  ! The line with the radii
        Logical                    :: lfound(100)   ! Was this element found

        Integer, Parameter         :: iufile = 100 ! File unit number

!       Open the file.
        if (isMom) then
            open(iufile, file=TRIM(filename), iostat=openstat, status='old')
            if (openstat /= 0) Then 
                call quit ("Error opening file '" // TRIM(filename) // "'.")
            end if

!           Read in the number of atoms
            read(iufile, *, iostat=readstat) g%nAtoms
            if (readstat /= 0) then
                call quit ("Error reading number of atoms from file '" // &
                           TRIM(filename) // "'.")
            end if

!           Read the line with the radii on it
            read(iufile,'(a160)', iostat=readstat) radline
            if (readstat /= 0) then
                call quit ("Error reading number of atomic radii from file '" // &
                           TRIM(filename) // "'.")
            else if (LEN_TRIM(radline) == 0) then
                call quit ("Missing atomic radii in file '" // TRIM(filename) // &
                           "'.")
            end if
        end if

#ifdef _MPI
        call MPI_BARRIER (MPI_COMM_WORLD, mpierr)
!       Broadcast the number of atoms and radline to all processors
        call MPI_BCAST (g%nAtoms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        if (mpierr /= MPI_SUCCESS) call quit ('Could not broadcast nAtoms')
        call MPI_BCAST (radline, 160, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
        if (mpierr /= MPI_SUCCESS) call quit ('Could not broadcast radline')
#endif

!       Allocate memory for the system
        allocate(g%name(g%nAtoms))
        allocate(g%xyz(3,g%nAtoms))
        allocate(g%rad(g%nAtoms))

        if (isMom) then
!           Read the coordinates.
            do i = 1, g%nAtoms
                read(iufile,'(a160)', iostat=readstat) line
!               Make sure there is more to read
                if (readstat /= 0) then
                    call quit ('Given no. of atoms larger than actual ' // &
                               "no. of atoms in '" // TRIM(filename) // "'.")
                end if
!               Read each atom's coordinates
                read(line,*) g%name(i), g%xyz(1:3,i)
            end do

!           Make sure there is no more to read
            read(iufile,'(a160)', iostat=readstat) line
            if (readstat == 0) then          
                call quit ('Actual no. of atoms larger than given no. of' // &
                           "atoms for '" // TRIM(filename) // "'.")
            end if

!           Close the file
            close(iufile)
        end if

#ifdef _MPI
        call MPI_BARRIER (MPI_COMM_WORLD, mpierr)
!       Broadcast the coordinates and atom names to all processers
        call MPI_BCAST (g%xyz, 3*g%nAtoms, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_WORLD, mpierr)
        if (mpierr /= MPI_SUCCESS) call quit ('Could not broadcast xyz')
        do iAtom = 1, g%nAtoms
            call MPI_BCAST (g%name(iAtom), 2, MPI_CHARACTER, 0, &
                            MPI_COMM_WORLD, mpierr)
            if (mpierr /= MPI_SUCCESS) call quit ('Could not broadcast name')
        end do
#endif

!       ----------------------------------------
!       Assign the radii in the proper locations
!       ----------------------------------------

!       Parse the radius line
        elements = ''
        radline  = TRIM(radline)
        radii    = ZERO
        nElem    = 0
        oldIndx  = 1
        do 
!           Find the next blank
            indx = INDEX(radline(oldIndx:), ' ') + ( oldIndx - 1 ) 
!           Quit if none were found
            if (indx == oldIndx) exit
!           Extract this parameter
            line = radline(oldIndx:indx-1)
!           Increment the old index
            oldIndx = indx + 1
!           From this line, find the '='
            indx = INDEX(line, '=')
!           It is an error if there is no '='
            if (indx == 0) then
                call quit ('Error reading atomic radii in file "' // &
                            TRIM(filename) // "'.")
            end if
!           Read the next element and the radii
            INCREMENT(nElem)
            elements(nElem) = line(1:indx-1)
            read(line(indx+1:), *, iostat=readstat) radii(nElem)
            if (readstat /= 0) then
                call quit ('Given radii not a float in file "' // &
                            TRIM(filename) // "'.")
            end if
        end do

!       If no elements were found, end in error
        if (nElem == 0) then
            call quit ('No atomic radii were found in file "' // &
                        TRIM(filename) // "'.")
        end if

!       Match the radii to the names found in the coordinates
        g%rad = ZERO
        atoms_: do iAtom = 1, g%nAtoms
            elements_: do iElem = 1, nElem
                if (g%name(iAtom) == elements(iElem)) then
                    g%rad(iAtom) = radii(iElem)
                    cycle atoms_
                end if
            end do elements_
        end do atoms_

!       Make sure all atoms got a radii
        if (ANY(g%rad == ZERO)) then
            call quit ('Not all atomic radii were defined in file "' // &
                        TRIM(filename) // "'.")
        end if

!       Convert to Angstroms if necessary
        if (.not. g%lbohr) then
            g%rad = g%rad * ANGSTROM2BOHR
            g%xyz = g%xyz * ANGSTROM2BOHR
        end if

    End Subroutine CollectCoordinates

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
!       largest posssible integer.
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

    Subroutine quit (msg)

!       =============================
!       Cleanly exits the calculation
!       =============================

        Character(*), Intent(In) :: msg
#ifdef _MPI
        Integer :: idummy
#endif

!       Deallocate the allocatables
        if (ALLOCATED(g%xyz)) deallocate(g%rad, g%name, g%xyz)

#ifdef _MPI
!       For MPI, use proper MPI routine to clean up.
        if (msg .eq. 'NORMAL TERMINATION') then
            call MPI_FINALIZE (mpierr)
        else
            if (isMom) write(*,*) msg
            call MPI_ABORT (MPI_COMM_WORLD, mpierr, idummy)
        end if
#else
!       If not MPI, force quit if not normally terminating
        if (msg .ne. 'NORMAL TERMINATION') then
            write(*,*) msg
        end if
#endif
        stop

    End Subroutine quit

    Subroutine Help ()

!       =========
!       Help info
!       =========

        if (isMom) then
            write(*,*)
            write(*,*) 'Program NPVol'
            write(*,*) 'Calcualte the volume of a nanoparticle using'
            write(*,*) 'Monte Carlo integration'
            write(*,*) 
            write(*,*) 'Simply give an XYZ file containing the coordinates'
            write(*,*) 'of your nanoparticle and the atomic radii and this '
            write(*,*) 'will determine the volume of your system.  The radii '
            write(*,*) 'are given on the second line of the file in the '
            write(*,*) 'comment section with the syntax "Ag=1.445 Au=1.445".'
            write(*,*) 'It is assumed that the units are Angstroms.'
            write(*,*)
            write(*,*) 'Options:'
            write(*,*) '-s NSAMPLES'
            write(*,*) '   Use this option to change the number of Monte Carlo'
            write(*,*) '   sample points that are used.  The default is based'
            write(*,*) '   on the number of atoms in the system'
            write(*,*)
            write(*,*) '-b'
            write(*,*) '   This changes the unit read in from Angstrom to Bohr'
            write(*,*)
            write(*,*) 'Author: Seth M. Morton'
            write(*,*)
        end if

        call quit ('NORMAL TERMINATION')

    End Subroutine Help

End Program NPVol
