Module constants

    Implicit None

!   Precision
    Integer,        Parameter :: KINDR = KIND(0d0)

    Real(KINDR),    Parameter :: ZERO        =  0.0e0_KINDR
    Real(KINDR),    Parameter :: ONE         =  1.0e0_KINDR
    Real(KINDR),    Parameter :: TWO         =  2.0e0_KINDR
    Real(KINDR),    Parameter :: THREE       =  3.0e0_KINDR
    Real(KINDR),    Parameter :: FOUR        =  4.0e0_KINDR
    Real(KINDR),    Parameter :: FIVE        =  5.0e0_KINDR
    Real(KINDR),    Parameter :: SIX         =  6.0e0_KINDR
    Real(KINDR),    Parameter :: SEVEN       =  7.0e0_KINDR
    Real(KINDR),    Parameter :: EIGHT       =  8.0e0_KINDR
    Real(KINDR),    Parameter :: NINE        =  9.0e0_KINDR
    Real(KINDR),    Parameter :: TEN         = 10.0e0_KINDR
    Real(KINDR),    Parameter :: ELEVEN      = 11.0e0_KINDR
    Real(KINDR),    Parameter :: TWELVE      = 12.0e0_KINDR
    Real(KINDR),    Parameter :: THIRTEEN    = 13.0e0_KINDR
    Real(KINDR),    Parameter :: FOURTEEN    = 14.0e0_KINDR
    Real(KINDR),    Parameter :: FIFTEEN     = 15.0e0_KINDR
    Real(KINDR),    Parameter :: SIXTEEN     = 16.0e0_KINDR
    Real(KINDR),    Parameter :: SEVENTEEN   = 17.0e0_KINDR
    Real(KINDR),    Parameter :: EIGHTEEN    = 18.0e0_KINDR
    Real(KINDR),    Parameter :: NINETEEN    = 19.0e0_KINDR
    Real(KINDR),    Parameter :: TWENTY      = 20.0e0_KINDR
    Real(KINDR),    Parameter :: TWENTYTWO   = 22.0e0_KINDR
    Real(KINDR),    Parameter :: TWENTYFOUR  = 24.0e0_KINDR
    Real(KINDR),    Parameter :: THIRTYTWO   = 32.0e0_KINDR
    Real(KINDR),    Parameter :: THIRTYSIX   = 36.0e0_KINDR
    Real(KINDR),    Parameter :: FORTYEIGHT  = 48.0e0_KINDR
    Real(KINDR),    Parameter :: SEVENTYTWO  = 72.0e0_KINDR
    Real(KINDR),    Parameter :: HUNDRED     =  1.0e2_KINDR
    Real(KINDR),    Parameter :: FOURTH      = 0.25e0_KINDR
    Real(KINDR),    Parameter :: HALF        =  0.5e0_KINDR
    Real(KINDR),    Parameter :: THREEFOURTH = 0.75e0_KINDR
    Real(KINDR),    Parameter :: THIRD       =  1.0e0_KINDR / 3.0e0_KINDR
    Real(KINDR),    Parameter :: TWOTHIRD    =  2.0e0_KINDR / 3.0e0_KINDR
    Real(KINDR),    Parameter :: SIXTH       =  1.0e0_KINDR / 6.0e0_KINDR
    Real(KINDR),    Parameter :: THREESECOND =  3.0e0_KINDR / 2.0e0_KINDR
    Real(KINDR),    Parameter :: FOURTHIRD   =  4.0e0_KINDR / 3.0e0_KINDR

    Real(KINDR),    Parameter :: PI            = 3.14159265358979323846_KINDR
    Real(KINDR),    Parameter :: SQRTPI        = 1.7724538509055159_KINDR
    Real(KINDR),    Parameter :: TWOINVPI      = TWO / PI
    Real(KINDR),    Parameter :: SQRTTWOINVPI  = 0.79788456080286541_KINDR
    Real(KINDR),    Parameter :: TWOINVSQRTPI  = TWO / SQRTPI
    Real(KINDR),    Parameter :: HART2NM       = 45.56335_KINDR
    Real(KINDR),    Parameter :: HART2WAVENUM  = 219474.629_KINDR
    Real(KINDR),    Parameter :: HART2EV       = 27.2113961_KINDR
    Real(KINDR),    Parameter :: HART2HZ       = 6.579683920729E15_KINDR
    Real(KINDR),    Parameter :: EV2HART       = 0.0367493088_KINDR
    Real(KINDR),    Parameter :: WAVENUM2HART  = 4.5563353E-6_KINDR
    Real(KINDR),    Parameter :: NM2HART       = 45.56335_KINDR
    Real(KINDR),    Parameter :: NM2EV         = 1239.0_KINDR
    Real(KINDR),    Parameter :: HBAR          = 6.5821189916e-16_KINDR !units = eV s
    Real(KINDR),    Parameter :: BOHR2ANGSTROM = 0.52917720859_KINDR
    Real(KINDR),    Parameter :: ANGSTROM2BOHR = 1.0E0_KINDR / BOHR2ANGSTROM
    Real(KINDR),    Parameter :: ANGSTROM2NM   = 1.0E-1_KINDR
    Real(KINDR),    Parameter :: NM2ANGSTROM   = 1.0E0_KINDR / ANGSTROM2NM
    Real(KINDR),    Parameter :: NM2M          = 1.0E-9_KINDR
    Real(KINDR),    Parameter :: BOHR2NM       = BOHR2ANGSTROM * ANGSTROM2NM
    Real(KINDR),    Parameter :: NM2BOHR       = 1.0E0_KINDR / BOHR2NM

    Complex(KINDR), Parameter :: I_C       = CMPLX(ZERO,    ONE,  KINDR)
    Complex(KINDR), Parameter :: ZERO_C    = CMPLX(ZERO,    ZERO, KINDR)
    Complex(KINDR), Parameter :: ONE_C     = CMPLX(ONE,     ZERO, KINDR)
    Complex(KINDR), Parameter :: TWO_C     = CMPLX(TWO,     ZERO, KINDR)
    Complex(KINDR), Parameter :: THREE_C   = CMPLX(THREE,   ZERO, KINDR)
    Complex(KINDR), Parameter :: FOUR_C    = CMPLX(FOUR,    ZERO, KINDR)
    Complex(KINDR), Parameter :: FIVE_C    = CMPLX(FIVE,    ZERO, KINDR)
    Complex(KINDR), Parameter :: SIX_C     = CMPLX(SIX,     ZERO, KINDR)
    Complex(KINDR), Parameter :: HALF_C    = CMPLX(HALF,    ZERO, KINDR)
    Complex(KINDR), Parameter :: HUNDRED_C = CMPLX(HUNDRED, ZERO, KINDR)

End Module constants
