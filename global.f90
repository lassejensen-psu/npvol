Module global

    use constants
    Implicit None
    Public
    Save

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

Contains

    Subroutine New (self)

!       ===================================
!       Allocate memory for the global type
!       ===================================

        Type(GlobalType), Intent(InOut) :: self

        allocate(self%name(self%nAtoms))
        allocate(self%xyz(3,self%nAtoms))
        allocate(self%rad(self%nAtoms))

    End Subroutine New

    Subroutine Delete (self)

!       =====================================
!       Deallocate memory for the global type
!       =====================================

        Type(GlobalType), Intent(InOut) :: self

        if (ALLOCATED(self%name)) deallocate(self%name)
        if (ALLOCATED(self%xyz))  deallocate(self%xyz)
        if (ALLOCATED(self%rad))  deallocate(self%rad)

    End Subroutine Delete

End Module global
