Subroutine CollectCoordinates (filename)

!   ======================================================
!   Subroutine to collect the coordinates in the XYZ file.
!   The results are saved in the global variable g.
!   ======================================================

    use constants
    use global,   ONLY : g, New
    use parallel, ONLY : IAmMom, broadcast, WaitForMom

    Implicit None

    Character(*),   Intent(In) :: filename      ! Name of the file
    Integer                    :: openstat      ! File open status
    Integer                    :: readstat      ! File read status
    Integer                    :: i             ! Counter
    Integer                    :: indx          ! Current index in a line
    Integer                    :: oldIndx       ! Current index in a line
    Integer                    :: nElem         ! Number of elements found
    Integer                    :: iElem         ! Element loop index
    Integer                    :: iAtom         ! Atom loop index
    Real(KINDR)                :: radii(100)    ! The radii detected
    Character(2)               :: elements(100) ! The elements detected
    Character(160)             :: line          ! A line of the file 
    Character(160)             :: radline       ! The line with the radii
    Logical                    :: lfound(100)   ! Was this element found

    Integer, Parameter         :: iufile = 100 ! File unit number

!   Open the file.
    if (IAmMom) then
        open(iufile, file=TRIM(filename), iostat=openstat, status='old')
        if (openstat /= 0) Then 
            call quit ("Error opening file '" // TRIM(filename) // "'.")
        end if

!       Read in the number of atoms
        read(iufile, *, iostat=readstat) g%nAtoms
        if (readstat /= 0) then
            call quit ("Error reading number of atoms from file '" // &
                       TRIM(filename) // "'.")
        end if

!       Read the line with the radii on it
        read(iufile,'(a160)', iostat=readstat) radline
        if (readstat /= 0) then
            call quit ("Error reading number of atomic radii from file '" // &
                       TRIM(filename) // "'.")
        else if (LEN_TRIM(radline) == 0) then
            call quit ("Missing atomic radii in file '" // TRIM(filename) // &
                       "'.")
        end if
    end if

!   Broadcast results to children then allocate arrays
    call WaitForMom
    call broadcast (g%nAtoms, 'nAtoms')
    call broadcast (160, radline, 'radline')
    call New (g)

    if (IAmMom) then
!       Read the coordinates.
        do i = 1, g%nAtoms
            read(iufile,'(a160)', iostat=readstat) line
!           Make sure there is more to read
            if (readstat /= 0) then
                call quit ('Given no. of atoms larger than actual ' // &
                           "no. of atoms in '" // TRIM(filename) // "'.")
            end if
!           Read each atom's coordinates
            read(line,*) g%name(i), g%xyz(1:3,i)
        end do

!       Make sure there is no more to read
        read(iufile,'(a160)', iostat=readstat) line
        if (readstat == 0) then          
            call quit ('Actual no. of atoms larger than given no. of' // &
                       "atoms for '" // TRIM(filename) // "'.")
        end if

!       Close the file
        close(iufile)
    end if

!   Broadcast results to children
    call WaitForMom
    call broadcast (3, g%nAtoms, g%xyz, 'xyz')
    do iAtom = 1, g%nAtoms
        call broadcast (2, g%name(iAtom), 'atom name')
    end do

!   ----------------------------------------
!   Assign the radii in the proper locations
!   ----------------------------------------

!   Parse the radius line
    elements = ''
    radline  = TRIM(radline)
    radii    = ZERO
    nElem    = 0
    oldIndx  = 1
    do 
!       Find the next blank
        indx = INDEX(radline(oldIndx:), ' ') + ( oldIndx - 1 ) 
!       Quit if none were found
        if (indx == oldIndx) exit
!       Extract this parameter
        line = radline(oldIndx:indx-1)
!       Increment the old index
        oldIndx = indx + 1
!       From this line, find the '='
        indx = INDEX(line, '=')
!       It is an error if there is no '='
        if (indx == 0) then
            call quit ('Error reading atomic radii in file "' // &
                        TRIM(filename) // "'.")
        end if
!       Read the next element and the radii
        nElem = nElem + 1
        elements(nElem) = line(1:indx-1)
        read(line(indx+1:), *, iostat=readstat) radii(nElem)
        if (readstat /= 0) then
            call quit ('Given radii not a float in file "' // &
                        TRIM(filename) // "'.")
        end if
    end do

!   If no elements were found, end in error
    if (nElem == 0) then
        call quit ('No atomic radii were found in file "' // &
                    TRIM(filename) // "'.")
    end if

!   Match the radii to the names found in the coordinates
    g%rad = ZERO
    atoms_: do iAtom = 1, g%nAtoms
        elements_: do iElem = 1, nElem
            if (g%name(iAtom) == elements(iElem)) then
                g%rad(iAtom) = radii(iElem)
                cycle atoms_
            end if
        end do elements_
    end do atoms_

!   Make sure all atoms got a radii
    if (ANY(g%rad == ZERO)) then
        call quit ('Not all atomic radii were defined in file "' // &
                    TRIM(filename) // "'.")
    end if

!   Convert to Angstroms if necessary
    if (.not. g%lbohr) then
        g%rad = g%rad * ANGSTROM2BOHR
        g%xyz = g%xyz * ANGSTROM2BOHR
    end if

End Subroutine CollectCoordinates
