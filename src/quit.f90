Subroutine quit (msg)

!   =============================
!   Cleanly exits the calculation
!   =============================

    use global,   ONLY : g, Delete
    use parallel, ONLY : IAmMom, endParallel

    Character(*), Intent(In) :: msg

!   Deallocate the allocatables
    call Delete (g)

!   Report possible error and quit
    if (msg .eq. 'NORMAL TERMINATION') then
        call endParallel (.false.)
    else
        if (IAmMom) write(*,*) msg
        call endParallel (.true.)
    end if
    stop

End Subroutine quit

Subroutine Help ()

!   =========
!   Help info
!   =========

    use parallel, ONLY : IAmMom

    if (IAmMom) then
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
