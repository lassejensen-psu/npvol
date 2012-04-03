Module GetOptions

!   ========================================================================
!   This is a slight modification from function to a subroutine of a
!   command-line option parser.  It is mostly in its original form
!
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

    Implicit None
    Private

!   Public Subroutine
    Public GetOpt

!   Saved variables
    Integer, Save :: optstr_ind
    Integer, Save :: opt_index = 1

Contains

    Subroutine GetOpt (optstr, option, optarg)

        Character(*), Intent(In)    :: optstr
        Character(1), Intent(InOut) :: option
        Character(*), Intent(InOut) :: optarg
        Character(80)               :: arg
        Character(1)                :: okey
        Integer                     :: argc
        Integer                     :: optstr_len
        Logical                     :: lfound

        Integer                     :: iargc

!       Grab the number of arguments on the command line
!       This is a Fortran builtin
        argc = iargc()

!       If there are zero arguments left, then stop
        if (opt_index .gt. argc) then
           option = '>'
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
                        optstr_ind = optstr_ind + 1
                        opt_index = opt_index + 1
                        Call GetArg(opt_index, optarg)
                    end if
                    exit
                end if
                optstr_ind = optstr_ind + 1
            end do

!           If found, return the option.  If not found return warning.
            if (lfound) then
                option = okey
            else
                option = '!'
                optarg = arg
            end if

!       If the argument is not an option, it is a regular command line argument
        else
            option = '.'
            optarg = arg
        end if

        opt_index = opt_index + 1

        return

    End Subroutine GetOpt

End Module GetOptions
