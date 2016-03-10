!  
!  Written by Leandro Martínez, 2009-2011.
!  Copyright (c) 2009-2011, Leandro Martínez, Jose Mario Martinez,
!  Ernesto G. Birgin.
!  
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  
!
! sizes.i: Define the maximum dimensions of the problems
!
!   maxatom:     Maximum number of atoms (total)
!   maxtype:     Maximum number of types of molecules (structures)
!   maxkeywords: Maximum number of keywords in input file
!   Maxlines:    Maximum number of lines of the input file
!   maxrest:     Maximum number of restrictions
!   mrperatom:   Maximum number of restrictions per atom
!   maxtry:      Number of tries for building the initial point  
!   nbp:         Maximum number of boxes for fast function evaluation (nbp**3)
!                FASTER FUNCTION EVALUATION IS OBTAINED WHEN
!                nbp**3 PARAMETER IS OF THE ORDER OF HALF OF THE NUMBER OF ATOMS
!                OF THE SYSTEM
!   nn:          Maximum number of variables 
!                (at least the number of molecules*6)
!

module sizes

  integer, parameter :: maxatom     =    500000
  integer, parameter :: maxtype     =        50
  integer, parameter :: maxkeywords =       200
  integer, parameter :: maxlines    =      1000
  integer, parameter :: maxrest     =       200
  integer, parameter :: mrperatom   =        10
  integer, parameter :: maxtry      =      1000
  integer, parameter :: nbp         =        84
  integer, parameter :: nn          = maxatom*6    

end module sizes

