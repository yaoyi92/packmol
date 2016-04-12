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

! Routine to print the title 

subroutine title()

  write(*,"( 62('#'), /,/&
             ' PACKMOL - Packing optimization for the automated', /&
             ' generation of starting configurations for',        /&
             ' molecular dynamics. ',/&
             ' ',/&
             t42,' Version 16.103 ',/&
             ,/,62('#'),               /,/)")

end subroutine title
