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
! Function that determines the length of a string
!
function length(string)

  implicit none
  integer :: length
  character(len=200) :: string
  
  length = 200
  do while(string(length:length) <= ' ')
    length = length - 1
    if ( length == 0 ) exit
  end do

end function length      
 
