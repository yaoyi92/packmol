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
! Subroutine resetboxes: Subroutine that resets the occupancy of
!                        linked cell boxes
!

subroutine resetboxes()
      
  use sizes
  use compute_data, only : latomfirst, latomfix, &
                           lboxfirst, lboxnext, hasfree
  implicit none
  integer :: i, j, k, ibox


  ! Reset boxes

  latomfirst = 0 ! Reset all boxes - three dimensional array

  ! Reset data for boxes that contain fixed atom

  ibox = lboxfirst
  do while( ibox > 0 ) 
    call ibox_to_ijk(ibox,i,j,k)
    latomfirst(i,j,k) = latomfix(i,j,k)
    hasfree(i,j,k) = .false.
    ibox = lboxnext(ibox)
  end do
  lboxfirst = 0

end subroutine resetboxes

