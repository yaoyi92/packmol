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
!  Subroutine that swaps indexes for packing molecules once
!  at a time
!

subroutine swaptype(n,x,itype,action)

  use sizes, only : nn
  use compute_data, only : ntype, comptype, nmols, ntotmol, init1
  use swaptypemod
  implicit none
  integer ::n, itype, ilubar, ilugan, i, action
  double precision :: x(nn)
  character(len=62), parameter :: dash1_line = "( 62('#') )" 

  ! Save original data

  if ( action == 0 ) then
    do i = 1, nn
      xfull(i) = x(i)
    end do
    ntemp = n
    ntottemp = ntotmol
    if(ntype.eq.1) then
      itype = 1
    else 
      itype = 0
    end if    
  end if

  ! Swapping data for packing this itype

  if ( action == 1 ) then
    if ( itype <= ntype ) then
      write(*,*)
      if ( .not. init1 ) write(*,dash1_line)
      write(*,*)
      write(*,*) ' Packing molecules of type ', itype
      write(*,*)
      if ( .not. init1 ) write(*,dash1_line)
      do i = 1, ntype
        if(i.eq.itype) then
          comptype(i) = .true.
        else
          comptype(i) = .false.
        end if
      end do
      n = nmols(itype) * 6
      ntotmol = nmols(itype)
      ilubar = 0
      do i = 1, itype - 1
        ilubar = ilubar + nmols(i) * 3
      end do
      ilubar = ilubar + 1
      ilugan = ntemp/2 + ilubar 
      do i = 1, n / 2
        x(i) = xfull(ilubar)
        x(i+n/2) = xfull(ilugan)
        ilubar = ilubar + 1
        ilugan = ilugan + 1
      end do
    else
      if(ntype.gt.1) then
        if ( .not. init1 ) then 
          write(*,*)
          write(*,dash1_line)
          write(*,*)
          write(*,*)' Solving the problem for all molecules together.'
          write(*,*)
          write(*,dash1_line)
          write(*,*)
        end if
      end if
    end if
  end if

  ! Restoring counters for packing next type of molecules

  if ( action == 2 ) then
    if ( itype <= ntype ) then
      ilubar = 0
      do i = 1, itype - 1
        ilubar = ilubar + nmols(i)*3
      end do
      ilubar = ilubar + 1
      ilugan = ntemp/2 + ilubar
      do i = 1, n/2
        xfull(ilubar) = x(i)
        xfull(ilugan) = x(i+n/2)
        ilubar = ilubar + 1
        ilugan = ilugan + 1
      end do
    end if
  end if

  ! If itype=ntype+1 restore original vectors and pack all molecules

  if ( action == 3 ) then
    if(itype.eq.ntype+1) then
      n = ntemp 
      ntotmol = ntottemp
      do i = 1, n
        x(i) = xfull(i)
      end do
      do itype = 1, ntype
        comptype(itype) = .true.
      end do
    end if                  
  end if

end subroutine swaptype





