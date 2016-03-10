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
module molpa

  use sizes
  double precision :: xcart(maxatom,3) 
  double precision :: coor(maxatom,3) 
  double precision :: restpars(maxrest,9) 
  double precision :: sizemin(3),sizemax(3)
  double precision :: boxl(3)
  double precision :: fdist, frest 
  double precision :: rot_bound(maxtype,3,2)
  double precision :: radius(maxatom), radius_ini(maxatom)
  
  double precision :: scale, scale2
  double precision :: fatom(maxatom)
  double precision :: dmax(maxtype)
  double precision :: cmxmin(maxtype),cmymin(maxtype),cmzmin(maxtype)
  double precision :: cmxmax(maxtype),cmymax(maxtype),cmzmax(maxtype)

  integer :: nmols(maxtype)    
  integer :: natoms(maxtype) 
  integer :: idfirst(maxtype)
  integer :: nratom(maxatom)   
  integer :: iratom(maxatom,mrperatom) 
  integer :: ityperest(maxrest)  
  integer :: nboxes(3)  
  integer :: ibmol(maxatom)  
  integer :: ibtype(maxatom)  

  integer :: ntotmol, ntype, natfix, ntotat

  logical :: constrain_rot(maxtype,3)
  logical :: comptype(maxtype)
  logical :: init1, move

  integer :: latomnext(maxatom),latomfirst(0:nbp+1,0:nbp+1,0:nbp+1),&
             latomfix(0:nbp+1,0:nbp+1,0:nbp+1)

end module molpa
