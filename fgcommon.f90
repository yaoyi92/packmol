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
! This file contains the parts of the evaluation of the function
! and gradient which are identical in the serial and parallel
! versions.
!
!
!
! Objective function evaluation
!
!
!
! Subroutine comprest: Compute the function value relative to
!                      to the restrictions for one atom
!
!

subroutine comprest(xcart,restpars,scale,scale2,nratom,ityperest,iratom,icart,f)

  use sizes
  implicit none

  integer :: iratcount, irest, icart 
  integer :: nratom(maxatom), ityperest(maxrest)
  integer :: iratom(maxatom,mrperatom)

  double precision :: xcart(maxatom,3), restpars(maxrest,9)
  double precision :: xmin, ymin, zmin, clength, a1, a2, a3, a4, w, b1, b2, b3, d, a5, a6
  double precision :: scale, scale2
  double precision :: f
  double precision :: xmax, ymax, zmax
  double precision :: v1, v2, v3
  double precision :: vnorm

  f = 0.d0
  do iratcount = 1, nratom(icart)
    irest = iratom(icart,iratcount)
    if(ityperest(irest).eq.2) then
      clength = restpars(irest,4)
      xmin = restpars(irest,1)
      ymin = restpars(irest,2)
      zmin = restpars(irest,3)
      xmax = restpars(irest,1) + clength
      ymax = restpars(irest,2) + clength
      zmax = restpars(irest,3) + clength
      a1 = dmin1(xcart(icart,1) - xmin, 0.d0)
      a2 = dmin1(xcart(icart,2) - ymin, 0.d0)
      a3 = dmin1(xcart(icart,3) - zmin, 0.d0)
      f = f + scale*(a1 * a1 + a2 * a2 + a3 * a3) 
      a1 = dmax1(xcart(icart,1) - xmax, 0.d0)
      a2 = dmax1(xcart(icart,2) - ymax, 0.d0)
      a3 = dmax1(xcart(icart,3) - zmax, 0.d0)
      f = f + scale*(a1 * a1 + a2 * a2 + a3 * a3)   
    else if(ityperest(irest).eq.3) then
      xmin = restpars(irest,1)
      ymin = restpars(irest,2)
      zmin = restpars(irest,3)
      xmax = restpars(irest,4)
      ymax = restpars(irest,5)
      zmax = restpars(irest,6)
      a1 = dmin1(xcart(icart,1) - xmin, 0.d0)
      a2 = dmin1(xcart(icart,2) - ymin, 0.d0)
      a3 = dmin1(xcart(icart,3) - zmin, 0.d0)
      f = f + scale*(a1 * a1 + a2 * a2 + a3 * a3) 
      a1 = dmax1(xcart(icart,1) - xmax, 0.d0)
      a2 = dmax1(xcart(icart,2) - ymax, 0.d0)
      a3 = dmax1(xcart(icart,3) - zmax, 0.d0)
      f = f + scale*(a1 * a1 + a2 * a2 + a3 * a3)   
    else if(ityperest(irest).eq.4) then     
      w = (xcart(icart,1)-restpars(irest,1))**2 + &
          (xcart(icart,2)-restpars(irest,2))**2 + &
          (xcart(icart,3)-restpars(irest,3))**2 - &
          restpars(irest,4)**2
      a1 = dmax1(w,0.d0)
      f = f + scale2*a1*a1
    else if(ityperest(irest).eq.5) then
      a1 = (xcart(icart,1)-restpars(irest,1))**2 / restpars(irest,4)**2
      a2 = (xcart(icart,2)-restpars(irest,2))**2 / restpars(irest,5)**2
      a3 = (xcart(icart,3)-restpars(irest,3))**2 / restpars(irest,6)**2
      a4 = restpars(irest,7)**2
      w = a1 + a2 + a3 - a4
      a1 = dmax1(w,0.d0)
      f = f + scale2*a1*a1
    else if(ityperest(irest).eq.6) then
      xmin = restpars(irest,1)
      ymin = restpars(irest,2)
      zmin = restpars(irest,3)
      xmax = restpars(irest,1) + restpars(irest,4)
      ymax = restpars(irest,2) + restpars(irest,4)
      zmax = restpars(irest,3) + restpars(irest,4)
      a1 = dmax1(xcart(icart,1) - xmin,0.d0)
      a2 = dmax1(xcart(icart,2) - ymin,0.d0)
      a3 = dmax1(xcart(icart,3) - zmin,0.d0)
      a4 = dmax1(xmax - xcart(icart,1),0.d0)
      a5 = dmax1(ymax - xcart(icart,2),0.d0)
      a6 = dmax1(zmax - xcart(icart,3),0.d0)
      f = f + a1*a2*a3*a4*a5*a6
    else if(ityperest(irest).eq.7) then
      xmin = restpars(irest,1)
      ymin = restpars(irest,2)
      zmin = restpars(irest,3)
      xmax = restpars(irest,4)
      ymax = restpars(irest,5)
      zmax = restpars(irest,6)
      a1 = dmax1(xcart(icart,1) - xmin,0.d0)
      a2 = dmax1(xcart(icart,2) - ymin,0.d0)
      a3 = dmax1(xcart(icart,3) - zmin,0.d0)
      a4 = dmax1(xmax - xcart(icart,1),0.d0)
      a5 = dmax1(ymax - xcart(icart,2),0.d0)
      a6 = dmax1(zmax - xcart(icart,3),0.d0)
      f = f + a1*a2*a3*a4*a5*a6
    else if(ityperest(irest).eq.8) then
      w = (xcart(icart,1)-restpars(irest,1))**2 + &
          (xcart(icart,2)-restpars(irest,2))**2 + &
          (xcart(icart,3)-restpars(irest,3))**2 - &
          restpars(irest,4)**2
      a1 = dmin1(w,0.d0)
      f = f + scale2*a1*a1
    else if(ityperest(irest).eq.9) then
      a1 = (xcart(icart,1)-restpars(irest,1))**2 / restpars(irest,4)**2
      a2 = (xcart(icart,2)-restpars(irest,2))**2 / restpars(irest,5)**2
      a3 = (xcart(icart,3)-restpars(irest,3))**2 / restpars(irest,6)**2
      a4 = restpars(irest,7)**2
      w = a1 + a2 + a3 - a4
      a1 = dmin1(w,0.d0)
      f = f + a1*a1
    else if(ityperest(irest).eq.10) then
      w = restpars(irest,1)*xcart(icart,1) + &
          restpars(irest,2)*xcart(icart,2) + &
          restpars(irest,3)*xcart(icart,3) - &
          restpars(irest,4)
      a1 = dmin1(w,0.d0)
      f = f + scale * a1*a1
    else if(ityperest(irest).eq.11) then
      w = restpars(irest,1)*xcart(icart,1) + &
          restpars(irest,2)*xcart(icart,2) + &
          restpars(irest,3)*xcart(icart,3) - &
          restpars(irest,4)
      a1 = dmax1(w,0.d0)
      f = f + scale * a1*a1
    else if(ityperest(irest).eq.12) then
      a1 = xcart(icart,1) - restpars(irest,1)
      a2 = xcart(icart,2) - restpars(irest,2)
      a3 = xcart(icart,3) - restpars(irest,3)
      vnorm = sqrt(restpars(irest,4)**2 + restpars(irest,5)**2 + restpars(irest,6)**2)
      v1 = restpars(irest,4)/vnorm
      v2 = restpars(irest,5)/vnorm
      v3 = restpars(irest,6)/vnorm
      b1 = v1 * a1
      b2 = v2 * a2
      b3 = v3 * a3
      w = b1 + b2 + b3
      d = ( a1 - v1*w )**2 + ( a2 - v2*w )**2 + ( a3 - v3*w )**2
      f = f + scale2 * ( &
           dmax1(-w , 0.d0)**2 + &
           dmax1(w - restpars(irest,9), 0.d0)**2 + &
           dmax1(d - restpars(irest,7)**2 , 0.d0 )**2 )
    else if(ityperest(irest).eq.13) then
      a1 = xcart(icart,1) - restpars(irest,1)
      a2 = xcart(icart,2) - restpars(irest,2)
      a3 = xcart(icart,3) - restpars(irest,3)
      vnorm = sqrt(restpars(irest,4)**2 + restpars(irest,5)**2 + restpars(irest,6)**2)
      v1 = restpars(irest,4)/vnorm
      v2 = restpars(irest,5)/vnorm
      v3 = restpars(irest,6)/vnorm
      b1 = v1 * a1
      b2 = v2 * a2
      b3 = v3 * a3
      w = b1 + b2 + b3
      d = ( a1 - v1*w )**2 +( a2 - v2*w )**2 + ( a3 - v3*w )**2
      f = f + scale2 * ( &
           dmin1(-w , 0.d0)**2 * &
           dmin1(w - restpars(irest,9), 0.d0)**2 * &
           dmin1(d - restpars(irest,7)**2 , 0.d0 )**2 )
    end if 
  end do 
  return
end subroutine comprest

!
! Gradient evaluation
!

!
! Gradient relative to restraints
!

subroutine gwalls(icart,irest,gxcar)
      
  use sizes
  use molpa
  implicit none

  integer :: icart, irest
  double precision :: a1, a2, a3, a4, a5, a6, xmin, ymin, zmin, &
                      xmax, ymax, zmax, &
                      clength, b1, b2, b3, c1, c2, w, d, rg(3), &
                      vnorm, vv1, vv2, vv3, frab, frac, frbc, &
                      dfra(3), dfrb(3), dfrc(3), gxcar(maxatom,3)

  if(ityperest(irest).eq.2) then
    clength = restpars(irest,4)
    xmin = restpars(irest,1)
    ymin = restpars(irest,2)
    zmin = restpars(irest,3)
    xmax = restpars(irest,1) + clength
    ymax = restpars(irest,2) + clength
    zmax = restpars(irest,3) + clength
    a1 = xcart(icart,1) - xmin
    a2 = xcart(icart,2) - ymin
    a3 = xcart(icart,3) - zmin
    if(a1.lt.0.d0) gxcar(icart,1) = gxcar(icart,1) + scale * 2.d0 * a1
    if(a2.lt.0.d0) gxcar(icart,2) = gxcar(icart,2) + scale * 2.d0 * a2
    if(a3.lt.0.d0) gxcar(icart,3) = gxcar(icart,3) + scale * 2.d0 * a3
    a1 = xcart(icart,1) - xmax
    a2 = xcart(icart,2) - ymax
    a3 = xcart(icart,3) - zmax
    if(a1.gt.0.d0) gxcar(icart,1) = gxcar(icart,1) + scale * 2.d0 * a1
    if(a2.gt.0.d0) gxcar(icart,2) = gxcar(icart,2) + scale * 2.d0 * a2
    if(a3.gt.0.d0) gxcar(icart,3) = gxcar(icart,3) + scale * 2.d0 * a3
  else if(ityperest(irest).eq.3) then
    xmin = restpars(irest,1) 
    ymin = restpars(irest,2) 
    zmin = restpars(irest,3) 
    xmax = restpars(irest,4) 
    ymax = restpars(irest,5) 
    zmax = restpars(irest,6) 
    a1 = xcart(icart,1) - xmin
    a2 = xcart(icart,2) - ymin
    a3 = xcart(icart,3) - zmin
    if(a1.lt.0.d0) gxcar(icart,1) = gxcar(icart,1) + scale * 2.d0 * a1
    if(a2.lt.0.d0) gxcar(icart,2) = gxcar(icart,2) + scale * 2.d0 * a2
    if(a3.lt.0.d0) gxcar(icart,3) = gxcar(icart,3) + scale * 2.d0 * a3
    a1 = xcart(icart,1) - xmax
    a2 = xcart(icart,2) - ymax
    a3 = xcart(icart,3) - zmax
    if(a1.gt.0.d0) gxcar(icart,1) = gxcar(icart,1) + scale * 2.d0 * a1
    if(a2.gt.0.d0) gxcar(icart,2) = gxcar(icart,2) + scale * 2.d0 * a2
    if(a3.gt.0.d0) gxcar(icart,3) = gxcar(icart,3) + scale * 2.d0 * a3
  else if(ityperest(irest).eq.4) then
    d = (xcart(icart,1)-restpars(irest,1))**2 + &
        (xcart(icart,2)-restpars(irest,2))**2 + &
        (xcart(icart,3)-restpars(irest,3))**2 - &
        restpars(irest,4)**2
    if(d.gt.0.d0) then
      gxcar(icart,1) = gxcar(icart,1) + 4.d0 * scale2 * &
                       (xcart(icart,1)-restpars(irest,1))*d
      gxcar(icart,2) = gxcar(icart,2) + 4.d0 * scale2 * &
                       (xcart(icart,2)-restpars(irest,2))*d
      gxcar(icart,3) = gxcar(icart,3) + 4.d0 * scale2 * &
                       (xcart(icart,3)-restpars(irest,3))*d
    end if
  else if(ityperest(irest).eq.5) then
    a1 = xcart(icart,1)-restpars(irest,1)
    b1 = xcart(icart,2)-restpars(irest,2)
    c1 = xcart(icart,3)-restpars(irest,3)
    a2 = restpars(irest,4)**2
    b2 = restpars(irest,5)**2
    c2 = restpars(irest,6)**2
    d = a1**2/a2+b1**2/b2+c1**2/c2-restpars(irest,7)**2
    if(d.gt.0) then
      gxcar(icart,1) = gxcar(icart,1) + scale2*4.d0*d*a1/a2 
      gxcar(icart,2) = gxcar(icart,2) + scale2*4.d0*d*b1/b2
      gxcar(icart,3) = gxcar(icart,3) + scale2*4.d0*d*c1/c2
    end if
  else if(ityperest(irest).eq.6) then
    xmin = restpars(irest,1)
    ymin = restpars(irest,2)
    zmin = restpars(irest,3)
    xmax = restpars(irest,1) + restpars(irest,4)
    ymax = restpars(irest,2) + restpars(irest,4)
    zmax = restpars(irest,3) + restpars(irest,4)
    a1 = dmax1(xcart(icart,1) - xmin,0.d0)
    a2 = dmax1(xcart(icart,2) - ymin,0.d0)
    a3 = dmax1(xcart(icart,3) - zmin,0.d0)
    a4 = dmax1(xmax - xcart(icart,1),0.d0)
    a5 = dmax1(ymax - xcart(icart,2),0.d0)
    a6 = dmax1(zmax - xcart(icart,3),0.d0)
    w = a1*a2*a3*a4*a5*a6
    if(w.gt.0.d0) then
      gxcar(icart,1) = gxcar(icart,1) + a2*a3*a5*a6*(a4-a1)
      gxcar(icart,2) = gxcar(icart,2) + a1*a3*a4*a6*(a5-a2)
      gxcar(icart,3) = gxcar(icart,3) + a1*a2*a4*a5*(a6-a3)
    end if
  else if(ityperest(irest).eq.7) then
    xmin = restpars(irest,1)
    ymin = restpars(irest,2)
    zmin = restpars(irest,3)
    xmax = restpars(irest,4)
    ymax = restpars(irest,5)
    zmax = restpars(irest,6)
    a1 = dmax1(xcart(icart,1) - xmin,0.d0)
    a2 = dmax1(xcart(icart,2) - ymin,0.d0)
    a3 = dmax1(xcart(icart,3) - zmin,0.d0)
    a4 = dmax1(xmax - xcart(icart,1),0.d0)
    a5 = dmax1(ymax - xcart(icart,2),0.d0)
    a6 = dmax1(zmax - xcart(icart,3),0.d0)
    w = a1*a2*a3*a4*a5*a6
    if(w.gt.0.d0) then
      gxcar(icart,1) = gxcar(icart,1) + a2*a3*a5*a6*(a4-a1)
      gxcar(icart,2) = gxcar(icart,2) + a1*a3*a4*a6*(a5-a2)
      gxcar(icart,3) = gxcar(icart,3) + a1*a2*a4*a5*(a6-a3)
    end if
  else if(ityperest(irest).eq.8) then
    d = (xcart(icart,1)-restpars(irest,1))**2 + &
        (xcart(icart,2)-restpars(irest,2))**2 + &
        (xcart(icart,3)-restpars(irest,3))**2 - &
        restpars(irest,4)**2
    if(d.lt.0.d0) then
      gxcar(icart,1) = gxcar(icart,1) + 4.d0 * scale2 * &
                       (xcart(icart,1)-restpars(irest,1))*d
      gxcar(icart,2) = gxcar(icart,2) + 4.d0 * scale2 * &
                       (xcart(icart,2)-restpars(irest,2))*d
      gxcar(icart,3) = gxcar(icart,3) + 4.d0 * scale2 * &
                       (xcart(icart,3)-restpars(irest,3))*d
    end if
  else if(ityperest(irest).eq.9) then
    a1 = xcart(icart,1)-restpars(irest,1)
    b1 = xcart(icart,2)-restpars(irest,2)
    c1 = xcart(icart,3)-restpars(irest,3)
    a2 = restpars(irest,4)**2
    b2 = restpars(irest,5)**2
    c2 = restpars(irest,6)**2
    d = a1**2/a2+b1**2/b2+c1**2/c2-restpars(irest,7)**2
    if(d.lt.0) then
      d = scale2 * d
      gxcar(icart,1) = gxcar(icart,1) + 4.d0*d*a1/a2 
      gxcar(icart,2) = gxcar(icart,2) + 4.d0*d*b1/b2
      gxcar(icart,3) = gxcar(icart,3) + 4.d0*d*c1/c2
    end if
  else if(ityperest(irest).eq.10) then
    d = restpars(irest,1)*xcart(icart,1) + &
        restpars(irest,2)*xcart(icart,2) + &
        restpars(irest,3)*xcart(icart,3) - &
        restpars(irest,4)
    if(d.lt.0.d0) then
      d = scale * d
      gxcar(icart,1) = gxcar(icart,1) + 2.d0*restpars(irest,1)*d
      gxcar(icart,2) = gxcar(icart,2) + 2.d0*restpars(irest,2)*d
      gxcar(icart,3) = gxcar(icart,3) + 2.d0*restpars(irest,3)*d
    end if
  else if(ityperest(irest).eq.11) then
    d = restpars(irest,1)*xcart(icart,1) + &
        restpars(irest,2)*xcart(icart,2) + &
        restpars(irest,3)*xcart(icart,3) - &
        restpars(irest,4)
    if(d.gt.0.d0) then
      d = scale * d 
      gxcar(icart,1) = gxcar(icart,1) + 2.d0*restpars(irest,1)*d
      gxcar(icart,2) = gxcar(icart,2) + 2.d0*restpars(irest,2)*d
      gxcar(icart,3) = gxcar(icart,3) + 2.d0*restpars(irest,3)*d
    end if 
  else if(ityperest(irest).eq.12) then
    rg(1) = 0.0d0
    rg(2) = 0.0d0
    rg(3) = 0.0d0
    a1 = xcart(icart,1) - restpars(irest,1)
    a2 = xcart(icart,2) - restpars(irest,2)
    a3 = xcart(icart,3) - restpars(irest,3)
    vnorm = sqrt(restpars(irest,4)**2 + restpars(irest,5)**2 &
         + restpars(irest,6)**2)
    vv1 = restpars(irest,4)/vnorm
    vv2 = restpars(irest,5)/vnorm
    vv3 = restpars(irest,6)/vnorm
    b1 = vv1 * a1
    b2 = vv2 * a2
    b3 = vv3 * a3
    w = b1 + b2 + b3
    d = (a1 - vv1*w)**2 + (a2 - vv2*w)**2 + (a3 - vv3*w)**2
    rg(1) = scale2 * ( &
         -2*dmax1(-w , 0.d0) * vv1 + &
         2*dmax1(w - restpars(irest,9), 0.d0) * vv1 + &
         2*dmax1(d - restpars(irest,7)**2 , 0.d0) * &
         (2*(a1 - vv1*w)*(1 - vv1**2)+ &
          2*(a2 - vv2*w)*(-vv2*vv1)+ &
          2*(a3 - vv3*w)*(-vv3*vv1) ))
    rg(2) = scale2 * ( &
         -2*dmax1(-w , 0.d0) * vv2 + &
         2*dmax1(w - restpars(irest,9), 0.d0) * vv2 + &
         2*dmax1(d - restpars(irest,7)**2 , 0.d0) * &
         (2*(a1 - vv1*w)*(-vv1*vv2)+ &
          2*(a2 - vv2*w)*(1 - vv2**2)+ &
          2*(a3 - vv3*w)*(-vv3*vv2) ))
    rg(3) = scale2 * ( &
         -2*dmax1(-w , 0.d0) * vv3 + &
         2*dmax1(w - restpars(irest,9), 0.d0) * vv3 + &
         2*dmax1(d - restpars(irest,7)**2 , 0.d0) * &
         (2*(a1 - vv1*w)*(-vv1*vv3)+ &
          2*(a2 - vv2*w)*(-vv2*vv3)+ &
          2*(a3 - vv3*w)*(1 - vv3**2) ))
    gxcar(icart,1) = gxcar(icart,1) + rg(1)
    gxcar(icart,2) = gxcar(icart,2) + rg(2)
    gxcar(icart,3) = gxcar(icart,3) + rg(3)
  else if(ityperest(irest).eq.13) then
    rg(1) = 0.0d0
    rg(2) = 0.0d0
    rg(3) = 0.0d0
    a1 = xcart(icart,1) - restpars(irest,1)
    a2 = xcart(icart,2) - restpars(irest,2)
    a3 = xcart(icart,3) - restpars(irest,3)
    vnorm = sqrt(restpars(irest,4)**2 + restpars(irest,5)**2 &
         + restpars(irest,6)**2)
    vv1 = restpars(irest,4)/vnorm
    vv2 = restpars(irest,5)/vnorm
    vv3 = restpars(irest,6)/vnorm
    b1 = vv1 * a1
    b2 = vv2 * a2
    b3 = vv3 * a3
    w = b1 + b2 + b3
    d = (a1 - vv1*w)**2 + (a2 - vv2*w)**2 + (a3 - vv3*w)**2
    frab = dmin1(-w , 0.d0)**2 * dmin1(w - restpars(irest,9), 0.d0)**2
    frac = dmin1(-w , 0.d0)**2 * dmin1(d - restpars(irest,7)**2 , 0.d0 )**2 
    frbc = dmin1(w - restpars(irest,9), 0.d0)**2 * &
           dmin1(d - restpars(irest,7)**2 , 0.d0 )**2
    dfra(1) = -2*dmin1(-w , 0.d0) * vv1
    dfrb(1) = 2*dmin1(w - restpars(irest,9), 0.d0) * vv1
    dfrc(1) = 2*dmin1(d - restpars(irest,7)**2 , 0.d0) * &
         (2*(a1 - vv1*w)*(1 - vv1**2)+ &
          2*(a2 - vv2*w)*(-vv2*vv1)+  &
          2*(a3 - vv3*w)*(-vv3*vv1) )
    dfra(2) = -2*dmin1(-w , 0.d0) * vv2
    dfrb(2) = 2*dmin1(w - restpars(irest,9), 0.d0) * vv2
    dfrc(2) = 2*dmin1(d - restpars(irest,7)**2 , 0.d0) * &
         (2*(a1 - vv1*w)*(-vv1*vv2)+  &
          2*(a2 - vv2*w)*(1 - vv2**2)+ & 
          2*(a3 - vv3*w)*(-vv3*vv2) )
    dfra(3) = -2*dmin1(-w , 0.d0) * vv3
    dfrb(3) = 2*dmin1(w - restpars(irest,9), 0.d0) * vv3
    dfrc(3) = 2*dmin1(d - restpars(irest,7)**2 , 0.d0) * &
         (2*(a1 - vv1*w)*(-vv1*vv3)+ &
          2*(a2 - vv2*w)*(-vv2*vv3)+ &
          2*(a3 - vv3*w)*(1 - vv3**2) )
    rg(1) = scale2 * ( dfra(1)*frbc + dfrb(1)*frac + dfrc(1)*frab)
    rg(2) = scale2 * ( dfra(2)*frbc + dfrb(2)*frac + dfrc(2)*frab)
    rg(3) = scale2 * ( dfra(3)*frbc + dfrb(3)*frac + dfrc(3)*frab)
    gxcar(icart,1) = gxcar(icart,1) + rg(1)
    gxcar(icart,2) = gxcar(icart,2) + rg(2)
    gxcar(icart,3) = gxcar(icart,3) + rg(3)
  end if
      
  return 
end subroutine gwalls

!
! Subroutine that performs finite difference and analytical gradient
! comparision. Used only for test purpouses
!

subroutine compgrad(n,x)

  use sizes
  implicit none

  integer :: n, i, iworst
  double precision :: x(*), g(nn), fx, step, gcomp, gbest, eworst, &
                      error, steperror, stepbest
  real :: time0, tarray(2), etime

  write(*,*)
  write(*,*) ' Comparing analytical and finite-difference '
  write(*,*) ' gradients... may take a while. '
  write(*,*)
  write(*,*) ' Five first center of masses and angles of tested point: '
  do i = 1, 15, 3
    write(*,"( i4,6(tr2,f8.3) )") (i+2)/3, x(i), x(i+1), x(i+2), x(n/2+i),&
                                  x(n/2+i+1),x(n/2+i+2)
  end do
  write(*,*) 
  write(*,*) ' Computing gradient ... ' 

  call feasy(x, fx) 
  write(*,*) ' Function value on test point: ', fx
  open(98, file = 'chkgrad.log',status='unknown') 
  write(98, *)'Function Value = ', fx 
  call geasy(n, x, g) 

  write(98,"( t2,'Component',t16,'Analytical',t33,'Discrete', &
              t51,'Error',t62,'Best step' )")
  time0 = etime(tarray)
  eworst = 0.d0
  do i = 1, n
    if(etime(tarray)-time0.gt.10.) then
      time0 = etime(tarray)
      write(*,*) ' Computing the ',i,'th of ',n,' components. '
    end if
    error = 1.d20
    step = 1.d-2
    do while(error.gt.1.d-6.and.step.ge.1.d-20)
      call discret(i,x,gcomp,step)
      if(dmin1(abs(g(i)),abs(gcomp)).gt.1.d-10) then
        steperror = abs( ( gcomp - g(i) ) / g(i) )
      else
        steperror = abs( gcomp - g(i) )
      end if
      if( steperror .lt. error ) then
        error = steperror
        gbest = gcomp
        stepbest = step
      end if
      step = step / 10.d0
    end do
    write(98,"(i10,4(tr2,d13.7))") i, g(i), gbest, error, stepbest 
    if(error.gt.eworst) then
      iworst = i
      eworst = error
    end if
  end do
  write(98,*) 'Maximum difference = ', iworst,' Error= ', eworst
  write(*,*) ' Done. '
  stop

end subroutine compgrad

subroutine discret(icomp,x,gcomp,step)

  implicit none
  integer :: icomp
  double precision :: save, step, x(*), fplus, fminus, gcomp

  save = x(icomp) 
  x(icomp) = save + step 
  call feasy(x, fplus) 
  x(icomp) = save - step 
  call feasy(x, fminus) 
  gcomp = (fplus - fminus) / (2.d0 * step) 
  x(icomp) = save 

  return      
end subroutine discret

!
! Subroutine that computes the function value
!
!
!
  
subroutine feasyseq(x,f) 
      
  use sizes
  use molpa
  implicit none

  double precision :: v1(3), v2(3), v3(3) 
  double precision :: x(nn)
  double precision :: f,fparcs,fplus
  double precision :: xtemp, ytemp, ztemp
  double precision :: xbar, ybar, zbar
  double precision :: beta, gama, teta
  double precision :: flast

  integer :: i, j, k
  integer :: ilugan, ilubar, icart, itype, imol, iatom, idatom, &
             iboxx, iboxy, iboxz 

  ! Reset function value

  f = 0.d0 
  frest = 0.d0
  fdist = 0.d0

  ! Reset boxes

  if(.not.init1) call resetboxes(nboxes,latomfirst,latomfix)

  ! Transform baricenter and angles into cartesian coordinates 
  ! Computes cartesian coordinates from vector x and coor 

  ilubar = 0 
  ilugan = ntotmol*3 
  icart = natfix

  do itype = 1, ntype 
    if(.not.comptype(itype)) then
      icart = icart + nmols(itype)*natoms(itype)
    else
    do imol = 1, nmols(itype) 

      xbar = x(ilubar+1) 
      ybar = x(ilubar+2) 
      zbar = x(ilubar+3) 
  
      ! Computing the rotation matrix

      beta = x(ilugan+1)
      gama = x(ilugan+2)
      teta = x(ilugan+3)

      call eulerrmat(beta,gama,teta,v1,v2,v3)  

      ! Looping over the atoms of this molecule
  
      idatom = idfirst(itype) - 1
      do iatom = 1, natoms(itype) 

        icart = icart + 1
        idatom = idatom + 1

        ! Computing the cartesian coordinates for this atom

        call compcart(icart,xcart, &
                      xbar,ybar,zbar, &
                      coor(idatom,1),coor(idatom,2),coor(idatom,3), &
                      v1,v2,v3)

        ! Adding to f the value relative to constraints for this atom

        call comprest(xcart,restpars, &
                      scale,scale2, &
                      nratom,ityperest,iratom, &
                      icart,fplus)
        f = f + fplus
        frest = dmax1(frest,fplus)
        if(move) fatom(icart) = fatom(icart) + fplus

        ! Putting atoms in their boxes

        if(.not.init1) then

          xtemp = xcart(icart,1) - sizemin(1) 
          ytemp = xcart(icart,2) - sizemin(2) 
          ztemp = xcart(icart,3) - sizemin(3) 

          iboxx = int(xtemp/boxl(1)) + 1
          iboxy = int(ytemp/boxl(2)) + 1
          iboxz = int(ztemp/boxl(3)) + 1

          if(xtemp.le.0) iboxx = 1
          if(ytemp.le.0) iboxy = 1
          if(ztemp.le.0) iboxz = 1 
          if(iboxx.gt.nboxes(1)) iboxx = nboxes(1)
          if(iboxy.gt.nboxes(2)) iboxy = nboxes(2)
          if(iboxz.gt.nboxes(3)) iboxz = nboxes(3)

          latomnext(icart) = latomfirst(iboxx,iboxy,iboxz)
          latomfirst(iboxx,iboxy,iboxz) = icart

          ibtype(icart) = itype
          ibmol(icart) = imol

        end if

      end do 
 
      ilugan = ilugan + 3 
      ilubar = ilubar + 3 

    end do
    end if
  end do            

  if(init1) return

  ! Minimum distance function evaluation

  do i = 1, nboxes(1)
    do j = 1, nboxes(2)
      do k = 1, nboxes(3)

        icart = latomfirst(i,j,k)
        do while ( icart .ne. 0 ) 

          if(comptype(ibtype(icart))) then

            ! Vector that keeps the value for this atom

            if(move) flast = f 

            ! Interactions inside box

            f = f + fparcs(icart,latomnext(icart))

            ! Interactions of boxes that share faces

            f = f + fparcs(icart,latomfirst(i+1,j,k))
            f = f + fparcs(icart,latomfirst(i,j+1,k))
            f = f + fparcs(icart,latomfirst(i,j,k+1))

            ! Interactions of boxes that share axes

            f = f + fparcs(icart,latomfirst(i+1,j+1,k))
            f = f + fparcs(icart,latomfirst(i+1,j,k+1))
            f = f + fparcs(icart,latomfirst(i+1,j-1,k))
            f = f + fparcs(icart,latomfirst(i+1,j,k-1))
            f = f + fparcs(icart,latomfirst(i,j+1,k+1))
            f = f + fparcs(icart,latomfirst(i,j+1,k-1))

            ! Interactions of boxes that share vertices

            f = f + fparcs(icart,latomfirst(i+1,j+1,k+1))
            f = f + fparcs(icart,latomfirst(i+1,j+1,k-1))
            f = f + fparcs(icart,latomfirst(i+1,j-1,k+1))
            f = f + fparcs(icart,latomfirst(i+1,j-1,k-1))

            ! If going to move bad molecules, update fatom

            if(move) fatom(icart) = fatom(icart) + f - flast

          end if

          icart = latomnext(icart)
        end do
  
      end do
    end do
  end do

  return
end subroutine feasyseq 

!
! Function that computes the main part of the objective function
!

double precision function fparcs(icart,firstjcart)

  use sizes
  use molpa

  ! SCALAR ARGUMENTS
  integer :: icart,firstjcart

  ! LOCAL SCALARS
  integer :: jcart
  double precision :: a1,a2,a3,datom,tol

  fparcs = 0.0d0
  jcart = firstjcart
  do while ( jcart .ne. 0 )
    if(comptype(ibtype(jcart))) then
      if(ibmol(icart).ne.ibmol(jcart).or. &
         ibtype(icart).ne.ibtype(jcart)) then
        a1 = xcart(icart, 1)-xcart(jcart, 1) 
        a2 = xcart(icart, 2)-xcart(jcart, 2) 
        a3 = xcart(icart, 3)-xcart(jcart, 3) 
        datom = a1 * a1 + a2 * a2 + a3 * a3
        tol = (radius(icart)+radius(jcart))**2
        a1 = dmin1(datom - tol, 0.d0)
        fparcs = fparcs + a1 * a1
        tol = (radius_ini(icart)+radius_ini(jcart))**2
        fdist = dmax1(tol-datom,fdist)
      end if
    end if
    jcart = latomnext(jcart)
  end do

end function fparcs

