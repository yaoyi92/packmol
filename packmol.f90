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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
!
! Packmol: A package for building initial configurations for
! molecular dynamics simulations, to be published, 2008.
!
! http://www.ime.unicamp.br/~martinez/packmol
!
! Usage (see the page above for further information):
!
! ./packmol < inputfile.inp
!
! References:
!
! L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez,
! PACKMOL: A package for building initial configurations for
! molecular dynamics simulations, J. Comp. Chem. 30:2157-2164, 2009.
!
! J. M. Martinez and L. Martinez, 
! Packing optimization for the automated generation of complex
! system's initial configurations for molcular dynamics and
! docking. J. Comp. Chem. 24:819-825, 2003.
!
! This version of Packmol uses the optimization method GENCAN which
! is a part of the TANGO (Trustable Algorithms for Nonlinear General
! Optimization) project.
! Reference:
! E. G. Birgin, J. M. Martinez, Comp. Opt. Appl. 23:101-125, 2002.
! http://www.ime.usp.br/~egbirgin/tango
!
!

program packmol

  use sizes
  use molpa
  use usegencan
  implicit none

  integer :: irestline(maxrest)
  integer :: linestrut(maxtype,2)
  integer :: itype, nrest, irest, idatom, iatom
  integer :: ntemp, ntfix, idtemp, nmtemp, natemp, nlines
  integer :: linesttmp1, linesttmp2, jtype
  integer :: ntmol, n, iftype, icart, imol, iicart, iline_atoms
  integer :: i, iline, iiatom, iat, iirest, iratcount, ival
  integer :: seed
  integer :: nloop, loop
  integer :: maxcon(maxatom)
  integer :: ntcon(9), nconnect(maxatom,8) 
  integer :: ntottemp, ilubar, ilugan  
  integer :: resnumbers(maxtype), resntemp
  integer :: charl, writeout, ioerr
  integer :: input_itype(maxtype)
      
  double precision :: v1(3),v2(3),v3(3)
  double precision :: x(nn), xfull(nn), xbest(nn)
  double precision :: rad, radiuswork(maxatom), radscale
  double precision :: cmx, cmy, cmz, beta, gama, teta
  double precision :: xtemp, ytemp, ztemp
  double precision :: fx, bestf, flast, fout
  double precision :: fimp, fimprov, precision
  double precision :: amass(maxatom), charge(maxatom)
  double precision :: discale, movefrac, dism
  double precision :: add_sides_fix
  double precision :: sidemax
  double precision, parameter :: pi=3.141592653589793d0

  real :: etime, tarray(2), time0
  
  character(len=200) :: keyword(maxlines,maxkeywords)
  character(len=200) :: record
  character(len=200) :: name(maxtype)
  character(len=80) :: pdbfile(maxtype), xyzfile
  character(len=3) :: ele(maxatom)
  character(len=200) :: xyzout        
  character(len=80) :: dash1_line

  logical :: fix,fixed(maxtype),fixtmp,randini,check,chkgrad
  logical :: pdb,tinker,xyz,moldy,rests,writebad,writexyz
  logical :: add_amber_ter, add_box_sides
  logical :: movebadprint, hasbad, movebadrandom
  logical :: thisisfixed(maxtype), changechains(maxtype)

  ! Start time computation

  time0 = etime(tarray)

  ! Printing title

  dash1_line = "( 62('#') )"
  write(*,"( 62('#'), /,/&
             ' PACKMOL - Packing optimization for the automated', /&
             ' generation of starting configurations for',        /&
             ' molecular dynamics. ',/&
             ' ',/&
             t42,' Version 16.070 ',/&
             ,/,62('#'),               /,/)")
      
  ! Reading input file

  call getinp(dism,precision,sidemax,&
              ntype,nlines,nrest,&
              natoms,idfirst,nconnect,maxcon,nmols,&
              seed,&
              discale,nloop,&
              irestline,ityperest,linestrut,&
              coor,amass,charge,restpars,&
              pdbfile,name,ele,keyword,&
              xyzout,writeout,writebad,&
              tinker,pdb,xyz,moldy,check,chkgrad,&
              randini,resnumbers,movefrac,movebadrandom,changechains,&
              add_amber_ter,add_box_sides,add_sides_fix)

  ! Put molecules in their center of mass

  call cenmass(coor,amass,ntype,nlines,idfirst,natoms,&
               keyword,linestrut)
 
  ! Computing the total number of atoms
     
  ntotat = 0
  do itype = 1, ntype
    ntotat = ntotat + natoms(itype) * nmols(itype)
  end do              
  write(*,*) ' Total number of atoms: ', ntotat

  if(ntotat.gt.maxatom) then
    write(*,*)' ERROR: Total number of atoms greater than maxatom.'
    write(*,*)'        Change the maxatom (sizes.f90 file) '
    stop
  end if

  ! Setting the vector that contains the default tolerance

  do i = 1, ntotat
    radius(i) = dism/2.d0
  end do

  ! Setting the radius defined for atoms of each molecule, 
  ! but not atom-specific, first

  icart = 0
  do itype = 1, ntype
    iline = linestrut(itype,1)
    iline_atoms = 0 
    do while( iline <= linestrut(itype,2) )
      if ( keyword(iline,1) == "atoms" ) then
        iline_atoms = iline
        iline = iline + 1
        cycle
      end if
      if ( keyword(iline,1) == "end" .and. &    
           keyword(iline,2) == "atoms" ) then
        iline_atoms = 0  
        iline = iline + 1
        cycle
      end if
      if ( iline_atoms == 0 ) then
        if ( keyword(iline,1) == "radius" ) then
          read(keyword(iline,2),*,iostat=ioerr) rad
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Could not read radius from keyword. '
            stop
          end if
          iicart = icart
          do imol = 1, nmols(itype)
            do iatom = 1, natoms(itype)
              iicart = iicart + 1
              radius(iicart) = rad 
            end do
          end do
        end if
      end if
      iline = iline + 1
    end do
    icart = icart + nmols(itype)*natoms(itype)
  end do
 
  ! If some radius was defined using atom-specific definitions, overwrite
  ! the general radius defined for the molecule

  icart = 0
  do itype = 1, ntype
    iline = linestrut(itype,1)
    iline_atoms = 0 
    do while( iline <= linestrut(itype,2) )
      if ( keyword(iline,1) == "atoms" ) then
        iline_atoms = iline
        iline = iline + 1
        cycle
      end if
      if ( keyword(iline,1) == "end" .and. &    
           keyword(iline,2) == "atoms" ) then
        iline_atoms = 0  
        iline = iline + 1
        cycle
      end if
      if ( iline_atoms /= 0 ) then
        if ( keyword(iline,1) == "radius" ) then
          read(keyword(iline,2),*,iostat=ioerr) rad
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Could not read radius from keyword. '
            stop
          end if
          ival = 2
          do
            read(keyword(iline_atoms,ival),*,iostat=ioerr) iat
            if ( ioerr /= 0 ) exit
            if ( iat > natoms(itype) ) then
              write(*,*) ' ERROR: atom selection with index greater than number of '
              write(*,*) '        atoms in structure ', itype
              stop
            end if
            radius(icart+iat) = rad
            ival = ival + 1
          end do
        end if
      end if
      iline = iline + 1
    end do
    iicart = icart
    icart = icart + natoms(itype)
    do imol = 2, nmols(itype)
      do iatom = 1, natoms(itype)
        icart = icart + 1
        radius(icart) = radius(iicart+iatom)
      end do
    end do
  end do

  ! Put fixed molecules in the specified position

  do itype = 1, ntype
    fixed(itype) = .false.
  end do

  do irest = 1, nrest
    if(ityperest(irest).eq.1) then
      do itype = 1, ntype
        if(irestline(irest).gt.linestrut(itype,1).and.&
           irestline(irest).lt.linestrut(itype,2)) then
          cmx = restpars(irest,1) 
          cmy = restpars(irest,2)
          cmz = restpars(irest,3)    
          beta = restpars(irest,4) 
          gama = restpars(irest,5) 
          teta = restpars(irest,6) 

          ! Compute rotation matrix from euler angles

          call eulerfixed(beta,gama,teta,v1,v2,v3)                 

          idatom = idfirst(itype) - 1
          do iatom = 1, natoms(itype)
            idatom = idatom + 1
            xtemp =   coor(idatom,1)*v1(1) &
                    + coor(idatom,2)*v2(1) &
                    + coor(idatom,3)*v3(1) 
            ytemp =   coor(idatom,1)*v1(2) &
                    + coor(idatom,2)*v2(2) &
                    + coor(idatom,3)*v3(2) 
            ztemp =   coor(idatom,1)*v1(3) &
                    + coor(idatom,2)*v2(3) &
                    + coor(idatom,3)*v3(3) 
            coor(idatom, 1) = xtemp + cmx
            coor(idatom, 2) = ytemp + cmy
            coor(idatom, 3) = ztemp + cmz 
          end do
          record = name(itype)
          write(*,*) ' Molecule ',record(1:charl(record)),&
                     '(',itype,') will be fixed.' 
          fixed(itype) = .true.
          if(nmols(itype).gt.1) then
            write(*,*)' ERROR: You cannot set number > 1',&
                      ' for fixed molecules. '
            stop
          end if
        end if
      end do
    end if
  end do 

  ! Reseting parameters for removing the fixed molecules

  fix = .false.
  ntemp = 0
  do itype = 1, ntype

  ! input_itype and thisisfixed vectors are used only to preserve the
  ! order of input in the output files

    input_itype(itype) = itype
    if(fixed(itype)) then
      fix = .true.
      thisisfixed(itype) = .true.
    else
      ntemp = ntemp + 1
      thisisfixed(itype) = .false.
    end if
  end do
  ntfix = ntype
  ntype = ntemp     

  do i = 1, ntfix - ntype 
    do itype = 1, ntfix - 1
      if(fixed(itype)) then
        record = name(itype)
        fixtmp = fixed(itype)
        idtemp = idfirst(itype)
        nmtemp = nmols(itype)
        natemp = natoms(itype)
        resntemp = resnumbers(itype)
        if(pdb) xyzfile = pdbfile(itype)
        linesttmp1 = linestrut(itype,1)
        linesttmp2 = linestrut(itype,2)
        jtype = itype + 1
        if(.not.fixed(jtype)) then
          name(itype) = name(jtype)
          name(jtype) = record(1:10)
          idfirst(itype) = idfirst(jtype)
          idfirst(jtype) = idtemp
          fixed(itype) = fixed(jtype)
          fixed(jtype) = fixtmp
          nmols(itype) = nmols(jtype)
          nmols(jtype) = nmtemp
          natoms(itype) = natoms(jtype)
          natoms(jtype) = natemp
          resnumbers(itype) = resnumbers(jtype)
          resnumbers(jtype) = resntemp
          if(pdb) then
            pdbfile(itype) = pdbfile(jtype) 
            pdbfile(jtype) = xyzfile
          end if
          linestrut(itype,1) = linestrut(jtype,1)
          linestrut(itype,2) = linestrut(jtype,2)
          linestrut(jtype,1) = linesttmp1
          linestrut(jtype,2) = linesttmp2
        end if
      end if
    end do
  end do

  ! Computing the number of variables
  !
  ! ntype: 1...ntype (counter for the number of free structures)
  !
  ! ntfix: 1...ntype...ntfix (counter for the total number of structures)
  !

  ntmol = 0
  do itype = 1, ntfix
    ntmol = ntmol + nmols(itype)
  end do
  ntotmol = 0 
  do itype = 1, ntype 
    ntotmol = ntotmol + nmols(itype)       
  end do     
  n = ntotmol * 6
  write(*,*) ' Total number of molecules: ', ntmol
  write(*,*) ' Number of fixed molecules: ', ntmol - ntotmol
  write(*,*) ' Number of free molecules: ', ntotmol
  write(*,*) ' Number of variables: ', n 

  ! Computing the total number of fixed atoms

  natfix = 0
  if(fix) then
    do iftype = ntype + 1, ntfix
      natfix = natfix + natoms(iftype)
    end do
  end if       
  write(*,*) ' Total number of fixed atoms: ', natfix

  ! Setting the array that contains the restrictions per atom

  icart = natfix
  do itype = 1, ntype
    rests = .false.
    do imol = 1, nmols(itype)
      idatom = idfirst(itype) - 1      
      do iatom = 1, natoms(itype) 
        icart = icart + 1
        idatom = idatom + 1
        nratom(icart) = 0
        iratcount = 0
        do i = 1, mrperatom
          iratom(icart,i) = 0
        end do
        iline = linestrut(itype,1)
        do while(iline.lt.linestrut(itype,2))
          iline = iline + 1
          if(keyword(iline,1).eq.'atoms') then
            iiatom = -1
            do iat = 2, maxkeywords
              read(keyword(iline,iat),*,iostat=ioerr) iiatom
              if ( ioerr /= 0 ) exit
              if(iatom.eq.iiatom) exit
            end do
            do while(keyword(iline,1).ne.'end'.and.&
                     keyword(iline,2).ne.'atoms')
              iline = iline + 1
              if(iatom.eq.iiatom) then
                if(keyword(iline,1).eq.'inside'.or.&
                   keyword(iline,1).eq.'outside'.or.&
                   keyword(iline,1).eq.'over'.or.&
                   keyword(iline,1).eq.'below') then
                  nratom(icart) = nratom(icart) + 1
                  iratcount = iratcount + 1
                  do irest = 1, nrest
                    if(irestline(irest).eq.iline) iirest = irest
                  end do
                  iratom(icart,iratcount) = iirest
                end if
              end if
            end do
            iline = iline - 1
          else if(keyword(iline,1).eq.'inside'.or.&
                  keyword(iline,1).eq.'outside'.or.&
                  keyword(iline,1).eq.'over'.or.&
                  keyword(iline,1).eq.'below') then
            nratom(icart) = nratom(icart) + 1    
            iratcount = iratcount + 1
            do irest = 1, nrest
              if(irestline(irest).eq.iline) iirest = irest
            end do
            iratom(icart,iratcount) = iirest
          end if
        end do
        if(nratom(icart).gt.0) rests = .true.
      end do 
      if(.not.rests) then
        write(*,*) ' ERROR: Some molecule has no geometrical',&
                   ' restriction defined: nothing to do.'
        stop
      end if
    end do
  end do

  ! Read the constraints to rotations about axis, if set

  do itype = 1, ntype
    constrain_rot(itype,1) = .false.
    constrain_rot(itype,2) = .false.
    constrain_rot(itype,3) = .false.
    iline = linestrut(itype,1)
    do while(iline.lt.linestrut(itype,2))
      iline = iline + 1
      if(keyword(iline,1).eq.'constrain_rotation') then
        if(iline.gt.linestrut(itype,1).and.&
           iline.lt.linestrut(itype,2)) then

           ! Note that for movable molecules, teta is a rotation on the x-axis,
           !                                  gama is a rotation on the z-axis,
           !                                  beta is a rotation on the y-axis
           !                                  (see eulerrmat routine)

          if(keyword(iline,2).eq.'x') then
            constrain_rot(itype,3) = .true.
            read(keyword(iline,3),*) rot_bound(itype,3,1)
            read(keyword(iline,4),*) rot_bound(itype,3,2)
            rot_bound(itype,3,1) = rot_bound(itype,3,1)*pi/180.d0
            rot_bound(itype,3,2) = rot_bound(itype,3,2)*pi/180.d0
  
            write(*,*) ' Rotations about x axis of molecules of ',&
                       ' type ', itype, ' will be constrained. '
          end if
          if(keyword(iline,2).eq.'y') then
            constrain_rot(itype,1) = .true.
            read(keyword(iline,3),*) rot_bound(itype,1,1)
            read(keyword(iline,4),*) rot_bound(itype,1,2)
            rot_bound(itype,1,1) = rot_bound(itype,1,1)*pi/180.d0
            rot_bound(itype,1,2) = rot_bound(itype,1,2)*pi/180.d0

            write(*,*) ' Rotations about y axis of molecules of ',&
                       ' type ', itype, ' will be constrained. '
          end if
          if(keyword(iline,2).eq.'z') then
            constrain_rot(itype,2) = .true.
            read(keyword(iline,3),*) rot_bound(itype,2,1)
            read(keyword(iline,4),*) rot_bound(itype,2,2)
            rot_bound(itype,2,1) = rot_bound(itype,2,1)*pi/180.d0
            rot_bound(itype,2,2) = rot_bound(itype,2,2)*pi/180.d0

            write(*,*) ' Rotations about z axis of molecules of ',&
                       ' type ', itype, ' will be constrained. '
          end if
          if ( keyword(iline,2) /= 'x' .and. &
               keyword(iline,2) /= 'y' .and. &
               keyword(iline,2) /= 'z' ) then
            write(*,*) ' ERROR: constrain_rotation option not properly defined (not x, y, or z) '
            stop
          end if
        end if
      end if
    end do
  end do
 
  ! If there are no variables (only fixed molecules, stop)

  if(n.eq.0) then
    call output(x,amass,&
                irestline,linestrut,maxcon,ntcon,nconnect,&
                nrest,ntfix,resnumbers,&
                ele,pdbfile,xyzout,name,&
                pdb,tinker,xyz,moldy,fix,&
                add_amber_ter,add_box_sides,add_sides_fix,&
                input_itype,thisisfixed,changechains)
    write(*,dash1_line)
    write(*,*) ' There are only fixed molecules, therefore '
    write(*,*) ' there is nothing to do. '
    write(*,*) ' The output file contains the fixed molecule '
    write(*,*) ' in the desired position. '
    write(*,dash1_line)
    stop
  end if
  
  !
  ! (Re)setting parameters and building initial point
  !

  call initial(seed,randini,x,n,ntfix,fix,moldy,&
               chkgrad,nloop,discale,precision,sidemax,&
               movefrac,movebadrandom,check)

  ! Computing the energy at the initial point

  radscale = 1.d0
  do i = 1, ntotat
    radius_ini(i) = radius(i)
  end do
  call feasy(x,fx)
  write(*,*) ' Objective function at initial point: ', fx
  bestf = fx
  flast = fx
  fout = fx
  do i = 1, n
    xbest(i) = x(i)
  end do

  if(check) then
    call output(x,amass,&
                irestline,linestrut,maxcon,ntcon,nconnect,&
                nrest,ntfix,resnumbers,&
                ele,pdbfile,xyzout,name,&
                pdb,tinker,xyz,moldy,fix,&
                add_amber_ter,add_box_sides,add_sides_fix,&
                input_itype,thisisfixed,changechains)
    write(*,*) ' Wrote initial point to output file: ',&
               xyzout(1:charl(xyzout)) 
    stop
  end if

  !
  ! Main loop: first pack types of molecules separately, then
  ! pack all molecules together
  !

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

  main : do while(itype.le.ntype+1)
    itype = itype + 1
 
    ! Use larger tolerance than required to improove separation

    radscale = discale
    do i = 1, ntotat
      radius(i) = discale*radius_ini(i)
    end do
    
    if(itype.le.ntype) then
      if(nmols(itype).eq.1) itype = itype + 1
    end if

    ! Adjusting parameters for packing only this type

    if(itype.le.ntype) then
      write(*,*)
      write(*,dash1_line)
      write(*,*)
      write(*,*) ' Packing molecules of type ', itype
      write(*,*)
      write(*,dash1_line)
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
    end if

    ! If itype=ntype+1 restore original vectors and pack all molecules

    if(itype.eq.ntype+1) then
      n = ntemp 
      ntotmol = ntottemp
      do i = 1, n
        x(i) = xfull(i)
      end do
      do itype = 1, ntype
        comptype(itype) = .true.
      end do
      if(ntype.gt.1) then
        write(*,*)
        write(*,dash1_line)
        write(*,*)
        write(*,*)' Solving the problem for all molecules together.'
        write(*,*)
        write(*,dash1_line)
        write(*,*)
      end if
    end if
 
    loop = -1
    gencanloop : do while(loop.le.nloop)
      loop = loop + 1

      ! Reseting the parameters relative to the improvement of the function
         
      if(loop.eq.0) then
        fimp = 1.d99
        fimprov = fimp
        do i = 1, ntotat
          radiuswork(i) = radius(i) 
          radius(i) = radius_ini(i)
        end do
        call feasy(x,fx)
        do i = 1, ntotat
          radius(i) = radiuswork(i)
        end do
        bestf = fx
        flast = fx
      end if

      ! Moving bad molecules

      if(radscale == 1.d0 .and. fimp.le.10.d0) then
        movebadprint = .true.
        call movebad(n,x,fx,movefrac,movebadrandom,precision,seed,hasbad,movebadprint)
        flast = fx
      end if

      if(loop.eq.nloop.and.itype.eq.ntype+1) then
        write(*,*)' STOP: Maximum number of GENCAN loops achieved.'
        call checkpoint(n,xbest,amass,&
                        nrest,ntfix,nloop,&
                        irestline,linestrut,maxcon,ntcon,nconnect,&
                        ele,pdbfile,xyzout,name,&
                        pdb,tinker,xyz,moldy,fix,&
                        movefrac,movebadrandom,precision,seed,resnumbers,&
                        add_amber_ter,add_box_sides,add_sides_fix,&
                        input_itype,thisisfixed,changechains)
        stop
      end if

      write(*,"( /, 17('-'),' Starting GENCAN loop(',i4,') ',17('-'),/&
      '          Scaling radii by:',f10.2 )") loop, radscale

      ! CALL GENCAN

      write(*,"( '  Packing:|0 ',tr39,'  10|' )")
      call pgencan(n,x,fx)

      ! Compute the statistics of the last optimization loop

      do i = 1, ntotat
        radiuswork(i) = radius(i)
        radius(i) = radius_ini(i)
      end do
      call feasy(x,fx)
      do i = 1, ntotat
        radius(i) = radiuswork(i)
      end do

      if(bestf.gt.0.d0) fimprov = -100.d0 * (fx - bestf) / bestf
      if(bestf.eq.0.d0) fimprov = 100.d0
      if(flast.gt.0.d0) fimp = -100.d0 * (fx - flast) / flast
      if(flast.eq.0.d0) fimp = 100.d0
      fimp = dmin1(99.99d0,dmax1(-99.99d0,fimp))
      fimprov = dmin1(99.99d0,dmax1(-99.99d0,fimprov))

      write(*,"( /&
             '  Function value from last GENCAN loop: f = ', e10.5, /&
             '  Best function value before: f = ', e10.5,           /&
             '  Improvement from best function value: ', f8.2, ' %',/&
             '  Improvement from last loop: ', f8.2, ' %',          /&
             '  Maximum violation of target distance: ', f12.6,/&
             '  Maximum violation of the constraints: ', e10.5,     /&
             62('-'),/ )")  fx, bestf, fimprov, fimp, fdist, frest
      flast = fx

      ! If the distance between molecules is satisfactory, restore the radii

      if ( radscale > 1.d0 ) then
        if( ( fdist < precision .and. fimp < 10.d0 ) .or. &
            fimp < 2.d0 ) then
          radscale = dmax1(0.9*radscale,1.d0)
          do i = 1, ntotat
            radius(i) = dmax1(radius_ini(i),0.9d0*radius(i))
          end do
        end if
      end if

      ! Updating best point

      if(fx.le.bestf) then
        bestf = fx 
        if(itype.eq.ntype+1) then
          do i = 1, n
            xbest(i) = x(i)
          end do
        end if
      end if

      ! Writing output file 

      writexyz = .false.
      if ( itype == ntype + 1 ) then

        ! If solution was found
        if ( ( fdist < precision .and. frest < precision ) .or. &
               bestf < precision ) then
          writexyz = .true.
          write(*,*) ' Solution written to file: ', xyzout(1:charl(xyzout))

        ! If this is the best structure so far
        else if( mod(loop+1,writeout) == 0 .and. bestf < fout ) then
          writexyz = .true.
          write(*,*) ' Best solution written to file: ', xyzout(1:charl(xyzout))
          fout = bestf

        ! If the user required printing even bad structures
        else if ( mod(loop+1,writeout) == 0 .and. writebad ) then
          writexyz = .true.
          write(*,*) ' Writing current (perhaps bad) structure to file: ', xyzout(1:charl(xyzout))
        end if

        if ( writexyz ) then
          call output(x,amass,&
                      irestline,linestrut,maxcon,ntcon,nconnect,&
                      nrest,ntfix,resnumbers,&
                      ele,pdbfile,xyzout,name,&
                      pdb,tinker,xyz,moldy,fix,&
                      add_amber_ter,add_box_sides,add_sides_fix,&
                      input_itype,thisisfixed,changechains)
        end if

      end if

      ! When the solution is found, print success information and stop

      if((fdist.lt.precision.and.&
          frest.lt.precision).or.&
          bestf.lt.precision) then

        if(itype.le.ntype) then
          write(*,dash1_line)
          write(*,*)' Packing solved for molecules of type', itype
          write(*,*)' Objective function value: ', bestf
          write(*,*)' Maximum violation of target distance: ',fdist
          write(*,*)' Max. constraint violation: ', frest
          write(*,dash1_line)
          loop = nloop + 1      
        else
          write(*,"( /, 62('#'),/,                           /,&
            t27, ' Success! ',                               /,&
            t10, ' Final objective function value: ', e10.5, /,&
            t10, ' Maximum violation of target distance: ', f10.6, /,&
            t10, ' Maximum violation of the constraints: ', e10.5,/,&
            62('-'), /,&
            ' Please cite this work if Packmol was useful: ',/,&
       ' L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez, ',/,&
       ' PACKMOL: A package for building initial configurations ',/,&
       ' for molecular dynamics simulations. ',/,&
       ' Journal of Computational Chemistry, 30:2157-2164,2009.',&
            /,/,62('#'),/ )" ) bestf, fdist, frest
          write(*,*) '  Running time: ', etime(tarray) - time0,' seconds. ' 
          stop 
        end if
      end if

    end do gencanloop

    if(itype.le.ntype) then
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

  end do main

end program packmol
