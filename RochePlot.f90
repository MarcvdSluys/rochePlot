!> \file  RochePlot.f90  Plots Roche lobes for given set of binaries
!!
!!  \mainpage RochePlot documentation
!!  RochePlot is a Fortran code using PGPlot to plot a series of binaries to illustrate the key moments in the evolution of a
!!  given binary.  The code was originally written by Frank Verbunt.
!!
!!  \par
!!  (c) 2003-2012, Frank Verbunt, Marc van der Sluys


!***********************************************************************************************************************************
!> \brief  Contains data from the input file and derived data

module input_data
  implicit none
  
  integer, parameter :: npl=100  ! Number of plotting points
  integer, parameter :: ng=10    ! Maximum number of binaries that can be plotted
  integer :: klabel, ktel
  integer :: blen, iscr
  
  real :: csep, rm1(ng),rm2(ng),pb(ng),rad1(ng),rad2(ng), age_mc(ng)
  real :: rsep(ng),rlag(ng),rlef(ng),rrig(ng), hei(ng)
  real :: xtl(5), xt,yt
  logical :: ce(ng)
  
  character :: txt(ng)*(50), label(5)*(50), text*(50), title*(50)
  
end module input_data
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief  Contains plot settings

module plot_settings
  use input_data, only: npl
  implicit none
  private :: npl
  
  integer :: lw
  real :: xpl(npl),ypl(npl),ypl2(npl)
  real :: xleft,xrigh, ymargin, ysize
  real :: xlen,yshift
  character :: outputfile*(50)
  logical :: use_colour
  
end module plot_settings
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief  Contains Roche-lobe data

module roche_data
  use input_data, only: ng
  implicit none
  private :: ng
  
  real :: q,q11, const1,const2,CEdiff(ng), xsq,onexsq
  
end module roche_data
!***********************************************************************************************************************************




!***********************************************************************************************************************************
!> \brief  Plots Roche lobes for given binaries
!!
!! For each graph, read m1, m2, a, r1, r2:
!! - m1, m2 = masses of left and right star, respectively
!! - a      = distance between stars (in solar radii)
!! - r1, r2 = radii of left and right stars in solar radii
!!
!! - if (r1,r2) > 1.e5 the Roche lobe is filled
!! - if both stars fill their Roche lobes, a common envelope is assumed
!! - if (r1,r2) < 0.   a circle with radius (r1,r2) + disc is drawn
!!
!! The plot is scaled automatically: 
!! - first all required parameters are read
!! - then lobe sizes and positions are estimated
!! - finally, the individual graphs are made

program rocheplot
  use input_data, only: label,csep,ktel, text,xt,yt
  use plot_settings, only: use_colour, outputfile
  
  implicit none
  integer :: i,itel, command_argument_count
  real :: gravc,sunm,sunr,twopi
  character :: inputfile*(50)
  
  
  !*** Initialise code:
  use_colour = .false.  ! B/W
  use_colour = .true.   ! Use colour
  
  ! Column headers:
  label = [character(len=50) :: 'M\d1\u(M\d\(2281)\u)','M\d2\u(M\d\(2281)\u)', 'P\dorb\u(d)','M\dc\u(M\d\(2281)\u)', '']
  
  ! Constant for orbital separation from mass and orbital period:
  gravc = 6.668e-8      ! G
  sunm  = 1.989e33      ! Mo
  sunr  = 6.96e10       ! Ro
  twopi = 8.*atan(1.0)  ! 2pi
  csep = ((24.*3600./twopi)**2 * gravc*sunm)**(1./3.) / sunr
  
  
  !*** Read command-line parameters, set input and output filenames:
  inputfile = 'input.dat'
  outputfile = 'rochelobes.eps'
  if(command_argument_count().eq.1) then
     call get_command_argument(1, inputfile)
     i = index(inputfile,'.dat', back=.true.)                      ! Find last dot in input file name
     if(i.gt.0.and.i.lt.50) outputfile = inputfile(1:i-1)//'.eps'  ! When file.dat is input, use file.eps as output
  end if
  
  
  !*** Read the input file:
  call read_input_file(trim(inputfile))
  
  
  !*** Initialise plot output; open output file, set page, create a white background, print plot title and column headers:
  call initialise_plot()
  
  
  !*** Plot the different binaries:
  do itel=1,ktel
     call plot_binary(itel)  ! Plot each binary; Roche lobes, stars and labels
  end do  ! do itel = 1,ktel
  
  
  
  !*** Finish plot:
  
  ! Plot scale bar:
  call plot_scale_bar()
  
  ! Plot axis of rotation:
  call plot_rotation_axis()
  
  ! Add texts, if necessary:
  if(xt.ne.0.) call pgtext(xt,yt,text)
  
  call pgend
  
end program rocheplot
!***********************************************************************************************************************************




!***********************************************************************************************************************************
!> \brief  Read the lines of the input file containting evolutionary states and compute positions of the Roche lobes

subroutine read_input_file(inputfile)
  use input_data
  use plot_settings, only: xleft,ysize,ymargin,xrigh
  use roche_data, only: q,q11, const1,CEdiff
  
  implicit none
  character, intent(in) :: inputfile*(*)
  
  integer :: io, itel, ki, nev
  real :: asep, dfx,dx,fx, rtsafe
  real :: x,x1,x2,xacc,xright,xshift, xmargin, xmin,xmax
  character :: tmpstr
  
  external :: rlimit
  
  
  ! Open input file:
  write(*,'(/,A,/)') ' Reading input file: '//trim(inputfile)
  open(unit=10,form='formatted',status='old',file=trim(inputfile))
  
  read(10,*) klabel            ! Number of labels per line - currently 3, 4 or 5
  read(10,*) nev               ! Number of evolutionary phases to plot = number of data lines in input file
  if(nev.gt.ng) write(0,'(A)') 'Increase the value of ng!'
  
  read(10,*) tmpstr
  tmpstr = tmpstr  ! Remove 'unused' compiler warnings
  
  
  
  ! Read the lines containting evolutionary states and compute positions/limits of the Roche lobes:
  xmin =  huge(xmin)
  xmax = -huge(xmax)
  do itel=1,nev
     select case(klabel)
     case(3)
        read(10,*, iostat=io) rm1(itel), rm2(itel), pb(itel), rad1(itel), rad2(itel)
     case(4)
        read(10,*, iostat=io) rm1(itel), rm2(itel), pb(itel), rad1(itel), rad2(itel), age_mc(itel)
     case(5)
        read(10,*, iostat=io) rm1(itel), rm2(itel), pb(itel), rad1(itel), rad2(itel), age_mc(itel), txt(itel)
     case default
        write(0,'(A,I3,A)') ' klabel =',klabel,' not supported, change the value in your input file.'
        stop
     end select
     
     if(io.lt.0) return  ! end of file
     
     rsep(itel) = csep * ((rm1(itel)+rm2(itel)) * pb(itel)**2)**(1./3.)  ! Kepler: P_orb -> a_orb
     ktel = itel
     
     
     q = rm1(ktel)/rm2(ktel)
     q11 = 1.0/(1.0+q)
     
     ! Calculate inner Lagrangian point, start with estimate:
     x = 0.5 + 0.2222222*log10(q)
     dx = huge(dx)
     do while(abs(dx).gt.1.e-6)
        fx = q/(x**2) - 1.0/(1.0-x)**2 - (1.0+q)*x + 1.0
        dfx = -2.0*q/x**3 - 2.0/(1.-x)**3 - (1.0+q)
        dx = -fx/(dfx*x)
        x = x*(1.0+dx)
     end do
     
     rlag(ktel) = x  ! Inner Lagrangian point
     
     
     ! Set vertical space for graph equal to max(x,1-x):
     if(q.gt.1.) then
        hei(ktel) = x*rsep(ktel)
     else
        hei(ktel) = (1.-x)*rsep(ktel)
     end if
     
     
     ! Compute the potential difference to subtract from the Roche potential for CE plots:
     CEdiff(itel) = 0.33/q
     if(q.lt.1) CEdiff(itel) = 0.33*q
     
     ce(itel) = .false.
     !ce(itel) = .true.
     if(min(rad1(itel),rad2(itel)).gt.1.e5) ce(itel) = .true.
     
     ! Calculate limits of lobes (before shift):
     const1 = q/x + 1./(1.-x) + 0.5*(1.+q)*(x-q11)**2
     if(ce(itel)) const1 = const1 - CEdiff(itel)    ! CE
     
     xacc = 1.e-4
     x1 = 1.5 - 0.5*x
     x2 = 2.0 - x
     rrig(ktel) = rtsafe(rlimit, x1,x2, xacc)  ! Right limit
     
     x1 = -0.5*x
     x2 = -x
     rlef(ktel) = rtsafe(rlimit, x1,x2, xacc)  ! Left limit
     
     
     write(*,'(A,I3,A1,4G12.3)') ' Roche limits binary',ktel,':', rlef(ktel), rlag(ktel), rrig(ktel), hei(ktel)
     
     
     ! Calculate limits after enlarging and shift, and keep track of minima and maxima:
     asep   = rsep(ktel)  ! Orbital separation
     xshift = -asep*rm2(ktel) / (rm1(ktel)+rm2(ktel))
     
     xleft = asep*rlef(ktel) + xshift
     xmin = min(xmin,xleft)
     
     xright = asep*rrig(ktel) + xshift
     xmax = max(xmax,xright)
  end do
  
  
  ! After all limits have been sampled, now calculate plot limits:
  ! - silly: if bar falls off plot, increase ysize
  
  xmargin = 0.2*(xmax-xmin)
  ysize = 0.
  do ki=1,ktel      
     ysize = ysize + hei(ki)
  end do
  
  ysize = 2.5*ysize*1.25
  ymargin = 0.02*ysize
  
  xleft = xmin - xmargin
  xrigh = xmax + xmargin*4.
  
  write(6,'(/,A,3F12.3)') ' Plot limits: ',xleft,xrigh,ysize
  
  
  
  ! Read the rest of the input file:
  read(10,*) iscr
  read(10,*) blen  ! Length of the scale bar
  do ki=1,klabel
     read(10,*) xtl(ki)  ! Column headers
  end do
  read(10,'(A50)') label(4)
  read(10,'(A50)') title  ! Plot title
  
  read(10,*) xt
  read(10,*) yt
  read(10,'(A)') text
  close(10)
  
  
end subroutine read_input_file
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Initialise plot output; open output file, set page, create a white background, print plot title and column headers

subroutine initialise_plot()
  use input_data, only: klabel, label,iscr,xtl,title
  use plot_settings, only: use_colour, xleft,xrigh,ysize,ymargin, outputfile, lw
  
  implicit none
  integer :: iaxis, kl
  
  ! the necessary parameters are read from file; all together,
  !  to enable calculation of the overall size of the graph.
  iaxis=-2  ! Draft: 0,  quality: -2
  
  
  if(iscr.eq.0) then
     write(6,'(/,A,/)')' Saving plot as '//trim(outputfile)
     if(use_colour) then
        call pgbegin(0,''//trim(outputfile)//'/cps',1,1)
     else
        call pgbegin(0,''//trim(outputfile)//'/ps',1,1)
     end if
     lw = 2
     call pgscf(1)
  else
     call pgbegin(1,'/xs',1,1)
     lw = 1
  end if
  
  
  if(iscr.eq.1.or.iscr.eq.2) call pgwhitebg()  ! Create a white background when plotting to screen; swap fg/bg colours
  call pgsfs(1)
  call pgslw(lw)
  
  call pgenv(xleft,xrigh,ysize,0., 1, iaxis)
  call pgsci(1)
  
  
  ! Print plot title:
  if(title(1:10).ne.'          ') then
     call pgsch(1.5)
     call pgslw(3*lw)
     call pgptxt(0.,-3*ymargin,0.,0.5,trim(title))
     call pgsch(1.)
  end if
  
  
  ! Print column headers:
  call pgslw(2*lw)
  do kl=1,klabel
     if(xtl(kl).ne.0.) call pgptxt(xtl(kl),0.,0.,0.5,trim(label(kl)))
  end do
  call pgslw(lw)
  
end subroutine initialise_plot
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief  Create a white background when plotting to screen; swap black (ci=0) and white (ci=1)

subroutine pgwhitebg()
  implicit none
  
  call pgsci(0)
  call pgscr(0,1.,1.,1.)
  call pgscr(1,0.,0.,0.)
  call pgsvp(0.,1.,0.,1.)
  call pgswin(-1.,1.,-1.,1.)
  call pgrect(-2.,2.,-2.,2.)
  call pgsci(1)

end subroutine pgwhitebg
!***********************************************************************************************************************************











!***********************************************************************************************************************************
!> \brief  Draw a circle
!!
!! \todo  Replace by pgcirc()?  -  Perhaps this looks nicer?

subroutine plot_circle(xc,yc,rad,n)
  implicit none
  real, intent(in) :: xc,yc,rad
  integer, intent(in) :: n
  
  integer :: i
  real :: x(n),y(n), step,phi
  
  step = 6.2831852/real(n-1)
  
  do i=1,n
     phi  = real(i)*step
     x(i) = xc + rad*cos(phi)
     y(i) = yc + rad*sin(phi)
  end do
  
  call pgpoly(n,x,y)
  
end subroutine plot_circle
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Draws an accretion disc centered on xc,yc between rad and rlen

subroutine plot_disc(xc,yc, rad,rlen)
  implicit none
  
  real, intent(in) :: xc,yc,rad,rlen
  real :: x(5),y(5), flare
  
  flare = 0.15  ! Disc's flare
  
  ! Draw right half:
  x(1) = xc + rad
  x(2) = x(1)
  x(3) = xc + rlen
  x(4) = x(3)
  x(5) = x(1)
  y(1) = yc + flare*rad
  y(2) = yc - flare*rad
  y(3) = yc - flare*rlen
  y(4) = yc + flare*rlen
  y(5) = y(1)
  call pgpoly(5,x,y)
  
  ! Draw left half:
  x(1) = xc - rad
  x(2) = x(1)
  x(3) = xc - rlen
  x(4) = x(3)
  x(5) = x(1)
  y(5) = y(1)
  call pgpoly(5,x,y)
  
end subroutine plot_disc
!***********************************************************************************************************************************










!***********************************************************************************************************************************
!> \brief  Plot each binary; Roche lobes, stars and labels
!!
!! \param itel  Number of the current binary/evolutionary state (1-ktel)

subroutine plot_binary(itel)
  use input_data, only: npl, age_mc, rm1,rm2,rsep,rlag,rlef,rrig,hei,rad1,rad2,klabel,label,ktel,txt, pb,xtl, ce
  use plot_settings, only: xpl,ypl,ypl2, use_colour, ysize,ymargin, yshift
  use roche_data, only: q,q11, const1,const2,CEdiff, xsq,onexsq
  
  implicit none
  integer, intent(in) :: itel
  integer :: il,pl,k,nl
  real :: asep, rad,radd,swap, rtsafe
  real :: x,xacc,xl,xm1,xm2, xmult,xshift
  real :: y1,y2,ysq, xtmp,ytmp, dy, xmap
  
  external :: rline
  
  xm1 = rm1(itel)     ! M1
  xm2 = rm2(itel)     ! M2
  asep = rsep(itel)   ! Orbital separation
  x = rlag(itel)      ! Inner Lagrangian point
  q = xm1/xm2         ! q1
  q11 = 1./(1.+q)     ! M2/Mtot
  
  const1 = q/x + 1./(1.-x) + 0.5*(1.+q)*(x-q11)**2
  if(ce(itel)) const1 = const1 - CEdiff(itel)  ! CE
  
  xpl(1)   = rlef(itel)  ! Left limit of Rl
  xpl(npl) = rrig(itel)  ! Right limit of Rl
  ypl(1)   = 0.            
  ypl(npl) = 0.          
  
  
  nl   = npl/2-1
  xacc = 1.e-4
  
  
  ! Compute left lobe:
  do il = 2,nl
     xl = xmap(il,nl, x,xpl(1), 1)  ! Map x-points more densely near outer Rl limit; 1=left
     
     xsq = xl*xl
     onexsq = (1.-xl)**2
     const2 = 0.5*(1.+q)*(xl-q11)**2 - const1
     
     y1 = 0.
     y2 = x**2
     ysq = rtsafe(rline,y1,y2,xacc)
     
     xpl(il) = xl
     ypl(il) = sqrt(ysq)
  end do
  xpl(nl+1) = x
  ypl(nl+1) = 0.
  
  ! Compute right lobe:
  do il = 2,nl+1
     xl = xmap(il,nl, x,xpl(npl), 2)  ! Map x-points more densely near outer Rl limit; 2=right
     
     xsq = xl*xl
     onexsq = (1.-xl)**2
     const2 = 0.5*(1.+q)*(xl-q11)**2 - const1
     
     y1 = 0.
     y2 = (1.0-x)**2
     ysq = rtsafe(rline,y1,y2,xacc)
     
     xpl(nl+il) = xl
     ypl(nl+il) = sqrt(ysq)
  end do
  
  
  ! Enlarge and shift lobes:
  xmult = asep
  xshift = -asep*xm2/(xm1+xm2)
  if(itel.eq.1) then
     yshift = hei(itel) + ymargin
  else
     yshift = yshift + hei(itel-1) + hei(itel) + ymargin
  end if
  do pl=1,npl
     xpl(pl)  = xpl(pl)*xmult + xshift
     swap     = ypl(pl)*xmult
     ypl(pl)  =  swap + yshift
     ypl2(pl) = -swap + yshift
  end do
  
  
  ! Plot left star/disc:
  if(.not.ce(itel)) then
     if(rad1(itel).gt.1.e5) then  ! Rl filling
        call pgsci(15)
        if(use_colour) call pgsci(2)  ! red
        call pgpoly(nl+1, xpl, ypl)    ! Bottom half
        call pgpoly(nl+1, xpl, ypl2)  ! Top half
        call pgsci(1)
     else
        rad = rad1(itel)
        if(rad2(itel).gt.1.e5.and.rad1(itel).gt.0.) then
           radd = 0.7*asep*x
           call pgsci(15)
           if(use_colour) call pgsci(5)  ! light blue
           call plot_disc(xshift,yshift,rad,radd)
           call pgsci(1)
        end if
        call plot_circle(xshift,yshift,max(abs(rad),ysize*0.002),40)
     end if
  end if
  
  
  ! Plot CE:
  if(ce(itel)) then
     ! Make sure Rl contour does not go through L1 point:
     xtmp = xpl(nl+1)
     xpl(nl+1)  = xpl(nl)
     
     ytmp = ypl(nl+1)
     ypl(nl+1)  = ypl(nl)
     ypl2(nl+1) = ypl2(nl)
     
     ! Plot both Roche lobe(s):
     call pgsci(2)
     dy = ysize*0.0005               ! Create some overlap between the two halves
     call pgpoly(npl,xpl,ypl  - dy)  ! RL Bottom
     call pgpoly(npl,xpl,ypl2 + dy)  ! RL Top
     
     call pgsci(1)
     call pgcirc(xshift,yshift,ysize*0.002)  ! Core
  end if
  
  call pgline(npl,xpl,ypl)   ! RL bottom contour
  call pgline(npl,xpl,ypl2)  ! RL top contour
  
  
  if(ce(itel)) then
     xpl(nl+1) = xtmp
     ypl(nl+1) = ytmp
  end if
  
  
  ! Plot right star/disc:
  if(.not.ce(itel)) then
     if(rad2(itel).gt.1.e5) then  ! Rl filling
        do il=1,nl+2
           xpl(il)  = xpl(il+nl)
           ypl(il)  = ypl(il+nl)
           ypl2(il) = ypl2(il+nl)
        end do
        call pgsci(15)
        if(use_colour) call pgsci(2)  ! red
        call pgpoly(nl+2, xpl, ypl)   ! Bottom half
        call pgpoly(nl+2, xpl, ypl2)  ! Top half
        call pgsci(1)
     else
        rad = rad2(itel)
        if(rad1(itel).gt.1.e5.and.rad2(itel).gt.0.) then
           radd = 0.7*asep*(1.-x)
           call pgsci(15)
           if(use_colour) call pgsci(5)  ! light blue
           call plot_disc(xshift+asep,yshift,rad,radd)
           call pgsci(1)
        end if
        
        call plot_circle(xshift+asep,yshift,max(abs(rad),ysize*0.002),40)
     end if
  end if
  
  
  if(ce(itel)) then
     call pgcirc(xshift+asep,yshift,ysize*0.002)  ! Plot core
     
     ! Make sure Rl contour doesn't go through L1 point: - not needed when not plotting right Rl (again)
     !ytmp = ypl(nl+1)
     !ypl(nl+1)  = ypl(nl)
     !ypl2(nl+1) = ypl2(nl)
  end if
  
  
  ! Plot right Roche lobe:
  if(.not.ce(itel)) then
     call pgline(nl+2,xpl,ypl)
     call pgline(nl+2,xpl,ypl2)
  end if
  
  
  
  ! Write labels:
  if(klabel.eq.3) then
     write(label(1),'(F7.3)') rm1(itel)
     write(label(2),'(F7.3)') rm2(itel)
  else
     write(label(1),'(F5.2)') rm1(itel)
     write(label(2),'(F5.2)') rm2(itel)
     if(maxval(age_mc(1:ktel)).lt.2.) then
        write(label(4),'(F7.3)') age_mc(itel)
     else if(maxval(age_mc(1:ktel)).lt.50.) then
        write(label(4),'(F6.2)') age_mc(itel)
     else
        write(label(4),'(I3)') nint(age_mc(itel))
     end if
     if(klabel.ge.5) write(label(5),'(A)') trim(txt(itel))
  end if
  write(label(3),'(F7.2)') pb(itel)
  
  do k=1,klabel
     if(k.eq.5) then
        call pgptxt(xtl(k),yshift, 0.,0.0, trim(label(k)))  ! Align left
     else
        call pgptxt(xtl(k),yshift, 0.,0.5, trim(label(k)))  ! Align centre
     endif
  end do
  
end subroutine plot_binary
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief  Map the x-array to plot a right-hand Roche lobe
!!
!! \param il     Current point number in array (1-nl)
!! \param nl     Total number of points in array
!! \param l1     L1 position
!! \param lim    Extreme extent of Roche lobe (farthest from L1)
!! \param lr     Left (1) or right (2) Roche lobe
!!
!! \retval xmap  Mapped position of il-th x value

function xmap(il,nl, l1,lim, lr)
  implicit none
  integer, intent(in) :: il, nl, lr
  real, intent(in) :: l1,lim
  real :: xmap, pio2
  
  pio2 = 2*atan(1.)   ! pi/2
  
  if(lr.eq.1) then
     xmap = (1.0-cos(real(il-1)/real(nl)*pio2))   * (l1-lim) + lim
  else if(lr.eq.2) then
     xmap = cos((1.0-real(il-1)/real(nl+1))*pio2) * (lim-l1) + l1
  else
     write(0,'(/,A,I4,/)') '*** xmap(): ERROR:  parameter lr should be 1 or 2, not',lr
     stop
  end if
  
end function xmap
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief  Find the root of a function bracketed by x1,x2 using a combination of a Newton-Raphson and bisection methods
!!
!! \param funcd  User-provided function
!! \param x1     Lower limit for solution
!! \param x2     Upper limit for solution
!! \param xacc   Desired accuracy for solution
!!
!! \see Numerical recipes, par.9.4 (p.258 / 359)

function rtsafe(funcd, x1,x2, xacc)
  implicit none
  integer, parameter :: maxit=100
  integer :: j
  real, intent(in) :: x1,x2,xacc
  real :: rtsafe, dx,dxold,xh,xl, f,df,fh,fl, swap,temp
  
  
  call funcd(x1,fl,df)
  call funcd(x2,fh,df)
  if(fl*fh.ge.0.) write(0,'(2(A,2ES12.3))') ' rtsafe(): root must be bracketed:  x1,x2:',x1,x2, '  fl,fh:',fl,fh
  if(fl.lt.0.) then
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
     swap=fl
     fl=fh
     fh=swap
  end if
  
  rtsafe = 0.5*(x1+x2)
  dxold  = abs(x2-x1)
  dx     = dxold
  
  call funcd(rtsafe,f,df)
  
  do j=1,maxit
     if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0. .or. abs(2.*f).gt.abs(dxold*df) ) then
        dxold = dx
        dx = 0.5*(xh-xl)
        rtsafe = xl+dx
        if(xl.eq.rtsafe) return
     else
        dxold = dx
        dx = f/df
        temp = rtsafe
        rtsafe = rtsafe-dx
        if(temp.eq.rtsafe) return
     end if
     if(abs(dx).lt.xacc) return
     call funcd(rtsafe,f,df)
     if(f.lt.0.) then
        xl = rtsafe
        fl = f
     else
        xh = rtsafe
        fh = f
     end if
  end do
  
  write(0,'(A)')' rtsafe() exceeded maximum number of iterations'
  
end function rtsafe
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief  Calculates outer limit of Roche lobe
!!
!! \param x   Position along binary axis
!! \param f   Roche potential
!! \param df  First derivative of f w.r.t. x

subroutine rlimit(x, f,df)
  use roche_data, only: q,q11, const1
  
  implicit none
  real, intent(in) :: x
  real, intent(out) :: f,df
  
  real :: r1,r2,r3
  
  r1 = abs(x)
  r2 = abs(1.-x)
  r3 = abs(x-q11)
  
  f  = q/r1 + 1./r2 + 0.5*(1.+q)*r3**2 - const1
  df = -q*x/r1**3 + (1.-x)/r2**3 + (1.+q)*(x-q11)
  
end subroutine rlimit
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Calculates value of y^2 for given x^2 value
!!
!! \param y   Position along binary axis
!! \param f   Position of Roche surface
!! \param df  First derivative of f w.r.t. x

subroutine rline(y, f,df)
  use roche_data, only: q, const2, xsq,onexsq
  implicit none
  real, intent(in) :: y
  real, intent(out) :: f,df
  
  real :: r1,r2
  
  r1 = sqrt(y + xsq)
  r2 = sqrt(y + onexsq)
  
  f  = q/r1 + 1./r2 + const2
  df = -0.5*q/r1**3 - 0.5/r2**3
  
end subroutine rline
!***********************************************************************************************************************************







!***********************************************************************************************************************************
!> \brief  Plot scale bar

subroutine plot_scale_bar()
  use input_data, only: blen, hei, ktel
  use plot_settings, only: yshift, ymargin
  
  implicit none
  real :: xlen,xpl(2),ypl(2)
  character :: text*(99)
  
  xlen = real(blen)
  xpl(2) = xlen/2.
  xpl(1) = -xpl(2)
  
  !yshift = yshift + hei(ktel) + 2*ymargin
  yshift = yshift + hei(ktel) + 5*ymargin
  ypl(1) = yshift
  ypl(2) = ypl(1)
  
  call pgline(2,xpl,ypl)
  
  ! Print label:
  write(text,'(I5,"R\d\(2281)")') blen
  call pgtext(xpl(2),ypl(2)+0.5*ymargin,text)
  
end subroutine plot_scale_bar
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief  Plot axis of rotation for the binaries

subroutine plot_rotation_axis()
  use plot_settings, only: yshift, ymargin
  
  implicit none
  real :: xpl(2),ypl(2)
  
  xpl(1) = 0.
  xpl(2) = 0.
  ypl(1) = 0.
  
  !ypl(2) = yshift+ymargin
  ypl(2) = yshift-ymargin
  
  call pgsls(4)
  call pgline(2,xpl,ypl)
  call pgsls(1)
  
end subroutine plot_rotation_axis
!***********************************************************************************************************************************

