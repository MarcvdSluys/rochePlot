program proche
  ! plots rochelobes for systems listed in fort.10
  !  for each graph: m1, m2, a, r1, r2
  ! m1, m2 = masses of left and right star, respectively
  ! a      = distance between stars (in solar radii)
  ! r1, r2 = radii of left and right stars in solar radii
  !  if (r1,r2) > 1.e5 the rochelobe is filled
  !  if (r1,r2) < 0.   a circle with radius (r1,r2) + disk is drawn
  !
  ! the plot is scaled automatically: to do this, first all
  !  required parameters are read, and lobe sizes and positions
  !  estimated
  !
  ! next, the individual graphs are made
  !
  parameter(npl=100,ng=10)
  integer :: nev
  dimension rm1(ng),rm2(ng),rsep(ng),rlag(ng),rlef(ng),rrig(ng), &
       hei(ng),rad1(ng),rad2(ng),xpl(npl),ypl(npl),ypl2(npl),xtl(4), &
       pb(ng),rmc(ng)
  integer :: narg,iargc
  character :: text*50,label(4)*50,yaa(8),bla,title*50,fname*50
  external rlimit,rline
  common/roche/ q,q11,const,const2,xsq,onexsq
  data label/'M\d1\u(M\d\(2281)\u)','M\d2\u(M\d\(2281)\u)', &
       'P\dorb\u(d)','M\dc\u(M\d\(2281)\u)'/
  data yaa/'c','d','e','f','g','h','g','h'/
  ! some physical constants
  data gravc,sunm,sunr,sunl,pi &
       /6.668e-8,1.989e33,6.96e10,3.846e33,3.1415926/
  !
  ! constant for orbital separation from mass and orbital period:
  csep=((24.*3600./2./pi)**2*gravc*sunm)**(1./3.)/sunr
  ! the necessary parameters are read from file; all together,
  !  to enable calculation of the overall size of the graph.
  !      write(6,*)'in test run, frame should be included: type 0'
  !      write(6,*)'otherwise: type -2'
  !      read(5,*)iaxis
  iaxis=-2
  !      write(6,*)'if M1, M2, and Pb are labels: type 3'
  !      write(6,*)'if to these the core mass is added: type 4'
  !      read(5,*)klabel
  !       klabel=3
  
  xmin = 0.
  xmax = 0.
  
  
  ! Read command-line variables
  narg = iargc()
  if(narg.eq.1) then
     call getarg(1,fname)
  else 
     fname = 'input.dat'
  end if
  
  
  print*,'Reading input file ',fname
  open(unit=10,form='formatted',status='old',file=fname)
  read(10,*)klabel
  read(10,*)nev
  if(nev.gt.ng) print*,'Increase the value of ng!'
  read(10,*)bla
  do itel=1,nev
     if(klabel.eq.3) then
        read(10,*,end=2)rm1(itel),rm2(itel),pb(itel),rad1(itel),rad2(itel)
     else
        read(10,*,end=2)rm1(itel),rm2(itel),pb(itel),rad1(itel),rad2(itel),rmc(itel)
     end if
     
     rsep(itel) = csep*((rm1(itel)+rm2(itel))*pb(itel)**2)**(1./3.)
     ktel = itel
     
     ! calculate inner lagrangian point, start with estimate
     q = rm1(ktel)/rm2(ktel)
     q11 = 1./(1.+q)
     x = 0.5+0.2222222*alog10(q)
1    fx = q/x/x-1./(1.-x)**2-(1.+q)*x+1.
     dfx = -2.*q/x**3-2./(1.-x)**3-(1.+q)
     dx = -fx/dfx/x
     x = x*(1.+dx)
     if(abs(dx).gt.1.e-6) goto 1
     rlag(ktel) = x
     ! set vertical space for graph equal to max(x,1-x)
     if(q.gt.1.)then
        hei(ktel) = x*rsep(ktel)
     else
        hei(ktel) = (1.-x)*rsep(ktel)
     end if
     
     ! calculate left limit of lobe (before shift)
     const = q/x+1./(1.-x)+0.5*(1.+q)*(x-q11)**2
     x1 = 1.5-0.5*x
     x2 = 2.0-x
     xacc = 1.e-4
     rrig(ktel) = rtsafe(rlimit,x1,x2,xacc)
     x1 = -0.5*x
     x2 = -x
     rlef(ktel) = rtsafe(rlimit,x1,x2,xacc)
     write(6,*)'Roche limits: ',rlef(ktel),rlag(ktel),rrig(ktel),hei(ktel)
     
     ! calculate limits after enlarging and shift, and keep track of maxima
     asep = rsep(ktel)
     xshift = -asep*rm2(ktel)/(rm1(ktel)+rm2(ktel))
     xleft = asep*rlef(ktel)+xshift
     if(xleft.lt.xmin) xmin = xleft
     xright = asep*rrig(ktel)+xshift
     if(xright.gt.xmax) xmax = xright
  end do
  
  
2 continue
  ! after all limits have been sampled, now calculate plot limits
  xmargin = 0.2*(xmax-xmin)
  ysize = 0.
  do i=1,ktel      
     ysize = ysize+hei(i)
  end do
  ysize = 2.5*ysize
  ymargin = 0.02*ysize
  xleft = xmin - xmargin
  xrigh = xmax + xmargin
  write(6,*)'Plot limits: ',xleft,xrigh,ysize
  
  read(10,*)iscr
  if(iscr.eq.0) then
     write(6,*)'Plot written to rochelobes.eps'
     call pgbegin(0,'rochelobes.eps/ps',1,1)
     call pgslw(2)
     call pgscf(1)
  else
     call pgbegin(1,'/xs',1,1)
  end if
  if(iscr.eq.1.or.iscr.eq.2) then     ! Create a white background; swap black (ci=0) and white (ci=1)
     call pgsci(0)
     call pgscr(0,1.,1.,1.)
     call pgscr(1,0.,0.,0.)
     call pgsvp(0.,1.,0.,1.)
     call pgswin(-1.,1.,-1.,1.)
     call pgrect(-2.,2.,-2.,2.)
     call pgsci(1)
  end if
  call pgsfs(1)
  
  
  
  call pgenv(xleft,xrigh,ysize,0.,1,iaxis) !iaxis: 0 for test, -2 for real (do or do not plot frame...?) 
  call pgscr(0,1.,1.,1.)            ! Repeat this, to make it work for AquaTerm, for which it was designed
  call pgscr(1,0.,0.,0.)
  call pgsci(0)
  call pgrect(xleft-10*abs(xleft-xrigh),xrigh+10*abs(xleft-xrigh),10*ysize,-10*ysize)
  call pgsci(1)
  
  !      write(6,*)'Length of size-bar? (integer in solar radii?)'
  !      read(5,*) ilen
  read(10,*) ilen
  !      if(klabel.eq.3) then
  !        write(6,*)'x-positions of: M1, M2, Pb,xaa'
  !      else
  !        write(6,*)'x-positions of: M1, M2, Pb, Mc'
  !      end if
  !      read(5,*)(xtl(k),k=1,klabel),xaa
  read(10,*) xtl(1)
  read(10,*) xtl(2)
  read(10,*) xtl(3)
  read(10,*) xtl(4)
  read(10,'(A50)') label(4)
  do kl=1,klabel
     if(xtl(kl).ne.0.) call pgptxt(xtl(kl),0.,0.,0.5,label(kl))
  end do
  read(10,'(A50)') title
  if(title(1:10).ne.'          ') then
     call pgsch(1.5)
     call pgptxt(0.,-3*ymargin,0.,0.5,title)
     call pgsch(1.)
  end if
  
  do itel=1,ktel
     xm1 = rm1(itel)
     xm2 = rm2(itel)
     asep = rsep(itel)
     x = rlag(itel)
     q = xm1/xm2
     q11 = 1./(1.+q)
     const = q/x+1./(1.-x)+0.5*(1.+q)*(x-q11)**2
     xpl(1) = rlef(itel)
     xpl(npl) = rrig(itel)
     ypl(1) = 0.
     ypl(npl) = 0.
     
     ! start on left lobe
     nl = npl/2-1
     dxl = (x-xpl(1))/nl
     do i = 2,nl
        xl = xpl(1)+(i-1)*dxl
        xsq = xl*xl
        onexsq = (1.-xl)**2
        const2 = 0.5*(1.+q)*(xl-q11)**2-const
        y1 = 0.
        y2 = x**2
        ysq = rtsafe(rline,y1,y2,xacc)
        xpl(i) = xl
        ypl(i) = sqrt(ysq)
     end do
     xpl(nl+1) = x
     ypl(nl+1) = 0.
     !        write(6,*)'left lobe done'
     
     ! and right lobe
     dxr = (xpl(npl)-x)/(nl+1)
     do i = 2,nl+1
        xl = xpl(nl+1)+(i-1)*dxr
        xsq = xl*xl
        onexsq = (1.-xl)**2
        const2 = 0.5*(1.+q)*(xl-q11)**2-const
        y1 = 0.
        y2 = (1-x)**2
        ysq = rtsafe(rline,y1,y2,xacc)
        xpl(nl+i) = xl
        ypl(nl+i) = sqrt(ysq)
     end do
     !        write(6,*)'right lobe done'
     
     ! now enlarge and shift lobes:
     xmult = asep
     xshift = -asep*xm2/(xm1+xm2)
     if(itel.eq.1) then
        yshift = hei(itel)+ymargin
     else
        yshift = yshift+hei(itel-1)+hei(itel)+ymargin
     end if
     do i=1,npl
        xpl(i) = xpl(i)*xmult+xshift
        swap = ypl(i)*xmult
        ypl(i) = swap + yshift
        ypl2(i) = -swap + yshift
     end do
     ! and plot them
     !        call pgline(npl,xpl,ypl)
     !        call pgline(npl,xpl,ypl2)
     
     ! start on stars, left first: (for use of rad1, see above, at begin)
     if(rad1(itel).gt.1.e5) then
        call pgsci(15)
        call pgpoly(nl+1,xpl,ypl)
        call pgpoly(nl+1,xpl,ypl2)
        call pgsci(1)
     else
        rad = rad1(itel)
        !          if(rad2(itel).gt.1.e5) then
        if(rad2(itel).gt.1.e5.and.rad1(itel).gt.0.) then
           radd = 0.7*asep*x
           call pgsci(15)
           call disk(xshift,yshift,rad,radd)
           call pgsci(1)
        end if
        !          if(abs(rad).lt.ysize*0.001) then
        !            call pgpoint(1,xshift,yshift,17)
        !          else
        !            call cirkel(xshift,yshift,abs(rad),40)
        !          end if
        call cirkel(xshift,yshift,max(abs(rad),ysize*0.002),40)
     end if
     
     ! Plot Roche lobe after star
     call pgline(npl,xpl,ypl)
     call pgline(npl,xpl,ypl2)
     
     ! right:
     if(rad2(itel).gt.1.e5) then
        do i=1,nl+2
           xpl(i) = xpl(i+nl)
           ypl(i) = ypl(i+nl)
           ypl2(i) = ypl2(i+nl)
        end do
        call pgsci(15)
        call pgpoly(nl+2,xpl,ypl)
        call pgpoly(nl+2,xpl,ypl2)
        call pgsci(1)
     else
        rad = rad2(itel)
        !          if(rad1(itel).gt.1.e5) then
        if(rad1(itel).gt.1.e5.and.rad2(itel).gt.0.) then
           radd = 0.7*asep*(1.-x)
           call pgsci(15)
           call disk(xshift+asep,yshift,rad,radd)
           call pgsci(1)
        end if
        !          if(abs(rad).lt.ysize*0.001) then
        !            call pgpoint(1,xshift+asep,yshift,17)
        !          else
        !            call cirkel(xshift+asep,yshift,abs(rad),40)
        !          end if
        call cirkel(xshift+asep,yshift,max(abs(rad),ysize*0.002),40)
     end if
     
     call pgline(nl+2,xpl,ypl)
     call pgline(nl+2,xpl,ypl2)
     
     ! and write labels
     if(klabel.eq.3) then
        write(label(1),102) rm1(itel)
        write(label(2),102) rm2(itel)
     else
        write(label(1),103) rm1(itel)
        write(label(2),103) rm2(itel)
        write(label(4),102) rmc(itel)
     end if
     write(label(3),104) pb(itel)
102  format(f7.3)
103  format(f5.2)
104  format(f7.2)
     do k=1,klabel
        call pgptxt(xtl(k),yshift,0.,0.5,label(k))
     end do
     !         call pgtext(xaa,yshift,yaa(itel))
  end do
  
  
  
  
  
  ! plot size bar
  xlen = ilen*1.
  xpl(2) = xlen/2.
  xpl(1) = -xpl(2)
  !      yshift = yshift+hei(ktel)+2*ymargin
  yshift = yshift+hei(ktel)+5*ymargin
  ypl(1) = yshift
  ypl(2) = ypl(1)
  call pgline(2,xpl,ypl)
  write(text,101)ilen
101 format(i5,'R\d\(2281)')
  !      call pgtext(xpl(2),ypl(2),text)
  call pgtext(xpl(2),ypl(2)+0.5*ymargin,text)
  
  ! and axis of rotation
  xpl(1) = 0.
  xpl(2) = 0.
  ypl(1) = 0.
  !      ypl(2) = yshift+ymargin
  ypl(2) = yshift-ymargin
  call pgsls(4)
  call pgline(2,xpl,ypl)
  call pgsls(1)
  
  ! add texts, if necessary
  !123   write(6,*)'give position (x,y) of text'
  !      write(6,*)'x=0. means: no text to be added'
  !      read(5,*)xt,yt
  !      if(xt.ne.0.) then
  !        write(6,*)'give text string'
  !        read(5,99)text
  !        call pgtext(xt,yt,text)
  !        goto 123
  !      end if
  read(10,*)xt
  read(10,*)yt
  read(10,99)text
  if(xt.ne.0.) call pgtext(xt,yt,text)
99 format(a)
  call pgend
end program proche











subroutine rlimit(x,f,df)
  ! calculates outer limit of roche-lobe
  common/roche/ q,q11,const,const2,xsq,onexsq
  r1=abs(x)
  r2=abs(1.-x)
  r3=abs(x-q11)
  f=q/r1+1./r2+0.5*(1.+q)*r3**2-const
  df=-q*x/r1**3+(1.-x)/r2**3+(1.+q)*(x-q11)
end subroutine rlimit


subroutine rline(y,f,df)
  ! calculates value of y^2 for x^2 value
  common/roche/ q,q11,const,const2,xsq,onexsq
  r1=sqrt(xsq+y)
  r2=sqrt(onexsq+y)
  f=q/r1+1./r2+const2
  df=-0.5*q/r1**3-0.5/r2**3
end subroutine rline


subroutine cirkel(xc,yc,rad,n)                                   
  dimension x(200),y(200)                                          
  step=6.2831852/(n-1)                                             
  do i=1,n                                                       
     phi=i*step                                                       
     x(i)=xc+rad*cos(phi)                                             
     y(i)=yc+rad*sin(phi)                                             
  end do
  call pgpoly(n,x,y)
end subroutine cirkel
  
  
subroutine disk(xc,yc,rad,rlen)                                 
  ! draws disk centered on xc,yc between rad and rlen
  real x(5),y(5)                                                  
  x(1)=xc+rad
  x(2)=x(1)                                                        
  x(3)=xc+rlen                                                     
  x(4)=x(3)                                                        
  x(5)=x(1)                                                        
  y(1)=yc+0.15*rad
  y(2)=yc-0.15*rad                                                
  y(3)=yc-0.15*rlen                                                     
  y(4)=yc+0.15*rlen                                                
  y(5)=y(1)                                                        
  call pgpoly(5,x,y)                                            
  x(1)=xc-rad
  x(2)=x(1)                                                        
  x(3)=xc-rlen                                                     
  x(4)=x(3)                                                        
  x(5)=x(1)                                                        
  y(5)=y(1)                                                        
  call pgpoly(5,x,y)                                            
end subroutine disk



function rtsafe(funcd,x1,x2,xacc)
  ! numerical recipes, p.258
  parameter(maxit=100)
  call funcd(x1,fl,df)
  call funcd(x2,fh,df)
  if(fl*fh.ge.0.) write(0,'(A)') 'root must be bracketed'
  if(fl.lt.0.)then
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
     swap=fl
     fl=fh
     fh=swap
  end if
  rtsafe=0.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call funcd(rtsafe,f,df)
  do j=1,maxit
     if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0..or. abs(2.*f).gt.abs(dxold*df) ) then
        dxold=dx
        dx=0.5*(xh-xl)
        rtsafe=xl+dx
        if(xl.eq.rtsafe) return
     else
        dxold=dx
        dx=f/df
        temp=rtsafe
        rtsafe=rtsafe-dx
        if(temp.eq.rtsafe) return
     end if
     if(abs(dx).lt.xacc) return
     call funcd(rtsafe,f,df)
     if(f.lt.0.)then
        xl=rtsafe
        fl=f
     else
        xh=rtsafe
        fh=f
     end if
  end do
  write(0,'(A)')' rtsafe exceeds maximum iterations'
end function rtsafe
