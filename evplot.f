       program proche
c plots rochelobes for systems listed in fort.10
c  for each graph: m1, m2, a, r1, r2
c m1, m2 = masses of left and right star, respectively
c a      = distance between stars (in solar radii)
c r1, r2 = radii of left and right stars in solar radii
c  if (r1,r2) > 1.e5 the rochelobe is filled
c  if (r1,r2) < 0.   a circle with radius (r1,r2) + disk is drawn
c
c the plot is scaled automatically: to do this, first all
c  required parameters are read, and lobe sizes and positions
c  estimated
c
c next, the individual graphs are made
c
       parameter(npl=100,ng=10)
       dimension rm1(ng),rm2(ng),rsep(ng),rlag(ng),rlef(ng),rrig(ng),
     c  hei(ng),rad1(ng),rad2(ng),xpl(npl),ypl(npl),ypl2(npl),xtl(4),
     c  pb(ng),rmc(ng)
       character*50 text,label(4)
       character*1 yaa(8)
       external rlimit,rline
       common/roche/ q,q11,const,const2,xsq,onexsq
       data label/'M\d1\u(M\d\(2281)\u)','M\d2\u(M\d\(2281)\u)',
     c  'P\db\u(d)','M\dc\u(M\d\(2281)\u)'/
       data yaa/'c','d','e','f','g','h','g','h'/
c some physical constants
       data gravc,sunm,sunr,sunl,pi
     c /6.668e-8,1.989e33,6.96e10,3.846e33,3.1415926/
c
c constant for orbital separation from mass and orbital period:
       csep=((24.*3600./2./pi)**2*gravc*sunm)**(1./3.)/sunr
c the necessary parameters are read from file; all together,
c  to enable calculation of the overall size of the graph.
c       write(6,*)'in test run, frame should be included: type 0'
c       write(6,*)'otherwise: type -2'
c       read(5,*)iaxis
        iaxis=-2
c       write(6,*)'if M1, M2, and Pb are labels: type 3'
c       write(6,*)'if to these the core mass is added: type 4'
c       read(5,*)klabel
        klabel=3
       xmin=0.
       xmax=0.
!       open (unit=10,form='formatted',status='old',file='input.dat')
       do itel=1,ng
        if(klabel.eq.3)then
        read(10,*,end=2)rm1(itel),rm2(itel),pb(itel),rad1(itel),
     c    rad2(itel)
        else
        read(10,*,end=2)rm1(itel),rm2(itel),pb(itel),rad1(itel),
     c    rad2(itel),rmc(itel)
        endif
        rsep(itel)=csep*((rm1(itel)+rm2(itel))*pb(itel)**2)**(1./3.) 
        ktel=itel
c calculate inner lagrangian point, start with estimate
        q=rm1(ktel)/rm2(ktel)
        q11=1./(1.+q)
        x=0.5+0.2222222*alog10(q)
1       fx=q/x/x-1./(1.-x)**2-(1.+q)*x+1.
        dfx=-2.*q/x**3-2./(1.-x)**3-(1.+q)
        dx=-fx/dfx/x
        x=x*(1.+dx)
        if(abs(dx).gt.1.e-6)goto1
        rlag(ktel)=x
c set vertical space for graph equal to max(x,1-x)
        if(q.gt.1.)then
         hei(ktel)=x*rsep(ktel)
        else
         hei(ktel)=(1.-x)*rsep(ktel)
        endif        
c calculate left limit of lobe (before shift)
        const=q/x+1./(1.-x)+0.5*(1.+q)*(x-q11)**2
        x1=1.5-0.5*x
        x2=2.0-x
        xacc=1.e-4
        rrig(ktel) = rtsafe(rlimit,x1,x2,xacc)
        x1=-0.5*x
        x2=-x
        rlef(ktel) = rtsafe(rlimit,x1,x2,xacc)
        write(6,*)'roche limits: ',rlef(ktel),rlag(ktel),rrig(ktel),
     c   hei(ktel)
c calculate limits after enlarging and shift, and keep track of maxima
        asep=rsep(ktel)
        xshift=-asep*rm2(ktel)/(rm1(ktel)+rm2(ktel))
        xleft=asep*rlef(ktel)+xshift
        if(xleft.lt.xmin)xmin=xleft
        xright=asep*rrig(ktel)+xshift
        if(xright.gt.xmax)xmax=xright
       enddo
2      continue
c after all limits have been sampled, now calculate plot limits
       xmargin=0.2*(xmax-xmin)
       ysize=0.
       do i=1,ktel      
        ysize=ysize+hei(i)
       enddo
       ysize=2.5*ysize
       ymargin=0.02*ysize
       xleft=xmin-xmargin
       xrigh=xmax+xmargin
       write(6,*)'plotlimits: ',xleft,xrigh,ysize
c start plotting
       write(6,*)'to plot on screen, type 1'
       read(5,*)iscr
       if(iscr.eq.1)then
        call pgbegin(0,'/xw',1,1)
       else
        write(6,*)'plot written to pgplot.ps'
        call pgbegin(0,'/ps',1,1)
        call pgslw(2)
       endif       
       call pgsci(1)
       call pgsfs(1)
       call pgenv(xleft,xrigh,ysize,0.,1,iaxis)     
       write(6,*)'length of size-bar? (integer in solar radii?)'
       read(5,*) ilen
       if(klabel.eq.3)then
       write(6,*)'x-positions of: M1, M2, Pb,xaa'
       else
       write(6,*)'x-positions of: M1, M2, Pb, Mc'
       endif
       read(5,*)(xtl(k),k=1,klabel),xaa
       do kl=1,klabel
        if(xtl(kl).ne.0.)then
        call pgtext(xtl(kl),0.,label(kl))
        endif
       enddo
       do itel=1,ktel
        xm1=rm1(itel)
        xm2=rm2(itel)
        asep=rsep(itel)
        x=rlag(itel)
        q=xm1/xm2
        q11=1./(1.+q)
        const=q/x+1./(1.-x)+0.5*(1.+q)*(x-q11)**2
        xpl(1)=rlef(itel)
        xpl(npl)=rrig(itel)
        ypl(1)=0.
        ypl(npl)=0.
c start on left lobe
        nl=npl/2-1
        dxl=(x-xpl(1))/nl
         do i=2,nl
          xl=xpl(1)+(i-1)*dxl
          xsq=xl*xl
          onexsq=(1.-xl)**2
          const2=0.5*(1.+q)*(xl-q11)**2-const
          y1=0.
          y2=x**2
          ysq = rtsafe(rline,y1,y2,xacc)
          xpl(i)=xl
          ypl(i)=sqrt(ysq)
         enddo
        xpl(nl+1)=x
        ypl(nl+1)=0.
        write(6,*)'left lobe done'
c and right lobe
        dxr=(xpl(npl)-x)/(nl+1)
         do i=2,nl+1
          xl=xpl(nl+1)+(i-1)*dxr
          xsq=xl*xl
          onexsq=(1.-xl)**2
          const2=0.5*(1.+q)*(xl-q11)**2-const
          y1=0.
          y2=(1-x)**2
          ysq = rtsafe(rline,y1,y2,xacc)
          xpl(nl+i)=xl
          ypl(nl+i)=sqrt(ysq)
         enddo
        write(6,*)'right lobe done'
c now enlarge and shift lobes:
         xmult=asep
         xshift=-asep*xm2/(xm1+xm2)
         if(itel.eq.1)then
          yshift=hei(itel)+ymargin
         else
          yshift=yshift+hei(itel-1)+hei(itel)+ymargin
         endif
         do i=1,npl
          xpl(i)=xpl(i)*xmult+xshift
          swap=ypl(i)*xmult
          ypl(i)=swap + yshift
          ypl2(i)=-swap + yshift
         enddo
c and plot them
         call pgline(npl,xpl,ypl)
         call pgline(npl,xpl,ypl2)
c start on stars, left first: (for use of rad1, see above, at begin)
         if(rad1(itel).gt.1.e5)then
          call pgpoly(nl+1,xpl,ypl)
          call pgpoly(nl+1,xpl,ypl2)
         else
          rad=rad1(itel)
          if(rad2(itel).gt.1.e5)then
           radd=0.7*asep*x
           call disk(xshift,yshift,rad,radd)
          endif
          if(rad.lt.ysize/500.)then
           call pgpoint(1,xshift,yshift,17)
          else
           call cirkel(xshift,yshift,rad,40)
          endif
         endif          
c right:
         if(rad2(itel).gt.1.e5)then
          do i=1,nl+2
           xpl(i)=xpl(i+nl)
           ypl(i)=ypl(i+nl)
           ypl2(i)=ypl2(i+nl)
          enddo
          call pgpoly(nl+2,xpl,ypl)
          call pgpoly(nl+2,xpl,ypl2)
         else
          rad=rad2(itel)
          if(rad1(itel).gt.1.e5)then
           radd=0.7*asep*(1.-x)
           call disk(xshift+asep,yshift,rad,radd)
          endif
          if(rad.lt.ysize/500.)then
           call pgpoint(1,xshift+asep,yshift,17)
          else
           call cirkel(xshift+asep,yshift,rad,40)
          endif
         endif          
c and write labels
         if(klabel.eq.3)then
          write(label(1),102)rm1(itel)
          write(label(2),102)rm2(itel)
         else
          write(label(1),103)rm1(itel)
          write(label(2),103)rm2(itel)
          write(label(4),103)rmc(itel)
         endif
         write(label(3),102)pb(itel)
102      format(f7.1)
103      format(f5.1)
         do k=1,klabel
           call pgtext(xtl(k),yshift,label(k))
         enddo
!         call pgtext(xaa,yshift,yaa(itel))
       enddo
c plot size bar
       xlen=ilen*1.
       xpl(2)=xlen/2.
       xpl(1)=-xpl(2)
       yshift=yshift+hei(ktel)+2.*ymargin
       ypl(1)=yshift
       ypl(2)=ypl(1)
       call pgline(2,xpl,ypl)
       write(text,101)ilen
101    format(i5,'R\d\(2281)')
       call pgtext(xpl(2),ypl(2),text)
c and axis of rotation
       xpl(1)=0.
       xpl(2)=0.
       ypl(1)=0.
       ypl(2)=yshift+ymargin
       call pgsls(4)
       call pgline(2,xpl,ypl)
       call pgsls(1)
c add texts, if necessary
123       write(6,*)'give position (x,y) of text'
       write(6,*)'x=0. means: no text to be added'
       read(5,*)xt,yt
       if(xt.ne.0.0)then
        write(6,*)'give text string'
        read(5,99)text
        call pgtext(xt,yt,text)
        goto 123
       endif
99     format(a)
       call pgend
       end

       function rlimit(x,f,df)
c calculates outer limit of roche-lobe
       common/roche/ q,q11,const,const2,xsq,onexsq
       r1=abs(x)
       r2=abs(1.-x)
       r3=abs(x-q11)
       f=q/r1+1./r2+0.5*(1.+q)*r3**2-const
       df=-q*x/r1**3+(1.-x)/r2**3+(1.+q)*(x-q11)
       return
       end
       function rline(y,f,df)
c calculates value of y^2 for x^2 value
       common/roche/ q,q11,const,const2,xsq,onexsq
       r1=sqrt(xsq+y)
       r2=sqrt(onexsq+y)
       f=q/r1+1./r2+const2
       df=-0.5*q/r1**3-0.5/r2**3
       return
       end
       SUBROUTINE CIRKEL(XC,YC,RAD,N)                                   
       DIMENSION X(200),Y(200)                                          
       STEP=6.2831852/(N-1)                                             
       DO 1 I=1,N                                                       
       PHI=I*STEP                                                       
       X(I)=XC+RAD*COS(PHI)                                             
       Y(I)=YC+RAD*SIN(PHI)                                             
1      CONTINUE                                                         
       call pgpoly(n,x,y)
        RETURN                                                          
       END                                                              
       SUBROUTINE DISK(XC,YC,rad,RLEN)                                 
C DRAWS DISK CENTERED ON XC,YC between rad and rlen
       REAL X(5),Y(5)                                                  
       X(1)=XC+rad
       X(2)=X(1)                                                        
       X(3)=XC+RLEN                                                     
       X(4)=X(3)                                                        
       X(5)=X(1)                                                        
       Y(1)=YC+0.15*rad
       Y(2)=YC-0.15*rad                                                
       Y(3)=YC-0.15*rlen                                                     
       Y(4)=YC+0.15*rlen                                                
       Y(5)=Y(1)                                                        
       CALL pgPOLY(5,X,Y)                                            
       X(1)=XC-rad
       X(2)=X(1)                                                        
       X(3)=XC-RLEN                                                     
       X(4)=X(3)                                                        
       X(5)=X(1)                                                        
       Y(5)=Y(1)                                                        
       CALL pgPOLY(5,X,Y)                                            
       RETURN                                                          
       END                                                             
       function rtsafe(funcd,x1,x2,xacc)
c numerical recipes, p.258
       parameter(maxit=100)
       call funcd(x1,fl,df)
       call funcd(x2,fh,df)
       if(fl*fh.ge.0.)pause 'root must be bracketed'
       if(fl.lt.0.)then
        xl=x1
        xh=x2
       else
        xh=x1
        xl=x2
        swap=fl
        fl=fh
        fh=swap
       endif
       rtsafe=0.5*(x1+x2)
       dxold=abs(x2-x1)
       dx=dxold
       call funcd(rtsafe,f,df)
       do j=1,maxit
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.
     c    .or. abs(2.*f).gt.abs(dxold*df) ) then
           dxold=dx
           dx=0.5*(xh-xl)
           rtsafe=xl+dx
           if(xl.eq.rtsafe)return
        else
           dxold=dx
           dx=f/df
           temp=rtsafe
           rtsafe=rtsafe-dx
           if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc)return
        call funcd(rtsafe,f,df)
        if(f.lt.0.)then
          xl=rtsafe
          fl=f
        else
          xh=rtsafe
          fh=f
        endif
       enddo
       pause 'rtsafe exceeds maximum iterations'
       return
       end
