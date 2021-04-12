c-----------------------------------------------------------------------
c..TSD SOR/SLOR SOLVER for NACA profiles                               |
c  Course     : AE443                                                  |
c  Instructor : Dr. Ismail H. TUNCER                                   |
c-----------------------------------------------------------------------
      program TSD
      parameter (imx=301,jmx=301)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/naca/naca,nair,dydx(2,imx)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)
      common/coef/Cr(jmx),Cs(jmx),Ct(jmx),Cv(jmx),Cw(jmx),Rhs(jmx)
      data itmx/10000/  io/1/

c..Generate the grid and initialize the solution
      call INIT
c..Start iterative solution loop
      DO iter = 1,itmx
c..Evaluate the mu(i,j) term
         call EVALMU
c..Apply explicit farfield BC
         call BCFAR
c..Evaluate circulation for lifting flows
         circ = phi(ite,jau)-phi(ite,jal)
c..Compute Dphi and update phi by JACOBI/SOR/SLOR
         call SOLVER(iter)
c..Output intermediate cp distribution and flow solution
         if (mod(iter,io).eq.0 ) call OUT(iter)
      ENDDO
      call OUT(iter)

      STOP
      END

c--------------------------------------------------------------------------
      subroutine INIT
      parameter (imx=301,jmx=301)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/naca/naca,nair,dydx(2,imx)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)

c..Read input parameters
      print*,' Enter NACA Profile, # of Surface Points, Mach# and AoA :'
      read(5,*) naca,nair,fsmach,alpha
c     naca    = 3
c     nair    = 51
c     fsmach  = 0.3
c     alpha   = 0.

      pi=4*atan(1.0)
      gamma   = 1.4
      gmp1    = gamma + 1.
      fsmach2 = fsmach**2

c..Set up the cartesian grid
      dx1 = 1./(nair-1)
      dy1 = 0.02
      call MAKEGRID(dx1,1.10,dy1,1.10, 15.)

c..Get the airfoil upper/lower surface slopes; dydx
      call NACAFOIL(naca,dydx)

c..Rotate dydx for alpha
      alpha=alpha*pi/180.
      do i=ile,ite,1
         dydx(1,i)=tan(atan(dydx(1,i))-alpha)
         dydx(2,i)=tan(atan(dydx(2,i))-alpha)
      end do
      
      open(41,file='rotated.dat')
      write(41,'(i3,3f8.4)')(i,xg(i),(dydx(k,i),k=1,2),i=ile,ite)
      close(41)
      
c..Initialize the solution arrays
      do j = 1,jg
      do i = 1,ig
         mu(i,j)   = 0
         phi(i,j)  = 0.
         dphi(i,j) = 0.
      enddo
      enddo

      return
      end

c--------------------------------------------------------------------------
      subroutine MAKEGRID(dx1,xratio,dy1,yratio,farbc)
      parameter (imx=301,jmx=301)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/naca/naca,nair,dydx(2,imx)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      real xstretch(imx), ystretch(jmx)
      character filename*32

c..Construct the stretching along x and y directions
      xstretch(1) = 0.
      is = 1
      dxi=dx1
      do while (xstretch(is) .lt. farbc )
         is = is+1
         xstretch(is) = xstretch(is-1) + dxi
         dxi = dxi*xratio
      enddo
      is = is -1
      ile = is
      ite = ile + nair - 1
      ig  = ite + is - 1

      ystretch(1) = 0.
      js = 1
      dyj=dy1
      do while (ystretch(js) .lt. farbc )
         js = js+1
         ystretch(js) = ystretch(js-1) + dyj
         dyj = dyj*yratio
      enddo
      js  = js - 1
      jal = js
      jau = jal + 1
      jg  = 2*js

      print*,  ' GRID:',ig,' x',jg,' - ile,ite,jal:',ile,ite,jal

c.....Construct x array
      do i=ile,1,-1
         xg(i) = xg(ile) - xstretch(ile+1-i)
      enddo

      do i=ile+1,ite
         xg(i) = xg(i-1) + dx1
      enddo

      xg(ite) = 1.0

      do i=ite+1,ig
         xg(i) = xg(ite) + xstretch(i+1-ite)
      enddo

c.....Construct y array
      yg(jal) = -0.01
      yg(jau) = 0.01

      do j=jal-1,1,-1
         yg(j) = yg(jal) - ystretch(jal+1-j)
      enddo

      do j=jau+1,jg
         yg(j) = yg(jau) + ystretch(j+1-jau)
      enddo

c.....Construct dx array
      do i=2,ig
         dx(i) = xg(i)-xg(i-1)
      enddo
c.....Construct dy array
      do j=2,jg
         dy(j) = yg(j)-yg(j-1)
      enddo
      
      dx(1) = dx(2)
      dy(1) = dy(2)

c..Output the grid for visualization
      filename='grid.tec'
      open(1,file=filename,form='formatted')
c..TECPLOT format
      write(1,*) ' variables="x","y"'
      write(1,*) ' zone i=',ig, ',j=',jg
      write(1,'(2e12.4)') ((xg(i),yg(j),i=1,ig),j=1,jg)
      close(1)
      
      return
      end
c--------------------------------------------------------------------------
      subroutine SOLVER(iter)
      parameter (imx=301,jmx=301)
      common/naca/naca,nair,dydx(2,imx)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)
      common/coef/Cr(jmx),Cs(jmx),Ct(jmx),Cv(jmx),Cw(jmx),Rhs(jmx)
      common/flow/fsmach,fsmach2,gamma,gmp1
      data omega/0.3/ resl21/0./
      save omega,resl21
      real Res3(jmx)

      DO i=2,ig-1
      Do j=2,jg-1
c..Evaluate the coefficients; Cr,Cs,Ct,Cv,Cw and Res
         call COEFF(i,j)
c..delta_phi for Point-Jacobi
c         dphi(i,j) = Rhs(j)/Ct(j)
c..delta_phi for SOR
c         Res   = Rhs(j)  - dphi(i-1,j)*Cs(j) - dphi(i-2,j)*Cr(j) -
c     &   dphi(i,j-1)*Cv(j)
c         dphi(i,j) = omega*Res/Ct(j)
c..residual for SLOR
          Res3(j) = Rhs(j)  - dphi(i-1,j)*Cs(j) - dphi(i-2,j)*Cr(j)

      ENDDO

c..In SLOR solve for delta_phi implicity along the i=constant lines
      ct(j) = ct(j)/omega
      call THOMAS(2,jg-1,cv,ct,cw,res3,jg-2)
c..Extract dphi(i,j) from the solution
      do j=2,jg-1
      dphi(i,j) = res3(j)
      enddo

      ENDDO

c.. Update the solution and evaluate the Max/L2 norm of the residual
      resmx = 0.0
      res2  = 0.0
      do i=2,ig-1
      do j=2,jg-1
         phi(i,j) = phi(i,j)+dphi(i,j)
         ares = abs(dphi(i,j))
         if(ares.gt.resmx) then
           resmx  = ares
           iresmx = i
           jresmx = j
         endif
         res2 = res2 + ares**2
      enddo
      enddo
      resl2 = sqrt(res2/((ig-2)*(jg-2)))
      if(resl21 .eq. 0.) resl21=resl2
      resl2log = alog10(resl2/resl21)
      print *, iter, resl2log, resmx,iresmx,jresmx

      return
      end

c--------------------------------------------------------------------------
      subroutine  COEFF(i,j)
c..Evaluate the matrix coefficients
      parameter (imx=301,jmx=301)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/naca/naca,nair,dydx(2,imx)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)
      common/coef/Cr(jmx),Cs(jmx),Ct(jmx),Cv(jmx),Cw(jmx),Rhs(jmx)

c....apply the Wall BC on the Ctq term
      if ( j .eq. jal .and. i .ge .ile .and. i .le .ite ) then
         Ctq = (1./dy(j))/(dy(j+1)+dy(j))
      elseif (j .eq. jau .and. i .ge .ile .and. i .le .ite ) then
         Ctq = (1./dy(j+1))/(dy(j+1)+dy(j))
      else
         Ctq = (1./dy(j+1)+1./dy(j))/(dy(j+1)+dy(j))
      endif

c.....Eliminate Cr term for i=2
      if (i.eq.2) then
      Cr(j) = 0.
      else
      Cr(j)  = (mu(i-1,j)*An(i-1,j)/dx(i-1))/(dx(i+1)+dx(i))
      endif

      if (i.eq.2) then
      Cs(j) = 0.
      else
      Cs(j)  = ((1-mu(i,j)-mu(i-1,j))*An(i,j)/dx(i) -
     +       mu(i-1,j)*An(i-1,j)/dx(i-1)) / (dx(i+1)+dx(i))
      endif
     
      Ct(j)  = -( (1-mu(i,j))*An(i+1,j)/dx(i+1)
     +           +(1-mu(i,j)-mu(i-1,j))*An(i,j)/dx(i))/(dx(i+1)+dx(i))
     +         - Ctq

c.....Eliminate j-1 derivative on upper surface
      if (j.eq.jau .and. i.ge.ile .and. i.le.ite) then
      Cv(j) = 0.
      else
      Cv(j)  = (1./dy(j))*(1./(dy(j+1)+dy(j)))
      endif

c.....Eliminate j+1 derivative on lower surface
      if (j.eq.jal .and. i.ge.ite .and. i.le.ite) then
      Cw(j) = 0.
      else
      Cw(j)  = (1./dy(j+1))*(1./(dy(j+1)+dy(j)))
      endif
      
      Rhs(j) = -(         ( (1-mu(i,j))*Pn(i+1,j)
     +           -(1-mu(i,j)-mu(i-1,j))*Pn(i,j)
     +           -            mu(i-1,j)*Pn(i-1,j) ) / (dx(i)+dx(i+1))
     +           +      (  Qn(i,j+1,j) - Qn(i,j,j)) / (dy(j)+dy(j+1)))


      return
      end

c--------------------------------------------------------------------------
      subroutine OUT(nio)
      parameter (imx=301,jmx=301)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)
      common/naca/naca,nair,dydx(2,imx)
      common/flow/fsmach,fsmach2,gamma,gmp1
      real cp(2,imx), u(imx,jmx),v(imx,jmx),mach(imx,jmx),t(imx,jmx),
     + rho(imx,jmx), pressure(imx,jmx), a(imx,jmx)
      character chd*10,cstep*5,filename*32,string*8
      data chd/'0123456789'/
      data Tinf/288.15/
      real uflow(imx,jmx), vflow(imx,jmx)

c..Evaluate u, v, cp, Mach no, Pressure, etc. first

      vinf = fsmach*sqrt(gamma*287*Tinf)

c..calculate the u and v
      do i=2,ig-1
      do j=2,jg-1

        u(i,j)=(phi(i+1,j)-phi(i-1,j))/(dx(i+1)+dx(i))

	 if(i.ge.ile.and.i.le.ite) then
		if(j.eq.jal) then
		v(i,j) = Qn(i,j+1,j)
		else
		v(i,j) = Qn(i,j,j)
		endif
	 else
	    v(i,j) = Qn(i,j,j)
	 endif

      end do
      end do

c..calculate temperature/density/pressure ratio
c..local speed of sound and local mach number
      do i=1,ig
      do j=1,jg
         uflow(i,j) = (1.0 + u(i,j))*vinf
         vflow(i,j) = v(i,j)*vinf
         vell = sqrt(uflow(i,j)**2+vflow(i,j)**2)
         t(i,j) = 1.0 + (gamma-1.0)/2.0*fsmach2*(1-vell**2/vinf**2)
         rho(i,j) = t(i,j)**(1.0/(gamma-1.0))
         pressure(i,j) = t(i,j)**(gamma/(gamma-1.0))
         a(i,j) = sqrt(gamma*287*t(i,j)*Tinf)
         mach(i,j) = vell/a(i,j)
      enddo
      enddo
      

c..calculate the pressure coefficient on both airfoil surfaces
      do i=ile,ite,1
          cp(1,i) = 2.0/(gamma*fsmach2)*(pressure(i,jal)-1.0)
          cp(2,i) = 2.0/(gamma*fsmach2)*(pressure(i,jau)-1.0)
      end do

c..Construct the filenames and write out the data
      write(string,'(f8.5)') float(nio)/100000
      read(string,'(3x,a5)') cstep
      filename = 'cp-'//cstep//'.dat'
      open(1,file=filename)
      write(1,'(2E13.5)') (xg(i),cp(1,i), i=ite,ile,-1)
      write(1,'(2E13.5)') (xg(i),cp(2,i), i=ile,ite)
      close(1)

      filename='vars-'//cstep//'.tec'
      open(1,file=filename,form='formatted')
c..Output Phi,u,v,pressure,density,Mach in TECPLOT format
c      write(1,*) ' variables="x","y","phi","mach", "pressure", "u", "v"'
      write(1,*) ' variables="x","y","phi"'
      write(1,*) ' zone i=',ig, ',j=',jg
      do j = 1,jg
      do i = 1,ig
c         write(1,'(8E12.4)') xg(i),yg(j),phi(i,j),mach(i,j),
c     +     pressure(i,j), uflow(i,j), vflow(i,j)
         write(1,'(8E12.4)') xg(i),yg(j),phi(i,j)
      enddo
      enddo

      close(1)
      return
      end

c--------------------------------------------------------------------------
      subroutine BCFAR
      parameter (imx=301,jmx=301)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)
      common/flow/fsmach,fsmach2,gamma,gmp1
      
      pi=4*atan(1.0)
c..First order farfield BC: Dphi/Dx=0, Dphi/Dy=0
      do i=2,ig-1
        phi(i,1)  =  phi(i,2)
        phi(i,jg) =  phi(i,jg-1)
      enddo
      do j=2,jg-1
        phi(1,j)  =  phi(2,j)
        phi(ig,j) =  phi(ig-1,j)
      enddo

c..Apply BC for lifting flows for more accuracy
c..The wake condition is applied in here

c      if ((abs(circ) .gt. 0.0001)) THEN
c      do i=2,ig-1
c         phi(i,1)=circ/(2*pi*(1-fsmach2))*atan(sqrt(1-fsmach2)*
c     >    yg(1)/xg(i))
c         phi(i,jg)=circ/(2*pi*(1-fsmach2))*atan(sqrt(1-fsmach2)*
c     >             yg(jg)/xg(i))
c      end do
c
c      do j=2,jg-1
c         phi(1,j)=circ/(2*pi*(1-fsmach2))*atan(sqrt(1-fsmach2)*
c     >            yg(j)/xg(1))
c         phi(ig,j)=circ/(2*pi*(1-fsmach2))*atan(sqrt(1-fsmach2)*
c     >             yg(j)/xg(ig))
c      end do
c      ENDIF


      return
      end



c--------------------------------------------------------------------------
      subroutine EVALMU
      parameter (imx=301,jmx=301)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)

      do i=1,ig
      do j=1,jg
         if (i.eq.1.or.j.eq.1.or.i.eq.ig.or.j.eq.jg) then
         mu(i,j) = 0.
         else
         phix = (phi(i,j)-phi(i-1,j))/dx(i)
         ss = fsmach2*(1.+gmp1*phix)
            if (ss.gt.1.0) then
            mu(i,j) = 1.
            else
            mu(i,j) = .0
            endif
         endif
      enddo
      enddo

      return
      end

c-------------------------------------------------------------------
      function An(i,j)
      parameter (imx=301,jmx=301)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)


      if (i.ne.1) then
      phix = (phi(i,j)-phi(i-1,j))/dx(i)
      An = (1-fsmach2)-(gmp1*fsmach2*phix)
      else
      An = 1-fsmach2
      endif

      return
      end

c-------------------------------------------------------------------
      function Pn(i,j)
      parameter (imx=301,jmx=301)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)

      if ( i .ne. 1) then
        phix = (phi(i,j)-phi(i-1,j))/dx(i)
        Pn   = (1.-fsmach2)*phix - 0.5*gmp1*fsmach2*phix**2
      else
        Pn   = 0.
      endif

      return
      end

c--------------------------------------------------------------------------
      function Qn(i,j,jp)                 !..apply wall/wake BC's inside
      parameter (imx=301,jmx=301)
      common/flow/fsmach,fsmach2,gamma,gmp1
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      common/sol/ circ,mu(imx,jmx),phi(imx,jmx),dphi(imx,jmx)
      common/naca/naca,nair,dydx(2,imx)
      
      if ((i.ge.ile).and.(i.le.ite)) then
         if((jp.eq.jal).and.(j.gt.jal)) then
            Qn=dydx(1,i)
         elseif ((jp.eq.jau) .and. (jp.eq.j)) then
            Qn=dydx(2,i)
         else
            Qn=(phi(i,j)-phi(i,j-1))/dy(j)
         endif

c..apply the wake condition here

      elseif (i.gt.ite.and.j.eq.jau.and.jp.eq.jal) then
              Qn=(phi(i,jau)-phi(i,jal)-circ)/dy(jau)
      elseif (i.gt.ite.and.j.eq.(jal+1).and.jp.eq.jau) then
             Qn=(phi(i,jau)-phi(i,jal)-circ)/dy(jau)
      else
         Qn=(phi(i,j)-phi(i,j-1))/dy(j)
      end if


      return
      end
      
c-----------------------------------------------------------------------
      subroutine NACAFOIL(naca,dydx)
      parameter (imx=301, jmx=301)
      common/grid/ig,jg,ile,ite,jal,jau,xg(imx),yg(jmx),dx(imx),dy(jmx)
      real dydx(2,imx)
      real xupper(imx),yupper(imx), xlower(imx),ylower(imx)
      pi = acos(-1.)

c..Decompose the NACA number to determine airfoil properties
       if (( naca .gt. 25999 ).or.( naca .lt. 1 )) naca = 12
       ieps = naca / 1000
       iptmax = naca / 100 - 10 * ieps
       itau = naca - 1000 * ieps - 100 * iptmax

c..Set the coefficients.
       epsmax = ieps   * 0.01
       ptmax  = iptmax * 0.1
       tau    = itau   * 0.01

c..Error correction for bogus NACA numbers.
       if (( naca .le. 9999 ) .and. ( epsmax .gt. 0 ) .and.
     >     ( ptmax .eq. 0 )) ptmax = 0.1

c..If NACA 5 digit coding is used, make neccessary changes.
       if ( ieps .ge. 10 ) then
         if ( ieps .eq. 21 ) then
           ptmax = 0.0580
           ak1 = 361.4
         elseif ( ieps .eq. 22 ) then
           ptmax = 0.1260
           ak1 = 51.64
         elseif ( ieps .eq. 23 ) then
           ptmax = 0.2025
           ak1 = 15.957
         elseif ( ieps .eq. 24 ) then
           ptmax = 0.2900
           ak1 = 6.643
         elseif ( ieps .eq. 25 ) then
           ptmax = 0.3910
           ak1 = 3.230
         endif
         epsmax = ak1 * ptmax**3 / 6
       endif

c.....Loop over the lower and upper surfaces
      do i=ite,ile,-1
         xc = xg(i)
         call NACA45(naca,tau,epsmax,ptmax,xc,thick,camber,beta)
         yupper(i) = camber + thick*cos(beta)
         ylower(i) = camber - thick*cos(beta)
      enddo
      
c.....Connect the points
         ylower(ite) = yupper(ite)
         ylower(ile) = yupper(ile)

c.....Calculate the slope dydx(1,imx) and dydx(2,imx)

      do i=ile,ite,1
         if ((i.gt.ile) .and. (i.lt.ite)) then
            dydx(1,i)=(ylower(i+1)-ylower(i-1))/(xg(i+1)-xg(i-1))
            dydx(2,i)=(yupper(i+1)-yupper(i-1))/(xg(i+1)-xg(i-1))
          elseif (i.eq.ile) then
            dydx(1,i)=(ylower(i+1)-ylower(i))/(xg(i+1)-xg(i))
            dydx(2,i)=(yupper(i+1)-yupper(i))/(xg(i+1)-xg(i))
          else
            dydx(1,i)=(ylower(i)-ylower(i-1))/(xg(i)-xg(i-1))
            dydx(2,i)=(yupper(i)-yupper(i-1))/(xg(i)-xg(i-1))
          endif
      enddo

c.....write the the airfoil/panel x and y coordinates
      open(2,file='airfoil.dat')
      write(2,'(e15.8,2x,e15.8,2x,e15.8)') (xg(n),yupper(n),dydx(2,n),
     c n=ite,ile,-1)
      write(2,'(e15.8,2x,e15.8,2x,e15.8)') (xg(n),ylower(n),dydx(1,n),
     c n=ile+1,ite)
      close(2)

      return
      end

c------------------------------------------------------------------------------
       subroutine NACA45(naca,tau,epsmax,ptmax,xc,thick,camber,beta)
c..Compute the thickness, camber, and angular location of an airfoil point.

c..Compute the thickness
       if ( xc .lt. 1.0E-10 ) then
         thick = 0.0    !..Thickness is corrected when xc is very small.
       else
         thick = tau * 5 * ( 0.2969 * SQRT(xc)
     >               - xc * ( 0.1260
     >               + xc * ( 0.3537
     >               - xc * ( 0.2843
     >               - xc * 0.1015))) )
       endif

c..Compute the camber
       if ( epsmax .eq. 0.0 ) then
c..For NACA 4-digit symmetrical arfoils.
         camber = 0.0
         beta = 0.0
       else
         if ( naca .gt. 9999 ) then
c..For NACA 5 digit airfoils
c..Ptmax = m and epsmax = (k_1*m^3)/6 from Abbott and Doenhoff.
           if ( xc .gt. ptmax ) then
             camber = epsmax * ( 1.0 - xc )
             dcamdx = - epsmax
           else
             w = xc / ptmax
             camber = epsmax * ( w**3 - 3*w**2 +(3.-ptmax)*w)
             dcamdx = epsmax/ptmax*(3*w**2 - 6*w + ( 3.0-ptmax))
           endif
         else

c..For NACA 4 digit airfoils.
           if ( xc .gt. ptmax ) then
             camber = epsmax / ( 1.0 - ptmax )**2
     >              * ( 1. + xc - ptmax * 2 ) * ( 1. - xc )
             dcamdx = epsmax * 2 / ( 1.0 - ptmax )**2
     >              * ( ptmax - xc )
           else
             camber = epsmax / ptmax**2 * ( ptmax*2 - xc ) * xc
             dcamdx = epsmax * 2 / ptmax**2  * ( ptmax - xc )
           endif
         endif

         beta = atan(dcamdx)
       endif

       return
       end
c-------------------------------------------------------------------
      SUBROUTINE THOMAS(il,iu, A,B,C,F,neq)
c............................................................
c Solution of a tridiagonal system of n equations of the form
c  A(k)*x(k-1) + B(k)*x(k) + C(k)*x(k+1) = F(i)  for k=il,iu
c  the solution X(i) is stored in F(i)
c  A(il-1) and C(iu+1) are not used.
c  A,B,C,F are arrays to be filled by the caller program
c............................................................
      real a(neq),b(neq),c(neq),f(neq),x(500)

      x(il)=c(il)/b(il)
      f(il)=f(il)/b(il)

      ilp1 = il+1
      do i=ilp1,iu
         z=1./(b(i)-a(i)*x(i-1))
         x(i)=c(i)*z
         f(i)=(f(i)-a(i)*f(i-1))*z
      enddo

      iupil=iu+il
      do ii=ilp1,iu
         i=iupil-ii
         f(i)=f(i)-x(i)*f(i+1)
      enddo

      return
      end
