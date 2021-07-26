Program MAIN
!***********************************************************************************************! 
! THIS PROGRAM SOLVER THE EULER TRANSPORT EQUATIONS FOR A COMPRESSIBLE FLOW IN 2D 
! USING ROE METHOD TO SOLVER THE HYPERBOLIC EQUATIONS
!***********************************************************************************************! 
	Implicit real*8(a-h, o-z)
	Real*8, Allocatable :: x(:,:), y(:,:), CL(:), CD(:), CofP(:,:), Ma(:,:)
	Real*8, Allocatable :: snx(:,:,:), sny(:,:,:), s(:,:,:), vol(:,:)
	Real*8, Allocatable :: r(:,:), u(:,:), v(:,:), e(:,:), p(:,:)
	Real*8, Allocatable :: rave(:,:), uave(:,:), vave(:,:), eave(:,:), pave(:,:)
	Real*8, Allocatable :: f(:,:,:), dt(:,:)

	integer, parameter :: dpp=8
	Real(kind=dpp) :: rtot,ptot,rspecific,ttot,aid,acr,radm,aux,tin,pin,rin,ein,cp,cv,Mach,q, eo

	print *,"Start"
	!*********************************** Grid **********************************************!
	Open(1,file='MESH.dat')
	Read(1,*) ni, nj ! <-------------- number of total i and j index

	Allocate(x(ni,nj), y(ni,nj), CL(ni), CD(ni), CofP(ni,nj), Ma(ni,nj))
	Allocate(snx(4,ni,nj), sny(4,ni,nj), s(4,ni,nj), vol(ni,nj))
	Allocate(r(ni,nj), u(ni,nj), v(ni,nj), e(ni,nj), p(ni,nj))
	Allocate(rave(ni,nj), uave(ni,nj), vave(ni,nj), eave(ni,nj), pave(ni,nj))
	Allocate(f(4,ni,nj), dt(ni,nj))

	Do j=1,nj
		Do i=1,ni
			Read(1,*) x(i,j), y(i,j)
		Enddo
	Enddo
	Close(1)	

	Open(1,file='spec.dat')
	Read(1,*) u0
	Read(1,*) v0
	Read(1,*) p0
	Read(1,*) T0
	Read(1,*) Cp
	Read(1,*) Cv
	Read(1,*) niter
	Close(1)

	gamma = Cp/Cv
	Rspecific = Cp - Cv
	Ptot = p0
	Ttot = T0
	rtot = Ptot/(Rspecific*Ttot)
	eo = Cv*Ttot + 0.5*(u0**2 + v0**2)
	q = 0.5*rtot*(u0**2.0 + v0**2.0)
	aid = sqrt(gamma*Ptot/rtot)
	uin = u0/aid
	vin = v0/aid
	v2 = uin**2 + vin**2
	acr = sqrt(2*gamma*(gamma-1)/(gamma+1)*Cv*Ttot)
	radm = Rspecific*Ttot/acr**2
	aux = (1-(gamma-1)/(gamma+1)*(1+tan(alpha)**2)*(uin/acr)**2)
	alpha = atan(v0/u0) 
	tin = (Ttot*aux)/Ttot
	Mach = sqrt(v2/tin)
	pin = Ptot/(rtot*(aid**2))
	!pin = (Ptot*aux**(gamma/(gamma-1)))/(rtot*acr**2)
	rin = (rtot*aux**(1/(gamma-1)))/rtot
	ein = eo/(aid**2)	
	!ein = (Ptot/(gamma-1) + 0.5D0*rtot*(uin**2 + vin**2))/(rtot*acr**2)

	Open(10,file='initial condition.dat')
	Write(10,*)"rhoin = ", rin
	Write(10,*)"uin = ", uin
	Write(10,*)"vin = ", vin
	Write(10,*)"Mach = ", Mach
	Write(10,*)"ein = ", ein
	Write(10,*)"Tin = ", tin
	Write(10,*)"pin = ", pin
	Write(10,*)"aid = ", aid
	Write(10,*)"acr = ", acr
	Write(10,*)"AoA = ", alpha
	Close(10)

	Print*,"rin = ", rin
	Print*,"uin = ", uin
	Print*,"vin = ", vin
	Print*,"ein = ", ein
	Print*,"tin = ", tin
	Print*,"pin = ", pin
	Print*,"aid = ", aid
	Print*,"acr = ", acr
	Print*,"AoA = ", alpha

	!************************************************** Volume ****************************************************************!
	Do j=1,nj-1
	   Do i=1,ni-1
	      s1x = y(i+1,j) - y(i,j)
	      s1y = x(i,j) - x(i+1,j)
	      s2x = y(i+1,j+1) - y(i+1,j)
	      s2y = x(i+1,j) - x(i+1,j+1)
	      s3x = y(i,j+1) - y(i+1,j+1)
	      s3y = x(i+1,j+1) - x(i,j+1)
	      s4x = y(i,j) - y(i,j+1)
	      s4y = x(i,j+1) - x(i,j)
	      s(1,i,j) = sqrt(s1x*s1x + s1y*s1y)
	      s(2,i,j) = sqrt(s2x*s2x + s2y*s2y)
	      s(3,i,j) = sqrt(s3x*s3x + s3y*s3y)
	      s(4,i,j) = sqrt(s4x*s4x + s4y*s4y)
	      snx(1,i,j) = s1x/s(1,i,j)
	      sny(1,i,j) = s1y/s(1,i,j)
	      snx(2,i,j) = s2x/s(2,i,j)
	      sny(2,i,j) = s2y/s(2,i,j)
	      snx(3,i,j) = s3x/s(3,i,j)
	      sny(3,i,j) = s3y/s(3,i,j)
	      snx(4,i,j) = s4x/s(4,i,j)
	      sny(4,i,j) = s4y/s(4,i,j)
	      v11 = (x(i,j) - x(i+1,j+1))*(y(i+1,j) - y(i,j+1))
	      v12 = (x(i,j+1) - x(i+1,j))*(y(i,j) - y(i+1,j+1))
	      vol(i,j) = 0.5*(v11 + v12)
	   Enddo
	Enddo

	!********************************************************* Initial condition **********************************************!
	Do j=1,nj-1
		Do i=1,ni-1
			r(i,j) = rin
			u(i,j) = uin
			v(i,j) = vin
			p(i,j) = pin
			e(i,j) = ein
		Enddo
	Enddo

	!********************************************************** Convective Flux ************************************************!
	print *,"Processing..."
	Open(10,file='Error.dat')
	Do iter=1,niter
		Do j=1,nj-1
			Do i=1,ni-1
				flux1 = 0.0 ! Continuity
				flux2 = 0.0 ! X-Momentum
				flux3 = 0.0 ! Y-Momentum
				flux4 = 0.0 ! Energy
				Do k=1,4
					rl = r(i,j)
					ul = snx(k,i,j)*u(i,j) + sny(k,i,j)*v(i,j)
					vl = -sny(k,i,j)*u(i,j) + snx(k,i,j)*v(i,j)
					pl = p(i,j)
					el = e(i,j)
					hl = (el + pl)/rl
					!***************************** 1st face *****************************************************!
					If(k==1) Then
						If(j==1) Then
							rr = rl
							ur = 0.0
							vr = 0.0
							pr = pl
							er = el
							hr = hl
						Else
							rr = r(i,j-1)
							ur = snx(k,i,j)*u(i,j-1) + sny(k,i,j)*v(i,j-1)
							vr = -sny(k,i,j)*u(i,j-1) + snx(k,i,j)*v(i,j-1)
							pr = p(i,j-1)
							er = e(i,j-1)
							hr = (er + pr)/rr
						Endif
					Elseif(k==2) Then
						If(i==ni-1) Then
							rr = r(1,j)
							ur = snx(k,i,j)*u(1,j) + sny(k,i,j)*v(1,j)
							vr = -sny(k,i,j)*u(1,j) + snx(k,i,j)*v(1,j)
							pr = p(1,j)
							er = e(1,j)
						Else
							rr = r(i+1,j)
							ur = snx(k,i,j)*u(i+1,j) + sny(k,i,j)*v(i+1,j)
							vr = -sny(k,i,j)*u(i+1,j) + snx(k,i,j)*v(i+1,j)
							pr = p(i+1,j)
							er = e(i+1,j)
						Endif
						hr = (er + pr)/rr
					Elseif(k==3) Then
						If(j==nj-1) Then
							rr = rl
							ur = ul
							vr = vl
							pr = pl
							er = el
							hr = hl
						Else
							rr = r(i,j+1)
							ur = snx(k,i,j)*u(i,j+1) + sny(k,i,j)*v(i,j+1)
							vr = -sny(k,i,j)*u(i,j+1) + snx(k,i,j)*v(i,j+1)
							pr = p(i,j+1)
							er = e(i,j+1)
							hr = (er + pr)/rr
						Endif
					Elseif(k==4) Then
						If(i==1) Then
							rr = r(ni-1,j)
							ur = snx(k,i,j)*u(ni-1,j) + sny(k,i,j)*v(ni-1,j)
							vr = -sny(k,i,j)*u(ni-1,j) + snx(k,i,j)*v(ni-1,j)
							pr = p(ni-1,j)
							er = e(ni-1,j)
						Else
							rr = r(i-1,j)
							ur = snx(k,i,j)*u(i-1,j) + sny(k,i,j)*v(i-1,j)
							vr = -sny(k,i,j)*u(i-1,j) + snx(k,i,j)*v(i-1,j)
							pr = p(i-1,j)
							er = e(i-1,j)
						Endif
						hr = (er + pr)/rr
					Endif
					rbar = sqrt(rl*rr)
					ubar = (sqrt(rl)*ul + sqrt(rr)*ur)/(sqrt(rl) + sqrt(rr))
					vbar = (sqrt(rl)*vl + sqrt(rr)*vr)/(sqrt(rl) + sqrt(rr))
					hbar = (sqrt(rl)*hl + sqrt(rr)*hr)/(sqrt(rl) + sqrt(rr))
					v2 = ubar*ubar + vbar*vbar
					abar = sqrt((gamma-1.0)*(hbar-0.5*v2))
					e1 = ubar - abar
					e2 = ubar
					e3 = e2
					e4 = ubar + abar
					ev11 = 1.0
					ev12 = e1
					ev13 = vbar
					ev14 = hbar - ubar*abar
					ev21 = 1.0
					ev22 = e2
					ev23 = vbar
					ev24 = 0.5*v2
					ev31 = 0.0
					ev32 = 0.0
					ev33 = 1.0
					ev34 = vbar
					ev41 = 1.0
					ev42 = e4
					ev43 = vbar
					ev44 = hbar + ubar*abar
					du1 = rr - rl
					du2 = rr*ur - rl*ul
					du3 = rr*vr - rl*vl
					du4 = er - el
					dp = pr - pl
					du = ur - ul
					dv = vr - vl
					dr = rr - rl
					!du4bar = du4 - vbar*(du3 - vbar*du1)
					c1 = 0.5*(dp - rbar*abar*du)/(abar*abar)
					c2 = dr - dp/(abar*abar)
					c3 = rbar*dv
					c4 = 0.5*(dp + rbar*abar*du)/(abar*abar)
					!c3 = du3 - vbar*du1
					!c2 = (gamma-1.0)*(du1*(hbar - ubar*ubar) + ubar*du2 - du4bar)/(abar*abar)
					!c1 = 0.5*(du1*(ubar + abar) - du2-abar*c2)/abar
					!c4 = du1 - (c1 + c2)
					f1s = c1*abs(e1)*ev11 + c2*abs(e2)*ev21 + c3*abs(e3)*ev31 + c4*abs(e4)*ev41
					f2s = c1*abs(e1)*ev12 + c2*abs(e2)*ev22 + c3*abs(e3)*ev32 + c4*abs(e4)*ev42
					f3s = c1*abs(e1)*ev13 + c2*abs(e2)*ev23 + c3*abs(e3)*ev33 + c4*abs(e4)*ev43
					f4s = c1*abs(e1)*ev14 + c2*abs(e2)*ev24 + c3*abs(e3)*ev34 + c4*abs(e4)*ev44
					f1l = rl*ul
					f2l = rl*ul*ul + pl
					f3l = rl*vl*ul
					f4l = ul*(el + pl)
					f1r = rr*ur
					f2r = rr*ur*ur + pr
					f3r = rr*vr*ur
					f4r = ur*(er + pr)
					f1 = 0.5*(f1l + f1r) - 0.5*f1s
					f2 = 0.5*(f2l + f2r) - 0.5*f2s
					f3 = 0.5*(f3l + f3r) - 0.5*f3s
					f4 = 0.5*(f4l + f4r) - 0.5*f4s
					flux1 = flux1 + f1*s(k,i,j)
					flux2 = flux2 + (snx(k,i,j)*f2 - sny(k,i,j)*f3)*s(k,i,j)
					flux3 = flux3 + (sny(k,i,j)*f2 + snx(k,i,j)*f3)*s(k,i,j)
					flux4 = flux4 + f4*s(k,i,j)
				Enddo
				f(1,i,j) = flux1
				f(2,i,j) = flux2
				f(3,i,j) = flux3
				f(4,i,j) = flux4
			Enddo
		Enddo

		!****************************************************************** Time step *************************************************************!
		sigma = 0.05 !<----------- SIGMA
		Do j=1,nj-1
			Do i=1,ni-1
				dsnxi = 0.5*(snx(2,i,j) - snx(4,i,j))
				dsnyi = 0.5*(sny(2,i,j) - sny(4,i,j))
				dsnxj = 0.5*(snx(3,i,j) - snx(1,i,j))
				dsnyj = 0.5*(sny(3,i,j) - sny(1,i,j))
				dsi = 0.5*(s(2,i,j) + s(4,i,j))
				dsj = 0.5*(s(1,i,j) + s(3,i,j))
				c = sqrt(gamma*p(i,j)/r(i,j))
			G1 = (abs(u(i,j)*dsnxi + v(i,j)*dsnyi) + c)*dsi
			G2 = (abs(u(i,j)*dsnxj + v(i,j)*dsnyj) + c)*dsj
			!G12 = sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j)) + c
			dt(i,j) = sigma*vol(i,j)/(G1 + G2)
			!dt(i,j) = sigma*vol(i,j)/(G12)
			time = dt(i,j) 
			Enddo
		Enddo

		resr = 0.0
		resu = 0.0
		resv = 0.0
		rese = 0.0
		resL2norm = 0.0

		!******************************************************************** Solver ****************************************************************!
		Do j=1,nj-1
			Do i=1,ni-1
				rnew = r(i,j) - dt(i,j)*f(1,i,j)/vol(i,j)
				runew = r(i,j)*u(i,j) - dt(i,j)*f(2,i,j)/vol(i,j)
				rvnew = r(i,j)*v(i,j) - dt(i,j)*f(3,i,j)/vol(i,j)
				enew = e(i,j) - dt(i,j)*f(4,i,j)/vol(i,j)
				dr = (rnew - r(i,j))*(rnew - r(i,j))
				dru = (runew - r(i,j)*u(i,j))*(runew - r(i,j)*u(i,j))
				drv = (rvnew - r(i,j)*v(i,j))*(rvnew - r(i,j)*v(i,j))
				de = (enew - e(i,j))*(enew - e(i,j))
				resL2norm = resL2norm + dr + dru + drv + de
				!resr = max(resr,dr)
				resr = resr + dr
				!resu = max(resu,du)
				resu = resu + dru
				!resv = max(resv,dv)
				resv = resv + drv
				!rese = max(rese,de)
				rese = rese + de 
				r(i,j) = rnew
				u(i,j) = runew/rnew
				v(i,j) = rvnew/rnew
				e(i,j) = enew
				p(i,j) = (enew - 0.5*rnew*(u(i,j)*u(i,j) + v(i,j)*v(i,j)))*(gamma-1.0)
			Enddo
		Enddo
		!********************************************************************* Residual Errors *******************************************************!
		res  = sqrt(resL2norm/((nj - 1)*(ni - 1)*4))
		resR = sqrt(resr/((nj-1)*(ni-1)))
		resU = sqrt(resu/((nj-1)*(ni-1)))
		resV = sqrt(resv/((nj-1)*(ni-1)))
		resE = sqrt(rese/((nj-1)*(ni-1)))
		Write(10,*)iter, resR, resU, resV, resE 
		Print*,"Iteration:",iter,"Residual=",res
		Print*,"==============================================================================="
	Enddo
	Close(10)

	!***************************************************************************** Post ******************************************************************!
	Do j=2,nj-1
		Do i=1,ni-1
			!************************************************************ Interface **************************************************************!
			If(i==1) Then
				rave(1,j) = 0.25*(r(1,j) + r(ni-1,j) + r(1,j-1) + r(ni-1,j-1))
				uave(1,j) = 0.25*(u(1,j) + u(ni-1,j) + u(1,j-1) + u(ni-1,j-1))
				vave(1,j) = 0.25*(v(1,j) + v(ni-1,j) + v(1,j-1) + v(ni-1,j-1))
				eave(1,j) = 0.25*(e(1,j) + e(ni-1,j) + e(1,j-1) + e(ni-1,j-1))
				rave(ni,j) = rave(1,j)
				uave(ni,j) = uave(1,j)
				vave(ni,j) = vave(1,j)
				eave(ni,j) = eave(1,j)
			Else
			!************************************************************ Domain *****************************************************************!
				rave(i,j) = 0.25*(r(i,j) + r(i-1,j) + r(i,j-1) + r(i-1,j-1))
				uave(i,j) = 0.25*(u(i,j) + u(i-1,j) + u(i,j-1) + u(i-1,j-1))
				vave(i,j) = 0.25*(v(i,j) + v(i-1,j) + v(i,j-1) + v(i-1,j-1))
				eave(i,j) = 0.25*(e(i,j) + e(i-1,j) + e(i,j-1) + e(i-1,j-1))
			Endif
		Enddo
	Enddo

	Do i=1,ni-1
		If(i==1) Then
		!******************************************************************** Interface **************************************************************!
			rave(1,1) = 0.5*(r(1,1) + r(ni-1,1))
			uave(1,1) = 0.5*(u(1,1) + u(ni-1,1))
			vave(1,1) = 0.5*(v(1,1) + v(ni-1,1))
			eave(1,1) = 0.5*(e(1,1) + e(ni-1,1))
			rave(ni,1) = rave(1,1)
			uave(ni,1) = uave(1,1)
			vave(ni,1) = vave(1,1)
			eave(ni,1) = eave(1,1)
		Else
		!********************************************************************* Wall ******************************************************************! 
			rave(i,1) = 0.5*(r(i,1) + r(i-1,1))
			uave(i,1) = 0.5*(u(i,1) + u(i-1,1))
			vave(i,1) = 0.5*(v(i,1) + v(i-1,1))
			eave(i,1) = 0.5*(e(i,1) + e(i-1,1))
		Endif
	Enddo

	!*************************************************************************** Far field ***************************************************************! 
	Do i=1,ni
		rave(i,nj) = rave(i,nj-1)
		uave(i,nj) = uave(i,nj-1)
		vave(i,nj) = vave(i,nj-1)
		eave(i,nj) = eave(i,nj-1)
	Enddo

	!******************************************************************** Cp, Forces, Cl and Cd  *********************************************************!
	Do j=1,nj
	   Do i=1,ni
	   	pave(i,j) = (eave(i,j) - 0.5*rave(i,j)*(uave(i,j)*uave(i,j) + vave(i,j)*vave(i,j)))*(gamma-1.0)
		CofP(i,j) = 2.0*(pave(i,j) - 1.0/gamma)/(Mach**2)
		Ma(i,j) = sqrt((uave(i,j)**2 + vave(i,j)**2)/tin)
	   Enddo
	Enddo

	Do i=1,ni-1
		press = 0.5*(pave(i,1) + pave(i+1,1))
		dfy = press*snx(1,i,1)*s(1,i,1)
		dfx = press*sny(1,i,1)*s(1,i,1)
		Fy = Fy + dfy
		Fx = Fx + dfx
		L = Fy*cos(alpha) - Fx*sin(alpha)
		D = Fx*cos(alpha) - Fy*sin(alpha)
		CL(i) = L/q
		CD(i) = D/q
	Enddo

	!************************************************************************** Date ********************************************************************!
	Open(1,file='cp_plot.dat')
	Do i=1,ni
	      Write(1,*) x(i,1), 2.0*(pave(i,1) - 1.0/gamma)/(Mach**2)
	Enddo

	Open(1,file='cl_plot.dat')
	Do i=1, ni
		Write(1,*) x(i,1), CL(i), CD(i)
	Enddo

	Open(1,file='solver.tec')
	Write(1,*) 'title="otype"'
	Write(1,*) 'variables= x, y, density, u, v, p, Cp, Mach, S'
	Write(1,*) 'zone t="grid_A", i=',ni,' j=',nj,'datapacking=point'
	Do j=1,nj
		Do i=1,ni
			Write(1,*) x(i,j), y(i,j), rave(i,j)*rtot, uave(i,j)*aid, vave(i,j)*aid, pave(i,j)*rtot*(aid**2), CofP(i,j), Ma(i,j), pave(i,j)*rtot*(aid**2)/((rave(i,j)*rtot)**gamma)
		Enddo
	Enddo
	Close(1)
	print *,"Finish"
	Stop
End Program Main
