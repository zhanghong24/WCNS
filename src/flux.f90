!_____________________________________________________________________!
subroutine flux_SW(priml,primr,nl,nx,ny,nz,nt,f,efix,gamal,gamar)
  implicit none
  integer :: m,nl
  real :: priml(nl),primr(nl),f(nl),fl(nl),fr(nl)
  real :: nx,ny,nz,nt,gamal,gamar,efix

  call flux_SW_pn(priml,nl,nx,ny,nz,nt,fl,efix,gamal, 1)
  call flux_SW_pn(primr,nl,nx,ny,nz,nt,fr,efix,gamar,-1)
  do m=1,nl
     f(m) = fl(m)+fr(m)
  enddo

  return
end subroutine flux_SW
!_____________________________________________________________________!
subroutine flux_SW_pn(prim,nl,nx,ny,nz,nt,f ,efix,gama,npn)
    use global_const,only:sml_sss
    implicit none
    integer :: npn
    integer :: m,nl
    real :: f(nl),prim(nl)
    real :: hint,gama
    real :: l1,l4,l5
    real :: efix,eps
    real :: rm,um,vm,wm,pm,cm,eint,em,hm
    real :: nx,ny,nz,nt,ct,nxa,nya,nza,cta,cgm,cgm1
    real :: x1,x2,u2,c2r

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    u2 = 0.5 * ( um*um + vm*vm + wm*wm )
    call gethmgm(prim,hint,gama)
    eint = hint - pm/rm          !内能
    em   = eint + u2             !总能
    hm   = hint + u2
    !求声速
    cm  = sqrt(gama*pm/rm)
    ct  = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)   !要考虑到cgm为零的情况

    l1 = ct
    l4 = ct + cm * cgm
    l5 = ct - cm * cgm
    eps = efix*efix*cgm*cgm

    l1 = 0.5*(l1 + npn*sqrt(l1*l1 + eps))
    l4 = 0.5*(l4 + npn*sqrt(l4*l4 + eps))
    l5 = 0.5*(l5 + npn*sqrt(l5*l5 + eps))

    cgm1 = 1.0/cgm
    nxa  = nx * cgm1
    nya  = ny * cgm1
    nza  = nz * cgm1
    cta  = (ct - nt )* cgm1 ! 减去nt ，在cta中不包含nt     
    c2r  = cm * cm / gama

    x1 = c2r * ( 2.0*l1 - l4 - l5 )/( 2.0 * cm * cm )
    x2 = c2r * ( l4 - l5 )/( 2.0 * cm )
    hm = em + pm/rm

    f(1) = (l1    - x1             ) * rm
    f(2) = (l1*um - x1*um + nxa*x2 ) * rm
    f(3) = (l1*vm - x1*vm + nya*x2 ) * rm
    f(4) = (l1*wm - x1*wm + nza*x2 ) * rm
    f(5) = (l1*em - x1*hm + cta*x2 ) * rm

    do m=6,nl
       f(m) = prim(m) * f(1)
    enddo

    return
end subroutine flux_SW_pn
!_____________________________________________________________________!
subroutine flux_SW_mod(priml,primr,nl,nx,ny,nz,nt,f,efix,gamal,gamar)
  implicit none
  integer :: m,nl
  real :: priml(nl),primr(nl),f(nl),fl(nl),fr(nl)
  real :: nx,ny,nz,nt,gamal,gamar,efix
  real :: ql(nl),qr(nl),ql_m(nl),qr_m(nl),priml_m(nl),primr_m(nl)
  real :: wp,rwp,pr,pl,gp,gykb(4)
  
  call p2q(nl,priml,ql,gamal)
  call p2q(nl,primr,qr,gamar)
  
  pl = priml(5)
  pr = primr(5)
  gp = 0.5*abs(pl - pr)/min(pl, pr)   !!!不够好,有待改进
  wp = 0.5/(gp*gp + 1.0)
  rwp = 1.0 - wp
  
  ql_m(:) = rwp*ql(:) + wp*qr(:)
  qr_m(:) = rwp*qr(:) + wp*ql(:)
  
  call q2p(nl,ql_m,priml_m,gamal)
  call q2p(nl,qr_m,primr_m,gamar)
  
  gykb(1) = nx 
  gykb(2) = ny 
  gykb(3) = nz 
  gykb(4) = nt 

  call mxdq(priml_m,gykb,ql,fl,efix, 1)
  call mxdq(primr_m,gykb,qr,fr,efix,-1)
  
  do m=1,nl
     f(m) = fl(m)+fr(m)
  enddo

end subroutine flux_SW_mod
!=============================================================================!
subroutine flux_Roe(ql,qr,nl,nx,ny,nz,nt,f,efix,gama)
    use global_const,only:sml_sss
    implicit none
!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine flux_Roe is used to calculate the convect flux by the use of !
!		Roe's flux difference scheme                                              !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			March,2006                                                              !                                                          !
! Reference                                                                   !
!     Liu Xin.PhD thesis                                                      ！
! Input                                                                       !
!     ql,qr							:						variables of both sides                   ！
!     nl    						:						number of flux components                 ！
!     efix		 					:						constant of flux                          ！
!     gama		 					:						gas parameter                             ！
!     nx,ny,nz,nt  	    :						geometry parameters                       ！
! Output                                                                      !
!      f								:						convect flux                              ！
!-----------------------------------------------------------------------------!
    integer :: m,nl
    real :: f(nl),df(nl),fl(nl),fr(nl),ql(nl),qr(nl),qroe(nl),dq(nl)
    real :: efix,gama
    real :: l1,l4,l5,rm,um,vm,wm,cm,hm,l12,l42,l52
    real :: nx,ny,nz,nt,cta,cgm,cgm1,rrl,rrl_1,ccgm
		real :: c2,c,dcta,cta1,gama1,gama2,deti,det2,rodctac
		real :: rodcta,dp_c2,a1,a2,a3,a4,a5,a6,a7,a8

!		to calcualte fluxes of right and left sides
		call flux(ql,gama,nl,nx,ny,nz,nt,fl)
		call flux(qr,gama,nl,nx,ny,nz,nt,fr)

!   the difference of variables between right and left sides
		do m=1,nl
			dq(m)=qr(m)-ql(m)  
		enddo

!   to calculate H at left and right sides
		gama1 = gama-1.0
		gama2 = gama/gama1
		ql(5)=gama2*ql(5)/ql(1) + 0.5*(ql(2)**2+ql(3)**2+ql(4)**2)
		qr(5)=gama2*qr(5)/qr(1) + 0.5*(qr(2)**2+qr(3)**2+qr(4)**2)

!   to calculate density of Roe average
		qroe(1)=sqrt(ql(1)*qr(1))

!   to calculate velocity and H of Roe average
		rrl = sqrt(qr(1)/ql(1))
		rrl_1 = rrl+1
		do m=2,nl
			qroe(m)=(ql(m)+qr(m)*rrl)/rrl_1
		enddo

		rm = qroe(1)
		um = qroe(2)
		vm = qroe(3)
		wm = qroe(4)
		hm = qroe(5)

!   to calculate the speed of sound
		c2 = (hm-0.5*(um*um+vm*vm+wm*wm))*gama1
		if(c2<=0.0)then
			write(*,*)'sqrt(c2) is a math error in subroutine flux_Roe'
			stop
		endif
		c = sqrt(c2)

		cta = um*nx+vm*ny+wm*nz+nt

		cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
		cgm1=1.0/cgm

		nx = nx*cgm1
		ny = ny*cgm1
		nz = nz*cgm1

		ccgm = c*cgm

		l1 = abs(cta-ccgm)
		l4 = abs(cta      )
		l5 = abs(cta+ccgm)

!		Harten's entropy modification
		l12   = l1*l1
		l42   = l4*l4
		l52   = l5*l5

		deti  = efix*cgm  
		det2  = deti*deti  

		l1 = sqrt(l12+det2)
		l4 = sqrt(l42+det2)
		l5 = sqrt(l52+det2)

!		ctal= ql(2)*nx+ql(3)*ny+ql(4)*nz+nt
!		ctar= qr(2)*nx+qr(3)*ny+qr(4)*nz+nt
!   dcta= ctar-ctal
		dcta= dq(2)*nx+dq(3)*ny+dq(4)*nz

    cta1= cta*cgm1

		rodcta = qroe(1)*dcta
    dp_c2  = dq(5)/c2

		rodctac = rodcta/c
		
		a1 = l4*(dq(1)-dp_c2)
		a2 = l5*(dp_c2+rodctac)*0.5
		a3 = l1*(dp_c2-rodctac)*0.5
		a4 = a1+a2+a3
		a5 = c*(a2-a3)
		a6 = l4*(rm*dq(2)-nx*rodcta)
		a7 = l4*(rm*dq(3)-ny*rodcta)
		a8 = l4*(rm*dq(4)-nz*rodcta)

		df(1) = a4
		df(2) = um*a4+nx*a5+a6
		df(3) = vm*a4+ny*a5+a7
		df(4) = wm*a4+nz*a5+a8
		df(5) = hm*a4+cta1*a5+um*a6+vm*a7+wm*a8-c2*a1/gama1

		do m=1,nl
			f(m) = 0.5*(fr(m)+fl(m)-df(m))
		enddo

    return
end subroutine flux_Roe
!=============================================================================!
subroutine flux_roemx(ql,qr,neqns,sx,sy,sz,st,flr,efix,gamma)
    use global_variables,only : rpmin,small,sml_sss
    implicit none
    real,parameter :: one=1.0,half=0.5,zero=0.0
    integer :: neqns
    real :: ql(neqns),qr(neqns),flr(neqns)
    real :: sx,sy,sz,st,efix,gamma
    real :: sn,osn,nx,ny,nz
    real :: rl,ul,vl,wl,pl,el,tl,gl,ael
    real :: rr,ur,vr,wr,pr,er,tr,gr,aer
    real :: vnl,vnr,vl2,vr2,hl,hr,rvnl,rvnr,plr
    real :: dq(neqns),mxdq(neqns),qi(neqns)
    real :: rls,rrs,oors,gami,vi2,vi2p2,rpmlr
    real :: mi,mia,ci2,ci,vni,b1,b2,oob12,b1b2,cmxd
    real :: dvn,dps,dvx,dvy,dvz,dh,rmp,fden,gcdc
    integer :: m,i
    common /rpind/ i

    ! Shock-stable Roe (M1,M2) scheme,Sung-soo Kim,Chongam Kim
    ! CMAME(2004), similiar to Harten-Lax-van Leer scheme in JCP(1991),273
    sn = max(sqrt(sx*sx + sy*sy + sz*sz),sml_sss)
    osn = one/sn
    nx = sx*osn
    ny = sy*osn
    nz = sz*osn

    rl = ql(1) + small
    ul = ql(2)
    vl = ql(3)
    wl = ql(4)
    pl = ql(5) + small

    rr = qr(1) + small
    ur = qr(2)
    vr = qr(3)
    wr = qr(4)
    pr = qr(5) + small

	vl2 = ul*ul + vl*vl + wl*wl
	vr2 = ur*ur + vr*vr + wr*wr

    gl = gamma
    ael = gl - one
    el = pl/(rl*ael)

    gr = gamma
    aer = gr - one
    er = pr/(rr*aer)

    el = el + half*vl2
    hl = el + pl/rl

    er = er + half*vr2
    hr = er + pr/rr

    vnl = nx*ul + ny*vl + nz*wl
    vnr = nx*ur + ny*vr + nz*wr

    dq(1) = rr - rl
    dq(2) = rr*ur - rl*ul
    dq(3) = rr*vr - rl*vl
    dq(4) = rr*wr - rl*wl
    dq(5) = rr*hr - rl*hl

    rls = sqrt(rl)
    rrs = sqrt(rr)
    qi(1) = rls*rrs

    oors = one/(rls + rrs)
    rls = rls*oors
    rrs = rrs*oors

    qi(2) = rls*ul + rrs*ur
    qi(3) = rls*vl + rrs*vr
    qi(4) = rls*wl + rrs*wr
    qi(5) = rls*hl + rrs*hr

    !gami = half*(gl + gr)
    gami = rls*gl + rrs*gr

    vi2 = qi(2)**2 + qi(3)**2 + qi(4)**2
    vi2p2 = half*vi2
    ci2 = (gami-one)*(qi(5) - vi2p2)
    ci = sqrt(ci2)

    vni = qi(2)*nx + qi(3)*ny + qi(4)*nz
    b1 = max(zero,vni+ci,vnr+ci)
    b2 = min(zero,vni-ci,vnl-ci)
    
    oob12 = one/(b1 - b2)
    b1 = b1*oob12
    b1b2 = b1*b2
    b2 = b2*oob12
    
    mi = vni/ci
    mia = abs(mi)
    if (vi2 <= small) then
       fden = one
    else
       fden = mia**(one - rpmin(i))
    end if

    rpmlr = min(pl/pr,pr/pl)
    if (mia <= small) then
       gcdc = one
    else
       gcdc = mia**(one - rpmlr)
    end if
    
    dvn = vnr - vnl
    dvx = ur - ul
    dvy = vr - vl
    dvz = wr - wl
    dps = pr - pl
    dh  = hr - hl

    rmp = dq(1) - fden*dps/ci2

    mxdq(1) = rmp
    mxdq(2) = rmp*qi(2) + qi(1)*(dvx - nx*dvn)
    mxdq(3) = rmp*qi(3) + qi(1)*(dvy - ny*dvn)
    mxdq(4) = rmp*qi(4) + qi(1)*(dvz - nz*dvn)
    mxdq(5) = rmp*qi(5) + qi(1)* dh

    !! rvnl = rl*vnl
    !! rvnr = rr*vnr
    !! flr(i,1) = b1*(rvnl           ) - b2*(rvnr           ) + b1b2*dq(1)
    !! flr(i,2) = b1*(rvnl*ul + nx*pl) - b2*(rvnr*ur + nx*pr) + b1b2*dq(2)
    !! flr(i,3) = b1*(rvnl*vl + ny*pl) - b2*(rvnr*vr + ny*pr) + b1b2*dq(3)
    !! flr(i,4) = b1*(rvnl*wl + nz*pl) - b2*(rvnr*wr + nz*pr) + b1b2*dq(4)
    !! flr(i,5) = b1*(rvnl*hl        ) - b2*(rvnr*hr        ) + b1b2*dq(5)

    rvnl = b1*rl*vnl
    rvnr = -b2*rr*vnr
    plr = b1*pl - b2*pr
    flr(1) = rvnl    + rvnr             + b1b2*dq(1)
    flr(2) = rvnl*ul + rvnr*ur + nx*plr + b1b2*dq(2)
    flr(3) = rvnl*vl + rvnr*vr + ny*plr + b1b2*dq(3)
    flr(4) = rvnl*wl + rvnr*wr + nz*plr + b1b2*dq(4)
    flr(5) = rvnl*hl + rvnr*hr          + b1b2*dq(5)

    cmxd = -gcdc*b1b2/(one + mia)

    do m=1,neqns
       flr(m) = (flr(m) + cmxd*mxdq(m))*sn
    end do

end subroutine flux_roemx
!=============================================================================!
subroutine flux_vanleer(priml,primr,nl,nx,ny,nz,nt,f,efix,gama)
  implicit none
  integer :: m,nl
  real :: priml(nl),primr(nl),f(nl),fl(nl),fr(nl)
  real :: nx,ny,nz,nt,gama,efix

  call flux_vanleer_pn(priml,nl,nx,ny,nz,nt,fl,efix,gama, 1)
  call flux_vanleer_pn(primr,nl,nx,ny,nz,nt,fr,efix,gama,-1)
  do m=1,nl
     f(m) = fl(m)+fr(m)
  enddo

  return
end subroutine flux_vanleer
!_____________________________________________________________________!
subroutine flux_vanleer_pn(prim,nl,nx,ny,nz,nt,f,efix,gama,npn)
    use global_const,only:sml_sss
    implicit none
    integer :: nl,npn,m,n
    real    :: f(nl),prim(nl+1),efix
    real    :: bmvl,ctcm,rmct
    real    :: rm,um,vm,wm,u2,pm,cm,c2,hm,gama,hint
    real    :: nx,ny,nz,nt,ct,nxa,nya,nza,cta,cgm,cgm1

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    u2 = 0.5 * ( um*um + vm*vm + wm*wm )

    call gethmgm(prim,hint,gama)

    hm   = hint + u2

    c2  = gama*pm/rm
    if ( c2 < 0.0 ) then
			write(*,*)'sqrt(c2)is a math error in subroutine flux_vanleer!'
      write(*,*)'c2=gama*pm/rm,gama,pm,rm:',c2,gama,pm,rm
			stop
    endif
    cm = sqrt(c2)

    ct   = nx*um + ny*vm + nz*wm + nt
    cgm  = sqrt(nx*nx + ny*ny + nz*nz) + sml_sss   
    cgm1 = 1.0/cgm
    nxa  = nx  * cgm1
    nya  = ny  * cgm1
    nza  = nz  * cgm1
    cta  = ct  * cgm1

    bmvl = cta/cm

    if ( npn == 1 ) then
       if ( abs(bmvl) < 1.0 ) then
          ctcm = ( -cta + 2.0*cm)/gama 
          f(1) = 0.25*cgm*rm*cm*(bmvl + 1.0)**2
          f(2) = f(1)*(um + nxa*ctcm)
          f(3) = f(1)*(vm + nya*ctcm)
          f(4) = f(1)*(wm + nza*ctcm)
          f(5) = f(1)* hm
       else
          if ( bmvl >= 1.0 ) then
             rmct = rm * cta
             f(1) = cgm * rmct
             f(2) = cgm *(rmct*um + pm*nxa)
             f(3) = cgm *(rmct*vm + pm*nya)
             f(4) = cgm *(rmct*wm + pm*nza)
             f(5) = cgm * rmct*hm
          else
             do n=1,5
                f(n) = 0.0
             enddo
          endif
       endif
    endif
    if ( npn == -1 ) then
       if ( abs(bmvl) < 1.0 ) then
          ctcm = ( -cta - 2.0*cm )/gama 
          f(1) = -0.25*cgm*rm*cm*(bmvl - 1.0)**2
          f(2) =  f(1) *(um + nxa*ctcm)
          f(3) =  f(1) *(vm + nya*ctcm)
          f(4) =  f(1) *(wm + nza*ctcm)
          f(5) =  f(1) * hm
       else
          if ( bmvl <= -1.0 ) then
             rmct = rm * cta
             f(1) = cgm * rmct
             f(2) = cgm *(rmct*um + pm*nxa)
             f(3) = cgm *(rmct*vm + pm*nya)
             f(4) = cgm *(rmct*wm + pm*nza)
             f(5) = cgm * rmct*hm
          else
             do n=1,5
                f(n) = 0.0
             enddo
          endif
       endif
    endif

    do m=6,nl
       f(m) = prim(m) * f(1)
    enddo

    return
end subroutine flux_vanleer_pn
!_____________________________________________________________________!
subroutine fluxsw_pn(prim,nl,nx,ny,nz,nt,f ,efix,gama,npn)
    use global_const,only:sml_sss
    implicit none
    integer :: npn
    integer :: m,nl
    real :: f(nl),prim(nl)
    real :: hint,gama
    real :: l1,l4,l5
    real :: efix,eps
    real :: rm,um,vm,wm,pm,cm,eint,em,hm
    real :: nx,ny,nz,nt,ct,nxa,nya,nza,cta,cgm,cgm1
    real :: x1,x2,u2,c2r

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    u2 = 0.5 * ( um*um + vm*vm + wm*wm )
    call gethmgm(prim,hint,gama)
    eint = hint - pm/rm          !内能
    em   = eint + u2             !总能
    hm   = hint + u2
    !求声速
    cm  = sqrt(gama*pm/rm)
    ct  = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)   !要考虑到cgm为零的情况

    l1 = ct
    l4 = ct + cm * cgm
    l5 = ct - cm * cgm
    eps = efix*efix*cgm*cgm

    l1 = 0.5*(l1 + npn*sqrt(l1*l1 + eps))
    l4 = 0.5*(l4 + npn*sqrt(l4*l4 + eps))
    l5 = 0.5*(l5 + npn*sqrt(l5*l5 + eps))

    cgm1 = 1.0/cgm
    nxa  = nx * cgm1
    nya  = ny * cgm1
    nza  = nz * cgm1
    cta  = (ct - nt )* cgm1 ! 减去nt ，在cta中不包含nt     
    c2r  = cm * cm / gama

    x1 = c2r * ( 2.0*l1 - l4 - l5 )/( 2.0 * cm * cm )
    x2 = c2r * ( l4 - l5 )/( 2.0 * cm )
    hm = em + pm/rm

    f(1) = (l1    - x1             ) * rm
    f(2) = (l1*um - x1*um + nxa*x2 ) * rm
    f(3) = (l1*vm - x1*vm + nya*x2 ) * rm
    f(4) = (l1*wm - x1*wm + nza*x2 ) * rm
    f(5) = (l1*em - x1*hm + cta*x2 ) * rm
!    f(5) = f(1)  * hm

    do m=6,nl
       f(m) = prim(m) * f(1)
    enddo

    return
end subroutine fluxsw_pn
!_____________________________________________________________________!
!=============================================================================!
subroutine flux_AUSM_plus(ql,qr,nl,nx,ny,nz,nt,f,efix,gama)
    use global_const,only:sml_sss
    implicit none
!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine flux_AUSM_plus is used to calculate the convect flux by the use of !
!		AUSM+ flux difference scheme                                              !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			March,2006                                                              !                                                          !
! Reference                                                                   !
!     Zhao huiyong.PhD thesis                                                      ！
! Input                                                                       !
!     ql,qr							:						variables of both sides                   ！
!     nl    						:						number of flux components                 ！
!     efix		 					:						constant of flux                          ！
!     gama		 					:						gas parameter                             ！
!     nx,ny,nz,nt  	    :						geometry parameters                       ！
! Output                                                                      !
!      f								:						convect flux                              ！
!-----------------------------------------------------------------------------!
    integer :: nl
    real :: f(nl),ql(nl),qr(nl)
    real :: efix,gama
    real :: nx,ny,nz,nt,cgm,cgmc,ctal,ctar
		real :: c2,c,gama1,gama2,cgmcmrl,cgmcmrr
		real :: rl,ul,vl,wl,pl,rr,ur,vr,wr,pr,hl,hr,el,er
		real :: cl2,cr2,cl,cr,ml,mr,mi,mi_p,mi_n,pi,ul2,ur2
		real,external :: mpn,ppn

		rl = ql(1)
		ul = ql(2)
		vl = ql(3)
		wl = ql(4)
		pl = ql(5)

		rr = qr(1)
		ur = qr(2)
		vr = qr(3)
		wr = qr(4)
		pr = qr(5)

		gama1 = gama-1.0

		cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

		ctal= (ul*nx+vl*ny+wl*nz+nt)/cgm
		ctar= (ur*nx+vr*ny+wr*nz+nt)/cgm

		ul2 = 0.5*(ul*ul+vl*vl+wl*wl)
		ur2 = 0.5*(ur*ur+vr*vr+wr*wr)

		hl  = gama/gama1*pl/rl+ul2
		hr  = gama/gama1*pr/rr+ur2

		el  = pl/rl/gama1+ul2
		er  = pr/rr/gama1+ur2

		gama2 = 2.0*gama1/(gama+1.0)
		cl2 = gama2*hl
		cr2 = gama2*hr
		cl  = cl2/amax1(abs(ctal),sqrt(cl2))
		cr  = cr2/amax1(abs(ctar),sqrt(cr2))

		c   = amin1(cl,cr)

		ml  = ctal/c
		mr  = ctar/c

		mi  = mpn(ml,1.0,0.125)+mpn(mr,-1.0,0.125)
		mi_p= 0.5*(mi+abs(mi))
		mi_n= 0.5*(mi-abs(mi))

		pi  = ppn(ml,1.0,0.1875)*pl+ppn(mr,-1.0,0.1875)*pr

		cgmc    = cgm*c
		cgmcmrl = cgmc*mi_p*rl
		cgmcmrr = cgmc*mi_n*rr

		f(1) = cgmcmrl    + cgmcmrr
		f(2) = cgmcmrl*ul + cgmcmrr*ur + nx*pi
		f(3) = cgmcmrl*vl + cgmcmrr*vr + ny*pi
		f(4) = cgmcmrl*wl + cgmcmrr*wr + nz*pi
		f(5) = cgmcmrl*hl + cgmcmrr*hr 

		return
end subroutine flux_AUSM_plus 
!=============================================================================!
real function mpn(x,npn,w)
    implicit none
    real :: x,npn,x2,w
		x2 = x*x
		if(-1.0<x .and. x<1.0)then
			mpn = npn*(0.25*(x+npn)*(x+npn)+w*(x2-1.0)*(x2-1.0))
		else
			mpn = 0.5*(x+npn*abs(x))
		endif
    return
end function mpn
!=============================================================================!
real function ppn(x,npn,w)
    implicit none
    real :: x,npn,x2,w
		x2 = x*x
		if(-1.0<x .and. x<1.0)then
			ppn = 0.25*(x+npn)*(x+npn)*(2.0-npn*x) &
			     +w*npn*x*(x2-1.0)*(x2-1.0)
		else
			ppn = 0.5*(1.0+sign(1.0,x))
		endif
    return
end function ppn
!=============================================================================!
!=============================================================================!
subroutine flux_AUSMPW(ql,qr,nl,nx,ny,nz,nt,f,efix,gama)
    use global_const,only:sml_sss
    implicit none
!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine flux_AUSMPW is used to calculate the convect flux            !
!     by the use of AUSMPW flux difference scheme                             !                                          !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			March,2006                                                              !                                                          !
! Reference                                                                   !
!     Zhao huiyong.PhD thesis                                                 ！
! Input                                                                       !
!     ql,qr							:						variables of both sides                   ！
!     nl    						:						number of flux components                 ！
!     efix		 					:						constant of flux                          ！
!     gama		 					:						gas parameter                             ！
!     nx,ny,nz,nt  	    :						geometry parameters                       ！
! Output                                                                      !
!      f								:						convect flux                              ！
!-----------------------------------------------------------------------------!
    integer :: nl
    real :: f(nl),ql(nl),qr(nl)
    real :: efix,gama
    real :: nx,ny,nz,nt,cgm,cgmc,ctal,ctar
		real :: c2,c,gama1,gama2,cgmcmrl,cgmcmrr,gama3
		real :: rl,ul,vl,wl,pl,rr,ur,vr,wr,pr,hl,hr
		real :: cl2,cr2,cl,cr,ml,mr,mi,pi,ul2,ur2,ps
		real :: wp,fl,fr,mlp,mrn,ml0,mr0,g1,g2,g3,g4,g5
		real,external :: mpn,ppn,pw,pli

		rl = ql(1)
		ul = ql(2)
		vl = ql(3)
		wl = ql(4)
		pl = ql(5)

		rr = qr(1)
		ur = qr(2)
		vr = qr(3)
		wr = qr(4)
		pr = qr(5)

		gama1 = gama-1.0

		cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

		ctal= (ul*nx+vl*ny+wl*nz+nt)/cgm
		ctar= (ur*nx+vr*ny+wr*nz+nt)/cgm

		ul2 = 0.5*(ul*ul+vl*vl+wl*wl)
		ur2 = 0.5*(ur*ur+vr*vr+wr*wr)

		gama3 = gama/gama1
		hl  = gama3*pl/rl+ul2
		hr  = gama3*pr/rr+ur2

		gama2 = 2.0*gama1/(gama+1.0)
		cl2 = gama2*hl
		cr2 = gama2*hr
		cl  = cl2/amax1(abs(ctal),sqrt(cl2))
		cr  = cr2/amax1(abs(ctar),sqrt(cr2))

		c   = amin1(cl,cr)

		ml  = ctal/c
		mr  = ctar/c

		mlp = mpn(ml, 1.0,0.125)
		mrn = mpn(mr,-1.0,0.125)


		mi  = mlp+mrn

		pi  = ppn(ml,1.0,0.1875)*pl+ppn(mr,-1.0,0.1875)*pr
!		ps  = ppn(ml,1.0,0.1875)*pl+ppn(mr, 1.0,0.1875)*pr ! this maybe a equation error

		wp  = 1.0-(amin1(pl/pr,pr/pl))**3.0

		fl  = 0.0
		fr  = 0.0
		
		if(abs(ml)<=1.0)then
			ml0 = mpn(ml,1.0,0.0)
			fl  = (pl/pi-1.0)*pli(pl,pr)*abs(ml0)*(amin1(1.0,abs(ml)**0.25))
		endif
		if(abs(mr)<=1.0)then
			mr0 = mpn(mr,-1.0,0.0)
			fr  = (pr/pi-1.0)*pli(pl,pr)*abs(mr0)*(amin1(1.0,abs(mr)**0.25))
		endif

		cgmc    = cgm*c

		if(mi>0.0)then

			g1 = pw(rl   ,rr   ,wp)
			g2 = pw(rl*ul,rr*ur,wp)
			g3 = pw(rl*vl,rr*vr,wp)
			g4 = pw(rl*wl,rr*wr,wp)
			g5 = g1*hl

			cgmcmrl = cgmc*mlp*(1.0+fl)*rl
			cgmcmrr = cgmc*mrn*(1.0+fr)

			f(1) = cgmcmrl    + cgmcmrr*g1
			f(2) = cgmcmrl*ul + cgmcmrr*g2 + nx*pi
			f(3) = cgmcmrl*vl + cgmcmrr*g3 + ny*pi
			f(4) = cgmcmrl*wl + cgmcmrr*g4 + nz*pi
			f(5) = cgmcmrl*hl + cgmcmrr*g5 

		else

			g1 = pw(rr   ,rl   ,wp)
			g2 = pw(rr*ur,rl*ul,wp)
			g3 = pw(rr*vr,rl*vl,wp)
			g4 = pw(rr*wr,rl*wl,wp)
			g5 = g1*hr

			cgmcmrl = cgmc*mlp*(1.0+fl)
			cgmcmrr = cgmc*mrn*(1.0+fr)*rr

			f(1) = cgmcmrl*g1 + cgmcmrr   
			f(2) = cgmcmrl*g2 + cgmcmrr*ur + nx*pi 
			f(3) = cgmcmrl*g3 + cgmcmrr*vr + ny*pi
			f(4) = cgmcmrl*g4 + cgmcmrr*wr + nz*pi 
			f(5) = cgmcmrl*g5 + cgmcmrr*hr 

		endif

		return
end subroutine flux_AUSMPW
!=============================================================================!
real function pw(x,y,w)
    implicit none
    real :: x,y,w

		pw = (1.0-w)*x + w*y

    return
end function pw
!=============================================================================!
real function pli(x,y)
    implicit none
    real :: x,y,w

		w = amin1(x/y,y/x)

		if(0<=w .and.w<0.75)then
			pli = 0.0
		else
			pli = 4.0*w-3.0
		endif

    return
end function pli
!=============================================================================!
subroutine flux_AUSMPW_plus(ql,qr,nl,nx,ny,nz,nt,f,efix,gama)
    use global_const,only:sml_sss
    implicit none
!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine flux_AUSMPW is used to calculate the convect flux            !
!     by the use of AUSMPW+ flux difference scheme                            !                                          !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			March,2006                                                              !                                                          !
! Reference                                                                   !
!     Zhao huiyong.PhD thesis                                                 ！
! Input                                                                       !
!     ql,qr							:						variables of both sides                   ！
!     nl    						:						number of flux components                 ！
!     efix		 					:						constant of flux                          ！
!     gama		 					:						gas parameter                             ！
!     nx,ny,nz,nt  	    :						geometry parameters                       ！
! Output                                                                      !
!      f								:						convect flux                              ！
!-----------------------------------------------------------------------------!
    integer :: nl
    real :: f(nl),ql(nl),qr(nl)
    real :: efix,gama
    real :: nx,ny,nz,nt,cgm,cgmc,ctal,ctar
		real :: c2,c,gama1,gama2,cgmcmrl,cgmcmrr,gama3
		real :: rl,ul,vl,wl,pl,rr,ur,vr,wr,pr,hl,hr,egamal,egamar
		real :: cl2,cl,cr,ml,mr,mi,pi,ul2,ur2
		real :: wp,fl,fr,mlp,mrn
		real,external :: mpn,ppn

		rl = ql(1)
		ul = ql(2)
		vl = ql(3)
		wl = ql(4)
		pl = ql(5)

		rr = qr(1)
		ur = qr(2)
		vr = qr(3)
		wr = qr(4)
		pr = qr(5)

		gama1 = gama-1.0

		cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

		ctal= (ul*nx+vl*ny+wl*nz+nt)/cgm
		ctar= (ur*nx+vr*ny+wr*nz+nt)/cgm

!		total entropy in normal direction
		gama3 = gama/gama1
		egamal= gama3*pl/rl
		egamar= gama3*pr/rr
		hl  = egamal+0.5*ctal*ctal
		hr  = egamar+0.5*ctar*ctar

		gama2 =gama1/(gama+1.0)

		cl2 = gama2*(hl+hr)

		if((ctal+ctar)>0.0)then
			c=cl2/amax1(abs(ctal),sqrt(cl2))
		else
			c=cl2/amax1(abs(ctar),sqrt(cl2))
		endif

		ml  = ctal/c
		mr  = ctar/c

		mlp = mpn(ml, 1.0,0.0)
		mrn = mpn(mr,-1.0,0.0)

		mi  = mlp+mrn

		pi  = ppn(ml,1.0,0.0)*pl+ppn(mr,-1.0,0.0)*pr

		wp  = 1.0-(amin1(pl/pr,pr/pl))**3.0

!		total entropies on both sides
		ul2 = 0.5*(ul*ul+vl*vl+wl*wl)
		ur2 = 0.5*(ur*ur+vr*vr+wr*wr)
		hl  = egamal+ul2
		hr  = egamar+ur2

		fl  = 0.0
		fr  = 0.0
		
		if(abs(ml)<=1.0)then
			fl  = pl/pi-1.0
		endif

		if(abs(mr)<=1.0)then
			fr  = pr/pi-1.0
		endif

		cgmc    = cgm*c

		if(mi>0.0)then

			mlp = mlp+mrn*((1.0-wp)*(1.0+fr)-fl)
			mrn = mrn*wp*(1.0+fr)

		else

			mrn = mrn+mlp*((1.0-wp)*(1.0+fl)-fr)
			mlp = mlp*wp*(1.0+fl)

		endif

		cgmcmrl = cgmc*mlp*rl
		cgmcmrr = cgmc*mrn*rr

		f(1) = cgmcmrl    + cgmcmrr   
		f(2) = cgmcmrl*ul + cgmcmrr*ur + nx*pi 
		f(3) = cgmcmrl*vl + cgmcmrr*vr + ny*pi
		f(4) = cgmcmrl*wl + cgmcmrr*wr + nz*pi 
		f(5) = cgmcmrl*hl + cgmcmrr*hr 

		return
end subroutine flux_AUSMPW_plus
!=============================================================================!
