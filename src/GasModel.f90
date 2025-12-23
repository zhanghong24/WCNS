!___________________________________________________________________!
subroutine air(h1,t,p,den,a)
    real :: h1,t,p,den,a
    h      = h1*1000.
    r      = 287.053
    g0     = 9.80665
    rp     = 6.37111e6
    g      = (rp/(rp+h))**2*g0
    t0     = 288.15
    p0     = 10.1325e2
    rho0   = 1.225
    t11    = 216.65
    p11    = 2.2632e2
    rho11  = 3.6392e-1
    t20    = t11
    p20    = 5.4747e1
    rho20  = 8.8035e-2
    t32    = 228.65
    p32    = 8.6789
    rho32  = 1.3225e-2
    t47    = 270.65
    p47    = 1.1090
    rho47  = 1.4275e-3
    t52    = t47
    p52    = 5.8997e-1
    rho52  = 7.5943e-4
    t61    = 252.65
    p61    = 1.8209e-1
    rho61  = 2.5109e-4
    t79    = 180.65
    p79    = 1.0376e-2
    rho79  = 2.0010e-5
    t90    = t79
    p90    = 1.6437e-3
    rho90  = 3.4165e-6
    t100   = 210.02
    p100   = 3.0070e-4
    rho100 = 5.6044e-7
    t110   = 257.00
    p110   = 7.3527e-5
    rho110 = 9.7081e-8
    t120   = 349.49
    p120   = 2.5209e-5
    rho120 = 2.2222e-8
    t150   = 892.79
    p150   = 5.0599e-6
    rho150 = 2.0752e-9
    t160   = 1022.20
    p160   = 3.6929e-6
    rho160 = 1.2336e-9
    t170   = 1103.40
    p170   = 2.7915e-6
    rho170 = 7.8155e-10
    t190   = 1205.40
    p190   = 1.6845e-6
    rho190 = 3.5807e-10
    t230   = 1322.30
    p230   = 6.7138e-7
    rho230 = 1.5640e-10
    t300   = 1432.10
    p300   = 1.8828e-7
    rho300 = 1.9159e-11
    t400   = 1487.40
    p400   = 4.0278e-8
    rho400 = 2.8028e-12
    t500   = 1499.20
    p500   = 1.0949e-8
    rho500 = 5.2148e-13
    t600   = 1506.10
    p600   = 3.4475e-9
    rho600 = 1.1367e-13
    t700   = 1507.60
    p700   = 1.1908e-9
    rho700 = 1.5270e-13
     
    if ( h <= 11019.0 ) then
       al1 = ( t11 - t0 )/11019.0
       t   = t0 + al1 * h
       p   = p0 * ( t/t0 )**( -g/(r * al1 ) )
       rho = rho0 * ( t/t0 )**(-1.0 - g/(r*al1) )
    elseif ( h <= 20063.0 ) then
       t   = t11
       p   = p11   * exp( -g*(h-11019.0)/(r*t11) )
       rho = rho11 * exp( -g*(h-11019.0)/(r*t11) )
    elseif ( h <= 32162.0 ) then
       al2 = (t32-t20)/(32162.0-20063.0)
       t   = t11+al2*(h-20063.0)
       p   = p20*(t/t11)**(-g/(r*al2))
       rho = rho20*(t/t11)**(-1.0-g/(r*al2))
    elseif ( h <= 47350.0 ) then
       al3 = (t47-t32)/(47350.0-32162.0)
       t   = t32+al3*(h-32162.0)
       p   = p32*(t/t32)**(-g/(r*al3))
       rho = rho32*(t/t32)**(-1.0-g/(r*al3))
    elseif ( h <= 52429.0 ) then
       t   = t47
       p   = p47*exp(-g*(h-47350.0)/(r*t47))
       rho = rho47*exp(-g*(h-47350.0)/(r*t47))
    else if(h.le.61591.0)then
         al4=(t61-t52)/(61591.0-52429.0)
         t=t47+al4*(h-52429.0)
         p=p52*(t/t47)**(-g/(r*al4))
         rho=rho52*(t/t47)**(-1.0-g/(r*al4))
    else if(h.le.79994.0)then
         al5=(t79-t61)/(79994.0-61591.0)
         t=t61+al5*(h-61591.0)
         p=p61*(t/t61)**(-g/(r*al5))
         rho=rho61*(t/t61)**(-1.0-g/(r*al5))
    else if(h.le.90000.0)then
         t=t79
         p=p79*exp(-g*(h-79994.0)/(r*t79))
         rho=rho79*exp(-g*(h-79994.0)/(r*t79))
    else if(h.le.100000.0)then
         al6=(t100-t90)/10000.0
         t=t79+al6*(h-90000.0)
         p=p90*(t/t79)**(-g/(r*al6))
         rho=rho90*(t/t79)**(-1.0-g/(r*al6))
    else if(h.le.110000.0)then
         al7=(t110-t100)/10000.0
         t=t100+al7*(h-100000.0)
         p=p100*(t/t100)**(-g/(r*al7))
         rho=rho100*(t/t100)**(-1.0-g/(r*al7))
    else if(h.le.120000.0)then
         al8=(t120-t110)/10000.0
         t=t110+al8*(h-110000.0)
         p=p110*(t/t110)**(-g/(r*al8))
         rho=rho110*(t/t110)**(-1.0-g/(r*al8))
    else if(h.le.150000.0)then
         al9=(t150-t120)/30000.0
         t=t120+al9*(h-120000.0)
         p=p120*(t/t120)**(-g/(r*al9))
         rho=rho120*(t/t120)**(-1.0-g/(r*al9))
    else if(h.le.160000.0)then
         al10=(t160-t150)/10000.0
         t=t150+al10*(h-150000.0)
         p=p150*(t/t150)**(-g/(r*al10))
         rho=rho150*(t/t150)**(-1.0-g/(r*al10))
    else if(h.le.170000.0)then
         al11=(t170-t160)/10000.0
         t=t160+al11*(h-160000.0)
         p=p160*(t/t160)**(-g/(r*al11))
         rho=rho160*(t/t160)**(-1.0-g/(r*al11))
    else if(h.le.190000.0)then
         al12=(t190-t170)/20000.0
         t=t170+al12*(h-170000.0)
         p=p170*(t/t170)**(-g/(r*al12))
         rho=rho170*(t/t170)**(-1.0-g/(r*al12))
    else if(h.le.230000.0)then
         al13=(t230-t190)/40000.0
         t=t190+al13*(h-190000.0)
         p=p190*(t/t190)**(-g/(r*al13))
         rho=rho190*(t/t190)**(-1.0-g/(r*al13))
    else if(h.le.300000.0)then
         al14=(t300-t230)/70000.0
         t=t230+al14*(h-230000.0)
         p=p230*(t/t230)**(-g/(r*al14))
         rho=rho230*(t/t230)**(-1.0-g/(r*al14))
    else if(h.le.400000.0)then
         al15=(t400-t300)/100000.0
         t=t300+al15*(h-300000.0)
         p=p300*(t/t300)**(-g/(r*al15))
         rho=rho300*(t/t300)**(-1.0-g/(r*al15))
    else if(h.le.500000.0)then
         al16=(t500-t400)/100000.0
         t=t400+al16*(h-400000.0)
         p=p400*(t/t400)**(-g/(r*al16))
         rho=rho400*(t/t400)**(-1.0-g/(r*al16))
    else if(h.le.600000.0)then
        al17=(t600-t500)/100000.0
        t=t500+al17*(h-500000.0)
        p=p500*(t/t500)**(-g/(r*al17))
        rho=rho500*(t/t500)**(-1.0-g/(r*al17))
    else if(h.le.700000.0)then
        al18=(t700-t600)/100000.0
        t=t600+al18*(h-600000.0)
        p=p600*(t/t600)**(-g/(r*al18))
        rho=rho600*(t/t600)**(-1.0-g/(r*al18))
    endif
    a   = sqrt(1.4*r*t)
    p   = p * 100.
    den = rho !/rho0
end subroutine air
!___________________________________________________________________!
subroutine read_perfect_gasmodel   !( nfile_par )
    use global_const
    implicit none
    integer :: nfile_par
    integer :: is,ir

    if( gasmodel == 'air.dat') then
	    ns = 2
		nf = 0
       allocate( varname(ns) )
       allocate( cn_init(ns) )
       allocate( ws(ns) , ms(ns) , ws1(ns) , ms1(ns) )
       !读各组元名称
       varname(1) = 'N2'
       varname(2) = 'O2'

       !读各组元分子量
       ws(1) = 28.0
       ws(2) = 32.0
       cn_init(1) = 0.79
       cn_init(2) = 0.21
       gama = 1.40
	   prl  = 0.72
	   prt  = 0.90
	else

	endif
    write(*,*)(trim(varname(is)),' ',is=1,ns)
    write(*,*)(ws(is),is=1,ns)
    write(*,*)(cn_init(is),is=1,ns)
!    call set_comp_par_nochem(nfile_par)

    return
end subroutine read_perfect_gasmodel
!_____________________________________________________________________!
subroutine compute_visl_ns
    use global_variables
    implicit none
    integer :: i,j,k,cmethod
    real :: tm
    cmethod = method - 1
    do k=1,nk + cmethod
       do j=1,nj + cmethod
          do i=1,ni + cmethod
             tm = t(i,j,k)
             visl(i,j,k) = tm*sqrt(tm)*(1.0+visc)/(tm+visc)
          enddo
       enddo
    enddo
    
    return
end subroutine compute_visl_ns

!_____________________________________________________________________!
subroutine get_c_t
    use global_variables
    implicit none
    integer :: i,j,k,cmethod
    integer :: nn(4)
    real :: prim(nl),q_q(nl)
    real :: a2,m2
    real :: rm,um,vm,wm,pm,em
    m2 = moo * moo
    cmethod = 1 - method 
    do k=1,nk - cmethod
       do j=1,nj - cmethod
          do i=1,ni - cmethod
             prim(1) = r(i,j,k)
             prim(2) = u(i,j,k)
             prim(3) = v(i,j,k)
             prim(4) = w(i,j,k)
             prim(5) = p(i,j,k)
             call prim_to_q(prim,q_q,gama)
             q(1,i,j,k) = q_q(1)
             q(2,i,j,k) = q_q(2)
             q(3,i,j,k) = q_q(3)
             q(4,i,j,k) = q_q(4)
             q(5,i,j,k) = q_q(5)

             a2 = gama * p(i,j,k) / r(i,j,k)
             t(i,j,k) = m2 * a2
             c(i,j,k) = sqrt( a2 )
             
          enddo
       enddo
    enddo

    nn(1) = -1 + method
    nn(2) = 0
    nn(3) = ni
    nn(4) = ni+1 

    do k=1,nk- cmethod
       do j=1,nj- cmethod
          do i=1,4
             prim(1) = r(nn(i),j,k)
             prim(2) = u(nn(i),j,k)
             prim(3) = v(nn(i),j,k)
             prim(4) = w(nn(i),j,k)
             prim(5) = p(nn(i),j,k)
             call prim_to_q_bc(prim,q_q,gama)   !prim_to_q(prim,q_q,gama)             
             q(1,nn(i),j,k) = q_q(1)
             q(2,nn(i),j,k) = q_q(2)
             q(3,nn(i),j,k) = q_q(3)
             q(4,nn(i),j,k) = q_q(4)
             q(5,nn(i),j,k) = q_q(5)

             r(nn(i),j,k) = prim(1)
             a2 = gama * p(nn(i),j,k) /  prim(1)   !r(nn(i),j,k)
             t(nn(i),j,k) = m2 * a2
             c(nn(i),j,k) = sqrt( a2 )
          enddo
       enddo
    enddo
    nn(3) = nj
    nn(4) = nj+1 
    do k=1,nk- cmethod
       do i=1,ni- cmethod
          do j=1,4
             prim(1) = r(i,nn(j),k)
             prim(2) = u(i,nn(j),k)
             prim(3) = v(i,nn(j),k)
             prim(4) = w(i,nn(j),k)
             prim(5) = p(i,nn(j),k)
             call prim_to_q_bc(prim,q_q,gama)   !prim_to_q(prim,q_q,gama)             
             q(1,i,nn(j),k) = q_q(1)
             q(2,i,nn(j),k) = q_q(2)
             q(3,i,nn(j),k) = q_q(3)
             q(4,i,nn(j),k) = q_q(4)
             q(5,i,nn(j),k) = q_q(5)

             r(i,nn(j),k) = prim(1)
             a2 = gama * p(i,nn(j),k) / prim(1)   !r(i,nn(j),k)
             if ( a2 <= 0.0 ) then
                write(*,*)i,j,k,' p=',p(i,nn(j),k),' r=',r(i,nn(j),k)
             endif
             t(i,nn(j),k) = m2 * a2
             c(i,nn(j),k) = sqrt( a2 )

          enddo
       enddo
    enddo

    nn(3) = nk
    nn(4) = nk+1
    do j=1,nj- cmethod
       do i=1,ni- cmethod
          do k=1,4
             prim(1) = r(i,j,nn(k))
             prim(2) = u(i,j,nn(k))
             prim(3) = v(i,j,nn(k))
             prim(4) = w(i,j,nn(k))
             prim(5) = p(i,j,nn(k))
             call prim_to_q_bc(prim,q_q,gama)   !prim_to_q(prim,q_q,gama)                          
             q(1,i,j,nn(k)) = q_q(1)
             q(2,i,j,nn(k)) = q_q(2)
             q(3,i,j,nn(k)) = q_q(3)
             q(4,i,j,nn(k)) = q_q(4)
             q(5,i,j,nn(k)) = q_q(5)

             r(i,j,nn(k)) = prim(1)
             a2 = gama * p(i,j,nn(k)) /  prim(1)   !r(i,j,nn(k))
             if ( a2 <= 0.0 ) then
                write(*,*)i,j,k,' p=',p(i,j,nn(k)),' r=',r(i,j,nn(k))
             endif
             t(i,j,nn(k)) = m2 * a2
             c(i,j,nn(k)) = sqrt( a2 )
          enddo
       enddo
    enddo

    return
end subroutine get_c_t
!_____________________________________________________________________!
subroutine getem(prim,em,gama)
    use global_const,only:nl,ns,nchem,nchem_source,tref,beta1,ms1,moo
    implicit none
    integer :: m
    real :: prim(nl),hs(ns)
    real :: rm,um,vm,wm,pm,em,hm,tm,gama
    real :: mav1,zf

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    em = pm/(gama-1.0) + 0.5* rm *(um*um +vm*vm +wm*wm)
    return
end subroutine getem
!_____________________________________________________________________!
subroutine gethmgm(prim,hm,gama)
    use global_const,only:nl
    implicit none
    real :: prim(nl)
    real :: hm,gama
          call gethmgm_ns(prim,hm,gama)
    return
end subroutine gethmgm
!_____________________________________________________________________!
subroutine gethmgm_ns(prim,hm,gamamm)
    use global_const,only : nl,gama
    implicit none
    real :: prim(nl)
    real :: hm,gamamm
    real :: rm,pm

    rm = prim(1)
    pm = prim(5)
    gamamm = gama

    hm = gama/(gama-1.0)*pm/rm

    return
end subroutine gethmgm_ns
!_____________________________________________________________________!
subroutine gettmhmgm(prim,tm,hm,gama)
    use global_const,only : nl
    implicit none
    real :: prim(nl)
    real :: tm,hm,gama
          call gettmhmgm_ns(prim,tm,hm,gama)
    return
end subroutine gettmhmgm

!_____________________________________________________________________!
subroutine gettmhmgm_ns(prim,tm,hm,gamamm)
    use global_const,only : nl,moo,gama
    implicit none
    real :: prim(nl)
    real :: tm,hm,gamamm
    real :: rm,pm

    rm = prim(1)
    pm = prim(5)

    gamamm= gama
    tm = gama*moo*moo*pm/rm
    hm = gama/(gama-1.0)*pm/rm

    return
end subroutine gettmhmgm_ns



