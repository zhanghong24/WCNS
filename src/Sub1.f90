!_____________________________________________________________________!
subroutine infinite_interplat1(f,ni,nj,nk,npdi,npdj,method)
    implicit none
    integer :: ns,ni,nj,nk,niend,njend,nkend,nmeth
    integer :: i,j,k,m,npdi,npdj,method
    real :: eps,fmin,fmax
    real :: cil,cir,cjl,cjr,ckl,ckr
    real :: f(npdi:ni+npdj,npdi:nj+npdj,npdi:nk+npdj)

    nmeth = 1 + method
!   对初场进行光滑
    do niend = 1,30
       do k=nmeth,nk-1 
          do j=nmeth,nj-1 
             do i=nmeth,ni-1 
                cir = f(i-1,j,k)-2.0*f(i,j,k)+f(i+1,j,k)
                cjr = f(i,j-1,k)-2.0*f(i,j,k)+f(i,j+1,k)
                ckr = f(i,j,k-1)-2.0*f(i,j,k)+f(i,j,k+1)
                f(i,j,k) = f(i,j,k)+ ( cir + cjr + ckr )/8.0
             enddo
          enddo
       enddo
    enddo
    return
end subroutine infinite_interplat1
!_____________________________________________________________________!
subroutine infinite_interplat2(f,ns,ni,nj,nk,npdi,npdj,method)
    implicit none
    integer :: ns,ni,nj,nk,niend,njend,nkend,nmeth
    integer :: i,j,k,m,npdi,npdj,method
    real :: f(ns,npdi:ni+npdj,npdi:nj+npdj,npdi:nk+npdj)
    real :: cil,cir,cjl,cjr,ckl,ckr

    nmeth = 1 + method
!   对初场进行光滑
    do niend = 1,30
       do k=nmeth,nk-1 
          do j=nmeth,nj-1 
             do i=nmeth,ni-1 
                do m=1,ns
                   cir = f(m,i-1,j,k) - 2.0 * f(m,i,j,k) + f(m,i+1,j,k)
                   cjr = f(m,i,j-1,k) - 2.0 * f(m,i,j,k) + f(m,i,j+1,k)
                   ckr = f(m,i,j,k-1) - 2.0 * f(m,i,j,k) + f(m,i,j,k+1)
                   f(m,i,j,k) = f(m,i,j,k)+ ( cir + cjr + ckr )/8.0
                enddo
             enddo
          enddo
       enddo
    enddo
    return
end subroutine infinite_interplat2
!__________________________________________!
subroutine nxyz(a,b,c,d,kxyz)
    !本子程序用于计算方向面积
    implicit none
    integer :: n
    real,dimension(1:3) :: a,b,c,d,kxyz
    real,dimension(1:3) :: v1,v2,p1,p2

    do n=1,3
       v1(n) = b(n) - a(n)
       v2(n) = c(n) - b(n)
    enddo

    call cross_product(v1,v2,p1)

    do n=1,3
       v1(n) = d(n) - c(n)
       v2(n) = a(n) - d(n)
    enddo

    call cross_product(v1,v2,p2)

    do n=1,3
       kxyz(n) = 0.5 * ( p1(n) + p2(n) )
    enddo

    return
end subroutine nxyz
!_____________________________________________________________________!
subroutine cross_product(a,b,cp)
!本子程序用于计算叉积, 向量积, 交叉乘积(cross product)
    implicit none
    real,dimension(1:3) :: a,b,cp
    cp(1) = a(2) * b(3) - a(3) * b(2)
    cp(2) = a(3) * b(1) - a(1) * b(3)
    cp(3) = a(1) * b(2) - a(2) * b(1)
    return;
end subroutine cross_product
!_____________________________________________________________________!
subroutine xyzcoor(nb,i,j,k,xc,yc,zc)
    use global_variables,only : mb_x,mb_y,mb_z
    implicit none
    integer :: nb,i,j,k
    real :: cc
    parameter(cc=0.125)
    real :: xc,yc,zc

    xc = mb_x(nb)%a3d(i,j,k)
    yc = mb_y(nb)%a3d(i,j,k)
    zc = mb_z(nb)%a3d(i,j,k)

    return
end subroutine xyzcoor
!_____________________________________________________________________!
subroutine prim_to_q(prim,q,gama)
    use global_const,only:nl
    implicit none
    integer :: m,n
    real :: prim(nl),q(nl)
    real :: rm,um,vm,wm,pm,em,gama

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)
    call getem(prim,em,gama)

    q(1) = rm
    q(2) = rm * um
    q(3) = rm * vm
    q(4) = rm * wm
    q(5) = em

    do m=6,nl
       q(m) = rm * prim(m)
    enddo

    return
end subroutine prim_to_q
!_____________________________________________________________________!
subroutine prim_to_q_bc(prim,q,gama)
    use global_const,only:nl
    implicit none
    integer :: m
    real :: prim(nl),q(nl)
    real :: rm,um,vm,wm,pm,em,gama

    q(1) = prim(1)
    if(prim(1) < 0 ) then
       prim(1) = - prim(1)
	  endif
	  rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)
    call getem(prim,em,gama)

    q(2) = rm * um
    q(3) = rm * vm
    q(4) = rm * wm
    q(5) = em

    do m=6,nl
       q(m) = rm * prim(m)
    enddo
    return
end subroutine prim_to_q_bc
!_____________________________________________________________________!
subroutine getprim(i,j,k,prim)
    use global_variables
    implicit none
    integer :: i,j,k,m
    real,dimension( 1:nl ) :: prim

    prim(1) = r(i,j,k)
    prim(2) = u(i,j,k)
    prim(3) = v(i,j,k)
    prim(4) = w(i,j,k)
    prim(5) = p(i,j,k)

    do m=6,nl
       prim(m) = fs(m-5,i,j,k)
    enddo

    return
end subroutine getprim
!_____________________________________________________________________!
subroutine stress(vis)
    use stress_module
    implicit none
    real :: vis,c23
    parameter(c23=0.66666666667)
    
    txx = vis * c23 * ( 2.0*dudx - dvdy - dwdz )
    tyy = vis * c23 * ( 2.0*dvdy - dwdz - dudx )
    tzz = vis * c23 * ( 2.0*dwdz - dudx - dvdy )
    txy = vis * ( dudy + dvdx )
    txz = vis * ( dudz + dwdx )
    tyz = vis * ( dvdz + dwdy )
    tyx = txy
    tzx = txz
    tzy = tyz

    return
end subroutine stress
!____________________________________________!
subroutine trisys(ldiag,diag,udiag,b,x,n)
    !求解三对角方程
    implicit none
    integer :: n,nn,j,k,kk
    real :: ldiag(n),diag(n),udiag(n),b(n),x(n)

    do j=1,n
       x(j) = b(j)
    enddo
 
    if ( n /= 1 ) goto 2
    x(1) = x(1)/diag(1)
2   do j=2,n
       ldiag(j) = ldiag(j)/diag(j-1)
       diag(j)  = diag(j) - ldiag(j) * udiag(j-1)
       x(j)     = x(j) - ldiag(j) * x(j-1)
    enddo

    x(n) = x(n)/diag(n)
    nn = n - 1
    do k=1,nn
       kk = n - k
       x(kk) = ( x(kk) - udiag(kk) * x(kk+1) )/diag(kk)
    enddo

    return
end subroutine trisys
!_____________________________________________________________________!
subroutine getgeo(i,j,k)
    use global_variables
    use geometry_module
    implicit none

    integer :: i,j,k

    kx = kcx(i,j,k)
    ky = kcy(i,j,k)
    kz = kcz(i,j,k)

    ex = etx(i,j,k)
    ey = ety(i,j,k)
    ez = etz(i,j,k)

    cx = ctx(i,j,k)
    cy = cty(i,j,k)
    cz = ctz(i,j,k)

    vjacob = vol(i,j,k)
    return
end subroutine getgeo
!_____________________________________________________________________!
subroutine getduvwdxyz
    use geometry_module
    use stress_module,only:dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    use duvwt_module,only:ukc,uet,uct,vkc,vet,vct,wkc,wet,wct
    implicit none

    dudx = (kx*ukc + ex*uet + cx*uct)/vjacob
    dudy = (ky*ukc + ey*uet + cy*uct)/vjacob
    dudz = (kz*ukc + ez*uet + cz*uct)/vjacob

    dvdx = (kx*vkc + ex*vet + cx*vct)/vjacob
    dvdy = (ky*vkc + ey*vet + cy*vct)/vjacob
    dvdz = (kz*vkc + ez*vet + cz*vct)/vjacob

    dwdx = (kx*wkc + ex*wet + cx*wct)/vjacob
    dwdy = (ky*wkc + ey*wet + cy*wct)/vjacob
    dwdz = (kz*wkc + ez*wet + cz*wct)/vjacob

    return
end subroutine getduvwdxyz
!_____________________________________________________________________!
subroutine getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)
    use geometry_module
    use duvwt_module
    implicit none

    real :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dtdx,dtdy,dtdz
    dudx = (kx*ukc + ex*uet + cx*uct)/vjacob
    dudy = (ky*ukc + ey*uet + cy*uct)/vjacob
    dudz = (kz*ukc + ez*uet + cz*uct)/vjacob

    dvdx = (kx*vkc + ex*vet + cx*vct)/vjacob
    dvdy = (ky*vkc + ey*vet + cy*vct)/vjacob
    dvdz = (kz*vkc + ez*vet + cz*vct)/vjacob

    dwdx = (kx*wkc + ex*wet + cx*wct)/vjacob
    dwdy = (ky*wkc + ey*wet + cy*wct)/vjacob
    dwdz = (kz*wkc + ez*wet + cz*wct)/vjacob

    dtdx = (kx*tkc + ex*tet + cx*tct)/vjacob
    dtdy = (ky*tkc + ey*tet + cy*tct)/vjacob
    dtdz = (kz*tkc + ez*tet + cz*tct)/vjacob

    return
end subroutine getdfdxyz
!_____________________________________________________________________!
subroutine getuvwtder0(i,j,k)
    use global_variables
    use duvwt_module
    implicit none

    integer :: i,j,k
    real :: aa,bb,cc
    real :: u1,v1,w1,t1,u2,v2,w2,t2,u3,v3,w3,t3,up,vp,wp,tp,um,vm,wm,tm

    aa =  3.0*0.5
    bb = -4.0*0.5
    cc =  1.0*0.5

    call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
    if ( i == 1 ) then
       call getuvwt( 2 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( 3 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       ukc = - ( aa * u1 + bb * u2 + cc * u3 )
       vkc = - ( aa * v1 + bb * v2 + cc * v3 )
       wkc = - ( aa * w1 + bb * w2 + cc * w3 )
       tkc = - ( aa * t1 + bb * t2 + cc * t3 )
    elseif ( i == ni ) then
       call getuvwt( i-1 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( i-2 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       ukc = aa * u1 + bb * u2 + cc * u3
       vkc = aa * v1 + bb * v2 + cc * v3
       wkc = aa * w1 + bb * w2 + cc * w3
       tkc = aa * t1 + bb * t2 + cc * t3
    else
       call getuvwt( i+1 ,j ,k ,up ,vp ,wp ,tp )
       call getuvwt( i-1 ,j ,k ,um ,vm ,wm ,tm )
       ukc = 0.5 * ( up - um )
       vkc = 0.5 * ( vp - vm )
       wkc = 0.5 * ( wp - wm )
       tkc = 0.5 * ( tp - tm )
    endif

    if ( j == 1 ) then
        call getuvwt( i ,2 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,3 ,k ,u3 ,v3 ,w3 ,t3 )
        uet = - ( aa * u1 + bb * u2 + cc * u3 )
        vet = - ( aa * v1 + bb * v2 + cc * v3 )
        wet = - ( aa * w1 + bb * w2 + cc * w3 )
        tet = - ( aa * t1 + bb * t2 + cc * t3 )
    elseif ( j == nj ) then
        call getuvwt( i ,j-1 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j-2 ,k ,u3 ,v3 ,w3 ,t3 )
        uet = aa * u1 + bb * u2 + cc * u3
        vet = aa * v1 + bb * v2 + cc * v3
        wet = aa * w1 + bb * w2 + cc * w3
        tet = aa * t1 + bb * t2 + cc * t3
    else
        call getuvwt( i ,j+1 ,k ,up ,vp ,wp ,tp )
        call getuvwt( i ,j-1 ,k ,um ,vm ,wm ,tm )
        uet = 0.5 * ( up - um )
        vet = 0.5 * ( vp - vm )
        wet = 0.5 * ( wp - wm )
        tet = 0.5 * ( tp - tm )
    endif

    if ( k == 1 ) then
        call getuvwt( i ,j ,2 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,3 ,u3 ,v3 ,w3 ,t3 )
        uct = - ( aa * u1 + bb * u2 + cc * u3 )
        vct = - ( aa * v1 + bb * v2 + cc * v3 )
        wct = - ( aa * w1 + bb * w2 + cc * w3 )
        tct = - ( aa * t1 + bb * t2 + cc * t3 )
    elseif ( k == nk ) then
        call getuvwt( i ,j ,k-1 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,k-2 ,u3 ,v3 ,w3 ,t3 )
        uct = aa * u1 + bb * u2 + cc * u3
        vct = aa * v1 + bb * v2 + cc * v3
        wct = aa * w1 + bb * w2 + cc * w3
        tct = aa * t1 + bb * t2 + cc * t3
    else
        call getuvwt( i ,j ,k+1 ,up ,vp ,wp ,tp )
        call getuvwt( i ,j ,k-1 ,um ,vm ,wm ,tm )
        uct = 0.5 * ( up - um )
        vct = 0.5 * ( vp - vm )
        wct = 0.5 * ( wp - wm )
        tct = 0.5 * ( tp - tm )
    endif

    return
end subroutine getuvwtder0
!_____________________________________________________________________!
subroutine getuvwtder(i,j,k)
    use global_variables
    use duvwt_module
    implicit none

    integer :: i,j,k
    real :: aa,bb,cc
    real :: u1,v1,w1,t1,u2,v2,w2,t2,u3,v3,w3,t3,up,vp,wp,tp,um,vm,wm,tm
!    real :: AC1,AC2,AC3,AC4,AC5,u4,v4,w4,t4,u5,v5,w5,t5

    aa =  1.5
    bb = -2.0
    cc =  0.5

    if ( i == 1 ) then
       call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1 )
       call getuvwt( 2 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( 3 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       ukc = - ( aa * u1 + bb * u2 + cc * u3 )
       vkc = - ( aa * v1 + bb * v2 + cc * v3 )
       wkc = - ( aa * w1 + bb * w2 + cc * w3 )
       tkc = - ( aa * t1 + bb * t2 + cc * t3 )

    elseif ( i == ni ) then
       call getuvwt( i   ,j ,k ,u1 ,v1 ,w1 ,t1 )
       call getuvwt( i-1 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( i-2 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       ukc = aa * u1 + bb * u2 + cc * u3
       vkc = aa * v1 + bb * v2 + cc * v3
       wkc = aa * w1 + bb * w2 + cc * w3
       tkc = aa * t1 + bb * t2 + cc * t3

    else
       call getuvwt( i+1 ,j ,k ,up ,vp ,wp ,tp )
       call getuvwt( i-1 ,j ,k ,um ,vm ,wm ,tm )
       ukc = 0.5 * ( up - um )
       vkc = 0.5 * ( vp - vm )
       wkc = 0.5 * ( wp - wm )
       tkc = 0.5 * ( tp - tm )
    endif

    if ( j == 1 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1 )
        call getuvwt( i ,2 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,3 ,k ,u3 ,v3 ,w3 ,t3 )
        uet = - ( aa * u1 + bb * u2 + cc * u3 )
        vet = - ( aa * v1 + bb * v2 + cc * v3 )
        wet = - ( aa * w1 + bb * w2 + cc * w3 )
        tet = - ( aa * t1 + bb * t2 + cc * t3 )

    elseif ( j == nj ) then
        call getuvwt( i ,j   ,k ,u1 ,v1 ,w1 ,t1 )
        call getuvwt( i ,j-1 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j-2 ,k ,u3 ,v3 ,w3 ,t3 )
        uet = aa * u1 + bb * u2 + cc * u3
        vet = aa * v1 + bb * v2 + cc * v3
        wet = aa * w1 + bb * w2 + cc * w3
        tet = aa * t1 + bb * t2 + cc * t3
    else
        call getuvwt( i ,j+1 ,k ,up ,vp ,wp ,tp )
        call getuvwt( i ,j-1 ,k ,um ,vm ,wm ,tm )
        uet = 0.5 * ( up - um )
        vet = 0.5 * ( vp - vm )
        wet = 0.5 * ( wp - wm )
        tet = 0.5 * ( tp - tm )
    endif

    if ( k == 1 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1 )
        call getuvwt( i ,j ,2 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,3 ,u3 ,v3 ,w3 ,t3 )
        uct = - ( aa * u1 + bb * u2 + cc * u3 )
        vct = - ( aa * v1 + bb * v2 + cc * v3 )
        wct = - ( aa * w1 + bb * w2 + cc * w3 )
        tct = - ( aa * t1 + bb * t2 + cc * t3 )

!			 endif

    elseif ( k == nk ) then
        call getuvwt( i ,j ,k   ,u1 ,v1 ,w1 ,t1 )
        call getuvwt( i ,j ,k-1 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,k-2 ,u3 ,v3 ,w3 ,t3 )
        uct = aa * u1 + bb * u2 + cc * u3
        vct = aa * v1 + bb * v2 + cc * v3
        wct = aa * w1 + bb * w2 + cc * w3
        tct = aa * t1 + bb * t2 + cc * t3

    else
        call getuvwt( i ,j ,k+1 ,up ,vp ,wp ,tp )
        call getuvwt( i ,j ,k-1 ,um ,vm ,wm ,tm )
        uct = 0.5 * ( up - um )
        vct = 0.5 * ( vp - vm )
        wct = 0.5 * ( wp - wm )
        tct = 0.5 * ( tp - tm )
    endif

    return
end subroutine getuvwtder
!_____________________________________________________________________!
subroutine getuvwt(i,j,k,um,vm,wm,tm)
    use global_variables
    implicit none

    integer :: i,j,k
    real :: um,vm,wm,tm
    um = u(i,j,k)
    vm = v(i,j,k)
    wm = w(i,j,k)
    tm = t(i,j,k)
    return
end subroutine getuvwt
!_____________________________________________________________________!
subroutine corner_point
    use global_variables
    implicit none
    integer :: m,mi,n1,n2
    integer :: i,j,k

    !先赋j,k面的角点
!    do i = 1,ni
    do i = 0,ni+1    
       do m=1,nl
          q(m,i,0 ,0 ) = 0.5 * ( q(m,i,1   ,0 ) + q(m,i,0 ,1   ) )
          q(m,i,0 ,nk) = 0.5 * ( q(m,i,1   ,nk) + q(m,i,0 ,nk-1) )
          q(m,i,nj,nk) = 0.5 * ( q(m,i,nj-1,nk) + q(m,i,nj,nk-1) )
          q(m,i,nj,0 ) = 0.5 * ( q(m,i,nj-1,0 ) + q(m,i,nj,1   ) )
       enddo
          
	   if(nl>5) then
       do m=1,nl-4
          fs(m,i,0 ,0 ) = 0.5 * ( fs(m,i,1   ,0 ) + fs(m,i,0 ,1   ) )
          fs(m,i,0 ,nk) = 0.5 * ( fs(m,i,1   ,nk) + fs(m,i,0 ,nk-1) )
          fs(m,i,nj,nk) = 0.5 * ( fs(m,i,nj-1,nk) + fs(m,i,nj,nk-1) )
          fs(m,i,nj,0 ) = 0.5 * ( fs(m,i,nj-1,0 ) + fs(m,i,nj,1   ) )
       enddo
	   endif

       r(i,0 ,0 ) = 0.5 * ( r(i,1   ,0 ) + r(i,0 ,1   ) )
       r(i,0 ,nk) = 0.5 * ( r(i,1   ,nk) + r(i,0 ,nk-1) )
       r(i,nj,nk) = 0.5 * ( r(i,nj-1,nk) + r(i,nj,nk-1) )
       r(i,nj,0 ) = 0.5 * ( r(i,nj-1,0 ) + r(i,nj,1   ) )

       u(i,0 ,0 ) = 0.5 * ( u(i,1   ,0 ) + u(i,0 ,1   ) )
       u(i,0 ,nk) = 0.5 * ( u(i,1   ,nk) + u(i,0 ,nk-1) )
       u(i,nj,nk) = 0.5 * ( u(i,nj-1,nk) + u(i,nj,nk-1) )
       u(i,nj,0 ) = 0.5 * ( u(i,nj-1,0 ) + u(i,nj,1   ) )

       v(i,0 ,0 ) = 0.5 * ( v(i,1   ,0 ) + v(i,0 ,1   ) )
       v(i,0 ,nk) = 0.5 * ( v(i,1   ,nk) + v(i,0 ,nk-1) )
       v(i,nj,nk) = 0.5 * ( v(i,nj-1,nk) + v(i,nj,nk-1) )
       v(i,nj,0 ) = 0.5 * ( v(i,nj-1,0 ) + v(i,nj,1   ) )

       w(i,0 ,0 ) = 0.5 * ( w(i,1   ,0 ) + w(i,0 ,1   ) )
       w(i,0 ,nk) = 0.5 * ( w(i,1   ,nk) + w(i,0 ,nk-1) )
       w(i,nj,nk) = 0.5 * ( w(i,nj-1,nk) + w(i,nj,nk-1) )
       w(i,nj,0 ) = 0.5 * ( w(i,nj-1,0 ) + w(i,nj,1   ) )

       p(i,0 ,0 ) = 0.5 * ( p(i,1   ,0 ) + p(i,0 ,1   ) )
       p(i,0 ,nk) = 0.5 * ( p(i,1   ,nk) + p(i,0 ,nk-1) )
       p(i,nj,nk) = 0.5 * ( p(i,nj-1,nk) + p(i,nj,nk-1) )
       p(i,nj,0 ) = 0.5 * ( p(i,nj-1,0 ) + p(i,nj,1   ) )

       t(i,0 ,0 ) = 0.5 * ( t(i,1   ,0 ) + t(i,0 ,1   ) )
       t(i,0 ,nk) = 0.5 * ( t(i,1   ,nk) + t(i,0 ,nk-1) )
       t(i,nj,nk) = 0.5 * ( t(i,nj-1,nk) + t(i,nj,nk-1) )
       t(i,nj,0 ) = 0.5 * ( t(i,nj-1,0 ) + t(i,nj,1   ) )

!	   gm(i,0 ,0 ) = 0.5 * ( gm(i,1   ,0 ) + gm(i,0 ,1   ) )
!	   gm(i,0 ,nk) = 0.5 * ( gm(i,1   ,nk) + gm(i,0 ,nk-1) )
!       gm(i,nj,nk) = 0.5 * ( gm(i,nj-1,nk) + gm(i,nj,nk-1) )
!       gm(i,nj,0 ) = 0.5 * ( gm(i,nj-1,0 ) + gm(i,nj,1   ) )

    enddo

    !赋i,k面的角点
    do j=0,nj+1
       do m=1,nm
          q(m,0 ,j,0 ) = 0.5 * ( q(m,1   ,j,0 ) + q(m,0 ,j,1   ) )
          q(m,0 ,j,nk) = 0.5 * ( q(m,1   ,j,nk) + q(m,0 ,j,nk-1) )
          q(m,ni,j,nk) = 0.5 * ( q(m,ni-1,j,nk) + q(m,ni,j,nk-1) )
          q(m,ni,j,0 ) = 0.5 * ( q(m,ni-1,j,0 ) + q(m,ni,j,1   ) )
       enddo
	   if(nl>5) then
       do m=1,nl-4
          fs(m,0 ,j,0 ) = 0.5 * ( fs(m,1   ,j,0 ) + fs(m,0 ,j,1   ) )
          fs(m,0 ,j,nk) = 0.5 * ( fs(m,1   ,j,nk) + fs(m,0 ,j,nk-1) )
          fs(m,ni,j,nk) = 0.5 * ( fs(m,ni-1,j,nk) + fs(m,ni,j,nk-1) )
          fs(m,ni,j,0 ) = 0.5 * ( fs(m,ni-1,j,0 ) + fs(m,ni,j,1   ) )
       enddo
	   endif
       r(0 ,j,0 ) = 0.5 * ( r(1   ,j,0 ) + r(0 ,j,1   ) )
       r(0 ,j,nk) = 0.5 * ( r(1   ,j,nk) + r(0 ,j,nk-1) )
       r(ni,j,nk) = 0.5 * ( r(ni-1,j,nk) + r(ni,j,nk-1) )
       r(ni,j,0 ) = 0.5 * ( r(ni-1,j,0 ) + r(ni,j,1   ) )

       u(0 ,j,0 ) = 0.5 * ( u(1   ,j,0 ) + u(0 ,j,1   ) )
       u(0 ,j,nk) = 0.5 * ( u(1   ,j,nk) + u(0 ,j,nk-1) )
       u(ni,j,nk) = 0.5 * ( u(ni-1,j,nk) + u(ni,j,nk-1) )
       u(ni,j,0 ) = 0.5 * ( u(ni-1,j,0 ) + u(ni,j,1   ) )

       v(0 ,j,0 ) = 0.5 * ( v(1   ,j,0 ) + v(0 ,j,1   ) )
       v(0 ,j,nk) = 0.5 * ( v(1   ,j,nk) + v(0 ,j,nk-1) )
       v(ni,j,nk) = 0.5 * ( v(ni-1,j,nk) + v(ni,j,nk-1) )
       v(ni,j,0 ) = 0.5 * ( v(ni-1,j,0 ) + v(ni,j,1   ) )

       w(0 ,j,0 ) = 0.5 * ( w(1   ,j,0 ) + w(0 ,j,1   ) )
       w(0 ,j,nk) = 0.5 * ( w(1   ,j,nk) + w(0 ,j,nk-1) )
       w(ni,j,nk) = 0.5 * ( w(ni-1,j,nk) + w(ni,j,nk-1) )
       w(ni,j,0 ) = 0.5 * ( w(ni-1,j,0 ) + w(ni,j,1   ) )

       p(0 ,j,0 ) = 0.5 * ( p(1   ,j,0 ) + p(0 ,j,1   ) )
       p(0 ,j,nk) = 0.5 * ( p(1   ,j,nk) + p(0 ,j,nk-1) )
       p(ni,j,nk) = 0.5 * ( p(ni-1,j,nk) + p(ni,j,nk-1) )
       p(ni,j,0 ) = 0.5 * ( p(ni-1,j,0 ) + p(ni,j,1   ) )

       t(0 ,j,0 ) = 0.5 * ( t(1   ,j,0 ) + t(0 ,j,1   ) )
       t(0 ,j,nk) = 0.5 * ( t(1   ,j,nk) + t(0 ,j,nk-1) )
       t(ni,j,nk) = 0.5 * ( t(ni-1,j,nk) + t(ni,j,nk-1) )
       t(ni,j,0 ) = 0.5 * ( t(ni-1,j,0 ) + t(ni,j,1   ) )

!       gm(0 ,j,0 ) = 0.5 * ( gm(1   ,j,0 ) + gm(0 ,j,1   ) )
!       gm(0 ,j,nk) = 0.5 * ( gm(1   ,j,nk) + gm(0 ,j,nk-1) )
!       gm(ni,j,nk) = 0.5 * ( gm(ni-1,j,nk) + gm(ni,j,nk-1) )
!       gm(ni,j,0 ) = 0.5 * ( gm(ni-1,j,0 ) + gm(ni,j,1   ) )
enddo

    !赋i,j面的角点
    do k=0,nk+1
       do m=1,nm
          q(m,0 ,0 ,k) = 0.5 * ( q(m,1   ,0 ,k) + q(m,0 ,1   ,k) )
          q(m,0 ,nj,k) = 0.5 * ( q(m,1   ,nj,k) + q(m,0 ,nj-1,k) )
          q(m,ni,nj,k) = 0.5 * ( q(m,ni-1,nj,k) + q(m,ni,nj-1,k) )
          q(m,ni,0 ,k) = 0.5 * ( q(m,ni-1,0 ,k) + q(m,ni,1   ,k) )
       enddo
	   if(nl>5) then
       do m=1,nl-4
          fs(m,0 ,0 ,k) = 0.5 * ( fs(m,1   ,0 ,k) + fs(m,0 ,1   ,k) )
          fs(m,0 ,nj,k) = 0.5 * ( fs(m,1   ,nj,k) + fs(m,0 ,nj-1,k) )
          fs(m,ni,nj,k) = 0.5 * ( fs(m,ni-1,nj,k) + fs(m,ni,nj-1,k) )
          fs(m,ni,0 ,k) = 0.5 * ( fs(m,ni-1,0 ,k) + fs(m,ni,1   ,k) )
       enddo
	   endif
       r(0 ,0 ,k) = 0.5 * ( r(1   ,0 ,k) + r(0 ,1   ,k) )
       r(0 ,nj,k) = 0.5 * ( r(1   ,nj,k) + r(0 ,nj-1,k) )
       r(ni,nj,k) = 0.5 * ( r(ni-1,nj,k) + r(ni,nj-1,k) )
       r(ni,0 ,k) = 0.5 * ( r(ni-1,0 ,k) + r(ni,1   ,k) )

       u(0 ,0 ,k) = 0.5 * ( u(1   ,0 ,k) + u(0 ,1   ,k) )
       u(0 ,nj,k) = 0.5 * ( u(1   ,nj,k) + u(0 ,nj-1,k) )
       u(ni,nj,k) = 0.5 * ( u(ni-1,nj,k) + u(ni,nj-1,k) )
       u(ni,0 ,k) = 0.5 * ( u(ni-1,0 ,k) + u(ni,1   ,k) )

       v(0 ,0 ,k) = 0.5 * ( v(1   ,0 ,k) + v(0 ,1   ,k) )
       v(0 ,nj,k) = 0.5 * ( v(1   ,nj,k) + v(0 ,nj-1,k) )
       v(ni,nj,k) = 0.5 * ( v(ni-1,nj,k) + v(ni,nj-1,k) )
       v(ni,0 ,k) = 0.5 * ( v(ni-1,0 ,k) + v(ni,1   ,k) )

       w(0 ,0 ,k) = 0.5 * ( w(1   ,0 ,k) + w(0 ,1   ,k) )
       w(0 ,nj,k) = 0.5 * ( w(1   ,nj,k) + w(0 ,nj-1,k) )
       w(ni,nj,k) = 0.5 * ( w(ni-1,nj,k) + w(ni,nj-1,k) )
       w(ni,0 ,k) = 0.5 * ( w(ni-1,0 ,k) + w(ni,1   ,k) )

       p(0 ,0 ,k) = 0.5 * ( p(1   ,0 ,k) + p(0 ,1   ,k) )
       p(0 ,nj,k) = 0.5 * ( p(1   ,nj,k) + p(0 ,nj-1,k) )
       p(ni,nj,k) = 0.5 * ( p(ni-1,nj,k) + p(ni,nj-1,k) )
       p(ni,0 ,k) = 0.5 * ( p(ni-1,0 ,k) + p(ni,1   ,k) )

       t(0 ,0 ,k) = 0.5 * ( t(1   ,0 ,k) + t(0 ,1   ,k) )
       t(0 ,nj,k) = 0.5 * ( t(1   ,nj,k) + t(0 ,nj-1,k) )
       t(ni,nj,k) = 0.5 * ( t(ni-1,nj,k) + t(ni,nj-1,k) )
       t(ni,0 ,k) = 0.5 * ( t(ni-1,0 ,k) + t(ni,1   ,k) )

!       gm(0 ,0 ,k) = 0.5 * ( gm(1   ,0 ,k) + gm(0 ,1   ,k) )
!      gm(0 ,nj,k) = 0.5 * ( gm(1   ,nj,k) + gm(0 ,nj-1,k) )
!       gm(ni,nj,k) = 0.5 * ( gm(ni-1,nj,k) + gm(ni,nj-1,k) )
!       gm(ni,0 ,k) = 0.5 * ( gm(ni-1,0 ,k) + gm(ni,1   ,k) )
    enddo

    return
end subroutine corner_point
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine readline(id,nline)
    implicit none
    integer :: id,nline,n

    do n=1,nline
       read(id,*)
    enddo

    return
end subroutine readline
!_____________________________________________________________________!
subroutine getrkec_mml(i,j,k,nx,ny,nz,nt,i_mml,j_mml,k_mml)
    use global_variables,only : kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,kct,ett,ctt,ni,nj,nk
    implicit none
    integer :: i,j,k,m,n,i_mml,j_mml,k_mml,i1,j1,k1
    real :: nx,ny,nz,nt


!___cic 隐式边界时间积分
    i1=i
	j1=j
	k1=k
	if(i1==0   )i1=1
    if(i1==ni+1)i1=ni
    if(j1==0   )j1=1
    if(j1==nj+1)j1=nj
    if(k1==0   )k1=1
    if(k1==nk+1)k1=nk

    do m=1,i_mml
       nx = kcx(i1,j1,k1)
       ny = kcy(i1,j1,k1)
       nz = kcz(i1,j1,k1)          
       nt = kct(i1,j1,k1)          
    enddo

    do m=1,j_mml
       nx = etx(i1,j1,k1)
       ny = ety(i1,j1,k1)
       nz = etz(i1,j1,k1)          
       nt = ett(i1,j1,k1)          
    enddo

    do m=1,k_mml
       nx = ctx(i1,j1,k1)
       ny = cty(i1,j1,k1)
       nz = ctz(i1,j1,k1)          
       nt = ctt(i1,j1,k1)          
    enddo
!___cic 隐式边界时间积分

    return
end subroutine getrkec_mml
!_____________________________________________________________________!
subroutine getr0kec_mml(i,j,k,nx,ny,nz,nt,i_mml,j_mml,k_mml)
    use global_variables, &
	only : method,kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,kct,ett,ctt
    implicit none
    integer :: i,j,k,m,n,i_mml,j_mml,k_mml
    real :: nx,ny,nz,nt

    do m=1,i_mml
          nx = kcx(i,j,k)
          ny = kcy(i,j,k)
          nz = kcz(i,j,k)          
          nt = kct(i,j,k)          
    enddo

    do m=1,j_mml
          nx = etx(i,j,k)
          ny = ety(i,j,k)
          nz = etz(i,j,k)          
          nt = ett(i,j,k)          
    enddo

    do m=1,k_mml
          nx = ctx(i,j,k)
          ny = cty(i,j,k)
          nz = ctz(i,j,k)          
          nt = ctt(i,j,k)          
    enddo

    return
end subroutine getr0kec_mml
!_____________________________________________________________________!
subroutine getnxyz_mml(nx,ny,nz,i,j,k,is,js,ks)     !计算单位法向矢量
    use global_variables
    implicit none
    integer :: is,js,ks,i,j,k,n
    real :: nx,ny,nz,sss1
    nx = 0
    ny = 0
    nz = 0
    do n=1,abs(is)
       nx = nx + kcx(i,j,k)
       ny = ny + kcy(i,j,k)
       nz = nz + kcz(i,j,k)
    enddo

    do n=1,abs(js)
       nx = nx + etx(i,j,k)
       ny = ny + ety(i,j,k)
       nz = nz + etz(i,j,k)
    enddo

    do n=1,abs(ks)
       nx = nx + ctx(i,j,k)
       ny = ny + cty(i,j,k)
       nz = nz + ctz(i,j,k)
    enddo
    sss1 = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
    nx = nx / sss1
    ny = ny / sss1
    nz = nz / sss1
    return
end subroutine getnxyz_mml
!_____________________________________________________________________!
subroutine getnxyz_tgh(nx,ny,nz,i,j,k,is,js,ks)
     !* 与getnxyz_mml不同的地方在与本程序考虑Jacobian
    use global_variables
    implicit none
    integer :: is,js,ks,i,j,k,n
    real :: nx,ny,nz,sss1
    nx = 0
    ny = 0
    nz = 0
    do n=1,abs(is)
       nx = nx + kcx(i,j,k)/vol(i,j,k)
       ny = ny + kcy(i,j,k)/vol(i,j,k)
       nz = nz + kcz(i,j,k)/vol(i,j,k)
    enddo

    do n=1,abs(js)
       nx = nx + etx(i,j,k)/vol(i,j,k)
       ny = ny + ety(i,j,k)/vol(i,j,k)
       nz = nz + etz(i,j,k)/vol(i,j,k)
    enddo

    do n=1,abs(ks)
       nx = nx + ctx(i,j,k)/vol(i,j,k)
       ny = ny + cty(i,j,k)/vol(i,j,k)
       nz = nz + ctz(i,j,k)/vol(i,j,k)
    enddo
    sss1 = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
    nx = nx / sss1
    ny = ny / sss1
    nz = nz / sss1
    return
end subroutine getnxyz_tgh
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine getrkec_tgh(i,j,k,nx,ny,nz,nt,i_mml,j_mml,k_mml,ni,nj,nk)
    use global_variables,only : kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,kct,ett,ctt
    implicit none
    integer :: i,j,k,m,n,i_mml,j_mml,k_mml,i1,j1,k1,ni,nj,nk
    real :: nx,ny,nz,nt

!___cic 隐式边界时间积分
    i1=i
	j1=j
	k1=k
	if(i1==0   )i1=1
    if(i1==ni+1)i1=ni
    if(j1==0   )j1=1
    if(j1==nj+1)j1=nj
    if(k1==0   )k1=1
    if(k1==nk+1)k1=nk

    do m=1,i_mml
       nx = kcx(i1,j1,k1)
       ny = kcy(i1,j1,k1)
       nz = kcz(i1,j1,k1)          
       nt = kct(i1,j1,k1)          
    enddo

    do m=1,j_mml
       nx = etx(i1,j1,k1)
       ny = ety(i1,j1,k1)
       nz = etz(i1,j1,k1)          
       nt = ett(i1,j1,k1)          
    enddo

    do m=1,k_mml
       nx = ctx(i1,j1,k1)
       ny = cty(i1,j1,k1)
       nz = ctz(i1,j1,k1)          
       nt = ctt(i1,j1,k1)          
    enddo
!___cic 隐式边界时间积分

    return
end subroutine getrkec_tgh

subroutine aver3(a,b,c,d)
    implicit none 
    real :: a,b,c,d
    d = ( a + b + c )/3.0
    return
end subroutine aver3
!_____________________________________________________________________!
subroutine halen(a,b,c,s)
    implicit none 
    real :: a,b,c,s,p
    p = 0.5 * ( a + b + c )
	S=p * ( p - a ) * ( p - b ) * ( p - c )
	IF(S <0 )THEN
	    WRITE(*,*)'求三角形面积出错，可能是网格的纵横比太大'
	    WRITE(*,*)'求三角形的三条边长：',A,B,C
		S=0.0
		RETURN
	ENDIF
    s = sqrt( S )
    return
end subroutine halen
!_____________________________________________________________________!
subroutine lenth(x1,y1,z1,x2,y2,z2,r)
    implicit none 
    real :: x1,y1,z1,x2,y2,z2,r
    real :: dx,dy,dz
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    r = sqrt( dx*dx + dy*dy + dz*dz )
    return
end subroutine lenth
!_____________________________________________________________________!
subroutine getuvwtder_vir(i,j,k)
    use global_variables,only : u,v,w,T,r,ni,nj,nk
    use duvwt_module
    implicit none

    integer :: i,j,k
!    real :: up,vp,wp,tp,um,vm,wm,tm
!    real :: u1,v1,w1,t1,u2,v2,w2,t2,u3,v3,w3,t3,up,vp,wp,tp,um,vm,wm,tm

		if(i==1 .and. r(-2,j,k)<0.0)then
			ukc = -1.5*u(i,j,k) + 2.0*u(i+1,j,k) - 0.5*u(i+2,j,k)
			vkc = -1.5*v(i,j,k) + 2.0*v(i+1,j,k) - 0.5*v(i+2,j,k)
			wkc = -1.5*w(i,j,k) + 2.0*w(i+1,j,k) - 0.5*w(i+2,j,k)
			tkc = -1.5*t(i,j,k) + 2.0*t(i+1,j,k) - 0.5*t(i+2,j,k)
		elseif(i==ni .and. r(ni+3,j,k)<0.0)then
			ukc =  1.5*u(i,j,k) - 2.0*u(i-1,j,k) + 0.5*u(i-2,j,k)
			vkc =  1.5*v(i,j,k) - 2.0*v(i-1,j,k) + 0.5*v(i-2,j,k)
			wkc =  1.5*w(i,j,k) - 2.0*w(i-1,j,k) + 0.5*w(i-2,j,k)
			tkc =  1.5*t(i,j,k) - 2.0*t(i-1,j,k) + 0.5*t(i-2,j,k)
		else
			ukc = 0.5 * ( u(i+1,j,k) - u(i-1,j,k) )
			vkc = 0.5 * ( v(i+1,j,k) - v(i-1,j,k) )
			wkc = 0.5 * ( w(i+1,j,k) - w(i-1,j,k) )
			tkc = 0.5 * ( t(i+1,j,k) - t(i-1,j,k) )
		endif

		if(j==1 .and. r(i,-2,k)<0.0)then
			uet = -1.5*u(i,j,k) + 2.0*u(i,j+1,k) - 0.5*u(i,j+2,k)
			vet = -1.5*v(i,j,k) + 2.0*v(i,j+1,k) - 0.5*v(i,j+2,k)
			wet = -1.5*w(i,j,k) + 2.0*w(i,j+1,k) - 0.5*w(i,j+2,k)
			tet = -1.5*t(i,j,k) + 2.0*t(i,j+1,k) - 0.5*t(i,j+2,k)
		elseif(j==nj .and. r(i,nj+3,k)<0.0)then
			uet =  1.5*u(i,j,k) - 2.0*u(i,j-1,k) + 0.5*u(i,j-2,k)
			vet =  1.5*v(i,j,k) - 2.0*v(i,j-1,k) + 0.5*v(i,j-2,k)
			wet =  1.5*w(i,j,k) - 2.0*w(i,j-1,k) + 0.5*w(i,j-2,k)
			tet =  1.5*t(i,j,k) - 2.0*t(i,j-1,k) + 0.5*t(i,j-2,k)
		else
			uet = 0.5 * ( u(i,j+1,k) - u(i,j-1,k) )
			vet = 0.5 * ( v(i,j+1,k) - v(i,j-1,k) )
			wet = 0.5 * ( w(i,j+1,k) - w(i,j-1,k) )
			tet = 0.5 * ( t(i,j+1,k) - t(i,j-1,k) )
		endif

		if(k==1 .and. r(i,j,-2)<0.0)then
			uct = -1.5*u(i,j,k) + 2.0*u(i,j,k+1) - 0.5*u(i,j,k+2)
			vct = -1.5*v(i,j,k) + 2.0*v(i,j,k+1) - 0.5*v(i,j,k+2)
			wct = -1.5*w(i,j,k) + 2.0*w(i,j,k+1) - 0.5*w(i,j,k+2)
			tct = -1.5*t(i,j,k) + 2.0*t(i,j,k+1) - 0.5*t(i,j,k+2)
		elseif(k==nk .and. r(i,j,nk+3)<0.0)then
			uct =  1.5*u(i,j,k) - 2.0*u(i,j,k-1) + 0.5*u(i,j,k-2)
			vct =  1.5*v(i,j,k) - 2.0*v(i,j,k-1) + 0.5*v(i,j,k-2)
			wct =  1.5*w(i,j,k) - 2.0*w(i,j,k-1) + 0.5*w(i,j,k-2)
			tct =  1.5*t(i,j,k) - 2.0*t(i,j,k-1) + 0.5*t(i,j,k-2)
		else
			uct = 0.5 * ( u(i,j,k+1) - u(i,j,k-1) )
			vct = 0.5 * ( v(i,j,k+1) - v(i,j,k-1) )
			wct = 0.5 * ( w(i,j,k+1) - w(i,j,k-1) )
			tct = 0.5 * ( t(i,j,k+1) - t(i,j,k-1) )
		endif
    return
end subroutine getuvwtder_vir
!_____________________________________________________________________!

