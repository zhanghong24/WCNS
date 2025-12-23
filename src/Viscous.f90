!_____________________________________________________________________!
subroutine spectvis
    use global_variables
    implicit none
    integer :: i,j,k,cmethod
    real :: rm_vol,ri,rj,rk
    real :: kx,ky,kz,ex,ey,ez,cx,cy,cz,nt
    real :: coef,vis,coef1

    do k=1,nk
       do j=1,nj
          do i=1,ni
             rm_vol = r(i,j,k)* vol(i,j,k) + small
             call getr0kec_mml(i,j,k,kx,ky,kz,nt,1,0,0)
             call getr0kec_mml(i,j,k,ex,ey,ez,nt,0,1,0)
             call getr0kec_mml(i,j,k,cx,cy,cz,nt,0,0,1)

             ri =  kx*kx + ky*ky + kz*kz
             rj =  ex*ex + ey*ey + ez*ez
             rk =  cx*cx + cy*cy + cz*cz
             vis = visl(i,j,k) + vist(i,j,k)
             coef = 2.0*vis/( reynolds * rm_vol  )
             coef1 = csrv*coef
             srva(i,j,k) = coef1*ri
             srvb(i,j,k) = coef1*rj
             srvc(i,j,k) = coef1*rk
          enddo
       enddo
    enddo

    return
end subroutine spectvis
!_____________________________________________________________________!
subroutine viscode
    use global_variables
    implicit none
       call vis_dif
    return
end subroutine viscode
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine vis_dif
!*TGH. 为了便于使用特征边界
!*TGH. 把求解结果区域从2 ~ max-1
!*TGH.  扩展为 1 ~ max
!*TGH. 计算结果保持在dq中
!
!*TGH. Modified by TU Guohua, 2009.2
!
!*TGH. 若某个方向的网格数小于等于nijk2d，则关闭该方向的离散（二维）
!*TGH. 还待完善!
!--------------------------------------------------------------
    use global_variables,only : dq,reynolds,nl,ni,nj,nk,vol,nmax,nijk2d
    implicit none
    integer :: i,j,k,m
    real :: ev_l(nl),re
    real,pointer,dimension( :,: ) :: fv

    allocate( fv(nl,nmax+1) )

 !          call set_dq_to_0

   re = 1.0 / reynolds

!----------------------------------------------
! Modified by TU Guohua
     if(ni <= nijk2d)goto 11  !*TGH. 二维
     do k=1,nk
       do j=1,nj
! end Modified
!----------------------------------------------
          !计算出一根网格线上的粘性通量


          do i =1,ni
             if ( vol(i,j,k) > 0.0 ) then
                call getev_dif(i,j,k,ev_l)
             else
                ev_l(:) = 0.0
             endif
             do m=1,nl
                fv(m,i) = ev_l(m)
             enddo
          enddo

          do i=2,ni-1
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) - 0.5 * re * ( fv(m,i+1) - fv(m,i-1) )
             enddo
          enddo
	      
         !*TGH. 补充边界上的单侧差分
		  do m=2,nl
			dq(m,1,j,k)  = dq(m,1,j,k)  - re* ( fv(m,2)  - fv(m,1)   )
			dq(m,ni,j,k) = dq(m,ni,j,k) - re* ( fv(m,ni) - fv(m,ni-1) )
!			dq(m,1,j,k)  = dq(m,1,j,k)  - 0.5*re* (-3.*fv(m,1) + 4.*fv(m,2) - fv(m,3))
!			dq(m,ni,j,k) = dq(m,ni,j,k) - 0.5*re* (3.*fv(m,ni) - 4.*fv(m,ni-1) + fv(m,ni-2))
		  enddo
		  !*TGH. end 补充边界上的单侧差分


       enddo
    enddo
11  continue

!----------------------------------------------
! Modified by TU Guohua

! Modified
    if(nj <= nijk2d)goto 12  !*TGH. 二维
    do k=1,nk
       do i=1,ni
! end Modified
!----------------------------------------------
          !计算出一根网格线上的粘性通量
 

          do j = 1,nj
             if ( vol(i,j,k) > 0.0 ) then
                call getfv_dif(i,j,k,ev_l)
             else
                ev_l(:) = 0.0
             endif
             do m=1,nl
                fv(m,j) = ev_l(m)
             enddo
          enddo

          do j=2,nj-1
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) - 0.5 * re * ( fv(m,j+1) - fv(m,j-1) )
             enddo
          enddo
         !*TGH. 补充边界上的单侧差分
		  do m=2,nl
			dq(m,i,1,k)  = dq(m,i,1,k) -  re* ( fv(m,2)  - fv(m,1)    )
			dq(m,i,nj,k) = dq(m,i,nj,k) - re* ( fv(m,nj) - fv(m,nj-1) )
!			dq(m,i,1,k)  = dq(m,i,1,k) -  0.5*re* (-3.*fv(m,1) + 4.*fv(m,2) - fv(m,3))
!			dq(m,i,nj,k) = dq(m,i,nj,k) - 0.5*re* ( 3.*fv(m,nj) - 4.*fv(m,nj-1) + fv(m,nj-2))
		  enddo
		  !*TGH. end 补充边界上的单侧差分


       enddo
    enddo
12  continue

!----------------------------------------------
! Modified by TU Guohua
! Modified
    if(nk <= nijk2d)goto 13  !*TGH. 二维

    do i=1,ni
       do j=1,nj
! end Modified
!----------------------------------------------
          !计算出一根网格线上的粘性通量

          do k = 1,nk
             if ( vol(i,j,k) > 0.0 ) then
                call getgv_dif(i,j,k,ev_l) !*TGH 计算一线上的粘性通量
             else
                ev_l(:) = 0.0
             endif
             do m=1,nl
                fv(m,k) = ev_l(m)
             enddo
          enddo

          do k=2,nk-1
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) - 0.5*re* (fv(m,k+1) - fv(m,k-1))
             enddo
          enddo
          !*TGH. 补充边界上的单侧差分
		  do m=2,nl
			dq(m,i,j,1)  = dq(m,i,j,1)  - re* ( fv(m,2)  - fv(m,1)    )
			dq(m,i,j,nk) = dq(m,i,j,nk) - re* ( fv(m,nk) - fv(m,nk-1) )
!			dq(m,i,j,1)  = dq(m,i,j,1)  - 0.5*re* (-3.*fv(m,1) + 4.*fv(m,2) - fv(m,3))
!			dq(m,i,j,nk) = dq(m,i,j,nk) - 0.5*re* ( 3.*fv(m,nk) - 4.*fv(m,nk-1) + fv(m,nk-2))
		  enddo
		  !*TGH. end 补充边界上的单侧差分

      enddo
    enddo
13  continue

    deallocate( fv )

    return
end subroutine vis_dif
!*TGH
!_____________________________________________________________________!
subroutine getev_dif(i,j,k,ev)
    use global_variables
    use stress_module
    use geometry_moduLE
    implicit none
    integer :: i,j,k,m
    real :: vis
    real :: cp,kcp
    real :: vis_l,vis_t,co_rdi
    real :: dtdx,dtdy,dtdz
    real :: um,vm,wm,qx,qy,qz
    real :: cs_init(ns),xi(ns),rdi(ns),hs(ns)
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: ev(nl)
    real :: temp

    !以下各程序的调用次序不可随意改动
    vis_l = visl(i,j,k)
    vis_t = vist(i,j,k)
    vis = vis_l + vis_t
    ev = 0.0

    call getgeo(i,j,k)
    call getuvwtder(i,j,k)
    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                   dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)
    call stress(vis)



    cp = 1.0/((gama-1.0)*moo*moo)

    kcp = ( vis_l/prl + vis_t/prt ) * cp

    !qx,qy,qz置零
!    qx = 0.0
!    qy = 0.0
!    qz = 0.0
!    qx = qx + kcp * dtdx
!    qy = qy + kcp * dtdy
!    qz = qz + kcp * dtdz

    qx = kcp * dtdx
    qy = kcp * dtdy
    qz = kcp * dtdz

    !这个地方要注意,不要忘掉了
    um = u(i,j,k)
    vm = v(i,j,k)
    wm = w(i,j,k)

    call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                       kx,ky,kz,qx,qy,qz,um,vm,wm,ev,nl)
    return
end subroutine getev_dif
!_____________________________________________________________________!
subroutine getfv_dif(i,j,k,fv)
    use global_variables
    use stress_module
    use geometry_moduLE
    implicit none
    integer :: i,j,k,m
    real :: vis
    real :: cp,kcp
    real :: vis_l,vis_t,co_rdi
    real :: dtdx,dtdy,dtdz
    real :: um,vm,wm,qx,qy,qz
    real :: cs_init(ns),xi(ns),rdi(ns),hs(ns)
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: fv(nl)
    real :: temp

    !以下各程序的调用次序不可随意改动
    vis_l = visl(i,j,k)
    vis_t = vist(i,j,k)
    vis = vis_l + vis_t

    call getgeo(i,j,k)
    call getuvwtder(i,j,k)
    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                   dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)
    call stress(vis)

    !qx,qy,qz置零
    qx = 0.0
    qy = 0.0
    qz = 0.0


    cp = 1.0/((gama-1.0)*moo*moo)

    kcp = ( vis_l/prl + vis_t/prt ) * cp

    qx = qx + kcp * dtdx
    qy = qy + kcp * dtdy
    qz = qz + kcp * dtdz

    !这个地方要注意,不要忘掉了
    um = u(i,j,k)
    vm = v(i,j,k)
    wm = w(i,j,k)

    call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                           ex,ey,ez,qx,qy,qz,um,vm,wm,fv,nl)

    return
end subroutine getfv_dif
!_____________________________________________________________________!
subroutine getgv_dif(i,j,k,gv)
    use global_variables
    use stress_module
    use geometry_moduLE
    implicit none
    integer :: i,j,k,m
    real :: vis
    real :: cp,kcp
    real :: vis_l,vis_t,co_rdi
    real :: dtdx,dtdy,dtdz
    real :: um,vm,wm,qx,qy,qz
    real :: cs_init(ns),xi(ns),rdi(ns),hs(ns)
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: gv(nl)
    real :: temp

    !以下各程序的调用次序不可随意改动
    vis_l = visl(i,j,k)
    vis_t = vist(i,j,k)
    vis = vis_l + vis_t

    call getgeo(i,j,k)
    call getuvwtder(i,j,k)
    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                    dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)
    call stress(vis)

    !qx,qy,qz置零
    qx = 0.0
    qy = 0.0
    qz = 0.0


    cp = 1.0/((gama-1.0)*moo*moo)
    kcp = ( vis_l/prl + vis_t/prt ) * cp

    qx = qx + kcp * dtdx
    qy = qy + kcp * dtdy
    qz = qz + kcp * dtdz

    !这个地方要注意,不要忘掉了
    um = u(i,j,k)
    vm = v(i,j,k)
    wm = w(i,j,k)

    call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                  cx,cy,cz,qx,qy,qz,um,vm,wm,gv,nl)
    return
end subroutine getgv_dif
!_____________________________________________________________________!
subroutine dfi_dxyz(i,j,k,dfidx,dfidy,dfidz)
    use global_const,only:ns
    use geometry_moduLE
    implicit none
    integer :: m,i,j,k
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: dfi_dkc(ns),dfi_det(ns),dfi_dct(ns)

    call getdfi_der(i,j,k,dfi_dkc,dfi_det,dfi_dct)

    do m=1,ns
       dfidx(m) = (kx*dfi_dkc(m) + ex*dfi_det(m) + cx*dfi_dct(m))/vjacob
       dfidy(m) = (ky*dfi_dkc(m) + ey*dfi_det(m) + cy*dfi_dct(m))/vjacob
       dfidz(m) = (kz*dfi_dkc(m) + ez*dfi_det(m) + cz*dfi_dct(m))/vjacob
    enddo

    return
end subroutine dfi_dxyz
!_____________________________________________________________________!
subroutine getdfi_der(i,j,k,dfi_dkc,dfi_det,dfi_dct)
    use global_variables
    use duvwt_module
    implicit none
    integer :: i,j,k
    real :: aa,bb,cc
    real :: u1,v1,w1,t1,u2,v2,w2,t2,u3,v3,w3,t3,up,vp,wp,tp,um,vm,wm,tm
    real :: fi_1(ns),fi_2(ns),fi_3(ns),fi_p(ns),fi_m(ns)
    real :: dfi_dkc(ns),dfi_det(ns),dfi_dct(ns)

    aa =  1.5
    bb = -2.0
    cc =  0.5

    fi_1(:) = fs(:,i,j,k)
    if ( i == 1 ) then
       fi_2(:) = fs(:,i+1,j,k)
       fi_3(:) = fs(:,i+2,j,k)
       dfi_dkc = - ( aa * fi_1 + bb * fi_2 + cc * fi_3 )
    elseif ( i == ni ) then
       fi_2(:) = fs(:,i-1,j,k)
       fi_3(:) = fs(:,i-2,j,k)
       dfi_dkc =   ( aa * fi_1 + bb * fi_2 + cc * fi_3 )
    else
       fi_p(:) = fs(:,i+1,j,k)
       fi_m(:) = fs(:,i-1,j,k)
       dfi_dkc = 0.5 * ( fi_p - fi_m )
    endif

    if ( j == 1 ) then
       fi_2(:) = fs(:,i,j+1,k)
       fi_3(:) = fs(:,i,j+2,k)
       dfi_det = - ( aa * fi_1 + bb * fi_2 + cc * fi_3 )
    elseif ( j == nj ) then
       fi_2(:) = fs(:,i,j-1,k)
       fi_3(:) = fs(:,i,j-2,k)
       dfi_det =   ( aa * fi_1 + bb * fi_2 + cc * fi_3 )
    else
       fi_p(:) = fs(:,i,j+1,k)
       fi_m(:) = fs(:,i,j-1,k)
       dfi_det = 0.5 * ( fi_p - fi_m )
    endif

    if ( k == 1 ) then
       fi_2(:) = fs(:,i,j,k+1)
       fi_3(:) = fs(:,i,j,k+2)
       dfi_dct = - ( aa * fi_1 + bb * fi_2 + cc * fi_3 )
    elseif ( k == nk ) then
       fi_2(:) = fs(:,i,j,k-1)
       fi_3(:) = fs(:,i,j,k-2)
       dfi_dct =   ( aa * fi_1 + bb * fi_2 + cc * fi_3 )
    else
       fi_p(:) = fs(:,i,j,k+1)
       fi_m(:) = fs(:,i,j,k-1)
       dfi_dct = 0.5 * ( fi_p - fi_m )
    endif

    return
end subroutine getdfi_der
!_____________________________________________________________________!
subroutine fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                   nx,ny,nz,qx,qy,qz,um,vm,wm,fv,nl)
    !--  此子程序计算粘性通量，未引入雷诺数
    !--  VCD 是viscid 的缩写
    integer :: nl
    real :: fv(nl)
    real :: txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
    real :: nx,ny,nz,qx,qy,qz,um,vm,wm

    fv(1) = 0.0
    fv(2) = nx*txx + ny*txy + nz*txz
    fv(3) = nx*tyx + ny*tyy + nz*tyz
    fv(4) = nx*tzx + ny*tzy + nz*tzz
    fv(5) = um*fv(2) + vm*fv(3) + wm*fv(4) + ( nx*qx + ny*qy + nz*qz )
            
    return
end subroutine fluxvcd_point
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine viscode_virtual
    use global_variables
    implicit none
       call vis_dif_virtual
    return
end subroutine viscode_virtual
!_____________________________________________________________________!
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine vis_dif_virtual
!*TGH. 为了便于使用特征边界
!*TGH. 把求解结果区域从2 ~ max-1
!*TGH.  扩展为 1 ~ max
!*TGH. 计算结果保持在dq中
!
!*TGH. Modified by TU Guohua, 2009.2
!
!*TGH. 若某个方向的网格数小于等于nijk2d，则关闭该方向的离散（二维）
!*TGH. 还待完善!
!--------------------------------------------------------------
    use global_variables,only : dq,reynolds,nl,ni,nj,nk,vol,nmax,nijk2d
    implicit none
    integer :: i,j,k,m
    real :: ev_l(nl),re
    real,pointer,dimension( :,: ) :: fv
    
    allocate( fv(nl,nmax+1) )

 !          call set_dq_to_0

   re = 1.0 / reynolds

!----------------------------------------------
! Modified by TU Guohua
    if(ni <= nijk2d)goto 11  !*TGH. 二维
     do k=1,nk
       do j=1,nj
! end Modified
!----------------------------------------------
          !计算出一根网格线上的粘性通量


          do i =1,ni
             if ( vol(i,j,k) > 0.0 ) then
                call getev_dif_vir(i,j,k,ev_l)
             else
                ev_l(:) = 0.0
             endif
             do m=1,nl
                fv(m,i) = ev_l(m)
             enddo
          enddo

          do i=2,ni-1
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) - 0.5 * re * ( fv(m,i+1) - fv(m,i-1) )
             enddo
          enddo
	      
         !*TGH. 补充边界上的单侧差分
		  do m=2,nl
			dq(m,1,j,k)  = dq(m,1,j,k)  - re* ( fv(m,2)  - fv(m,1)   )
			dq(m,ni,j,k) = dq(m,ni,j,k) - re* ( fv(m,ni) - fv(m,ni-1) )
!			dq(m,1,j,k)  = dq(m,1,j,k)  - 0.5*re* (-3.*fv(m,1) + 4.*fv(m,2) - fv(m,3))
!			dq(m,ni,j,k) = dq(m,ni,j,k) - 0.5*re* (3.*fv(m,ni) - 4.*fv(m,ni-1) + fv(m,ni-2))
		  enddo
		  !*TGH. end 补充边界上的单侧差分


       enddo
    enddo
11  continue

!----------------------------------------------
! Modified by TU Guohua

! Modified
    if(nj <= nijk2d)goto 12  !*TGH. 二维
    do k=1,nk
       do i=1,ni
! end Modified
!----------------------------------------------
          !计算出一根网格线上的粘性通量
 

          do j = 1,nj
             if ( vol(i,j,k) > 0.0 ) then
                call getfv_dif_vir(i,j,k,ev_l)
             else
                ev_l(:) = 0.0
             endif
             do m=1,nl
                fv(m,j) = ev_l(m)
             enddo
          enddo

          do j=2,nj-1
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) - 0.5 * re * ( fv(m,j+1) - fv(m,j-1) )
             enddo
          enddo
         !*TGH. 补充边界上的单侧差分
		  do m=2,nl
			dq(m,i,1,k)  = dq(m,i,1,k) -  re* ( fv(m,2)  - fv(m,1)    )
			dq(m,i,nj,k) = dq(m,i,nj,k) - re* ( fv(m,nj) - fv(m,nj-1) )
!			dq(m,i,1,k)  = dq(m,i,1,k) -  0.5*re* (-3.*fv(m,1) + 4.*fv(m,2) - fv(m,3))
!			dq(m,i,nj,k) = dq(m,i,nj,k) - 0.5*re* ( 3.*fv(m,nj) - 4.*fv(m,nj-1) + fv(m,nj-2))
		  enddo
		  !*TGH. end 补充边界上的单侧差分


       enddo
    enddo
12  continue

!----------------------------------------------
! Modified by TU Guohua
! Modified
    if(nk <= nijk2d)goto 13  !*TGH. 二维

    do i=1,ni
       do j=1,nj
! end Modified
!----------------------------------------------
          !计算出一根网格线上的粘性通量

          do k = 1,nk
             if ( vol(i,j,k) > 0.0 ) then
                call getgv_dif_vir(i,j,k,ev_l) !*TGH 计算一线上的粘性通量
             else
                ev_l(:) = 0.0
             endif
             do m=1,nl
                fv(m,k) = ev_l(m)
             enddo
          enddo

          do k=2,nk-1
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) - 0.5*re* (fv(m,k+1) - fv(m,k-1))
             enddo
          enddo
          !*TGH. 补充边界上的单侧差分
		  do m=2,nl
			dq(m,i,j,1)  = dq(m,i,j,1)  - re* ( fv(m,2)  - fv(m,1)    )
			dq(m,i,j,nk) = dq(m,i,j,nk) - re* ( fv(m,nk) - fv(m,nk-1) )
!			dq(m,i,j,1)  = dq(m,i,j,1)  - 0.5*re* (-3.*fv(m,1) + 4.*fv(m,2) - fv(m,3))
!			dq(m,i,j,nk) = dq(m,i,j,nk) - 0.5*re* ( 3.*fv(m,nk) - 4.*fv(m,nk-1) + fv(m,nk-2))
		  enddo
		  !*TGH. end 补充边界上的单侧差分

      enddo
    enddo
13  continue

    deallocate( fv )

    return
end subroutine vis_dif_virtual
!_____________________________________________________________________!
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine getev_dif_vir(i,j,k,ev)
    use global_variables
    use stress_module
    use geometry_moduLE
    implicit none
    integer :: i,j,k,m
    real :: vis
    real :: cp,kcp
    real :: vis_l,vis_t,co_rdi
    real :: dtdx,dtdy,dtdz
    real :: um,vm,wm,qx,qy,qz
    real :: cs_init(ns),xi(ns),rdi(ns),hs(ns)
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: ev(nl)
    real :: temp

    !以下各程序的调用次序不可随意改动
    vis_l = visl(i,j,k)
    vis_t = vist(i,j,k)
    vis = vis_l + vis_t
    ev = 0.0

    call getgeo(i,j,k)
    call getuvwtder_vir(i,j,k)
    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                   dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)
    call stress(vis)



    cp = 1.0/((gama-1.0)*moo*moo)

    kcp = ( vis_l/prl + vis_t/prt ) * cp

    !qx,qy,qz置零
!    qx = 0.0
!    qy = 0.0
!    qz = 0.0
!    qx = qx + kcp * dtdx
!    qy = qy + kcp * dtdy
!    qz = qz + kcp * dtdz

    qx = kcp * dtdx
    qy = kcp * dtdy
    qz = kcp * dtdz

    !这个地方要注意,不要忘掉了
    um = u(i,j,k)
    vm = v(i,j,k)
    wm = w(i,j,k)

    call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                       kx,ky,kz,qx,qy,qz,um,vm,wm,ev,nl)
    return
end subroutine getev_dif_vir
!_____________________________________________________________________!
subroutine getfv_dif_vir(i,j,k,fv)
    use global_variables
    use stress_module
    use geometry_moduLE
    implicit none
    integer :: i,j,k,m
    real :: vis
    real :: cp,kcp
    real :: vis_l,vis_t,co_rdi
    real :: dtdx,dtdy,dtdz
    real :: um,vm,wm,qx,qy,qz
    real :: cs_init(ns),xi(ns),rdi(ns),hs(ns)
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: fv(nl)
    real :: temp

    !以下各程序的调用次序不可随意改动
    vis_l = visl(i,j,k)
    vis_t = vist(i,j,k)
    vis = vis_l + vis_t

    call getgeo(i,j,k)
    call getuvwtder_vir(i,j,k)
    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                   dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)
    call stress(vis)

    !qx,qy,qz置零
    qx = 0.0
    qy = 0.0
    qz = 0.0


    cp = 1.0/((gama-1.0)*moo*moo)

    kcp = ( vis_l/prl + vis_t/prt ) * cp

    qx = qx + kcp * dtdx
    qy = qy + kcp * dtdy
    qz = qz + kcp * dtdz

    !这个地方要注意,不要忘掉了
    um = u(i,j,k)
    vm = v(i,j,k)
    wm = w(i,j,k)

    call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                           ex,ey,ez,qx,qy,qz,um,vm,wm,fv,nl)

    return
end subroutine getfv_dif_vir
!_____________________________________________________________________!
subroutine getgv_dif_vir(i,j,k,gv)
    use global_variables
    use stress_module
    use geometry_moduLE
    implicit none
    integer :: i,j,k,m
    real :: vis
    real :: cp,kcp
    real :: vis_l,vis_t,co_rdi
    real :: dtdx,dtdy,dtdz
    real :: um,vm,wm,qx,qy,qz
    real :: cs_init(ns),xi(ns),rdi(ns),hs(ns)
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: gv(nl)
    real :: temp

    !以下各程序的调用次序不可随意改动
    vis_l = visl(i,j,k)
    vis_t = vist(i,j,k)
    vis = vis_l + vis_t

    call getgeo(i,j,k)
    call getuvwtder_vir(i,j,k)
    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                    dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)
    call stress(vis)

    !qx,qy,qz置零
    qx = 0.0
    qy = 0.0
    qz = 0.0


    cp = 1.0/((gama-1.0)*moo*moo)
    kcp = ( vis_l/prl + vis_t/prt ) * cp

    qx = qx + kcp * dtdx
    qy = qy + kcp * dtdy
    qz = qz + kcp * dtdz

    !这个地方要注意,不要忘掉了
    um = u(i,j,k)
    vm = v(i,j,k)
    wm = w(i,j,k)

    call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                  cx,cy,cz,qx,qy,qz,um,vm,wm,gv,nl)
    return
end subroutine getgv_dif_vir


subroutine WCNSE5_VIS_weighted
    use define_precision_mod
    use global_variables,only : r,u,v,w,p,t,reynolds,visl,vist,dq,nvis,     &
                                kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol, &
                                nl,nmax,ni,nj,nk,gama,moo,prl,prt
    implicit none
    integer :: i,j,k,m
    real :: re,cp,cp_prl,cp_prt,vis,kcp,vis2p3
    real :: kx,ky,kz,ex,ey,ez,cx,cy,cz,oovol
    real :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real :: dtdx,dtdy,dtdz,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    real,pointer,dimension(:,:,:,:) :: duijk,dvijk,dwijk,dtijk
    real,pointer :: fv(:,:),dfv(:)
    
    re = 1.0/reynolds

    cp = 1.0/((gama-1.0)*moo*moo)
    cp_prl = cp/prl
    cp_prt = cp/prt
    
    allocate( duijk(ni,nj,nk,3) )
    call get_dvnode(ni,nj,nk,r,u,duijk)
    
    allocate( dvijk(ni,nj,nk,3) )
    call get_dvnode(ni,nj,nk,r,v,dvijk)
    
    allocate( dwijk(ni,nj,nk,3) )
    call get_dvnode(ni,nj,nk,r,w,dwijk)
    
    allocate( dtijk(ni,nj,nk,3) )
    call get_dtnode(ni,nj,nk,r,p,t,dtijk)
    
    do k=1,nk
    do j=1,nj
    do i=1,ni
       kx = kcx(i,j,k)
       ky = kcy(i,j,k)
       kz = kcz(i,j,k)
         
       ex = etx(i,j,k)
       ey = ety(i,j,k)
       ez = etz(i,j,k)
       
       cx = ctx(i,j,k)
       cy = cty(i,j,k)
       cz = ctz(i,j,k)
       
       oovol = 1.0/vol(i,j,k)
       
       dudx = ( kx*duijk(i,j,k,1) + ex*duijk(i,j,k,2) + cx*duijk(i,j,k,3) ) * oovol
       dudy = ( ky*duijk(i,j,k,1) + ey*duijk(i,j,k,2) + cy*duijk(i,j,k,3) ) * oovol
       dudz = ( kz*duijk(i,j,k,1) + ez*duijk(i,j,k,2) + cz*duijk(i,j,k,3) ) * oovol

       dvdx = ( kx*dvijk(i,j,k,1) + ex*dvijk(i,j,k,2) + cx*dvijk(i,j,k,3) ) * oovol
       dvdy = ( ky*dvijk(i,j,k,1) + ey*dvijk(i,j,k,2) + cy*dvijk(i,j,k,3) ) * oovol
       dvdz = ( kz*dvijk(i,j,k,1) + ez*dvijk(i,j,k,2) + cz*dvijk(i,j,k,3) ) * oovol
       
       dwdx = ( kx*dwijk(i,j,k,1) + ex*dwijk(i,j,k,2) + cx*dwijk(i,j,k,3) ) * oovol
       dwdy = ( ky*dwijk(i,j,k,1) + ey*dwijk(i,j,k,2) + cy*dwijk(i,j,k,3) ) * oovol
       dwdz = ( kz*dwijk(i,j,k,1) + ez*dwijk(i,j,k,2) + cz*dwijk(i,j,k,3) ) * oovol
       
       dtdx = ( kx*dtijk(i,j,k,1) + ex*dtijk(i,j,k,2) + cx*dtijk(i,j,k,3) ) * oovol
       dtdy = ( ky*dtijk(i,j,k,1) + ey*dtijk(i,j,k,2) + cy*dtijk(i,j,k,3) ) * oovol
       dtdz = ( kz*dtijk(i,j,k,1) + ez*dtijk(i,j,k,2) + cz*dtijk(i,j,k,3) ) * oovol    
       
       duijk(i,j,k,1) = dudx
       duijk(i,j,k,2) = dudy
       duijk(i,j,k,3) = dudz
       
       dvijk(i,j,k,1) = dvdx
       dvijk(i,j,k,2) = dvdy
       dvijk(i,j,k,3) = dvdz
       
       dwijk(i,j,k,1) = dwdx
       dwijk(i,j,k,2) = dwdy
       dwijk(i,j,k,3) = dwdz
       
       dtijk(i,j,k,1) = dtdx
       dtijk(i,j,k,2) = dtdy
       dtijk(i,j,k,3) = dtdz
    end do
    end do
    end do
    
    allocate( fv(ni,5), dfv(ni))
    do k=1,nk
    do j=1,nj
       do i=1,ni
          vis = visl(i,j,k)
          kcp = vis*cp_prl
          
          if (nvis > 0) then
             vis = vis + vist(i,j,k)
             kcp = kcp + vist(i,j,k)*cp_prt
          end if

          kx = kcx(i,j,k)
          ky = kcy(i,j,k)
          kz = kcz(i,j,k)
            
          dudx = duijk(i,j,k,1)
          dudy = duijk(i,j,k,2)
          dudz = duijk(i,j,k,3)
                 
          dvdx = dvijk(i,j,k,1)
          dvdy = dvijk(i,j,k,2)
          dvdz = dvijk(i,j,k,3)
                 
          dwdx = dwijk(i,j,k,1)
          dwdy = dwijk(i,j,k,2)
          dwdz = dwijk(i,j,k,3)
                 
          dtdx = dtijk(i,j,k,1)
          dtdy = dtijk(i,j,k,2)
          dtdz = dtijk(i,j,k,3)
          
          vis2p3 = 2.0*vis/3.0
          tauxx = vis2p3 * ( 2.0*dudx - dvdy - dwdz )
          tauyy = vis2p3 * ( 2.0*dvdy - dwdz - dudx )
		  tauzz = vis2p3 * ( 2.0*dwdz - dudx - dvdy )
          tauxy = vis * ( dudy + dvdx )
          tauxz = vis * ( dudz + dwdx )
          tauyz = vis * ( dvdz + dwdy )
          
          fv(i,1) = 0.0
          fv(i,2) = kx*tauxx + ky*tauxy + kz*tauxz
          fv(i,3) = kx*tauxy + ky*tauyy + kz*tauyz
          fv(i,4) = kx*tauxz + ky*tauyz + kz*tauzz
          fv(i,5) = u(i,j,k)*fv(i,2) + v(i,j,k)*fv(i,3) + w(i,j,k)*fv(i,4) &
                  + kcp*( kx*dtdx + ky*dtdy + kz*dtdz )
       end do
       
       do m=1,5
          call dfnode_via_node(ni,fv(:,m),dfv)

          do i=1,ni
             dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(i)
          end do
       end do
    end do
    end do
    deallocate( dfv, fv )
    
    allocate( fv(nj,5), dfv(nj))
    do k=1,nk
    do i=1,ni
       do j=1,nj
          vis = visl(i,j,k)
          kcp = vis*cp_prl
          
          if (nvis > 0) then
             vis = vis + vist(i,j,k)
             kcp = kcp + vist(i,j,k)*cp_prt
          end if

          ex = etx(i,j,k)
          ey = ety(i,j,k)
          ez = etz(i,j,k)
          
          dudx = duijk(i,j,k,1)
          dudy = duijk(i,j,k,2)
          dudz = duijk(i,j,k,3)
                 
          dvdx = dvijk(i,j,k,1)
          dvdy = dvijk(i,j,k,2)
          dvdz = dvijk(i,j,k,3)
                 
          dwdx = dwijk(i,j,k,1)
          dwdy = dwijk(i,j,k,2)
          dwdz = dwijk(i,j,k,3)
                 
          dtdx = dtijk(i,j,k,1)
          dtdy = dtijk(i,j,k,2)
          dtdz = dtijk(i,j,k,3)
          
          vis2p3 = 2.0*vis/3.0
          tauxx = vis2p3 * ( 2.0*dudx - dvdy - dwdz )
          tauyy = vis2p3 * ( 2.0*dvdy - dwdz - dudx )
		  tauzz = vis2p3 * ( 2.0*dwdz - dudx - dvdy )
          tauxy = vis * ( dudy + dvdx )
          tauxz = vis * ( dudz + dwdx )
          tauyz = vis * ( dvdz + dwdy )
          
          fv(j,1) = 0.0
          fv(j,2) = ex*tauxx + ey*tauxy + ez*tauxz
          fv(j,3) = ex*tauxy + ey*tauyy + ez*tauyz
          fv(j,4) = ex*tauxz + ey*tauyz + ez*tauzz
          fv(j,5) = u(i,j,k)*fv(j,2) + v(i,j,k)*fv(j,3) + w(i,j,k)*fv(j,4) &
                  + kcp*( ex*dtdx + ey*dtdy + ez*dtdz )
       end do
       
       do m=1,5
          call dfnode_via_node(nj,fv(:,m),dfv)

          do j=1,nj
             dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(j)
          end do
       end do
    end do
    end do
    deallocate( dfv, fv )
    
    allocate( fv(nk,5), dfv(nk))
    do j=1,nj
    do i=1,ni
       do k=1,nk
          vis = visl(i,j,k)
          kcp = vis*cp_prl
          
          if (nvis > 0) then
             vis = vis + vist(i,j,k)
             kcp = kcp + vist(i,j,k)*cp_prt
          end if

          cx = ctx(i,j,k)
          cy = cty(i,j,k)
          cz = ctz(i,j,k)
          
          dudx = duijk(i,j,k,1)
          dudy = duijk(i,j,k,2)
          dudz = duijk(i,j,k,3)
                 
          dvdx = dvijk(i,j,k,1)
          dvdy = dvijk(i,j,k,2)
          dvdz = dvijk(i,j,k,3)
                 
          dwdx = dwijk(i,j,k,1)
          dwdy = dwijk(i,j,k,2)
          dwdz = dwijk(i,j,k,3)
                 
          dtdx = dtijk(i,j,k,1)
          dtdy = dtijk(i,j,k,2)
          dtdz = dtijk(i,j,k,3)
          
          vis2p3 = 2.0*vis/3.0
          tauxx = vis2p3 * ( 2.0*dudx - dvdy - dwdz )
          tauyy = vis2p3 * ( 2.0*dvdy - dwdz - dudx )
		  tauzz = vis2p3 * ( 2.0*dwdz - dudx - dvdy )
          tauxy = vis * ( dudy + dvdx )
          tauxz = vis * ( dudz + dwdx )
          tauyz = vis * ( dvdz + dwdy )
          
          fv(k,1) = 0.0
          fv(k,2) = cx*tauxx + cy*tauxy + cz*tauxz
          fv(k,3) = cx*tauxy + cy*tauyy + cz*tauyz
          fv(k,4) = cx*tauxz + cy*tauyz + cz*tauzz
          fv(k,5) = u(i,j,k)*fv(k,2) + v(i,j,k)*fv(k,3) + w(i,j,k)*fv(k,4) &
                  + kcp*( cx*dtdx + cy*dtdy + cz*dtdz )
       end do
       
       do m=1,5
          call dfnode_via_node(nk,fv(:,m),dfv)

          do k=1,nk
             dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(k)
          end do
       end do
    end do
    end do
    deallocate( dfv, fv )

    deallocate( dtijk, dwijk, dvijk, duijk )

end subroutine WCNSE5_VIS_weighted


subroutine get_dvnode(ni,nj,nk,ro,vn,dvn)
    use global_variables,only: nijk2nd,small
    implicit none
    integer :: ni,nj,nk
    real :: ro(-2:ni+3,-2:nj+3,-2:nk+3)
    real :: vn(-2:ni+3,-2:nj+3,-2:nk+3)
    real :: dvn(1:ni,1:nj,1:nk,3)
    integer :: i,j,k,sbcl,sbcr
    real :: vni(-2:ni+3),vnj(-2:nj+3),vnk(-2:nk+3)
    real :: vci(0:ni),vcj(0:nj),vck(0:nk)
    real :: dvni(1:ni),dvnj(1:nj),dvnk(1:nk)

    if ( ni > nijk2nd ) then
    
    do k=1,nk
    do j=1,nj
       do i=-2,ni+3
          vni(i) = vn(i,j,k)
       end do
       
       if ( ro(-2,j,k) < small ) then
          sbcl = 1
       else
          sbcl = 0
       end if
       
       if ( ro(ni+3,j,k) < small ) then
          sbcr = 1
       else
          sbcr = 0
       end if
       
       call vcenter_via_node(ni,vni,sbcl,sbcr,vci)
       call dvnode_via_center(ni,vci,dvni)
       do i=1,ni
          dvn(i,j,k,1) = dvni(i)
       end do
    end do
    end do
    
    else
       do k=1,nk
       do j=1,nj
          do i=1,ni
             dvn(i,j,k,1) = 0.0
          end do
       end do
       end do
    end if    
    
    
    if ( nj > nijk2nd ) then
    
    do k=1,nk
    do i=1,ni
       do j=-2,nj+3
          vnj(j) = vn(i,j,k)
       end do
       
       if ( ro(i,-2,k) < small ) then
          sbcl = 1
       else
          sbcl = 0
       end if
       
       if ( ro(i,nj+3,k) < small ) then
          sbcr = 1
       else
          sbcr = 0
       end if  
            
       call vcenter_via_node(nj,vnj,sbcl,sbcr,vcj)
       call dvnode_via_center(nj,vcj,dvnj)
       do j=1,nj
          dvn(i,j,k,2) = dvnj(j)
       end do
    end do
    end do
    
    else
       do k=1,nk
       do i=1,ni
          do j=1,nj
             dvn(i,j,k,2) = 0.0
          end do
       end do
       end do
    end if    
   
    
    if ( nk > nijk2nd ) then
    
    do j=1,nj
    do i=1,ni
       do k=-2,nk+3
          vnk(k) = vn(i,j,k)
       end do
       
       if ( ro(i,j,-2) < small ) then
          sbcl = 1
       else
          sbcl = 0
       end if
       
       if ( ro(i,j,nk+3) < small ) then
          sbcr = 1
       else
          sbcr = 0
       end if  
              
       call vcenter_via_node(nk,vnk,sbcl,sbcr,vck)
       call dvnode_via_center(nk,vck,dvnk)
       do k=1,nk
          dvn(i,j,k,3) = dvnk(k)
       end do
    end do
    end do

    else
       do j=1,nj
       do i=1,ni
          do k=1,nk
             dvn(i,j,k,3) = 0.0
          end do
       end do
       end do
    end if    
        
end subroutine get_dvnode

subroutine get_dtnode(ni,nj,nk,ro,ps,vn,dvn)
    use global_variables,only : moo,gama,nijk2nd,small
    implicit none
    integer :: ni,nj,nk
    real :: ro(-2:ni+3,-2:nj+3,-2:nk+3)
    real :: ps(-2:ni+3,-2:nj+3,-2:nk+3)
    real :: vn(-2:ni+3,-2:nj+3,-2:nk+3)
    real :: dvn(1:ni,1:nj,1:nk,3)
    integer :: i,j,k,sbcl,sbcr
    real :: vni(-2:ni+3),vnj(-2:nj+3),vnk(-2:nk+3)
    real :: vci(0:ni),vcj(0:nj),vck(0:nk)
    real :: dvni(1:ni),dvnj(1:nj),dvnk(1:nk)
    real :: gmmoo2
    
    gmmoo2 = gama*moo*moo
    
    
    if ( ni > nijk2nd ) then
    
    do k=1,nk
    do j=1,nj
       do i=-2,ni+3
          vni(i) = vn(i,j,k)
       end do
       
       if ( ro(-2,j,k) < small ) then
          sbcl = 1
       else
          sbcl = 0
          
          do i=-2,0
             vni(i) = gmmoo2*ps(i,j,k)/ro(i,j,k)
          end do
       end if
       
       if ( ro(ni+3,j,k) < small ) then
          sbcr = 1
       else
          sbcr = 0
          
          do i=ni+1,ni+3
             vni(i) = gmmoo2*ps(i,j,k)/ro(i,j,k)
          end do
       end if
       
       call vcenter_via_node(ni,vni,sbcl,sbcr,vci)
       call dvnode_via_center(ni,vci,dvni)
       do i=1,ni
          dvn(i,j,k,1) = dvni(i)
       end do
    end do
    end do
    
    else
       do k=1,nk
       do j=1,nj
          do i=1,ni
             dvn(i,j,k,1) = 0.0
          end do
       end do
       end do
    end if    
    
    
    if ( nj > nijk2nd ) then
    
    do k=1,nk
    do i=1,ni
       do j=-2,nj+3
          vnj(j) = vn(i,j,k)
       end do
       
       if ( ro(i,-2,k) < small ) then
          sbcl = 1
       else
          sbcl = 0
          
          do j=-2,0
             vnj(j) = gmmoo2*ps(i,j,k)/ro(i,j,k)
          end do
       end if
       
       if ( ro(i,nj+3,k) < small ) then
          sbcr = 1
       else
          sbcr = 0
          
          do j=nj+1,nj+3
             vnj(j) = gmmoo2*ps(i,j,k)/ro(i,j,k)
          end do
       end if  
            
       call vcenter_via_node(nj,vnj,sbcl,sbcr,vcj)
       call dvnode_via_center(nj,vcj,dvnj)
       do j=1,nj
          dvn(i,j,k,2) = dvnj(j)
       end do
    end do
    end do
    
    else
       do k=1,nk
       do i=1,ni
          do j=1,nj
             dvn(i,j,k,2) = 0.0
          end do
       end do
       end do
    end if    
    
    
    if ( nk > nijk2nd ) then
    
    do j=1,nj
    do i=1,ni
       do k=-2,nk+3
          vnk(k) = vn(i,j,k)
       end do
       
       if ( ro(i,j,-2) < small ) then
          sbcl = 1
       else
          sbcl = 0
          
          do k=-2,0
             vnk(k) = gmmoo2*ps(i,j,k)/ro(i,j,k)
          end do
       end if
       
       if ( ro(i,j,nk+3) < small ) then
          sbcr = 1
       else
          sbcr = 0

          do k=nk+1,nk+3
             vnk(k) = gmmoo2*ps(i,j,k)/ro(i,j,k)
          end do
       end if  
              
       call vcenter_via_node(nk,vnk,sbcl,sbcr,vck)
       call dvnode_via_center(nk,vck,dvnk)
       do k=1,nk
          dvn(i,j,k,3) = dvnk(k)
       end do
    end do
    end do
    
    else
       do j=1,nj
       do i=1,ni
          do k=1,nk
             dvn(i,j,k,3) = 0.0
          end do
       end do
       end do
    end if 
            
end subroutine get_dtnode

subroutine dfnode_via_node4(nj,fn,dfn)     !! fourth order approximation to the first derivative of the flux
    implicit none
    integer :: nj
    real :: fn(1:nj),dfn(1:nj)
    integer :: j
    
    do j=3,nj-2
       dfn(j) = ( 8.0*( fn(j+1) - fn(j-1) ) - ( fn(j+2) - fn(j-2) ) ) / 12.0                                ! fourth-order
    end do
    
    dfn(1) = ( -11.0*fn(1) + 18.0*fn(2) -  9.0*fn(3) +  2.0*fn(4) ) / 6.0                                   ! third-order
    dfn(2) = (  -3.0*fn(1) - 10.0*fn(2) + 18.0*fn(3) -  6.0*fn(4) + fn(5) ) / 12.0                          ! fourth-order

    dfn(nj-1) = -(  -3.0*fn(nj) - 10.0*fn(nj-1) + 18.0*fn(nj-2) -  6.0*fn(nj-3) + fn(nj-4) ) / 12.0         ! fourth-order
    dfn(nj)   = -( -11.0*fn(nj) + 18.0*fn(nj-1) -  9.0*fn(nj-2) +  2.0*fn(nj-3) ) / 6.0                     ! third-order
    
end subroutine dfnode_via_node4

subroutine dfnode_via_node(nj,fn,dfn)     !! sixth order approximation to the first derivative of the flux
    use global_variables,only : nijk2nd
    implicit none
    integer :: nj
    real :: fn(1:nj),dfn(1:nj)
    integer :: j
    
    if ( nj > nijk2nd ) then
    do j=4,nj-3
       dfn(j) = ( 45.0*( fn(j+1) - fn(j-1) ) - 9.0*( fn(j+2) - fn(j-2) ) + ( fn(j+3) - fn(j-3) ) ) / 60.0                   ! sixth-order
    end do
    
    dfn(1) = ( -11.0*fn(1) + 18.0*fn(2) -  9.0*fn(3) +  2.0*fn(4) ) / 6.0                                                   ! third-order
    dfn(2) = (  -3.0*fn(1) - 10.0*fn(2) + 18.0*fn(3) -  6.0*fn(4) +      fn(5) ) / 12.0                                     ! fourth-order
    dfn(3) = (   3.0*fn(1) - 30.0*fn(2) - 20.0*fn(3) + 60.0*fn(4) - 15.0*fn(5) + 2.0*fn(6) ) / 60.0                         ! fifth-order
!!    dfn(1) = ( -25.0*fn(1) + 48.0*fn(2) - 36.0*fn(3) + 16.0*fn(4) -  3.0*fn(5) ) / 12.0                                     ! fourth-order
!!    dfn(3) = (       fn(1) -  8.0*fn(2)              +  8.0*fn(4) -      fn(5) ) / 12.0                                     ! fourth-order
    
    dfn(nj-2) = -(   3.0*fn(nj) - 30.0*fn(nj-1) - 20.0*fn(nj-2) + 60.0*fn(nj-3) - 15.0*fn(nj-4) + 2.0*fn(nj-5) ) / 60.0     ! fifth-order
    dfn(nj-1) = -(  -3.0*fn(nj) - 10.0*fn(nj-1) + 18.0*fn(nj-2) -  6.0*fn(nj-3) +      fn(nj-4) ) / 12.0                    ! fourth-order
    dfn(nj)   = -( -11.0*fn(nj) + 18.0*fn(nj-1) -  9.0*fn(nj-2) +  2.0*fn(nj-3) ) / 6.0                                     ! third-order
!!    dfn(nj-2) = -(       fn(nj) -  8.0*fn(nj-1)                 +  8.0*fn(nj-3) -      fn(nj-4) ) / 12.0                    ! fourth-order
!!    dfn(nj)   = -( -25.0*fn(nj) + 48.0*fn(nj-1) - 36.0*fn(nj-2) + 16.0*fn(nj-3) -  3.0*fn(nj-4) ) / 12.0                    ! fourth-order
    else
       do j=1,nj
          dfn(j) = 0.0
       end do
    end if
    
end subroutine dfnode_via_node

subroutine dvnode_via_center(nj,vc,dvn)   !! fourth order approximation to the first derivative of the variable
    implicit none
    integer :: nj
    real :: vc(0:nj),dvn(1:nj)
    integer :: j
    
    do j=2,nj-1
       dvn(j) = ( vc(j-2) - 27.0*vc(j-1) + 27.0*vc(j) - vc(j+1) ) / 24.0           ! fourth-order
    end do
    
    dvn(1)  =  ( -23.0*vc(0)  + 21.0*vc(1)    + 3.0*vc(2)    - vc(3)    ) / 24.0   ! third-order
    dvn(nj) = -( -23.0*vc(nj) + 21.0*vc(nj-1) + 3.0*vc(nj-2) - vc(nj-3) ) / 24.0   ! third-order
       
end subroutine dvnode_via_center

subroutine dvnode_via_center6(nj,vc,dvn)    !! sixth order approximation to the first derivative of the variable
    implicit none
    integer :: nj
    real :: vc(0:nj),dvn(1:nj)
    integer :: j

    do j=3,nj-2
       dvn(j) = ( 2250.0 * ( vc(j) - vc(j-1) )  - 125.0 * ( vc(j+1) - vc(j-2) ) + 9.0 * ( vc(j+2) - vc(j-3) ) ) / 1920.0          ! sixth-order
    end do
    
    dvn(2)    =  ( 71.0*vc(0)  - 2115.0*vc(1)    + 2070.0*vc(2)    + 10.0*vc(3)    - 45.0*vc(4)    + 9.0*vc(5)    ) / 1920.0      ! fifth-order
    dvn(nj-1) = -( 71.0*vc(nj) - 2115.0*vc(nj-1) + 2070.0*vc(nj-2) + 10.0*vc(nj-3) - 45.0*vc(nj-4) + 9.0*vc(nj-5) ) / 1920.0      ! fifth-order
    
    dvn(1)  =  ( -22.0*vc(0)  + 17.0*vc(1)    + 9.0*vc(2)    - 5.0*vc(3)    + vc(4)    ) / 24.0                                   ! fourth-order
    dvn(nj) = -( -22.0*vc(nj) + 17.0*vc(nj-1) + 9.0*vc(nj-2) - 5.0*vc(nj-3) + vc(nj-4) ) / 24.0                                   ! fourth-order
       
end subroutine dvnode_via_center6

subroutine vcenter_via_node(nj,vn,sbcl,sbcr,vc)  !! sixth order approximation to the centered variable
    implicit none
    integer :: nj,sbcl,sbcr
    real :: vn(-2:nj+3),vc(0:nj)
    integer :: j,k
    real :: f(3,0:nj),s(3,0:nj),t(3,0:nj)
    real :: ck(3),eps,is,bk(3),bksum,wk(3),rvk(3)
    
    ck(1) =  3.0D0/16.0D0
    ck(2) = 10.0D0/16.0D0
    ck(3) = ck(1)
    
    eps = 1.0D-6
    
    do j=0,nj
       f(1,j) = (       vn(j-2) -  3.0*vn(j-1) - 21.0*vn(j)   + 23.0*vn(j+1) ) / 24.0        ! third-order
       f(2,j) = (       vn(j-1) - 27.0*vn(j)   + 27.0*vn(j+1) -      vn(j+2) ) / 24.0
       f(3,j) = ( -23.0*vn(j)   + 21.0*vn(j+1) +  3.0*vn(j+2) -      vn(j+3) ) / 24.0
       
!!       s(1,j) = (    -vn(j-2) + 5.0*vn(j-1) - 7.0*vn(j)   + 3.0*vn(j+1) ) / 2.0              ! second-order
!!       s(2,j) = (     vn(j-1) -     vn(j)   -     vn(j+1) +     vn(j+2) ) / 2.0
!!       s(3,j) = ( 3.0*vn(j)   - 7.0*vn(j+1) + 5.0*vn(j+2) -     vn(j+3) ) / 2.0
       
       t(1,j) = -vn(j-2) + 3.0*vn(j-1) - 3.0*vn(j)   + vn(j+1)                               ! first-order
       !!t(2,j) = -vn(j-1) + 3.0*vn(j)   - 3.0*vn(j+1) + vn(j+2)
       !!t(3,j) = -vn(j)   + 3.0*vn(j+1) - 3.0*vn(j+2) + vn(j+3)
    end do
    
    do j=0,nj-2
       t(2,j) = t(1,j+1)
       t(3,j) = t(1,j+2)
    end do
    
    t(2,nj-1) = -vn(nj-2) + 3.0*vn(nj-1) - 3.0*vn(nj) + vn(nj+1)
    t(2,nj)   = -vn(nj-1) + 3.0*vn(nj) - 3.0*vn(nj+1) + vn(nj+2)
    
    t(3,nj-1) = t(2,nj)
    t(3,nj  ) = -vn(nj) + 3.0*vn(nj+1) - 3.0*vn(nj+2) + vn(nj+3)
    
    do j=0,nj
       do k=1,3
          is = f(k,j)*f(k,j) + t(k,j)*t(k,j)
          is = is + eps
          bk(k) = ck(k)/(is*is)
       end do
    
       bksum = bk(1) + bk(2) + bk(3)
       do k=1,3
          !!wk(k) = ck(k)
          wk(k) = bk(k)/bksum
       end do
       rvk(1) = (     vn(j-2) -  5.0*vn(j-1) + 15.0*vn(j)   + 5.0*vn(j+1) ) / 16.0
       rvk(2) = (    -vn(j-1) +  9.0*vn(j)   +  9.0*vn(j+1) -     vn(j+2) ) / 16.0
       rvk(3) = ( 5.0*vn(j)   + 15.0*vn(j+1) -  5.0*vn(j+2) +     vn(j+3) ) / 16.0
       
       vc(j) = wk(1)*rvk(1) + wk(2)*rvk(2) + wk(3)*rvk(3)                                       ! sixth-order
    end do
    
    if (sbcl == 1) then
       vc(2) = (  -5.0*vn(1) +  60.0*vn(2) +  90.0*vn(3) -  20.0*vn(4) +  3.0*vn(5) ) / 128.0   ! fifth-order
       vc(1) = (  35.0*vn(1) + 140.0*vn(2) -  70.0*vn(3) +  28.0*vn(4) -  5.0*vn(5) ) / 128.0   ! fifth-order
       vc(0) = ( 315.0*vn(1) - 420.0*vn(2) + 378.0*vn(3) - 180.0*vn(4) + 35.0*vn(5) ) / 128.0   ! fifth-order
!!       vc(1) = (   5.0*vn(1) +  15.0*vn(2) -   5.0*vn(3) +       vn(4) ) / 16.0                 ! fourth-order
!!       vc(0) = (  35.0*vn(1) -  35.0*vn(2) +  21.0*vn(3) -   5.0*vn(4) ) / 16.0                 ! fourth-order
       
!!       vc(2) = (     -vn(1) +  9.0*vn(2) +  9.0*vn(3) -      vn(4) ) / 16.0                     ! fourth-order
!!       vc(0) = ( 15.0*vn(1) - 10.0*vn(2) +  3.0*vn(3) ) / 8.0                                   ! third-order
    end if
    
    if (sbcr == 1) then
       vc(nj-2) = (  -5.0*vn(nj) +  60.0*vn(nj-1) +  90.0*vn(nj-2) -  20.0*vn(nj-3) +  3.0*vn(nj-4) ) / 128.0   ! fifth-order
       vc(nj-1) = (  35.0*vn(nj) + 140.0*vn(nj-1) -  70.0*vn(nj-2) +  28.0*vn(nj-3) -  5.0*vn(nj-4) ) / 128.0   ! fifth-order
       vc(nj)   = ( 315.0*vn(nj) - 420.0*vn(nj-1) + 378.0*vn(nj-2) - 180.0*vn(nj-3) + 35.0*vn(nj-4) ) / 128.0   ! fifth-order
!!       vc(nj-1) = (   5.0*vn(nj) +  15.0*vn(nj-1) -   5.0*vn(nj-2) +       vn(nj-3) ) / 16.0                    ! fourth-order
!!       vc(nj)   = (  35.0*vn(nj) -  35.0*vn(nj-1) +  21.0*vn(nj-2) -   5.0*vn(nj-3) ) / 16.0                    ! fourth-order

!!       vc(nj-2) = (      -vn(nj) +   9.0*vn(nj-1) +   9.0*vn(nj-2) -       vn(nj-3) ) / 16.0                    ! fourth-order
!!       vc(nj)   = (  15.0*vn(nj) -  10.0*vn(nj-1) +   3.0*vn(nj-2) ) / 8.0                                      ! third-order
    end if
    
end subroutine vcenter_via_node


subroutine vcenter_via_node0(nj,vn,sbcl,sbcr,vc)
    implicit none
    integer :: nj,sbcl,sbcr
    real :: vn(-2:nj+3),vc(0:nj)
    integer :: j,k
    real :: f(3,0:nj),s(3,0:nj),t(3,0:nj)
    real :: ck(3),eps,is,bk(3),bksum,wk(3),vtc,rvk(3)
    
    ck(1) =  3.0D0/16.0D0
    ck(2) = 10.0D0/16.0D0
    ck(3) = ck(1)
    
    eps = 1.0D-6
    
    do j=0,nj
       f(1,j) = (       vn(j-2) -  6.0*vn(j-1) + 3.0*vn(j  ) + 2.0*vn(j+1) ) / 6.0
       f(2,j) = ( - 2.0*vn(j-1) -  3.0*vn(j  ) + 6.0*vn(j+1) -     vn(j+2) ) / 6.0
       f(3,j) = ( -11.0*vn(j  ) + 18.0*vn(j+1) - 9.0*vn(j+2) + 2.0*vn(j+3) ) / 6.0
       
       s(1,j) =     vn(j-1) - 2.0*vn(j  ) +     vn(j+1)
       s(2,j) =      s(1,j)
       s(3,j) = 2.0*vn(j  ) - 5.0*vn(j+1) + 4.0*vn(j+2) - vn(j+3)
       
       t(1,j) =    -vn(j-2) + 3.0*vn(j-1) - 3.0*vn(j  ) + vn(j+1)
    end do
    
    do j=0,nj-2
       t(2,j) = t(1,j+1)
       t(3,j) = t(1,j+2)
    end do
    
    do j=nj-1,nj
       t(2,j) = -vn(j-1) + 3.0*vn(j) - 3.0*vn(j+1) + vn(j+2)
    end do
    
    t(3,nj-1) = t(2,nj)
    t(3,nj  ) = -vn(nj) + 3.0*vn(nj+1) - 3.0*vn(nj+2) + vn(nj+3)
    
    do j=0,nj
       do k=1,3
          is = f(k,j)*f(k,j) + t(k,j)*t(k,j)
          is = is + eps
          bk(k) = ck(k)/(is*is)
       enddo
    
       bksum = bk(1) + bk(2) + bk(3)
       do k=1,3
          wk(k) = bk(k)/bksum
          rvk(k) = 0.5*f(k,j) + 0.125*s(k,j) + t(k,j)/48.0
       enddo
       
       vtc = vn(j)
       vc(j) = vtc + wk(1)*rvk(1) + wk(2)*rvk(2) + wk(3)*rvk(3)
    end do  
    
    if (sbcl == 1) then
       vc(2) = (  -5.0*vn(1) +  60.0*vn(2) +  90.0*vn(3) -  20.0*vn(4) +  3.0*vn(5) ) / 128.0   ! fifth-order
       vc(1) = (  35.0*vn(1) + 140.0*vn(2) -  70.0*vn(3) +  28.0*vn(4) -  5.0*vn(5) ) / 128.0   ! fifth-order
       vc(0) = ( 315.0*vn(1) - 420.0*vn(2) + 378.0*vn(3) - 180.0*vn(4) + 35.0*vn(5) ) / 128.0   ! fifth-order
       
!!       vc(1) = (  5.0*vn(1) + 15.0*vn(2) -  5.0*vn(3) +      vn(4) ) / 16.0                     ! fourth-order
!!       vc(0) = ( 35.0*vn(1) - 35.0*vn(2) + 21.0*vn(3) -  5.0*vn(4) ) / 16.0                     ! fourth-order
!!       vc(2) = (     -vn(1) +  9.0*vn(2) +  9.0*vn(3) -      vn(4) ) / 16.0                     ! fourth-order
!!       vc(0) = ( 15.0*vn(1) - 10.0*vn(2) +  3.0*vn(3) ) / 8.0                                   ! third-order
    end if
    
    if (sbcr == 1) then
       vc(nj-2) = (  -5.0*vn(nj) +  60.0*vn(nj-1) +  90.0*vn(nj-2) -  20.0*vn(nj-3) +  3.0*vn(nj-4) ) / 128.0   ! fifth-order
       vc(nj-1) = (  35.0*vn(nj) + 140.0*vn(nj-1) -  70.0*vn(nj-2) +  28.0*vn(nj-3) -  5.0*vn(nj-4) ) / 128.0   ! fifth-order
       vc(nj)   = ( 315.0*vn(nj) - 420.0*vn(nj-1) + 378.0*vn(nj-2) - 180.0*vn(nj-3) + 35.0*vn(nj-4) ) / 128.0   ! fifth-order
       
!!       vc(nj-1) = (  5.0*vn(nj) + 15.0*vn(nj-1) -  5.0*vn(nj-2) +      vn(nj-3) ) / 16.0                        ! fourth-order
!!       vc(nj)   = ( 35.0*vn(nj) - 35.0*vn(nj-1) + 21.0*vn(nj-2) -  5.0*vn(nj-3) ) / 16.0                        ! fourth-order
!!       vc(nj-2) = (     -vn(nj) +  9.0*vn(nj-1) +  9.0*vn(nj-2) -      vn(nj-3) ) / 16.0                        ! fourth-order
!!       vc(nj)   = ( 15.0*vn(nj) - 10.0*vn(nj-1) +  3.0*vn(nj-2) ) / 8.0                                         ! third-order
    end if
    
end subroutine vcenter_via_node0


