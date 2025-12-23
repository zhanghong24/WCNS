!_____________________________________________________________________!
subroutine boundary(nb)
    use global_variables, only : mb_bc,connect_order,cic1,cic2,cic3,cic4,cic5
#ifdef PARALLEL
    use mod_parallels,only : boundary_n1_parallel,boundary_n1_vir1_parallel,boundary_n1_vir2_parallel
#endif
    implicit none
    integer :: nb,nr,nrmax,bctype,m,k

    nrmax = mb_bc(nb)%nregions             !本块共有nrmax个边界需要处理
    !具体处理边界条件
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype    !* 对接类型 *!
       if( bctype < 0 ) then               !* 对接边界
          if (connect_order == 1) then
#ifdef PARALLEL
             call boundary_n1_vir1_parallel(nb,nr,bctype)
#else
             call boundary_n1_vir1(nb,nr,bctype)   !* 1阶精度导数法求虚点值
#endif
          else if (connect_order == 2) then
#ifdef PARALLEL
             call boundary_n1_vir2_parallel(nb,nr,bctype)
#else
             call boundary_n1_vir2(nb,nr,bctype)   !* 2阶精度导数法求虚点值
#endif
          else if (connect_order == 0) then
#ifdef PARALLEL
             call boundary_n1_parallel(nb,nr,bctype)
#else
             call boundary_n1(nb,nr,bctype)        !* 直接取对接面另一侧对应点的值
                                                   !*TGH. 一般对接和特征对接都适用，且都赋了虚点的值
                                                   !*TGH. 一般对接还要计算边界上的值，湍流模型也计算了边界值
                                                   !*TGH. 特征对接不计算边界上的值，但湍流模型计算了边界值
#endif
          else
             write(*,*)'确定对接虚点方案出错，建议connect_order=0，1，或2'
             write(*,*)'Please input connect_order ='
             read(*,*)connect_order
          endif

       elseif( bctype == 2 ) then          !完全气体/完全非催化 绝热壁面
          call boundary_2(nb,nr,bctype)
       elseif( bctype == 3 ) then          !对称条件
          call boundary_3(nb,nr,bctype)
       elseif( bctype == 4 ) then          !远场边界
          call boundary_4(nb,nr,bctype)
       elseif( bctype == 5 ) then          !超声速入流边界
          call boundary_5(nb,nr,bctype)
       elseif( bctype == 6 ) then          !超声速出口边界
          call boundary_6(nb,nr,bctype)
       elseif( bctype == 7 ) then          !奇性面边界
          call boundary_7(nb,nr,bctype)
       elseif( bctype == 8 ) then          !周期边界
          call boundary_8(nb,nr,bctype)
       elseif( bctype == 20 ) then         !完全气体 等温壁面
          call boundary_20(nb,nr,bctype)
       elseif( bctype == 21 ) then         !完全气体 恒温壁面
          call boundary_21(nb,nr,bctype)
       elseif( bctype == 71 .or. bctype == 72 .or. bctype == 73) then
          call boundary_7123_old(nb,nr,bctype)  !奇性轴
       else
          write(*,*)' 程序不能处理边界条件 bctype = ',bctype
          stop
       endif
    enddo

    return
end subroutine boundary
!_____________________________________________________________________!
subroutine boundary_a(nb)
    use global_variables, only : mb_bc,cic1,cic2,cic3,cic4,cic5
    implicit none
    integer :: nb,nr,nrmax,bctype,m,k
! 处理物面，远场，奇性轴和超声速出入口

    nrmax = mb_bc(nb)%nregions             !本块共有nrmax个边界需要处理
    !具体处理边界条件
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype    !* 对接类型 *!
       if( bctype < 0 ) then               !* 对接边界

!         call boundary_n1_vir3(nb,nr,bctype)  !* 直接取对接面另一侧对应点的值（取3个点）
                                              !*TGH. 一般对接和特征对接都适用，且都赋了虚点的值
                                              !*TGH. 一般对接还要计算边界上的值，湍流模型也计算了边界值

       elseif( bctype == 2 ) then          !完全气体/完全非催化 绝热壁面
          call boundary_2(nb,nr,bctype)
       elseif( bctype == 3 ) then          !对称条件
 !         call boundary_3(nb,nr,bctype)
       elseif( bctype == 4 ) then          !远场边界
          call boundary_4(nb,nr,bctype)
       elseif( bctype == 5 ) then          !超声速入流边界
          call boundary_5(nb,nr,bctype)
       elseif( bctype == 6 ) then          !超声速出口边界
          call boundary_6(nb,nr,bctype)
       elseif( bctype == 7 ) then          !奇性面边界
          call boundary_7(nb,nr,bctype)
       elseif( bctype == 8 ) then          !周期边界
!          call boundary_8_TGH(nb,nr,bctype)   !该程序从网格上实现多点重叠周期
!          call boundary_8(nb,nr,bctype)       !该程序在网格上只有一点周期,且在本块上，但是在程序上再加上3点周期
!       elseif( bctype == 20 ) then         !完全气体 等温壁面
!          call boundary_20(nb,nr,bctype)
!       elseif( bctype == 21 ) then         !完全气体 恒温壁面
!          call boundary_21(nb,nr,bctype)
       elseif( bctype == 71 .or. bctype == 72 .or. bctype == 73) then
          call boundary_7123_old(nb,nr,bctype)  !奇性轴
       else
          write(*,*)' 程序不能处理边界条件 bctype = ',bctype
          stop
       endif
    enddo
    return
end subroutine boundary_a
!_____________________________________________________________________!
subroutine boundary_B(nb)
    use global_variables, only : mb_bc,cic1,cic2,cic3,cic4,cic5
#ifdef PARALLEL
    use mod_parallels,only : boundary_n1_vir3_parallel
#endif
    implicit none
    integer :: nb,nr,nrmax,bctype,m,k
!	处理周期边界和对称边界，以及对接边界

    nrmax = mb_bc(nb)%nregions             !本块共有nrmax个边界需要处理
    !具体处理边界条件
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype    !* 对接类型 *!
       if( bctype < 0 ) then               !* 对接边界

#ifdef PARALLEL
         call boundary_n1_vir3_parallel(nb,nr,bctype)
#else
         call boundary_n1_vir3(nb,nr,bctype)  !* 直接取对接面另一侧对应点的值（取3个点）
                                               !*TGH. 一般对接和特征对接都适用，且都赋了虚点的值
                                               !*TGH. 一般对接还要计算边界上的值，湍流模型也计算了边界值
#endif

       elseif( bctype == 2 ) then          !完全气体/完全非催化 绝热壁面
!          call boundary_2(nb,nr,bctype)
       elseif( bctype == 3 ) then          !对称条件
          call boundary_3(nb,nr,bctype)
       elseif( bctype == 4 ) then          !远场边界
!          call boundary_4(nb,nr,bctype)
       elseif( bctype == 5 ) then          !超声速入流边界
!          call boundary_5(nb,nr,bctype)
       elseif( bctype == 6 ) then          !超声速出口边界
!          call boundary_6(nb,nr,bctype)
       elseif( bctype == 7 ) then          !奇性面边界
!          call boundary_7(nb,nr,bctype)
       elseif( bctype == 8 ) then          !周期边界
!          call boundary_8_TGH(nb,nr,bctype)  !该程序从网格上实现多点重叠周期
          call boundary_8(nb,nr,bctype)       !该程序在网格上只有一点周期,且在本块上，但是在程序上再加上3点周期
!       elseif( bctype == 20 ) then         !完全气体 等温壁面
!          call boundary_20(nb,nr,bctype)
!       elseif( bctype == 21 ) then         !完全气体 恒温壁面
!          call boundary_21(nb,nr,bctype)
       elseif( bctype == 71 .or. bctype == 72 .or. bctype == 73) then
!          call boundary_7123_old(nb,nr,bctype)  !奇性轴
       else
          write(*,*)' 程序不能处理边界条件 bctype = ',bctype
          stop
       endif
    enddo

	return
end subroutine boundary_B
!_____________________________________________________________________!
subroutine boundary_sequence(nb)
    use global_variables, only : mb_bc,cic1,cic2,cic3,cic4,cic5
#ifdef PARALLEL
    use mod_parallels,only : boundary_n1_vir3_parallel
#endif
    implicit none
    integer :: nb,nr,nrmax,bctype,ib
!	处理周期边界和对称边界，以及对接边界

    nrmax = mb_bc(nb)%nregions             !本块共有nrmax个边界需要处理
    !具体处理边界条件
    do ib=1,nrmax
       nr = mb_bc(nb)%bcindexs(ib)
       bctype = mb_bc(nb)%bc(nr)%bctype    !* 对接类型 *!
       if( bctype < 0 ) then               !* 对接边界

#ifdef PARALLEL
         call boundary_n1_vir3_parallel(nb,nr,bctype)
#else
         call boundary_n1_vir3(nb,nr,bctype)  !* 直接取对接面另一侧对应点的值（取3个点）
                                              !*TGH. 一般对接和特征对接都适用，且都赋了虚点的值
                                              !*TGH. 一般对接还要计算边界上的值，湍流模型也计算了边界值
#endif

       elseif( bctype == 2 ) then          !完全气体/完全非催化 绝热壁面
          call boundary_2(nb,nr,bctype)
       elseif( bctype == 3 ) then          !对称条件
          call boundary_3(nb,nr,bctype)
       elseif( bctype == 4 ) then          !远场边界
          call boundary_4(nb,nr,bctype)
       elseif( bctype == 5 ) then          !超声速入流边界
          call boundary_5(nb,nr,bctype)
       elseif( bctype == 6 ) then          !超声速出口边界
          call boundary_6(nb,nr,bctype)
       else
          write(*,*)' 程序不能处理边界条件 bctype = ',bctype
          stop
       endif
    enddo

	return
end subroutine boundary_sequence
!_____________________________________________________________________!
subroutine boundary_2(nb,nr,bctype)   !完全气体 绝热边界
    use global_variables,only : nvis
    implicit none
    integer :: nb,nr,bctype
    
    if ( nvis == 1 ) then
       call viscid_wall(nb,nr,bctype)
    else
       call inviscid_wall(nb,nr,bctype)
    endif


    return
end subroutine boundary_2
!_____________________________________________________________________!
subroutine boundary_3(nb,nr,bctype)      !对称边界条件
    use global_variables
    implicit none
    integer :: nb,nr,bctype,m
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,i,j,k,n,nt
    integer :: it,jt,kt,it0,jt0,kt0
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: prim_s1(nl),q_1(nl),zf1
    real :: nx,ny,nz,sss1,nx1,ny1,nz1
	real :: rm,um,vm,wm,pm
	integer :: nbi,nbj,nbk

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
    enddo
    s_nd  = mb_bc(nb)%bc(nr)%s_nd           !边界面方向:1,2,3对应于i,j,k

	nbi = mb_dim(nb,1)
	nbj = mb_dim(nb,2)
	nbk = mb_dim(nb,3)

	if(mb_dim(nb,s_nd) <= nijk2d)then  !二维情况,且展向必须为Z方向
        do k = s_st(3),s_ed(3)
        do j = s_st(2),s_ed(2)
        do i = s_st(1),s_ed(1)
		    do n=1,1 !2-method
				is = i + s_lr3d(1)
				js = j + s_lr3d(2)
				ks = k + s_lr3d(3)
				it = i - s_lr3d(1)
				jt = j - s_lr3d(2)
				kt = k - s_lr3d(3)

				it = min(nbi,max(1,it))
				jt = min(nbj,max(1,jt))
				kt = min(nbk,max(1,kt))


					call BC_symmetry_turbulence(nb,is,js,ks,it,jt,kt,i,j,k)

				rm =  mb_r(nb)%a3d(it,jt,kt)
				um =  mb_u(nb)%a3d(it,jt,kt)
				vm =  mb_v(nb)%a3d(it,jt,kt)
				wm = -mb_w(nb)%a3d(it,jt,kt)
				pm =  mb_p(nb)%a3d(it,jt,kt)

                mb_r(nb)%a3d(is,js,ks) =  rm
                mb_u(nb)%a3d(is,js,ks) =  um
                mb_v(nb)%a3d(is,js,ks) =  vm
                mb_w(nb)%a3d(is,js,ks) =  wm
                mb_p(nb)%a3d(is,js,ks) =  pm
				
!				mb_vist(nb)%a3d(i,j,k) = mb_vist(nb)%a3d(it,jt,kt)


                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
!!                it = it - s_lr3d(1)
!!                jt = jt - s_lr3d(2)
!!                kt = kt - s_lr3d(3)
!!				it = min(nbi,max(1,it))
!!				jt = min(nbj,max(1,jt))
!!				kt = min(nbk,max(1,kt))

!!				rm =  mb_r(nb)%a3d(it,jt,kt)
!!				um =  mb_u(nb)%a3d(it,jt,kt)
!!				vm =  mb_v(nb)%a3d(it,jt,kt)
				wm = 0.0 !-mb_w(nb)%a3d(it,jt,kt)
!!				pm =  mb_p(nb)%a3d(it,jt,kt)
                mb_r(nb)%a3d(is,js,ks) =  rm
                mb_u(nb)%a3d(is,js,ks) =  um
                mb_v(nb)%a3d(is,js,ks) =  vm
                mb_w(nb)%a3d(is,js,ks) =  wm
                mb_p(nb)%a3d(is,js,ks) =  pm


                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
!!                it = it - s_lr3d(1)
!!                jt = jt - s_lr3d(2)
!!                kt = kt - s_lr3d(3)
!!				it = min(nbi,max(1,it))
!!				jt = min(nbj,max(1,jt))
!!				kt = min(nbk,max(1,kt))
!!
!!
!!				rm =  mb_r(nb)%a3d(it,jt,kt)
!!				um =  mb_u(nb)%a3d(it,jt,kt)
!!				vm =  mb_v(nb)%a3d(it,jt,kt)
				wm = 0.0 !-mb_w(nb)%a3d(it,jt,kt)
!!				pm =  mb_p(nb)%a3d(it,jt,kt)
                mb_r(nb)%a3d(is,js,ks) =  rm
                mb_u(nb)%a3d(is,js,ks) =  um
                mb_v(nb)%a3d(is,js,ks) =  vm
                mb_w(nb)%a3d(is,js,ks) =  wm
                mb_p(nb)%a3d(is,js,ks) =  pm

			enddo
		enddo
		enddo
		enddo
        if ( method == 1 ) then !有限差分 给定边界上的值
          call dif_average(nb,nr,bctype,s_st,s_ed)
        endif    
		return
	endif

    nx1 = 0.0
    ny1 = 0.0
    nz1 = 0.0

    it = s_lr3d(1) !*tgh. 边界面时=1或-1，非边界面时=0 *!
    jt = s_lr3d(2)
    kt = s_lr3d(3)

    !确定对称面法向
    if ( (it+jt+kt) == 1 .and. method == 0 ) then !* tgh. 右边界面+有限体积 *!
       is = it
       js = jt
       ks = kt
    else
       is = 0
       js = 0
       ks = 0
    endif

!   it0 = (1 -abs(it))*method !* tgh.这三个量在非边界面且为有限差分时取1，其他情况为零 *!
!   jt0 = (1 -abs(jt))*method
!   kt0 = (1 -abs(kt))*method

    it0 = min(ni/2,2*(1 -abs(it)))  !* tgh.这三个量在非边界面且为有限差分时取2，其他情况为零 *!
    jt0 = min(nj/2,2*(1 -abs(jt)))
    kt0 = min(nk/2,2*(1 -abs(kt)))

!*TGH. Modified by TU Guohua
    if(.false.)then !tgh. 原求对称面法向的方法：所有点平均
       do k = ks+s_st(3)+kt0, ks+s_ed(3)-kt0
       do j = js+s_st(2)+jt0, js+s_ed(2)-jt0  !* tgh. 即去掉了4条边
       do i = is+s_st(1)+it0, is+s_ed(1)-it0     !* tgh. 在边界面内的坐标活动范围为S_ST+1:S_ED-1
             call getnxyz_mml(nx,ny,nz,i,j,k,it,jt,kt)   !计算单位法向矢量 
!             call getnxyz_tgh(nx,ny,nz,i,j,k,it,jt,kt)   !计算单位法向矢量       
             nx1 = nx + nx1
             ny1 = ny + ny1
             nz1 = nz + nz1
       enddo
       enddo
       enddo
    else !tgh. 新的求对称面的法向方向只取中心点
       is=(s_st(1)+s_ed(1))/2
	   js=(s_st(2)+s_ed(2))/2
	   ks=(s_st(3)+s_ed(3))/2
       call getnxyz_mml(nx,ny,nz,is,js,ks,it,jt,kt)   !计算单位法向矢量          
!      call getnxyz_tgh(nx,ny,nz,is,js,ks,it,jt,kt)   !计算单位法向矢量          
      nx1 = nx
      ny1 = ny
      nz1 = nz
    endif
!*TGH. end Modified by TU Guohua

    sss1 = sqrt( nx1*nx1 + ny1*ny1 + nz1*nz1 )
    nx = nx1 / sss1   !*tgh. 这好像求得的时边界面上所有点法向矢量的平均值*!
    ny = ny1 / sss1
    nz = nz1 / sss1

    do k = s_st(3),s_ed(3)
    do j = s_st(2),s_ed(2)
    do i = s_st(1),s_ed(1)
        do n=1,2-method
            is = i + s_lr3d(1)*n
            js = j + s_lr3d(2)*n
			ks = k + s_lr3d(3)*n
			nt = n - 1 + method
			it = i - s_lr3d(1)*nt
			jt = j - s_lr3d(2)*nt
			kt = k - s_lr3d(3)*nt
			sss1 = 2.0*( nx*mb_u(nb)%a3d(it,jt,kt) + &
			             ny*mb_v(nb)%a3d(it,jt,kt) + nz*mb_w(nb)%a3d(it,jt,kt) )
			mb_r(nb)%a3d(is,js,ks) =  mb_r(nb)%a3d(it,jt,kt)
			mb_u(nb)%a3d(is,js,ks) =  mb_u(nb)%a3d(it,jt,kt) - sss1*nx
			mb_v(nb)%a3d(is,js,ks) =  mb_v(nb)%a3d(it,jt,kt) - sss1*ny
			mb_w(nb)%a3d(is,js,ks) =  mb_w(nb)%a3d(it,jt,kt) - sss1*nz
			mb_p(nb)%a3d(is,js,ks) =  mb_p(nb)%a3d(it,jt,kt)

				call BC_symmetry_turbulence(nb,is,js,ks,it,jt,kt,i,j,k)

		enddo

        !* 以下另外再赋两层虚点值，特别注意，第三层密度为正
        is = is + s_lr3d(1)
        js = js + s_lr3d(2)
        ks = ks + s_lr3d(3)
        it = it - s_lr3d(1)
        jt = jt - s_lr3d(2)
        kt = kt - s_lr3d(3)
        sss1 = 2.0*( nx*mb_u(nb)%a3d(it,jt,kt) + &
                     ny*mb_v(nb)%a3d(it,jt,kt) + nz*mb_w(nb)%a3d(it,jt,kt) )
        mb_r(nb)%a3d(is,js,ks) =  mb_r(nb)%a3d(it,jt,kt)
        mb_u(nb)%a3d(is,js,ks) =  mb_u(nb)%a3d(it,jt,kt) - sss1*nx
        mb_v(nb)%a3d(is,js,ks) =  mb_v(nb)%a3d(it,jt,kt) - sss1*ny
        mb_w(nb)%a3d(is,js,ks) =  mb_w(nb)%a3d(it,jt,kt) - sss1*nz
        mb_p(nb)%a3d(is,js,ks) =  mb_p(nb)%a3d(it,jt,kt)

        is = is + s_lr3d(1)
        js = js + s_lr3d(2)
        ks = ks + s_lr3d(3)
        it = it - s_lr3d(1)
        jt = jt - s_lr3d(2)
        kt = kt - s_lr3d(3)
        sss1 = 2.0*( nx*mb_u(nb)%a3d(it,jt,kt) + &
                     ny*mb_v(nb)%a3d(it,jt,kt) + nz*mb_w(nb)%a3d(it,jt,kt) )
        mb_r(nb)%a3d(is,js,ks) =  mb_r(nb)%a3d(it,jt,kt)          !第三层密度给正值
        mb_u(nb)%a3d(is,js,ks) =  mb_u(nb)%a3d(it,jt,kt) - sss1*nx
        mb_v(nb)%a3d(is,js,ks) =  mb_v(nb)%a3d(it,jt,kt) - sss1*ny
        mb_w(nb)%a3d(is,js,ks) =  mb_w(nb)%a3d(it,jt,kt) - sss1*nz
        mb_p(nb)%a3d(is,js,ks) =  mb_p(nb)%a3d(it,jt,kt)

    enddo
    enddo
    enddo

    if ( method == 1 ) then !有限差分 给定边界上的值
       call dif_average(nb,nr,bctype,s_st,s_ed)
    endif    
    return
end subroutine boundary_3
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine boundary_4(nb,nr,bctype)   !远场边界
    use global_variables, &
    only : nl,ns,roo,uoo,voo,woo,poo,coo,gama,method,mb_bc, &
           mb_r,mb_u,mb_v,mb_w,mb_p,mb_fs,sml_sss
    implicit none
    integer :: nb,nr,bctype,m
    integer :: nbs,s_nd,s_fix,s_lr,lr
    integer :: i,j,k,is,js,ks,it,jt,kt,is1,js1,ks1,n
    integer :: s_st(3),s_ed(3),s_lr3d(3),s0(3)
    real :: prim_s1(nl),q_1(nl)
    real :: zf1,zf2
    real :: nx,ny,nz,nxa,nya,nza,cgm1
    real :: rin,uin,vin,win,pin
    real :: rb,ub,vb,wb,pb,sb,ab,vnb
    real :: uref,vref,wref,vnref
    real :: cin,s_in,vnin,mach
    real :: s_oo,vnoo,rp,rn

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
       s0  (m) = 0
       if ( s_lr3d(m) == 1 .and. method == 0 ) s0  (m) = 1
    enddo

    s_nd  = mb_bc(nb)%bc(nr)%s_nd
    lr    = mb_bc(nb)%bc(nr)%s_lr

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             n = 1 - method
             is = i + s_lr3d(1)*n !*TGH. 有限差分时在边界面上，体积时在虚层上
             js = j + s_lr3d(2)*n
             ks = k + s_lr3d(3)*n
             it = i - s_lr3d(1) + s0(1) !*TGH. 在边界面的内一层上
             jt = j - s_lr3d(2) + s0(2)
             kt = k - s_lr3d(3) + s0(3)
             call getnxyz_mml ( nx,ny,nz,i+s0(1),j+s0(2),k+s0(3),s_lr3d(1),s_lr3d(2),s_lr3d(3) ) 
			 !计算单位法向矢量 
			           
             rin = mb_r(nb)%a3d(it,jt,kt)
             uin = mb_u(nb)%a3d(it,jt,kt)
             vin = mb_v(nb)%a3d(it,jt,kt)
             win = mb_w(nb)%a3d(it,jt,kt)
             pin = mb_p(nb)%a3d(it,jt,kt)

             s_in = pin/( rin ** gama )
             s_oo = poo/( roo ** gama )

             cin = sqrt(gama*pin/rin)

             cgm1 = lr/max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
             nxa = nx * cgm1
             nya = ny * cgm1
             nza = nz * cgm1
             !nxa,nya,nza为边界的单位外法向

             vnoo = nxa*uoo + nya*voo + nza*woo
             vnin = nxa*uin + nya*vin + nza*win

             !mach = abs(vnin)/cin
             mach = sqrt(uin*uin +vin*vin +win*win)/cin
             if ( mach < 1.0 ) then
                rp = vnin + 2.0 * cin / (gama-1.0)
                rn = vnoo - 2.0 * coo / (gama-1.0)
                vnb = 0.5  * ( rp + rn )
                ab  = 0.25 * ( gama - 1.0 ) * ( rp - rn )
                if ( vnb > 0.0 ) then
                   sb    = s_in
                   uref  = uin
                   vref  = vin
                   wref  = win
                   vnref = vnin
                else
                   sb    = s_oo
                   uref  = uoo
                   vref  = voo
                   wref  = woo
                   vnref = vnoo
                endif
                rb = ( ab*ab/gama/sb ) ** ( 1.0/(gama-1.0) )
                ub = uref + nxa * ( vnb - vnref )
                vb = vref + nya * ( vnb - vnref )
                wb = wref + nza * ( vnb - vnref )
                pb = sb * ( rb ** gama )
            else
                if ( vnin > 0.0 ) then
                   rb = rin
                   ub = uin
                   vb = vin
                   wb = win
                   pb = pin
                else
                   rb = roo
                   ub = uoo
                   vb = voo
                   wb = woo
                   pb = poo
                endif
             ENDIF

             mb_r(nb)%a3d(is,js,ks) =   rb
             mb_u(nb)%a3d(is,js,ks) =   ub
             mb_v(nb)%a3d(is,js,ks) =   vb
             mb_w(nb)%a3d(is,js,ks) =   wb
             mb_p(nb)%a3d(is,js,ks) =   pb

!						 call BC_farfield_turbulence(nb,is,js,ks)
						 call BC_farfield_turbulence_tgh(nb,is,js,ks,it,jt,kt,s_lr3d,s_nd)
	           !*tgh s_nd边界面方向:1,2,3对应于i,j,k


             n = 2 - method
             is1 = i + s_lr3d(1)*n
             js1 = j + s_lr3d(2)*n
             ks1 = k + s_lr3d(3)*n

             mb_r(nb)%a3d(is1,js1,ks1) =   rb
             mb_u(nb)%a3d(is1,js1,ks1) =   ub
             mb_v(nb)%a3d(is1,js1,ks1) =   vb
             mb_w(nb)%a3d(is1,js1,ks1) =   wb
             mb_p(nb)%a3d(is1,js1,ks1) =   pb

        !* 以下另外再赋两层虚点值，特别注意，第三层密度给的是负值
             is1 = is1 + s_lr3d(1)
             js1 = js1 + s_lr3d(2)
             ks1 = ks1 + s_lr3d(3)
             mb_r(nb)%a3d(is1,js1,ks1) =   rb
             mb_u(nb)%a3d(is1,js1,ks1) =   ub
             mb_v(nb)%a3d(is1,js1,ks1) =   vb
             mb_w(nb)%a3d(is1,js1,ks1) =   wb
             mb_p(nb)%a3d(is1,js1,ks1) =   pb

             is1 = is1 + s_lr3d(1)
             js1 = js1 + s_lr3d(2)
             ks1 = ks1 + s_lr3d(3)
             mb_r(nb)%a3d(is1,js1,ks1) =  -rb    !第三层密度给负值
             mb_u(nb)%a3d(is1,js1,ks1) =   ub
             mb_v(nb)%a3d(is1,js1,ks1) =   vb
             mb_w(nb)%a3d(is1,js1,ks1) =   wb
             mb_p(nb)%a3d(is1,js1,ks1) =   pb

          enddo
       enddo
    enddo
end subroutine boundary_4
!_____________________________________________________________________!
subroutine boundary_5(nb,nr,bctype)   !超声速来流边界
    use global_variables
    implicit none
    integer :: nb,nr,bctype,m
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: i,j,k,is,js,ks,n
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: prim_s1(nl),q_1(nl)
    real :: zf1,zf2

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
    enddo
    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             do n=1-method,2-method
                is = i + s_lr3d(1)*n
                js = j + s_lr3d(2)*n
                ks = k + s_lr3d(3)*n

                mb_r(nb)%a3d(is,js,ks) =   roo
                mb_u(nb)%a3d(is,js,ks) =   uoo
                mb_v(nb)%a3d(is,js,ks) =   voo
                mb_w(nb)%a3d(is,js,ks) =   woo
                mb_p(nb)%a3d(is,js,ks) =   poo

								call BC_farfield_turbulence(nb,is,js,ks)


             enddo

        !* 以下另外再赋两层虚点值，特别注意，第三层密度为正
             is = is + s_lr3d(1)
             js = js + s_lr3d(2)
             ks = ks + s_lr3d(3)
             mb_r(nb)%a3d(is,js,ks) =   roo
             mb_u(nb)%a3d(is,js,ks) =   uoo
             mb_v(nb)%a3d(is,js,ks) =   voo
             mb_w(nb)%a3d(is,js,ks) =   woo
             mb_p(nb)%a3d(is,js,ks) =   poo

             is = is + s_lr3d(1)
             js = js + s_lr3d(2)
             ks = ks + s_lr3d(3)
             mb_r(nb)%a3d(is,js,ks) =   roo
             mb_u(nb)%a3d(is,js,ks) =   uoo
             mb_v(nb)%a3d(is,js,ks) =   voo
             mb_w(nb)%a3d(is,js,ks) =   woo
             mb_p(nb)%a3d(is,js,ks) =   poo
          enddo
       enddo
    enddo
    return
end subroutine boundary_5
!_____________________________________________________________________!
subroutine boundary_6(nb,nr,bctype) !超声速出口边界
    use global_variables
    implicit none
    integer :: nb,nr,bctype,m,n
    integer :: is,js,ks,it,jt,kt,i,j,k
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: prim_s1(nl), q_1(nl), zf1
	real :: rm,um,vm,wm,pm
    integer :: it1,jt1,kt1,it2,jt2,kt2

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
    enddo

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             it = i - s_lr3d(1)*method
             jt = j - s_lr3d(2)*method
             kt = k - s_lr3d(3)*method
             do n=1-method,2-method
                is = i + s_lr3d(1)*n
                js = j + s_lr3d(2)*n
                ks = k + s_lr3d(3)*n
!___2nd
                it1 = is - s_lr3d(1)
                jt1 = js - s_lr3d(2)
                kt1 = ks - s_lr3d(3)
                it2 = is - s_lr3d(1)*2
                jt2 = js - s_lr3d(2)*2
                kt2 = ks - s_lr3d(3)*2
                
!!                rm = 2.d0*mb_r(nb)%a3d(it1,jt1,kt1) - mb_r(nb)%a3d(it2,jt2,kt2)
!!                um = 2.d0*mb_u(nb)%a3d(it1,jt1,kt1) - mb_u(nb)%a3d(it2,jt2,kt2)
!!                vm = 2.d0*mb_v(nb)%a3d(it1,jt1,kt1) - mb_v(nb)%a3d(it2,jt2,kt2)
!!                wm = 2.d0*mb_w(nb)%a3d(it1,jt1,kt1) - mb_w(nb)%a3d(it2,jt2,kt2)
!!                pm = 2.d0*mb_p(nb)%a3d(it1,jt1,kt1) - mb_p(nb)%a3d(it2,jt2,kt2)


!!				if( rm < rmin_limit .or. pm < pmin_limit)then
                  rm =  mb_r(nb)%a3d(it,jt,kt)
                  um =  mb_u(nb)%a3d(it,jt,kt) 
                  vm =  mb_v(nb)%a3d(it,jt,kt) 
                  wm =  mb_w(nb)%a3d(it,jt,kt) 
                  pm =  mb_p(nb)%a3d(it,jt,kt)
!!				endif
                mb_r(nb)%a3d(is,js,ks) = rm
                mb_u(nb)%a3d(is,js,ks) = um
                mb_v(nb)%a3d(is,js,ks) = vm
                mb_w(nb)%a3d(is,js,ks) = wm
                mb_p(nb)%a3d(is,js,ks) = pm

                call BC_outflow_turbulence(nb,is,js,ks,it,jt,kt) ! turbulence model


             enddo

        !* 以下另外再赋两层虚点值，特别注意，第三层密度为正
             is = is + s_lr3d(1)
             js = js + s_lr3d(2)
             ks = ks + s_lr3d(3)
             mb_r(nb)%a3d(is,js,ks) =   rm
             mb_u(nb)%a3d(is,js,ks) =   um
             mb_v(nb)%a3d(is,js,ks) =   vm
             mb_w(nb)%a3d(is,js,ks) =   wm
             mb_p(nb)%a3d(is,js,ks) =   pm

             is = is + s_lr3d(1)
             js = js + s_lr3d(2)
             ks = ks + s_lr3d(3)
             mb_r(nb)%a3d(is,js,ks) =   rm  ! 第3层密度为正
             mb_u(nb)%a3d(is,js,ks) =   um
             mb_v(nb)%a3d(is,js,ks) =   vm
             mb_w(nb)%a3d(is,js,ks) =   wm
             mb_p(nb)%a3d(is,js,ks) =   pm

         enddo
       enddo
    enddo
    return
end subroutine boundary_6
!_____________________________________________________________________!
subroutine boundary_8(nb,nr,bctype) !周期边界
!*tgh. 该程序只能处理周期边界的对应点在本块的情况

    use global_variables
    implicit none
    integer :: nwhole,nd,half_dim,whole_dim
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: nb,nr,bctype,m,n,nt(3),nsijk(3)
    integer :: is,js,ks,it,jt,kt,i,j,k
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: prim_s1(nl), q_1(nl), zf1

    t_lr = 0
    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
!      nsijk(m) = abs(mb_bc(nb)%bc(nr)%s_lr3d(m))
       nsijk(m) = abs( s_lr3d(m) )
       nt  (m) = 1 - nsijk(m) 
       if(s_lr3d(m) == 1) t_lr = 1
    enddo
    s_nd  = mb_bc(nb)%bc(nr)%s_nd     !边界面方向:1,2,3对应于i,j,k
    s_lr  = mb_bc(nb)%bc(nr)%s_lr     !左右边界-1,1对应于左右边界
    s_fix = mb_bc(nb)%bc(nr)%s_fix    !固定坐标(fixed_coor)
    nbs   = nb ! mb_bc(nb)%bc(nr)%nbs      !块号 
    if ( t_lr == 1 ) then
       whole_dim = s_st(s_nd) - 1 + method
    else
       whole_dim = s_ed(s_nd) + 1 - method
    endif
    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             do n=1,2-method
                is = i + s_lr3d(1) * n
                js = j + s_lr3d(2) * n
                ks = k + s_lr3d(3) * n
                it = nt(1) * is + nsijk(1) * ( whole_dim + s_lr3d(1) * n )
                jt = nt(2) * js + nsijk(2) * ( whole_dim + s_lr3d(2) * n )
                kt = nt(3) * ks + nsijk(3) * ( whole_dim + s_lr3d(3) * n )

                mb_r(nbs)%a3d(is,js,ks) =  mb_r(nbs)%a3d(it,jt,kt)
                mb_u(nbs)%a3d(is,js,ks) =  mb_u(nbs)%a3d(it,jt,kt) 
                mb_v(nbs)%a3d(is,js,ks) =  mb_v(nbs)%a3d(it,jt,kt) 
                mb_w(nbs)%a3d(is,js,ks) =  mb_w(nbs)%a3d(it,jt,kt) 
                mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)

					call BC_peroidic_turbulence(nbs,nbs,is,js,ks,it,jt,kt,i,j,k) ! turbulence model

             enddo

        !* 以下另外再赋两层虚点值，特别注意，第三层密度为正
             is = is + s_lr3d(1) 
             js = js + s_lr3d(2)
             ks = ks + s_lr3d(3)
             it = it + nsijk(1) * s_lr3d(1)
             jt = jt + nsijk(2) * s_lr3d(2)
             kt = kt + nsijk(3) * s_lr3d(3)
             mb_r(nbs)%a3d(is,js,ks) =  mb_r(nbs)%a3d(it,jt,kt)
             mb_u(nbs)%a3d(is,js,ks) =  mb_u(nbs)%a3d(it,jt,kt) 
             mb_v(nbs)%a3d(is,js,ks) =  mb_v(nbs)%a3d(it,jt,kt) 
             mb_w(nbs)%a3d(is,js,ks) =  mb_w(nbs)%a3d(it,jt,kt) 
             mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)
            
			 is = is + s_lr3d(1) 
             js = js + s_lr3d(2)
             ks = ks + s_lr3d(3)
             it = it + nsijk(1) * s_lr3d(1)
             jt = jt + nsijk(2) * s_lr3d(2)
             kt = kt + nsijk(3) * s_lr3d(3)
             mb_r(nbs)%a3d(is,js,ks) =  mb_r(nbs)%a3d(it,jt,kt)
             mb_u(nbs)%a3d(is,js,ks) =  mb_u(nbs)%a3d(it,jt,kt) 
             mb_v(nbs)%a3d(is,js,ks) =  mb_v(nbs)%a3d(it,jt,kt) 
             mb_w(nbs)%a3d(is,js,ks) =  mb_w(nbs)%a3d(it,jt,kt) 
             mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)

          enddo
       enddo
    enddo
    if ( method == 1 ) then !有限差分 给定边界上的值
       call dif_average(nb,nr,bctype,s_st,s_ed)
    endif
    return    
end subroutine boundary_8
!_____________________________________________________________________!
subroutine boundary_20_mml(nb,nr,bctype) !完全气体/完全非催化 等温边界
    use global_variables
    implicit none
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: nb,nr,bctype,m,n,nt
    integer :: is,js,ks,it,jt,kt,i,j,k
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: prim_s1(nl), q_1(nl), zf1
    real :: tw1,pres,dens,mav1,ff,tw2,rw,rw1,rw2

     if(nvis == 1 ) then
       do m=1,3
          s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
          s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
          s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
       enddo
       s_nd  = mb_bc(nb)%bc(nr)%s_nd     !边界面方向:1,2,3对应于i,j,k
       s_lr  = mb_bc(nb)%bc(nr)%s_lr     !左右边界-1,1对应于左右边界
       s_fix = mb_bc(nb)%bc(nr)%s_fix    !固定坐标(fixed_coor)
       nbs   = mb_bc(nb)%bc(nr)%nbs      !块号 
       tw2   = 0.50*mb_bc(nb)%bc(nr)%bc_par(1)/tref
       if( method == 0 ) then !有限体积
	       do i = s_st(1),s_ed(1)
              do j = s_st(2),s_ed(2)
                 do k = s_st(3),s_ed(3)
!   第一排给按对称条件延拓
                    is = i + s_lr3d(1)
                    js = j + s_lr3d(2)
                    ks = k + s_lr3d(3)

                    it = i 
                    jt = j 
                    kt = k 
                
                    do m=6,nl
                       prim_s1(m) = mb_fs(nb)%a4d(m-5,it,jt,kt)
                    enddo
                    pres = mb_p(nbs)%a3d(it,jt,kt)

					rw = 2.0*gama*moo*moo*pres/tw2

					rw1 = 2.0*rw - mb_r(nbs)%a3d(it,jt,kt)

                    mb_r(nbs)%a3d(is,js,ks) =  rw1
                    mb_u(nbs)%a3d(is,js,ks) = -mb_u(nbs)%a3d(it,jt,kt)
                    mb_v(nbs)%a3d(is,js,ks) = -mb_v(nbs)%a3d(it,jt,kt)
                    mb_w(nbs)%a3d(is,js,ks) = -mb_w(nbs)%a3d(it,jt,kt)
                    mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)
									  
										call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                    if ( nl > 5 ) then
                       zf1 = 0.0
                       do m=6,nl
                          zf1 = zf1 + prim_s1(m)
                       enddo
                       zf1 = min(1.0,zf1)
                       mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                    endif

                    do m=6,nl
                       mb_fs(nb)%a4d(m-5,is,js,ks) = prim_s1(m)  
                    enddo
                    is = i + 2*s_lr3d(1)
                    js = j + 2*s_lr3d(2)
                    ks = k + 2*s_lr3d(3)

                    it = i - s_lr3d(1)
                    jt = j - s_lr3d(2)
                    kt = k - s_lr3d(3)
                
                    do m=6,nl
                       prim_s1(m) = mb_fs(nb)%a4d(m-5,it,jt,kt)
                    enddo

					pres = mb_p(nbs)%a3d(it,jt,kt)

					rw = 2.0*gama*moo*moo*pres/tw2

                    rw2 = 2.0*rw - mb_r(nbs)%a3d(it,jt,kt)



                    mb_r(nbs)%a3d(is,js,ks) =  rw2
                    mb_u(nbs)%a3d(is,js,ks) = -mb_u(nbs)%a3d(it,jt,kt)
                    mb_v(nbs)%a3d(is,js,ks) = -mb_v(nbs)%a3d(it,jt,kt)
                    mb_w(nbs)%a3d(is,js,ks) = -mb_w(nbs)%a3d(it,jt,kt)
                    mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)
									  
										call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)
                    
!                给密度反号，表示为等温固壁
!                    mb_r(nbs)%a3d(is,js,ks) = -dens
                    if ( nl > 5 ) then
                       zf1 = 0.0
                       do m=6,nl
                          zf1 = zf1 + prim_s1(m)
                       enddo
                       zf1 = min(1.0,zf1)
                       mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                    endif

                    do m=6,nl
                       mb_fs(nb)%a4d(m-5,is,js,ks) = prim_s1(m)  
                    enddo
                 enddo
              enddo
           enddo
	   else  !有限差分 给定边界上的值
          do i = s_st(1),s_ed(1)
             do j = s_st(2),s_ed(2)
                do k = s_st(3),s_ed(3)
                   it = i - s_lr3d(1)
                   jt = j - s_lr3d(2)
                   kt = k - s_lr3d(3)

                   is = i 
                   js = j 
                   ks = k 

!                   pres = mb_p(nbs)%a3d(it,jt,kt)
                   pres = (4.D0*mb_p(nbs)%a3d(it,jt,kt) - mb_p(nbs)%a3d(it-s_lr3d(1),jt-s_lr3d(2),kt-s_lr3d(3)))/3.d0 
				           tw1  = 2.0*tw2 

                      if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                      else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                      endif

                   mb_r(nbs)%a3d(is,js,ks) =  dens
                   mb_u(nbs)%a3d(is,js,ks) =  0.0
                   mb_v(nbs)%a3d(is,js,ks) =  0.0
                   mb_w(nbs)%a3d(is,js,ks) =  0.0
                   mb_p(nbs)%a3d(is,js,ks) =  pres

									 call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                   do m=6,nl
                      mb_fs(nb)%a4d(m-5,is,js,ks) = mb_fs(nb)%a4d(m-5,it,jt,kt)  
                   enddo

                   if ( nl > 5 ) then
                      zf1 = 0.0
                      do m=6,nl
                         zf1 = zf1 + mb_fs(nb)%a4d(m-5,is,js,ks)
                      enddo
                      zf1 = min(1.0,zf1)
                      mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                   endif

                   is = i + s_lr3d(1) 
                   js = j + s_lr3d(2) 
                   ks = k + s_lr3d(3) 

				    tw1  = 4.0*tw2 - mb_t(nbs)%a3d(it,jt,kt)
				    if(tw1 < tw2 ) tw1 = tw2

                       if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       endif

                   mb_r(nbs)%a3d(is,js,ks) =  dens
                   mb_u(nbs)%a3d(is,js,ks) = -mb_u(nbs)%a3d(it,jt,kt)
                   mb_v(nbs)%a3d(is,js,ks) = -mb_v(nbs)%a3d(it,jt,kt)
                   mb_w(nbs)%a3d(is,js,ks) = -mb_w(nbs)%a3d(it,jt,kt)
                   mb_p(nbs)%a3d(is,js,ks) =  pres

									 call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                   do m=6,nl
                      mb_fs(nb)%a4d(m-5,is,js,ks) = mb_fs(nb)%a4d(m-5,it,jt,kt)  
                   enddo

                   if ( nl > 5 ) then
                      zf1 = 0.0
                      do m=6,nl
                         zf1 = zf1 + mb_fs(nb)%a4d(m-5,is,js,ks)
                      enddo
                      zf1 = min(1.0,zf1)
                      mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                   endif
!                给密度反号，表示为等温固壁
!                   mb_r(nbs)%a3d(i+2*s_lr3d(1),j+2*s_lr3d(2),k+2*s_lr3d(3) ) = -dens
                enddo
             enddo
          enddo
	   endif
	else
	   call boundary_2(nb,nr,bctype)
	endif

    return
end subroutine boundary_20_mml
!_____________________________________________________________________!
subroutine boundary_20(nb,nr,bctype) !完全气体/完全非催化 等温边界
    use global_variables
    implicit none
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: nb,nr,bctype,m,n,nt
    integer :: is,js,ks,it,jt,kt,i,j,k
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: prim_s1(nl), q_1(nl), zf1
    real :: tw1,pres,dens,mav1,ff,tw2

    if(nvis == 1 ) then
       do m=1,3
          s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
          s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
          s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
       enddo
       s_nd  = mb_bc(nb)%bc(nr)%s_nd     !边界面方向:1,2,3对应于i,j,k
       s_lr  = mb_bc(nb)%bc(nr)%s_lr     !左右边界-1,1对应于左右边界
       s_fix = mb_bc(nb)%bc(nr)%s_fix    !固定坐标(fixed_coor)
       nbs   = mb_bc(nb)%bc(nr)%nbs      !块号 
       tw2   = 0.50*mb_bc(nb)%bc(nr)%bc_par(1)/tref
       if( method == 0 ) then !有限体积
           do i = s_st(1),s_ed(1)
              do j = s_st(2),s_ed(2)
                 do k = s_st(3),s_ed(3)
!   第一排给按对称条件延拓
                    is = i + s_lr3d(1)
                    js = j + s_lr3d(2)
                    ks = k + s_lr3d(3)

                    it = i 
                    jt = j 
                    kt = k 
                
                    do m=6,nl
                       prim_s1(m) = mb_fs(nb)%a4d(m-5,it,jt,kt)
                    enddo
                    pres = mb_p(nbs)%a3d(it,jt,kt)
                    tw1  = 4.0*tw2 - mb_t(nbs)%a3d(it,jt,kt)
                    if(tw1 < tw2 ) tw1 = tw2

                       if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       endif

                    mb_r(nbs)%a3d(is,js,ks) =  dens
                    mb_u(nbs)%a3d(is,js,ks) = -mb_u(nbs)%a3d(it,jt,kt)
                    mb_v(nbs)%a3d(is,js,ks) = -mb_v(nbs)%a3d(it,jt,kt)
                    mb_w(nbs)%a3d(is,js,ks) = -mb_w(nbs)%a3d(it,jt,kt)
                    mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)

                    call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                    if ( nl > 5 ) then
                       zf1 = 0.0
                       do m=6,nl
                          zf1 = zf1 + prim_s1(m)
                       enddo
                       zf1 = min(1.0,zf1)
                       mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                    endif

                    do m=6,nl
                       mb_fs(nb)%a4d(m-5,is,js,ks) = prim_s1(m)
                    enddo
                    is = i + 2*s_lr3d(1)
                    js = j + 2*s_lr3d(2)
                    ks = k + 2*s_lr3d(3)

                    it = i - s_lr3d(1)
                    jt = j - s_lr3d(2)
                    kt = k - s_lr3d(3)
                
                    do m=6,nl
                       prim_s1(m) = mb_fs(nb)%a4d(m-5,it,jt,kt)
                    enddo
                    pres = mb_p(nbs)%a3d(it,jt,kt)
                    tw1  = 4.0*tw2 - mb_t(nbs)%a3d(it,jt,kt)
                    if(tw1 < tw2 ) tw1 = tw2
                       if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       endif

                    mb_r(nbs)%a3d(is,js,ks) =  dens
                    mb_u(nbs)%a3d(is,js,ks) = -mb_u(nbs)%a3d(it,jt,kt)
                    mb_v(nbs)%a3d(is,js,ks) = -mb_v(nbs)%a3d(it,jt,kt)
                    mb_w(nbs)%a3d(is,js,ks) = -mb_w(nbs)%a3d(it,jt,kt)
                    mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)

                    call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)
                    
!                给密度反号，表示为等温固壁
!                    mb_r(nbs)%a3d(is,js,ks) = -dens
                    if ( nl > 5 ) then
                       zf1 = 0.0
                       do m=6,nl
                          zf1 = zf1 + prim_s1(m)
                       enddo
                       zf1 = min(1.0,zf1)
                       mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                    endif

                    do m=6,nl
                       mb_fs(nb)%a4d(m-5,is,js,ks) = prim_s1(m)  
                    enddo
                 enddo
              enddo
           enddo
       else  !有限差分 给定边界上的值
          do i = s_st(1),s_ed(1)
             do j = s_st(2),s_ed(2)
                do k = s_st(3),s_ed(3)
                   it = i - s_lr3d(1)
                   jt = j - s_lr3d(2)
                   kt = k - s_lr3d(3)

                   is = i 
                   js = j 
                   ks = k 

!                   pres = mb_p(nbs)%a3d(it,jt,kt)
                   pres = (4.D0*mb_p(nbs)%a3d(it,jt,kt) - mb_p(nbs)%a3d(it-s_lr3d(1),jt-s_lr3d(2),kt-s_lr3d(3)))/3.d0 
                   tw1  = 2.0*tw2 
                      if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                      else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                      endif
                   mb_r(nbs)%a3d(is,js,ks) =  dens
                   mb_u(nbs)%a3d(is,js,ks) =  0.0
                   mb_v(nbs)%a3d(is,js,ks) =  0.0
                   mb_w(nbs)%a3d(is,js,ks) =  0.0
                   mb_p(nbs)%a3d(is,js,ks) =  pres		                           !*TGH old
!                   mb_p(nbs)%a3d(is,js,ks) =  0.5*(pres+dens*tw1/(gama*moo*moo) )  !*TGH 临时修改

                   call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                   do m=6,nl
                      mb_fs(nb)%a4d(m-5,is,js,ks) = mb_fs(nb)%a4d(m-5,it,jt,kt)
                   enddo

                   if ( nl > 5 ) then
                      zf1 = 0.0
                      do m=6,nl
                         zf1 = zf1 + mb_fs(nb)%a4d(m-5,is,js,ks)
                      enddo
                      zf1 = min(1.0,zf1)
                      mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                   endif

                   is = i + s_lr3d(1) 
                   js = j + s_lr3d(2) 
                   ks = k + s_lr3d(3) 

                   tw1  = 4.0*tw2 - mb_t(nbs)%a3d(it,jt,kt)
                   if(tw1 < tw2 ) tw1 = tw2

                       if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       endif

                   mb_r(nbs)%a3d(is,js,ks) =  dens
                   mb_u(nbs)%a3d(is,js,ks) = -mb_u(nbs)%a3d(it,jt,kt)
                   mb_v(nbs)%a3d(is,js,ks) = -mb_v(nbs)%a3d(it,jt,kt)
                   mb_w(nbs)%a3d(is,js,ks) = -mb_w(nbs)%a3d(it,jt,kt)
                   mb_p(nbs)%a3d(is,js,ks) =  pres

                   call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                   do m=6,nl
                      mb_fs(nb)%a4d(m-5,is,js,ks) = mb_fs(nb)%a4d(m-5,it,jt,kt)  
                   enddo

                   if ( nl > 5 ) then
                      zf1 = 0.0
                      do m=6,nl
                         zf1 = zf1 + mb_fs(nb)%a4d(m-5,is,js,ks)
                      enddo
                      zf1 = min(1.0,zf1)
                      mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                   endif
!                给密度反号，表示为等温固壁
!                   mb_r(nbs)%a3d(i+2*s_lr3d(1),j+2*s_lr3d(2),k+2*s_lr3d(3) ) = -dens
                enddo
             enddo
          enddo
       endif
    else
       call inviscid_wall(nb,nr,bctype)
    endif

    return
end subroutine boundary_20
!_____________________________________________________________________!
subroutine boundary_21(nb,nr,bctype) !完全气体恒温壁
    implicit none
    integer :: nb,nr,bctype
    return
end subroutine boundary_21
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine viscid_wall(nb,nr,bctype) !粘性壁面边界
    use define_precision_mod
    use global_const,only : pmin_limit,tmin_limit
    use global_variables,only : mb_bc,mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,moo,gama,twall,tref
	implicit none
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nb,nr,bctype,m,n
    integer :: is,js,ks,it,jt,kt,i,j,k
    integer :: s_st(3),s_ed(3),s_lr3d(3)
	real(prec) :: moocp,dres,pres,tres,twall_present

	moocp=moo*moo*gama

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)  
    enddo

	if(twall > 0.0)then !等温壁
		twall_present = twall/tref
        do k = s_st(3),s_ed(3)
        do j = s_st(2),s_ed(2)
        do i = s_st(1),s_ed(1)
                is = i + s_lr3d(1)
                js = j + s_lr3d(2)
                ks = k + s_lr3d(3)
                it = i - s_lr3d(1)
                jt = j - s_lr3d(2)
                kt = k - s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) =   mb_r(nb)%a3d(it,jt,kt)
                mb_u(nb)%a3d(is,js,ks) =  -mb_u(nb)%a3d(it,jt,kt)
                mb_v(nb)%a3d(is,js,ks) =  -mb_v(nb)%a3d(it,jt,kt)
                mb_w(nb)%a3d(is,js,ks) =  -mb_w(nb)%a3d(it,jt,kt) 
                mb_p(nb)%a3d(is,js,ks) =   mb_p(nb)%a3d(it,jt,kt)

							 call BC_wall_turbulence(nb,is,js,ks,it,jt,kt,i,j,k)!*TGH. 湍流模型赋虚点值


!*tgh.
!*tgh. 2阶
               pres = (4.D0*mb_p(nb)%a3d(it,jt,kt) - mb_p(nb)%a3d(it-s_lr3d(1),jt-s_lr3d(2),kt-s_lr3d(3)))/3.d0
			   if(pres < pmin_limit) pres = mb_p(nb)%a3d(it,jt,kt)
!*tgh. 以上 2 阶
!*tgh. 以下 1 阶
!				pres= mb_p(nb)%a3d(it,jt,kt)
!*TGH. END 1 阶

!			    pres = max(pres,1.0d-5)
			    tres= twall_present
			    dres= moocp*pres/tres

                   mb_r(nb)%a3d(i,j,k) =  dres
                   mb_u(nb)%a3d(i,j,k) =  0.0
                   mb_v(nb)%a3d(i,j,k) =  0.0
                   mb_w(nb)%a3d(i,j,k) =  0.0
                   mb_p(nb)%a3d(i,j,k) =  pres
								call BC_face_dif_turbulence(nb,nb,is,js,ks,it,jt,kt,i,j,k) !*TGH. 求边界上的值
!*tgh. end

            !* 以下另外再赋两层虚点值，特别注意，第三层密度给的是负值，在使用时需要乘“-1”
                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
                it = it - s_lr3d(1)
                jt = jt - s_lr3d(2)
                kt = kt - s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) =   mb_r(nb)%a3d(it,jt,kt)
                mb_u(nb)%a3d(is,js,ks) =  -mb_u(nb)%a3d(it,jt,kt)
                mb_v(nb)%a3d(is,js,ks) =  -mb_v(nb)%a3d(it,jt,kt)
                mb_w(nb)%a3d(is,js,ks) =  -mb_w(nb)%a3d(it,jt,kt) 
                mb_p(nb)%a3d(is,js,ks) =   mb_p(nb)%a3d(it,jt,kt)

                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
                it = it - s_lr3d(1)
                jt = jt - s_lr3d(2)
                kt = kt - s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) =  -mb_r(nb)%a3d(it,jt,kt)   !第三层密度给负值
                mb_u(nb)%a3d(is,js,ks) =  -mb_u(nb)%a3d(it,jt,kt)
                mb_v(nb)%a3d(is,js,ks) =  -mb_v(nb)%a3d(it,jt,kt)
                mb_w(nb)%a3d(is,js,ks) =  -mb_w(nb)%a3d(it,jt,kt) 
                mb_p(nb)%a3d(is,js,ks) =   mb_p(nb)%a3d(it,jt,kt)
        enddo
        enddo
        enddo


	else !绝热壁

        do k = s_st(3),s_ed(3)
        do j = s_st(2),s_ed(2)
        do i = s_st(1),s_ed(1)
                is = i + s_lr3d(1)
                js = j + s_lr3d(2)
                ks = k + s_lr3d(3)
                it = i - s_lr3d(1)
                jt = j - s_lr3d(2)
                kt = k - s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) =   mb_r(nb)%a3d(it,jt,kt)
                mb_u(nb)%a3d(is,js,ks) =  -mb_u(nb)%a3d(it,jt,kt)
                mb_v(nb)%a3d(is,js,ks) =  -mb_v(nb)%a3d(it,jt,kt)
                mb_w(nb)%a3d(is,js,ks) =  -mb_w(nb)%a3d(it,jt,kt) 
                mb_p(nb)%a3d(is,js,ks) =   mb_p(nb)%a3d(it,jt,kt)

							 call BC_wall_turbulence(nb,is,js,ks,it,jt,kt,i,j,k)!*TGH. 湍流模型赋虚点值


!*tgh.
!*tgh. 2阶
               pres = (4.D0*mb_p(nb)%a3d(it,jt,kt) - mb_p(nb)%a3d(it-s_lr3d(1),jt-s_lr3d(2),kt-s_lr3d(3)))/3.d0
               tres = (4.D0*mb_t(nb)%a3d(it,jt,kt) - mb_t(nb)%a3d(it-s_lr3d(1),jt-s_lr3d(2),kt-s_lr3d(3)))/3.d0
			    if(pres < pmin_limit) pres = mb_p(nb)%a3d(it,jt,kt)
			    if(tres < tmin_limit) tres = mb_t(nb)%a3d(it,jt,kt)
!*tgh. 以上 2 阶
!*tgh. 以下 1 阶
!				pres= mb_p(nb)%a3d(it,jt,kt)
!				dres= mb_r(nb)%a3d(it,jt,kt)
!				tres= moocp*pres/dres
!*TGH. END 1 阶
!			   pres = max(pres,1.0d-5)
!			   tres = max(tres,1.0d-3)			   
			   dres = moocp*pres/tres

                   mb_r(nb)%a3d(i,j,k) =  dres
                   mb_u(nb)%a3d(i,j,k) =  0.0
                   mb_v(nb)%a3d(i,j,k) =  0.0
                   mb_w(nb)%a3d(i,j,k) =  0.0
                   mb_p(nb)%a3d(i,j,k) =  pres
								call BC_face_dif_turbulence(nb,nb,is,js,ks,it,jt,kt,i,j,k) !*TGH. 求边界上的值
!*tgh. end

            !* 以下另外再赋两层虚点值，特别注意，第三层密度给的是负值，在使用时需要乘“-1”
                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
                it = it - s_lr3d(1)
                jt = jt - s_lr3d(2)
                kt = kt - s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) =   mb_r(nb)%a3d(it,jt,kt)
                mb_u(nb)%a3d(is,js,ks) =  -mb_u(nb)%a3d(it,jt,kt)
                mb_v(nb)%a3d(is,js,ks) =  -mb_v(nb)%a3d(it,jt,kt)
                mb_w(nb)%a3d(is,js,ks) =  -mb_w(nb)%a3d(it,jt,kt) 
                mb_p(nb)%a3d(is,js,ks) =   mb_p(nb)%a3d(it,jt,kt)
                
                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
                it = it - s_lr3d(1)
                jt = jt - s_lr3d(2)
                kt = kt - s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) =  -mb_r(nb)%a3d(it,jt,kt)   !第三层密度给负值
                mb_u(nb)%a3d(is,js,ks) =  -mb_u(nb)%a3d(it,jt,kt)
                mb_v(nb)%a3d(is,js,ks) =  -mb_v(nb)%a3d(it,jt,kt)
                mb_w(nb)%a3d(is,js,ks) =  -mb_w(nb)%a3d(it,jt,kt) 
                mb_p(nb)%a3d(is,js,ks) =   mb_p(nb)%a3d(it,jt,kt)

        enddo
        enddo
        enddo

	endif

    return
end subroutine viscid_wall
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine inviscid_wall(nb,nr,bctype) !无粘壁面边界
    use define_precision_mod
    use global_variables,only : mb_bc,mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,moo,gama &
	                          & ,U,V,W,METHOD
	implicit none
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: nb,nr,bctype,m,n,nt
    integer :: is,js,ks,it,jt,kt,i,j,k,iq,jq,kq
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real(PREC) :: nx,ny,nz,SSS1

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)  
    enddo

    t_nd  = s_lr3d(1)
    t_fix = s_lr3d(2)
    t_lr  = s_lr3d(3)

    do k = s_st(3),s_ed(3)
       do j = s_st(2),s_ed(2)
          do i = s_st(1),s_ed(1)
             call getnxyz_mml(nx,ny,nz,i,j,k,t_nd,t_fix,t_lr)             
                is = i + s_lr3d(1)
                js = j + s_lr3d(2)
                ks = k + s_lr3d(3)

                it = i - s_lr3d(1)
                jt = j - s_lr3d(2)
                kt = k - s_lr3d(3)
                sss1 =2.0* ( u(it,jt,kt)*nx + v(it,jt,kt)*ny + w(it,jt,kt)*nz )
                mb_r(nb)%a3d(is,js,ks) =  mb_r(nb)%a3d(it,jt,kt)
                mb_u(nb)%a3d(is,js,ks) =  mb_u(nb)%a3d(it,jt,kt) - sss1*nx
                mb_v(nb)%a3d(is,js,ks) =  mb_v(nb)%a3d(it,jt,kt) - sss1*ny
                mb_w(nb)%a3d(is,js,ks) =  mb_w(nb)%a3d(it,jt,kt) - sss1*nz
                mb_p(nb)%a3d(is,js,ks) =  mb_p(nb)%a3d(it,jt,kt)
 
            !* 以下另外再赋两层虚点值，特别注意，第三层密度给的是负值，在使用时需要乘“-1”
                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
                it = it - s_lr3d(1)
                jt = jt - s_lr3d(2)
                kt = kt - s_lr3d(3)
                sss1 =2.0* ( u(it,jt,kt)*nx + v(it,jt,kt)*ny + w(it,jt,kt)*nz )
                mb_r(nb)%a3d(is,js,ks) =  mb_r(nb)%a3d(it,jt,kt)
                mb_u(nb)%a3d(is,js,ks) =  mb_u(nb)%a3d(it,jt,kt) - sss1*nx
                mb_v(nb)%a3d(is,js,ks) =  mb_v(nb)%a3d(it,jt,kt) - sss1*ny
                mb_w(nb)%a3d(is,js,ks) =  mb_w(nb)%a3d(it,jt,kt) - sss1*nz
                mb_p(nb)%a3d(is,js,ks) =  mb_p(nb)%a3d(it,jt,kt)

                is = is + s_lr3d(1)
                js = js + s_lr3d(2)
                ks = ks + s_lr3d(3)
                it = it - s_lr3d(1)
                jt = jt - s_lr3d(2)
                kt = kt - s_lr3d(3)
                sss1 =2.0* ( u(it,jt,kt)*nx + v(it,jt,kt)*ny + w(it,jt,kt)*nz )
                mb_r(nb)%a3d(is,js,ks) = -mb_r(nb)%a3d(it,jt,kt)   !第三层密度给负值
                mb_u(nb)%a3d(is,js,ks) =  mb_u(nb)%a3d(it,jt,kt) - sss1*nx
                mb_v(nb)%a3d(is,js,ks) =  mb_v(nb)%a3d(it,jt,kt) - sss1*ny
                mb_w(nb)%a3d(is,js,ks) =  mb_w(nb)%a3d(it,jt,kt) - sss1*nz
                mb_p(nb)%a3d(is,js,ks) =  mb_p(nb)%a3d(it,jt,kt)
          enddo
       enddo
    enddo

    if ( method == 1 ) then !有限差分 给定边界上的值
       call dif_average(nb,nr,bctype,s_st,s_ed)
    endif    
    return
end subroutine inviscid_wall
!_____________________________________________________________________!
subroutine boundary_7(nb,nr,bctype) !奇性面
    use global_variables
    implicit none
    integer :: nwhole,half_dim,whole_dim
    integer :: nb,nr,naxir,bcaxi,m,n
    integer :: bctype,s_nd,nd,nd3
    integer :: i,j,k,is1,it1,jj,is,js,ks,it,jt,kt
    integer :: s_lr3d(3),s_st(3),s_ed(3)
    real :: nx,ny,nz,prim_s1(nl),zf1,q_1(nl),cn,prim_t1(nl),sss1,un

    call  boundary_6(nb,nr,bctype)
    return
end subroutine boundary_7
!_____________________________________________________________________!
subroutine boundary_7123_old(nb,nr,bctype) !奇性轴
    use global_variables
    implicit none
    integer :: nwhole,half_dim,whole_dim
    integer :: nb,nr,naxir,bcaxi,m,n
    integer :: bctype,s_nd,nd,nd3
    integer :: i,j,k,is1,it1,jj,is,js,ks,it,jt,kt
    integer :: s_lr3d(3),s_st(3),s_ed(3)
    real :: nx,ny,nz,prim_s1(nl),zf1,q_1(nl),cn,prim_t1(nl),sss1,un
    real :: turb_s1(nlamtur),turb_t1(nlamtur)
	real :: rm,um,vm,wm,pm

    !寻找算奇轴为半奇轴还是全奇轴的判据
    nd = 0                   !赋初值，假定为i向
    nwhole = 0
    s_nd  = mb_bc(nb)%bc(nr)%s_nd
    do naxir = 1,mb_bc(nb)%nregions
       bcaxi = mb_bc(nb)%bc(naxir)%bctype
       if ( bcaxi == 3 ) then
          nwhole = 0         !存在对称边界条件，改算半场
          nd = mb_bc(nb)%bc(naxir)%s_nd
          if ( s_nd /= nd ) goto 10
       endif
       if ( bcaxi == 8 ) then
          nwhole = 1         !存在周期边界条件，改算全场
          nd = mb_bc(nb)%bc(naxir)%s_nd
          if ( s_nd /= nd ) goto 10
       endif
       if ( bcaxi < 0 ) then   !块交界面边界
          if ( nb == mb_bc(nb)%bc(naxir)%nbt ) then
             if ( mb_bc(nb)%bc(naxir)%t_nd == mb_bc(nb)%bc(naxir)%s_nd ) then
                nwhole = 1         !存在周期边界条件，改算全场
                nd = mb_bc(nb)%bc(naxir)%s_nd
                if ( s_nd /= nd ) goto 10
             endif
          endif
       endif
    enddo
10  continue

    if ( nd == 0 ) then
        call boundary_7(nb,nr,bctype)
        goto 20
    endif
    
    nd3 = 6 - s_nd - nd
    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)  = mb_bc(nb)%bc(nr)%s_lr3d(m)
    enddo

    if ( nwhole == 0 ) then
       !确定对称面的法向
       if ( nd == 1 ) then
          i = s_st(1)
          j = (s_ed(2) + s_st(2))/2 - 3*s_lr3d(2)
          k = (s_ed(3) + s_st(3))/2 - 3*s_lr3d(3)
          call  getnxyz_mml(nx,ny,nz,i,j,k,1,0,0)      !计算单位法向矢量
       elseif ( nd == 2 ) then
          i = (s_ed(1) + s_st(1))/2 - 3*s_lr3d(1)
          j = s_st(2)
          k = (s_ed(3) + s_st(3))/2 - 3*s_lr3d(3)
          call  getnxyz_mml(nx,ny,nz,i,j,k,0,1,0)      !计算单位法向矢量
       else
          i = (s_ed(1) + s_st(1))/2 - 3*s_lr3d(1)
          j = (s_ed(2) + s_st(2))/2 - 3*s_lr3d(2)
          k =  s_st(3)
          call  getnxyz_mml(nx,ny,nz,i,j,k,0,0,1)      !计算单位法向矢量
       endif
    else
       nx = 0.0
       ny = 0.0
       nz = 0.0
    endif

    half_dim  = s_st(nd) + ( s_ed(nd) - s_st(nd) )/( 1 + nwhole ) - method 
    whole_dim = s_ed(nd) 

    do i = s_st(s_nd),s_ed(s_nd)
       do n=1,2-method
          is1 = i + n * mb_bc(nb)%bc(nr)%s_lr
          it1 = i - ( n - 1 + method ) * mb_bc(nb)%bc(nr)%s_lr
          do j = s_st(nd),s_ed(nd)
             if ( nwhole == 0 ) then
                jj = whole_dim - j + 1 
             else
                jj = half_dim + j 
                if( jj > whole_dim ) jj = jj - whole_dim + method 
             endif
             do k = s_st(nd3),s_ed(nd3)
                is = iii(s_nd,1)*is1 + iii(nd,1)*j  + iii(nd3,1)*k
                js = iii(s_nd,2)*is1 + iii(nd,2)*j  + iii(nd3,2)*k
                ks = iii(s_nd,3)*is1 + iii(nd,3)*j  + iii(nd3,3)*k
                it = iii(s_nd,1)*it1 + iii(nd,1)*jj + iii(nd3,1)*k
                jt = iii(s_nd,2)*it1 + iii(nd,2)*jj + iii(nd3,2)*k
                kt = iii(s_nd,3)*it1 + iii(nd,3)*jj + iii(nd3,3)*k
                sss1 = 2.0 * ( mb_u(nb)%a3d(it,jt,kt)*nx + &
                               mb_v(nb)%a3d(it,jt,kt)*ny + &
                               mb_w(nb)%a3d(it,jt,kt)*nz )
                mb_r(nb)%a3d(is,js,ks) =  mb_r(nb)%a3d(it,jt,kt)
                mb_u(nb)%a3d(is,js,ks) =  mb_u(nb)%a3d(it,jt,kt) - sss1*nx
                mb_v(nb)%a3d(is,js,ks) =  mb_v(nb)%a3d(it,jt,kt) - sss1*ny
                mb_w(nb)%a3d(is,js,ks) =  mb_w(nb)%a3d(it,jt,kt) - sss1*nz
                mb_p(nb)%a3d(is,js,ks) =  mb_p(nb)%a3d(it,jt,kt)
!                给密度反号，表示为奇性轴
!				        if( n==2 ) mb_r(nb)%a3d(is,js,ks) = mb_r(nb)%a3d(it,jt,kt)
                do m = 1,nl-5    !ns
                   mb_fs(nb)%a4d( m,is,js,ks) = mb_fs(nb)%a4d(m,it,jt,kt)
                enddo

				do m=1,nlamtur
				   mb_qke(nb)%a4d(is,js,ks,m) =  mb_qke(nb)%a4d(it,jt,kt,m)								
				enddo
                do m=6,nl
                   prim_s1(m) = mb_fs(nb)%a4d(m-5,is,js,ks)
                enddo


                if ( nl > 5 ) then
                   zf1 = 0.0
                   do m=6,nl
                      zf1 = zf1 + prim_s1(m)
                   enddo
                   zf1 = min(1.0,zf1)
                   mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                endif
             enddo
          enddo
       enddo
    enddo

!    call boundary_7(nb,nr,bctype)


    if ( method == 1 ) then !有限差分 给定边界上的值
       cn = 1.0/(s_ed(nd) - s_st(nd) + 1)
       do i = s_st(s_nd),s_ed(s_nd)
          is1 = i + mb_bc(nb)%bc(nr)%s_lr
          it1 = i - mb_bc(nb)%bc(nr)%s_lr
          do j = s_st(nd3),s_ed(nd3)
             do m=1,nl
                prim_s1(m) = 0.0
                prim_t1(m) = 0.0
             enddo

             do m=1,nlamtur
                turb_s1(m) = 0.0
                turb_t1(m) = 0.0
             enddo

             do k = s_st(nd),s_ed(nd)
                is = iii(s_nd,1)*is1 + iii(nd,1)*k  + iii(nd3,1)*j
                js = iii(s_nd,2)*is1 + iii(nd,2)*k  + iii(nd3,2)*j
                ks = iii(s_nd,3)*is1 + iii(nd,3)*k  + iii(nd3,3)*j
             
                do m=1,nl
                   prim_t1(m) = prim_t1(m) + mb_q(nb)%a4d(m,is,js,ks)
                enddo
             
                do m=1,nlamtur
                   turb_t1(m) = turb_t1(m) + mb_qke(nb)%a4d(is,js,ks,m)
                enddo
             
                prim_s1(1) = prim_s1(1) +  mb_r(nb)%a3d(is,js,ks)
                prim_s1(2) = prim_s1(2) +  mb_u(nb)%a3d(is,js,ks)
                prim_s1(3) = prim_s1(3) +  mb_v(nb)%a3d(is,js,ks)
                prim_s1(4) = prim_s1(4) +  mb_w(nb)%a3d(is,js,ks)
                prim_s1(5) = prim_s1(5) +  mb_p(nb)%a3d(is,js,ks)

                do m=1,nlamtur
                   turb_s1(m) = turb_s1(m) + mb_qke(nb)%a4d(is,js,ks,m)
                enddo
             
            
                do m=6,nl
                   prim_s1(m) = prim_s1(m) + mb_fs(nb)%a4d(m-5,is,js,ks) 
                enddo
             enddo
             do m=1,nl
              prim_s1(m) = prim_s1(m) * cn
             enddo

             do m=1,nlamtur
                turb_s1(m) = turb_s1(m) * cn
             enddo
             
             un = nx*prim_s1(2) + ny*prim_s1(3) + nz*prim_s1(4)
             prim_s1(2) = prim_s1(2) - nx * un
             prim_s1(3) = prim_s1(3) - ny * un
             prim_s1(4) = prim_s1(4) - nz * un

             do m=1,nl
                prim_t1(m) = prim_t1(m) * cn
             enddo

             do m=1,nlamtur
                turb_t1(m) = turb_t1(m) * cn
             enddo
             
             un = nx*prim_t1(2) + ny*prim_t1(3) + nz*prim_t1(4)
             prim_t1(2) = prim_t1(2) - nx * un
             prim_t1(3) = prim_t1(3) - ny * un
             prim_t1(4) = prim_t1(4) - nz * un
             
             do k = s_st(nd),s_ed(nd)
                is = iii(s_nd,1)*i + iii(nd,1)*k  + iii(nd3,1)*j
                js = iii(s_nd,2)*i + iii(nd,2)*k  + iii(nd3,2)*j
                ks = iii(s_nd,3)*i + iii(nd,3)*k  + iii(nd3,3)*j
                do m=1,nl
                   mb_q(nb)%a4d(m,is,js,ks) = prim_t1(m) 
                enddo
				
				if(prim_s1(1) < rmin_limit)prim_s1(1) = rmin_limit             
				if(prim_s1(5) < pmin_limit)prim_s1(5) = pmin_limit
				          
                mb_r(nb)%a3d(is,js,ks) = prim_s1(1)  
                mb_u(nb)%a3d(is,js,ks) = prim_s1(2)  
                mb_v(nb)%a3d(is,js,ks) = prim_s1(3)  
                mb_w(nb)%a3d(is,js,ks) = prim_s1(4)  
                mb_p(nb)%a3d(is,js,ks) = prim_s1(5)  
             
								do m=1,nlamtur
									mb_qke(nb)%a4d(is,js,ks,m) = turb_s1(m)
								enddo
             
                do m=6,nl
                   mb_fs(nb)%a4d(m-5,is,js,ks) = prim_s1(m)  
                enddo


                if ( nl > 5 ) then
                    zf1 = 0.0
                    do m=6,nl
                       mb_fs(nb)%a4d(m-5,is,js,ks) = prim_s1(m)  
                       zf1 = zf1 + prim_s1(m)
                    enddo
                    zf1 = min(1.0,zf1)
                    mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                endif

!********************
!* 以下另外再赋两层虚点值，特别注意，第三层密度给的是负值，在使用时需要乘“-1”
! 还没有检验虚层上值的正确性
! 
				is = is + s_lr3d(1)
				js = js + s_lr3d(2)
				ks = ks + s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) = prim_s1(1)  
                mb_u(nb)%a3d(is,js,ks) = prim_s1(2)  
                mb_v(nb)%a3d(is,js,ks) = prim_s1(3)  
                mb_w(nb)%a3d(is,js,ks) = prim_s1(4)  
                mb_p(nb)%a3d(is,js,ks) = prim_s1(5)

				is = is + s_lr3d(1)
				js = js + s_lr3d(2)
				ks = ks + s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) = prim_s1(1)  
                mb_u(nb)%a3d(is,js,ks) = prim_s1(2)  
                mb_v(nb)%a3d(is,js,ks) = prim_s1(3)  
                mb_w(nb)%a3d(is,js,ks) = prim_s1(4)  
                mb_p(nb)%a3d(is,js,ks) = prim_s1(5)
				 
				is = is + s_lr3d(1)
				js = js + s_lr3d(2)
				ks = ks + s_lr3d(3)
                mb_r(nb)%a3d(is,js,ks) =-prim_s1(1)  
                mb_u(nb)%a3d(is,js,ks) = prim_s1(2)  
                mb_v(nb)%a3d(is,js,ks) = prim_s1(3)  
                mb_w(nb)%a3d(is,js,ks) = prim_s1(4)  
                mb_p(nb)%a3d(is,js,ks) = prim_s1(5)  

             enddo
          enddo
       enddo
    endif
20  continue
    return
end subroutine boundary_7123_old
!_____________________________________________________________________!
subroutine boundary_7123(nb,nr,bctype) !奇性轴
    use global_variables
    implicit none
    integer :: nwhole,half_dim,whole_dim
    integer :: nb,nr,naxir,bcaxi,m,n
    integer :: bctype,s_nd,nd,nd3
    integer :: i,j,k,is1,is2,it1,it2,jj,is,js,ks,it,jt,kt
    integer :: s_lr3d(3),s_st(3),s_ed(3)
    real :: nx,ny,nz,cn,sss1,un,zf1,c_q1,c_q2
    real :: prim_s1(nl),prim_s2(nl),q_1(nl),prim_t1(nl),prim_t2(nl)

    !寻找算奇轴为半奇轴还是全奇轴的判据
    nd = 0                   !赋初值，假定为i向
    nwhole = 0
    s_nd  = mb_bc(nb)%bc(nr)%s_nd
    do naxir = 1,mb_bc(nb)%nregions
       bcaxi = mb_bc(nb)%bc(naxir)%bctype
       if ( bcaxi == 3 ) then
          nwhole = 0         !存在对称边界条件，改算半场
          nd = mb_bc(nb)%bc(naxir)%s_nd
          if ( s_nd /= nd ) goto 10
       endif
       if ( bcaxi == 8 ) then
          nwhole = 1         !存在周期边界条件，改算全场
          nd = mb_bc(nb)%bc(naxir)%s_nd
          if ( s_nd /= nd ) goto 10
       endif
       if ( bcaxi < 0 ) then   !块交界面边界
          if ( nb == mb_bc(nb)%bc(naxir)%nbt ) then
             if ( mb_bc(nb)%bc(naxir)%t_nd == mb_bc(nb)%bc(naxir)%s_nd ) then
                nwhole = 1         !存在周期边界条件，改算全场
                nd = mb_bc(nb)%bc(naxir)%s_nd
                if ( s_nd /= nd ) goto 10
             endif
          endif
       endif
    enddo
10  continue

    if ( nd == 0 ) then
        call boundary_7(nb,nr,bctype)
        goto 20
    endif
    
    nd3 = 6 - s_nd - nd
    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)  = mb_bc(nb)%bc(nr)%s_lr3d(m)
    enddo

    if ( nwhole == 0 ) then
       !确定对称面的法向
       if ( nd == 1 ) then
          i = s_st(1)
          j = (s_ed(2) + s_st(2))/2 - 3*s_lr3d(2)
          k = (s_ed(3) + s_st(3))/2 - 3*s_lr3d(3)
          call  getnxyz_mml(nx,ny,nz,i,j,k,1,0,0)      !计算单位法向矢量
       elseif ( nd == 2 ) then
          i = (s_ed(1) + s_st(1))/2 - 3*s_lr3d(1)
          j = s_st(2)
          k = (s_ed(3) + s_st(3))/2 - 3*s_lr3d(3)
          call  getnxyz_mml(nx,ny,nz,i,j,k,0,1,0)      !计算单位法向矢量
       else
          i = (s_ed(1) + s_st(1))/2 - 3*s_lr3d(1)
          j = (s_ed(2) + s_st(2))/2 - 3*s_lr3d(2)
          k =  s_st(3)
          call  getnxyz_mml(nx,ny,nz,i,j,k,0,0,1)      !计算单位法向矢量
       endif
    else
       nx = 0.0
       ny = 0.0
       nz = 0.0
    endif

    half_dim  = s_st(nd) + ( s_ed(nd) - s_st(nd) )/( 1 + nwhole ) - method 
    whole_dim = s_ed(nd) 
!!!!  确定奇性轴上流场
       if(method == 1 ) then
	      c_q1 = 4.0/3.0
		  c_q2 =-1.0/3.0
	   else
	      c_q1 = 9.0/8.0
		  c_q2 =-1.0/8.0
	   endif
       cn = 1.0/(s_ed(nd) - s_st(nd) + 1)
       do i = s_st(s_nd),s_ed(s_nd)
!          is1 = i   + mb_bc(nb)%bc(nr)%s_lr
          it1 = i   - mb_bc(nb)%bc(nr)%s_lr
          it2 = it1 - mb_bc(nb)%bc(nr)%s_lr
          do j = s_st(nd3),s_ed(nd3)
             do m=1,nl
                prim_s1(m) = 0.0
                prim_s2(m) = 0.0 ! prim_s1(nl),prim_s2(nl)
             enddo
             do k = s_st(nd),s_ed(nd)
                is = iii(s_nd,1)*it1 + iii(nd,1)*k  + iii(nd3,1)*j
                js = iii(s_nd,2)*it1 + iii(nd,2)*k  + iii(nd3,2)*j
                ks = iii(s_nd,3)*it1 + iii(nd,3)*k  + iii(nd3,3)*j
                prim_s1(1) = prim_s1(1) +  mb_r(nb)%a3d(is,js,ks)
                prim_s1(2) = prim_s1(2) +  mb_u(nb)%a3d(is,js,ks)
                prim_s1(3) = prim_s1(3) +  mb_v(nb)%a3d(is,js,ks)
                prim_s1(4) = prim_s1(4) +  mb_w(nb)%a3d(is,js,ks)
                prim_s1(5) = prim_s1(5) +  mb_p(nb)%a3d(is,js,ks)
                do m=6,nl
                   prim_s1(m) = prim_s1(m) + mb_fs(nb)%a4d(m-5,is,js,ks) 
                enddo

                is = iii(s_nd,1)*it2 + iii(nd,1)*k  + iii(nd3,1)*j
                js = iii(s_nd,2)*it2 + iii(nd,2)*k  + iii(nd3,2)*j
                ks = iii(s_nd,3)*it2 + iii(nd,3)*k  + iii(nd3,3)*j
                prim_s2(1) = prim_s2(1) +  mb_r(nb)%a3d(is,js,ks)
                prim_s2(2) = prim_s2(2) +  mb_u(nb)%a3d(is,js,ks)
                prim_s2(3) = prim_s2(3) +  mb_v(nb)%a3d(is,js,ks)
                prim_s2(4) = prim_s2(4) +  mb_w(nb)%a3d(is,js,ks)
                prim_s2(5) = prim_s2(5) +  mb_p(nb)%a3d(is,js,ks)
                do m=6,nl
                   prim_s2(m) = prim_s2(m) + mb_fs(nb)%a4d(m-5,is,js,ks) 
                enddo
             enddo
             do m=1,nl
                prim_s1(m) = prim_s1(m) * cn
                prim_s2(m) = prim_s2(m) * cn
             enddo
             un = nx*prim_s1(2) + ny*prim_s1(3) + nz*prim_s1(4)
             prim_s1(2) = prim_s1(2) - nx * un
             prim_s1(3) = prim_s1(3) - ny * un
             prim_s1(4) = prim_s1(4) - nz * un

             un = nx*prim_s2(2) + ny*prim_s2(3) + nz*prim_s2(4)
             prim_s2(2) = prim_s2(2) - nx * un
             prim_s2(3) = prim_s2(3) - ny * un
             prim_s2(4) = prim_s2(4) - nz * un

             do m=1,nl
			    q_1(m) = c_q1*prim_s1(m) + c_q2*prim_s2(m)
             enddo
             do n=1,2-method
                is1 = i + n * mb_bc(nb)%bc(nr)%s_lr
                it1 = i - ( n - 1 + method ) * mb_bc(nb)%bc(nr)%s_lr
                do k = s_st(nd),s_ed(nd)
                   is = iii(s_nd,1)*is1 + iii(nd,1)*k  + iii(nd3,1)*j
                   js = iii(s_nd,2)*is1 + iii(nd,2)*k  + iii(nd3,2)*j
                   ks = iii(s_nd,3)*is1 + iii(nd,3)*k  + iii(nd3,3)*j
                   it = iii(s_nd,1)*i   + iii(nd,1)*k  + iii(nd3,1)*j
                   jt = iii(s_nd,2)*i   + iii(nd,2)*k  + iii(nd3,2)*j
                   kt = iii(s_nd,3)*i   + iii(nd,3)*k  + iii(nd3,3)*j
                   mb_r(nb)%a3d(is,js,ks) =  2.0*q_1(1) - mb_r(nb)%a3d(it,jt,kt)
                   mb_u(nb)%a3d(is,js,ks) =  2.0*q_1(2) - mb_u(nb)%a3d(it,jt,kt) 
                   mb_v(nb)%a3d(is,js,ks) =  2.0*q_1(3) - mb_v(nb)%a3d(it,jt,kt) 
                   mb_w(nb)%a3d(is,js,ks) =  2.0*q_1(4) - mb_w(nb)%a3d(it,jt,kt) 
                   mb_p(nb)%a3d(is,js,ks) =  2.0*q_1(5) - mb_p(nb)%a3d(it,jt,kt)
                   do m = 1,nl-5    !ns
                      mb_fs(nb)%a4d( m,is,js,ks) = q_1(m) 
                   enddo
                   if(mb_r(nb)%a3d(is,js,ks)< rmin_limit .or. &
				      mb_p(nb)%a3d(is,js,ks)< pmin_limit ) then
                      mb_r(nb)%a3d(is,js,ks) =  q_1(1) 
                      mb_u(nb)%a3d(is,js,ks) =  q_1(2) 
                      mb_v(nb)%a3d(is,js,ks) =  q_1(3) 
                      mb_w(nb)%a3d(is,js,ks) =  q_1(4) 
                      mb_p(nb)%a3d(is,js,ks) =  q_1(5) 
                      do m = 1,nl-5    !ns
                         mb_fs(nb)%a4d( m,is,js,ks) = q_1(m) 
                      enddo
				   endif
                   do m=6,nl
                      prim_s1(m) = mb_fs(nb)%a4d(m-5,is,js,ks)
                   enddo

                   if ( nl > 5 ) then
                      zf1 = 0.0
                      do m=6,nl
                         zf1 = zf1 + prim_s1(m)
                      enddo
                      zf1 = min(1.0,zf1)
                      mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                   endif
                enddo
             enddo


          enddo
       enddo

20  continue

    return
end subroutine boundary_7123
!_____________________________________________________________________!
subroutine dif_average(nb,nr,bctype,s_st,s_ed) !有限差分边界值的给定
    use global_variables,only : mb_r,mb_u,mb_v,mb_w,mb_p,mb_bc &
	                          &, r,u,v,w,p
    implicit none
    integer :: s_st(3),s_ed(3),s_lr3d(3),nb,nr,bctype
    integer :: i,j,k,is,js,ks,it,jt,kt,m
!    real :: prim_s1(nl),q_1(nl)
!    real :: zf1

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m)  = mb_bc(nb)%bc(nr)%s_lr3d(m)
    enddo
    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             is = i + s_lr3d(1)
             js = j + s_lr3d(2)
             ks = k + s_lr3d(3)
             it = i - s_lr3d(1)
             jt = j - s_lr3d(2)
             kt = k - s_lr3d(3)  !*TGH. 边界上的值等于两侧值的平均 *!

             mb_r(nb)%a3d(i,j,k) = 0.5*( r(is,js,ks) + r(it,jt,kt) )
             mb_u(nb)%a3d(i,j,k) = 0.5*( u(is,js,ks) + u(it,jt,kt) )
             mb_v(nb)%a3d(i,j,k) = 0.5*( v(is,js,ks) + v(it,jt,kt) )
             mb_w(nb)%a3d(i,j,k) = 0.5*( w(is,js,ks) + w(it,jt,kt) )
             mb_p(nb)%a3d(i,j,k) = 0.5*( p(is,js,ks) + p(it,jt,kt) )
						 
						 call BC_face_dif_turbulence(nb,nb,is,js,ks,it,jt,kt,i,j,k)

          enddo
       enddo
    enddo
    return
end subroutine dif_average
!_____________________________________________________________________!
subroutine boundary_wall(nb)
    use global_variables
    implicit none
    integer :: nb,nr,nrmax,bctype,m,k

    nrmax = mb_bc(nb)%nregions             !本块共有nrmax个边界需要处理
    !具体处理边界条件
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if( bctype < 0 ) then               !对接边界
!          call boundary_n1(nb,nr,bctype)
       elseif( bctype == 20 ) then         !完全气体/完全非催化 等温壁面
          call boundary_wall_20(nb,nr,bctype)
       endif
    enddo
		
    return
end subroutine boundary_wall
!_____________________________________________________________________!
subroutine boundary_wall_20(nb,nr,bctype) !完全气体/完全非催化 等温边界
    use global_variables
    implicit none
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: nb,nr,bctype,m,n,nt
    integer :: is,js,ks,it,jt,kt,i,j,k
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: prim_s1(nl), q_1(nl), zf1
    real :: tw1,pres,dens,mav1,ff,tw2

     if(nvis == 1 ) then !* tgh. 粘性物面
       do m=1,3
          s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
          s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
          s_lr3d(m)= mb_bc(nb)%bc(nr)%s_lr3d(m)
       enddo
       s_nd  = mb_bc(nb)%bc(nr)%s_nd     !边界面方向:1,2,3对应于i,j,k
       s_lr  = mb_bc(nb)%bc(nr)%s_lr     !左右边界-1,1对应于左右边界
       s_fix = mb_bc(nb)%bc(nr)%s_fix    !固定坐标(fixed_coor)
       nbs   = mb_bc(nb)%bc(nr)%nbs      !块号 
       tw2   = 0.50*mb_bc(nb)%bc(nr)%bc_par(1)/tref
       if( method == 0 ) then !有限体积
	       do i = s_st(1),s_ed(1)
              do j = s_st(2),s_ed(2)
                 do k = s_st(3),s_ed(3)
!   第一排给按对称条件延拓
                    is = i + s_lr3d(1)
                    js = j + s_lr3d(2)
                    ks = k + s_lr3d(3)

                    it = i 
                    jt = j 
                    kt = k 
                
                    do m=6,nl
                       prim_s1(m) = mb_fs(nb)%a4d(m-5,it,jt,kt)
                    enddo
                    pres = mb_p(nbs)%a3d(it,jt,kt)
				    tw1  = 4.0*tw2 - mb_t(nbs)%a3d(it,jt,kt)
                    mb_t(nbs)%a3d(is,js,ks) =  tw1

                    is = i + 2*s_lr3d(1)
                    js = j + 2*s_lr3d(2)
                    ks = k + 2*s_lr3d(3)

                    it = i - s_lr3d(1)
                    jt = j - s_lr3d(2)
                    kt = k - s_lr3d(3)
				    tw1  = 4.0*tw2 - mb_t(nbs)%a3d(it,jt,kt)
                    mb_t(nbs)%a3d(is,js,ks) =  tw1
                 enddo
              enddo
           enddo
	   else  !有限差分 给定边界上的值
          do i = s_st(1),s_ed(1)
             do j = s_st(2),s_ed(2)
                do k = s_st(3),s_ed(3)
                   it = i - s_lr3d(1)
                   jt = j - s_lr3d(2)
                   kt = k - s_lr3d(3)

                   is = i 
                   js = j 
                   ks = k 

                   pres = mb_p(nbs)%a3d(it,jt,kt)
				   tw1  = 2.0*tw2 
                      if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                      else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                      endif
                   mb_r(nbs)%a3d(is,js,ks) =  dens
                   mb_u(nbs)%a3d(is,js,ks) =  0.0
                   mb_v(nbs)%a3d(is,js,ks) =  0.0
                   mb_w(nbs)%a3d(is,js,ks) =  0.0
                   mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)

									 call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                   do m=6,nl
                      mb_fs(nb)%a4d(m-5,is,js,ks) = mb_fs(nb)%a4d(m-5,it,jt,kt)  
                   enddo

                   if ( nl > 5 ) then
                      zf1 = 0.0
                      do m=6,nl
                         zf1 = zf1 + mb_fs(nb)%a4d(m-5,is,js,ks)
                      enddo
                      zf1 = min(1.0,zf1)
                      mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                   endif

                   is = i + s_lr3d(1) 
                   js = j + s_lr3d(2) 
                   ks = k + s_lr3d(3) 

				    tw1  = 4.0*tw2 - mb_t(nbs)%a3d(it,jt,kt)
				    if(tw1 < tw2 ) tw1 = tw2

                       if ( nchem_source == 0 ) then
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       else
                          dens = gama*moo*moo*pres/tw1  !注意温度用的是壁面的温度
                       endif

                   mb_r(nbs)%a3d(is,js,ks) =  dens
                   mb_u(nbs)%a3d(is,js,ks) = -mb_u(nbs)%a3d(it,jt,kt)
                   mb_v(nbs)%a3d(is,js,ks) = -mb_v(nbs)%a3d(it,jt,kt)
                   mb_w(nbs)%a3d(is,js,ks) = -mb_w(nbs)%a3d(it,jt,kt)
                   mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbs)%a3d(it,jt,kt)

									 call BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)

                   do m=6,nl
                      mb_fs(nb)%a4d(m-5,is,js,ks) = mb_fs(nb)%a4d(m-5,it,jt,kt)  
                   enddo

                   if ( nl > 5 ) then
                      zf1 = 0.0
                      do m=6,nl
                         zf1 = zf1 + mb_fs(nb)%a4d(m-5,is,js,ks)
                      enddo
                      zf1 = min(1.0,zf1)
                      mb_fs(nb)%a4d(ns,is,js,ks) = 1.0 - zf1
                   endif
!                给密度反号，表示为等温固壁
!                   mb_r(nbs)%a3d(i+2*s_lr3d(1),j+2*s_lr3d(2),k+2*s_lr3d(3) ) = -dens
                enddo
             enddo
          enddo
	   endif
	else  !* tgh. 无粘物面
	   call boundary_2(nb,nr,bctype)
	endif

    return
end subroutine boundary_wall_20

!===========================================================================================!
!===========================================================================================!
subroutine boundary_n1_vir1(nb,nr,bctype) !对接边界条件
!_____________________________________________________________________!
!* 1阶精度 导数法 求解 对接面上的虚网格上的原始变量
!* 通过原始变量在笛卡尔坐标系中的导数来求虚点上的值
!* 湍流mb_qke采用简单对接的方法直接取对于块上的值；
!*     然后call dif_average （含湍流call BC_face_dif_turbulence），
!*     边界上的值等于边界两侧的值的简单平均
!* 特征对接时只求虚点上的值，一般对接时还要求对接面上的值，但是还不适合化学反应的情况
!* 设计：毛枚良
!* 调试：涂国华
!_____________________________________________________________________!
	use define_precision_mod
    use global_variables,only : mb_r,mb_u,mb_v,mb_w,mb_p,mb_bc,method,cic1,mb_vol &
	                          , mb_kcx,mb_kcy,mb_kcz,mb_etx,mb_ety,mb_etz,mb_ctx,mb_cty,mb_ctz &
							  , mb_dim,nl,nmax,mb_x,mb_y,mb_z,dis_tao
    implicit none
    integer :: nb,nr,bctype,m,n
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,i,j,k,nt,ntarg
    integer :: it,jt,kt,ist,jst,kst,it0,jt0,kt0
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real(prec) :: prim_s1(nl),q_1(nl),zf1,ttt1
	real(prec) :: d1i,d1j,d1k,d1x,d1y,d1z
	real(prec) :: volp,nxyz(9),nijk(9),druvwp(15),druvwpdxyz(15)
	real(prec) :: tk1,tk2,tk3
	real(prec) :: rm,um,vm,wm,pm
	real(prec) :: dis1,dis2,dis12 ! 对接两侧网格间距
	integer :: nit,njt,nkt,nis,njs,nks
	integer :: nerr,key_scm

!	dis_tao = 5.0 !间距比超过dis_tao时，修正虚点法求得的差量
    key_scm = 1

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !起始点坐标(可动坐标)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !终止点坐标(可动坐标)
    enddo
    s_nd  = mb_bc(nb)%bc(nr)%s_nd            !边界面方向:1,2,3对应于i,j,k
    s_lr  = mb_bc(nb)%bc(nr)%s_lr            !左右边界-1,1对应于左右边界
    s_fix = mb_bc(nb)%bc(nr)%s_fix           !固定坐标(fixed_coor)
    nbs   = mb_bc(nb)%bc(nr)%nbs             !块号

    nbt   = mb_bc(nb)%bc(nr)%nbt             !对应于块与块边界条件,指出对应目标区域信息在第几个块上
    t_nd  = mb_bc(nb)%bc(nr)%t_nd            !边界面方向:1,2,3对应于i,j,k
    t_lr  = mb_bc(nb)%bc(nr)%t_lr            !左右边界-1,1对应于左右边界
    t_fix = mb_bc(nb)%bc(nr)%t_fix           !固定坐标(fixed_coor)

    nis = mb_dim(nb,1)
    njs = mb_dim(nb,2)
    nks = mb_dim(nb,3)

    nit = mb_dim(nbt,1)
    njt = mb_dim(nbt,2)
    nkt = mb_dim(nbt,3)

	s_lr3d(1) = mb_bc(nb)%bc(nr)%s_lr3d(1)
	s_lr3d(2) = mb_bc(nb)%bc(nr)%s_lr3d(2)
	s_lr3d(3) = mb_bc(nb)%bc(nr)%s_lr3d(3)

	   nerr = 0
       do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
       do k = s_st(3),s_ed(3)


            is = i + s_lr3d(1)
            js = j + s_lr3d(2)
            ks = k + s_lr3d(3)
            ist = i - s_lr3d(1)
            jst = j - s_lr3d(2)
            kst = k - s_lr3d(3)

			it0 = mb_bc(nb)%bc(nr)%image(i,j,k)
			jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,k)
			kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,k)
            it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)
            jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)
            kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)

			dis1 = sqrt( (mb_x(nbs)%a3d(i , j , k ) - mb_x(nbs)%a3d(ist,jst,kst))**2 + &
						 (mb_y(nbs)%a3d(i , j , k ) - mb_y(nbs)%a3d(ist,jst,kst))**2 + &
						 (mb_z(nbs)%a3d(i , j , k ) - mb_z(nbs)%a3d(ist,jst,kst))**2 )

			dis2 = sqrt( (mb_x(nbt)%a3d(it0,jt0,kt0) - mb_x(nbt)%a3d(it ,jt ,kt ))**2 + &
						 (mb_y(nbt)%a3d(it0,jt0,kt0) - mb_y(nbt)%a3d(it ,jt ,kt ))**2 + &
						 (mb_z(nbt)%a3d(it0,jt0,kt0) - mb_z(nbt)%a3d(it ,jt ,kt ))**2 )

			dis12 = dis1/dis2

            !----------------------------------------------------------------------!
            !* 以下开始计算对应块在对接面上原始变量在直角坐标系中的一阶导数        !
			volp    = mb_vol(nbt)%a3d(it0,jt0,kt0)
			nxyz(1) = mb_kcx(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(2) = mb_kcy(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(3) = mb_kcz(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(4) = mb_etx(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(5) = mb_ety(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(6) = mb_etz(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(7) = mb_ctx(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(8) = mb_cty(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(9) = mb_ctz(nbt)%a3d(it0,jt0,kt0) / volp


			if(it0 == 1) then

			  if(key_scm == 1)then
			  druvwp(1) =  mb_r(nbt)%a3d(2,jt0,kt0) &
			              -mb_r(nbt)%a3d(1,jt0,kt0)
			  druvwp(2) =  mb_u(nbt)%a3d(2,jt0,kt0) &
			              -mb_u(nbt)%a3d(1,jt0,kt0)
			  druvwp(3) =  mb_v(nbt)%a3d(2,jt0,kt0) &
			              -mb_v(nbt)%a3d(1,jt0,kt0)
			  druvwp(4) =  mb_w(nbt)%a3d(2,jt0,kt0) &
			              -mb_w(nbt)%a3d(1,jt0,kt0)
			  druvwp(5) =  mb_p(nbt)%a3d(2,jt0,kt0) &
			              -mb_p(nbt)%a3d(1,jt0,kt0)

			  else

			  druvwp(1) = -1.5_prec*mb_r(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_r(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_r(nbt)%a3d(3,jt0,kt0)
			  druvwp(2) = -1.5_prec*mb_u(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_u(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_u(nbt)%a3d(3,jt0,kt0)
			  druvwp(3) = -1.5_prec*mb_v(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_v(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_v(nbt)%a3d(3,jt0,kt0)
			  druvwp(4) = -1.5_prec*mb_w(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_w(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_w(nbt)%a3d(3,jt0,kt0)
			  druvwp(5) = -1.5_prec*mb_p(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_p(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_p(nbt)%a3d(3,jt0,kt0)
			  endif

			elseif(it0 == nit)then

			  if(key_scm == 1)then

			  druvwp(1) =  mb_r(nbt)%a3d(it0  ,jt0,kt0) &
			              -mb_r(nbt)%a3d(it0-1,jt0,kt0)
			  druvwp(2) =  mb_u(nbt)%a3d(it0  ,jt0,kt0) &
			              -mb_u(nbt)%a3d(it0-1,jt0,kt0)
			  druvwp(3) =  mb_v(nbt)%a3d(it0  ,jt0,kt0) &
			              -mb_v(nbt)%a3d(it0-1,jt0,kt0)
			  druvwp(4) =  mb_w(nbt)%a3d(it0  ,jt0,kt0) &
			              -mb_w(nbt)%a3d(it0-1,jt0,kt0)
			  druvwp(5) =  mb_p(nbt)%a3d(it0  ,jt0,kt0) &
			              -mb_p(nbt)%a3d(it0-1,jt0,kt0)

			  else

			  druvwp(1) =  1.5_prec*mb_r(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_r(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_r(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(2) =  1.5_prec*mb_u(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_u(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_u(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(3) =  1.5_prec*mb_v(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_v(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_v(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(4) =  1.5_prec*mb_w(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_w(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_w(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(5) =  1.5_prec*mb_p(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_p(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_p(nbt)%a3d(it0-2,jt0,kt0)
			  endif

			else
			  druvwp(1) =  0.5_prec*( mb_r(nbt)%a3d(it0+1,jt0,kt0)-mb_r(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(2) =  0.5_prec*( mb_u(nbt)%a3d(it0+1,jt0,kt0)-mb_u(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(3) =  0.5_prec*( mb_v(nbt)%a3d(it0+1,jt0,kt0)-mb_v(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(4) =  0.5_prec*( mb_w(nbt)%a3d(it0+1,jt0,kt0)-mb_w(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(5) =  0.5_prec*( mb_p(nbt)%a3d(it0+1,jt0,kt0)-mb_p(nbt)%a3d(it0-1,jt0,kt0) )
			endif

			if(jt0 == 1)then

  			  if(key_scm == 1)then

			  druvwp(6) =  mb_r(nbt)%a3d(it0,2,kt0) &
			              -mb_r(nbt)%a3d(it0,1,kt0)
			  druvwp(7) =  mb_u(nbt)%a3d(it0,2,kt0) &
			              -mb_u(nbt)%a3d(it0,1,kt0)
			  druvwp(8) =  mb_v(nbt)%a3d(it0,2,kt0) &
			              -mb_v(nbt)%a3d(it0,1,kt0)
			  druvwp(9) =  mb_w(nbt)%a3d(it0,2,kt0) &
			              -mb_w(nbt)%a3d(it0,1,kt0)
			  druvwp(10)=  mb_p(nbt)%a3d(it0,2,kt0) &
			              -mb_p(nbt)%a3d(it0,1,kt0)

			  else

			  druvwp(6) = -1.5_prec*mb_r(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_r(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_r(nbt)%a3d(it0,3,kt0)
			  druvwp(7) = -1.5_prec*mb_u(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_u(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_u(nbt)%a3d(it0,3,kt0)
			  druvwp(8) = -1.5_prec*mb_v(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_v(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_v(nbt)%a3d(it0,3,kt0)
			  druvwp(9) = -1.5_prec*mb_w(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_w(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_w(nbt)%a3d(it0,3,kt0)
			  druvwp(10)= -1.5_prec*mb_p(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_p(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_p(nbt)%a3d(it0,3,kt0)

			  endif

			elseif(jt0 == njt)then

			  if(key_scm == 1)then
			  druvwp(6) =  mb_r(nbt)%a3d(it0,jt0  ,kt0) &
			              -mb_r(nbt)%a3d(it0,jt0-1,kt0)
			  druvwp(7) =  mb_u(nbt)%a3d(it0,jt0  ,kt0) &
			              -mb_u(nbt)%a3d(it0,jt0-1,kt0)
			  druvwp(8) =  mb_v(nbt)%a3d(it0,jt0  ,kt0) &
			              -mb_v(nbt)%a3d(it0,jt0-1,kt0)
			  druvwp(9) =  mb_w(nbt)%a3d(it0,jt0  ,kt0) &
			              -mb_w(nbt)%a3d(it0,jt0-1,kt0)
			  druvwp(10)=  mb_p(nbt)%a3d(it0,jt0  ,kt0) &
			              -mb_p(nbt)%a3d(it0,jt0-1,kt0)

			  else

			  druvwp(6) =  1.5_prec*mb_r(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(7) =  1.5_prec*mb_u(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(8) =  1.5_prec*mb_v(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(9) =  1.5_prec*mb_w(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(10)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0-2,kt0)
			  endif

			else
			  druvwp(6) =  0.5_prec*( mb_r(nbt)%a3d(it0,jt0+1,kt0)-mb_r(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(7) =  0.5_prec*( mb_u(nbt)%a3d(it0,jt0+1,kt0)-mb_u(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(8) =  0.5_prec*( mb_v(nbt)%a3d(it0,jt0+1,kt0)-mb_v(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(9) =  0.5_prec*( mb_w(nbt)%a3d(it0,jt0+1,kt0)-mb_w(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(10)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0+1,kt0)-mb_p(nbt)%a3d(it0,jt0-1,kt0) )
			endif

			if(kt0 == 1)then

			  if(key_scm == 1)then

			  druvwp(11)=  mb_r(nbt)%a3d(it0,jt0,2) &
			              -mb_r(nbt)%a3d(it0,jt0,1)
			  druvwp(12)=  mb_u(nbt)%a3d(it0,jt0,2) &
			              -mb_u(nbt)%a3d(it0,jt0,1)
			  druvwp(13)=  mb_v(nbt)%a3d(it0,jt0,2) &
			              -mb_v(nbt)%a3d(it0,jt0,1)
			  druvwp(14)=  mb_w(nbt)%a3d(it0,jt0,2) &
			              -mb_w(nbt)%a3d(it0,jt0,1)
			  druvwp(15)=  mb_p(nbt)%a3d(it0,jt0,2) &
			              -mb_p(nbt)%a3d(it0,jt0,1)

			  else

			  druvwp(11)= -1.5_prec*mb_r(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_r(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_r(nbt)%a3d(it0,jt0,3)
			  druvwp(12)= -1.5_prec*mb_u(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_u(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_u(nbt)%a3d(it0,jt0,3)
			  druvwp(13)= -1.5_prec*mb_v(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_v(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_v(nbt)%a3d(it0,jt0,3)
			  druvwp(14)= -1.5_prec*mb_w(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_w(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_w(nbt)%a3d(it0,jt0,3)
			  druvwp(15)= -1.5_prec*mb_p(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_p(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_p(nbt)%a3d(it0,jt0,3)
			  endif

			elseif(kt0 == nkt)then

			  if(key_scm == 1)then

			  druvwp(11)=  mb_r(nbt)%a3d(it0,jt0,kt0  ) &
			              -mb_r(nbt)%a3d(it0,jt0,kt0-1)
			  druvwp(12)=  mb_u(nbt)%a3d(it0,jt0,kt0  ) &
			              -mb_u(nbt)%a3d(it0,jt0,kt0-1)
			  druvwp(13)=  mb_v(nbt)%a3d(it0,jt0,kt0  ) &
			              -mb_v(nbt)%a3d(it0,jt0,kt0-1)
			  druvwp(14)=  mb_w(nbt)%a3d(it0,jt0,kt0  ) &
			              -mb_w(nbt)%a3d(it0,jt0,kt0-1)
			  druvwp(15)=  mb_p(nbt)%a3d(it0,jt0,kt0  ) &
			              -mb_p(nbt)%a3d(it0,jt0,kt0-1)

			  else

			  druvwp(11)=  1.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(12)=  1.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(13)=  1.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(14)=  1.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(15)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0-2)

			  endif

			else
			  druvwp(11)=  0.5_prec*( mb_r(nbt)%a3d(it0,jt0,kt0+1)-mb_r(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(12)=  0.5_prec*( mb_u(nbt)%a3d(it0,jt0,kt0+1)-mb_u(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(13)=  0.5_prec*( mb_v(nbt)%a3d(it0,jt0,kt0+1)-mb_v(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(14)=  0.5_prec*( mb_w(nbt)%a3d(it0,jt0,kt0+1)-mb_w(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(15)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0,kt0+1)-mb_p(nbt)%a3d(it0,jt0,kt0-1) )
			endif

			call druvwp_xyz(nxyz,druvwp,druvwpdxyz)
                 !* 把计算系下原始变量的一阶导数转换成直角系的一阶导数
                 !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                 !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk
                 !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz

            !------------------------------------------------------------------------------------!
			!* 以下把对接面上原始变量在直角坐标系中的一阶导数化成在本块网格中计算系下的一阶导数  !

			volp    = mb_vol(nb)%a3d(i,j,k)
			nxyz(1) = mb_kcx(nb)%a3d(i,j,k) / volp
			nxyz(2) = mb_kcy(nb)%a3d(i,j,k) / volp
			nxyz(3) = mb_kcz(nb)%a3d(i,j,k) / volp
			nxyz(4) = mb_etx(nb)%a3d(i,j,k) / volp
			nxyz(5) = mb_ety(nb)%a3d(i,j,k) / volp
			nxyz(6) = mb_etz(nb)%a3d(i,j,k) / volp
			nxyz(7) = mb_ctx(nb)%a3d(i,j,k) / volp
			nxyz(8) = mb_cty(nb)%a3d(i,j,k) / volp
			nxyz(9) = mb_ctz(nb)%a3d(i,j,k) / volp

!			call get_nijk(nb,nis,njs,nks,i,j,k,nijk) !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
			call nxyz_nijk(nxyz,nijk,volp) !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                                           !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
			call druvwp_xyz(nijk,druvwpdxyz,druvwp)
                 !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz
			     !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk

           !---------------------------------------------------------------------------------!
            !* 以下根据梯度求 虚点与边界点之间的差量                 ------------------------!
			rm =  s_lr3d(1)*druvwp(1) + s_lr3d(2)*druvwp(6)  + s_lr3d(3)*druvwp(11)
			um =  s_lr3d(1)*druvwp(2) + s_lr3d(2)*druvwp(7)  + s_lr3d(3)*druvwp(12)
			vm =  s_lr3d(1)*druvwp(3) + s_lr3d(2)*druvwp(8)  + s_lr3d(3)*druvwp(13)
			wm =  s_lr3d(1)*druvwp(4) + s_lr3d(2)*druvwp(9)  + s_lr3d(3)*druvwp(14)
			pm =  s_lr3d(1)*druvwp(5) + s_lr3d(2)*druvwp(10) + s_lr3d(3)*druvwp(15)

			if(dis12 > dis_tao)then
				dis12 = dis_tao/dis12

				rm = rm * dis12
				um = um * dis12
				vm = vm * dis12
				wm = wm * dis12
				pm = pm * dis12
			endif

          !* 以下根据 "差量" 赋虚点上的值                               ------------------------!
			rm = mb_r(nb)%a3d(i,j,k) + rm
			um = mb_u(nb)%a3d(i,j,k) + um
			vm = mb_v(nb)%a3d(i,j,k) + vm
			wm = mb_w(nb)%a3d(i,j,k) + wm
			pm = mb_p(nb)%a3d(i,j,k) + pm

901	format(1x,'第',i3,' 块的第',i2,' 边界采用导数法求虚点值失败，密度：',e12.4,3i4)
902	format(1x,'第',i3,' 块的第',i2,' 边界采用导数法求虚点值失败，压力：',e12.4,3i4)

 			if(pm <= 0._prec) then
			    nerr = nerr + 1
!				if(nerr<3)write(*,902)nb,nr,pm,i,j,k

                rm = mb_r(nbt)%a3d(it,jt,kt)
                um = mb_u(nbt)%a3d(it,jt,kt)
                vm = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. 虚层的原始变量 赋 对应块上向内收缩一层后的值
                wm = mb_w(nbt)%a3d(it,jt,kt)
                pm = mb_p(nbt)%a3d(it,jt,kt)

			endif

			if(rm <= 0._prec) then
			    nerr = nerr + 1
!				if(nerr<3)write(*,901)nb,nr,rm,i,j,k

                rm = mb_r(nbt)%a3d(it,jt,kt)
                um = mb_u(nbt)%a3d(it,jt,kt)
                vm = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. 虚层的原始变量 赋 对应块上向内收缩一层后的值
                wm = mb_w(nbt)%a3d(it,jt,kt)
                pm = mb_p(nbt)%a3d(it,jt,kt)

			endif

            mb_r(nb)%a3d(is,js,ks) = rm
            mb_u(nb)%a3d(is,js,ks) = um
            mb_v(nb)%a3d(is,js,ks) = vm
            mb_w(nb)%a3d(is,js,ks) = wm
            mb_p(nb)%a3d(is,js,ks) = pm
           !* 虚点赋值完毕  ------------------------!

			call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC

			do m=1,cic1
			    call BC_face_dif_turbulence(nb,nb,is,js,ks,ist,jst,kst,i,j,k)
			enddo

       enddo
       enddo
       enddo

903	format(1x,'第',i3,' 块的第',i2,'边界共', i4 ' 点采用导数法求虚点值失败')

	   if(nerr > 0)write(*,903)nb,nr,nerr

	   if(cic1 /= 1) call dif_average(nb,nr,bctype,s_st,s_ed) !* TGH. 边界上的值等于两侧的平均，同时还有处理湍流模式边界 *!


    return
end subroutine boundary_n1_vir1

!===========================================================================================!
!===========================================================================================!
!===========================================================================================!

subroutine boundary_n1_vir2(nb,nr,bctype) !对接边界条件
!_____________________________________________________________________!
!* 2阶精度 导数法 求解 对接面上的虚网格上的原始变量
!* 通过原始变量在笛卡尔坐标系中的导数来求虚点上的值
!* 湍流mb_qke采用简单对接的方法直接取对于块上的值；
!*     然后call dif_average （含湍流call BC_face_dif_turbulence），
!*     边界上的值等于边界两侧的值的简单平均
!* 特征对接时只求虚点上的值，一般对接时还要求对接面上的值，但是还不适合化学反应的情况
!* 设计：毛枚良
!* 调试：涂国华
!_____________________________________________________________________!
	use define_precision_mod
    use global_variables,only : mb_r,mb_u,mb_v,mb_w,mb_p,mb_bc,method,cic1,mb_vol &
	                          , mb_kcx,mb_kcy,mb_kcz,mb_etx,mb_ety,mb_etz,mb_ctx,mb_cty,mb_ctz &
							  , mb_dim,nl,nmax,mb_x,mb_y,mb_z,dis_tao
    implicit none
    integer :: nb,nr,bctype,m,n
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,i,j,k,nt,ntarg
    integer :: it,jt,kt,ist,jst,kst,it0,jt0,kt0
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real(prec) :: prim_s1(nl),q_1(nl),zf1,ttt1
	real(prec) :: d1i,d1j,d1k,d1x,d1y,d1z
	real(prec) :: volp,nxyz(9),nijk(9),druvwp(15),druvwpdxyz(15)
	real(prec) :: tk1,tk2,tk3
	real(prec) :: rm,um,vm,wm,pm
	real(prec) :: dis1,dis2,dis12  ! 对接两侧网格间距
	integer :: nit,njt,nkt,nis,njs,nks
	integer :: nerr

!	dis_tao = 2.0 !间距比超过dis_tao时，修正虚点法求得的差量
    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !起始点坐标(可动坐标)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !终止点坐标(可动坐标)
    enddo
    s_nd  = mb_bc(nb)%bc(nr)%s_nd            !边界面方向:1,2,3对应于i,j,k
    s_lr  = mb_bc(nb)%bc(nr)%s_lr            !左右边界-1,1对应于左右边界
    s_fix = mb_bc(nb)%bc(nr)%s_fix           !固定坐标(fixed_coor)
    nbs   = mb_bc(nb)%bc(nr)%nbs             !块号

    nbt   = mb_bc(nb)%bc(nr)%nbt             !对应于块与块边界条件,指出对应目标区域信息在第几个块上
    t_nd  = mb_bc(nb)%bc(nr)%t_nd            !边界面方向:1,2,3对应于i,j,k
    t_lr  = mb_bc(nb)%bc(nr)%t_lr            !左右边界-1,1对应于左右边界
    t_fix = mb_bc(nb)%bc(nr)%t_fix           !固定坐标(fixed_coor)

    nis = mb_dim(nb,1)
    njs = mb_dim(nb,2)
    nks = mb_dim(nb,3)

    nit = mb_dim(nbt,1)
    njt = mb_dim(nbt,2)
    nkt = mb_dim(nbt,3)

	s_lr3d(1) = mb_bc(nb)%bc(nr)%s_lr3d(1)
	s_lr3d(2) = mb_bc(nb)%bc(nr)%s_lr3d(2)
	s_lr3d(3) = mb_bc(nb)%bc(nr)%s_lr3d(3)

	   nerr = 0
       do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
       do k = s_st(3),s_ed(3)


            is = i + s_lr3d(1)
            js = j + s_lr3d(2)
            ks = k + s_lr3d(3)
            ist = i - s_lr3d(1)
            jst = j - s_lr3d(2)
            kst = k - s_lr3d(3)

			it0 = mb_bc(nb)%bc(nr)%image(i,j,k)
			jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,k)
			kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,k)
            it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)
            jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)
            kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)

			dis1 = sqrt( (mb_x(nbs)%a3d(i , j , k ) - mb_x(nbs)%a3d(ist,jst,kst))**2 + &
						 (mb_y(nbs)%a3d(i , j , k ) - mb_y(nbs)%a3d(ist,jst,kst))**2 + &
						 (mb_z(nbs)%a3d(i , j , k ) - mb_z(nbs)%a3d(ist,jst,kst))**2 )

			dis2 = sqrt( (mb_x(nbt)%a3d(it0,jt0,kt0) - mb_x(nbt)%a3d(it ,jt ,kt ))**2 + &
						 (mb_y(nbt)%a3d(it0,jt0,kt0) - mb_y(nbt)%a3d(it ,jt ,kt ))**2 + &
						 (mb_z(nbt)%a3d(it0,jt0,kt0) - mb_z(nbt)%a3d(it ,jt ,kt ))**2 )

			dis12 = dis1/dis2

            !----------------------------------------------------------------------!
            !* 以下开始计算对应块在对接面上原始变量在直角坐标系中的一阶导数        !
			volp    = mb_vol(nbt)%a3d(it0,jt0,kt0)
			nxyz(1) = mb_kcx(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(2) = mb_kcy(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(3) = mb_kcz(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(4) = mb_etx(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(5) = mb_ety(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(6) = mb_etz(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(7) = mb_ctx(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(8) = mb_cty(nbt)%a3d(it0,jt0,kt0) / volp
			nxyz(9) = mb_ctz(nbt)%a3d(it0,jt0,kt0) / volp


			if(it0 == 1) then
			  druvwp(1) = -1.5_prec*mb_r(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_r(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_r(nbt)%a3d(3,jt0,kt0)
			  druvwp(2) = -1.5_prec*mb_u(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_u(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_u(nbt)%a3d(3,jt0,kt0)
			  druvwp(3) = -1.5_prec*mb_v(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_v(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_v(nbt)%a3d(3,jt0,kt0)
			  druvwp(4) = -1.5_prec*mb_w(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_w(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_w(nbt)%a3d(3,jt0,kt0)
			  druvwp(5) = -1.5_prec*mb_p(nbt)%a3d(1,jt0,kt0) &
			              +2.0_prec*mb_p(nbt)%a3d(2,jt0,kt0) &
			              -0.5_prec*mb_p(nbt)%a3d(3,jt0,kt0)
			elseif(it0 == nit)then
			  druvwp(1) =  1.5_prec*mb_r(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_r(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_r(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(2) =  1.5_prec*mb_u(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_u(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_u(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(3) =  1.5_prec*mb_v(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_v(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_v(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(4) =  1.5_prec*mb_w(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_w(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_w(nbt)%a3d(it0-2,jt0,kt0)
			  druvwp(5) =  1.5_prec*mb_p(nbt)%a3d(it0  ,jt0,kt0) &
			              -2.0_prec*mb_p(nbt)%a3d(it0-1,jt0,kt0) &
			              +0.5_prec*mb_p(nbt)%a3d(it0-2,jt0,kt0)
			else
			  druvwp(1) =  0.5_prec*( mb_r(nbt)%a3d(it0+1,jt0,kt0)-mb_r(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(2) =  0.5_prec*( mb_u(nbt)%a3d(it0+1,jt0,kt0)-mb_u(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(3) =  0.5_prec*( mb_v(nbt)%a3d(it0+1,jt0,kt0)-mb_v(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(4) =  0.5_prec*( mb_w(nbt)%a3d(it0+1,jt0,kt0)-mb_w(nbt)%a3d(it0-1,jt0,kt0) )
			  druvwp(5) =  0.5_prec*( mb_p(nbt)%a3d(it0+1,jt0,kt0)-mb_p(nbt)%a3d(it0-1,jt0,kt0) )
			endif

			if(jt0 == 1)then
			  druvwp(6) = -1.5_prec*mb_r(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_r(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_r(nbt)%a3d(it0,3,kt0)
			  druvwp(7) = -1.5_prec*mb_u(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_u(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_u(nbt)%a3d(it0,3,kt0)
			  druvwp(8) = -1.5_prec*mb_v(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_v(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_v(nbt)%a3d(it0,3,kt0)
			  druvwp(9) = -1.5_prec*mb_w(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_w(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_w(nbt)%a3d(it0,3,kt0)
			  druvwp(10)= -1.5_prec*mb_p(nbt)%a3d(it0,1,kt0) &
			              +2.0_prec*mb_p(nbt)%a3d(it0,2,kt0) &
			              -0.5_prec*mb_p(nbt)%a3d(it0,3,kt0)
			elseif(jt0 == njt)then
			  druvwp(6) =  1.5_prec*mb_r(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(7) =  1.5_prec*mb_u(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(8) =  1.5_prec*mb_v(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(9) =  1.5_prec*mb_w(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0-2,kt0)
			  druvwp(10)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0  ,kt0) &
			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0-1,kt0) &
			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0-2,kt0)
			else
			  druvwp(6) =  0.5_prec*( mb_r(nbt)%a3d(it0,jt0+1,kt0)-mb_r(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(7) =  0.5_prec*( mb_u(nbt)%a3d(it0,jt0+1,kt0)-mb_u(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(8) =  0.5_prec*( mb_v(nbt)%a3d(it0,jt0+1,kt0)-mb_v(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(9) =  0.5_prec*( mb_w(nbt)%a3d(it0,jt0+1,kt0)-mb_w(nbt)%a3d(it0,jt0-1,kt0) )
			  druvwp(10)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0+1,kt0)-mb_p(nbt)%a3d(it0,jt0-1,kt0) )
			endif

			if(kt0 == 1)then
			  druvwp(11)= -1.5_prec*mb_r(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_r(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_r(nbt)%a3d(it0,jt0,3)
			  druvwp(12)= -1.5_prec*mb_u(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_u(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_u(nbt)%a3d(it0,jt0,3)
			  druvwp(13)= -1.5_prec*mb_v(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_v(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_v(nbt)%a3d(it0,jt0,3)
			  druvwp(14)= -1.5_prec*mb_w(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_w(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_w(nbt)%a3d(it0,jt0,3)
			  druvwp(15)= -1.5_prec*mb_p(nbt)%a3d(it0,jt0,1) &
			              +2.0_prec*mb_p(nbt)%a3d(it0,jt0,2) &
			              -0.5_prec*mb_p(nbt)%a3d(it0,jt0,3)
			elseif(kt0 == nkt)then
			  druvwp(11)=  1.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(12)=  1.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(13)=  1.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(14)=  1.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0-2)
			  druvwp(15)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0  ) &
			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0,kt0-1) &
			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0-2)
			else
			  druvwp(11)=  0.5_prec*( mb_r(nbt)%a3d(it0,jt0,kt0+1)-mb_r(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(12)=  0.5_prec*( mb_u(nbt)%a3d(it0,jt0,kt0+1)-mb_u(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(13)=  0.5_prec*( mb_v(nbt)%a3d(it0,jt0,kt0+1)-mb_v(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(14)=  0.5_prec*( mb_w(nbt)%a3d(it0,jt0,kt0+1)-mb_w(nbt)%a3d(it0,jt0,kt0-1) )
			  druvwp(15)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0,kt0+1)-mb_p(nbt)%a3d(it0,jt0,kt0-1) )
			endif

			call druvwp_xyz(nxyz,druvwp,druvwpdxyz)
                 !* 把计算系下原始变量的一阶导数转换成直角系的一阶导数
                 !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                 !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk
                 !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz

            !------------------------------------------------------------------------------------!
			!* 以下把对接面上原始变量在直角坐标系中的一阶导数化成在本块网格中计算系下的一阶导数  !

			volp    = mb_vol(nb)%a3d(i,j,k)
			nxyz(1) = mb_kcx(nb)%a3d(i,j,k) / volp
			nxyz(2) = mb_kcy(nb)%a3d(i,j,k) / volp
			nxyz(3) = mb_kcz(nb)%a3d(i,j,k) / volp
			nxyz(4) = mb_etx(nb)%a3d(i,j,k) / volp
			nxyz(5) = mb_ety(nb)%a3d(i,j,k) / volp
			nxyz(6) = mb_etz(nb)%a3d(i,j,k) / volp
			nxyz(7) = mb_ctx(nb)%a3d(i,j,k) / volp
			nxyz(8) = mb_cty(nb)%a3d(i,j,k) / volp
			nxyz(9) = mb_ctz(nb)%a3d(i,j,k) / volp

!			call get_nijk(nb,nis,njs,nks,i,j,k,nijk) !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
			call nxyz_nijk(nxyz,nijk,volp) !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                                           !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
			call druvwp_xyz(nijk,druvwpdxyz,druvwp)
                 !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz
			     !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk

           !---------------------------------------------------------------------------------!
            !* 以下根据梯度求 虚点与边界点之间的差量                 ------------------------!
			rm =  s_lr3d(1)*druvwp(1) + s_lr3d(2)*druvwp(6)  + s_lr3d(3)*druvwp(11)
			um =  s_lr3d(1)*druvwp(2) + s_lr3d(2)*druvwp(7)  + s_lr3d(3)*druvwp(12)
			vm =  s_lr3d(1)*druvwp(3) + s_lr3d(2)*druvwp(8)  + s_lr3d(3)*druvwp(13)
			wm =  s_lr3d(1)*druvwp(4) + s_lr3d(2)*druvwp(9)  + s_lr3d(3)*druvwp(14)
			pm =  s_lr3d(1)*druvwp(5) + s_lr3d(2)*druvwp(10) + s_lr3d(3)*druvwp(15)

			if(dis12 > dis_tao)then
				dis12 = dis_tao/dis12

				rm = rm * dis12
				um = um * dis12
				vm = vm * dis12
				wm = wm * dis12
				pm = pm * dis12
			endif

          !* 以下根据 "差量" 赋虚点上的值                               ------------------------!
			rm = mb_r(nb)%a3d(ist,jst,kst) + rm*2.0
			um = mb_u(nb)%a3d(ist,jst,kst) + um*2.0
			vm = mb_v(nb)%a3d(ist,jst,kst) + vm*2.0
			wm = mb_w(nb)%a3d(ist,jst,kst) + wm*2.0
			pm = mb_p(nb)%a3d(ist,jst,kst) + pm*2.0

901	format(1x,'第',i3,' 块的第',i2,' 边界采用导数法求虚点值失败，密度：',e12.4,3i4)
902	format(1x,'第',i3,' 块的第',i2,' 边界采用导数法求虚点值失败，压力：',e12.4,3i4)
903	format(1x,'第',i3,' 块的第',i2,'边界共', i4 ' 点采用导数法求虚点值失败')

 			if(pm <= 0._prec) then
			    nerr = nerr + 1
!				if(nerr<3)write(*,902)nb,nr,pm,i,j,k

                rm = mb_r(nbt)%a3d(it,jt,kt)
                um = mb_u(nbt)%a3d(it,jt,kt)
                vm = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. 虚层的原始变量 赋 对应块上向内收缩一层后的值
                wm = mb_w(nbt)%a3d(it,jt,kt)
                pm = mb_p(nbt)%a3d(it,jt,kt)

			endif

			if(rm <= 0._prec) then
			    nerr = nerr + 1
!				if(nerr<3)write(*,901)nb,nr,rm,i,j,k

                rm = mb_r(nbt)%a3d(it,jt,kt)
                um = mb_u(nbt)%a3d(it,jt,kt)
                vm = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. 虚层的原始变量 赋 对应块上向内收缩一层后的值
                wm = mb_w(nbt)%a3d(it,jt,kt)
                pm = mb_p(nbt)%a3d(it,jt,kt)

			endif

            mb_r(nb)%a3d(is,js,ks) = rm
            mb_u(nb)%a3d(is,js,ks) = um
            mb_v(nb)%a3d(is,js,ks) = vm
            mb_w(nb)%a3d(is,js,ks) = wm
            mb_p(nb)%a3d(is,js,ks) = pm
           !* 虚点赋值完毕  ------------------------!

			call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC

			do m=1,cic1
			    call BC_face_dif_turbulence(nb,nb,is,js,ks,ist,jst,kst,i,j,k)
			enddo

       enddo
       enddo
       enddo

	   if(nerr > 0)write(*,903)nb,nr,nerr

	   if(cic1 /= 1) call dif_average(nb,nr,bctype,s_st,s_ed) !* TGH. 边界上的值等于两侧的平均，同时还有处理湍流模式边界 *!


    return
end subroutine boundary_n1_vir2

!===========================================================================================!
!===========================================================================================!

!===========================================================================================!

subroutine druvwp_xyz(nxyz,druvwp,druvwpdxyz)
!_____________________________________________________________________!
!* 把计算系下原始变量的一阶导数转换成直角系的一阶导数
!*    或反变换
!* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
!* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk
!* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz
!_____________________________________________________________________!
	use define_precision_mod
	implicit none

	real(prec) :: nxyz(9),druvwp(15),druvwpdxyz(15)
	integer :: n,n1,n2,n3

	do n=1,5
		n1 = n
		n2 = 5+n
		n3 = 10+n

		druvwpdxyz(n1) = nxyz(1)*druvwp(n1)+nxyz(4)*druvwp(n2) + nxyz(7)*druvwp(n3)

		druvwpdxyz(n2) = nxyz(2)*druvwp(n1)+nxyz(5)*druvwp(n2) + nxyz(8)*druvwp(n3)

		druvwpdxyz(n3) = nxyz(3)*druvwp(n1)+nxyz(6)*druvwp(n2) + nxyz(9)*druvwp(n3)

	enddo

end subroutine	druvwp_xyz


!===========================================================================================!
!===========================================================================================!

subroutine nxyz_nijk(nxyz,nijk,vol)
!* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
!* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
    use define_precision_mod
	implicit none
	real :: nxyz(9),nijk(9),vol
	integer :: i

	nijk(1) = nxyz(5)*nxyz(9) - nxyz(8)*nxyz(6)
	nijk(2) = nxyz(8)*nxyz(3) - nxyz(2)*nxyz(9)
	nijk(3) = nxyz(2)*nxyz(6) - nxyz(5)*nxyz(3)

	nijk(4) = nxyz(6)*nxyz(7) - nxyz(9)*nxyz(4)
	nijk(5) = nxyz(9)*nxyz(1) - nxyz(3)*nxyz(7)
	nijk(6) = nxyz(3)*nxyz(4) - nxyz(6)*nxyz(1)

	nijk(7) = nxyz(4)*nxyz(8) - nxyz(7)*nxyz(5)
	nijk(8) = nxyz(7)*nxyz(2) - nxyz(1)*nxyz(8)
	nijk(9) = nxyz(1)*nxyz(5) - nxyz(4)*nxyz(2)

	do i=1,9
		nijk(i) = nijk(i)*vol
	enddo

end subroutine nxyz_nijk


!===========================================================================================!
!===========================================================================================!
subroutine get_nijk(nbt,nit,njt,nkt,it0,jt0,kt0,nijk)
!* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
	use define_precision_mod
    use global_variables,only : mb_x,mb_y,mb_z
	integer :: nbt,nit,njt,nkt,it0,jt0,kt0
	real(prec) :: nijk(9)

	if(it0 == 1) then
	    nijk(1) = -1.5_prec*mb_x(nbt)%a3d(1,jt0,kt0) &
			      +2.0_prec*mb_x(nbt)%a3d(2,jt0,kt0) &
			      -0.5_prec*mb_x(nbt)%a3d(3,jt0,kt0)
	    nijk(4) = -1.5_prec*mb_y(nbt)%a3d(1,jt0,kt0) &
			      +2.0_prec*mb_y(nbt)%a3d(2,jt0,kt0) &
			      -0.5_prec*mb_y(nbt)%a3d(3,jt0,kt0)
		nijk(7) = -1.5_prec*mb_z(nbt)%a3d(1,jt0,kt0) &
		          +2.0_prec*mb_z(nbt)%a3d(2,jt0,kt0) &
			      -0.5_prec*mb_z(nbt)%a3d(3,jt0,kt0)
	elseif(it0 == nit)then
		nijk(1) =  1.5_prec*mb_x(nbt)%a3d(it0  ,jt0,kt0) &
		          -2.0_prec*mb_x(nbt)%a3d(it0-1,jt0,kt0) &
		          +0.5_prec*mb_x(nbt)%a3d(it0-2,jt0,kt0)
		nijk(4) =  1.5_prec*mb_y(nbt)%a3d(it0  ,jt0,kt0) &
		          -2.0_prec*mb_y(nbt)%a3d(it0-1,jt0,kt0) &
		          +0.5_prec*mb_y(nbt)%a3d(it0-2,jt0,kt0)
		nijk(7) =  1.5_prec*mb_z(nbt)%a3d(it0  ,jt0,kt0) &
		          -2.0_prec*mb_z(nbt)%a3d(it0-1,jt0,kt0) &
		          +0.5_prec*mb_z(nbt)%a3d(it0-2,jt0,kt0)
	else
	    nijk(1) =  0.5_prec*( mb_x(nbt)%a3d(it0+1,jt0,kt0)-mb_x(nbt)%a3d(it0-1,jt0,kt0) )
	    nijk(4) =  0.5_prec*( mb_y(nbt)%a3d(it0+1,jt0,kt0)-mb_y(nbt)%a3d(it0-1,jt0,kt0) )
	    nijk(7) =  0.5_prec*( mb_z(nbt)%a3d(it0+1,jt0,kt0)-mb_z(nbt)%a3d(it0-1,jt0,kt0) )
    endif

	if(jt0 == 1)then
	    nijk(2) = -1.5_prec*mb_x(nbt)%a3d(it0,1,kt0) &
		          +2.0_prec*mb_x(nbt)%a3d(it0,2,kt0) &
		          -0.5_prec*mb_x(nbt)%a3d(it0,3,kt0)
		nijk(5) = -1.5_prec*mb_y(nbt)%a3d(it0,1,kt0) &
		          +2.0_prec*mb_y(nbt)%a3d(it0,2,kt0) &
		          -0.5_prec*mb_y(nbt)%a3d(it0,3,kt0)
		nijk(8) = -1.5_prec*mb_z(nbt)%a3d(it0,1,kt0) &
		          +2.0_prec*mb_z(nbt)%a3d(it0,2,kt0) &
		          -0.5_prec*mb_z(nbt)%a3d(it0,3,kt0)
	elseif(jt0 == njt)then
		nijk(2) =  1.5_prec*mb_x(nbt)%a3d(it0,jt0  ,kt0) &
		          -2.0_prec*mb_x(nbt)%a3d(it0,jt0-1,kt0) &
		          +0.5_prec*mb_x(nbt)%a3d(it0,jt0-2,kt0)
		nijk(5) =  1.5_prec*mb_y(nbt)%a3d(it0,jt0  ,kt0) &
		          -2.0_prec*mb_y(nbt)%a3d(it0,jt0-1,kt0) &
		          +0.5_prec*mb_y(nbt)%a3d(it0,jt0-2,kt0)
		nijk(8) =  1.5_prec*mb_z(nbt)%a3d(it0,jt0  ,kt0) &
		          -2.0_prec*mb_z(nbt)%a3d(it0,jt0-1,kt0) &
		          +0.5_prec*mb_z(nbt)%a3d(it0,jt0-2,kt0)

	else
	    nijk(2) =  0.5_prec*( mb_x(nbt)%a3d(it0,jt0+1,kt0)-mb_x(nbt)%a3d(it0,jt0-1,kt0) )
	    nijk(5) =  0.5_prec*( mb_y(nbt)%a3d(it0,jt0+1,kt0)-mb_y(nbt)%a3d(it0,jt0-1,kt0) )
	    nijk(8) =  0.5_prec*( mb_z(nbt)%a3d(it0,jt0+1,kt0)-mb_z(nbt)%a3d(it0,jt0-1,kt0) )
    endif

	if(kt0 == 1)then
		nijk(3) = -1.5_prec*mb_x(nbt)%a3d(it0,jt0,1) &
		          +2.0_prec*mb_x(nbt)%a3d(it0,jt0,2) &
		          -0.5_prec*mb_x(nbt)%a3d(it0,jt0,3)
		nijk(6) = -1.5_prec*mb_y(nbt)%a3d(it0,jt0,1) &
		          +2.0_prec*mb_y(nbt)%a3d(it0,jt0,2) &
		          -0.5_prec*mb_y(nbt)%a3d(it0,jt0,3)
		nijk(9) = -1.5_prec*mb_z(nbt)%a3d(it0,jt0,1) &
		          +2.0_prec*mb_z(nbt)%a3d(it0,jt0,2) &
		          -0.5_prec*mb_z(nbt)%a3d(it0,jt0,3)
	elseif(kt0 == nkt)then
		nijk(3) =  1.5_prec*mb_x(nbt)%a3d(it0,jt0,kt0  ) &
		          -2.0_prec*mb_x(nbt)%a3d(it0,jt0,kt0-1) &
		          +0.5_prec*mb_x(nbt)%a3d(it0,jt0,kt0-2)
		nijk(6) =  1.5_prec*mb_y(nbt)%a3d(it0,jt0,kt0  ) &
		          -2.0_prec*mb_y(nbt)%a3d(it0,jt0,kt0-1) &
		          +0.5_prec*mb_y(nbt)%a3d(it0,jt0,kt0-2)
		nijk(8) =  1.5_prec*mb_z(nbt)%a3d(it0,jt0,kt0  ) &
		          -2.0_prec*mb_z(nbt)%a3d(it0,jt0,kt0-1) &
		          +0.5_prec*mb_z(nbt)%a3d(it0,jt0,kt0-2)
	else
		nijk(3) =  0.5_prec*( mb_x(nbt)%a3d(it0,jt0,kt0+1)-mb_x(nbt)%a3d(it0,jt0,kt0-1) )
		nijk(6) =  0.5_prec*( mb_y(nbt)%a3d(it0,jt0,kt0+1)-mb_y(nbt)%a3d(it0,jt0,kt0-1) )
		nijk(9) =  0.5_prec*( mb_z(nbt)%a3d(it0,jt0,kt0+1)-mb_z(nbt)%a3d(it0,jt0,kt0-1) )
	endif

end subroutine get_nijk

!===========================================================================================!
!===========================================================================================!
!===========================================================================================!
