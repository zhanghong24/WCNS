subroutine WCNSE5_VIS_FLUX_new
  use global_variables,only : u,v,w,t,reynolds,visl,vist,dq,nvis     &
                            , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol &
							, nl,nmax,ni,nj,nk,gama,moo,prl,prt
  use duvwt_all_field
  implicit none
!-----------------------------------------------------------------------------!
!	在原 WCNSE5_VIS_FLUX 基础上做了部分程序结构调整                           !
!	改进地方：                                                                !
!   （1） 把求剪切应力的的子程序挪到求粘性通量的模块中，减少了一次数据传递    !
!   （2） 程序中考虑到了层流和湍流情况，计算速度更快                          !
!   （3） 在层流时，原程序必须对湍流涡粘系数是给0初值，修改后不用             !
!   （4） 原程序把“节点上计算系下的一阶导数插值到半结点”时有1/3的重复计算   !
!         修改后消除了重复计算                                                !
!   （5） 本程序比WCNSE5_VIS_FLUX节约计算时间18% ，同时还更节约内存           !
!                                                                             !
!		设计：毛枚良，涂国华                                                  !
!		调试：涂国华 2009.03.14                                               !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
!	real,pointer,dimension(:,:,:,:) :: duvwt_mid,duvwt
	real    :: duvwtdxyz (12,0:nmax),nxyz_line(3,0:nmax)
	real    :: duvwt_line(12,nmax),uvwt_line(4,nmax)
	real    :: duvwt_half(12,0:nmax),uvwt_half(4,0:nmax)
	real    :: kxyz_line ( 9,nmax),kxyz_half(9,0:nmax)
	real    :: vol_line(nmax),fv(nl,0:nmax),dfv(nl,nmax)
	real    :: vol_half(0:nmax) !,txyz_half(9,0:nmax)
	real    :: vslt1_half(0:nmax),vslt2_half(0:nmax),cp,cp_prl,cp_prt
	real    :: vslt1_line(nmax),vslt2_line(nmax),re
	integer :: mvist

	allocate(duvwt(12,ni,nj,nk))   ! overcome stack over problem
	allocate(duvwt_mid(12,0:ni,0:nj,0:nk))   !___ 直接计算半结点1阶导数  2009.2.1

	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)
    cp_prl = cp/prl
    cp_prt = cp/prt

! to calculate the deriative of u,v,w,t in computing	coordinate

	call UVWT_DER_4th_ORDER (12) !*TGH. 直接求节点1阶导数，暂时没有用虚点值
    call UVWT_DER_4th_ORDER_half(12) !*TGH. 直接求半结点1阶导数
         !*TGH. 修改前要用虚点上的值，修改后不用虚点上的值
!**************** I direction ****************
	do k=1,nk
	do j=1,nj

			do i=1,ni

				uvwt_line(1,i) = u(i,j,k)
				uvwt_line(2,i) = v(i,j,k)
				uvwt_line(3,i) = w(i,j,k)
				uvwt_line(4,i) = t(i,j,k)

				vslt1_line (i) = visl(i,j,k)
				vslt2_line (i) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line (i) = vslt1_line(i) + vist(i,j,k)
				  vslt2_line (i) = vslt2_line(i) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,i) = kcx(i,j,k)
				kxyz_line(2,i) = etx(i,j,k)
				kxyz_line(3,i) = ctx(i,j,k)

				kxyz_line(4,i) = kcy(i,j,k)
				kxyz_line(5,i) = ety(i,j,k)
				kxyz_line(6,i) = cty(i,j,k)

				kxyz_line(7,i) = kcz(i,j,k)
				kxyz_line(8,i) = etz(i,j,k)
				kxyz_line(9,i) = ctz(i,j,k)

				vol_line(   i) = vol(i,j,k)

!				do m=1,12
				do m=5,12
					duvwt_line(m,i) = duvwt(m,i,j,k)
				enddo

			enddo

			call VALUE_HALF_NODE(4 ,nmax,ni,uvwt_line ,uvwt_half ) !求原始变量在半结点上的值
			call VALUE_HALF_NODE_IJK(0,1,1,12,nmax,ni ,duvwt_line,duvwt_half) !把节点上的一阶导数插值到半结点
			call VALUE_HALF_NODE(9 ,nmax,ni,kxyz_line ,kxyz_half )   !把网格导数插值到半结点
			call VALUE_HALF_NODE(1 ,nmax,ni,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,ni,vslt1_line,vslt1_half  ) !*坐标：1:ni --> 0:ni
			call VALUE_HALF_NODE(1 ,nmax,ni,vslt2_line,vslt2_half  )

			do i=0,ni
				do m=1,4
					duvwt_half(m,i) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo


			call DUVWT_DXYZ(nmax,ni,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz) 
!*TGH. 计算空间半结点上的一阶导数换算成物理空间的导数, 传递数值范围0:nmax
			
			do i=0,ni
					nxyz_line(1,i) = kxyz_half(1,i)   !半结点上稳定在物理空间中的一阶导数
					nxyz_line(2,i) = kxyz_half(4,i)
					nxyz_line(3,i) = kxyz_half(7,i)
			enddo


			call FLUX_VIS_LINE_new(nl,nmax,ni,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    )
!			                    12,duvwtdxyz,vslt1_half,vslt2_half,prl,prt,cp,fv    )
		  !*TGH. 按线计算粘性通量在计算空间半结点上的通量

			call FLUX_DXYZ(nl,nmax,ni,fv,dfv)

			do i = 1,ni
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
				enddo

			enddo

	enddo
	enddo
!
!**************** J direction ****************
	do k=1,nk
	do i=1,ni
			do j=1,nj

				uvwt_line(1,j) = u(i,j,k)
				uvwt_line(2,j) = v(i,j,k)
				uvwt_line(3,j) = w(i,j,k)
				uvwt_line(4,j) = t(i,j,k)

				vslt1_line (j) = visl(i,j,k)
				vslt2_line (j) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line (j) = vslt1_line(j) + vist(i,j,k)
				  vslt2_line (j) = vslt2_line(j) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,j) = kcx(i,j,k)
				kxyz_line(2,j) = etx(i,j,k)
				kxyz_line(3,j) = ctx(i,j,k)

				kxyz_line(4,j) = kcy(i,j,k)
				kxyz_line(5,j) = ety(i,j,k)
				kxyz_line(6,j) = cty(i,j,k)

				kxyz_line(7,j) = kcz(i,j,k)
				kxyz_line(8,j) = etz(i,j,k)
				kxyz_line(9,j) = ctz(i,j,k)

				vol_line(   j) = vol(i,j,k)

!				do m=1,12
!					duvwt_line(m,j) = duvwt(m,i,j,k)
!				enddo
				do m=1,4
					duvwt_line(m,j) = duvwt(m,i,j,k)
				enddo
				do m=9,12
					duvwt_line(m,j) = duvwt(m,i,j,k)
				enddo

			enddo

			call VALUE_HALF_NODE(4 ,nmax,nj,uvwt_line ,uvwt_half )
			call VALUE_HALF_NODE_IJK(1,0,1,12,nmax,nj,duvwt_line,duvwt_half)
			call VALUE_HALF_NODE(9 ,nmax,nj,kxyz_line ,kxyz_half )
			call VALUE_HALF_NODE(1 ,nmax,nj,vol_line  ,vol_half  )

			call VALUE_HALF_NODE(1 ,nmax,nj,vslt1_line,vslt1_half  ) !*坐标：1:nj --> 0:nj
			call VALUE_HALF_NODE(1 ,nmax,nj,vslt2_line,vslt2_half  )

			do j=0,nj
				do m=5,8
					duvwt_half(m,j) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo

			call DUVWT_DXYZ(nmax,nj,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)
			
			do j=0,nj
					nxyz_line(1,j) = kxyz_half(2,j)
					nxyz_line(2,j) = kxyz_half(5,j)
					nxyz_line(3,j) = kxyz_half(8,j)
			enddo 

			call  FLUX_VIS_LINE_new(nl,nmax,nj,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    )

			call FLUX_DXYZ(nl,nmax,nj,fv,dfv)

			do j = 1,nj
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j)
				enddo
			enddo

	enddo
	enddo
! 
!**************** K direction ****************
	do j=1,nj 
	do i=1,ni 
			do k=1,nk

				uvwt_line(1,k) = u(i,j,k)
				uvwt_line(2,k) = v(i,j,k)
				uvwt_line(3,k) = w(i,j,k)
				uvwt_line(4,k) = t(i,j,k)

				vslt1_line (k) = visl(i,j,k)
				vslt2_line (k) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line (k) = vslt1_line(k) + vist(i,j,k)
				  vslt2_line (k) = vslt2_line(k) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,k) = kcx(i,j,k)
				kxyz_line(2,k) = etx(i,j,k)
				kxyz_line(3,k) = ctx(i,j,k)

				kxyz_line(4,k) = kcy(i,j,k)
				kxyz_line(5,k) = ety(i,j,k)
				kxyz_line(6,k) = cty(i,j,k)

				kxyz_line(7,k) = kcz(i,j,k)
				kxyz_line(8,k) = etz(i,j,k)
				kxyz_line(9,k) = ctz(i,j,k)

				vol_line(   k) = vol(i,j,k)

!				do m=1,12
				do m=1,8
					duvwt_line(m,k) = duvwt(m,i,j,k)
				enddo
			enddo

			call VALUE_HALF_NODE(4 ,nmax,nk,uvwt_line ,uvwt_half )
			call VALUE_HALF_NODE_IJK(1,1,0,12,nmax,nk,duvwt_line,duvwt_half)
			call VALUE_HALF_NODE(9 ,nmax,nk,kxyz_line ,kxyz_half )
			call VALUE_HALF_NODE(1 ,nmax,nk,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,nk,vslt1_line,vslt1_half  ) !*坐标：1:nk --> 0:nk
			call VALUE_HALF_NODE(1 ,nmax,nk,vslt2_line,vslt2_half  )

			do k=0,nk
				do m=9,12
					duvwt_half(m,k) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo

			call DUVWT_DXYZ(nmax,nk,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)
			
			do k=0,nk 
					nxyz_line(1,k) = kxyz_half(3,k)
					nxyz_line(2,k) = kxyz_half(6,k)
					nxyz_line(3,k) = kxyz_half(9,k)
			enddo 

			call FLUX_VIS_LINE_new(nl,nmax,nk,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    )

			call FLUX_DXYZ(nl,nmax,nk,fv,dfv)

			do k = 1,nk
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
				enddo
			enddo

	enddo
	enddo
!
	deallocate(duvwt)
	deallocate(duvwt_mid)
  return
end subroutine WCNSE5_VIS_FLUX_new
!=============================================================================!
!=============================================================================!
subroutine FLUX_VIS_LINE_new(nl,nmax,ni,n1,uvwt,n3,kxyz,n4,duvwtdxyz, &
                         vslt1,vslt2,fv)
!--------------------------------------------------------------------------!
!* 从0点开始计算
!* 该模块求粘性通量前，同时求剪切应力
!* 该模块传入的数组都是   0 : nmax
!* 该模块传入的工作数组是 0 : ni
!--------------------------------------------------------------------------!

	implicit none
	
	integer :: nmax,ni,i,m,nl,n1,n2,n3,n4
	real    :: uvwt(n1,0:nmax),kxyz(n3,0:nmax),fv(nl,0:nmax)
	real    :: duvwtdxyz(n4,0:nmax),vslt1(0:nmax),vslt2(0:nmax),kcp !,prl,prt
    real    :: txyz(9),dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	real    :: cc,vs,vscc,nx,ny,nz
			
    CC=2.0/3.0
	do i=0,ni      !___cic

		vs   = vslt1(i)

		vscc = vs*CC

		dudx = duvwtdxyz(1,i)
		dudy = duvwtdxyz(2,i)
		dudz = duvwtdxyz(3,i)

		dvdx = duvwtdxyz(4,i)
		dvdy = duvwtdxyz(5,i)
		dvdz = duvwtdxyz(6,i)

		dwdx = duvwtdxyz(7,i)
		dwdy = duvwtdxyz(8,i)
		dwdz = duvwtdxyz(9,i)
			

        txyz(1) = vscc * ( 2.0*dudx - dvdy - dwdz )
        txyz(5) = vscc * ( 2.0*dvdy - dwdz - dudx )
		txyz(9) = vscc * ( 2.0*dwdz - dudx - dvdy )
        txyz(2) = vs * ( dudy + dvdx )
        txyz(3) = vs * ( dudz + dwdx )
        txyz(6) = vs * ( dvdz + dwdy )
        txyz(4) = txyz(2)
        txyz(7) = txyz(3)
        txyz(8) = txyz(6)
		
		kcp = vslt2(i)

		nx  = kxyz(1,i)
		ny  = kxyz(2,i)
		nz  = kxyz(3,i)

		fv(1,i) = 0.0

		fv(2,i) = txyz(1)*nx+txyz(2)*ny+txyz(3)*nz
		fv(3,i) = txyz(4)*nx+txyz(5)*ny+txyz(6)*nz
		fv(4,i) = txyz(7)*nx+txyz(8)*ny+txyz(9)*nz

		fv(5,i) = uvwt(1,I)*fv(2,i)+uvwt(2,i)*fv(3,i)+uvwt(3,i)*fv(4,i)    +  &
							kcp*(duvwtdxyz(10,i)*nx+duvwtdxyz(11,i)*ny +  &
							duvwtdxyz(12,i)*nz )
		
	enddo

	return
end subroutine FLUX_VIS_LINE_new
!=============================================================================!
!=============================================================================!
subroutine FLUX_VIS_LINE_new1(nl,nmax,ni,n1,uvwt,n3,kxyz,n4,duvwtdxyz, &
                         vslt1,vslt2,fv)
!--------------------------------------------------------------------------!
!* 从1点开始计算
!* FLUX_VIS_LINE_new1与FLUX_VIS_LINE_new的数值起点不同
!* 传入前数组为：uvwt(n1,1:nmax),kxyz(n3,0:nmax),duvwtdxyz(n4,0:nmax)
!*               vslt1(1:nmax),wslt2(1:nmax),fv(nl,0:nmax)
!* 该模块传入的工作数组是 1 : ni
!--------------------------------------------------------------------------!

	implicit none
	
	integer :: nmax,ni,i,m,nl,n1,n2,n3,n4
	real    :: uvwt(n1,1:nmax),kxyz(n3,0:nmax),fv(nl,0:nmax)
	real    :: duvwtdxyz(n4,0:nmax),vslt1(1:nmax),vslt2(1:nmax),kcp !,prl,prt
    real    :: txyz(9),dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	real    :: cc,vs,vscc,nx,ny,nz
			
    CC=2.0/3.0
	do i=1,ni      !* 从1点开始计

		vs   = vslt1(i)

		vscc = vs*CC

		dudx = duvwtdxyz(1,i)
		dudy = duvwtdxyz(2,i)
		dudz = duvwtdxyz(3,i)

		dvdx = duvwtdxyz(4,i)
		dvdy = duvwtdxyz(5,i)
		dvdz = duvwtdxyz(6,i)

		dwdx = duvwtdxyz(7,i)
		dwdy = duvwtdxyz(8,i)
		dwdz = duvwtdxyz(9,i)
		
		nx = kxyz(1,i)	
		ny = kxyz(2,i)	
		nz = kxyz(3,i)	

        txyz(1) = vscc * ( 2.0*dudx - dvdy - dwdz )
        txyz(5) = vscc * ( 2.0*dvdy - dwdz - dudx )
		txyz(9) = vscc * ( 2.0*dwdz - dudx - dvdy )
        txyz(2) = vs * ( dudy + dvdx )
        txyz(3) = vs * ( dudz + dwdx )
        txyz(6) = vs * ( dvdz + dwdy )
        txyz(4) = txyz(2)
        txyz(7) = txyz(3)
        txyz(8) = txyz(6)
		
		kcp = vslt2(i)

		fv(1,i) = 0.0

		fv(2,i) = txyz(1)*nx+txyz(2)*ny+txyz(3)*nz
		fv(3,i) = txyz(4)*nx+txyz(5)*ny+txyz(6)*nz
		fv(4,i) = txyz(7)*nx+txyz(8)*ny+txyz(9)*nz

		fv(5,i) = uvwt(1,I)*fv(2,i)+uvwt(2,i)*fv(3,i)+uvwt(3,i)*fv(4,i)    +  &
							kcp*(duvwtdxyz(10,i)*nx+duvwtdxyz(11,i)*ny +  &
							duvwtdxyz(12,i)*nz )
		
	enddo

	return
end subroutine FLUX_VIS_LINE_new1
!=============================================================================!
!=============================================================================!
subroutine VALUE_HALF_NODE_IJK(IP,JP,KP,n,nmax,ni,q,q_half) !___INTERV2
!----------------------------------------------------------------------------!
!*TGH 从节点插值半结点 (N,1:NI) --> (N,0:NI)                                 !
!*TGH 1     : N/3   对应 i 方向，用ip标识是否求解                            !
!*TGH N/3+1 : 2*N/3 对应 j 方向，用jp标识是否求解                            !
!*TGH 2*N/3+1 : N   对应 k 方向，用kp标识是否求解                            !
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
	use global_variables,only : nijk2nd
	implicit none

	real,parameter :: A1=9.0,B1=-1.0
	real,parameter :: A2=5.0,B2=15.0,C2=-5.0,D2=1.0

	integer :: nmax,n,ni,m,i,IP,JP,KP,N3,NTEM,MST,MED
	real    :: q(n,nmax),q_half(n,0:nmax)

	N3 = N/3	
	
	if(mod(n,n3) /= 0 )then
		write(*,*)'在调用VALUE_HALF_NODE_IJK时，物理量个数n不对'
		stop
	endif

	if( ni <= nijk2nd )then  !!*TGH. 注意，可能是二维，降为2阶

		DO NTEM=1,IP
		    MST=1
			MED=N3
		    do i=1,ni-1
		    do m=MST,MED
				q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
			enddo
		    enddo
		    do m=MST,MED
			    q_half(m, 0) = 1.5*q(m, 1)-0.5*q(m,min(2,ni))  !___cic !*tgh. Modified by TU Guohua
			    q_half(m,ni) = 1.5*q(m,ni)-0.5*q(m,max(ni-1,1))  !___cic !*tgh. Modified by TU Guohua
			ENDDO
		ENDDO

		DO NTEM=1,jP
		    MST=N3+1
			MED=N3+N3
		    do i=1,ni-1
		    do m=MST,MED
				q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
			enddo
		    enddo
		    do m=MST,MED
			    q_half(m, 0) = 1.5*q(m, 1)-0.5*q(m,min(2,ni))  !___cic !*tgh. Modified by TU Guohua
			    q_half(m,ni) = 1.5*q(m,ni)-0.5*q(m,max(ni-1,1))  !___cic !*tgh. Modified by TU Guohua
			ENDDO
		ENDDO

		DO NTEM=1,kP
		    MST=N3+N3+1
			MED=N
		    do i=1,ni-1
		    do m=MST,MED
				q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
			enddo
		    enddo
		    do m=MST,MED
			    q_half(m, 0) = 1.5*q(m, 1)-0.5*q(m,min(2,ni))  !___cic !*tgh. Modified by TU Guohua
			    q_half(m,ni) = 1.5*q(m,ni)-0.5*q(m,max(ni-1,1))  !___cic !*tgh. Modified by TU Guohua
			ENDDO
		ENDDO

	else

		DO NTEM=1,IP
		    MST=1
			MED=N3
	        do i=2,ni-2  !
		    do m=MST,MED
				q_half(m,i) = (A1*(q(m,i) + q(m,i+1)) + B1*(q(m,i+2) + q(m,i-1)))/16.0
			enddo
			enddo

		    do m=MST,MED
			  q_half(m,1   ) = (A2*q(m,1 )+B2*q(m,2   )+C2*q(m,3   )+D2*q(m,4   ))/16.0
			  q_half(m,ni-1) = (A2*q(m,ni)+B2*q(m,ni-1)+C2*q(m,ni-2)+D2*q(m,ni-3))/16.0
			  q_half(m,0   ) = (35.d0*q(m,1 )-35.d0*q(m,2   )+21.d0*q(m,3   )-5.d0*q(m,4   ))/16.0 !cic
			  q_half(m,ni  ) = (35.d0*q(m,ni)-35.d0*q(m,ni-1)+21.d0*q(m,ni-2)-5.d0*q(m,ni-3))/16.0 !cic
!			  q_half(m,0   ) = (15.d0*q(m,1 )-10.d0*q(m,2   )+3.d0*q(m,3   ))/8.0 !*TGH. NEW WCNS_E_5_BORDER
!			  q_half(m,ni  ) = (15.d0*q(m,ni)-10.d0*q(m,ni-1)+3.d0*q(m,ni-2))/8.0 !*TGH. NEW WCNS_E_5_BORDER
			enddo

		enddo

		DO NTEM=1,jP
		    MST=N3+1
			MED=N3+N3
	        do i=2,ni-2  !
		    do m=MST,MED
				q_half(m,i) = (A1*(q(m,i) + q(m,i+1)) + B1*(q(m,i+2) + q(m,i-1)))/16.0
			enddo
			enddo

		    do m=MST,MED
			  q_half(m,1   ) = (A2*q(m,1 )+B2*q(m,2   )+C2*q(m,3   )+D2*q(m,4   ))/16.0
			  q_half(m,ni-1) = (A2*q(m,ni)+B2*q(m,ni-1)+C2*q(m,ni-2)+D2*q(m,ni-3))/16.0
			  q_half(m,0   ) = (35.d0*q(m,1 )-35.d0*q(m,2   )+21.d0*q(m,3   )-5.d0*q(m,4   ))/16.0 !cic
			  q_half(m,ni  ) = (35.d0*q(m,ni)-35.d0*q(m,ni-1)+21.d0*q(m,ni-2)-5.d0*q(m,ni-3))/16.0 !cic
!			  q_half(m,0   ) = (15.d0*q(m,1 )-10.d0*q(m,2   )+3.d0*q(m,3   ))/8.0 !*TGH. NEW WCNS_E_5_BORDER
!			  q_half(m,ni  ) = (15.d0*q(m,ni)-10.d0*q(m,ni-1)+3.d0*q(m,ni-2))/8.0 !*TGH. NEW WCNS_E_5_BORDER
			enddo

		enddo

		DO NTEM=1,kP
		    MST=N3+N3+1
			MED=N
	        do i=2,ni-2  !
		    do m=MST,MED
				q_half(m,i) = (A1*(q(m,i) + q(m,i+1)) + B1*(q(m,i+2) + q(m,i-1)))/16.0
			enddo
			enddo

		    do m=MST,MED
			  q_half(m,1   ) = (A2*q(m,1 )+B2*q(m,2   )+C2*q(m,3   )+D2*q(m,4   ))/16.0
			  q_half(m,ni-1) = (A2*q(m,ni)+B2*q(m,ni-1)+C2*q(m,ni-2)+D2*q(m,ni-3))/16.0
			  q_half(m,0   ) = (35.d0*q(m,1 )-35.d0*q(m,2   )+21.d0*q(m,3   )-5.d0*q(m,4   ))/16.0 !cic
			  q_half(m,ni  ) = (35.d0*q(m,ni)-35.d0*q(m,ni-1)+21.d0*q(m,ni-2)-5.d0*q(m,ni-3))/16.0 !cic
!			  q_half(m,0   ) = (15.d0*q(m,1 )-10.d0*q(m,2   )+3.d0*q(m,3   ))/8.0 !*TGH. NEW WCNS_E_5_BORDER
!			  q_half(m,ni  ) = (15.d0*q(m,ni)-10.d0*q(m,ni-1)+3.d0*q(m,ni-2))/8.0 !*TGH. NEW WCNS_E_5_BORDER
			enddo

		enddo

	endif

	return
end subroutine VALUE_HALF_NODE_IJK
!=============================================================================!
!=============================================================================!



!=============================================================================!
subroutine VIS_conservation
  use global_variables,only : u,v,w,t,reynolds,visl,vist,dq,nvis     &
                            , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol &
							, nl,nmax,ni,nj,nk,gama,moo,prl,prt
  use duvwt_all_field
  implicit none
!-----------------------------------------------------------------------------!
!                                                                             !
!  把同方向导数和交叉导数分离开，采用不同格式的求解方法                       !
!                                                                             !
!                                                                             !
! 特征注意：数值的起始坐标                                                    !
!                                                                             !
!  比WCNSE5_VIS_FLUX （1）多耗时18% —— 二阶中心计算交叉                     !
!                    （2）少耗时27% —— 不计算交叉                           !
!                                                                             !
!  设计：毛枚良，涂国华                                                       !
!  调试：涂国华 2009.03.15                                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
!	real,pointer,dimension(:,:,:,:) :: duvwt,duvwt_mid
	real    :: duvwtdxyz (12,0:nmax),nxyz_line(3,0:nmax)
	real    :: duvwt_line(12,nmax),uvwt_line(4,nmax)
	real    :: duvwt_half(12,0:nmax),uvwt_half(4,0:nmax)
	real    :: kxyz_line ( 9,1:nmax),kxyz_half(9,0:nmax)
	real    :: vol_line(nmax),fv(nl,0:nmax),dfv(nl,nmax)
	real    :: vol_half(0:nmax),txyz_half(9,0:nmax)
	real    :: vslt1_half(0:nmax),vslt2_half(0:nmax),cp,cp_prl,cp_prt
	real    :: vslt1_line(nmax),vslt2_line(nmax),re
	integer :: mvist,ICROSS

	ICROSS = 0 !特别注意，是否计算交叉导数

	allocate(duvwt(12,ni,nj,nk))   ! overcome stack over problem
	allocate(duvwt_mid(12,0:ni,0:nj,0:nk))   !___ 直接计算半结点1阶导数  2009.2.1

	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)

    cp_prl = cp/prl
    cp_prt = cp/prt

! to calculate the deriative of u,v,w,t in computing	coordinate

  call UVWT_DER_4th_ORDER (12) !*TGH. 直接求节点1阶导数，暂时没有用虚点值
  call UVWT_DER_4th_ORDER_half(12) !*TGH. 直接求半结点1阶导数

!**************** I direction ****************
	do k=1,nk 
    do j=1,nj

			do i=1,ni

				uvwt_line(1,i) = u(i,j,k)
				uvwt_line(2,i) = v(i,j,k)
				uvwt_line(3,i) = w(i,j,k)
				uvwt_line(4,i) = t(i,j,k)

				vslt1_line (i) = visl(i,j,k)
				vslt2_line (i) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line (i) = vslt1_line(i) + vist(i,j,k)
				  vslt2_line (i) = vslt2_line(i) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,i) = kcx(i,j,k)
				kxyz_line(2,i) = etx(i,j,k)
				kxyz_line(3,i) = ctx(i,j,k)

				kxyz_line(4,i) = kcy(i,j,k)
				kxyz_line(5,i) = ety(i,j,k)
				kxyz_line(6,i) = cty(i,j,k)

				kxyz_line(7,i) = kcz(i,j,k)
				kxyz_line(8,i) = etz(i,j,k)
				kxyz_line(9,i) = ctz(i,j,k)

				vol_line(   i) = vol(i,j,k)

				do m=5,12
					duvwt_line(m,i) = duvwt(m,i,j,k)
				enddo
			enddo
!=======================================================================================
!           求交叉导数形成的在整节点上粘性通量
!            DUVWT_DXYZ_CROSS(0,1,1..........)计算j和k方向导数形成的变量梯度
!            DUVWT_DXYZ_CROSS(1,0,1..........)计算i和k方向导数形成的变量梯度
!            DUVWT_DXYZ_CROSS(1,1,0..........)计算i和j方向导数形成的变量梯度
!            为了少开一个数组，此处txyz_half存储整节点的应力张量
!      注：其实可以将应力张量的计算可以放入粘性通量的计算子程序中更合理
!          可以，去掉 STRESS_LINE，FLUX_VIS_LINE也可以少传txyz_half

	                   IF(ICROSS == 1)THEN !* 是否计算交叉项

            call  DUVWT_DXYZ_CROSS(0,1,1,nmax,ni,12,duvwt_line,9,kxyz_line,vol_line,duvwtdxyz) !从1点开始计算

			do i=1,ni
			   nxyz_line(1,i) = kxyz_line(1,i)   
			   nxyz_line(2,i) = kxyz_line(4,i)
			   nxyz_line(3,i) = kxyz_line(7,i)
			enddo
			call  FLUX_VIS_LINE_new1(nl,nmax,ni,4,uvwt_line,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_line,vslt2_line,fv    ) !从1点开始计算

!			CALL DFLUX_VIS_CROSS_2nd(NL,NMAX,NI,FV,DFV) !* 求结点上交叉通量部分的一阶导数
			CALL DFLUX_VIS_CROSS_4th(NL,NMAX,NI,FV,DFV) !* 求结点上交叉通量部分的一阶导数

			do i = 1,ni
			do m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
			enddo
			enddo

			             ENDIF ! 求交叉导数形成的在整节点上粘性通量的贡献计算完毕 (I)
!==========================================================================================
!          (I) 求交叉导数形成的在整节点上粘性通量的贡献计算完毕
!==========================================================================================


!=========================================================================================
!             现计算同向导数对粘性项的贡献
!            DUVWT_DXYZ_single(1,0,0..........)计算i方向导数形成的变量梯度
!            DUVWT_DXYZ_single(0,1,0..........)计算j方向导数形成的变量梯度
!            DUVWT_DXYZ_single(0,0,1..........)计算k方向导数形成的变量梯度
!

			call VALUE_HALF_NODE(4 ,nmax,ni,uvwt_line ,uvwt_half ) !求原始变量在半结点上的值
	    	call VALUE_HALF_NODE(9 ,nmax,ni,kxyz_line ,kxyz_half ) !把网格导数插值到半结点
			call VALUE_HALF_NODE(1 ,nmax,ni,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,ni,vslT1_line,vslT1_half  )
			call VALUE_HALF_NODE(1 ,nmax,ni,vsLt2_line,vslt2_half  )

			do i=0,ni
				do m=1,4
					duvwt_half(m,i) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo


!			call DUVWT_DXYZ_single(1,0,0,nmax,ni,12,duvwt_mid,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
			call DUVWT_DXYZ_single(1,0,0,nmax,ni,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
!*TGH. 计算空间半结点上的一阶导数换算成物理空间的导数
		
			do i=0,ni
					nxyz_line(1,i) = kxyz_half(1,i)   !半结点上稳定在物理空间中的一阶导数
					nxyz_line(2,i) = kxyz_half(4,i)
					nxyz_line(3,i) = kxyz_half(7,i)
			enddo


			call  FLUX_VIS_LINE_NEW(nl,nmax,ni,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    ) !从0点开始计算
		  !*TGH. 按线计算粘性通量在计算空间半结点上的通量

			call	FLUX_DXYZ(nl,nmax,ni,fv,dfv)

			do i = 1,ni 
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
				enddo

			enddo

		enddo
	enddo
!
!**************** J direction ****************

	do k=1,nk 
		do i=1,ni 
			do j=1,nj

				uvwt_line(1,j) = u(i,j,k)
				uvwt_line(2,j) = v(i,j,k)
				uvwt_line(3,j) = w(i,j,k)
				uvwt_line(4,j) = t(i,j,k)

				vslt1_line(j) = visl(i,j,k)
				vslt2_line(j) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line(j) = vslt1_line(j) + vist(i,j,k)
				  vslt2_line(j) = vslt2_line(j) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,j) = kcx(i,j,k)
				kxyz_line(2,j) = etx(i,j,k)
				kxyz_line(3,j) = ctx(i,j,k)

				kxyz_line(4,j) = kcy(i,j,k)
				kxyz_line(5,j) = ety(i,j,k)
				kxyz_line(6,j) = cty(i,j,k)

				kxyz_line(7,j) = kcz(i,j,k)
				kxyz_line(8,j) = etz(i,j,k)
				kxyz_line(9,j) = ctz(i,j,k)

				vol_line(   j) = vol(i,j,k)

				do m=1,4
					duvwt_line(m,j) = duvwt(m,i,j,k)
				enddo
				do m=9,12
					duvwt_line(m,j) = duvwt(m,i,j,k)
				enddo

			enddo

!=======================================================================================
!           现在处理交叉导数部分 (J)

	                   IF(ICROSS == 1)THEN

            call  DUVWT_DXYZ_CROSS(1,0,1,nmax,nj,12,duvwt_line,9,kxyz_line,vol_line,duvwtdxyz) !从1点开始计算
			do j=1,nj
			   nxyz_line(1,j) = kxyz_line(2,j)   
			   nxyz_line(2,j) = kxyz_line(5,j)
			   nxyz_line(3,j) = kxyz_line(8,j)
			enddo
			call  FLUX_VIS_LINE_new1(nl,nmax,nj,4,uvwt_line,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_line,vslt2_line,fv    ) !从1点开始计算

!			CALL DFLUX_VIS_CROSS_2nd(NL,NMAX,NJ,FV,DFV) !* 求结点上交叉通量部分的一阶导数
			CALL DFLUX_VIS_CROSS_4th(NL,NMAX,NJ,FV,DFV) !* 求结点上交叉通量部分的一阶导数

			do J = 1,nJ
			do m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,J)
			enddo
			enddo

			             ENDIF ! 求交叉导数形成的在整节点上粘性通量的贡献计算完毕 (J)
!==========================================================================================
!          (J) 求交叉导数形成的在整节点上粘性通量的贡献计算完毕
!==========================================================================================

!             现计算同向导数对粘性项的贡献
				
			call VALUE_HALF_NODE(4 ,nmax,nj,uvwt_line ,uvwt_half )
			call VALUE_HALF_NODE(9 ,nmax,nj,kxyz_line ,kxyz_half )
			call VALUE_HALF_NODE(1 ,nmax,nj,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,nj,vslt1_line  ,vslt1_half  )
			call VALUE_HALF_NODE(1 ,nmax,nj,vslt2_line  ,vslt2_half  )

			do j=0,nj
				do m=5,8
					duvwt_half(m,j) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo

			call  DUVWT_DXYZ_SINGLE(0,1,0,nmax,nj,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
			
			do j=0,nj    !___cic
					nxyz_line(1,j) = kxyz_half(2,j)
					nxyz_line(2,j) = kxyz_half(5,j)
					nxyz_line(3,j) = kxyz_half(8,j)
			enddo 

			call  FLUX_VIS_LINE_new(nl,nmax,nj,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    ) !从0点开始计算

			call	FLUX_DXYZ(nl,nmax,nj,fv,dfv)

			do j = 1,nj
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j)
				enddo
			enddo

		enddo
	enddo
! 
!**************** K direction ****************
	do j=1,nj
		do i=1,ni
			do k=1,nk

				uvwt_line(1,k) = u(i,j,k)
				uvwt_line(2,k) = v(i,j,k)
				uvwt_line(3,k) = w(i,j,k)
				uvwt_line(4,k) = t(i,j,k)

				vslt1_line(k) = visl(i,j,k)
				vslt2_line(k) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line (k) = vslt1_line(k) + vist(i,j,k)
				  vslt2_line (k) = vslt2_line(k) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,k) = kcx(i,j,k)
				kxyz_line(2,k) = etx(i,j,k)
				kxyz_line(3,k) = ctx(i,j,k)

				kxyz_line(4,k) = kcy(i,j,k)
				kxyz_line(5,k) = ety(i,j,k)
				kxyz_line(6,k) = cty(i,j,k)

				kxyz_line(7,k) = kcz(i,j,k)
				kxyz_line(8,k) = etz(i,j,k)
				kxyz_line(9,k) = ctz(i,j,k)

				vol_line(   k) = vol(i,j,k)

				do m=1,8
					duvwt_line(m,k) = duvwt(m,i,j,k)
				enddo
			enddo

!=======================================================================================
!           现在处理交叉导数部分 (K)
	                   IF(ICROSS == 1)THEN
            call  DUVWT_DXYZ_CROSS(1,1,0,nmax,nK,12,duvwt_line,9,kxyz_line,vol_line,duvwtdxyz) !从1点开始计算
			do k=1,nk
			   nxyz_line(1,k) = kxyz_line(2,k)   
			   nxyz_line(2,k) = kxyz_line(5,k)
			   nxyz_line(3,k) = kxyz_line(8,k)
			enddo
			call  FLUX_VIS_LINE_new1(nl,nmax,nk,4,uvwt_line,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_line,vslt2_line,fv    ) !从1点开始计算

!			CALL DFLUX_VIS_CROSS_2nd(NL,NMAX,NK,FV,DFV) !* 求结点上交叉通量部分的一阶导数
			CALL DFLUX_VIS_CROSS_4th(NL,NMAX,NK,FV,DFV) !* 求结点上交叉通量部分的一阶导数

			do k = 1,nk
			do m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
			enddo
			enddo

			             ENDIF ! 求交叉导数形成的在整节点上粘性通量的贡献计算完毕 (k)
!==========================================================================================
!          (k) 求交叉导数形成的在整节点上粘性通量的贡献计算完毕
!==========================================================================================

!             现计算同向导数对粘性项的贡献

			call VALUE_HALF_NODE(4 ,nmax,nk,uvwt_line ,uvwt_half )
			call VALUE_HALF_NODE(9 ,nmax,nk,kxyz_line ,kxyz_half )
			call VALUE_HALF_NODE(1 ,nmax,nk,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,nk,vslt1_line  ,vslt1_half  )
			call VALUE_HALF_NODE(1 ,nmax,nk,vslt2_line  ,vslt2_half  )

			do k=0,nk
				do m=9,12
					duvwt_half(m,k) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo

			call  DUVWT_DXYZ_SINGLE(0,0,1,nmax,nk,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
			
!			do k=1,nk-1
			do k=0,nk    !___cic
					nxyz_line(1,k) = kxyz_half(3,k)
					nxyz_line(2,k) = kxyz_half(6,k)
					nxyz_line(3,k) = kxyz_half(9,k)
			enddo 

			call  FLUX_VIS_LINE_new(nl,nmax,nk,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    ) !从0点开始计算

			call	FLUX_DXYZ(nl,nmax,nk,fv,dfv)

			do k = 1,nk
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
				enddo
			enddo

		enddo
	enddo
!
	deallocate(duvwt)
	deallocate(duvwt_mid)
  return
end subroutine VIS_conservation
!=============================================================================!
!=============================================================================!
subroutine DUVWT_DXYZ_CROSS(r1,r2,r3,nmax,ni,n1,duvwt,n2,kxyz,vol,duvwtdxyz)
!*TGH. DUVWT_DXYZ_CROSS与DUVWT_DXYZ_single的功能相同
!*TGH.     但是数值的起始坐标不同，本程序除duvwtdxyz为(:,0:nmax)外
!*TGH.     其他都为为(:,1:nmax)
!*TGH. 半计算空间半结点上的一阶导数换算成物理空间的导数
!*TGH. 从1点开始计算，但是从duvwtdxyz点开始传数据
	implicit none

	integer :: nmax,ni,i,m,m1,m2,m3,n1,n2,n3,r1,r2,r3,n
	real,dimension(n1,1:nmax) :: duvwt
	real,dimension(n2,1:nmax) :: kxyz
	real,dimension(   1:nmax) :: vol
	real,dimension(n1,0:nmax) :: duvwtdxyz

	real,dimension(12) :: temp

	do i=1,ni  !0 点不用计算
		do m=1,3

			m1=M+M+M-2
			m2=m1+1
			m3=m2+1
			temp(m  ) = 0.0
			temp(m+3) = 0.0
			temp(m+6) = 0.0
			temp(m+9) = 0.0

            do n=1,r1   !m1表示i方向,r1=0剔除该方向对变量梯度的关系
			   temp(m  ) = temp(m  ) + kxyz(m1,i)*duvwt(1 ,i) 
			   temp(m+3) = temp(m+3) + kxyz(m1,i)*duvwt(2 ,i)
			   temp(m+6) = temp(m+6) + kxyz(m1,i)*duvwt(3 ,i) 
			   temp(m+9) = temp(m+9) + kxyz(m1,i)*duvwt(4 ,i) 
		    enddo

            do n=1,r2  !m2表示j方向,r2=0剔除该方向对变量梯度的关系
			   temp(m  ) = temp(m  ) + kxyz(m2,i)*duvwt(5 ,i) 
			   temp(m+3) = temp(m+3) + kxyz(m2,i)*duvwt(6 ,i)
			   temp(m+6) = temp(m+6) + kxyz(m2,i)*duvwt(7 ,i) 
			   temp(m+9) = temp(m+9) + kxyz(m2,i)*duvwt(8 ,i) 
		    enddo

            do n=1,r3  !m3表示k方向,r3=0剔除该方向对变量梯度的关系
			   temp(m  ) = temp(m  ) + kxyz(m3,i)*duvwt(9 ,i) 
			   temp(m+3) = temp(m+3) + kxyz(m3,i)*duvwt(10 ,i)
			   temp(m+6) = temp(m+6) + kxyz(m3,i)*duvwt(11 ,i) 
			   temp(m+9) = temp(m+9) + kxyz(m3,i)*duvwt(12 ,i) 
		    enddo

		   	duvwtdxyz(m  ,i) = temp(m  )/vol(i)
			duvwtdxyz(m+3,i) = temp(m+3)/vol(i)
			duvwtdxyz(m+6,i) = temp(m+6)/vol(i)
			duvwtdxyz(m+9,i) = temp(m+9)/vol(i)

		enddo
	enddo

	return
end subroutine DUVWT_DXYZ_CROSS   

!=============================================================================!
subroutine DUVWT_DXYZ_single(r1,r2,r3,nmax,ni,n1,duvwt,n2,kxyz,vol,duvwtdxyz)
!*TGH. DUVWT_DXYZ_CROSS与DUVWT_DXYZ_single的功能相同
!*TGH. 但是数值的起始坐标不同，本程序数值都为(:,0:nmax)
!*TGH. 半计算空间半结点上的一阶导数换算成物理空间的导数
!*TGH. 从0点开始计算

	implicit none

	integer :: nmax,ni,i,m,m1,m2,m3,n1,n2,n3,r1,r2,r3,n
	real,dimension(n1,0:nmax) :: duvwtdxyz
	real,dimension(n1,0:nmax) :: duvwt
	real,dimension(n2,0:nmax) :: kxyz
	real,dimension(   0:nmax) :: vol
	real,dimension(12) :: temp

	do i=0,ni
		do m=1,3

			m1=m+m+m-2
			m2=m1+1
			m3=m2+1
			temp(m  ) = 0.0
			temp(m+3) = 0.0
			temp(m+6) = 0.0
			temp(m+9) = 0.0

            do n=1,r1   !m1表示i方向,r1=0剔除该方向对变量梯度的关系
			   temp(m  ) = temp(m  ) + kxyz(m1,i)*duvwt(1 ,i) 
			   temp(m+3) = temp(m+3) + kxyz(m1,i)*duvwt(2 ,i)
			   temp(m+6) = temp(m+6) + kxyz(m1,i)*duvwt(3 ,i) 
			   temp(m+9) = temp(m+9) + kxyz(m1,i)*duvwt(4 ,i) 
		    enddo

            do n=1,r2  !m2表示j方向,r2=0剔除该方向对变量梯度的关系
			   temp(m  ) = temp(m  ) + kxyz(m2,i)*duvwt(5 ,i) 
			   temp(m+3) = temp(m+3) + kxyz(m2,i)*duvwt(6 ,i)
			   temp(m+6) = temp(m+6) + kxyz(m2,i)*duvwt(7 ,i) 
			   temp(m+9) = temp(m+9) + kxyz(m2,i)*duvwt(8 ,i) 
		    enddo

            do n=1,r3  !m3表示k方向,r3=0剔除该方向对变量梯度的关系
			   temp(m  ) = temp(m  ) + kxyz(m3,i)*duvwt(9 ,i) 
			   temp(m+3) = temp(m+3) + kxyz(m3,i)*duvwt(10 ,i)
			   temp(m+6) = temp(m+6) + kxyz(m3,i)*duvwt(11 ,i) 
			   temp(m+9) = temp(m+9) + kxyz(m3,i)*duvwt(12 ,i) 
		    enddo

		   	duvwtdxyz(m  ,i) = temp(m  )/vol(i)
			duvwtdxyz(m+3,i) = temp(m+3)/vol(i)
			duvwtdxyz(m+6,i) = temp(m+6)/vol(i)
			duvwtdxyz(m+9,i) = temp(m+9)/vol(i)

		enddo
	enddo

	return
end subroutine DUVWT_DXYZ_single


SUBROUTINE DFLUX_VIS_CROSS_2nd(NL,NMAX,NI,FV,DFV)
!------------------------------------------------------------------!
!* 已知结点fv，求结点dfv的一阶导数（计算空间）                     !
!* 暂时用二阶                                                      !
!* 特别注意：考虑到调用该程序时的FV(NL,0:NMAX),DFV(NL,1:NMAX)      !
!------------------------------------------------------------------!
	use global_variables,only : nijk2d
	USE DEFINE_PRECISION_MOD
	IMPLICIT NONE
	INTEGER :: I,NL,NMAX,NI,M
	REAL(PREC) :: FV(NL,0:NMAX),DFV(NL,1:NMAX)
	REAL(PREC) :: B1,B2,B3,B4,COE

	COE =  0.5_PREC
	B1  = -1.5_PREC
	B2  =  2.0_PREC
	B3  = -0.5_PREC 
	
	IF( NI <= nijk2d ) THEN !*TGH. 注意，可能是二维，将为2阶
         DFV(:,:) = 0.0
		RETURN
	ENDIF

	DO I=2,NI-1
	DO M=1,NL
		DFV(M,I) = COE * ( FV(M,I+1) - FV(M,I-1) )
	ENDDO
	ENDDO

	DO M=1,NL
		DFV(M,1) = B1*FV(M,1) + B2*FV(M,2) + B3*FV(M,3)
		DFV(M,NI) = -B1*FV(M,NI) - B2*FV(M,NI-1) - B3*FV(M,NI-2)
	ENDDO

END SUBROUTINE DFLUX_VIS_CROSS_2nd

!=============================================================================!
SUBROUTINE DFLUX_VIS_CROSS_4th(NL,NMAX,NI,FV,DFV)
!------------------------------------------------------------------!
!* 已知结点fv，求结点dfv的一阶导数（计算空间）                     !
!* 暂时用二阶                                                      !
!* 特别注意：考虑到调用该程序时的FV(NL,0:NMAX),DFV(NL,1:NMAX)      !
!------------------------------------------------------------------!
	use global_variables,only : nijk2d
	USE DEFINE_PRECISION_MOD
	IMPLICIT NONE
	INTEGER :: I,NL,NMAX,NI,M
	REAL(PREC) :: FV(NL,0:NMAX),DFV(NL,1:NMAX)
	REAL(PREC) :: B11,B12,B13,B14,b21,b22,b23,b24,coe1,coe2,dd

	COE1 =  8.0_prec/12._prec
	coe2 =   -1._prec/12._prec

	b11 = -11._prec
	b12 =  18._prec
	b13 =   9._prec
	b14 =  -2._prec

	b21 =  -2._prec
	b22 =  -3._prec
	b23 =   6._prec
	b24 =  -1._prec

	dd  = 6.0_prec

	IF( NI <= nijk2d) THEN !*TGH. 注意，可能是二维，将为2阶
         DFV(:,:) =  0.0
		RETURN
	ENDIF

	DO I=3,NI-2
	DO M=1,NL
		DFV(M,I) = COE1 * ( FV(M,I+1) - FV(M,I-1) ) + COE2 * ( FV(M,I+2) - FV(M,I-2) )
	ENDDO
	ENDDO

	DO M=1,NL
		DFV(M,1)    = ( B11*FV(M,1 ) + B12*FV(M,2   ) + B13*FV(M,3   ) + b14*FV(M,4)    )/dd
		DFV(M,NI)   =-( B11*FV(M,NI) + B12*FV(M,NI-1) + B13*FV(M,NI-2) + B14*FV(M,NI-3) )/dd
		DFV(M,2)    = ( B21*FV(M,1 ) + B22*FV(M,2   ) + B23*FV(M,3   ) + b24*FV(M,4)    )/dd
		DFV(M,NI-1) =-( B21*FV(M,NI) + B22*FV(M,NI-1) + B23*FV(M,NI-2) + B24*FV(M,NI-3) )/dd
	ENDDO


END SUBROUTINE DFLUX_VIS_CROSS_4th

!=============================================================================!!=============================================================================!

!=============================================================================!
subroutine VIS_conservation_virtual
  use global_variables,only : u,v,w,t,reynolds,visl,vist,dq,nvis     &
                            , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol &
							, nl,nmax,ni,nj,nk,gama,moo,prl,prt
  use duvwt_all_field
  implicit none
!-----------------------------------------------------------------------------!
!                                                                             !
!  把同方向导数和交叉导数分离开，采用不同格式的求解方法                       !
!  需要使用虚点上的值                                                         !
!                                                                             !
! 特征注意：数值的起始坐标                                                    !
!                                                                             !
!  比WCNSE5_VIS_FLUX （1）多耗时18% —— 二阶中心计算交叉                     !
!                    （2）少耗时27% —— 不计算交叉                           !
!                                                                             !
!  设计：毛枚良，涂国华                                                       !
!  调试：涂国华 2009.03.15                                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
!	real,pointer,dimension(:,:,:,:) :: duvwt,duvwt_mid
	real    :: duvwtdxyz (12,0:nmax),nxyz_line(3,0:nmax)
	real    :: duvwt_line(12,nmax),uvwt_line(4,nmax),uvwt_line_vir(4,0:nmax+1)
	real    :: duvwt_half(12,0:nmax),uvwt_half(4,0:nmax)
	real    :: kxyz_line ( 9,1:nmax),kxyz_half(9,0:nmax)
	real    :: vol_line(nmax),fv(nl,0:nmax),dfv(nl,nmax)
	real    :: vol_half(0:nmax),txyz_half(9,0:nmax)
	real    :: vslt1_half(0:nmax),vslt2_half(0:nmax),cp,cp_prl,cp_prt
	real    :: vslt1_line(nmax),vslt2_line(nmax),re
	integer :: mvist,ICROSS

	ICROSS = 1 !特别注意，是否计算交叉导数

	allocate(duvwt(12,ni,nj,nk))   ! overcome stack over problem
	allocate(duvwt_mid(12,0:ni,0:nj,0:nk))   !___ 直接计算半结点1阶导数  2009.2.1

	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)

    cp_prl = cp/prl
    cp_prt = cp/prt

! to calculate the deriative of u,v,w,t in computing	coordinate

  call UVWT_DER_4th_vir !*TGH. 直接求节点1阶导数，暂时没有用虚点值
  call UVWT_DER_4th_half_vir !*TGH. 直接求半结点1阶导数

!**************** I direction ****************
	do k=1,nk 
    do j=1,nj

			do i=1,ni

				uvwt_line(1,i) = u(i,j,k)
				uvwt_line(2,i) = v(i,j,k)
				uvwt_line(3,i) = w(i,j,k)
				uvwt_line(4,i) = t(i,j,k)

				do m=1,4
					uvwt_line_vir(m,i) = uvwt_line(m,i)
				enddo

				vslt1_line (i) = visl(i,j,k)
				vslt2_line (i) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line (i) = vslt1_line(i) + vist(i,j,k)
				  vslt2_line (i) = vslt2_line(i) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,i) = kcx(i,j,k)
				kxyz_line(2,i) = etx(i,j,k)
				kxyz_line(3,i) = ctx(i,j,k)

				kxyz_line(4,i) = kcy(i,j,k)
				kxyz_line(5,i) = ety(i,j,k)
				kxyz_line(6,i) = cty(i,j,k)

				kxyz_line(7,i) = kcz(i,j,k)
				kxyz_line(8,i) = etz(i,j,k)
				kxyz_line(9,i) = ctz(i,j,k)

				vol_line(   i) = vol(i,j,k)

				do m=5,12
					duvwt_line(m,i) = duvwt(m,i,j,k)
				enddo
			enddo
!=======================================================================================
!           求交叉导数形成的在整节点上粘性通量
!            DUVWT_DXYZ_CROSS(0,1,1..........)计算j和k方向导数形成的变量梯度
!            DUVWT_DXYZ_CROSS(1,0,1..........)计算i和k方向导数形成的变量梯度
!            DUVWT_DXYZ_CROSS(1,1,0..........)计算i和j方向导数形成的变量梯度
!            为了少开一个数组，此处txyz_half存储整节点的应力张量
!      注：其实可以将应力张量的计算可以放入粘性通量的计算子程序中更合理
!          可以，去掉 STRESS_LINE，FLUX_VIS_LINE也可以少传txyz_half

	                   IF(ICROSS == 1)THEN !* 是否计算交叉项

            call  DUVWT_DXYZ_CROSS(0,1,1,nmax,ni,12,duvwt_line,9,kxyz_line,vol_line,duvwtdxyz) !从1点开始计算

			do i=1,ni
			   nxyz_line(1,i) = kxyz_line(1,i)   
			   nxyz_line(2,i) = kxyz_line(4,i)
			   nxyz_line(3,i) = kxyz_line(7,i)
			enddo
			call  FLUX_VIS_LINE_new1(nl,nmax,ni,4,uvwt_line,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_line,vslt2_line,fv    ) !从1点开始计算

!			CALL DFLUX_VIS_CROSS_2nd(NL,NMAX,NI,FV,DFV) !* 求结点上交叉通量部分的一阶导数
			CALL DFLUX_VIS_CROSS_4th(NL,NMAX,NI,FV,DFV) !* 求结点上交叉通量部分的一阶导数

			do i = 1,ni
			do m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
			enddo
			enddo

			             ENDIF ! 求交叉导数形成的在整节点上粘性通量的贡献计算完毕 (I)
!==========================================================================================
!          (I) 求交叉导数形成的在整节点上粘性通量的贡献计算完毕
!==========================================================================================


!=========================================================================================
!             现计算同向导数对粘性项的贡献
!            DUVWT_DXYZ_single(1,0,0..........)计算i方向导数形成的变量梯度
!            DUVWT_DXYZ_single(0,1,0..........)计算j方向导数形成的变量梯度
!            DUVWT_DXYZ_single(0,0,1..........)计算k方向导数形成的变量梯度
!
			i=0
			uvwt_line_vir(1,i) = u(i,j,k)
			uvwt_line_vir(2,i) = v(i,j,k)
			uvwt_line_vir(3,i) = w(i,j,k)
			uvwt_line_vir(4,i) = t(i,j,k)
			i=ni+1
			uvwt_line_vir(1,i) = u(i,j,k)
			uvwt_line_vir(2,i) = v(i,j,k)
			uvwt_line_vir(3,i) = w(i,j,k)
			uvwt_line_vir(4,i) = t(i,j,k)

			call VALUE_HALF_NODE_vir(4 ,nmax,ni,uvwt_line_vir ,uvwt_half ) !求原始变量在半结点上的值
	    	call VALUE_HALF_NODE(9 ,nmax,ni,kxyz_line ,kxyz_half ) !把网格导数插值到半结点
			call VALUE_HALF_NODE(1 ,nmax,ni,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,ni,vslT1_line,vslT1_half  )
			call VALUE_HALF_NODE(1 ,nmax,ni,vsLt2_line,vslt2_half  )

			do i=0,ni
				do m=1,4
					duvwt_half(m,i) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo


!			call DUVWT_DXYZ_single(1,0,0,nmax,ni,12,duvwt_mid,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
			call DUVWT_DXYZ_single(1,0,0,nmax,ni,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
!*TGH. 计算空间半结点上的一阶导数换算成物理空间的导数
		
			do i=0,ni
					nxyz_line(1,i) = kxyz_half(1,i)   !半结点上稳定在物理空间中的一阶导数
					nxyz_line(2,i) = kxyz_half(4,i)
					nxyz_line(3,i) = kxyz_half(7,i)
			enddo


			call  FLUX_VIS_LINE_NEW(nl,nmax,ni,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    ) !从0点开始计算
		  !*TGH. 按线计算粘性通量在计算空间半结点上的通量

			call	FLUX_DXYZ(nl,nmax,ni,fv,dfv)

			do i = 1,ni 
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
				enddo

			enddo

		enddo
	enddo
!
!**************** J direction ****************

	do k=1,nk 
		do i=1,ni 
			do j=1,nj

				uvwt_line(1,j) = u(i,j,k)
				uvwt_line(2,j) = v(i,j,k)
				uvwt_line(3,j) = w(i,j,k)
				uvwt_line(4,j) = t(i,j,k)

				do m=1,4
					uvwt_line_vir(m,j) = uvwt_line(m,j)
				enddo

				vslt1_line(j) = visl(i,j,k)
				vslt2_line(j) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line(j) = vslt1_line(j) + vist(i,j,k)
				  vslt2_line(j) = vslt2_line(j) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,j) = kcx(i,j,k)
				kxyz_line(2,j) = etx(i,j,k)
				kxyz_line(3,j) = ctx(i,j,k)

				kxyz_line(4,j) = kcy(i,j,k)
				kxyz_line(5,j) = ety(i,j,k)
				kxyz_line(6,j) = cty(i,j,k)

				kxyz_line(7,j) = kcz(i,j,k)
				kxyz_line(8,j) = etz(i,j,k)
				kxyz_line(9,j) = ctz(i,j,k)

				vol_line(   j) = vol(i,j,k)

				do m=1,4
					duvwt_line(m,j) = duvwt(m,i,j,k)
				enddo
				do m=9,12
					duvwt_line(m,j) = duvwt(m,i,j,k)
				enddo

			enddo

!=======================================================================================
!           现在处理交叉导数部分 (J)

	                   IF(ICROSS == 1)THEN

            call  DUVWT_DXYZ_CROSS(1,0,1,nmax,nj,12,duvwt_line,9,kxyz_line,vol_line,duvwtdxyz) !从1点开始计算
			do j=1,nj
			   nxyz_line(1,j) = kxyz_line(2,j)   
			   nxyz_line(2,j) = kxyz_line(5,j)
			   nxyz_line(3,j) = kxyz_line(8,j)
			enddo
			call  FLUX_VIS_LINE_new1(nl,nmax,nj,4,uvwt_line,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_line,vslt2_line,fv    ) !从1点开始计算

!			CALL DFLUX_VIS_CROSS_2nd(NL,NMAX,NJ,FV,DFV) !* 求结点上交叉通量部分的一阶导数
			CALL DFLUX_VIS_CROSS_4th(NL,NMAX,NJ,FV,DFV) !* 求结点上交叉通量部分的一阶导数

			do J = 1,nJ
			do m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,J)
			enddo
			enddo

			             ENDIF ! 求交叉导数形成的在整节点上粘性通量的贡献计算完毕 (J)
!==========================================================================================
!          (J) 求交叉导数形成的在整节点上粘性通量的贡献计算完毕
!==========================================================================================

!             现计算同向导数对粘性项的贡献
				
			j=0
			uvwt_line_vir(1,j) = u(i,j,k)
			uvwt_line_vir(2,j) = v(i,j,k)
			uvwt_line_vir(3,j) = w(i,j,k)
			uvwt_line_vir(4,j) = t(i,j,k)
			j=nj+1
			uvwt_line_vir(1,j) = u(i,j,k)
			uvwt_line_vir(2,j) = v(i,j,k)
			uvwt_line_vir(3,j) = w(i,j,k)
			uvwt_line_vir(4,j) = t(i,j,k)

			call VALUE_HALF_NODE_vir(4 ,nmax,nj,uvwt_line_vir ,uvwt_half )
			call VALUE_HALF_NODE(9 ,nmax,nj,kxyz_line ,kxyz_half )
			call VALUE_HALF_NODE(1 ,nmax,nj,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,nj,vslt1_line  ,vslt1_half  )
			call VALUE_HALF_NODE(1 ,nmax,nj,vslt2_line  ,vslt2_half  )

			do j=0,nj
				do m=5,8
					duvwt_half(m,j) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo

			call  DUVWT_DXYZ_SINGLE(0,1,0,nmax,nj,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
			
			do j=0,nj    !___cic
					nxyz_line(1,j) = kxyz_half(2,j)
					nxyz_line(2,j) = kxyz_half(5,j)
					nxyz_line(3,j) = kxyz_half(8,j)
			enddo 

			call  FLUX_VIS_LINE_new(nl,nmax,nj,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    ) !从0点开始计算

			call	FLUX_DXYZ(nl,nmax,nj,fv,dfv)

			do j = 1,nj
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j)
				enddo
			enddo

		enddo
	enddo
! 
!**************** K direction ****************
	do j=1,nj
		do i=1,ni
			do k=1,nk

				uvwt_line(1,k) = u(i,j,k)
				uvwt_line(2,k) = v(i,j,k)
				uvwt_line(3,k) = w(i,j,k)
				uvwt_line(4,k) = t(i,j,k)

				do m=1,4
					uvwt_line_vir(m,k) = uvwt_line(m,k)
				enddo

				vslt1_line(k) = visl(i,j,k)
				vslt2_line(k) = visl(i,j,k)*cp_prl
			    do mvist=1,nvis !*湍流时
				  vslt1_line (k) = vslt1_line(k) + vist(i,j,k)
				  vslt2_line (k) = vslt2_line(k) + vist(i,j,k)*cp_prt
				enddo

				kxyz_line(1,k) = kcx(i,j,k)
				kxyz_line(2,k) = etx(i,j,k)
				kxyz_line(3,k) = ctx(i,j,k)

				kxyz_line(4,k) = kcy(i,j,k)
				kxyz_line(5,k) = ety(i,j,k)
				kxyz_line(6,k) = cty(i,j,k)

				kxyz_line(7,k) = kcz(i,j,k)
				kxyz_line(8,k) = etz(i,j,k)
				kxyz_line(9,k) = ctz(i,j,k)

				vol_line(   k) = vol(i,j,k)

				do m=1,8
					duvwt_line(m,k) = duvwt(m,i,j,k)
				enddo
			enddo

!=======================================================================================
!           现在处理交叉导数部分 (K)
	                   IF(ICROSS == 1)THEN
            call  DUVWT_DXYZ_CROSS(1,1,0,nmax,nK,12,duvwt_line,9,kxyz_line,vol_line,duvwtdxyz) !从1点开始计算
			do k=1,nk
			   nxyz_line(1,k) = kxyz_line(2,k)   
			   nxyz_line(2,k) = kxyz_line(5,k)
			   nxyz_line(3,k) = kxyz_line(8,k)
			enddo
			call  FLUX_VIS_LINE_new1(nl,nmax,nk,4,uvwt_line,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_line,vslt2_line,fv    ) !从1点开始计算

!			CALL DFLUX_VIS_CROSS_2nd(NL,NMAX,NK,FV,DFV) !* 求结点上交叉通量部分的一阶导数
			CALL DFLUX_VIS_CROSS_4th(NL,NMAX,NK,FV,DFV) !* 求结点上交叉通量部分的一阶导数

			do k = 1,nk
			do m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
			enddo
			enddo

			             ENDIF ! 求交叉导数形成的在整节点上粘性通量的贡献计算完毕 (k)
!==========================================================================================
!          (k) 求交叉导数形成的在整节点上粘性通量的贡献计算完毕
!==========================================================================================

!             现计算同向导数对粘性项的贡献

			k=0
			uvwt_line_vir(1,k) = u(i,j,k)
			uvwt_line_vir(2,k) = v(i,j,k)
			uvwt_line_vir(3,k) = w(i,j,k)
			uvwt_line_vir(4,k) = t(i,j,k)
			k=nk+1
			uvwt_line_vir(1,k) = u(i,j,k)
			uvwt_line_vir(2,k) = v(i,j,k)
			uvwt_line_vir(3,k) = w(i,j,k)
			uvwt_line_vir(4,k) = t(i,j,k)

			call VALUE_HALF_NODE_vir(4 ,nmax,nk,uvwt_line_vir ,uvwt_half )
			call VALUE_HALF_NODE(9 ,nmax,nk,kxyz_line ,kxyz_half )
			call VALUE_HALF_NODE(1 ,nmax,nk,vol_line  ,vol_half  )
			call VALUE_HALF_NODE(1 ,nmax,nk,vslt1_line  ,vslt1_half  )
			call VALUE_HALF_NODE(1 ,nmax,nk,vslt2_line  ,vslt2_half  )

			do k=0,nk
				do m=9,12
					duvwt_half(m,k) = duvwt_mid(m,i,j,k)    !___同向导数
				enddo
			enddo

			call  DUVWT_DXYZ_SINGLE(0,0,1,nmax,nk,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz) !从0点开始计算
			
!			do k=1,nk-1
			do k=0,nk    !___cic
					nxyz_line(1,k) = kxyz_half(3,k)
					nxyz_line(2,k) = kxyz_half(6,k)
					nxyz_line(3,k) = kxyz_half(9,k)
			enddo 

			call  FLUX_VIS_LINE_new(nl,nmax,nk,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    ) !从0点开始计算

			call	FLUX_DXYZ(nl,nmax,nk,fv,dfv)

			do k = 1,nk
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
				enddo
			enddo

		enddo
	enddo
!
	deallocate(duvwt)
	deallocate(duvwt_mid)
  return
end subroutine VIS_conservation_virtual
!=============================================================================!
subroutine UVWT_DER_4th_half_vir    !____直接求半结点1阶导数 2009.2.1
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt_mid
  implicit none
	integer :: i,j,k,m,nv
	real    :: uvwt(4,-1:nmax+1),duvwt(4,0:nmax)
!	real    :: duvwtdkc(nv,0:ni,0:nj,0:nk)

! I direction

	do k=1,nk
		do j=1,nj
			do i=0,ni+1

				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)

			enddo

			call DUVWT_halfNODE(nmax,ni,4,uvwt,duvwt)     !*TGH. 要用虚点上的值

			do i=0,ni
				do m=1,4

					duvwt_mid(m,i,j,k)=duvwt(m,i)

				enddo
			enddo

		enddo
	enddo
!
! J direction

	do k=1,nk
		do i=1,ni
			do j=0,nj+1

				uvwt(1,j) = u(i,j,k)
				uvwt(2,j) = v(i,j,k)
				uvwt(3,j) = w(i,j,k)
				uvwt(4,j) = t(i,j,k)

			enddo

			call DUVWT_halfNODE(nmax,nj,4,uvwt,duvwt)     !*TGH. 要用虚点上的值

			do j=0,nj
				do m=1,4

					duvwt_mid(m+4,i,j,k)=duvwt(m,j)

				enddo
			enddo

		enddo
	enddo

! K direction

	do j=1,nj
		do i=1,ni
			do k=0,nk+1

				uvwt(1,k) = u(i,j,k)
				uvwt(2,k) = v(i,j,k)
				uvwt(3,k) = w(i,j,k)
				uvwt(4,k) = t(i,j,k)

			enddo

			call DUVWT_halfNODE(nmax,nk,4,uvwt,duvwt)     !*TGH. 要用虚点上的值

			do k=0,nk
				do m=1,4

					duvwt_mid(m+8,i,j,k)=duvwt(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_4th_half_vir
!=============================================================================!
!=============================================================================!
subroutine UVWT_DER_4th_vir
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt
  implicit none
	integer :: i,j,k,m,nv
	real    :: uvwt(4,-1:nmax+1),dtem(4,nmax)
!	real    :: duvwtdkc(nv,ni,nj,nk)

! I direction

	do k=1,nk
		do j=1,nj
			do i=0,ni+1

				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)

			enddo

			call DUVWT_NODE_vir(nmax,ni,4,uvwt,dtem) !*TGH. 暂时没有用虚点值

			do i=1,ni
				do m=1,4

					duvwt(m,i,j,k)=dtem(m,i)

				enddo
			enddo

		enddo
	enddo
!
! J direction

	do k=1,nk
		do i=1,ni
			do j=0,nj+1

				uvwt(1,j) = u(i,j,k)
				uvwt(2,j) = v(i,j,k)
				uvwt(3,j) = w(i,j,k)
				uvwt(4,j) = t(i,j,k)

			enddo

			call DUVWT_NODE_vir(nmax,nj,4,uvwt,dtem)

			do j=1,nj
				do m=1,4

					duvwt(m+4,i,j,k)=dtem(m,j)

				enddo
			enddo

		enddo
	enddo

! K direction

	do j=1,nj
		do i=1,ni
			do k=0,nk+1

				uvwt(1,k) = u(i,j,k)
				uvwt(2,k) = v(i,j,k)
				uvwt(3,k) = w(i,j,k)
				uvwt(4,k) = t(i,j,k)

			enddo

			call DUVWT_NODE_vir(nmax,nk,4,uvwt,dtem)

			do k=1,nk
				do m=1,4

					duvwt(m+8,i,j,k)=dtem(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_4th_vir
!=============================================================================!
!=============================================================================!
subroutine DUVWT_NODE_vir(nmax,ni,n1,uvwt,duvwt)  !___DERINODE
	use global_variables,only : nijk2nd
    use define_precision_mod
	implicit none

	real(prec) :: A1,B1
	real(prec) :: A2,B2,C2,D2,E2

	integer :: nmax,i,m,ni,n1,n2
	real(prec) :: uvwt(n1,-1:nmax+1),duvwt(n1,nmax)

	A1=8._prec ; B1=-1._prec
	A2=-3._prec ; B2=-10._prec ; C2=18._prec ; D2=-6._prec ; E2=1._prec

	if( ni <= nijk2nd ) then   !*TGH. 注意，可能是二维，将为2阶
		do m=1,n1
			do i=2,ni-1
				duvwt(m,i) = 0.5*(uvwt(m,i+1) - uvwt(m,i-1))
			enddo
				duvwt(m,1) = uvwt(m,2 ) - uvwt(m,1)
				duvwt(m,ni)= uvwt(m,ni) - uvwt(m,ni-1)
		enddo
	
	else
		do m=1,n1
			do i=3,ni-2
				duvwt(m,i) = (A1*(uvwt(m,i+1) - uvwt(m,i-1)) + &
										  B1*(uvwt(m,i+2) - uvwt(m,i-2)))/12.d0
			enddo
			duvwt(m,2   ) =  (A2*uvwt(m,1) + B2*uvwt(m,2) + C2*uvwt(m,3) + &
			                  D2*uvwt(m,4) + E2*uvwt(m,5))/12.d0
			duvwt(m,ni-1) = -(A2*uvwt(m,ni  ) + B2*uvwt(m,ni-1) + C2*uvwt(m,ni-2) + &
			                  D2*uvwt(m,ni-3) + E2*uvwt(m,ni-4))/12.d0

			duvwt(m,1   ) =  (A2*uvwt(m,0) + B2*uvwt(m,1) + C2*uvwt(m,2) + &
			                  D2*uvwt(m,3) + E2*uvwt(m,4))/12._prec
			duvwt(m,ni  ) = -(A2*uvwt(m,ni+1) + B2*uvwt(m,ni  ) + C2*uvwt(m,ni-1) + &
			                  D2*uvwt(m,ni-2) + E2*uvwt(m,ni-3))/12._prec

		enddo
	endif
!
  return
end subroutine DUVWT_NODE_vir
!=============================================================================!
!=============================================================================!
subroutine VALUE_HALF_NODE_vir(n,nmax,ni,q,q_half) !___INTERV2
	use global_variables,only : nijk2nd
	use define_precision_mod
	implicit none
	
	real(prec) :: A1,B1
	real(prec) :: A2,B2,C2,D2

	integer :: nmax,n,ni,m,i
	real(prec) :: q(n,0:nmax+1),q_half(n,0:nmax)
!
! deal with the symmetric boundary condiction
!
	if( ni <= nijk2nd)then !*TGH. 注意，可能是二维，将为2阶

		do m=1,n
			do i=1,ni-1
				q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
			enddo
			  q_half(m, 0) = 1.5*q(m, 1)-0.5*q(m,min(2,ni))  !___cic !*tgh. Modified by TU Guohua
			  q_half(m,ni) = 1.5*q(m,ni)-0.5*q(m,max(ni-1,1))  !___cic !*tgh. Modified by TU Guohua
		enddo

	else

	    A1=9.0_prec ; B1=-1.0_prec
		A2=5.0_prec ; B2=15.0_prec ; C2=-5.0_prec ; D2=1.0_prec

		do m=1,n
			do i=2,ni-2  !
				q_half(m,i) = (A1*(q(m,i) + q(m,i+1)) + B1*(q(m,i+2) + q(m,i-1)))/16.0_prec
			enddo

			q_half(m,1   ) = (A2*q(m,1 )+B2*q(m,2   )+C2*q(m,3   )+D2*q(m,4   ))/16.0_prec
			q_half(m,ni-1) = (A2*q(m,ni)+B2*q(m,ni-1)+C2*q(m,ni-2)+D2*q(m,ni-3))/16.0_prec

			q_half(m,0   ) = (A2*q(m,0   )+B2*q(m,1   )+C2*q(m,2   )+D2*q(m,3   ))/16.0_prec
			q_half(m,ni  ) = (A2*q(m,ni+1)+B2*q(m,ni  )+C2*q(m,ni-1)+D2*q(m,ni-2))/16.0_prec

		enddo

	endif

	return
end subroutine VALUE_HALF_NODE_vir
!=============================================================================!
!=============================================================================!
!=============================================================================!
subroutine VIS_nonconservation
  use define_precision_mod
  use global_variables,only : u,v,w,t,reynolds,visl,vist,dq,nvis     &
                            , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol &
							, nl,nmax,ni,nj,nk,gama,moo,prl,prt
  use duvwt_all_field,only : duvwt
  implicit none
!-----------------------------------------------------------------------------!
!                                                                             !
!  把同方向导数和交叉导数分离开，采用不同格式的求解方法                       !
!  注意：求得的粘性项贡献对应的时间离散为：(dQ/Dt)/J                          !
!                                                                             !
!                                                                             !
! 特别注意：数组的起始坐标                                                    !
!                                                                             !
!  比WCNSE5_VIS_FLUX （1）少耗时3.7% —— 二阶中心计算交叉(还需要测试正确性） !
!                    （2）少耗时42%  —— 不计算交叉                           !
!                                                                             !
!  设计：毛枚良，涂国华                                                       !
!  调试：涂国华 2009.03.20                                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
	real(prec) :: duvwtdxyz(12,nmax),nxyz_line(3,nmax)
	real(prec) :: uvwt_line(4,nmax),duvwt_line(12,nmax)
	real(prec) :: d1uvwt(4,nmax),d2uvwt(4,nmax) !tgh
	real(prec) :: kxyz_line(9,nmax) !,kxyz_half(9,0:nmax)
	real(prec) :: vol_line(nmax),fv(nl,nmax),dfv(nl,nmax)
	real(prec) :: vslt1_line(nmax),vslt2_line(nmax)
	real(prec) :: cp,cp_prl,cp_prt,re,vslt1,vslt2
	integer :: mvist,ICROSS,mm
	real(prec) :: volp !,vslt1p,vslt2p
	real(prec) :: nx,ny,nz !,nxyz_line(3,nmax) !tgh
!	real(prec) :: d1uvwt(4,nmax),d2uvwt(4,nmax)

	ICROSS = 1 !特别注意，是否计算交叉导数

	allocate(duvwt(12,ni,nj,nk))

    call D1UVWT_4ORDER    !*TGH. 直接求结点1阶导数



	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)

    cp_prl = cp/prl
    cp_prt = cp/prt

!**************** I direction  (含或不含交叉项) ****************
	do k=1,nk 
    do j=1,nj

		do i=1,ni

			uvwt_line(1,i) = u(i,j,k)
			uvwt_line(2,i) = v(i,j,k)
			uvwt_line(3,i) = w(i,j,k)
			uvwt_line(4,i) = t(i,j,k)

			nxyz_line(1,i) = kcx(i,j,k)
			nxyz_line(2,i) = kcy(i,j,k)
			nxyz_line(3,i) = kcz(i,j,k)

			volp        = vol(i,j,k)

			if(icross == 1) then
			    kxyz_line(1,i) = kcx(i,j,k)
			    kxyz_line(2,i) = etx(i,j,k)
			    kxyz_line(3,i) = ctx(i,j,k)

			    kxyz_line(4,i) = kcy(i,j,k)
			    kxyz_line(5,i) = ety(i,j,k)
			    kxyz_line(6,i) = cty(i,j,k)

			    kxyz_line(7,i) = kcz(i,j,k)
			    kxyz_line(8,i) = etz(i,j,k)
			    kxyz_line(9,i) = ctz(i,j,k)

			    vol_line(i) = volp
			endif


			vslt1 = visl(i,j,k)  
			vslt2 = visl(i,j,k)*cp_prl  
		    do mvist=1,nvis !*湍流时
			    vslt1 = vslt1 + vist(i,j,k)  
			    vslt2 = vslt2 + vist(i,j,k)*cp_prt
			enddo
            vslt1_line(i)  =  vslt1/ volp
            vslt2_line(i)  =  vslt2/ volp
			
			do m=1,4
				mm = m +0
				d1uvwt(m,i) = duvwt(mm,i,j,k)
			enddo
		
		enddo

	    call d2nd_4(4,nmax,ni,uvwt_line,d2uvwt)

		call vis_single(nl,nmax,ni,nxyz_line,uvwt_line,d1uvwt,d2uvwt,vslt1_line,vslt2_line,dfv)

		do i=1,ni 
		   do m=2,nl !*注意，有源项时要m=1,nl
			  dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i) !* vol(i,j,k)
		  enddo
	    enddo

!mml------------------------------
!mml------------------------------
        IF(ICROSS == 1)THEN     !* 是否计算交叉项（i方向）

		    do i=1,ni
			    do m=5,12
				    duvwt_line(m,i) = duvwt(m,i,j,k)
			    enddo
			enddo

            call  DFLUX_VIS_LINE_NEW(0,1,1,nmax,ni,4,uvwt_line,12,duvwt_line,9,kxyz_line &
			                           ,vslt1_line,vslt2_line,dfv,1)         !从1点开始计算

			!*注: 前面3个常数控制是否计算该方向对粘性的贡献
			!*注: 最后常数取1，2或3 ，分别表示i，j或k方向

			do i=1,ni
			do m=2,nl !*注意，有源项时要m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
			enddo
			enddo

		endif !* 结束i方向交叉导数项
!mml------------------------------
!mml------------------------------
	enddo
	enddo


!**************** J direction  (含或不含交叉项) ****************
	do k=1,nk 
    do i=1,ni

		do j=1,nj

			uvwt_line(1,j) = u(i,j,k)
			uvwt_line(2,j) = v(i,j,k)
			uvwt_line(3,j) = w(i,j,k)
			uvwt_line(4,j) = t(i,j,k)

			nxyz_line(1,j) = etx(i,j,k)
			nxyz_line(2,j) = ety(i,j,k)
			nxyz_line(3,j) = etz(i,j,k)

			volp        = vol(i,j,k)

			if(icross == 1) then

			    kxyz_line(1,j) = kcx(i,j,k)
			    kxyz_line(2,j) = etx(i,j,k)
			    kxyz_line(3,j) = ctx(i,j,k)

			    kxyz_line(4,j) = kcy(i,j,k)
			    kxyz_line(5,j) = ety(i,j,k)
			    kxyz_line(6,j) = cty(i,j,k)

			    kxyz_line(7,j) = kcz(i,j,k)
			    kxyz_line(8,j) = etz(i,j,k)
			    kxyz_line(9,j) = ctz(i,j,k)

			    vol_line(j) = volp
			endif

			vslt1 = visl(i,j,k)  
			vslt2 = visl(i,j,k)*cp_prl  
		    do mvist=1,nvis !*湍流时
			    vslt1 = vslt1 + vist(i,j,k)  
			    vslt2 = vslt2 + vist(i,j,k)*cp_prt
			enddo
            vslt1_line(j)  =  vslt1/ volp
            vslt2_line(j)  =  vslt2/ volp

			do m=1,4
			    mm = 4+m
				d1uvwt(m,j) = duvwt(mm,i,j,k)
			enddo
		
		enddo

	    call d2nd_4(4,nmax,nj,uvwt_line,d2uvwt)

		call vis_single(nl,nmax,nj,nxyz_line,uvwt_line,d1uvwt,d2uvwt,vslt1_line,vslt2_line,dfv)

		do j=1,nj 
		do m=2,nl !*注意，有源项时要m=1,nl
			dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j) !* vol(i,j,k)
		enddo
		enddo

!mml------------------------------
!mml------------------------------
        IF(ICROSS == 1)THEN     !* 是否计算交叉项（j方向）

		    do j=1,nj
			    do m=1,4
				    duvwt_line(m,j) = duvwt(m,i,j,k)
			    enddo
			    do m=9,12
				    duvwt_line(m,j) = duvwt(m,i,j,k)
			    enddo
			enddo

            call  DFLUX_VIS_LINE_NEW(1,0,1,nmax,nj,4,uvwt_line,12,duvwt_line,9,kxyz_line &
			                           ,vslt1_line,vslt2_line,dfv,2)         !从1点开始计算
			!*注: 前面3个常数控制是否计算该方向对粘性的贡献
			!*注: 最后常数取1，2或3 ，分别表示i，j或k方向

			do j=1,nj
			do m=2,nl !*注意，有源项时要m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j)
			enddo
			enddo

		endif !* 结束j方向交叉导数项
!mml------------------------------
!mml------------------------------

	enddo
	enddo

!**************** k direction (含或不含交叉项) ****************
	do j=1,nj 
    do i=1,ni
		do k=1,nk

			uvwt_line(1,k) = u(i,j,k)
			uvwt_line(2,k) = v(i,j,k)
			uvwt_line(3,k) = w(i,j,k)
			uvwt_line(4,k) = t(i,j,k)

			nxyz_line(1,k) = ctx(i,j,k)
			nxyz_line(2,k) = cty(i,j,k)
			nxyz_line(3,k) = ctz(i,j,k)

			volp        = vol(i,j,k)

			if(icross == 1) then
			    kxyz_line(1,k) = kcx(i,j,k)
			    kxyz_line(2,k) = etx(i,j,k)
			    kxyz_line(3,k) = ctx(i,j,k)

			    kxyz_line(4,k) = kcy(i,j,k)
			    kxyz_line(5,k) = ety(i,j,k)
			    kxyz_line(6,k) = cty(i,j,k)

			    kxyz_line(7,k) = kcz(i,j,k)
			    kxyz_line(8,k) = etz(i,j,k)
			    kxyz_line(9,k) = ctz(i,j,k)

			    vol_line(k) = volp
			endif

			vslt1 = visl(i,j,k)  
			vslt2 = visl(i,j,k)*cp_prl  
		    do mvist=1,nvis !*湍流时
			    vslt1 = vslt1 + vist(i,j,k)  
			    vslt2 = vslt2 + vist(i,j,k)*cp_prt
			enddo
            vslt1_line(k)  =  vslt1/ volp
            vslt2_line(k)  =  vslt2/ volp

		    do m=1,4
				mm = 8+m
				d1uvwt(m,k) = duvwt(mm,i,j,k)
			enddo

		enddo

	    call d2nd_4(4,nmax,nk,uvwt_line,d2uvwt)

		call vis_single(nl,nmax,nk,nxyz_line,uvwt_line,d1uvwt,d2uvwt,vslt1_line,vslt2_line,dfv)

		do k=1,nk 
		do m=2,nl   !*注意，有源项时要m=1,nl
			dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k) !* vol(i,j,k)
		enddo
		enddo

!mml------------------------------
!mml------------------------------
        IF(ICROSS == 1)THEN     !* 是否计算交叉项（k方向）

		    do k=1,nk
			    do m=1,8
				    duvwt_line(m,k) = duvwt(m,i,j,k)
			    enddo
			enddo

            call  DFLUX_VIS_LINE_NEW(1,1,0,nmax,nk,4,uvwt_line,12,duvwt_line,9,kxyz_line &
			                           ,vslt1_line,vslt2_line,dfv,3)         !从1点开始计算
			!*注: 前面3个常数控制是否计算该方向对粘性的贡献
			!*注: 最后常数取1，2或3 ，分别表示i，j或k方向

			do k=1,nk
			do m=2,nl !*注意，有源项时要m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
			enddo
			enddo

		endif !* 结束k方向交叉导数项
!mml------------------------------
!mml------------------------------
	enddo
	enddo


    deallocate(duvwt)


!
  return
end subroutine VIS_nonconservation
!=============================================================================!

!=============================================================================!
subroutine d1st_4(nl,nmax,ni,uvwt,duvwt)
!-----------------------------------------------------------------------------!
! 4阶中心格式求一阶导数
! 输入：uvwt(nl,nmax)，计算后不能被改写
! 输出：duvwt(nl,nmax)
! 计算范围(NL,1:NI)
!-----------------------------------------------------------------------------!
    use define_precision_mod
	use global_variables,only : nijk2nd
	implicit none
	integer :: nl,ni,nmax,i,m
    real(prec),dimension(nl,nmax) :: uvwt,duvwt
    real(prec) :: cc1,cc2,dc
    real(prec) :: b11,b12,b13,b14,b21,b22,b23,b24,db

	if(ni <= nijk2nd )then  !*TGH. 注意，可能是二维，将为2阶
        do i = 2,ni-1
  	    do m=1,nl
		    duvwt(m,i) = 0.5_prec*(uvwt(m,i+1) - uvwt(m,i-1))
	    enddo
	    enddo

  	    do m=1,nl
		    duvwt(m,1)  = -1.5_prec*uvwt(m,1 ) - 2._prec*uvwt(m,min(2,ni))   + 0.5_prec*uvwt(m,min(3,ni))
		    duvwt(m,ni) = -1.5_prec*uvwt(m,ni) + 2._prec*uvwt(m,max(ni-1,1)) - 0.5_prec*uvwt(m,max(ni-2,1))
		enddo

		return		
	endif

	cc1=8._prec; cc2=-1._prec; dc=12._prec
	b11=-11._prec; b12=18._prec; b13=-9._prec; b14= 2._prec !1边界
	b21= -2._prec; b22=-3._prec; b23= 6._prec; b24=-1._prec !2边界
	db =  6._prec
	do m=1,nl
		duvwt(m,1) = ( b11*uvwt(m,1) + b12*uvwt(m,2) + b13*uvwt(m,3) + b14*uvwt(m,4) )/db
		duvwt(m,2) = ( b21*uvwt(m,1) + b22*uvwt(m,2) + b23*uvwt(m,3) + b24*uvwt(m,4) )/db
		duvwt(m,ni  ) = ( -b11*uvwt(m,ni) - b12*uvwt(m,ni-1) - b13*uvwt(m,ni-2) - b14*uvwt(m,ni-3) )/db
		duvwt(m,ni-1) = ( -b21*uvwt(m,ni) - b22*uvwt(m,ni-1) - b23*uvwt(m,ni-2) - b24*uvwt(m,ni-3) )/db
	enddo

    do i = 3,ni-2
	do m=1,nl
		duvwt(m,i) = ( cc1*(uvwt(m,i+1) - uvwt(m,i-1)) + cc2*(uvwt(m,i+2) - uvwt(m,i-2)) )/dc
	enddo
	enddo

end subroutine d1st_4
!=============================================================================!
!=============================================================================!
subroutine d2ND_4(nl,nmax,ni,uvwt,d2uvwt)
!-----------------------------------------------------------------------------!
! 4阶中心格式求2阶导数
! 输入：uvwt(nl,nmax)，计算后不能被改写
! 输出：d2uvwt(nl,nmax)
! 计算范围(NL,1:NI)
!-----------------------------------------------------------------------------!
	use global_variables,only : nijk2d
    use define_precision_mod
 	implicit none
	integer :: nl,ni,nmax,i,m
    real(prec),dimension(nl,nmax) :: uvwt,d2uvwt
    real(prec) :: c0,c1,c2,dd
	REAL(PREC) :: B11,B12,B13,B14,b15,b21,b22,b23,b24,b25

	if(ni <= nijk2d)then  !*TGH. 注意，可能是二维，将为2阶
        do i = 1,ni
  	    do m=1,nl
		    d2uvwt(m,i) = 0.0
	    enddo
	    enddo

		return		
	endif

	c0=-30._prec; c1=16._prec; c2=-1._prec; dd=12._prec

	b11 =  35._prec/12._prec
	b12 =-104._prec/12._prec
	b13 = 114._prec/12._prec
	b14 = -56._prec/12._prec
	b15 =  11._prec/12._prec

	b21 =  11._prec/12._prec
	b22 = -20._prec/12._prec
	b23 =   6._prec/12._prec
	b24 =   4._prec/12._prec
	b25 =  -1._prec/12._prec

	do m=1,nl
		d2uvwt(m,1) = b11*uvwt(m,1)+b12*uvwt(m,2)+b13*uvwt(m,3)+b14*uvwt(m,4)+b15*uvwt(m,5)
		d2uvwt(m,2) = b21*uvwt(m,1)+b22*uvwt(m,2)+b23*uvwt(m,3)+b24*uvwt(m,4)+b25*uvwt(m,5)
		d2uvwt(m,ni  ) = b11*uvwt(m,ni)+b12*uvwt(m,ni-1)+b13*uvwt(m,ni-2)+b14*uvwt(m,ni-3)+b15*uvwt(m,ni-4)
		d2uvwt(m,ni-1) = b21*uvwt(m,ni)+b22*uvwt(m,ni-1)+b23*uvwt(m,ni-2)+b24*uvwt(m,ni-3)+b25*uvwt(m,ni-4)
	enddo

    do i = 3,ni-2
	do m=1,nl
		d2uvwt(m,i) = ( c0*uvwt(m,i) + c1*( uvwt(m,i+1) + uvwt(m,i-1) ) + c2*( uvwt(m,i+2) + uvwt(m,i-2) ) )/dd
	enddo
	enddo

end subroutine d2ND_4
!=============================================================================!
!=============================================================================!
subroutine vis_single(nl,nmax,ni,nxyz,uvwt,d1uvwt,d2uvwt,vslt1,vslt2,dfv)
!-----------------------------------------------------------------------------!
!* 计算nxyz_line方向上的粘性通量贡献                                          !
!* nxyz -- 除了Jacobian的网格导数                                             !
!* d1uvwt，d2uvwt -- u,v,w,t的一二阶导数                                      !
!* vsjac -- 粘性系数*网格导致/JACOBIAN 后的一阶导数 -- 与粘性应力相关         !
!* tkjac -- (粘性系数/PR)*网格导致/JACOBIAN 后的一阶导数 -- 与热流相关        !
!* fv -- 计算结果                                                             !
!  注意：求得的粘性项贡献对应的时间离散为：(dQ/Dt)/J                          !
!-----------------------------------------------------------------------------!
	use define_precision_mod
	implicit none
	integer :: nl,nmax,ni,i,m
	real(prec) :: nxyz(3,nmax),vslt1(nmax),vslt2(nmax)
	real(prec) :: uvwt(4,nmax),d1uvwt(4,nmax),d2uvwt(4,nmax),dfv(5,nmax)
	real(prec) :: vsn(3,nmax),tkn(3,nmax),dvsn(3,nmax),dtkn(3,nmax)
	real(prec) :: vsp,tkp,dvs,dtk,nxyz_2,fxyz
	real(prec) :: nx,ny,nz,volp,nxv,nyv,nzv,c43,c13,c23
	real(prec) :: taoxx,taoxy,taoxz,taoyy,taoyz,taozz,tem_tao,visp

	c43 = 4._prec/3._prec
	c13 = 1._prec/3._prec
	c23 = 2._prec/3._prec

	do i=1,ni
		nx = nxyz(1,i)
		ny = nxyz(2,i)
		nz = nxyz(3,i)
		
		vsn(1,i) = nx*vslt1(i)
		vsn(2,i) = ny*vslt1(i)
		vsn(3,i) = nz*vslt1(i)
		
		tkn(1,i) = nx*vslt2(i)
		tkn(2,i) = ny*vslt2(i)
		tkn(3,i) = nz*vslt2(i)
	enddo


	call d1st_4(3,nmax,ni,vsn,dvsn)
	call d1st_4(3,nmax,ni,tkn,dtkn)


	do i=1,ni

		nx = nxyz(1,i)
		ny = nxyz(2,i)
		nz = nxyz(3,i)

		visp = vslt1(i)

		tem_tao = c23*( nx*d1uvwt(1,i)+ny*d1uvwt(2,i)+nz*d1uvwt(3,i) )
		
		taoxx = visp * ( 2._prec*nx*d1uvwt(1,i) - tem_tao )
		taoyy = visp * ( 2._prec*ny*d1uvwt(2,i) - tem_tao )
		taozz = visp * ( 2._prec*nz*d1uvwt(3,i) - tem_tao )

		taoxy = visp * ( ny*d1uvwt(1,i) + nx*d1uvwt(2,i) )
		taoxz = visp * ( nz*d1uvwt(1,i) + nx*d1uvwt(3,i) )
		taoyz = visp * ( ny*d1uvwt(3,i) + nz*d1uvwt(2,i) )
		
		dfv(1,i) = 0._prec

		nxyz_2 = nx*nx + ny*ny + nz*nz
        fxyz   = c13*( nx*d2uvwt(1,i) + ny*d2uvwt(2,i) + nz*d2uvwt(3,i) )
		dfv(2,i) = nxyz_2*d2uvwt(1,i) + nx* fxyz
		dfv(3,i) = nxyz_2*d2uvwt(2,i) + ny* fxyz
		dfv(4,i) = nxyz_2*d2uvwt(3,i) + nz* fxyz


		dfv(2,i) = ( c43*nx*nx + ny*ny + nz*nz )*d2uvwt(1,i)  &
		         +   c13*( nx*ny*d2uvwt(2,i) + nx*nz*d2uvwt(3,i) )

		dfv(2,i) = visp*dfv(2,i) + (c43*nx*dvsn(1,i) + ny*dvsn(2,i) + nz*dvsn(3,i) )* d1uvwt(1,i) &
		         + ( ny*dvsn(1,i) - c23*nx*dvsn(2,i) )* d1uvwt(2,i)  &
				 + ( nz*dvsn(1,i) - c23*nx*dvsn(3,i) )* d1uvwt(3,i)

		dfv(3,i) = ( nx*nx + c43*ny*ny + nz*nz )*d2uvwt(2,i)  &
		         +   c13*( nx*ny*d2uvwt(1,i) + ny*nz*d2uvwt(3,i) )
		dfv(3,i) = visp*dfv(3,i) + (nx*dvsn(1,i) + c43*ny*dvsn(2,i) + nz*dvsn(3,i) )* d1uvwt(2,i) &
		         + ( nx*dvsn(2,i) - c23*ny*dvsn(1,i) )* d1uvwt(1,i)  &
				 + ( nz*dvsn(2,i) - c23*ny*dvsn(3,i) )* d1uvwt(3,i)

		dfv(4,i) = ( nx*nx + ny*ny + c43*nz*nz )*d2uvwt(3,i)  &
		         +   c13*( nx*nz*d2uvwt(1,i) + ny*nz*d2uvwt(2,i) )
		dfv(4,i) = visp*dfv(4,i) + (nx*dvsn(1,i) + ny*dvsn(2,i) + c43*nz*dvsn(3,i) )* d1uvwt(3,i) &
		         + ( nx*dvsn(3,i) - c23*nz*dvsn(1,i) )* d1uvwt(1,i)  &
				 + ( ny*dvsn(3,i) - c23*nz*dvsn(2,i) )* d1uvwt(2,i)

		dfv(5,i) = uvwt(1,i)*dfv(2,i) + uvwt(2,i)*dfv(3,i) +uvwt(3,i)*dfv(4,i)  &
		         + (nx*taoxx + ny*taoxy + nz*taoxz)*d1uvwt(1,i)                 &
			 	 + (nx*taoxy + ny*taoyy + nz*taoyz)*d1uvwt(2,i)                 &
				 + (nx*taoxz + ny*taoyz + nz*taozz)*d1uvwt(3,i)
		dfv(5,i) = dfv(5,i)   +    vslt2(i)*(nx*nx+ny*ny+nz*nz) * d2uvwt(4,i)   &
		         + ( nx*dtkn(1,i)+ny*dtkn(2,i)+nz*dtkn(3,i) ) * d1uvwt(4,i)

	enddo

end subroutine vis_single
!=============================================================================!
!=============================================================================!
subroutine DFLUX_VIS_LINE_NEW(r1,r2,r3,nmax,ni,n1,uvwt,n2,duvwt,n3,kxyz, &
                         vslt1,vslt2,dfv,IJK)
!--------------------------------------------------------------------------!
!* 从1点开始计算
!* 先把计算空间的一阶导数duvwt转换成笛卡尔系中的duvwtdxtz
!* 再构造粘性通量
!* r1,r2,r3分别控制是否计算该方向导数对粘性项的贡献
!* IJK = 1,2,3 分别表示计算i，j或k方向上的粘性通量
!* nl = 5 为三维情况下的粘性项
!* 该模块传入的工作数组是 1 : ni
!--------------------------------------------------------------------------!

	implicit none
	
	integer :: r1,r2,r3,IJK,nl,nmax,ni,i,m,m1,m2,m3,n1,n2,n3,n

	real :: uvwt(n1,nmax),duvwt(n2,nmax),kxyz(n3,nmax),dfv(5,nmax)
	real :: vslt1(nmax),vslt2(nmax),fv(5,nmax)
    real :: duvwtdxyz(12),txyz(9),dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	real :: cc,vs,vscc,nx,ny,nz,kcp,temp(12)

    CC=2.0/3.0

	do i=1,ni      !* 从1点开始计

        !------- 计算系下的一阶导数转换到笛卡尔系 ------!
		do m=1,3   !* 分别对应x，y，z方向
			m1=M+M+M-2
			m2=m1+1
			m3=m2+1
			duvwtdxyz(m  ) = 0.0
			duvwtdxyz(m+3) = 0.0
			duvwtdxyz(m+6) = 0.0
			duvwtdxyz(m+9) = 0.0
!			temp(m+9) = 0.0

            do n=1,r1   !m1表示i方向,r1=0剔除该方向对变量梯度的关系
			   duvwtdxyz(m  ) = duvwtdxyz(m  ) + kxyz(m1,i)*duvwt(1 ,i) 
			   duvwtdxyz(m+3) = duvwtdxyz(m+3) + kxyz(m1,i)*duvwt(2 ,i)
			   duvwtdxyz(m+6) = duvwtdxyz(m+6) + kxyz(m1,i)*duvwt(3 ,i) 
			   duvwtdxyz(m+9) = duvwtdxyz(m+9) + kxyz(m1,i)*duvwt(4 ,i) 
		    enddo

            do n=1,r2  !m2表示j方向,r2=0剔除该方向对变量梯度的关系
			   duvwtdxyz(m  ) = duvwtdxyz(m  ) + kxyz(m2,i)*duvwt(5 ,i) 
			   duvwtdxyz(m+3) = duvwtdxyz(m+3) + kxyz(m2,i)*duvwt(6 ,i)
			   duvwtdxyz(m+6) = duvwtdxyz(m+6) + kxyz(m2,i)*duvwt(7 ,i) 
			   duvwtdxyz(m+9) = duvwtdxyz(m+9) + kxyz(m2,i)*duvwt(8 ,i) 
		    enddo

            do n=1,r3  !m3表示k方向,r3=0剔除该方向对变量梯度的关系
			   duvwtdxyz(m  ) = duvwtdxyz(m  ) + kxyz(m3,i)*duvwt(9 ,i) 
			   duvwtdxyz(m+3) = duvwtdxyz(m+3) + kxyz(m3,i)*duvwt(10 ,i)
			   duvwtdxyz(m+6) = duvwtdxyz(m+6) + kxyz(m3,i)*duvwt(11 ,i) 
			   duvwtdxyz(m+9) = duvwtdxyz(m+9) + kxyz(m3,i)*duvwt(12 ,i) 
		    enddo
		enddo
        !------- 完毕. 计算系下的一阶导数转换到笛卡尔系 ------!


		vs   = vslt1(i)

		vscc = vs*CC

		dudx = duvwtdxyz(1)
		dudy = duvwtdxyz(2)
		dudz = duvwtdxyz(3)

		dvdx = duvwtdxyz(4)
		dvdy = duvwtdxyz(5)
		dvdz = duvwtdxyz(6)

		dwdx = duvwtdxyz(7)
		dwdy = duvwtdxyz(8)
		dwdz = duvwtdxyz(9)
		
        txyz(1) = vscc * ( 2.0*dudx - dvdy - dwdz )
        txyz(5) = vscc * ( 2.0*dvdy - dwdz - dudx )
		txyz(9) = vscc * ( 2.0*dwdz - dudx - dvdy )
        txyz(2) = vs * ( dudy + dvdx )
        txyz(3) = vs * ( dudz + dwdx )
        txyz(6) = vs * ( dvdz + dwdy )
        txyz(4) = txyz(2)
        txyz(7) = txyz(3)
        txyz(8) = txyz(6)
		

		nx = kxyz(IJK  ,i)	
		ny = kxyz(IJK+3,i)	
		nz = kxyz(IJK+6,i)	

		kcp = vslt2(i)

		fv(1,i) = 0.0

		fv(2,i) = txyz(1)*nx+txyz(2)*ny+txyz(3)*nz
		fv(3,i) = txyz(4)*nx+txyz(5)*ny+txyz(6)*nz
		fv(4,i) = txyz(7)*nx+txyz(8)*ny+txyz(9)*nz

		fv(5,i) = uvwt(1,I)*fv(2,i)+uvwt(2,i)*fv(3,i)+uvwt(3,i)*fv(4,i)    +  &
							kcp*(duvwtdxyz(10)*nx+duvwtdxyz(11)*ny +  &
							duvwtdxyz(12)*nz )
		
	enddo

	CALL DFLUX_VIS_4th_1toni(5,NMAX,NI,FV,DFV)

	return
end subroutine DFLUX_VIS_LINE_NEW

!=============================================================================!
!=============================================================================!
SUBROUTINE DFLUX_VIS_2nd_1toni(NL,NMAX,NI,FV,DFV)
!------------------------------------------------------------------!
!* 已知结点fv，求结点dfv的一阶导数（计算空间）                     !
!* 暂时用二阶                                                      !
!* 特别注意：考虑到调用该程序时的FV(NL,0:NMAX),DFV(NL,1:NMAX)      !
!------------------------------------------------------------------!
	USE DEFINE_PRECISION_MOD
	use global_variables,only : nijk2d
	IMPLICIT NONE
	INTEGER :: I,NL,NMAX,NI,M
	REAL(PREC) :: FV(NL,NMAX),DFV(NL,NMAX)
	REAL(PREC) :: B1,B2,B3,B4,COE

	COE =  0.5_PREC
	B1  = -1.5_PREC
	B2  =  2.0_PREC
	B3  = -0.5_PREC 
	
	IF(NI <= nijk2d) THEN !*TGH. 注意，可能是二维，将为2阶
         DFV(:,:) = 0.0
		RETURN
	ENDIF

	DO I=2,NI-1
	DO M=1,NL
		DFV(M,I) = COE * ( FV(M,I+1) - FV(M,I-1) )
	ENDDO
	ENDDO

	DO M=1,NL
		DFV(M,1) = B1*FV(M,1) + B2*FV(M,2) + B3*FV(M,3)
		DFV(M,NI) = -B1*FV(M,NI) - B2*FV(M,NI-1) - B3*FV(M,NI-2)
	ENDDO

END SUBROUTINE DFLUX_VIS_2nd_1toni
!=============================================================================!
!=============================================================================!
SUBROUTINE DFLUX_VIS_4TH_1toni(NL,NMAX,NI,FV,DFV)
!------------------------------------------------------------------!
!* 已知结点fv，求结点dfv的一阶导数（计算空间）                     !
!* 暂时用4阶                                                      !
!* 特别注意：考虑到调用该程序时的FV(NL,0:NMAX),DFV(NL,1:NMAX)      !
!------------------------------------------------------------------!
	USE DEFINE_PRECISION_MOD
	use global_variables,only : nijk2nd
	IMPLICIT NONE
	INTEGER :: I,NL,NMAX,NI,M
	REAL(PREC) :: FV(NL,NMAX),DFV(NL,NMAX)
	REAL(PREC) :: B11,B12,B13,B14,b21,b22,b23,b24,coe1,coe2,dd

	COE1 =  8.0_prec/12._prec
	coe2 =   -1._prec/12._prec

	b11 = -11._prec
	b12 =  18._prec
	b13 =   9._prec
	b14 =  -2._prec

	b21 =  -2._prec
	b22 =  -3._prec
	b23 =   6._prec
	b24 =  -1._prec

	dd  = 6.0_prec

	IF( NI <= nijk2nd) THEN !*TGH. 注意，可能是二维，将为2阶
	  DO I=2,NI-1
	  DO M=1,NL
		DFV(M,I) = 0.5_prec * ( FV(M,I+1) - FV(M,I-1) )
	  ENDDO
	  ENDDO

	  DO M=1,NL
		DFV(M,1) = -1.5_prec*FV(M,1) + 2._prec*FV(M,2) - 0.5_prec*FV(M,3)
		DFV(M,NI) = 1.5_prec*FV(M,NI) - 2._prec*FV(M,NI-1) + 0.5_prec*FV(M,NI-2)
	  ENDDO

	  return
	ENDIF

	DO I=3,NI-2
	DO M=1,NL
		DFV(M,I) = COE1 * ( FV(M,I+1) - FV(M,I-1) ) + COE2 * ( FV(M,I+2) - FV(M,I-2) )
	ENDDO
	ENDDO

	DO M=1,NL
		DFV(M,1)    = ( B11*FV(M,1 ) + B12*FV(M,2   ) + B13*FV(M,3   ) + b14*FV(M,4)    )/dd
		DFV(M,NI)   =-( B11*FV(M,NI) + B12*FV(M,NI-1) + B13*FV(M,NI-2) + B14*FV(M,NI-3) )/dd
		DFV(M,2)    = ( B21*FV(M,1 ) + B22*FV(M,2   ) + B23*FV(M,3   ) + b24*FV(M,4)    )/dd
		DFV(M,NI-1) =-( B21*FV(M,NI) + B22*FV(M,NI-1) + B23*FV(M,NI-2) + B24*FV(M,NI-3) )/dd
	ENDDO

END SUBROUTINE DFLUX_VIS_4TH_1toni
!=============================================================================!
!=============================================================================!
subroutine D1UVWT_4ORDER
!-------------------------------------------------------------------------!
!* 求u v w T 的一阶导数（4阶精度）
!* 通过 mod global_variables 和 mod duvwt_all_field 来传递数据
!
!--------------------------------------------------------------------------
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt
  implicit none
	integer :: i,j,k,m
	real    :: uvwt(4,nmax),d1uvwt(4,nmax)

! I direction

	do k=1,nk
	do j=1,nj
			do i=1,ni

				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)

			enddo

	        call d1st_4(4,nmax,ni,uvwt,d1uvwt)

			do i=1,ni
			do m=1,4
				duvwt(m,i,j,k)=d1uvwt(m,i)
			enddo
			enddo

	enddo
	enddo
!
! J direction

	do k=1,nk
	do i=1,ni
			do j=1,nj

				uvwt(1,j) = u(i,j,k)
				uvwt(2,j) = v(i,j,k)
				uvwt(3,j) = w(i,j,k)
				uvwt(4,j) = t(i,j,k)

			enddo

			call d1st_4(4,nmax,nj,uvwt,d1uvwt)

			do j=1,nj
			do m=1,4
				duvwt(m+4,i,j,k)=d1uvwt(m,j)
			enddo
			enddo

	enddo
	enddo

! K direction

	do j=1,nj
	do i=1,ni
			do k=1,nk

				uvwt(1,k) = u(i,j,k)
				uvwt(2,k) = v(i,j,k)
				uvwt(3,k) = w(i,j,k)
				uvwt(4,k) = t(i,j,k)

			enddo

			call d1st_4(4,nmax,nk,uvwt,d1uvwt)

			do k=1,nk
			do m=1,4

					duvwt(m+8,i,j,k)=d1uvwt(m,k)

			enddo
			enddo

	enddo
	enddo
!
  return
end subroutine D1UVWT_4ORDER
!=============================================================================!
!=============================================================================!
!=============================================================================!
subroutine VIS_nonconservation_virtual
  use define_precision_mod
  use global_variables,only : u,v,w,t,reynolds,visl,vist,dq,nvis     &
                            , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol &
							, nl,nmax,ni,nj,nk,gama,moo,prl,prt
  use duvwt_all_field,only : duvwt
  implicit none
!-----------------------------------------------------------------------------!
!   要用到一个虚点值                                                          !
!                                                                             !
!  把同方向导数和交叉导数分离开，采用不同格式的求解方法                       !
!  注意：求得的粘性项贡献对应的时间离散为：(dQ/Dt)/J                          !
!                                                                             !
!                                                                             !
! 特别注意：数组的起始坐标                                                    !
!                                                                             !
!  比WCNSE5_VIS_FLUX （1）少耗时3.7% —— 二阶中心计算交叉(还需要测试正确性） !
!                    （2）少耗时42%  —— 不计算交叉                           !
!                                                                             !
!  设计：毛枚良，涂国华                                                       !
!  调试：涂国华 2009.03.20                                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
	real(prec) :: duvwtdxyz(12,nmax),nxyz_line(3,nmax)
	real(prec) :: uvwt_line(4,0:nmax+1),duvwt_line(12,nmax)
	real(prec) :: d1uvwt(4,nmax),d2uvwt(4,nmax) !tgh
	real(prec) :: kxyz_line(9,nmax) !,kxyz_half(9,0:nmax)
	real(prec) :: vol_line(nmax),fv(nl,nmax),dfv(nl,nmax)
	real(prec) :: vslt1_line(nmax),vslt2_line(nmax)
	real(prec) :: cp,cp_prl,cp_prt,re,vslt1,vslt2
	integer :: mvist,ICROSS,mm
	real(prec) :: volp !,vslt1p,vslt2p
	real(prec) :: nx,ny,nz !,nxyz_line(3,nmax) !tgh
!	real(prec) :: d1uvwt(4,nmax),d2uvwt(4,nmax)

	ICROSS = 1 !特别注意，是否计算交叉导数

	allocate(duvwt(12,ni,nj,nk))

    call D1UVWT_4ORDER_virtual    !*TGH. 直接求结点1阶导数



	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)

    cp_prl = cp/prl
    cp_prt = cp/prt

!**************** I direction  (含或不含交叉项) ****************
	do k=1,nk 
    do j=1,nj

		do i=1,ni

			uvwt_line(1,i) = u(i,j,k)
			uvwt_line(2,i) = v(i,j,k)
			uvwt_line(3,i) = w(i,j,k)
			uvwt_line(4,i) = t(i,j,k)

			nxyz_line(1,i) = kcx(i,j,k)
			nxyz_line(2,i) = kcy(i,j,k)
			nxyz_line(3,i) = kcz(i,j,k)

			volp        = vol(i,j,k)

			if(icross == 1) then
			    kxyz_line(1,i) = kcx(i,j,k)
			    kxyz_line(2,i) = etx(i,j,k)
			    kxyz_line(3,i) = ctx(i,j,k)

			    kxyz_line(4,i) = kcy(i,j,k)
			    kxyz_line(5,i) = ety(i,j,k)
			    kxyz_line(6,i) = cty(i,j,k)

			    kxyz_line(7,i) = kcz(i,j,k)
			    kxyz_line(8,i) = etz(i,j,k)
			    kxyz_line(9,i) = ctz(i,j,k)

			    vol_line(i) = volp
			endif


			vslt1 = visl(i,j,k)  
			vslt2 = visl(i,j,k)*cp_prl  
		    do mvist=1,nvis !*湍流时
			    vslt1 = vslt1 + vist(i,j,k)  
			    vslt2 = vslt2 + vist(i,j,k)*cp_prt
			enddo
            vslt1_line(i)  =  vslt1/ volp
            vslt2_line(i)  =  vslt2/ volp
			
			do m=1,4
				mm = m +0
				d1uvwt(m,i) = duvwt(mm,i,j,k)
			enddo
		
		enddo

		i=0
		uvwt_line(1,i) = u(i,j,k)
		uvwt_line(2,i) = v(i,j,k)
		uvwt_line(3,i) = w(i,j,k)
		uvwt_line(4,i) = t(i,j,k)
		i=ni+1
		uvwt_line(1,i) = u(i,j,k)
		uvwt_line(2,i) = v(i,j,k)
		uvwt_line(3,i) = w(i,j,k)
		uvwt_line(4,i) = t(i,j,k)

	    call d2nd_4_virtual(4,nmax,ni,uvwt_line,d2uvwt)

		call vis_single_vir(nl,nmax,ni,nxyz_line,uvwt_line,d1uvwt,d2uvwt,vslt1_line,vslt2_line,dfv)

		do i=1,ni 
		   do m=2,nl !*注意，有源项时要m=1,nl
			  dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i) !* vol(i,j,k)
		  enddo
	    enddo

!mml------------------------------
!mml------------------------------
        IF(ICROSS == 1)THEN     !* 是否计算交叉项（i方向）

		    do i=1,ni
			    do m=5,12
				    duvwt_line(m,i) = duvwt(m,i,j,k)
			    enddo
			enddo

            call  DFLUX_VIS_LINE_NEW_vir(0,1,1,nmax,ni,4,uvwt_line,12,duvwt_line,9,kxyz_line &
			                           ,vslt1_line,vslt2_line,dfv,1)         !从1点开始计算

			!*注: 前面3个常数控制是否计算该方向对粘性的贡献
			!*注: 最后常数取1，2或3 ，分别表示i，j或k方向

			do i=1,ni
			do m=2,nl !*注意，有源项时要m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
			enddo
			enddo

		endif !* 结束i方向交叉导数项
!mml------------------------------
!mml------------------------------
	enddo
	enddo


!**************** J direction  (含或不含交叉项) ****************
	do k=1,nk 
    do i=1,ni

		do j=1,nj

			uvwt_line(1,j) = u(i,j,k)
			uvwt_line(2,j) = v(i,j,k)
			uvwt_line(3,j) = w(i,j,k)
			uvwt_line(4,j) = t(i,j,k)

			nxyz_line(1,j) = etx(i,j,k)
			nxyz_line(2,j) = ety(i,j,k)
			nxyz_line(3,j) = etz(i,j,k)

			volp        = vol(i,j,k)

			if(icross == 1) then

			    kxyz_line(1,j) = kcx(i,j,k)
			    kxyz_line(2,j) = etx(i,j,k)
			    kxyz_line(3,j) = ctx(i,j,k)

			    kxyz_line(4,j) = kcy(i,j,k)
			    kxyz_line(5,j) = ety(i,j,k)
			    kxyz_line(6,j) = cty(i,j,k)

			    kxyz_line(7,j) = kcz(i,j,k)
			    kxyz_line(8,j) = etz(i,j,k)
			    kxyz_line(9,j) = ctz(i,j,k)

			    vol_line(j) = volp
			endif

			vslt1 = visl(i,j,k)  
			vslt2 = visl(i,j,k)*cp_prl  
		    do mvist=1,nvis !*湍流时
			    vslt1 = vslt1 + vist(i,j,k)  
			    vslt2 = vslt2 + vist(i,j,k)*cp_prt
			enddo
            vslt1_line(j)  =  vslt1/ volp
            vslt2_line(j)  =  vslt2/ volp

			do m=1,4
			    mm = 4+m
				d1uvwt(m,j) = duvwt(mm,i,j,k)
			enddo
		
		enddo

		j=0
		uvwt_line(1,j) = u(i,j,k)
		uvwt_line(2,j) = v(i,j,k)
		uvwt_line(3,j) = w(i,j,k)
		uvwt_line(4,j) = t(i,j,k)
		j=nj+1
		uvwt_line(1,j) = u(i,j,k)
		uvwt_line(2,j) = v(i,j,k)
		uvwt_line(3,j) = w(i,j,k)
		uvwt_line(4,j) = t(i,j,k)

	    call d2nd_4_virtual(4,nmax,nj,uvwt_line,d2uvwt)

		call vis_single_vir(nl,nmax,nj,nxyz_line,uvwt_line,d1uvwt,d2uvwt,vslt1_line,vslt2_line,dfv)

		do j=1,nj 
		do m=2,nl !*注意，有源项时要m=1,nl
			dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j) !* vol(i,j,k)
		enddo
		enddo

!mml------------------------------
!mml------------------------------
        IF(ICROSS == 1)THEN     !* 是否计算交叉项（j方向）

		    do j=1,nj
			    do m=1,4
				    duvwt_line(m,j) = duvwt(m,i,j,k)
			    enddo
			    do m=9,12
				    duvwt_line(m,j) = duvwt(m,i,j,k)
			    enddo
			enddo

            call  DFLUX_VIS_LINE_NEW_vir(1,0,1,nmax,nj,4,uvwt_line,12,duvwt_line,9,kxyz_line &
			                           ,vslt1_line,vslt2_line,dfv,2)         !从1点开始计算
			!*注: 前面3个常数控制是否计算该方向对粘性的贡献
			!*注: 最后常数取1，2或3 ，分别表示i，j或k方向

			do j=1,nj
			do m=2,nl !*注意，有源项时要m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j)
			enddo
			enddo

		endif !* 结束j方向交叉导数项
!mml------------------------------
!mml------------------------------

	enddo
	enddo

!**************** k direction (含或不含交叉项) ****************
	do j=1,nj 
    do i=1,ni
		do k=1,nk

			uvwt_line(1,k) = u(i,j,k)
			uvwt_line(2,k) = v(i,j,k)
			uvwt_line(3,k) = w(i,j,k)
			uvwt_line(4,k) = t(i,j,k)

			nxyz_line(1,k) = ctx(i,j,k)
			nxyz_line(2,k) = cty(i,j,k)
			nxyz_line(3,k) = ctz(i,j,k)

			volp        = vol(i,j,k)

			if(icross == 1) then
			    kxyz_line(1,k) = kcx(i,j,k)
			    kxyz_line(2,k) = etx(i,j,k)
			    kxyz_line(3,k) = ctx(i,j,k)

			    kxyz_line(4,k) = kcy(i,j,k)
			    kxyz_line(5,k) = ety(i,j,k)
			    kxyz_line(6,k) = cty(i,j,k)

			    kxyz_line(7,k) = kcz(i,j,k)
			    kxyz_line(8,k) = etz(i,j,k)
			    kxyz_line(9,k) = ctz(i,j,k)

			    vol_line(k) = volp
			endif

			vslt1 = visl(i,j,k)  
			vslt2 = visl(i,j,k)*cp_prl  
		    do mvist=1,nvis !*湍流时
			    vslt1 = vslt1 + vist(i,j,k)  
			    vslt2 = vslt2 + vist(i,j,k)*cp_prt
			enddo
            vslt1_line(k)  =  vslt1/ volp
            vslt2_line(k)  =  vslt2/ volp

		    do m=1,4
				mm = 8+m
				d1uvwt(m,k) = duvwt(mm,i,j,k)
			enddo

		enddo

		k=0
		uvwt_line(1,k) = u(i,j,k)
		uvwt_line(2,k) = v(i,j,k)
		uvwt_line(3,k) = w(i,j,k)
		uvwt_line(4,k) = t(i,j,k)
		k=nk+1
		uvwt_line(1,k) = u(i,j,k)
		uvwt_line(2,k) = v(i,j,k)
		uvwt_line(3,k) = w(i,j,k)
		uvwt_line(4,k) = t(i,j,k)
	    call d2nd_4_virtual(4,nmax,nk,uvwt_line,d2uvwt)

		call vis_single_vir(nl,nmax,nk,nxyz_line,uvwt_line,d1uvwt,d2uvwt,vslt1_line,vslt2_line,dfv)

		do k=1,nk 
		do m=2,nl   !*注意，有源项时要m=1,nl
			dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k) !* vol(i,j,k)
		enddo
		enddo

!mml------------------------------
!mml------------------------------
        IF(ICROSS == 1)THEN     !* 是否计算交叉项（k方向）

		    do k=1,nk
			    do m=1,8
				    duvwt_line(m,k) = duvwt(m,i,j,k)
			    enddo
			enddo

            call  DFLUX_VIS_LINE_NEW_vir(1,1,0,nmax,nk,4,uvwt_line,12,duvwt_line,9,kxyz_line &
			                           ,vslt1_line,vslt2_line,dfv,3)         !从1点开始计算
			!*注: 前面3个常数控制是否计算该方向对粘性的贡献
			!*注: 最后常数取1，2或3 ，分别表示i，j或k方向

			do k=1,nk
			do m=2,nl !*注意，有源项时要m=1,nl
				dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
			enddo
			enddo

		endif !* 结束k方向交叉导数项
!mml------------------------------
!mml------------------------------
	enddo
	enddo

    deallocate(duvwt)
!
  return
end subroutine VIS_nonconservation_virtual
!=============================================================================!

!=============================================================================!
subroutine d1st_4_virtual(nl,nmax,ni,uvwt,duvwt)
!-----------------------------------------------------------------------------!
! 4阶中心格式求一阶导数
! 输入：uvwt(nl,0:nmax+1)，计算后不能被改写
! 输出：duvwt(nl,nmax)
! 计算范围(NL,1:NI)
!-----------------------------------------------------------------------------!
    use define_precision_mod
	use global_variables,only : nijk2nd
	implicit none
	integer :: nl,ni,nmax,i,m
    real(prec) :: uvwt(nl,0:nmax+1),duvwt(nl,nmax)
    real(prec) :: cc1,cc2,dc
    real(prec) :: b11,b12,b13,b14,b21,b22,b23,b24,db

	if( ni <= nijk2nd )then  !*TGH. 注意，可能是二维，将为2阶
        do i = 1,ni
  	    do m=1,nl
		    duvwt(m,i) = (uvwt(m,ni+1) - uvwt(m,ni-1) )*0.5_PREC
	    enddo
	    enddo

		return		
	endif

	cc1=8._prec; cc2=-1._prec; dc=12._prec
	b11=-11._prec; b12=18._prec; b13=-9._prec; b14= 2._prec !1边界
	b21= -2._prec; b22=-3._prec; b23= 6._prec; b24=-1._prec !2边界
	db =  6._prec
	do m=1,nl
!		duvwt(m,1) = ( b11*uvwt(m,1) + b12*uvwt(m,2) + b13*uvwt(m,3) + b14*uvwt(m,4) )/db
		duvwt(m,1) = ( b21*uvwt(m,0) + b22*uvwt(m,1) + b23*uvwt(m,2) + b24*uvwt(m,3) )/db
		duvwt(m,2) = ( b21*uvwt(m,1) + b22*uvwt(m,2) + b23*uvwt(m,3) + b24*uvwt(m,4) )/db
!		duvwt(m,ni  ) = ( -b11*uvwt(m,ni  ) - b12*uvwt(m,ni-1) - b13*uvwt(m,ni-2) - b14*uvwt(m,ni-3) )/db
		duvwt(m,ni  ) = ( -b21*uvwt(m,ni+1) - b22*uvwt(m,ni  ) - b23*uvwt(m,ni-1) - b24*uvwt(m,ni-2) )/db
		duvwt(m,ni-1) = ( -b21*uvwt(m,ni  ) - b22*uvwt(m,ni-1) - b23*uvwt(m,ni-2) - b24*uvwt(m,ni-3) )/db
	enddo

    do i = 3,ni-2
	do m=1,nl
		duvwt(m,i) = ( cc1*(uvwt(m,i+1) - uvwt(m,i-1)) + cc2*(uvwt(m,i+2) - uvwt(m,i-2)) )/dc
	enddo
	enddo

end subroutine d1st_4_virtual
!=============================================================================!
!=============================================================================!
subroutine d2ND_4_virtual(nl,nmax,ni,uvwt,d2uvwt)
!-----------------------------------------------------------------------------!
! 4阶中心格式求2阶导数
! 输入：uvwt(nl,0:nmax+1)，计算后不能被改写
! 输出：d2uvwt(nl,nmax)
! 计算范围(NL,1:NI)
!-----------------------------------------------------------------------------!
    use define_precision_mod
	use global_variables,only : nijk2d
 	implicit none
	integer :: nl,ni,nmax,i,m
    real(prec) :: uvwt(nl,0:nmax+1),d2uvwt(nl,nmax)
    real(prec) :: c0,c1,c2,dd
	REAL(PREC) :: B11,B12,B13,B14,b15,b21,b22,b23,b24,b25


	if(ni <= nijk2d)then  !*TGH. 注意，可能是二维，将为2阶
        do i = 1,ni
  	    do m=1,nl
		  d2uvwt(m,ni  ) = 0._prec
	    enddo
	    enddo

		return		
	endif

	c0=-30._prec; c1=16._prec; c2=-1._prec; dd=12._prec

!	b11 =  35._prec/12._prec
!	b12 =-104._prec/12._prec
!	b13 = 114._prec/12._prec
!	b14 = -56._prec/12._prec
!	b15 =  11._prec/12._prec !1边界

	b21 =  11._prec/12._prec
	b22 = -20._prec/12._prec
	b23 =   6._prec/12._prec
	b24 =   4._prec/12._prec
	b25 =  -1._prec/12._prec !2边界



	do m=1,nl
!		d2uvwt(m,1) = b11*uvwt(m,1)+b12*uvwt(m,2)+b13*uvwt(m,3)+b14*uvwt(m,4)+b15*uvwt(m,5)
		d2uvwt(m,1) = b21*uvwt(m,0)+b22*uvwt(m,1)+b23*uvwt(m,2)+b24*uvwt(m,3)+b25*uvwt(m,4)
		d2uvwt(m,2) = b21*uvwt(m,1)+b22*uvwt(m,2)+b23*uvwt(m,3)+b24*uvwt(m,4)+b25*uvwt(m,5)
!		d2uvwt(m,ni  ) = b11*uvwt(m,ni  )+b12*uvwt(m,ni-1)+b13*uvwt(m,ni-2)+b14*uvwt(m,ni-3)+b15*uvwt(m,ni-4)
		d2uvwt(m,ni  ) = b21*uvwt(m,ni+1)+b22*uvwt(m,ni  )+b23*uvwt(m,ni-1)+b24*uvwt(m,ni-2)+b25*uvwt(m,ni-3)
		d2uvwt(m,ni-1) = b21*uvwt(m,ni  )+b22*uvwt(m,ni-1)+b23*uvwt(m,ni-2)+b24*uvwt(m,ni-3)+b25*uvwt(m,ni-4)
	enddo

    do i = 3,ni-2
	do m=1,nl
		d2uvwt(m,i) = ( c0*uvwt(m,i) + c1*( uvwt(m,i+1) + uvwt(m,i-1) ) + c2*( uvwt(m,i+2) + uvwt(m,i-2) ) )/dd
	enddo
	enddo

end subroutine d2ND_4_virtual

!=============================================================================!
!=============================================================================!
!=============================================================================!

subroutine DFLUX_VIS_LINE_NEW_vir(r1,r2,r3,nmax,ni,n1,uvwt,n2,duvwt,n3,kxyz, &
                         vslt1,vslt2,dfv,IJK)
!--------------------------------------------------------------------------!
!* 从1点开始计算
!* 先把计算空间的一阶导数duvwt转换成笛卡尔系中的duvwtdxtz
!* 再构造粘性通量
!* r1,r2,r3分别控制是否计算该方向导数对粘性项的贡献
!* IJK = 1,2,3 分别表示计算i，j或k方向上的粘性通量
!* nl = 5 为三维情况下的粘性项
!* 该模块传入的工作数组是 1 : ni
!--------------------------------------------------------------------------!

	implicit none
	
	integer :: r1,r2,r3,IJK,nl,nmax,ni,i,m,m1,m2,m3,n1,n2,n3,n

	real :: uvwt(n1,0:nmax+1),duvwt(n2,nmax),kxyz(n3,nmax),dfv(5,nmax)
	real :: vslt1(nmax),vslt2(nmax),fv(5,nmax)
    real :: duvwtdxyz(12),txyz(9),dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	real :: cc,vs,vscc,nx,ny,nz,kcp,temp(12)

    CC=2.0/3.0

	do i=1,ni      !* 从1点开始计

        !------- 计算系下的一阶导数转换到笛卡尔系 ------!
		do m=1,3   !* 分别对应x，y，z方向
			m1=M+M+M-2
			m2=m1+1
			m3=m2+1
			duvwtdxyz(m  ) = 0.0
			duvwtdxyz(m+3) = 0.0
			duvwtdxyz(m+6) = 0.0
			duvwtdxyz(m+9) = 0.0
!			temp(m+9) = 0.0

            do n=1,r1   !m1表示i方向,r1=0剔除该方向对变量梯度的关系
			   duvwtdxyz(m  ) = duvwtdxyz(m  ) + kxyz(m1,i)*duvwt(1 ,i) 
			   duvwtdxyz(m+3) = duvwtdxyz(m+3) + kxyz(m1,i)*duvwt(2 ,i)
			   duvwtdxyz(m+6) = duvwtdxyz(m+6) + kxyz(m1,i)*duvwt(3 ,i) 
			   duvwtdxyz(m+9) = duvwtdxyz(m+9) + kxyz(m1,i)*duvwt(4 ,i) 
		    enddo

            do n=1,r2  !m2表示j方向,r2=0剔除该方向对变量梯度的关系
			   duvwtdxyz(m  ) = duvwtdxyz(m  ) + kxyz(m2,i)*duvwt(5 ,i) 
			   duvwtdxyz(m+3) = duvwtdxyz(m+3) + kxyz(m2,i)*duvwt(6 ,i)
			   duvwtdxyz(m+6) = duvwtdxyz(m+6) + kxyz(m2,i)*duvwt(7 ,i) 
			   duvwtdxyz(m+9) = duvwtdxyz(m+9) + kxyz(m2,i)*duvwt(8 ,i) 
		    enddo

            do n=1,r3  !m3表示k方向,r3=0剔除该方向对变量梯度的关系
			   duvwtdxyz(m  ) = duvwtdxyz(m  ) + kxyz(m3,i)*duvwt(9 ,i) 
			   duvwtdxyz(m+3) = duvwtdxyz(m+3) + kxyz(m3,i)*duvwt(10 ,i)
			   duvwtdxyz(m+6) = duvwtdxyz(m+6) + kxyz(m3,i)*duvwt(11 ,i) 
			   duvwtdxyz(m+9) = duvwtdxyz(m+9) + kxyz(m3,i)*duvwt(12 ,i) 
		    enddo
		enddo
        !------- 完毕. 计算系下的一阶导数转换到笛卡尔系 ------!


		vs   = vslt1(i)

		vscc = vs*CC

		dudx = duvwtdxyz(1)
		dudy = duvwtdxyz(2)
		dudz = duvwtdxyz(3)

		dvdx = duvwtdxyz(4)
		dvdy = duvwtdxyz(5)
		dvdz = duvwtdxyz(6)

		dwdx = duvwtdxyz(7)
		dwdy = duvwtdxyz(8)
		dwdz = duvwtdxyz(9)
		
        txyz(1) = vscc * ( 2.0*dudx - dvdy - dwdz )
        txyz(5) = vscc * ( 2.0*dvdy - dwdz - dudx )
		txyz(9) = vscc * ( 2.0*dwdz - dudx - dvdy )
        txyz(2) = vs * ( dudy + dvdx )
        txyz(3) = vs * ( dudz + dwdx )
        txyz(6) = vs * ( dvdz + dwdy )
        txyz(4) = txyz(2)
        txyz(7) = txyz(3)
        txyz(8) = txyz(6)
		

		nx = kxyz(IJK  ,i)	
		ny = kxyz(IJK+3,i)	
		nz = kxyz(IJK+6,i)	

		kcp = vslt2(i)

		fv(1,i) = 0.0

		fv(2,i) = txyz(1)*nx+txyz(2)*ny+txyz(3)*nz
		fv(3,i) = txyz(4)*nx+txyz(5)*ny+txyz(6)*nz
		fv(4,i) = txyz(7)*nx+txyz(8)*ny+txyz(9)*nz

		fv(5,i) = uvwt(1,I)*fv(2,i)+uvwt(2,i)*fv(3,i)+uvwt(3,i)*fv(4,i)    +  &
							kcp*(duvwtdxyz(10)*nx+duvwtdxyz(11)*ny +  &
							duvwtdxyz(12)*nz )
		
	enddo

	CALL DFLUX_VIS_4th_1toni(5,NMAX,NI,FV,DFV)

	return
end subroutine DFLUX_VIS_LINE_NEW_vir

!=============================================================================!
!=============================================================================!
subroutine D1UVWT_4ORDER_virtual
!-------------------------------------------------------------------------!
!* 求u v w T 的一阶导数（4阶精度）
!* 通过 mod global_variables 和 mod duvwt_all_field 来传递数据
!
!--------------------------------------------------------------------------
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt
  implicit none
	integer :: i,j,k,m
	real    :: uvwt(4,0:nmax+1),d1uvwt(4,nmax)

! I direction

	do k=1,nk
	do j=1,nj
			do i=0,ni+1

				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)

			enddo

	        call d1st_4_virtual(4,nmax,ni,uvwt,d1uvwt)

			do i=1,ni
			do m=1,4
				duvwt(m,i,j,k)=d1uvwt(m,i)
			enddo
			enddo

	enddo
	enddo
!
! J direction

	do k=1,nk
	do i=1,ni
			do j=0,nj+1

				uvwt(1,j) = u(i,j,k)
				uvwt(2,j) = v(i,j,k)
				uvwt(3,j) = w(i,j,k)
				uvwt(4,j) = t(i,j,k)

			enddo

			call d1st_4_virtual(4,nmax,nj,uvwt,d1uvwt)

			do j=1,nj
			do m=1,4
				duvwt(m+4,i,j,k)=d1uvwt(m,j)
			enddo
			enddo

	enddo
	enddo

! K direction

	do j=1,nj
	do i=1,ni
			do k=0,nk+1

				uvwt(1,k) = u(i,j,k)
				uvwt(2,k) = v(i,j,k)
				uvwt(3,k) = w(i,j,k)
				uvwt(4,k) = t(i,j,k)

			enddo

			call d1st_4_virtual(4,nmax,nk,uvwt,d1uvwt)

			do k=1,nk
			do m=1,4

					duvwt(m+8,i,j,k)=d1uvwt(m,k)

			enddo
			enddo

	enddo
	enddo
!
  return
end subroutine D1UVWT_4ORDER_virtual
!=============================================================================!
!=============================================================================!
!=============================================================================!
!=============================================================================!
subroutine vis_single_vir(nl,nmax,ni,nxyz,uvwt,d1uvwt,d2uvwt,vslt1,vslt2,dfv)
!-----------------------------------------------------------------------------!
!* 计算nxyz_line方向上的粘性通量贡献                                          !
!* nxyz -- 除了Jacobian的网格导数                                             !
!* d1uvwt，d2uvwt -- u,v,w,t的一二阶导数                                      !
!* vsjac -- 粘性系数*网格导致/JACOBIAN 后的一阶导数 -- 与粘性应力相关         !
!* tkjac -- (粘性系数/PR)*网格导致/JACOBIAN 后的一阶导数 -- 与热流相关        !
!* fv -- 计算结果                                                             !
!  注意：求得的粘性项贡献对应的时间离散为：(dQ/Dt)/J                          !
!-----------------------------------------------------------------------------!
	use define_precision_mod
	implicit none
	integer :: nl,nmax,ni,i,m
	real(prec) :: nxyz(3,nmax),vslt1(nmax),vslt2(nmax)
	real(prec) :: uvwt(4,0:nmax+1),d1uvwt(4,nmax),d2uvwt(4,nmax),dfv(5,nmax)
	real(prec) :: vsn(3,nmax),tkn(3,nmax),dvsn(3,nmax),dtkn(3,nmax)
	real(prec) :: vsp,tkp,dvs,dtk,nxyz_2,fxyz
	real(prec) :: nx,ny,nz,volp,nxv,nyv,nzv,c43,c13,c23
	real(prec) :: taoxx,taoxy,taoxz,taoyy,taoyz,taozz,tem_tao,visp

	c43 = 4._prec/3._prec
	c13 = 1._prec/3._prec
	c23 = 2._prec/3._prec

	do i=1,ni
		nx = nxyz(1,i)
		ny = nxyz(2,i)
		nz = nxyz(3,i)
		
		vsn(1,i) = nx*vslt1(i)
		vsn(2,i) = ny*vslt1(i)
		vsn(3,i) = nz*vslt1(i)
		
		tkn(1,i) = nx*vslt2(i)
		tkn(2,i) = ny*vslt2(i)
		tkn(3,i) = nz*vslt2(i)
	enddo


	call d1st_4(3,nmax,ni,vsn,dvsn)
	call d1st_4(3,nmax,ni,tkn,dtkn)


	do i=1,ni

		nx = nxyz(1,i)
		ny = nxyz(2,i)
		nz = nxyz(3,i)

		visp = vslt1(i)

		tem_tao = c23*( nx*d1uvwt(1,i)+ny*d1uvwt(2,i)+nz*d1uvwt(3,i) )
		
		taoxx = visp * ( 2._prec*nx*d1uvwt(1,i) - tem_tao )
		taoyy = visp * ( 2._prec*ny*d1uvwt(2,i) - tem_tao )
		taozz = visp * ( 2._prec*nz*d1uvwt(3,i) - tem_tao )

		taoxy = visp * ( ny*d1uvwt(1,i) + nx*d1uvwt(2,i) )
		taoxz = visp * ( nz*d1uvwt(1,i) + nx*d1uvwt(3,i) )
		taoyz = visp * ( ny*d1uvwt(3,i) + nz*d1uvwt(2,i) )
		
		dfv(1,i) = 0._prec

		nxyz_2 = nx*nx + ny*ny + nz*nz
        fxyz   = c13*( nx*d2uvwt(1,i) + ny*d2uvwt(2,i) + nz*d2uvwt(3,i) )
		dfv(2,i) = nxyz_2*d2uvwt(1,i) + nx* fxyz
		dfv(3,i) = nxyz_2*d2uvwt(2,i) + ny* fxyz
		dfv(4,i) = nxyz_2*d2uvwt(3,i) + nz* fxyz


		dfv(2,i) = ( c43*nx*nx + ny*ny + nz*nz )*d2uvwt(1,i)  &
		         +   c13*( nx*ny*d2uvwt(2,i) + nx*nz*d2uvwt(3,i) )

		dfv(2,i) = visp*dfv(2,i) + (c43*nx*dvsn(1,i) + ny*dvsn(2,i) + nz*dvsn(3,i) )* d1uvwt(1,i) &
		         + ( ny*dvsn(1,i) - c23*nx*dvsn(2,i) )* d1uvwt(2,i)  &
				 + ( nz*dvsn(1,i) - c23*nx*dvsn(3,i) )* d1uvwt(3,i)

		dfv(3,i) = ( nx*nx + c43*ny*ny + nz*nz )*d2uvwt(2,i)  &
		         +   c13*( nx*ny*d2uvwt(1,i) + ny*nz*d2uvwt(3,i) )
		dfv(3,i) = visp*dfv(3,i) + (nx*dvsn(1,i) + c43*ny*dvsn(2,i) + nz*dvsn(3,i) )* d1uvwt(2,i) &
		         + ( nx*dvsn(2,i) - c23*ny*dvsn(1,i) )* d1uvwt(1,i)  &
				 + ( nz*dvsn(2,i) - c23*ny*dvsn(3,i) )* d1uvwt(3,i)

		dfv(4,i) = ( nx*nx + ny*ny + c43*nz*nz )*d2uvwt(3,i)  &
		         +   c13*( nx*nz*d2uvwt(1,i) + ny*nz*d2uvwt(2,i) )
		dfv(4,i) = visp*dfv(4,i) + (nx*dvsn(1,i) + ny*dvsn(2,i) + c43*nz*dvsn(3,i) )* d1uvwt(3,i) &
		         + ( nx*dvsn(3,i) - c23*nz*dvsn(1,i) )* d1uvwt(1,i)  &
				 + ( ny*dvsn(3,i) - c23*nz*dvsn(2,i) )* d1uvwt(2,i)

		dfv(5,i) = uvwt(1,i)*dfv(2,i) + uvwt(2,i)*dfv(3,i) +uvwt(3,i)*dfv(4,i)  &
		         + (nx*taoxx + ny*taoxy + nz*taoxz)*d1uvwt(1,i)                 &
			 	 + (nx*taoxy + ny*taoyy + nz*taoyz)*d1uvwt(2,i)                 &
				 + (nx*taoxz + ny*taoyz + nz*taozz)*d1uvwt(3,i)
		dfv(5,i) = dfv(5,i)   +    vslt2(i)*(nx*nx+ny*ny+nz*nz) * d2uvwt(4,i)   &
		         + ( nx*dtkn(1,i)+ny*dtkn(2,i)+nz*dtkn(3,i) ) * d1uvwt(4,i)

	enddo

end subroutine vis_single_vir
!=============================================================================!
!=============================================================================!
