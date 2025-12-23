!=============================================================================!
subroutine WCNS_E_5(cc1,cc2,ni,nmax,nl,method,cmethd,flux_type, &
                    limiter,efix,trxyz,q,df)
    use define_precision_mod
    use global_const,only:small,rmin_limit,pmin_limit
	implicit none
!-------------------------------------------------------------------------------!
!	Function                                                                    !
!		Subroutine WCNS_E_5 is used to calculate the convect flux by the use of !
!		a type of fifth order accurate weighted compact nonlinear difference    !
!		schemes(WCNS-E-5).                                                      !
!	Called by                                                                   !
!			Subroutine inviscd3d in Convect.f90                                 !
!	Editor                                                                      !
!			TU Guohua                                                           !
!	Date                                                                        !
!			28-Oct-2009                                                         !
!                                                                               !
!   ������Ҫ�õ�3�����ֵ                                                       !
!   �����߽磨����3����㣩�ϵ��ܶ�Ϊ������ʾ�ö˵����ֵ��������������Բ�ֵ !
!                                          �ö˱�Ϊ�߽׵�����                 !
!  !
!  !
!  !
!  !
! Reference                                                                     !
!                                                                               !
!                                                                               !
! Input                                                                         !
!     cc1,cc2		:						constants of scheme(not used)       !
!     ni		 	:						dimension of calculating direction  !
!     nmax		    :						max dimension of three direction    !
!     nl    	    :						number of flux components           !
!     method,cmethd :						constants of method(not used)       !
!     flux_type 	:						the type used to calculate the flux !
!     limiter 		:						limiter                             !
!     efix		 	:						constant of flux                    !
!     q		 	    :						primitive variables and gama        !
!											 q(1:5,i)=ro,u,v,w,e;  q(6,i)=gama  !
!     trxyz		    :						kcx,kcy,kcz,kct,jacobian            !
! Output                                                                        !
!     df			:						derivative of convect flux          !
!-------------------------------------------------------------------------------!
	real,parameter:: CL1=1.0/16.0,CL2=10.0/16.0,CL3=5.0/16.0
	real,parameter:: CR1=CL3     ,CR2=CL2      ,CR3=CL1

	real,parameter:: EPS=1.0e-6

	real,external :: limiter,minmod
	external flux_type
	integer :: ni,nmax,nl,method,cmethd
	real  :: cc1,cc2,efix,trxyz(5,nmax)
	real    :: q(1:nl+1,-2:nmax+3,2),df(nl,nmax)
	real    :: g(nl,3,-1:ni+1) ,s(nl,3,-1:ni+1) ,bl(3),br(3)
	real    :: wl(nl,3,-1:ni+1),wr(nl,3,-1:ni+1),CL(3),CR(3)
	real    :: ql(nl),qr(nl),gamaeq
	real    :: f(nl,0:nmax),flr(nl),qwl(nl,0:nmax),qwr(nl,0:nmax)
	real    :: EIS2,IS,nx(5,0:nmax),kx,ky,kz,kt
	integer :: i,m,n,i1,ist,ied
	integer :: nerror

	CL(1)=CL1;CL(2)=CL2;CL(3)=CL3
	CR(1)=CR1;CR(2)=CR2;CR(3)=CR3
!
! average geometry derivative by 4th order
!
	call VALUE_HALF_NODE(5,nmax,ni,trxyz,nx)

	if(q(1,-2,1) < small)then
		ist = 2
	else
		ist =0
	endif

	if(q(1,ni+3,1) < small)then
		ied = ni-2
	else
		ied = ni
	endif

!
!	calculating three 1th and 2th derivatives & corresponding weighted constants
!
	do m=1,nl
		i = ist-1
	    s(m,2,i) = q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1)
		s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)

		do i=ist,ied+1

			g(m,1,i) = 0.5*(     q(m,i-2,1) - 4.0*q(m,i-1,1) + 3.0*q(m,i,  1))
			g(m,2,i) = 0.5*(     q(m,i+1,1) -     q(m,i-1,1)                 )
			g(m,3,i) = 0.5*(-3.0*q(m,i,  1) + 4.0*q(m,i+1,1) -     q(m,i+2,1))

			s(m,1,i) = s(m,2,i-1) !q(m,i-2,1) - 2.0*q(m,i-1,1) + q(m,i,  1)
			s(m,2,i) = s(m,3,i-1) !q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1) 
			s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)

			do n=1,3
				IS = g(m,n,i)*g(m,n,i) + s(m,n,i)*s(m,n,i)
				EIS2 = (EPS+IS)**2
				bl(n) = CL(n)/EIS2
				br(n) = CR(n)/EIS2
			enddo

			IS = bl(1) + bl(2) + bl(3)
			do n=1,3
				wl(m,n,i) = bl(n)/IS
			enddo

			IS = br(1) + br(2) + br(3)
			do n=1,3
				wr(m,n,i) = br(n)/IS
			enddo

        enddo
	enddo

!
!	calculating ql & qr
!
	do i=ist,ied
		i1 = i+1
		do m=1,nl
			qwl(m,i) = q(m,i ,1) + 0.125*(wl(m,1,i )*(s(m,1,i )+4.0*g(m,1,i )) + &
								          wl(m,2,i )*(s(m,2,i )+4.0*g(m,2,i )) + &
								          wl(m,3,i )*(s(m,3,i )+4.0*g(m,3,i )) )
			qwr(m,i) = q(m,i1,1) + 0.125*(wr(m,1,i1)*(s(m,1,i1)-4.0*g(m,1,i1)) + &
			                              wr(m,2,i1)*(s(m,2,i1)-4.0*g(m,2,i1)) + &
							              wr(m,3,i1)*(s(m,3,i1)-4.0*g(m,3,i1)) )  
		enddo
	enddo

! *********************************************
! �����õ�����ϵ�ֵ����߽��ϲ������Ը߽ײ�ֵ 
	if( ist > 1 )then
	  do m=1,nl
		qwl(m,0) = (35.d0*q(m,1,1) - 35.d0*q(m,2,1) + 21.d0*q(m,3,1) - 5.d0*q(m,4,1))/16.d0  !cic
		qwl(m,1) = (5.d0*q(m,1,1)  + 15.d0*q(m,2,1) -  5.d0*q(m,3,1) +  q(m,4,1))/16.d0
		qwl(m,2) = (    -q(m,1,1)  +  9.d0*q(m,2,1) +  9.d0*q(m,3,1) -  q(m,4,1))/16.d0
		qwr(m,0) = qwl(m,0)
		qwr(m,1) = qwl(m,1)
			
	  enddo
	endif

	if( ied < ni )then
	  do m=1,nl
		qwr(m,ni)= (35.d0*q(m,ni,1)-35.d0*q(m,ni-1,1)+21.d0*q(m,ni-2,1)-5.d0*q(m,ni-3,1))/16.d0 !cic
		qwr(m,ni-1)= (5.d0*q(m,ni,1) + 15.d0*q(m,ni-1,1) - 5.d0*q(m,ni-2,1) + q(m,ni-3,1))/16.d0
		qwr(m,ni-2)= (    -q(m,ni,1) +  9.d0*q(m,ni-1,1) + 9.d0*q(m,ni-2,1) - q(m,ni-3,1))/16.d0
		qwl(m,ni)= qwr(m,ni)
		qwl(m,ni-1)= qwr(m,ni-1)
	  enddo
	endif

	nerror=0
!____________________
!  do i=1,ni-1
   do i=0,ni   !___cic

		do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror=nerror+1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		endif

		if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror=nerror+1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
			enddo
		endif

		kx = nx(1,i)
		ky = nx(2,i)
		kz = nx(3,i)
		kt = nx(4,i)

		gamaeq = q(nl+1,i,1) !0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))    

		call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		do m=1,nl
			f(m,i) = flr(m)
		enddo

	enddo

	call FLUX_DXYZ(nl,nmax,ni,f,df)


	if(nerror > 8)then
		!write(*,*) 'WCNSE5 failed points along a line:', nerror
	endif

	return
end subroutine WCNS_E_5
!=============================================================================!
!=============================================================================!
subroutine FLUX_DXYZ(nl,nmax,ni,f,df)
	use define_precision_mod
	use global_variables,only: nijk2nd
	implicit none

	real(prec) :: AC,BC,CC,DC
	real(prec) :: A1,B1,C1,D1,E1
	real(prec) :: A2,B2,C2,D2,DD

	integer :: nmax,nl,ni,i,m
	real    :: f(nl,0:nmax),df(nl,nmax)


    AC= 2250.0_prec;BC=-125.0_prec;CC=9.0_prec;  DC=1920.0_prec
    A1=-22.0_prec;  B1= 17.0_prec; C1= 9.0_prec; D1=-5.0_prec; E1=1.0_prec
!	A1=-23.0_prec;  B1= 21.0_prec; C1= 3.0_prec; d1=-1.0_prec; E1=0.0_prec !*TGH NEW WCNS_E_5_BORDER
	A2=  1.0_prec;  B2=-27.0_prec; C2= 27.0_prec;D2=-1.0_prec; DD=24.0_prec;

	if( ni <= nijk2nd )then  !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��

!		do i=2,ni-1
		do i=1,ni   !___cic
			do m=1,nl
				df(m,i) = f(m,i)  - f(m,i-1)  ! degrade to 2th order
			enddo      
		enddo

	else

       do i=3,ni-2
       do m=1,nl
          df(m,i) = ( 2250.0 * ( f(m,i) - f(m,i-1) )  - 125.0 * ( f(m,i+1) - f(m,i-2) ) + 9.0 * ( f(m,i+2) - f(m,i-3) ) ) / 1920.0          ! sixth-order
       end do
       end do
       
       do m=1,nl
          df(m,2)    =  ( f(m,0)  - 27.0*f(m,1)    + 27.0*f(m,2)    - f(m,3)    ) / 24.0      ! fourth-order
          df(m,ni-1) = -( f(m,ni) - 27.0*f(m,ni-1) + 27.0*f(m,ni-2) - f(m,ni-3) ) / 24.0      ! fourth-order
          df(m,1)  =  ( -23.0*f(m,0)  + 21.0*f(m,1)    + 3.0*f(m,2)    - f(m,3)    ) / 24.0   ! third-order
          df(m,ni) = -( -23.0*f(m,ni) + 21.0*f(m,ni-1) + 3.0*f(m,ni-2) - f(m,ni-3) ) / 24.0   ! third-order
       end do
       
	endif

	return
end subroutine FLUX_DXYZ
!=============================================================================!
!=============================================================================!
subroutine WCNSE5_VIS_FLUX  
  use global_variables
  use duvwt_all_field
  implicit none
!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine WCNSE5_VIS_FLUX is used to calculate the viscous flux        !
!		by the use of a type of fifth order accurate weighted compact nonlinear   !
!		difference schemes(WCNS-E-5).                                             !                                                        !
!	Called by                                                                   !
!			Subroutine r_h_s in RHS.f90                                             !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			February,2006                                                           !                                                          !
! Reference                                                                   !
!     Liu Xin.PhD thesis                                                      ��
! Input                                                                       !
!			global_variables                                                        !
! Output                                                                      !
!     derivative of viscous flux                                              ��
!	modified by                                                                 !
!			Zhang Yifeng                                                            !
!	Date                                                                        !
!			August,2008                                                             !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
	real    :: duvwtdxyz (12,0:nmax),nxyz_line(3,0:nmax)
	real    :: duvwt_line(12,nmax),uvwt_line(4,nmax)
	real    :: duvwt_half(12,0:nmax),uvwt_half(4,0:nmax)
	real    :: kxyz_line ( 9,nmax),kxyz_half(9,0:nmax)
	real    :: vis_line(nmax),vol_line(nmax),fv(nl,0:nmax),dfv(nl,nmax)
	real    :: vis_half(0:nmax),vol_half(0:nmax),txyz_half(9,0:nmax)
	real    :: vsl_half(0:nmax),vst_half(0:nmax),cp
	real    :: vsl_line(nmax),vst_line(nmax)

	allocate(duvwt(12,ni,nj,nk))   ! overcome stack over problem

	allocate(duvwt_mid(12,0:ni,0:nj,0:nk))   !___ ֱ�Ӽ������1�׵���  2009.2.1

	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)

! to calculate the deriative of u,v,w,t in computing	coordinate

	call UVWT_DER_4th_ORDER (12) !*TGH. ֱ����ڵ�1�׵�������ʱû�������ֵ
    call UVWT_DER_4th_ORDER_half(12) !*TGH. ֱ�������1�׵���
         !*TGH. �޸�ǰҪ������ϵ�ֵ���޸ĺ�������ϵ�ֵ
!**************** I direction ****************
	do k=1,nk       !___cic
		do j=1,nj     !___cic

			do i=1,ni

				uvwt_line(1,i) = u(i,j,k)
				uvwt_line(2,i) = v(i,j,k)
				uvwt_line(3,i) = w(i,j,k)
				uvwt_line(4,i) = t(i,j,k)

				vis_line (  i) = visl(i,j,k)+vist(i,j,k)
				vsl_line (  i) = visl(i,j,k)
				vst_line (  i) = vist(i,j,k)

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

				do m=1,12
					duvwt_line(m,i) = duvwt(m,i,j,k)
				enddo

			enddo

			call	VALUE_HALF_NODE(4 ,nmax,ni,uvwt_line ,uvwt_half ) !��ԭʼ�����ڰ����ϵ�ֵ

			call	VALUE_HALF_NODE(12,nmax,ni,duvwt_line,duvwt_half) !�ѽڵ��ϵ�һ�׵�����ֵ������
	                                             !*TGH. ����������ظ����㣬i�����ò�ֵ��ʹ��duvwt_mid����

			call	VALUE_HALF_NODE(9 ,nmax,ni,kxyz_line ,kxyz_half ) !����������ֵ������
			call	VALUE_HALF_NODE(1 ,nmax,ni,vol_line  ,vol_half  )
			call	VALUE_HALF_NODE(1 ,nmax,ni,vis_line  ,vis_half  )
			call	VALUE_HALF_NODE(1 ,nmax,ni,vsl_line  ,vsl_half  )
			call	VALUE_HALF_NODE(1 ,nmax,ni,vst_line  ,vst_half  )

			do i=0,ni
				do m=1,4
					duvwt_half(m,i) = duvwt_mid(m,i,j,k)    !___ͬ����
				enddo
			enddo


			call	DUVWT_DXYZ(nmax,ni,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)
!*TGH. ����ռ�����ϵ�һ�׵�������������ռ�ĵ���
			
			call  STRESS_LINE(nmax,ni,12,duvwtdxyz,vis_half,9,txyz_half)
!*TGH. ��ð�������������Ӧ��
			
			do i=0,ni    !___cic
					nxyz_line(1,i) = kxyz_half(1,i)   !�������ȶ��������ռ��е�һ�׵���
					nxyz_line(2,i) = kxyz_half(4,i)
					nxyz_line(3,i) = kxyz_half(7,i)
			enddo


			call  FLUX_VIS_LINE(nl,nmax,ni,4,uvwt_half,9,txyz_half,3,nxyz_line, &
			                    12,duvwtdxyz,vsl_half,vst_half,prl,prt,cp,fv    )
		  !*TGH. ���߼���ճ��ͨ���ڼ���ռ�����ϵ�ͨ��

			call FLUX_DXYZ(nl,nmax,ni,fv,dfv)

			do i = 1,ni     !___cic
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
				enddo

			enddo

		enddo
	enddo
!
!**************** J direction ****************
	do k=1,nk       !___cic
		do i=1,ni     !___cic
			do j=1,nj

				uvwt_line(1,j) = u(i,j,k)
				uvwt_line(2,j) = v(i,j,k)
				uvwt_line(3,j) = w(i,j,k)
				uvwt_line(4,j) = t(i,j,k)

				vis_line (  j) = visl(i,j,k)+vist(i,j,k)
				vsl_line (  j) = visl(i,j,k)
				vst_line (  j) = vist(i,j,k)

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

				do m=1,12
					duvwt_line(m,j) = duvwt(m,i,j,k)
				enddo

			enddo

			call	VALUE_HALF_NODE(4 ,nmax,nj,uvwt_line ,uvwt_half )
			call	VALUE_HALF_NODE(12,nmax,nj,duvwt_line,duvwt_half)
			call	VALUE_HALF_NODE(9 ,nmax,nj,kxyz_line ,kxyz_half )
			call	VALUE_HALF_NODE(1 ,nmax,nj,vol_line  ,vol_half  )
			call	VALUE_HALF_NODE(1 ,nmax,nj,vis_line  ,vis_half  )
			call	VALUE_HALF_NODE(1 ,nmax,nj,vsl_line  ,vsl_half  )
			call	VALUE_HALF_NODE(1 ,nmax,nj,vst_line  ,vst_half  )

			do j=0,nj
				do m=5,8
					duvwt_half(m,j) = duvwt_mid(m,i,j,k)    !___ͬ����
				enddo
			enddo

			call	DUVWT_DXYZ(nmax,nj,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)
			
			call  STRESS_LINE(nmax,nj,12,duvwtdxyz,vis_half,9,txyz_half)

			do j=0,nj    !___cic
					nxyz_line(1,j) = kxyz_half(2,j)
					nxyz_line(2,j) = kxyz_half(5,j)
					nxyz_line(3,j) = kxyz_half(8,j)
			enddo 

			call  FLUX_VIS_LINE(nl,nmax,nj,4,uvwt_half,9,txyz_half,3,nxyz_line, &
			                    12,duvwtdxyz,vsl_half,vst_half,prl,prt,cp,fv    )

			call	FLUX_DXYZ(nl,nmax,nj,fv,dfv)

			do j = 1,nj     !___cic
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j)
				enddo
			enddo

		enddo
	enddo
! 
!**************** K direction ****************
	do j=1,nj        !___cic
		do i=1,ni      !___cic
			do k=1,nk

				uvwt_line(1,k) = u(i,j,k)
				uvwt_line(2,k) = v(i,j,k)
				uvwt_line(3,k) = w(i,j,k)
				uvwt_line(4,k) = t(i,j,k)

				vis_line (  k) = visl(i,j,k)+vist(i,j,k)
				vsl_line (  k) = visl(i,j,k)
				vst_line (  k) = vist(i,j,k)

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

				do m=1,12
					duvwt_line(m,k) = duvwt(m,i,j,k)
				enddo
			enddo

			call	VALUE_HALF_NODE(4 ,nmax,nk,uvwt_line ,uvwt_half )
			call	VALUE_HALF_NODE(12,nmax,nk,duvwt_line,duvwt_half)
			call	VALUE_HALF_NODE(9 ,nmax,nk,kxyz_line ,kxyz_half )
			call	VALUE_HALF_NODE(1 ,nmax,nk,vol_line  ,vol_half  )
			call	VALUE_HALF_NODE(1 ,nmax,nk,vis_line  ,vis_half  )
			call	VALUE_HALF_NODE(1 ,nmax,nk,vsl_line  ,vsl_half  )
			call	VALUE_HALF_NODE(1 ,nmax,nk,vst_line  ,vst_half  )

			do k=0,nk
				do m=9,12
					duvwt_half(m,k) = duvwt_mid(m,i,j,k)    !___ͬ����
				enddo
			enddo

			call	DUVWT_DXYZ(nmax,nk,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)
			
			call  STRESS_LINE(nmax,nk,12,duvwtdxyz,vis_half,9,txyz_half)

			do k=0,nk    !___cic
					nxyz_line(1,k) = kxyz_half(3,k)
					nxyz_line(2,k) = kxyz_half(6,k)
					nxyz_line(3,k) = kxyz_half(9,k)
			enddo 

			call  FLUX_VIS_LINE(nl,nmax,nk,4,uvwt_half,9,txyz_half,3,nxyz_line, &
			                    12,duvwtdxyz,vsl_half,vst_half,prl,prt,cp,fv    )

			call	FLUX_DXYZ(nl,nmax,nk,fv,dfv)

			do k = 1,nk     !___cic
				do m=1,nl
					dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
				enddo
			enddo

		enddo
	enddo

	deallocate(duvwt)
	deallocate(duvwt_mid)
  return
end subroutine WCNSE5_VIS_FLUX
!=============================================================================!
subroutine UVWT_DER_4th_ORDER_half(nv)    !____ֱ�������1�׵��� 2009.2.1
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt_mid
  implicit none
	integer :: i,j,k,m,nv
	real    :: uvwt(4,-1:nmax+1),duvwt(4,0:nmax)

! I direction

	do k=1,nk
		do j=1,nj
			do i=0,ni+1

				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)

			enddo

			call DUVWT_halfNODE_tgh(nmax,ni,4,uvwt,duvwt) !*TGH. �޸ĺ�������ϵ�ֵ

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

			call DUVWT_halfNODE_tgh(nmax,nj,4,uvwt,duvwt) !*TGH. �޸ĺ�������ϵ�ֵ

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

			call DUVWT_halfNODE_tgh(nmax,nk,4,uvwt,duvwt) !*TGH. �޸ĺ�������ϵ�ֵ

			do k=0,nk
				do m=1,4

					duvwt_mid(m+8,i,j,k)=duvwt(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_4th_ORDER_half
!=============================================================================!
subroutine UVWT_DER_4th_ORDER (nv)
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt
  implicit none
	integer :: i,j,k,m,nv
	real    :: uvwt(4,-1:nmax+1),dtem(4,nmax)

! I direction

	do k=1,nk
		do j=1,nj
			do i=0,ni+1

				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)

			enddo

			call DUVWT_NODE(nmax,ni,4,uvwt,dtem) !*TGH. ��ʱû�������ֵ

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

			call DUVWT_NODE(nmax,nj,4,uvwt,dtem)

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

			call DUVWT_NODE(nmax,nk,4,uvwt,dtem)

			do k=1,nk
				do m=1,4

					duvwt(m+8,i,j,k)=dtem(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_4th_ORDER
!=============================================================================!
subroutine DUVWT_NODE(nmax,ni,n1,uvwt,duvwt)  !___DERINODE
	use global_variables,only: nijk2nd
    implicit none

	real,parameter :: A1=8.d0,B1=-1.d0
	real,parameter :: A2=-3.d0,B2=-10.d0,C2=18.d0,D2=-6.d0,E2=1.d0 
	real,parameter :: AC2=-25.,BC2= 48.0,CC2=-36.,DC2=16.0,EC2=-3. ! 2008.8.1 

	integer :: nmax,i,m,ni,n1,n2
	real    :: uvwt(n1,-1:nmax+1),duvwt(n1,nmax)

	if( ni <= nijk2nd ) then  !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
		do i=2,ni-1
		do m=1,n1
				duvwt(m,i) = 0.5*(uvwt(m,i+1) - uvwt(m,i-1))
		enddo
		enddo
		
		do m=1,n1
				duvwt(m,1) = uvwt(m,2 ) - uvwt(m,1)
				duvwt(m,ni)= uvwt(m,ni) - uvwt(m,ni-1)
		enddo
	
	else
		do i=3,ni-2
		do m=1,n1
				duvwt(m,i) = (A1*(uvwt(m,i+1) - uvwt(m,i-1)) + &
										  B1*(uvwt(m,i+2) - uvwt(m,i-2)))/12.d0
		enddo
		enddo

		do m=1,n1
			duvwt(m,2   ) =  (A2*uvwt(m,1) + B2*uvwt(m,2) + C2*uvwt(m,3) + &
			                  D2*uvwt(m,4) + E2*uvwt(m,5))/12.d0
			duvwt(m,ni-1) = -(A2*uvwt(m,ni  ) + B2*uvwt(m,ni-1) + C2*uvwt(m,ni-2) + &
			                  D2*uvwt(m,ni-3) + E2*uvwt(m,ni-4))/12.d0


			duvwt(m,1   ) =  (-11.d0*uvwt(m,1) +18.d0*uvwt(m,2)   -9*uvwt(m,3)   +2*uvwt(m,4))/6.d0  !___cic ����� 
			duvwt(m,ni  ) = -(-11.d0*uvwt(m,ni)+18.d0*uvwt(m,ni-1)-9*uvwt(m,ni-2)+2*uvwt(m,ni-3))/6.d0

		enddo
	endif
!
  return
end subroutine DUVWT_NODE
!=============================================================================!
subroutine DUVWT_halfNODE(nmax,ni,n1,uvwt,duvwt)  !___DERIEGE   2009.2.1
	use global_variables,only: nijk2nd
    implicit none

	real,parameter :: A1=27.d0,B1=-1.d0
	real,parameter :: A2=-22.d0,B2=17.d0,C2=9.d0,D2=-5.d0,E2=1.d0  !*TGH. OLD

	integer :: nmax,i,m,ni,n1,n2
	real    :: uvwt(n1,-1:nmax+1),duvwt(n1,0:nmax)

	if( ni <= nijk2nd )then  !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
!
! deal with the symmetric boundary condiction 
		do i=0,ni
		do m=1,n1
				duvwt(m,i) = uvwt(m,i+1) - uvwt(m,i)
		enddo
		enddo
	else  
		do i=2,ni-2
		do m=1,n1
				duvwt(m,i) = (A1*(uvwt(m,i+1) - uvwt(m,i)) + &
										  B1*(uvwt(m,i+2) - uvwt(m,i-1)))/24.d0
		enddo
		enddo

		do m=1,n1
			duvwt(m,1   ) =  (A2*uvwt(m,1)  + B2*uvwt(m,2)    + C2*uvwt(m,3)    + D2*uvwt(m,4)    + E2*uvwt(m,5)   )/24.d0
			duvwt(m,ni-1) = -(A2*uvwt(m,ni) + B2*uvwt(m,ni-1) + C2*uvwt(m,ni-2) + D2*uvwt(m,ni-3) + E2*uvwt(m,ni-4))/24.d0

			duvwt(m,0   ) =  (A2*uvwt(m,0)    + B2*uvwt(m,1)  + C2*uvwt(m,2)   + D2*uvwt(m,3)   + E2*uvwt(m,4)   )/24.d0 !��1��
			duvwt(m,ni  ) = -(A2*uvwt(m,ni+1) + B2*uvwt(m,ni) + C2*uvwt(m,ni-1)+ D2*uvwt(m,ni-2)+ E2*uvwt(m,ni-3))/24.d0 !��1��

		enddo

	endif
!
  return
end subroutine DUVWT_halfNODE
!=============================================================================!
subroutine VALUE_HALF_NODE(n,nmax,ni,q,q_half) !___INTERV2
	use global_variables,only: nijk2nd
	use define_precision_mod
	implicit none

	real(prec) :: A1,B1,A2,B2,C2,D2

	integer :: nmax,n,ni,m,i
	real    :: q(n,nmax),q_half(n,0:nmax)

	A1=9.0_prec
	B1=-1.0_prec
    A2=5.0_prec
	B2=15.0_prec
	C2=-5.0_prec
	D2=1.0_prec

!
! deal with the symmetric boundary condiction
!
	if( ni <= nijk2nd )then !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��

		do i=1,ni-1
		do m=1,n
				q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
		enddo
		enddo

		do m=1,n
		  q_half(m, 0) = ( 15._prec*q(m, 1)-10._prec*q(m,min(2,ni))  + 3._prec*q(m,min(3,ni))   )/8._prec
		  q_half(m,ni) = ( 15._prec*q(m,ni)-10._prec*q(m,max(ni-1,1))+ 3._prec*q(m,max(ni-2,1)) )/8._prec
		enddo

	else

		do i=2,ni-2  !
		do m=1,n
				q_half(m,i) = (A1*(q(m,i) + q(m,i+1)) + B1*(q(m,i+2) + q(m,i-1)))/16.0
		enddo
		enddo

		do m=1,n
			q_half(m,1   ) = (A2*q(m,1 )+B2*q(m,2   )+C2*q(m,3   )+D2*q(m,4   ))/16.0
			q_half(m,ni-1) = (A2*q(m,ni)+B2*q(m,ni-1)+C2*q(m,ni-2)+D2*q(m,ni-3))/16.0

			q_half(m,0   ) = (35.d0*q(m,1 )-35.d0*q(m,2   )+21.d0*q(m,3   )-5.d0*q(m,4   ))/16.0 !cic
			q_half(m,ni  ) = (35.d0*q(m,ni)-35.d0*q(m,ni-1)+21.d0*q(m,ni-2)-5.d0*q(m,ni-3))/16.0 !cic

		enddo

	endif

	return
end subroutine VALUE_HALF_NODE
!=============================================================================!
subroutine DUVWT_DXYZ(nmax,ni,n1,duvwt,n2,kxyz,vol,duvwtdxyz)
!*TGH. �����ռ�����ϵ�һ�׵�������������ռ�ĵ���
	implicit none

	integer :: nmax,ni,i,m,m1,m2,m3,n1,n2
	real,dimension(n1,0:nmax) :: duvwtdxyz
	real,dimension(n1,0:nmax) :: duvwt
	real,dimension(n2,0:nmax) :: kxyz
	real,dimension(   0:nmax) :: vol

	do i=0,ni
		do m=1,3

			m1=m+m+m-2
			m2=m1+1
			m3=m2+1

			duvwtdxyz(m  ,i) =(kxyz(m1,i)*duvwt(1 ,i) + &
			                   kxyz(m2,i)*duvwt(5 ,i) + &
											   kxyz(m3,i)*duvwt(9 ,i))/vol(i)

			duvwtdxyz(m+3,i) =(kxyz(m1,i)*duvwt(2 ,i) + &
			                   kxyz(m2,i)*duvwt(6 ,i) + &
											   kxyz(m3,i)*duvwt(10,i))/vol(i)

			duvwtdxyz(m+6,i) =(kxyz(m1,i)*duvwt(3 ,i) + &
			                   kxyz(m2,i)*duvwt(7 ,i) + &
											   kxyz(m3,i)*duvwt(11,i))/vol(i)

			duvwtdxyz(m+9,i) =(kxyz(m1,i)*duvwt(4 ,i) + &
			                   kxyz(m2,i)*duvwt(8 ,i) + &
											   kxyz(m3,i)*duvwt(12,i))/vol(i)
		enddo
	enddo

	return
end subroutine DUVWT_DXYZ
!=============================================================================!
subroutine STRESS_LINE(nmax,ni,n1,duvwtdxyz,vis,n2,txyz)
!*TGH. �����������Ӧ��
	implicit none
	real,parameter :: CC=2.0/3.0
	integer :: nmax,i,ni,m,n1,n2
	real    :: vis(0:nmax),txyz(n2,0:nmax),duvwtdxyz(n1,0:nmax),vs,vscc
	real    :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

	do i=0,ni     !___cic

		vs   = vis(i)

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
			

    txyz(1,i) = vscc * ( 2.0*dudx - dvdy - dwdz )
    txyz(5,i) = vscc * ( 2.0*dvdy - dwdz - dudx )
		txyz(9,i) = vscc * ( 2.0*dwdz - dudx - dvdy )
    txyz(2,i) = vs * ( dudy + dvdx )
    txyz(3,i) = vs * ( dudz + dwdx )
    txyz(6,i) = vs * ( dvdz + dwdy )
    txyz(4,i) = txyz(2,i)
    txyz(7,i) = txyz(3,i)
    txyz(8,i) = txyz(6,i)

	enddo

	return
end subroutine STRESS_LINE
!=============================================================================!
subroutine FLUX_VIS_LINE(nl,nmax,ni,n1,uvwt,n2,txyz,n3,kxyz,n4,duvwtdxyz, &
                         vsl,vst,prl,prt,cp,fv)
	implicit none
	
	integer :: nmax,ni,i,m,nl,n1,n2,n3,n4
	real    :: uvwt(n1,0:nmax),txyz(n2,0:nmax),kxyz(n3,0:nmax),fv(nl,0:nmax)
	real    :: duvwtdxyz(n4,0:nmax),vsl(0:nmax),vst(0:nmax),cp,prl,prt,kcp

	do i=0,ni      !___cic

		kcp = cp*(vsl(i)/prl+vst(i)/prt)

		fv(1,i) = 0.0

		fv(2,i) = txyz(1,i)*kxyz(1,i)+txyz(2,i)*kxyz(2,i)+txyz(3,i)*kxyz(3,i)

		fv(3,i) = txyz(4,i)*kxyz(1,i)+txyz(5,i)*kxyz(2,i)+txyz(6,i)*kxyz(3,i)

		fv(4,i) = txyz(7,i)*kxyz(1,i)+txyz(8,i)*kxyz(2,i)+txyz(9,i)*kxyz(3,i)

		fv(5,i) = uvwt(1,i)*fv(2,i)+uvwt(2,i)*fv(3,i)+uvwt(3,i)*fv(4,i)    +  &
							kcp*(duvwtdxyz(10,i)*kxyz(1,i)+duvwtdxyz(11,i)*kxyz(2,i) +  &
							duvwtdxyz(12,i)*kxyz(3,i) )
		          
		
	enddo

	return
end subroutine FLUX_VIS_LINE
!=============================================================================!
subroutine GRID_DERIVATIVE1_4th_ORDER(nb)
	use define_precision_mod
	use global_variables,only:ni,nj,nk,x,y,z,kcx,kcy,kcz, &
                            etx,ety,etz,ctx,cty,ctz,vol,nmax,gcl,nijk2nd
	implicit none
!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine GRID_DERIVATIVE1_4th_ORDER is used to calculate the mesh     !
!   derivatives flux by the use of 4th orders                                 !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			February,2006                                                           !                                                          !
! Reference                                                                   !
!     Liu Xin.PhD thesis & AIAA 99-0557 & AIAA 2003-620                       ��
!
!   Modified by TU Guohua, February,2009
!
!-----------------------------------------------------------------------------!
	real(prec) :: dd,half,aa,bb,cc
	real(prec) :: xkc,xet,xct,ykc,yet,yct,zkc,zet,zct,xm,ym,zm
	real(prec) :: temp(6,nmax),dtemp(6,nmax),temp1(6,0:nmax) !*tgh. new
	integer :: nb,i,j,k,m,MG
	real(prec),pointer,dimension(:,:,:,:) :: dxyzdkecxyz,xkcetct
	integer :: order

	order = 2  !* order = 2 -- ���ü����غ���ʱ���ڲ��������ľ���Ϊ2��
	           !*              ��������������ü����غ���ʱ,����������Ϊ4��

	dd=1._prec/12._prec ; half=0.5_prec ; aa=1.5_prec ; bb=-2.0_prec ; cc=0.5_prec


    do MG=1,GCL !*tgh. ���� Visbal conservative metrix
	    allocate( dxyzdkecxyz(18,ni,nj,nk), xkcetct(3,ni,nj,nk) )
	enddo

    if(order == 2) then
	    do k = 1,nk	
	    do j = 1,nj	
	    do i = 1,ni

		   if(ni==2)then
			  if(i==1)then
		  		 xkc = -x(1,j,k) + x(2,j,k)
				 ykc = -y(1,j,k) + y(2,j,k)
				 zkc = -z(1,j,k) + z(2,j,k)
			  else if(i==ni)then
			    	xkc = x(i,j,k) - x(i-1,j,k)
					ykc = y(i,j,k) - y(i-1,j,k)
					zkc = z(i,j,k) - z(i-1,j,k)
			  endif
		    else

			    if(i==1)then
					xkc = -( aa*x(1,j,k) + bb*x(2,j,k) + cc*x(3,j,k) )
					ykc = -( aa*y(1,j,k) + bb*y(2,j,k) + cc*y(3,j,k) )
					zkc = -( aa*z(1,j,k) + bb*z(2,j,k) + cc*z(3,j,k) )
				else if(i==ni)then
					xkc = aa*x(i,j,k) + bb*x(i-1,j,k) + cc*x(i-2,j,k)
					ykc = aa*y(i,j,k) + bb*y(i-1,j,k) + cc*y(i-2,j,k)
					zkc = aa*z(i,j,k) + bb*z(i-1,j,k) + cc*z(i-2,j,k)
				else
					xkc = half * ( x(i+1,j,k) - x(i-1,j,k) )
					ykc = half * ( y(i+1,j,k) - y(i-1,j,k) )
					zkc = half * ( z(i+1,j,k) - z(i-1,j,k) )
				endif
			endif


			if(nj==2)then
					if(j==1)then
						xet = -x(i,1,k) + x(i,2,k)
						yet = -y(i,1,k) + y(i,2,k)
						zet = -z(i,1,k) + z(i,2,k)
					else if(j==nj)then
						xet = x(i,j,k) - x(i,j-1,k)
						yet = y(i,j,k) - y(i,j-1,k)
						zet = z(i,j,k) - z(i,j-1,k)
					endif
		 	else
									
					if(j==1)then
						xet = -( aa*x(i,1,k) + bb*x(i,2,k) + cc*x(i,3,k) )
						yet = -( aa*y(i,1,k) + bb*y(i,2,k) + cc*y(i,3,k) )
						zet = -( aa*z(i,1,k) + bb*z(i,2,k) + cc*z(i,3,k) )
					else if(j==nj)then
						xet = aa*x(i,j,k) + bb*x(i,j-1,k) + cc*x(i,j-2,k)
						yet = aa*y(i,j,k) + bb*y(i,j-1,k) + cc*y(i,j-2,k)
						zet = aa*z(i,j,k) + bb*z(i,j-1,k) + cc*z(i,j-2,k)
					else
						xet = half * ( x(i,j+1,k) - x(i,j-1,k) )
						yet = half * ( y(i,j+1,k) - y(i,j-1,k) )
						zet = half * ( z(i,j+1,k) - z(i,j-1,k) )
					endif
			  endif

			  if(nk==2)then
					if(k==1)then
						xct = -x(i,j,1) + x(i,j,2)
						yct = -y(i,j,1) + y(i,j,2)
						zct = -z(i,j,1) + z(i,j,2)
					else if(k==nk)then
						xct = x(i,j,k) - x(i,j,k-1)
						yct = y(i,j,k) - y(i,j,k-1)
						zct = z(i,j,k) - z(i,j,k-1)
					endif

			  else
					if(k==1)then
						xct = -( aa*x(i,j,1) + bb*x(i,j,2) + cc*x(i,j,3) )
						yct = -( aa*y(i,j,1) + bb*y(i,j,2) + cc*y(i,j,3) )
						zct = -( aa*z(i,j,1) + bb*z(i,j,2) + cc*z(i,j,3) )
					else if(k==nk)then
						xct = aa*x(i,j,k) + bb*x(i,j,k-1) + cc*x(i,j,k-2)
						yct = aa*y(i,j,k) + bb*y(i,j,k-1) + cc*y(i,j,k-2)
						zct = aa*z(i,j,k) + bb*z(i,j,k-1) + cc*z(i,j,k-2)
					else
						xct = half * ( x(i,j,k+1) - x(i,j,k-1) )
						yct = half * ( y(i,j,k+1) - y(i,j,k-1) )
						zct = half * ( z(i,j,k+1) - z(i,j,k-1) )
					endif

			  endif

            vol(i,j,k) = xkc*yet*zct + xet*yct*zkc + xct*ykc*zet &
                        -xkc*yct*zet - xet*ykc*zct - xct*yet*zkc

            kcx(i,j,k) = yet*zct - zet*yct
            kcy(i,j,k) = zet*xct - xet*zct
            kcz(i,j,k) = xet*yct - yet*xct

            etx(i,j,k) = yct*zkc - zct*ykc
            ety(i,j,k) = zct*xkc - xct*zkc
            etz(i,j,k) = xct*ykc - yct*xkc

            ctx(i,j,k) = ykc*zet - zkc*yet
            cty(i,j,k) = zkc*xet - xkc*zet
            ctz(i,j,k) = xkc*yet - ykc*xet

!_____________ Visbal conservative metrix
            do MG=1,GCL !*tgh. Visbal conservative metrix

			xkcetct(1,i,j,k) = xkc
			xkcetct(2,i,j,k) = xet
			xkcetct(3,i,j,k) = xct

			xm = x(i,j,k)
			ym = y(i,j,k)
			zm = z(i,j,k)

			kcx(i,j,k)= xkc*ym
			kcy(i,j,k)= xet*ym
			kcz(i,j,k)= xct*ym

			etx(i,j,k)= ykc*zm
			ety(i,j,k)= yet*zm
			etz(i,j,k)= yct*zm

			ctx(i,j,k)= zkc*xm
			cty(i,j,k)= zet*xm
			ctz(i,j,k)= zct*xm

		    ENDDO

	    enddo
	    enddo
	    enddo
	 
    else !�߽��㷨

    	do k = 1,nk	
	    do j = 1,nj	
	    do i = 1,ni

		   if(ni==2)then
			  if(i==1)then
		  		 xkc = -x(1,j,k) + x(2,j,k)
				 ykc = -y(1,j,k) + y(2,j,k)
				 zkc = -z(1,j,k) + z(2,j,k)
			  else if(i==ni)then
			    	xkc = x(i,j,k) - x(i-1,j,k)
					ykc = y(i,j,k) - y(i-1,j,k)
					zkc = z(i,j,k) - z(i-1,j,k)
			  endif

		    elseif(ni <= nijk2nd)then !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
				if(i==1)then
					xkc = -( aa*x(1,j,k) + bb*x(2,j,k) + cc*x(3,j,k) )
					ykc = -( aa*y(1,j,k) + bb*y(2,j,k) + cc*y(3,j,k) )
					zkc = -( aa*z(1,j,k) + bb*z(2,j,k) + cc*z(3,j,k) )
				else if(i==ni)then
					xkc = aa*x(i,j,k) + bb*x(i-1,j,k) + cc*x(i-2,j,k)
					ykc = aa*y(i,j,k) + bb*y(i-1,j,k) + cc*y(i-2,j,k)
					zkc = aa*z(i,j,k) + bb*z(i-1,j,k) + cc*z(i-2,j,k)
				else
					xkc = half * ( x(i+1,j,k) - x(i-1,j,k) )
					ykc = half * ( y(i+1,j,k) - y(i-1,j,k) )
					zkc = half * ( z(i+1,j,k) - z(i-1,j,k) )
				endif

		     else

				if(i==1)then
						xkc = dd*(-25.0*x(1,j,k)+48.0*x(2,j,k)-36.0*x(3,j,k) &
						          +16.0*x(4,j,k)- 3.0*x(5,j,k))
						ykc = dd*(-25.0*y(1,j,k)+48.0*y(2,j,k)-36.0*y(3,j,k) &
						          +16.0*y(4,j,k)- 3.0*y(5,j,k))
						zkc = dd*(-25.0*z(1,j,k)+48.0*z(2,j,k)-36.0*z(3,j,k) &
						          +16.0*z(4,j,k)- 3.0*z(5,j,k))
				elseif(i==2)then
						xkc = dd*(- 3.0*x(1,j,k)-10.0*x(2,j,k)+18.0*x(3,j,k) &
						          - 6.0*x(4,j,k)+     x(5,j,k))
						ykc = dd*(- 3.0*y(1,j,k)-10.0*y(2,j,k)+18.0*y(3,j,k) &
						          - 6.0*y(4,j,k)+     y(5,j,k))
						zkc = dd*(- 3.0*z(1,j,k)-10.0*z(2,j,k)+18.0*z(3,j,k) &
						          - 6.0*z(4,j,k)+     z(5,j,k))
				elseif(i==ni-1)then
						xkc = -dd*(- 3.0*x(ni  ,j,k)-10.0*x(ni-1,j,k)+18.0*x(ni-2,j,k) &
						           - 6.0*x(ni-3,j,k)+     x(ni-4,j,k))
						ykc = -dd*(- 3.0*y(ni  ,j,k)-10.0*y(ni-1,j,k)+18.0*y(ni-2,j,k) &
						           - 6.0*y(ni-3,j,k)+     y(ni-4,j,k))
						zkc = -dd*(- 3.0*z(ni  ,j,k)-10.0*z(ni-1,j,k)+18.0*z(ni-2,j,k) &
						           - 6.0*z(ni-3,j,k)+     z(ni-4,j,k))
				elseif(i==ni)then
						xkc = -dd*(-25.0*x(ni  ,j,k)+48.0*x(ni-1,j,k)-36.0*x(ni-2,j,k) &
						           +16.0*x(ni-3,j,k)- 3.0*x(ni-4,j,k))
						ykc = -dd*(-25.0*y(ni  ,j,k)+48.0*y(ni-1,j,k)-36.0*y(ni-2,j,k) &
						           +16.0*y(ni-3,j,k)- 3.0*y(ni-4,j,k))
						zkc = -dd*(-25.0*z(ni  ,j,k)+48.0*z(ni-1,j,k)-36.0*z(ni-2,j,k) &
						           +16.0*z(ni-3,j,k)- 3.0*z(ni-4,j,k))
				else
						xkc = dd*(8.0*(x(i+1,j,k)-x(i-1,j,k))-(x(i+2,j,k)-x(i-2,j,k)))
						ykc = dd*(8.0*(y(i+1,j,k)-y(i-1,j,k))-(y(i+2,j,k)-y(i-2,j,k)))
						zkc = dd*(8.0*(z(i+1,j,k)-z(i-1,j,k))-(z(i+2,j,k)-z(i-2,j,k)))
				endif
			endif

			if(nj==2)then 
					if(j==1)then
						xet = -x(i,1,k) + x(i,2,k)
						yet = -y(i,1,k) + y(i,2,k)
						zet = -z(i,1,k) + z(i,2,k)
					else if(j==nj)then
						xet = x(i,j,k) - x(i,j-1,k)
						yet = y(i,j,k) - y(i,j-1,k)
						zet = z(i,j,k) - z(i,j-1,k)
					endif
									
			elseif(nj <= nijk2nd)then  !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
					if(j==1)then
						xet = -( aa*x(i,1,k) + bb*x(i,2,k) + cc*x(i,3,k) )
						yet = -( aa*y(i,1,k) + bb*y(i,2,k) + cc*y(i,3,k) )
						zet = -( aa*z(i,1,k) + bb*z(i,2,k) + cc*z(i,3,k) )
					else if(j==nj)then
						xet = aa*x(i,j,k) + bb*x(i,j-1,k) + cc*x(i,j-2,k)
						yet = aa*y(i,j,k) + bb*y(i,j-1,k) + cc*y(i,j-2,k)
						zet = aa*z(i,j,k) + bb*z(i,j-1,k) + cc*z(i,j-2,k)
					else
						xet = half * ( x(i,j+1,k) - x(i,j-1,k) )
						yet = half * ( y(i,j+1,k) - y(i,j-1,k) )
						zet = half * ( z(i,j+1,k) - z(i,j-1,k) )
					endif

			else
					if(j==1)then
						xet = dd*(-25.0*x(i,1,k)+48.0*x(i,2,k)-36.0*x(i,3,k) &
						          +16.0*x(i,4,k)- 3.0*x(i,5,k))
						yet = dd*(-25.0*y(i,1,k)+48.0*y(i,2,k)-36.0*y(i,3,k) &
						          +16.0*y(i,4,k)- 3.0*y(i,5,k))
						zet = dd*(-25.0*z(i,1,k)+48.0*z(i,2,k)-36.0*z(i,3,k) &
						          +16.0*z(i,4,k)- 3.0*z(i,5,k))
					elseif(j==2)then
						xet = dd*(- 3.0*x(i,1,k)-10.0*x(i,2,k)+18.0*x(i,3,k) &
						          - 6.0*x(i,4,k)+     x(i,5,k))
						yet = dd*(- 3.0*y(i,1,k)-10.0*y(i,2,k)+18.0*y(i,3,k) &
						          - 6.0*y(i,4,k)+     y(i,5,k))
						zet = dd*(- 3.0*z(i,1,k)-10.0*z(i,2,k)+18.0*z(i,3,k) &
						          - 6.0*z(i,4,k)+     z(i,5,k))
					elseif(j==nj-1)then
						xet = -dd*(- 3.0*x(i,nj  ,k)-10.0*x(i,nj-1,k)+18.0*x(i,nj-2,k) &
						           - 6.0*x(i,nj-3,k)+     x(i,nj-4,k))
						yet = -dd*(- 3.0*y(i,nj  ,k)-10.0*y(i,nj-1,k)+18.0*y(i,nj-2,k) &
						           - 6.0*y(i,nj-3,k)+     y(i,nj-4,k))
						zet = -dd*(- 3.0*z(i,nj  ,k)-10.0*z(i,nj-1,k)+18.0*z(i,nj-2,k) &
						           - 6.0*z(i,nj-3,k)+     z(i,nj-4,k))
					elseif(j==nj)then
						xet = -dd*(-25.0*x(i,nj  ,k)+48.0*x(i,nj-1,k)-36.0*x(i,nj-2,k) &
						           +16.0*x(i,nj-3,k)- 3.0*x(i,nj-4,k))
						yet = -dd*(-25.0*y(i,nj  ,k)+48.0*y(i,nj-1,k)-36.0*y(i,nj-2,k) &
						           +16.0*y(i,nj-3,k)- 3.0*y(i,nj-4,k))
						zet = -dd*(-25.0*z(i,nj  ,k)+48.0*z(i,nj-1,k)-36.0*z(i,nj-2,k) &
						           +16.0*z(i,nj-3,k)- 3.0*z(i,nj-4,k))
					else
						xet = dd*(8.0*(x(i,j+1,k)-x(i,j-1,k))-(x(i,j+2,k)-x(i,j-2,k)))
						yet = dd*(8.0*(y(i,j+1,k)-y(i,j-1,k))-(y(i,j+2,k)-y(i,j-2,k)))
						zet = dd*(8.0*(z(i,j+1,k)-z(i,j-1,k))-(z(i,j+2,k)-z(i,j-2,k)))
					endif

			endif

			if(nk==2)then
					if(k==1)then
						xct = -x(i,j,1) + x(i,j,2)
						yct = -y(i,j,1) + y(i,j,2)
						zct = -z(i,j,1) + z(i,j,2)
					else if(k==nk)then
						xct = x(i,j,k) - x(i,j,k-1)
						yct = y(i,j,k) - y(i,j,k-1)
						zct = z(i,j,k) - z(i,j,k-1)
					endif

			elseif(nk <= nijk2nd)then  !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
					if(k==1)then
						xct = -( aa*x(i,j,1) + bb*x(i,j,2) + cc*x(i,j,3) )
						yct = -( aa*y(i,j,1) + bb*y(i,j,2) + cc*y(i,j,3) )
						zct = -( aa*z(i,j,1) + bb*z(i,j,2) + cc*z(i,j,3) )
					else if(k==nk)then
						xct = aa*x(i,j,k) + bb*x(i,j,k-1) + cc*x(i,j,k-2)
						yct = aa*y(i,j,k) + bb*y(i,j,k-1) + cc*y(i,j,k-2)
						zct = aa*z(i,j,k) + bb*z(i,j,k-1) + cc*z(i,j,k-2)
					else
						xct = half * ( x(i,j,k+1) - x(i,j,k-1) )
						yct = half * ( y(i,j,k+1) - y(i,j,k-1) )
						zct = half * ( z(i,j,k+1) - z(i,j,k-1) )
					endif

			else

					if(k==1)then
						xct = dd*(-25.0*x(i,j,1)+48.0*x(i,j,2)-36.0*x(i,j,3) &
						          +16.0*x(i,j,4)- 3.0*x(i,j,5))
						yct = dd*(-25.0*y(i,j,1)+48.0*y(i,j,2)-36.0*y(i,j,3) &
						          +16.0*y(i,j,4)- 3.0*y(i,j,5))
						zct = dd*(-25.0*z(i,j,1)+48.0*z(i,j,2)-36.0*z(i,j,3) &
						          +16.0*z(i,j,4)- 3.0*z(i,j,5))
					elseif(k==2)then
						xct = dd*(- 3.0*x(i,j,1)-10.0*x(i,j,2)+18.0*x(i,j,3) &
						          - 6.0*x(i,j,4)+     x(i,j,5))
						yct = dd*(- 3.0*y(i,j,1)-10.0*y(i,j,2)+18.0*y(i,j,3) &
						          - 6.0*y(i,j,4)+     y(i,j,5))
						zct = dd*(- 3.0*z(i,j,1)-10.0*z(i,j,2)+18.0*z(i,j,3) &
						          - 6.0*z(i,j,4)+     z(i,j,5))
					elseif(k==nk-1)then
						xct = -dd*(- 3.0*x(i,j,nk  )-10.0*x(i,j,nk-1)+18.0*x(i,j,nk-2) &
						           - 6.0*x(i,j,nk-3)+     x(i,j,nk-4))
						yct = -dd*(- 3.0*y(i,j,nk  )-10.0*y(i,j,nk-1)+18.0*y(i,j,nk-2) &
						           - 6.0*y(i,j,nk-3)+     y(i,j,nk-4))
						zct = -dd*(- 3.0*z(i,j,nk  )-10.0*z(i,j,nk-1)+18.0*z(i,j,nk-2) &
						           - 6.0*z(i,j,nk-3)+     z(i,j,nk-4))
					elseif(k==nk)then
						xct = -dd*(-25.0*x(i,j,nk  )+48.0*x(i,j,nk-1)-36.0*x(i,j,nk-2) &
						           +16.0*x(i,j,nk-3)- 3.0*x(i,j,nk-4))
						yct = -dd*(-25.0*y(i,j,nk  )+48.0*y(i,j,nk-1)-36.0*y(i,j,nk-2) &
						           +16.0*y(i,j,nk-3)- 3.0*y(i,j,nk-4))
						zct = -dd*(-25.0*z(i,j,nk  )+48.0*z(i,j,nk-1)-36.0*z(i,j,nk-2) &
						           +16.0*z(i,j,nk-3)- 3.0*z(i,j,nk-4))
					else
						xct = dd*(8.0*(x(i,j,k+1)-x(i,j,k-1))-(x(i,j,k+2)-x(i,j,k-2)))
						yct = dd*(8.0*(y(i,j,k+1)-y(i,j,k-1))-(y(i,j,k+2)-y(i,j,k-2)))
						zct = dd*(8.0*(z(i,j,k+1)-z(i,j,k-1))-(z(i,j,k+2)-z(i,j,k-2)))
					endif

			endif

            vol(i,j,k) = xkc*yet*zct + xet*yct*zkc + xct*ykc*zet &
                        -xkc*yct*zet - xet*ykc*zct - xct*yet*zkc

            kcx(i,j,k) = yet*zct - zet*yct
            kcy(i,j,k) = zet*xct - xet*zct
            kcz(i,j,k) = xet*yct - yet*xct

            etx(i,j,k) = yct*zkc - zct*ykc
            ety(i,j,k) = zct*xkc - xct*zkc
            etz(i,j,k) = xct*ykc - yct*xkc

            ctx(i,j,k) = ykc*zet - zkc*yet
            cty(i,j,k) = zkc*xet - xkc*zet
            ctz(i,j,k) = xkc*yet - ykc*xet

!_____________ Visbal conservative metrix
            do MG=1,GCL !*tgh. Visbal conservative metrix

			xkcetct(1,i,j,k) = xkc
			xkcetct(2,i,j,k) = xet
			xkcetct(3,i,j,k) = xct

			xm = x(i,j,k)
			ym = y(i,j,k)
			zm = z(i,j,k)

			kcx(i,j,k)= xkc*ym
			kcy(i,j,k)= xet*ym
			kcz(i,j,k)= xct*ym

			etx(i,j,k)= ykc*zm
			ety(i,j,k)= yet*zm
			etz(i,j,k)= yct*zm

			ctx(i,j,k)= zkc*xm
			cty(i,j,k)= zet*xm
			ctz(i,j,k)= zct*xm

		    ENDDO

	    enddo
	    enddo
	    enddo
    endif

    do MG=1,GCL !*tgh. ������ Visbal conservative metrix

	    do k=1,nk
	    do j=1,nj

			do i=1,ni
				temp(1,i) = cty(i,j,k) 
				temp(2,i) = ctz(i,j,k) 
				temp(3,i) = kcy(i,j,k) 
				temp(4,i) = kcz(i,j,k) 
				temp(5,i) = ety(i,j,k) 
				temp(6,i) = etz(i,j,k) 
			enddo
			call VALUE_HALF_NODE(6,nmax,ni,temp,temp1)
	        call FLUX_DXYZ(6,nmax,ni,temp1,dtemp)

			do i=1,ni
				do m=1,6
					dxyzdkecxyz(m,i,j,k) = dtemp(m,i)
				enddo
			enddo

	    enddo
	    enddo
  
	    do k=1,nk
	    do i=1,ni

			do j=1,nj
				temp(1,j) = ctx(i,j,k) 
				temp(2,j) = ctz(i,j,k) 
				temp(3,j) = kcx(i,j,k) 
				temp(4,j) = kcz(i,j,k) 
				temp(5,j) = etx(i,j,k) 
				temp(6,j) = etz(i,j,k) 
			enddo
			call VALUE_HALF_NODE(6,nmax,nj,temp,temp1)
	        call FLUX_DXYZ(6,nmax,nj,temp1,dtemp)

			do j=1,nj
				do m=1,6
					dxyzdkecxyz(m+6,i,j,k) = dtemp(m,j)
				enddo
			enddo

	    enddo
	    enddo

	    do j=1,nj
	    do i=1,ni

			do k=1,nk
				temp(1,k) = ctx(i,j,k) 
				temp(2,k) = cty(i,j,k) 
				temp(3,k) = kcx(i,j,k) 
				temp(4,k) = kcy(i,j,k) 
				temp(5,k) = etx(i,j,k) 
				temp(6,k) = ety(i,j,k) 
			enddo
			call VALUE_HALF_NODE(6,nmax,nk,temp,temp1)
	        call FLUX_DXYZ(6,nmax,nk,temp1,dtemp)

			do k=1,nk
				do m=1,6
					dxyzdkecxyz(m+12,i,j,k) = dtemp(m,k)
				enddo
			enddo

	    enddo
	    enddo


	    do i=1,ni
		do j=1,nj
		do k=1,nk

				kcx(i,j,k) = dxyzdkecxyz(18,i,j,k)-dxyzdkecxyz(12,i,j,k)
				kcy(i,j,k) = dxyzdkecxyz(14,i,j,k)-dxyzdkecxyz( 8,i,j,k)
				kcz(i,j,k) = dxyzdkecxyz(16,i,j,k)-dxyzdkecxyz(10,i,j,k)
												
				etx(i,j,k) = dxyzdkecxyz( 6,i,j,k)-dxyzdkecxyz(17,i,j,k)
				ety(i,j,k) = dxyzdkecxyz( 2,i,j,k)-dxyzdkecxyz(13,i,j,k)
				etz(i,j,k) = dxyzdkecxyz( 4,i,j,k)-dxyzdkecxyz(15,i,j,k)
												
				ctx(i,j,k) = dxyzdkecxyz(11,i,j,k)-dxyzdkecxyz( 5,i,j,k)
				cty(i,j,k) = dxyzdkecxyz( 7,i,j,k)-dxyzdkecxyz( 1,i,j,k)
				ctz(i,j,k) = dxyzdkecxyz( 9,i,j,k)-dxyzdkecxyz( 3,i,j,k)

			vol(i,j,k) = xkcetct(1,i,j,k)*kcx(i,j,k) & !( dxyzkecxyz(18,i,j,k) - dxyzkecxyz(12,i,j,k) ) &
			           + xkcetct(2,i,j,k)*etx(i,j,k) & !( dxyzkecxyz(6 ,i,j,k) - dxyzkecxyz(17,i,j,k) ) &
					   + xkcetct(3,i,j,k)*ctx(i,j,k)   !( dxyzkecxyz(11,i,j,k) - dxyzkecxyz(5 ,i,j,k) )                 
												
		enddo
		enddo
	    enddo
!_____________ Visbal conservative metrix  -- end

	    deallocate(dxyzdkecxyz,xkcetct )!*tgh. end Visbal conservative metrix

     ENDDO !* �������ü����غ���
end subroutine GRID_DERIVATIVE1_4th_ORDER
!=============================================================================!


!=============================================================================!
subroutine GRID_DERIVATIVE_gcl(nb)
	use define_precision_mod
	use global_variables,only:ni,nj,nk,x,y,z,kcx,kcy,kcz, &
                            etx,ety,etz,ctx,cty,ctz,vol,nmax,gcl,nscheme
	implicit none
!
!-----------------------------------------------------------------------------!
	real(prec) :: dd,half,aa,bb,cc
	real(prec) :: xkc,xet,xct,ykc,yet,yct,zkc,zet,zct,xm,ym,zm
	real(prec) :: temp(6,nmax),dtemp(6,nmax),temp1(6,0:nmax) !*tgh. new
	real(prec) :: tem3(3,nmax),temp3(3,0:nmax) !*tgh. new
	integer :: nb,i,j,k,m,MG
	real(prec),pointer,dimension(:,:,:,:) :: dxyzdkecxyz,xkcetct
	integer :: nschemegrid

	dd=1._prec/12._prec ; half=0.5_prec ; aa=1.5_prec ; bb=-2.0_prec ; cc=0.5_prec

	nschemegrid = nscheme  !* ������ɢ����

	if(nschemegrid /= nscheme)then
		write(*,*)'������ɢ���Ⱥ�������ɢ���ȷֱ�Ϊ��',nschemegrid, nscheme
		write(*,*)'�����غ����Զ��ر�!!!'
		gcl = 0
	endif

	allocate( dxyzdkecxyz(18,ni,nj,nk), xkcetct(3,ni,nj,nk) )

	do k = 1,nk	
	do j = 1,nj	
	    do i = 1,ni
				tem3(1,i) = x(i,j,k)
				tem3(2,i) = y(i,j,k)
				tem3(3,i) = z(i,j,k)
		enddo

		call value_half_node(3,nmax,ni,tem3,temp3)
		call flux_dxyz(3,nmax,ni,temp3,tem3)

		do i=1,ni
				dxyzdkecxyz(1 ,i,j,k) = tem3(1,i) !
				dxyzdkecxyz(4 ,i,j,k) = tem3(2,i) !
				dxyzdkecxyz(7 ,i,j,k) = tem3(3,i) !
		enddo
		

	enddo
	enddo

	do k=1,nk
	do i=1,ni
	    do j=1,nj
				tem3(1,j) = x(i,j,k)
				tem3(2,j) = y(i,j,k)
				tem3(3,j) = z(i,j,k)
		enddo
		
		call value_half_node(3,nmax,nj,tem3,temp3)
		call flux_dxyz(3,nmax,nj,temp3,tem3)

		do j=1,nj
				dxyzdkecxyz(2 ,i,j,k) = tem3(1,j) !
				dxyzdkecxyz(5 ,i,j,k) = tem3(2,j) !
				dxyzdkecxyz(8 ,i,j,k) = tem3(3,j) !
		enddo

	enddo
	enddo

	do j=1,nj
	do i=1,ni
	    do k=1,nk
				tem3(1,k) = x(i,j,k)
				tem3(2,k) = y(i,j,k)
				tem3(3,k) = z(i,j,k)
		enddo

		call value_half_node(3,nmax,nk,tem3,temp3)
		call flux_dxyz(3,nmax,nk,temp3,tem3)

		do k=1,nk
				dxyzdkecxyz(3 ,i,j,k) = tem3(1,k) !
				dxyzdkecxyz(6 ,i,j,k) = tem3(2,k) !
				dxyzdkecxyz(9 ,i,j,k) = tem3(3,k) !
		enddo
			
	enddo
	enddo


    do k = 1,nk	
    do j = 1,nj	
    do i = 1,ni

        kcx(i,j,k) = dxyzdkecxyz(5 ,i,j,k)*dxyzdkecxyz(9 ,i,j,k) - dxyzdkecxyz(8 ,i,j,k)*dxyzdkecxyz(6 ,i,j,k) !yet*zct - zet*yct
        kcy(i,j,k) = dxyzdkecxyz(8 ,i,j,k)*dxyzdkecxyz(3 ,i,j,k) - dxyzdkecxyz(2 ,i,j,k)*dxyzdkecxyz(9 ,i,j,k) !zet*xct - xet*zct
        kcz(i,j,k) = dxyzdkecxyz(2 ,i,j,k)*dxyzdkecxyz(6 ,i,j,k) - dxyzdkecxyz(5 ,i,j,k)*dxyzdkecxyz(3 ,i,j,k) !xet*yct - yet*xct

        etx(i,j,k) = dxyzdkecxyz(6 ,i,j,k)*dxyzdkecxyz(7 ,i,j,k) - dxyzdkecxyz(9 ,i,j,k)*dxyzdkecxyz(4 ,i,j,k) !yct*zkc - zct*ykc
        ety(i,j,k) = dxyzdkecxyz(9 ,i,j,k)*dxyzdkecxyz(1 ,i,j,k) - dxyzdkecxyz(3 ,i,j,k)*dxyzdkecxyz(7 ,i,j,k) !zct*xkc - xct*zkc
        etz(i,j,k) = dxyzdkecxyz(3 ,i,j,k)*dxyzdkecxyz(4 ,i,j,k) - dxyzdkecxyz(6 ,i,j,k)*dxyzdkecxyz(1 ,i,j,k) !xct*ykc - yct*xkc

        ctx(i,j,k) = dxyzdkecxyz(4 ,i,j,k)*dxyzdkecxyz(8 ,i,j,k) - dxyzdkecxyz(7 ,i,j,k)*dxyzdkecxyz(5 ,i,j,k) !ykc*zet - zkc*yet
        cty(i,j,k) = dxyzdkecxyz(7 ,i,j,k)*dxyzdkecxyz(2 ,i,j,k) - dxyzdkecxyz(1 ,i,j,k)*dxyzdkecxyz(8 ,i,j,k) !zkc*xet - xkc*zet
        ctz(i,j,k) = dxyzdkecxyz(1 ,i,j,k)*dxyzdkecxyz(5 ,i,j,k) - dxyzdkecxyz(4 ,i,j,k)*dxyzdkecxyz(2 ,i,j,k) !xkc*yet - ykc*xet

		vol(i,j,k) = dxyzdkecxyz(1,i,j,k)*kcx(i,j,k) &
		           + dxyzdkecxyz(2,i,j,k)*etx(i,j,k) &
		           + dxyzdkecxyz(3,i,j,k)*ctx(i,j,k)
	enddo
	enddo
	enddo

	do k=1,nk
	do j=1,nj
	do i=1,ni
!_____________ Visbal conservative metrix
        do MG=1,GCL !*tgh. Visbal conservative metrix

			xkcetct(1,i,j,k) = dxyzdkecxyz(1,i,j,k)
			xkcetct(2,i,j,k) = dxyzdkecxyz(2,i,j,k)
			xkcetct(3,i,j,k) = dxyzdkecxyz(3,i,j,k)

			xm = x(i,j,k)
			ym = y(i,j,k)
			zm = z(i,j,k)

			kcx(i,j,k)= dxyzdkecxyz(1 ,i,j,k)*ym
			kcy(i,j,k)= dxyzdkecxyz(2 ,i,j,k)*ym
			kcz(i,j,k)= dxyzdkecxyz(3 ,i,j,k)*ym

			etx(i,j,k)= dxyzdkecxyz(4 ,i,j,k)*zm
			ety(i,j,k)= dxyzdkecxyz(5 ,i,j,k)*zm
			etz(i,j,k)= dxyzdkecxyz(6 ,i,j,k)*zm

			ctx(i,j,k)= dxyzdkecxyz(7 ,i,j,k)*xm
			cty(i,j,k)= dxyzdkecxyz(8 ,i,j,k)*xm
			ctz(i,j,k)= dxyzdkecxyz(9 ,i,j,k)*xm

	    ENDDO

	enddo
	enddo
	enddo
	 


    do MG=1,GCL !*tgh. ���� Visbal conservative metrix

	    do k=1,nk
	    do j=1,nj

			do i=1,ni
				temp(1,i) = cty(i,j,k) 
				temp(2,i) = ctz(i,j,k) 
				temp(3,i) = kcy(i,j,k) 
				temp(4,i) = kcz(i,j,k) 
				temp(5,i) = ety(i,j,k) 
				temp(6,i) = etz(i,j,k) 
			enddo

		    call VALUE_HALF_NODE(6,nmax,ni,temp,temp1)
	        call FLUX_DXYZ(6,nmax,ni,temp1,dtemp)

			do i=1,ni
				do m=1,6
					dxyzdkecxyz(m,i,j,k) = dtemp(m,i)
				enddo
			enddo

	    enddo
	    enddo
  
	    do k=1,nk
	    do i=1,ni

			do j=1,nj
				temp(1,j) = ctx(i,j,k) 
				temp(2,j) = ctz(i,j,k) 
				temp(3,j) = kcx(i,j,k) 
				temp(4,j) = kcz(i,j,k) 
				temp(5,j) = etx(i,j,k) 
				temp(6,j) = etz(i,j,k) 
			enddo

		    call VALUE_HALF_NODE(6,nmax,nj,temp,temp1)
	        call FLUX_DXYZ(6,nmax,nj,temp1,dtemp)

			do j=1,nj
				do m=1,6
					dxyzdkecxyz(m+6,i,j,k) = dtemp(m,j)
				enddo
			enddo

	    enddo
	    enddo

	    do j=1,nj
	    do i=1,ni

			do k=1,nk
				temp(1,k) = ctx(i,j,k) 
				temp(2,k) = cty(i,j,k) 
				temp(3,k) = kcx(i,j,k) 
				temp(4,k) = kcy(i,j,k) 
				temp(5,k) = etx(i,j,k) 
				temp(6,k) = ety(i,j,k) 
			enddo

		    call VALUE_HALF_NODE(6,nmax,nk,temp,temp1)
	        call FLUX_DXYZ(6,nmax,nk,temp1,dtemp)

			do k=1,nk
				do m=1,6
					dxyzdkecxyz(m+12,i,j,k) = dtemp(m,k)
				enddo
			enddo

	    enddo
	    enddo


	    do i=1,ni
		do j=1,nj
		do k=1,nk

				kcx(i,j,k) = dxyzdkecxyz(18,i,j,k)-dxyzdkecxyz(12,i,j,k)
				kcy(i,j,k) = dxyzdkecxyz(14,i,j,k)-dxyzdkecxyz( 8,i,j,k)
				kcz(i,j,k) = dxyzdkecxyz(16,i,j,k)-dxyzdkecxyz(10,i,j,k)
												
				etx(i,j,k) = dxyzdkecxyz( 6,i,j,k)-dxyzdkecxyz(17,i,j,k)
				ety(i,j,k) = dxyzdkecxyz( 2,i,j,k)-dxyzdkecxyz(13,i,j,k)
				etz(i,j,k) = dxyzdkecxyz( 4,i,j,k)-dxyzdkecxyz(15,i,j,k)
												
				ctx(i,j,k) = dxyzdkecxyz(11,i,j,k)-dxyzdkecxyz( 5,i,j,k)
				cty(i,j,k) = dxyzdkecxyz( 7,i,j,k)-dxyzdkecxyz( 1,i,j,k)
				ctz(i,j,k) = dxyzdkecxyz( 9,i,j,k)-dxyzdkecxyz( 3,i,j,k)

			vol(i,j,k) = xkcetct(1,i,j,k)*kcx(i,j,k) & !( dxyzkecxyz(18,i,j,k) - dxyzkecxyz(12,i,j,k) ) &
			           + xkcetct(2,i,j,k)*etx(i,j,k) & !( dxyzkecxyz(6 ,i,j,k) - dxyzkecxyz(17,i,j,k) ) &
					   + xkcetct(3,i,j,k)*ctx(i,j,k)   !( dxyzkecxyz(11,i,j,k) - dxyzkecxyz(5 ,i,j,k) )                 
		enddo
		enddo
	    enddo
!_____________ Visbal conservative metrix  -- end


     ENDDO !* �������ü����غ���

	 deallocate(dxyzdkecxyz,xkcetct )

end subroutine GRID_DERIVATIVE_gcl
!=============================================================================!
!=============================================================================!
!=============================================================================!
subroutine VALUE_HALF_2nd(n,nmax,ni,q,q_half) !___INTERV2
	use define_precision_mod
	implicit none

	real,parameter :: A1=9.0,B1=-1.0
	real,parameter :: A2=5.0,B2=15.0,C2=-5.0,D2=1.0

	integer :: nmax,n,ni,m,i
	real    :: q(n,nmax),q_half(n,0:nmax)
!
! deal with the symmetric boundary condiction
!
		
	
	do i=1,ni-1
	do m=1,n
		q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
	enddo
	enddo

	do m=1,n
!		q_half(m, 0) = 1.5*q(m, 1)-0.5*q(m,min(2,ni))  
!		q_half(m,ni) = 1.5*q(m,ni)-0.5*q(m,max(ni-1,1)) 
		q_half(m, 0) = ( 15._prec*q(m, 1)-10._prec*q(m,min(2,ni))  + 3._prec*q(m,min(3,ni))   )/8._prec
		q_half(m,ni) = ( 15._prec*q(m,ni)-10._prec*q(m,max(ni-1,1))+ 3._prec*q(m,max(ni-2,1)) )/8._prec
	enddo

	return
end subroutine VALUE_HALF_2nd
!=============================================================================!
!=============================================================================!

!=============================================================================!
subroutine FLUX_DXYZ_2nd(nl,nmax,ni,f,df)
	use define_precision_mod
	implicit none

	integer :: nmax,nl,ni,i,m
	real    :: f(nl,0:nmax),df(nl,nmax)

	do i=1,ni    !
	do m=1,nl
		df(m,i) = f(m,i  ) - f(m,i-1)
	enddo
	enddo

	return
end subroutine FLUX_DXYZ_2nd

!=============================================================================!
!=============================================================================!
subroutine WCNS_E_5_J(cc1,cc2,ni,nmax,nl,method,cmethd,flux_type, &
                    limiter,efix,trxyz,q,df)
    use define_precision_mod
    use global_const,only:rmin_limit,pmin_limit
	implicit none
!-----------------------------------------------------------------------------!
! ����Jacobian��Ȩ��WCNS-E-5��ʽ
!	Called by                                                                   !
!			Subroutine inviscd3d in Convect.f90                                     !
!	Editor                                                                      !    
!			TU Guohua                                                          !
!	Date                                                                        !
!			February,2006                                                           !                                                          !
! Reference                                                                   !
!     Liu Xin.PhD thesis                                                      ��
! Input                                                                       !
!     cc1,cc2						:						constants of scheme(not used)             ��
!     ni								:						dimension of calculating direction        ��
!     nmax							:						max dimension of three direction          ��
!     nl    						:						number of flux components                 ��
!     method,cmethd 		:						constants of method(not used)             ��
!     flux_type 		    :						the type used to calculate the flux       ��
!     limiter 					:						limiter                                   ��
!     efix		 					:						constant of flux                          ��
!     q		 					    :						primitive variables and gama              ��
!																		q(1:5,i)=ro,u,v,w,e;  q(6,i)=gama         ��
!     trxyz					    :						kcx,kcy,kcz,kct,jacobian                  ��
! Output                                                                      !
!     df								:						derivative of convect flux                ��
!	modified by                                                                 !
!			TU Guohua                                                            !
!	Date                                                                        !
!			August,2008                                                             !
!-----------------------------------------------------------------------------!
	real,parameter:: CL1=1.0/16.0,CL2=10.0/16.0,CL3=5.0/16.0
	real,parameter:: CR1=CL3     ,CR2=CL2      ,CR3=CL1

	real,parameter:: EPS=1.0e-6

	real,external :: limiter,minmod
	external flux_type
	integer :: ni,nmax,nl,method,cmethd
	real  :: cc1,cc2,efix,trxyz(5,nmax)
	real    :: q(1:nl+1,-1:nmax+1,2),df(nl,nmax)

	real    :: g(nl,3,-1:ni+2) ,s(nl,3,-1:ni+2) ,bl(3),br(3)
	real    :: wl(nl,3,-1:ni+2),wr(nl,3,-1:ni+2),CL(3),CR(3)
	real    :: dq(nl,-3:nmax+3),ql(nl),qr(nl),gamaeq
	real    :: f(nl,0:nmax),flr(nl),qwl(nl,0:nmax+1),qwr(nl,0:nmax+1)
	real    :: EIS2,IS,nx(5,0:nmax),kx,ky,kz,kt
	integer :: i,m,n,i1
	real :: qj(nl,0:ni+1),gj(3,-1:ni+2),sj(3,-1:ni+2),jac_l,jac_r

	CL(1)=CL1;CL(2)=CL2;CL(3)=CL3
	CR(1)=CR1;CR(2)=CR2;CR(3)=CR3
!
! average geometry derivative by 4th order
!
	call VALUE_HALF_NODE(5,nmax,ni,trxyz,nx)
!
!	calculating three 1th and 2th derivatives & corresponding weighted constants
!
	do i=1,ni
    do m=1,nl
		qj(m,i)=trxyz(5,i)*q(m,i,1)
	enddo
	enddo
!    do m=1,nl
!	    qj(m,0   ) = qj(m,1 )
!	    qj(m,ni+1) = qj(m,ni)
!	enddo


    do i=3,ni-2

		m=5
		gj(1,i) = 0.5*(     trxyz(m,i-2) - 4.0*trxyz(m,i-1) + 3.0*trxyz(m,i  )) !Jacobian
		gj(2,i) = 0.5*(     trxyz(m,i+1) -     trxyz(m,i-1)                   )
		gj(3,i) = 0.5*(-3.0*trxyz(m,i  ) + 4.0*trxyz(m,i+1) -     trxyz(m,i+2))

		sj(1,i) = trxyz(m,i-2) - 2.0*trxyz(m,i-1) + trxyz(m,i  )
		sj(2,i) = trxyz(m,i-1) - 2.0*trxyz(m,i  ) + trxyz(m,i+1) 
		sj(3,i) = trxyz(m,i  ) - 2.0*trxyz(m,i+1) + trxyz(m,i+2)

	    do m=1,nl

			g(m,1,i) = 0.5*(     qj(m,i-2) - 4.0*qj(m,i-1) + 3.0*qj(m,i  ))
			g(m,2,i) = 0.5*(     qj(m,i+1) -     qj(m,i-1)                )
			g(m,3,i) = 0.5*(-3.0*qj(m,i  ) + 4.0*qj(m,i+1) -     qj(m,i+2))

			s(m,1,i) = qj(m,i-2) - 2.0*qj(m,i-1) + qj(m,i  )
			s(m,2,i) = qj(m,i-1) - 2.0*qj(m,i  ) + qj(m,i+1) 
			s(m,3,i) = qj(m,i  ) - 2.0*qj(m,i+1) + qj(m,i+2)

			do n=1,3
				IS = g(m,n,i)*g(m,n,i) + s(m,n,i)*s(m,n,i)
				EIS2 = (EPS+IS)**2
				bl(n) = CL(n)/EIS2
				br(n) = CR(n)/EIS2
			enddo

			IS = bl(1) + bl(2) + bl(3)
			do n=1,3
				wl(m,n,i) = bl(n)/IS
			enddo

			IS = br(1) + br(2) + br(3)
			do n=1,3
				wr(m,n,i) = br(n)/IS
			enddo
      enddo


  enddo

!
!	calculating ql & qr
!

    do i=3,ni-2 
		i1=i-1
	    do m=1,nl

		    jac_l = trxyz(5,i)  + 0.125*(wl(m,1,i )*(sj(1,i )+4.0*gj(1,i )) + &
								         wl(m,2,i )*(sj(2,i )+4.0*gj(2,i )) + &
							             wl(m,3,i )*(sj(3,i )+4.0*gj(3,i )) )

			qwl(m,i)= qj(m,i) + 0.125*(wl(m,1,i )*(s(m,1,i )+4.0*g(m,1,i )) + &
								       wl(m,2,i )*(s(m,2,i )+4.0*g(m,2,i )) + &
							           wl(m,3,i )*(s(m,3,i )+4.0*g(m,3,i )) )
			if(jac_l <= 0.0)then
				write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',i
!				jac_l = 1.d-15
				jac_l   = 0.5*(trxyz(5,i) + trxyz(5,i+1))
				qwl(m,i)= 0.5*(qj(m,i) + qj(m,i+1))
			endif
			qwl(m,i) =qwl(m,i)/jac_l


		    jac_r = trxyz(5,i) + 0.125*(wr(m,1,i)*(sj(1,i)-4.0*gj(1,i)) + &
			                            wr(m,2,i)*(sj(2,i)-4.0*gj(2,i)) + &
								        wr(m,3,i)*(sj(3,i)-4.0*gj(3,i)) )  

			qwr(m,i1)= qj(m,i) + 0.125*(wr(m,1,i)*(s(m,1,i)-4.0*g(m,1,i)) + &
			                            wr(m,2,i)*(s(m,2,i)-4.0*g(m,2,i)) + &
								        wr(m,3,i)*(s(m,3,i)-4.0*g(m,3,i)) )  
			if(jac_r <= 0.0)then
				write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',i
!				jac_r = 1.d-15
				jac_r = 0.5*(trxyz(5,i) + trxyz(5,i-1))
				qwr(m,i1)= 0.5*(qj(m,i) + qj(m,i-1))
			endif
			qwr(m,i1) =qwr(m,i1)/jac_r
	    enddo
	enddo
															   
!
!	to deal with near boundary edges
!
if(.false.)then !�߽�Ҳ����Jacpobian��Ȩ
	do m=1,nl

		jac_l = (5.d0*trxyz(5,1)+15.d0*trxyz(5,2)-5.d0*trxyz(5,3)+trxyz(5,4))/16.d0
		if(jac_l <= 0.0)then
			write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',1
!			jac_l = 1.d-15
			jac_l = 0.5*(trxyz(5,1)+trxyz(5,2))
		endif
		qwl(m,1)   = (5.d0*qj(m,1)  + 15.d0*qj(m,2) -  5.d0*qj(m,3) +  qj(m,4))/16.d0
		qwl(m,1)   = qwl(m,1)/jac_l
		qwr(m,1)   = qwl(m,1)

		jac_l = (-trxyz(5,1)+9.d0*trxyz(5,2)+9.d0*trxyz(5,3)-trxyz(5,4))/16.d0
		if(jac_l <= 0.0)then
			write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',2
!			jac_l = 1.d-15
			jac_l = 0.5*(trxyz(5,2)+trxyz(5,3))
		endif
		qwl(m,2)   = (    -qj(m,1)  +  9.d0*qj(m,2) +  9.d0*qj(m,3) -  qj(m,4))/16.d0
		qwl(m,2)   = qwl(m,2)/jac_l

		jac_r = (5.d0*trxyz(5,ni)+15.d0*trxyz(5,ni-1)-5.d0*trxyz(5,ni-2)+trxyz(5,ni-3))/16.d0
		if(jac_r <= 0.0)then
			write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',ni-1
!			jac_r = 1.d-15
			jac_r = 0.5*(trxyz(5,ni)+trxyz(5,ni-1))
		endif
		qwr(m,ni-1)= (5.d0*qj(m,ni) + 15.d0*qj(m,ni-1) - 5.d0*qj(m,ni-2) + qj(m,ni-3))/16.d0
		qwr(m,ni-1)= qwr(m,ni-1)/jac_r
		qwl(m,ni-1)= qwr(m,ni-1)

		jac_r = (-trxyz(5,ni)+9.d0*trxyz(5,ni-1)+9.d0*trxyz(5,ni-2)-trxyz(5,ni-3))/16.d0
		if(jac_r <= 0.0)then
			write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',ni-2
!			jac_r = 1.d-15
			jac_r = 0.5*(trxyz(5,ni-2)+trxyz(5,ni-1))
		endif
		qwr(m,ni-2)= (    -qj(m,ni) +  9.d0*qj(m,ni-1) + 9.d0*qj(m,ni-2) - qj(m,ni-3))/16.d0
		qwr(m,ni-2)= qwr(m,ni-2)/jac_r


		jac_l = (35.d0*trxyz(5,1)-35.d0*trxyz(5,2)+21.d0*trxyz(5,3)-5.d0*trxyz(5,4))/16.d0
		if(jac_l <= 0.0)then
			write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',0
!			jac_l = 1.d-15
			jac_l = trxyz(5,1)
		endif
		qwl(m,0) = (35.d0*qj(m,1) - 35.d0*qj(m,2) + 21.d0*qj(m,3) - 5.d0*qj(m,4))/16.d0  !cic
		qwl(m,0) = qwl(m,0)/jac_l


		jac_r = (35.d0*trxyz(5,ni)-35.d0*trxyz(5,ni-1)+21.d0*trxyz(5,ni-1)-5.d0*trxyz(5,ni-3))/16.d0
		if(jac_r <= 0.0)then
			write(*,*)' WCNS_E_5_J, �ſ˱Ȳ�ֵ������i=',ni
!		    jac_r = 1.d-15
		    jac_r = trxyz(5,ni)
		endif
		qwr(m,ni)= (35.d0*qj(m,ni)-35.d0*qj(m,ni-1)+21.d0*qj(m,ni-2)-5.d0*qj(m,ni-3))/16.d0 !cic
		qwr(m,ni) = qwr(m,ni)/jac_r

		qwr(m,0) = qwl(m,0)
		qwl(m,ni)= qwr(m,ni)

!		qwl(m,0) = (15.d0*q(m,1,1) - 10.d0*q(m,2,1   ) + 3.d0*q(m,3,1)    )/8.d0  !*TGH. NEW WCNS_E_5_BORDER
!		qwr(m,ni)= (15.d0*q(m,ni,1)- 10.d0*q(m,ni-1,1) + 3.d0*q(m,ni-2,1) )/8.d0  !*TGH. NEW WCNS_E_5_BORDER

	enddo

else	!�߽粻����Jacobian��Ȩ

	do m=1,nl

		qwl(m,1)   = (5.d0*q(m,1,1)  + 15.d0*q(m,2,1) -  5.d0*q(m,3,1) +  q(m,4,1))/16.d0
		qwl(m,2)   = (    -q(m,1,1)  +  9.d0*q(m,2,1) +  9.d0*q(m,3,1) -  q(m,4,1))/16.d0
		qwr(m,1)   = qwl(m,1)

		qwr(m,ni-1)= (5.d0*q(m,ni,1) + 15.d0*q(m,ni-1,1) - 5.d0*q(m,ni-2,1) + q(m,ni-3,1))/16.d0
		qwr(m,ni-2)= (    -q(m,ni,1) +  9.d0*q(m,ni-1,1) + 9.d0*q(m,ni-2,1) - q(m,ni-3,1))/16.d0
		qwl(m,ni-1)= qwr(m,ni-1)


		qwl(m,0) = (35.d0*q(m,1,1) - 35.d0*q(m,2,1) + 21.d0*q(m,3,1) - 5.d0*q(m,4,1))/16.d0  !cic
		qwr(m,ni)= (35.d0*q(m,ni,1)-35.d0*q(m,ni-1,1)+21.d0*q(m,ni-2,1)-5.d0*q(m,ni-3,1))/16.d0 !cic

!		qwl(m,0) = (15.d0*q(m,1,1) - 10.d0*q(m,2,1   ) + 3.d0*q(m,3,1)    )/8.d0  !*TGH. NEW WCNS_E_5_BORDER
!		qwr(m,ni)= (15.d0*q(m,ni,1)- 10.d0*q(m,ni-1,1) + 3.d0*q(m,ni-2,1) )/8.d0  !*TGH. NEW WCNS_E_5_BORDER

		qwr(m,0) = qwl(m,0)
		qwl(m,ni)= qwr(m,ni)

	enddo

endif

!____________________
!  do i=1,ni-1
  do i=0,ni   !___cic

		do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			write(*,*)'interplation failed on the left side!',i
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		endif

		if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			write(*,*)'interplation failed on the right side!',i
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		endif

		kx = nx(1,i)
		ky = nx(2,i)
		kz = nx(3,i)
		kt = nx(4,i)

		gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))    

		call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		do m=1,nl
			f(m,i) = flr(m)
		enddo

	enddo

	call FLUX_DXYZ(nl,nmax,ni,f,df)

	return
end subroutine WCNS_E_5_J
!=============================================================================!


!=============================================================================!
subroutine WCNS_E_5_41(cc1,cc2,ni,nmax,nl,method,cmethd,flux_type, &
                    limiter,efix,trxyz,q,df)
    use define_precision_mod
    use global_const,only:small,rmin_limit,pmin_limit
	implicit none
!-----------------------------------------------------------------------------!
!	�߽罵Ϊ1�׾���                                                           !
!	By                                                                        !
!			TU Guohua                                                         !
!	Date                                                                      !
!			2009.6                                                            !
!-----------------------------------------------------------------------------!
	real,parameter:: CL1=1.0/16.0,CL2=10.0/16.0,CL3=5.0/16.0
	real,parameter:: CR1=CL3     ,CR2=CL2      ,CR3=CL1
	real,parameter:: EPS=1.0e-6
	real,external :: limiter,minmod
	external flux_type
	integer :: ni,nmax,nl,method,cmethd
	real  :: cc1,cc2,efix,trxyz(5,nmax)
	real    :: q(1:nl+1,-2:nmax+3,2),df(nl,nmax)
	real    :: g(nl,3,-1:ni+1) ,s(nl,3,-1:ni+1) ,bl(3),br(3)
	real    :: wl(nl,3,-1:ni+2),wr(nl,3,-1:ni+2),CL(3),CR(3)
	real    :: ql(nl),qr(nl),gamaeq
	real    :: f(nl,0:nmax),flr(nl),qwl(nl,0:nmax),qwr(nl,0:nmax)
	real    :: EIS2,IS,nx(5,0:nmax),kx,ky,kz,kt
	integer :: i,m,n,i1,ist,ied
	integer :: nerror
	real :: u_l_old(nl,8),u_r_old(nl,8)  !���ڱ߽總������ʱʹ��

	CL(1)=CL1;CL(2)=CL2;CL(3)=CL3
	CR(1)=CR1;CR(2)=CR2;CR(3)=CR3
!
! average geometry derivative by 4th order
!
	call VALUE_HALF_NODE(5,nmax,ni,trxyz,nx)

	if(q(1,-2,1) < small)then
		ist = 2
	else
		ist =0
	endif

	if(q(1,ni+3,1) < small)then
		ied = ni-2
	else
		ied = ni
	endif

!
!	calculating three 1th and 2th derivatives & corresponding weighted constants
!
	do m=1,nl
		i = ist-1
	    s(m,2,i) = q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1)
		s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)

		do i=ist,ied+1

			g(m,1,i) = 0.5*(     q(m,i-2,1) - 4.0*q(m,i-1,1) + 3.0*q(m,i,  1))
			g(m,2,i) = 0.5*(     q(m,i+1,1) -     q(m,i-1,1)                 )
			g(m,3,i) = 0.5*(-3.0*q(m,i,  1) + 4.0*q(m,i+1,1) -     q(m,i+2,1))

			s(m,1,i) = s(m,2,i-1) !q(m,i-2,1) - 2.0*q(m,i-1,1) + q(m,i,  1)
			s(m,2,i) = s(m,3,i-1) !q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1) 
			s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)

			do n=1,3
				IS = g(m,n,i)*g(m,n,i) + s(m,n,i)*s(m,n,i)
				EIS2 = (EPS+IS)**2
				bl(n) = CL(n)/EIS2
				br(n) = CR(n)/EIS2
			enddo

			IS = bl(1) + bl(2) + bl(3)
			do n=1,3
				wl(m,n,i) = bl(n)/IS
			enddo

			IS = br(1) + br(2) + br(3)
			do n=1,3
				wr(m,n,i) = br(n)/IS
			enddo

        enddo
	enddo

!
!	calculating ql & qr
!

	do i=ist,ied
		i1 = i+1
		do m=1,nl
			qwl(m,i) = q(m,i ,1) + 0.125*(wl(m,1,i )*(s(m,1,i )+4.0*g(m,1,i )) + &
								          wl(m,2,i )*(s(m,2,i )+4.0*g(m,2,i )) + &
								          wl(m,3,i )*(s(m,3,i )+4.0*g(m,3,i )) )
			qwr(m,i) = q(m,i1,1) + 0.125*(wr(m,1,i1)*(s(m,1,i1)-4.0*g(m,1,i1)) + &
			                              wr(m,2,i1)*(s(m,2,i1)-4.0*g(m,2,i1)) + &
							              wr(m,3,i1)*(s(m,3,i1)-4.0*g(m,3,i1)) )  
		enddo
	enddo

! *********************************************
! �����õ�����ϵ�ֵ����߽��ϲ������Ը߽ײ�ֵ 
	if( ist > 1 )then
	  do m=1,nl
		qwl(m,0) = (35.d0*q(m,1,1) - 35.d0*q(m,2,1) + 21.d0*q(m,3,1) - 5.d0*q(m,4,1))/16.d0  !cic
		qwl(m,1) = (5.d0*q(m,1,1)  + 15.d0*q(m,2,1) -  5.d0*q(m,3,1) +  q(m,4,1))/16.d0
		qwl(m,2) = (    -q(m,1,1)  +  9.d0*q(m,2,1) +  9.d0*q(m,3,1) -  q(m,4,1))/16.d0
		qwr(m,0) = qwl(m,0)
		qwr(m,1) = qwl(m,1)
			
	  enddo
	endif

	if( ied < ni )then
	  do m=1,nl
		qwr(m,ni)= (35.d0*q(m,ni,1)-35.d0*q(m,ni-1,1)+21.d0*q(m,ni-2,1)-5.d0*q(m,ni-3,1))/16.d0 !cic
		qwr(m,ni-1)= (5.d0*q(m,ni,1) + 15.d0*q(m,ni-1,1) - 5.d0*q(m,ni-2,1) + q(m,ni-3,1))/16.d0
		qwr(m,ni-2)= (    -q(m,ni,1) +  9.d0*q(m,ni-1,1) + 9.d0*q(m,ni-2,1) - q(m,ni-3,1))/16.d0
		qwl(m,ni)= qwr(m,ni)
		qwl(m,ni-1)= qwr(m,ni-1)
	  enddo
	endif

	nerror=0
!____________________
!  do i=1,ni-1
   do i=0,ni   !___cic

		do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the left side!',i
!			endif

			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		endif

		if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the right side!',i
!			endif
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		endif

		kx = nx(1,i)
		ky = nx(2,i)
		kz = nx(3,i)
		kt = nx(4,i)

		gamaeq = q(nl+1,i,1) !0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))    

		call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		do m=1,nl
			f(m,i) = flr(m)
		enddo

	enddo

	call FLUX_DXYZ(nl,nmax,ni,f,df)

!!!!------------------------------------------------------------!!!!!!
!!!!------------------------------------------------------------!!!!!!

! �˶γ������߽��ʽ

	if( ist > 1 )then !��߽罵�ף��������ֵ

!		�Ȱѱ߽總��ԭ��ʽ��������ҷ�����������
	    do i=1,4
		do m=1,nl
			u_l_old(m,i) = qwl(m,i)
			u_r_old(m,i) = qwr(m,i)
		enddo 
	    enddo

!		�߽�ԭ��ʽ�ݴ����
!-------------------------------------------------------

! (1)  �߽�1��������

		do m=1,nl

!------------------------------------------
! 1/2��
		  qwl(m,0) = (3.d0*q(m,1,1) - 1.d0*q(m,2,1)   - 0.*q(m,3,1)    )/2.d0
		  qwr(m,0) = qwl(m,0)
!------------------------------------------

!------------------------------------------
! 3/2��
		  qwl(m,1)   = (19.d0*q(m,1,1)  + 25.d0*q(m,2,1)    - 2.*q(m,3,1)     )/42.d0 
		  qwr(m,1)   = qwl(m,1)
!------------------------------------------

!------------------------------------------
! 5/2��
		  qwl(m,2)   = ( 4.d0*q(m,1,1) + 3.d0*q(m,2,1)  + 9.d0*q(m,3,1)    + 2.*q(m,4,1)   )/18.d0 
		  qwr(m,2)   = qwl(m,2)  
!------------------------------------------

!------------------------------------------
! 7/2��
		  qwl(m,3) =    (-2.d0*q(m,1,1) +  3.d0*q(m,2,1)  + 3.d0*q(m,3,1)    + 2.*q(m,4,1)   )/6.d0
		  qwr(m,3) =  qwl(m,3)
	    enddo

!------------------------------------------
! ��ͨ��
		do i=0,3

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

	    enddo

	    do m=1,nl
! !������һ��
		  df(m,1) = (-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
	    enddo

! (2)  �߽�2��������
	
		do m=1,nl

!------------------------------------------
! 3/2��
		  qwl(m,1)   = (3.d0*q(m,1,1)  + 6.d0*q(m,2,1)   -  q(m,3,1)    )/8.d0 
		  qwr(m,1)   = qwl(m,1)
!------------------------------------------

!------------------------------------------
! 5/2��
		  qwl(m,2)   = ( -30.d0*q(m,1,1) +145.d0*q(m,2,1)  +  29.d0*q(m,3,1)    - 8.d0*q(m,4,1)   )/136.d0 
		  qwr(m,2)   = (   2.d0*q(m,1,1) + 49.d0*q(m,2,1)  + 125.d0*q(m,3,1)    -40.d0*q(m,4,1)   )/136.d0 
!------------------------------------------

!------------------------------------------
! 7/2�� 
		  qwl(m,3) = (-q(m,2,1) + 6.d0*q(m,3,1) + 3.d0*q(m,4,1))/8.d0
		  qwr(m,3) =  qwl(m,3)
!------------------------------------------
! 9/2��
		  qwl(m,4) = (q(m,2,1) + 2.d0*q(m,3,1) + 5.d0*q(m,4,1))/8.d0
		  qwr(m,4) =  qwl(m,4)

!------------------------------------------
! 11/2��
		  qwl(m,5) = (q(m,2,1) + q(m,3,1) + 6.d0*q(m,4,1))/8.d0
		  qwr(m,5) =  qwl(m,5)
		enddo

!------------------------------------------
! ��ͨ��
		do i=1,5

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������2��
		  df(m,2) = (-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!	(3)  �߽�3����������
	
		do m=1,nl

!------------------------------------------
! 3/2�� 
!------------------------------------------
		  ! qwl(m,1)    = (ǰ���Ѿ����)
		  qwr(m,1)      = u_r_old(m,1)

!------------------------------------------
! 5/2��
		  qwl(m,2)   = ( -29.d0*q(m,1,1) + 170.d0*q(m,2,1) +  63.d0*q(m,3,1)    + 12.d0*q(m,4,1)   )/216.d0 
		  qwr(m,2)   = u_r_old(m,2)
!------------------------------------------

!------------------------------------------
! 7/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,3)    = (ǰ���Ѿ����)
		  qwr(m,3)      = u_r_old(m,3)

!------------------------------------------
! 9/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,4)    = (ǰ���Ѿ����)
		  qwr(m,4)      = u_r_old(m,4)

		enddo

!------------------------------------------
! ��ͨ��
		do i=1,4

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������3��
		  df(m,3) = (f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!

	if( ied < ni )then !�ұ߽罵�״������������ֵ

!		�߽��ֵ�ݴ�
		do i=1,4
		  i1 = ni-i
		  do m=1,nl
			u_l_old(m,9-i) = qwl(m,i1)
			u_r_old(m,9-i) = qwr(m,i1)
		  enddo 
		enddo
!		�߽�ԭ��ʽ�ݴ����

! (4)  �߽�N����������



		do m=1,nl

!------------------------------------------
! N + 1/2��
		  qwr(m,ni)= (3.d0*q(m,ni,1)- 1.d0*q(m,ni-1,1)- 0.*q(m,ni-2,1) )/2.d0
		  qwl(m,ni)= qwr(m,ni)
!------------------------------------------

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1)= (19.d0*q(m,ni,1) + 25.d0*q(m,ni-1,1) - 2.* q(m,ni-2,1) )/42.d0 
		  qwl(m,ni-1)= qwr(m,ni-1)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( 4.d0*q(m,ni,1)+ 3.d0*q(m,ni-1,1)+9.d0*q(m,ni-2,1) + 2.*q(m,ni-3,1))/18.d0 
	      qwl(m,ni-2)= qwr(m,ni-2)
!------------------------------------------

!------------------------------------------
! N - 5/2��
		  qwr(m,ni-3) = (-2.d0*q(m,ni,1) + 3.d0*q(m,ni-1,1)+3.d0*q(m,ni-2,1) + 2.*q(m,ni-3,1))/6.d0
		  qwl(m,ni-3) =  qwr(m,ni-3)
		enddo

!------------------------------------------
! ��ͨ��
		do i=ni,ni-3,-1

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,min(i+1,1),1) + q(m,i,1) )/2.0
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) =(q(m,min(i+1,1),1) + q(m,i,1)  )/2.0
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo
		enddo

		do m=1,nl
! ��N������
		  df(m,ni) = -(-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
		enddo

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! (5)  �߽�N-1����������
	
		do m=1,nl

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1)= (3.d0*q(m,ni,1) + 6.d0*q(m,ni-1,1) - q(m,ni-2,1) )/8.d0 
		  qwl(m,ni-1)= qwr(m,ni-1)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -30.d0*q(m,ni,1)+145.d0*q(m,ni-1,1) +29.d0*q(m,ni-2,1) - 8.d0*q(m,ni-3,1))/136.d0 
	      qwl(m,ni-2)= (   2.d0*q(m,ni,1)+ 49.d0*q(m,ni-1,1)+125.d0*q(m,ni-2,1) -40.d0*q(m,ni-3,1))/136.d0
!------------------------------------------

!------------------------------------------
! N - 5/2�� 
		  qwr(m,ni-3) = (-q(m,ni-1,1) + 6.d0*q(m,ni-2,1) + 3.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-3) =  qwr(m,ni-3)

!------------------------------------------
! N - 7/2��
		  qwr(m,ni-4) = (q(m,ni-1,1) + 2.d0*q(m,ni-2,1) + 5.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-4) =  qwr(m,ni-4)

!------------------------------------------
! N - 9/2��
		  qwr(m,ni-5) = (q(m,ni-1,1) + q(m,ni-2,1) + 6.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-5) =  qwr(m,ni-5)

		enddo

!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-5,-1
		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ������n-1��
		  df(m,ni-1) = -(-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!-----------------------------------------------------------------------
! (6)  �߽�N-2����������
	
		do m=1,nl

!------------------------------------------
! N - 1/2�� 
!------------------------------------------
		  ! qwr(m,ni-1) = (ǰ���Ѿ����)
		  qwl(m,ni-1)   = u_l_old(m,8)
!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -29.d0*q(m,ni,1)+ 170.d0*q(m,ni-1,1)+63.d0*q(m,ni-2,1) + 12.d0*q(m,ni-3,1))/216.d0 
		  qwl(m,ni-2)= u_l_old(m,7)
!------------------------------------------

!------------------------------------------
! N - 5/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-3) = (ǰ���Ѿ����)
		  qwl(m,ni-3)   = u_l_old(m,6)

!------------------------------------------
! N - 7/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-4) = (ǰ���Ѿ����)
		  qwl(m,ni-4)   = u_l_old(m,5)
		enddo

!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-4,-1

		  do m=1,nl
			qr(m) = qwr(m,i)
			ql(m) = qwl(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ! ��N-2������
		  df(m,ni-2) = -(f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif

! ����. �˶γ������߽��ʽ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(nerror > 8)then
		!write(*,*) 'WCNSE5 failed points along a line:', nerror
	endif

	return
end subroutine WCNS_E_5_41

!=============================================================================!
!=============================================================================!


!=============================================================================!
subroutine WCNS_E_5_42(cc1,cc2,ni,nmax,nl,method,cmethd,flux_type, &
                    limiter,efix,trxyz,q,df)
    use define_precision_mod
    use global_const,only:small,rmin_limit,pmin_limit
	implicit none
!-----------------------------------------------------------------------------!
!	�߽罵Ϊ2�׾���                                                           !
!	By                                                                        !
!			TU Guohua                                                         !
!	Date                                                                      !
!			2009.6                                                            !
!-----------------------------------------------------------------------------!
	real,parameter:: CL1=1.0/16.0,CL2=10.0/16.0,CL3=5.0/16.0
	real,parameter:: CR1=CL3     ,CR2=CL2      ,CR3=CL1
	real,parameter:: EPS=1.0e-6
	real,external :: limiter,minmod
	external flux_type
	integer :: ni,nmax,nl,method,cmethd
	real  :: cc1,cc2,efix,trxyz(5,nmax)
	real    :: q(1:nl+1,-2:nmax+3,2),df(nl,nmax)
	real    :: g(nl,3,-1:ni+1) ,s(nl,3,-1:ni+1) ,bl(3),br(3)
	real    :: wl(nl,3,-1:ni+2),wr(nl,3,-1:ni+2),CL(3),CR(3)
	real    :: ql(nl),qr(nl),gamaeq
	real    :: f(nl,0:nmax),flr(nl),qwl(nl,0:nmax),qwr(nl,0:nmax)
	real    :: EIS2,IS,nx(5,0:nmax),kx,ky,kz,kt
	integer :: i,m,n,i1,ist,ied
	integer :: nerror
	real :: u_l_old(nl,8),u_r_old(nl,8)  !���ڱ߽總������ʱʹ��
    common /rpind/ i

	CL(1)=CL1;CL(2)=CL2;CL(3)=CL3
	CR(1)=CR1;CR(2)=CR2;CR(3)=CR3
!
! average geometry derivative by 4th order
!
	call VALUE_HALF_NODE(5,nmax,ni,trxyz,nx)

	if(q(1,-2,1) < small)then
		ist = 2
	else
		ist =0
	endif

	if(q(1,ni+3,1) < small)then
		ied = ni-2
	else
		ied = ni
	endif

!
!	calculating three 1th and 2th derivatives & corresponding weighted constants
!
	do m=1,nl
		i = ist-1
	    s(m,2,i) = q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1)
		s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)

		do i=ist,ied+1
			g(m,1,i) = 0.5*(     q(m,i-2,1) - 4.0*q(m,i-1,1) + 3.0*q(m,i,  1))
			g(m,2,i) = 0.5*(     q(m,i+1,1) -     q(m,i-1,1)                 )
			g(m,3,i) = 0.5*(-3.0*q(m,i,  1) + 4.0*q(m,i+1,1) -     q(m,i+2,1))

			s(m,1,i) = s(m,2,i-1) !q(m,i-2,1) - 2.0*q(m,i-1,1) + q(m,i,  1)
			s(m,2,i) = s(m,3,i-1) !q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1) 
			s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)
			
			do n=1,3
				IS = g(m,n,i)*g(m,n,i) + s(m,n,i)*s(m,n,i)
				EIS2 = (EPS+IS)**2
				bl(n) = CL(n)/EIS2
				br(n) = CR(n)/EIS2
			enddo

			IS = bl(1) + bl(2) + bl(3)
			do n=1,3
				wl(m,n,i) = bl(n)/IS
			enddo

			IS = br(1) + br(2) + br(3)
			do n=1,3
				wr(m,n,i) = br(n)/IS
			enddo

        enddo
	enddo

!
!	calculating ql & qr
!

	do i=ist,ied
		i1 = i+1
		do m=1,nl
			qwl(m,i) = q(m,i ,1) + 0.125*(wl(m,1,i )*(s(m,1,i )+4.0*g(m,1,i )) + &
								          wl(m,2,i )*(s(m,2,i )+4.0*g(m,2,i )) + &
								          wl(m,3,i )*(s(m,3,i )+4.0*g(m,3,i )) )
			qwr(m,i) = q(m,i1,1) + 0.125*(wr(m,1,i1)*(s(m,1,i1)-4.0*g(m,1,i1)) + &
			                              wr(m,2,i1)*(s(m,2,i1)-4.0*g(m,2,i1)) + &
							              wr(m,3,i1)*(s(m,3,i1)-4.0*g(m,3,i1)) )
		enddo
	enddo

! *********************************************
! �����õ�����ϵ�ֵ����߽��ϲ������Ը߽ײ�ֵ 
	if( ist > 1 )then
	  do m=1,nl
		qwl(m,0) = (35.d0*q(m,1,1) - 35.d0*q(m,2,1) + 21.d0*q(m,3,1) - 5.d0*q(m,4,1))/16.d0  !cic
		qwl(m,1) = (5.d0*q(m,1,1)  + 15.d0*q(m,2,1) -  5.d0*q(m,3,1) +  q(m,4,1))/16.d0
		qwl(m,2) = (    -q(m,1,1)  +  9.d0*q(m,2,1) +  9.d0*q(m,3,1) -  q(m,4,1))/16.d0
		qwr(m,0) = qwl(m,0)
		qwr(m,1) = qwl(m,1)
			
	  enddo
	endif

	if( ied < ni )then
	  do m=1,nl
		qwr(m,ni)= (35.d0*q(m,ni,1)-35.d0*q(m,ni-1,1)+21.d0*q(m,ni-2,1)-5.d0*q(m,ni-3,1))/16.d0 !cic
		qwr(m,ni-1)= (5.d0*q(m,ni,1) + 15.d0*q(m,ni-1,1) - 5.d0*q(m,ni-2,1) + q(m,ni-3,1))/16.d0
		qwr(m,ni-2)= (    -q(m,ni,1) +  9.d0*q(m,ni-1,1) + 9.d0*q(m,ni-2,1) - q(m,ni-3,1))/16.d0
		qwl(m,ni)= qwr(m,ni)
		qwl(m,ni-1)= qwr(m,ni-1)
	  enddo
	endif

	nerror=0
!____________________
!  do i=1,ni-1
   do i=0,ni   !___cic

		do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the left side!',i
!			endif

			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		endif

		if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the right side!',i
!			endif
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		endif

		kx = nx(1,i)
		ky = nx(2,i)
		kz = nx(3,i)
		kt = nx(4,i)

		gamaeq = q(nl+1,i,1) !0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))    

		call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		do m=1,nl
			f(m,i) = flr(m)
		enddo

	enddo

	call FLUX_DXYZ(nl,nmax,ni,f,df)

!!!!---------------------------------------------------------------!!!
!!!!---------------------------------------------------------------!!!

! �˶γ������߽��ʽ

	if( ist > 1 )then !��߽罵�ף��������ֵ

!		�Ȱѱ߽總��ԭ��ʽ��������ҷ�����������
	    do i=1,4
		do m=1,nl
			u_l_old(m,i) = qwl(m,i)
			u_r_old(m,i) = qwr(m,i)
		enddo 
	    enddo

!		�߽�ԭ��ʽ�ݴ����
!-------------------------------------------------------

! (1)  �߽�1��������

		do m=1,nl

!------------------------------------------
! 1/2��
		  qwl(m,0) = (677.d0*q(m,1,1) - 421.d0*q(m,2,1)   + 99.d0*q(m,3,1)   + 13.d0*q(m,4,1)    )/368.d0
		  qwr(m,0) = qwl(m,0)
!------------------------------------------

!------------------------------------------
! 3/2��
		qwl(m,1)   = (5.d0*q(m,1,1)  + 15.d0*q(m,2,1)    - 5.d0*q(m,3,1)   + q(m,4,1)    )/16.d0 
		qwr(m,1)   = qwl(m,1)
!------------------------------------------

!------------------------------------------
! 5/2��
		qwl(m,2)   = ( -q(m,1,1) + 9.d0*q(m,2,1)  + 9.d0*q(m,3,1)    - q(m,4,1)   )/16.d0 
		qwr(m,2)   = qwl(m,2)  
!------------------------------------------

!------------------------------------------
! 7/2��
		qwl(m,3)   = (q(m,1,1) - 5.d0*q(m,2,1)   + 15.d0*q(m,3,1)    + 5.d0*q(m,4,1)    )/16.d0
		qwr(m,3)   =  qwl(m,3)
	    enddo

!------------------------------------------
! ��ͨ��
		do i=0,3

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

	    enddo

	    do m=1,nl
! !������һ��
		  df(m,1) = (-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
	    enddo

! (2)  �߽�2��������
	
		do m=1,nl

!------------------------------------------
! 3/2��
		  qwl(m,1)   = (3.d0*q(m,1,1)  + 6.d0*q(m,2,1)   -  q(m,3,1)    )/8.d0 
		  qwr(m,1)   = qwl(m,1)
!------------------------------------------

!------------------------------------------
! 5/2��
		  qwl(m,2)   = ( -30.d0*q(m,1,1) +145.d0*q(m,2,1)  +  29.d0*q(m,3,1)    - 8.d0*q(m,4,1)   )/136.d0 
		  qwr(m,2)   = (   2.d0*q(m,1,1) + 49.d0*q(m,2,1)  + 125.d0*q(m,3,1)    -40.d0*q(m,4,1)   )/136.d0 
!------------------------------------------

!------------------------------------------
! 7/2�� 
		  qwl(m,3) = (-q(m,2,1) + 6.d0*q(m,3,1) + 3.d0*q(m,4,1))/8.d0
		  qwr(m,3) =  qwl(m,3)
!------------------------------------------
! 9/2��
		  qwl(m,4) = (q(m,2,1) + 2.d0*q(m,3,1) + 5.d0*q(m,4,1))/8.d0
		  qwr(m,4) =  qwl(m,4)

!------------------------------------------
! 11/2��
		  qwl(m,5) = (q(m,2,1) + q(m,3,1) + 6.d0*q(m,4,1))/8.d0
		  qwr(m,5) =  qwl(m,5)
		enddo

!------------------------------------------
! ��ͨ��
		do i=1,5

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������2��
		  df(m,2) = (-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!	(3)  �߽�3����������
	
		do m=1,nl

!------------------------------------------
! 3/2�� 
!------------------------------------------
		  ! qwl(m,1)    = (ǰ���Ѿ����)
		  qwr(m,1)      = u_r_old(m,1)

!------------------------------------------
! 5/2��
		  qwl(m,2)   = ( -29.d0*q(m,1,1) + 170.d0*q(m,2,1) +  63.d0*q(m,3,1)    + 12.d0*q(m,4,1)   )/216.d0 
		  qwr(m,2)   = u_r_old(m,2)
!------------------------------------------

!------------------------------------------
! 7/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,3)    = (ǰ���Ѿ����)
		  qwr(m,3)      = u_r_old(m,3)

!------------------------------------------
! 9/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,4)    = (ǰ���Ѿ����)
		  qwr(m,4)      = u_r_old(m,4)

		enddo

!------------------------------------------
! ��ͨ��
		do i=1,4

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������3��
		  df(m,3) = (f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!

	if( ied < ni )then !�ұ߽罵�״������������ֵ

!		�߽��ֵ�ݴ�
		do i=1,4
		  i1 = ni-i
		  do m=1,nl
			u_l_old(m,9-i) = qwl(m,i1)
			u_r_old(m,9-i) = qwr(m,i1)
		  enddo 
		enddo
!		�߽�ԭ��ʽ�ݴ����

! (4)  �߽�N����������



		do m=1,nl

!------------------------------------------
! N + 1/2��
		  qwr(m,ni)= (677.d0*q(m,ni,1)- 421.d0*q(m,ni-1,1)+ 99.d0*q(m,ni-2,1)+ 13.d0*q(m,ni-3,1) )/368.d0
		  qwl(m,ni)= qwr(m,ni)
!------------------------------------------

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1)= (5.d0*q(m,ni,1) + 15.d0*q(m,ni-1,1) - 5.d0*q(m,ni-2,1)+ q(m,ni-3,1) )/16.d0 
		  qwl(m,ni-1)= qwr(m,ni-1)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -q(m,ni,1)+ 9.d0*q(m,ni-1,1)+9.d0*q(m,ni-2,1) - q(m,ni-3,1))/16.d0 
	      qwl(m,ni-2)= qwr(m,ni-2)
!------------------------------------------

!------------------------------------------
! N - 5/2��
		  qwr(m,ni-3)= (q(m,ni,1)- 5.d0*q(m,ni-1,1)+ 15.d0*q(m,ni-2,1) + 5.d0*q(m,ni-3,1) )/16.d0
		  qwl(m,ni-3)=  qwr(m,ni-3)
		enddo

!------------------------------------------
! ��ͨ��
		do i=ni,ni-3,-1

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,min(i+1,1),1) + q(m,i,1) )/2.0
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) =(q(m,min(i+1,1),1) + q(m,i,1)  )/2.0
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo
		enddo

		do m=1,nl
! ��N������
		  df(m,ni) = -(-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
		enddo

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! (5)  �߽�N-1����������
	
		do m=1,nl

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1)= (3.d0*q(m,ni,1) + 6.d0*q(m,ni-1,1) - q(m,ni-2,1) )/8.d0 
		  qwl(m,ni-1)= qwr(m,ni-1)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -30.d0*q(m,ni,1)+145.d0*q(m,ni-1,1) +29.d0*q(m,ni-2,1) - 8.d0*q(m,ni-3,1))/136.d0 
	      qwl(m,ni-2)= (   2.d0*q(m,ni,1)+ 49.d0*q(m,ni-1,1)+125.d0*q(m,ni-2,1) -40.d0*q(m,ni-3,1))/136.d0
!------------------------------------------

!------------------------------------------
! N - 5/2�� 
		  qwr(m,ni-3) = (-q(m,ni-1,1) + 6.d0*q(m,ni-2,1) + 3.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-3) =  qwr(m,ni-3)

!------------------------------------------
! N - 7/2��
		  qwr(m,ni-4) = (q(m,ni-1,1) + 2.d0*q(m,ni-2,1) + 5.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-4) =  qwr(m,ni-4)

!------------------------------------------
! N - 9/2��
		  qwr(m,ni-5) = (q(m,ni-1,1) + q(m,ni-2,1) + 6.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-5) =  qwr(m,ni-5)

		enddo

!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-5,-1
		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ������n-1��
		  df(m,ni-1) = -(-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!-----------------------------------------------------------------------
! (6)  �߽�N-2����������
	
		do m=1,nl

!------------------------------------------
! N - 1/2�� 
!------------------------------------------
		  ! qwr(m,ni-1) = (ǰ���Ѿ����)
		  qwl(m,ni-1)   = u_l_old(m,8)
!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -29.d0*q(m,ni,1)+ 170.d0*q(m,ni-1,1)+63.d0*q(m,ni-2,1) + 12.d0*q(m,ni-3,1))/216.d0 
		  qwl(m,ni-2)= u_l_old(m,7)
!------------------------------------------

!------------------------------------------
! N - 5/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-3) = (ǰ���Ѿ����)
		  qwl(m,ni-3)   = u_l_old(m,6)

!------------------------------------------
! N - 7/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-4) = (ǰ���Ѿ����)
		  qwl(m,ni-4)   = u_l_old(m,5)
		enddo

!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-4,-1

		  do m=1,nl
			qr(m) = qwr(m,i)
			ql(m) = qwl(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ! ��N-2������
		  df(m,ni-2) = -(f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif

! ����. �˶γ������߽��ʽ

	if(nerror > 8)then
		!write(*,*) 'WCNSE5 failed points along a line:', nerror
	endif

	return

end subroutine WCNS_E_5_42

!=============================================================================!

subroutine WCNS_E_5_45(cc1,cc2,ni,nmax,nl,method,cmethd,flux_type, &
                    limiter,efix,trxyz,q,df)
    use define_precision_mod
    use global_const,only:small,rmin_limit,pmin_limit
	implicit none
!-----------------------------------------------------------------------------!
!	�߽罵Ϊ2�׾���                                                           !
!	By                                                                        !
!			Fox Liu,TU Guohua                                                 !
!	Date                                                                      !
!			2010.1                                                            !
!-----------------------------------------------------------------------------!
	real,parameter:: CL1=1.0/16.0,CL2=10.0/16.0,CL3=5.0/16.0
	real,parameter:: CR1=CL3     ,CR2=CL2      ,CR3=CL1
	real,parameter:: EPS=1.0e-6
	real,external :: limiter,minmod
	external flux_type
	integer :: ni,nmax,nl,method,cmethd
	real  :: cc1,cc2,efix,trxyz(5,-3:nmax+4)
	real    :: q(1:nl+1,-2:nmax+3,2),df(nl,nmax)
	real    :: g(nl,3,-1:ni+1) ,s(nl,3,-1:ni+1) ,bl(3),br(3)
	real    :: wl(nl,3,-1:ni+2),wr(nl,3,-1:ni+2),CL(3),CR(3)
	real    :: ql(nl),qr(nl),gamaeq
	real    :: f(nl,0:nmax),flr(nl),qwl(nl,0:nmax,1),qwr(nl,0:nmax,1)
	real    :: EIS2,IS,kx,ky,kz,kt,nx(5,-2:nmax+2)
    real    :: cq(nl,-2:nmax+3),norm(4),qtmp(nl),qave(nl),qal(nl),qar(nl)
	integer :: i,m,n,i1,i2,ist,ied,k,kc
	integer :: nerror
	real :: u_l_old(nl,8),u_r_old(nl,8)  !���ڱ߽總������ʱʹ��

	CL(1)=CL1;CL(2)=CL2;CL(3)=CL3
	CR(1)=CR1;CR(2)=CR2;CR(3)=CR3
!
! average geometry derivative by 4th order
!
	call VALUE_HALF_NODE_45(5,nmax,ni,trxyz,nx)

	if(q(1,-2,1) < small)then
		ist = 2
	else
		ist =0
	endif

	if(q(1,ni+3,1) < small)then
		ied = ni-2
	else
		ied = ni
	endif
	
!	
!   calculating special characteristic varialbes from primary varialbes.
!
    do i=ist-2,ied+2
        i1 = i+1
        qal = q(1:nl,i ,1)
        qar = q(1:nl,i1,1)
        
        gamaeq = q(nl+1,i,1)
        
        qal(1) = abs( qal(1) )
        qar(1) = abs( qar(1) )
 		do m=1,nl
           qtmp(m) = qar(m) - qal(m)
        end do
        
        call pave(nl,qal,qar,qave,gamaeq)
        
        !!call LxQ(qave,trxyz(1:4,i),qtmp,cq(1:nl,i),gamaeq)       
        call LxQ(qave,nx(1:4,i),qtmp,cq(1:nl,i),gamaeq)       
    end do
    
!
!	calculating three 1th and 2th derivatives & corresponding weighted constants
!
	do i=ist,ied+1
		do m=1,nl
            g(m,1,i) = 0.5*( 3.0*cq(m,i-1) - cq(m,i-2) )
            g(m,2,i) = 0.5*(     cq(m,i  ) + cq(m,i-1) )             
            g(m,3,i) = 0.5*( 3.0*cq(m,i  ) - cq(m,i+1) )

            s(m,1,i) = cq(m,i-1) - cq(m,i-2)
            s(m,2,i) = cq(m,i  ) - cq(m,i-1) 
            s(m,3,i) = cq(m,i+1) - cq(m,i  )

			do n=1,3
				IS = g(m,n,i)*g(m,n,i) + s(m,n,i)*s(m,n,i)
				EIS2 = (EPS+IS)**2
				bl(n) = CL(n)/EIS2
				br(n) = CR(n)/EIS2
			enddo

			IS = bl(1) + bl(2) + bl(3)
			do n=1,3
				wl(m,n,i) = bl(n)/IS
			enddo

			IS = br(1) + br(2) + br(3)
			do n=1,3
				wr(m,n,i) = br(n)/IS
			enddo

        enddo
	enddo

	k = 1      !! very important(lhy)
	
!
!	calculating ql & qr
!
	do i=ist,ied
        i1 = i+1
        do m=1,nl
            qal(m) = 0.125*( wl(m,1,i )*(s(m,1,i )+4.0*g(m,1,i )) + &
                             wl(m,2,i )*(s(m,2,i )+4.0*g(m,2,i )) + &
                             wl(m,3,i )*(s(m,3,i )+4.0*g(m,3,i )) )
            qar(m) = 0.125*( wr(m,1,i1)*(s(m,1,i1)-4.0*g(m,1,i1)) + &
                             wr(m,2,i1)*(s(m,2,i1)-4.0*g(m,2,i1)) + &
                             wr(m,3,i1)*(s(m,3,i1)-4.0*g(m,3,i1)) )  
        enddo
        
        gamaeq = q(nl+1,i,1)
        
!!        qtmp = qal
!!        call RxQ(q(1:nl,i ,1),trxyz(1:4,i ),qtmp,qal,gamaeq)
!!        qtmp = qar
!!        call RxQ(q(1:nl,i1,1),trxyz(1:4,i1),qtmp,qar,gamaeq)

        call pave(nl,q(1:nl,i ,1),q(1:nl,i1,1),qave,gamaeq)
        qtmp = qal
        call RxQ(qave,nx(1:4,i ),qtmp,qal,gamaeq)
        qtmp = qar
        call RxQ(qave,nx(1:4,i),qtmp,qar,gamaeq)
        
        do m=1,nl
            qwl(m,i,k) = q(m,i ,k) + qal(m)
            qwr(m,i,k) = q(m,i1,k) + qar(m)  
        enddo
	enddo


! *********************************************
! �����õ�����ϵ�ֵ����߽��ϲ������Ը߽ײ�ֵ 
	if( ist > 1 )then
	  do m=1,nl
		qwl(m,0,k) = (35.d0*q(m,1,k) - 35.d0*q(m,2,k) + 21.d0*q(m,3,k) - 5.d0*q(m,4,k))/16.d0  !cic
		qwl(m,1,k) = ( 5.d0*q(m,1,k) + 15.d0*q(m,2,k) -  5.d0*q(m,3,k) +      q(m,4,k))/16.d0
		qwl(m,2,k) = (     -q(m,1,k) +  9.d0*q(m,2,k) +  9.d0*q(m,3,k) -      q(m,4,k))/16.d0
		qwr(m,0,k) = qwl(m,0,k)
		qwr(m,1,k) = qwl(m,1,k)
	  enddo
	endif

	if( ied < ni )then
	  do m=1,nl
		qwr(m,ni  ,k)= (35.d0*q(m,ni,k) - 35.d0*q(m,ni-1,k) + 21.d0*q(m,ni-2,k) - 5.d0*q(m,ni-3,k))/16.d0 !cic
		qwr(m,ni-1,k)= ( 5.d0*q(m,ni,k) + 15.d0*q(m,ni-1,k) -  5.d0*q(m,ni-2,k) +      q(m,ni-3,k))/16.d0
		qwr(m,ni-2,k)= (     -q(m,ni,k) +  9.d0*q(m,ni-1,k) +  9.d0*q(m,ni-2,k) -      q(m,ni-3,k))/16.d0
		qwl(m,ni  ,k)= qwr(m,ni  ,k)
		qwl(m,ni-1,k)= qwr(m,ni-1,k)
	  enddo
	endif
	
    

	kc = 1      !! very important(lhy)
	

	nerror=0
!____________________
!    do i=1,ni-1
    do i=0,ni   !___cic

		do m=1,nl
			ql(m) = qwl(m,i,kc)
			qr(m) = qwr(m,i,kc)
		enddo 

!
!		to prevent the advent of negtive values of density and pressure
!
		if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the left side!',i
!			endif

			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		endif

		if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the right side!',i
!			endif
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		endif

		kx = nx(1,i)
		ky = nx(2,i)
		kz = nx(3,i)
		kt = nx(4,i)

		gamaeq = q(nl+1,i,1) !0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))    

		call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		do m=1,nl
			f(m,i) = flr(m)
		enddo

	enddo

	call FLUX_DXYZ(nl,nmax,ni,f,df)

!!!!---------------------------------------------------------------!!!
!!!!---------------------------------------------------------------!!!

! �˶γ������߽��ʽ

	if( ist > 1 )then !��߽罵�ף��������ֵ

!		�Ȱѱ߽總��ԭ��ʽ��������ҷ�����������
	    do i=1,4
		do m=1,nl
			u_l_old(m,i) = qwl(m,i,kc)
			u_r_old(m,i) = qwr(m,i,kc)
		enddo 
	    enddo

!		�߽�ԭ��ʽ�ݴ����
!-------------------------------------------------------
       
! (1)  �߽�1��������

		do m=1,nl

!------------------------------------------
! 1/2��
		  qwl(m,0,kc) = ( 677.d0*q(m,1,kc) - 421.d0*q(m,2,kc) + 99.d0*q(m,3,kc) + 13.d0*q(m,4,kc) )/368.d0
		  qwr(m,0,kc) = qwl(m,0,kc)
!------------------------------------------

!------------------------------------------
! 3/2��
		  qwl(m,1,kc) = ( 5.d0*q(m,1,kc)  + 15.d0*q(m,2,kc) - 5.d0*q(m,3,kc) + q(m,4,kc) )/16.d0 
		  qwr(m,1,kc) = qwl(m,1,kc)
!------------------------------------------

!------------------------------------------
! 5/2��
		  qwl(m,2,kc) = ( -q(m,1,kc) + 9.d0*q(m,2,kc) + 9.d0*q(m,3,kc) - q(m,4,kc) )/16.d0 
		  qwr(m,2,kc) = qwl(m,2,kc)  
!------------------------------------------

!------------------------------------------
! 7/2��
		  qwl(m,3,kc) = ( q(m,1,kc) - 5.d0*q(m,2,kc) + 15.d0*q(m,3,kc) + 5.d0*q(m,4,kc) )/16.d0
		  qwr(m,3,kc) = qwl(m,3,kc)
	    enddo

!------------------------------------------
! ��ͨ��
		do i=0,3

		  do m=1,nl
			ql(m) = qwl(m,i,kc)
			qr(m) = qwr(m,i,kc)
		  enddo

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

	    enddo

	    do m=1,nl
! !������һ��
		  df(m,1) = (-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
	    enddo

! (2)  �߽�2��������
		do m=1,nl

!------------------------------------------
! 3/2��
		  qwl(m,1,kc) = ( 3.d0*q(m,1,kc) + 6.d0*q(m,2,kc) - q(m,3,kc) )/8.d0 
		  qwr(m,1,kc) = qwl(m,1,kc)
!------------------------------------------

!------------------------------------------
! 5/2��
		  qwl(m,2,kc) = (-30.d0*q(m,1,kc) + 145.d0*q(m,2,kc) +  29.d0*q(m,3,kc) -  8.d0*q(m,4,kc) )/136.d0 
		  qwr(m,2,kc) = (  2.d0*q(m,1,kc) +  49.d0*q(m,2,kc) + 125.d0*q(m,3,kc) - 40.d0*q(m,4,kc) )/136.d0 
!------------------------------------------

!------------------------------------------
! 7/2�� 
		  qwl(m,3,kc) = (-q(m,2,kc) + 6.d0*q(m,3,kc) + 3.d0*q(m,4,kc) )/8.d0
		  qwr(m,3,kc) = qwl(m,3,kc)
!------------------------------------------
! 9/2��
		  qwl(m,4,kc) = ( q(m,2,kc) + 2.d0*q(m,3,kc) + 5.d0*q(m,4,kc) )/8.d0
		  qwr(m,4,kc) = qwl(m,4,kc)

!------------------------------------------
! 11/2��
		  qwl(m,5,kc) = ( q(m,2,kc) + q(m,3,kc) + 6.d0*q(m,4,kc) )/8.d0
		  qwr(m,5,kc) = qwl(m,5,kc)
		enddo

!------------------------------------------
! ��ͨ��
		do i=1,5

		  do m=1,nl
			ql(m) = qwl(m,i,kc)
			qr(m) = qwr(m,i,kc)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������2��
		  df(m,2) = (-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!	(3)  �߽�3����������

		do m=1,nl

!------------------------------------------
! 3/2�� 
!------------------------------------------
		  ! qwl(m,1)    = (ǰ���Ѿ����)
		  qwr(m,1,kc) = u_r_old(m,1)

!------------------------------------------
! 5/2��
		  qwl(m,2,kc) = (-29.d0*q(m,1,kc) + 170.d0*q(m,2,kc) + 63.d0*q(m,3,kc) + 12.d0*q(m,4,kc) )/216.d0 
		  qwr(m,2,kc) = u_r_old(m,2)
!------------------------------------------

!------------------------------------------
! 7/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,3)    = (ǰ���Ѿ����)
		  qwr(m,3,kc) = u_r_old(m,3)

!------------------------------------------
! 9/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,4)    = (ǰ���Ѿ����)
		  qwr(m,4,kc) = u_r_old(m,4)

		enddo
		

!------------------------------------------
! ��ͨ��
		do i=1,4

		  do m=1,nl
			ql(m) = qwl(m,i,kc)
			qr(m) = qwr(m,i,kc)
		  enddo

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������3��
		  df(m,3) = (f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!

	if( ied < ni )then !�ұ߽罵�״������������ֵ

!		�߽��ֵ�ݴ�
		do i=1,4
		  i1 = ni-i
		  do m=1,nl
			u_l_old(m,9-i) = qwl(m,i1,kc)
			u_r_old(m,9-i) = qwr(m,i1,kc)
		  enddo 
		enddo
!		�߽�ԭ��ʽ�ݴ����

! (4)  �߽�N����������

		do m=1,nl

!------------------------------------------
! N + 1/2��
		  qwr(m,ni  ,kc) = ( 677.d0*q(m,ni,kc) - 421.d0*q(m,ni-1,kc) + 99.d0*q(m,ni-2,kc) + 13.d0*q(m,ni-3,kc) )/368.d0
		  qwl(m,ni  ,kc) = qwr(m,ni,kc)
!------------------------------------------

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1,kc) = ( 5.d0*q(m,ni,kc) + 15.d0*q(m,ni-1,kc) - 5.d0*q(m,ni-2,kc) + q(m,ni-3,kc) )/16.d0 
		  qwl(m,ni-1,kc) = qwr(m,ni-1,kc)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2,kc) = ( -q(m,ni,kc) + 9.d0*q(m,ni-1,kc) + 9.d0*q(m,ni-2,kc) - q(m,ni-3,kc) )/16.d0 
	      qwl(m,ni-2,kc) = qwr(m,ni-2,kc)
!------------------------------------------

!------------------------------------------
! N - 5/2��
		  qwr(m,ni-3,kc) = ( q(m,ni,kc) - 5.d0*q(m,ni-1,kc) + 15.d0*q(m,ni-2,kc) + 5.d0*q(m,ni-3,kc) )/16.d0
		  qwl(m,ni-3,kc) = qwr(m,ni-3,kc)
		enddo


!------------------------------------------
! ��ͨ��
		do i=ni,ni-3,-1

		  do m=1,nl
			ql(m) = qwl(m,i,kc)
			qr(m) = qwr(m,i,kc)
		  enddo 
		  
!
!		to prevent the advent of negtive values of density and pressure
!
		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,min(i+1,1),1) + q(m,i,1) )/2.0
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) =(q(m,min(i+1,1),1) + q(m,i,1)  )/2.0
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo
		enddo

		do m=1,nl
! ��N������
		  df(m,ni) = -(-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
		enddo

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! (5)  �߽�N-1����������
		do m=1,nl

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1,kc) = ( 3.d0*q(m,ni,kc) + 6.d0*q(m,ni-1,kc) - q(m,ni-2,kc) )/8.d0 
		  qwl(m,ni-1,kc) = qwr(m,ni-1,kc)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2,kc) = ( -30.d0*q(m,ni,kc) + 145.d0*q(m,ni-1,kc) +  29.d0*q(m,ni-2,kc) -  8.d0*q(m,ni-3,kc) )/136.d0 
	      qwl(m,ni-2,kc) = (   2.d0*q(m,ni,kc) +  49.d0*q(m,ni-1,kc) + 125.d0*q(m,ni-2,kc) - 40.d0*q(m,ni-3,kc) )/136.d0
!------------------------------------------

!------------------------------------------
! N - 5/2�� 
		  qwr(m,ni-3,kc) = ( -q(m,ni-1,kc) + 6.d0*q(m,ni-2,kc) + 3.d0*q(m,ni-3,kc) )/8.d0
		  qwl(m,ni-3,kc) = qwr(m,ni-3,kc)

!------------------------------------------
! N - 7/2��
		  qwr(m,ni-4,kc) = ( q(m,ni-1,kc) + 2.d0*q(m,ni-2,kc) + 5.d0*q(m,ni-3,kc) )/8.d0
		  qwl(m,ni-4,kc) = qwr(m,ni-4,kc)

!------------------------------------------
! N - 9/2��
		  qwr(m,ni-5,kc) = ( q(m,ni-1,kc) + q(m,ni-2,kc) + 6.d0*q(m,ni-3,kc) )/8.d0
		  qwl(m,ni-5,kc) = qwr(m,ni-5,kc)

		enddo

!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-5,-1
		  do m=1,nl
			ql(m) = qwl(m,i,kc)
			qr(m) = qwr(m,i,kc)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ������n-1��
		  df(m,ni-1) = -(-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!-----------------------------------------------------------------------
! (6)  �߽�N-2����������
		do m=1,nl

!------------------------------------------
! N - 1/2�� 
!------------------------------------------
		  ! qwr(m,ni-1) = (ǰ���Ѿ����)
		  qwl(m,ni-1,kc)   = u_l_old(m,8)
!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2,kc)= ( -29.d0*q(m,ni,kc) + 170.d0*q(m,ni-1,kc) + 63.d0*q(m,ni-2,kc) + 12.d0*q(m,ni-3,kc) )/216.d0 
		  qwl(m,ni-2,kc)= u_l_old(m,7)
!------------------------------------------

!------------------------------------------
! N - 5/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-3) = (ǰ���Ѿ����)
		  qwl(m,ni-3,kc)   = u_l_old(m,6)

!------------------------------------------
! N - 7/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-4) = (ǰ���Ѿ����)
		  qwl(m,ni-4,kc)   = u_l_old(m,5)
		enddo
		
!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-4,-1

		  do m=1,nl
			qr(m) = qwr(m,i,kc)
			ql(m) = qwl(m,i,kc)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ! ��N-2������
		  df(m,ni-2) = -(f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif

! ����. �˶γ������߽��ʽ

	if(nerror > 8)then
		!write(*,*) 'WCNSE5 failed points along a line:', nerror
	endif

	return

end subroutine WCNS_E_5_45

!=============================================================================!
subroutine WCNS_E_5_47(cc1,cc2,ni,nmax,nl,method,cmethd,flux_type, &
                    limiter,efix,trxyz,q,df)
    use define_precision_mod
    use global_const,only:small,rmin_limit,pmin_limit
	implicit none
!-----------------------------------------------------------------------------!
!	�߽罵Ϊ2�׾���                                                           !
!	By                                                                        !
!			Fox Liu,TU Guohua                                                 !
!	Date                                                                      !
!			2010.1                                                            !
!-----------------------------------------------------------------------------!
	real,parameter:: CL1=1.0/16.0,CL2=10.0/16.0,CL3=5.0/16.0
	real,parameter:: CR1=CL3     ,CR2=CL2      ,CR3=CL1
	real,parameter:: EPS=1.0e-6
	real,external :: limiter,minmod
	external flux_type
	integer :: ni,nmax,nl,method,cmethd
	real  :: cc1,cc2,efix,trxyz(5,-3:nmax+4)
	real    :: q(1:nl+1,-2:nmax+3,2),df(nl,nmax)
	real    :: g(nl,3,-1:ni+1) ,s(nl,3,-1:ni+1) ,bl(3),br(3)
	real    :: wl(nl,3,-1:ni+2),wr(nl,3,-1:ni+2),CL(3),CR(3)
	real    :: ql(nl),qr(nl),gamaeq,IS0(nl,3),gch(nl),sch(nl)
	real    :: f(nl,0:nmax),flr(nl),qwl(nl,0:nmax),qwr(nl,0:nmax)
	real    :: EIS2,IS,nx(5,-2:nmax+2),kx,ky,kz,kt
	integer :: i,m,n,i1,ist,ied
	integer :: nerror
	real :: u_l_old(nl,8),u_r_old(nl,8)  !���ڱ߽總������ʱʹ��
    common /rpind/ i

	CL(1)=CL1;CL(2)=CL2;CL(3)=CL3
	CR(1)=CR1;CR(2)=CR2;CR(3)=CR3
!
! average geometry derivative by 4th order
!
	call VALUE_HALF_NODE_45(5,nmax,ni,trxyz,nx)

	if(q(1,-2,1) < small)then
		ist = 2
	else
		ist =0
	endif

	if(q(1,ni+3,1) < small)then
		ied = ni-2
	else
		ied = ni
	endif

!
!	calculating three 1th and 2th derivatives & corresponding weighted constants
!
	i = ist-1
	do m=1,nl
	    s(m,2,i) = q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1)
		s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)
	enddo

	do i=ist,ied+1
	    do m=1,nl
			g(m,1,i) = 0.5*(     q(m,i-2,1) - 4.0*q(m,i-1,1) + 3.0*q(m,i,  1))
			g(m,2,i) = 0.5*(     q(m,i+1,1) -     q(m,i-1,1)                 )
			g(m,3,i) = 0.5*(-3.0*q(m,i,  1) + 4.0*q(m,i+1,1) -     q(m,i+2,1))

			s(m,1,i) = s(m,2,i-1) !q(m,i-2,1) - 2.0*q(m,i-1,1) + q(m,i,  1)
			s(m,2,i) = s(m,3,i-1) !q(m,i-1,1) - 2.0*q(m,i  ,1) + q(m,i+1,1) 
			s(m,3,i) = q(m,i  ,1) - 2.0*q(m,i+1,1) + q(m,i+2,1)
        enddo
	enddo

	do i=ist,ied+1
        call LxQ(q(1:nl,i,1),trxyz(1:4,i),q(1:nl,i,1),q(1:nl,i,2),q(nl+1,i,1)) 
           
	    do n=1,3
            call LxQ(q(1:nl,i,1),trxyz(1:4,i),g(1:nl,n,i),gch(1:nl),q(nl+1,i,1))    
            call LxQ(q(1:nl,i,1),trxyz(1:4,i),s(1:nl,n,i),sch(1:nl),q(nl+1,i,1))
            
            g(1:nl,n,i) = gch(1:nl)     
            s(1:nl,n,i) = sch(1:nl)     

			do m=1,nl
				IS0(m,n) = gch(m)*gch(m) + sch(m)*sch(m)
			enddo
        enddo
        
	    do m=1,nl
			do n=1,3
				EIS2 = (EPS+IS0(m,n))**2
				bl(n) = CL(n)/EIS2
				br(n) = CR(n)/EIS2
			enddo	
				
			IS = bl(1) + bl(2) + bl(3)
			do n=1,3
				wl(m,n,i) = bl(n)/IS
			enddo

			IS = br(1) + br(2) + br(3)
			do n=1,3
				wr(m,n,i) = br(n)/IS
			enddo

        enddo
	enddo

!
!	calculating ql & qr
!
	do i=ist,ied
		i1 = i+1
		do m=1,nl
			qwl(m,i) = q(m,i ,2) + 0.125*(wl(m,1,i )*(s(m,1,i )+4.0*g(m,1,i )) + &
								          wl(m,2,i )*(s(m,2,i )+4.0*g(m,2,i )) + &
								          wl(m,3,i )*(s(m,3,i )+4.0*g(m,3,i )) )
			qwr(m,i) = q(m,i1,2) + 0.125*(wr(m,1,i1)*(s(m,1,i1)-4.0*g(m,1,i1)) + &
			                              wr(m,2,i1)*(s(m,2,i1)-4.0*g(m,2,i1)) + &
							              wr(m,3,i1)*(s(m,3,i1)-4.0*g(m,3,i1)) )
							              
		enddo
		
        call RxQ(q(1:nl,i ,1),trxyz(1:4 ,i),qwl(1:nl,i),gch(1:nl),q(nl+1,i ,1))    
        call RxQ(q(1:nl,i1,1),trxyz(1:4,i1),qwr(1:nl,i),sch(1:nl),q(nl+1,i1,1))
        qwl(1:nl,i) = gch(1:nl)
        qwr(1:nl,i) = sch(1:nl)
	enddo

! *********************************************
! �����õ�����ϵ�ֵ����߽��ϲ������Ը߽ײ�ֵ 
	if( ist > 1 )then
	  do m=1,nl
		qwl(m,0) = (35.d0*q(m,1,1) - 35.d0*q(m,2,1) + 21.d0*q(m,3,1) - 5.d0*q(m,4,1))/16.d0  !cic
		qwl(m,1) = (5.d0*q(m,1,1)  + 15.d0*q(m,2,1) -  5.d0*q(m,3,1) +  q(m,4,1))/16.d0
		qwl(m,2) = (    -q(m,1,1)  +  9.d0*q(m,2,1) +  9.d0*q(m,3,1) -  q(m,4,1))/16.d0
		qwr(m,0) = qwl(m,0)
		qwr(m,1) = qwl(m,1)
			
	  enddo
	endif

	if( ied < ni )then
	  do m=1,nl
		qwr(m,ni)= (35.d0*q(m,ni,1)-35.d0*q(m,ni-1,1)+21.d0*q(m,ni-2,1)-5.d0*q(m,ni-3,1))/16.d0 !cic
		qwr(m,ni-1)= (5.d0*q(m,ni,1) + 15.d0*q(m,ni-1,1) - 5.d0*q(m,ni-2,1) + q(m,ni-3,1))/16.d0
		qwr(m,ni-2)= (    -q(m,ni,1) +  9.d0*q(m,ni-1,1) + 9.d0*q(m,ni-2,1) - q(m,ni-3,1))/16.d0
		qwl(m,ni)= qwr(m,ni)
		qwl(m,ni-1)= qwr(m,ni-1)
	  enddo
	endif

	nerror=0
!____________________
!  do i=1,ni-1
   do i=0,ni   !___cic

		do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the left side!',i
!			endif

			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		endif

		if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror=nerror+1
!		    if(nerror <= 3)then
!			   write(*,*)'interplation failed on the right side!',i
!			endif
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		endif

		kx = nx(1,i)
		ky = nx(2,i)
		kz = nx(3,i)
		kt = nx(4,i)

		gamaeq = q(nl+1,i,1) !0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))    

		call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		do m=1,nl
			f(m,i) = flr(m)
		enddo

	enddo

	call FLUX_DXYZ(nl,nmax,ni,f,df)

!!!!---------------------------------------------------------------!!!
!!!!---------------------------------------------------------------!!!

! �˶γ������߽��ʽ

	if( ist > 1 )then !��߽罵�ף��������ֵ

!		�Ȱѱ߽總��ԭ��ʽ��������ҷ�����������
	    do i=1,4
		do m=1,nl
			u_l_old(m,i) = qwl(m,i)
			u_r_old(m,i) = qwr(m,i)
		enddo 
	    enddo

!		�߽�ԭ��ʽ�ݴ����
!-------------------------------------------------------

! (1)  �߽�1��������

		do m=1,nl

!------------------------------------------
! 1/2��
		  qwl(m,0) = (677.d0*q(m,1,1) - 421.d0*q(m,2,1)   + 99.d0*q(m,3,1)   + 13.d0*q(m,4,1)    )/368.d0
		  qwr(m,0) = qwl(m,0)
!------------------------------------------

!------------------------------------------
! 3/2��
		qwl(m,1)   = (5.d0*q(m,1,1)  + 15.d0*q(m,2,1)    - 5.d0*q(m,3,1)   + q(m,4,1)    )/16.d0 
		qwr(m,1)   = qwl(m,1)
!------------------------------------------

!------------------------------------------
! 5/2��
		qwl(m,2)   = ( -q(m,1,1) + 9.d0*q(m,2,1)  + 9.d0*q(m,3,1)    - q(m,4,1)   )/16.d0 
		qwr(m,2)   = qwl(m,2)  
!------------------------------------------

!------------------------------------------
! 7/2��
		qwl(m,3)   = (q(m,1,1) - 5.d0*q(m,2,1)   + 15.d0*q(m,3,1)    + 5.d0*q(m,4,1)    )/16.d0
		qwr(m,3)   =  qwl(m,3)
	    enddo

!------------------------------------------
! ��ͨ��
		do i=0,3

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = (q(m,max(i,1),1)+q(m,i+1,1))/2.
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

	    enddo

	    do m=1,nl
! !������һ��
		  df(m,1) = (-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
	    enddo

! (2)  �߽�2��������
	
		do m=1,nl

!------------------------------------------
! 3/2��
		  qwl(m,1)   = (3.d0*q(m,1,1)  + 6.d0*q(m,2,1)   -  q(m,3,1)    )/8.d0 
		  qwr(m,1)   = qwl(m,1)
!------------------------------------------

!------------------------------------------
! 5/2��
		  qwl(m,2)   = ( -30.d0*q(m,1,1) +145.d0*q(m,2,1)  +  29.d0*q(m,3,1)    - 8.d0*q(m,4,1)   )/136.d0 
		  qwr(m,2)   = (   2.d0*q(m,1,1) + 49.d0*q(m,2,1)  + 125.d0*q(m,3,1)    -40.d0*q(m,4,1)   )/136.d0 
!------------------------------------------

!------------------------------------------
! 7/2�� 
		  qwl(m,3) = (-q(m,2,1) + 6.d0*q(m,3,1) + 3.d0*q(m,4,1))/8.d0
		  qwr(m,3) =  qwl(m,3)
!------------------------------------------
! 9/2��
		  qwl(m,4) = (q(m,2,1) + 2.d0*q(m,3,1) + 5.d0*q(m,4,1))/8.d0
		  qwr(m,4) =  qwl(m,4)

!------------------------------------------
! 11/2��
		  qwl(m,5) = (q(m,2,1) + q(m,3,1) + 6.d0*q(m,4,1))/8.d0
		  qwr(m,5) =  qwl(m,5)
		enddo

!------------------------------------------
! ��ͨ��
		do i=1,5

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������2��
		  df(m,2) = (-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!	(3)  �߽�3����������
	
		do m=1,nl

!------------------------------------------
! 3/2�� 
!------------------------------------------
		  ! qwl(m,1)    = (ǰ���Ѿ����)
		  qwr(m,1)      = u_r_old(m,1)

!------------------------------------------
! 5/2��
		  qwl(m,2)   = ( -29.d0*q(m,1,1) + 170.d0*q(m,2,1) +  63.d0*q(m,3,1)    + 12.d0*q(m,4,1)   )/216.d0 
		  qwr(m,2)   = u_r_old(m,2)
!------------------------------------------

!------------------------------------------
! 7/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,3)    = (ǰ���Ѿ����)
		  qwr(m,3)      = u_r_old(m,3)

!------------------------------------------
! 9/2�� (ǰ���Ѿ����)
!------------------------------------------
		  ! qwl(m,4)    = (ǰ���Ѿ����)
		  qwr(m,4)      = u_r_old(m,4)

		enddo

!------------------------------------------
! ��ͨ��
		do i=1,4

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! !������3��
		  df(m,3) = (f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!
!!!!!!!!--------------------------------------------------------------!!!!!!!!!!!!

	if( ied < ni )then !�ұ߽罵�״������������ֵ

!		�߽��ֵ�ݴ�
		do i=1,4
		  i1 = ni-i
		  do m=1,nl
			u_l_old(m,9-i) = qwl(m,i1)
			u_r_old(m,9-i) = qwr(m,i1)
		  enddo 
		enddo
!		�߽�ԭ��ʽ�ݴ����

! (4)  �߽�N����������



		do m=1,nl

!------------------------------------------
! N + 1/2��
		  qwr(m,ni)= (677.d0*q(m,ni,1)- 421.d0*q(m,ni-1,1)+ 99.d0*q(m,ni-2,1)+ 13.d0*q(m,ni-3,1) )/368.d0
		  qwl(m,ni)= qwr(m,ni)
!------------------------------------------

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1)= (5.d0*q(m,ni,1) + 15.d0*q(m,ni-1,1) - 5.d0*q(m,ni-2,1)+ q(m,ni-3,1) )/16.d0 
		  qwl(m,ni-1)= qwr(m,ni-1)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -q(m,ni,1)+ 9.d0*q(m,ni-1,1)+9.d0*q(m,ni-2,1) - q(m,ni-3,1))/16.d0 
	      qwl(m,ni-2)= qwr(m,ni-2)
!------------------------------------------

!------------------------------------------
! N - 5/2��
		  qwr(m,ni-3)= (q(m,ni,1)- 5.d0*q(m,ni-1,1)+ 15.d0*q(m,ni-2,1) + 5.d0*q(m,ni-3,1) )/16.d0
		  qwl(m,ni-3)=  qwr(m,ni-3)
		enddo

!------------------------------------------
! ��ͨ��
		do i=ni,ni-3,-1

		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 
!
!		to prevent the advent of negtive values of density and pressure
!
		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = (q(m,min(i+1,1),1) + q(m,i,1) )/2.0
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) =(q(m,min(i+1,1),1) + q(m,i,1)  )/2.0
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo
		enddo

		do m=1,nl
! ��N������
		  df(m,ni) = -(-23.d0*f(m,0) + 21.d0*f(m,1) + 3.d0*f(m,2) - 1.d0*f(m,3))/24.d0 
		enddo

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! (5)  �߽�N-1����������
	
		do m=1,nl

!------------------------------------------
! N - 1/2��
		  qwr(m,ni-1)= (3.d0*q(m,ni,1) + 6.d0*q(m,ni-1,1) - q(m,ni-2,1) )/8.d0 
		  qwl(m,ni-1)= qwr(m,ni-1)
!------------------------------------------

!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -30.d0*q(m,ni,1)+145.d0*q(m,ni-1,1) +29.d0*q(m,ni-2,1) - 8.d0*q(m,ni-3,1))/136.d0 
	      qwl(m,ni-2)= (   2.d0*q(m,ni,1)+ 49.d0*q(m,ni-1,1)+125.d0*q(m,ni-2,1) -40.d0*q(m,ni-3,1))/136.d0
!------------------------------------------

!------------------------------------------
! N - 5/2�� 
		  qwr(m,ni-3) = (-q(m,ni-1,1) + 6.d0*q(m,ni-2,1) + 3.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-3) =  qwr(m,ni-3)

!------------------------------------------
! N - 7/2��
		  qwr(m,ni-4) = (q(m,ni-1,1) + 2.d0*q(m,ni-2,1) + 5.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-4) =  qwr(m,ni-4)

!------------------------------------------
! N - 9/2��
		  qwr(m,ni-5) = (q(m,ni-1,1) + q(m,ni-2,1) + 6.d0*q(m,ni-3,1))/8.d0
		  qwl(m,ni-5) =  qwr(m,ni-5)

		enddo

!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-5,-1
		  do m=1,nl
			ql(m) = qwl(m,i)
			qr(m) = qwr(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
!				ql(m) = q(m,i,1)
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
!				qr(m) = q(m,i+1,1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ������n-1��
		  df(m,ni-1) = -(-22.d0*f(m,1) + 17.d0*f(m,2) + 9.d0*f(m,3) - 5.d0*f(m,4) + f(m,5))/24.d0 
		enddo

!-----------------------------------------------------------------------
! (6)  �߽�N-2����������
	
		do m=1,nl

!------------------------------------------
! N - 1/2�� 
!------------------------------------------
		  ! qwr(m,ni-1) = (ǰ���Ѿ����)
		  qwl(m,ni-1)   = u_l_old(m,8)
!------------------------------------------
! N - 3/2��
		  qwr(m,ni-2)= ( -29.d0*q(m,ni,1)+ 170.d0*q(m,ni-1,1)+63.d0*q(m,ni-2,1) + 12.d0*q(m,ni-3,1))/216.d0 
		  qwl(m,ni-2)= u_l_old(m,7)
!------------------------------------------

!------------------------------------------
! N - 5/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-3) = (ǰ���Ѿ����)
		  qwl(m,ni-3)   = u_l_old(m,6)

!------------------------------------------
! N - 7/2�� (ǰ���Ѿ����)
!------------------------------------------

		  ! qwr(m,ni-4) = (ǰ���Ѿ����)
		  qwl(m,ni-4)   = u_l_old(m,5)
		enddo

!------------------------------------------
! ��ͨ��
		do i=ni-1,ni-4,-1

		  do m=1,nl
			qr(m) = qwr(m,i)
			ql(m) = qwl(m,i)
		  enddo 

		  if(ql(1)<=rmin_limit .or. ql(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				ql(m) = q(m,max(i,1),1)
			enddo
		  endif

		  if(qr(1)<=rmin_limit .or. qr(5)<=pmin_limit)then
			nerror = nerror + 1
			do m=1,nl
				qr(m) = q(m,min(i+1,ni),1)
			enddo
		  endif

		  kx = nx(1,i)
		  ky = nx(2,i)
		  kz = nx(3,i)
		  kt = nx(4,i)

		  gamaeq = 0.5*(q(nl+1,i,1) + q(nl+1,i  ,1))

		  call flux_type(ql,qr,nl,kx,ky,kz,kt,flr,efix,gamaeq,gamaeq)

		  do m=1,nl
			f(m,ni-i) = flr(m)
		  enddo

		enddo

		do m=1,nl
! ! ��N-2������
		  df(m,ni-2) = -(f(m,1) - 27.d0*f(m,2) + 27.d0*f(m,3) - f(m,4) )/24.d0 
		enddo
	endif

! ����. �˶γ������߽��ʽ

	if(nerror > 8)then
		!write(*,*) 'WCNSE5 failed points along a line:', nerror
	endif

	return

end subroutine WCNS_E_5_47

!================================================================!

subroutine pave(nl,ql,qr,qa,gama)
    implicit none
    integer :: nl
    real :: ql(nl),qr(nl),qa(nl),gama
    real :: rl,ul,vl,wl,pl,hl,vl2
    real :: rr,ur,vr,wr,pr,hr,vr2
    real :: gam1,gam2,rls,rrs,oors
    real :: ra,ua,va,wa,ha,pa,va2
    
    gam1 = gama-1.0
    gam2 = gama/gam1

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

	vl2 = ul*ul + vl*vl + wl*wl
	vr2 = ur*ur + vr*vr + wr*wr

    hl = gam2*pl/rl + 0.5*vl2
    hr = gam2*pr/rr + 0.5*vr2

    rls = sqrt(rl)
    rrs = sqrt(rr)
    ra = rls*rrs

    oors = 1.0/(rls + rrs)
    rls = rls*oors
    rrs = rrs*oors

    ua = rls*ul + rrs*ur
    va = rls*vl + rrs*vr
    wa = rls*wl + rrs*wr
    ha = rls*hl + rrs*hr
    
    va2 = ua*ua + va*va + wa*wa
    pa = (ha - 0.5*va2)*ra/gam2
    
    qa(1) = ra
    qa(2) = ua
    qa(3) = va
    qa(4) = wa
    qa(5) = pa
   
!!    qa(:) = 0.5*(ql(:) + qr(:))
     
end subroutine pave


subroutine p2q(nl,prim,q,gama)
    implicit none
    integer :: nl
    real :: prim(nl),q(nl)
    real :: rm,um,vm,wm,pm,em,gama
  
    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)
    
    em = pm/(gama-1.0) + 0.5*rm*( um*um + vm*vm + wm*wm )

    q(1) = rm
    q(2) = rm * um
    q(3) = rm * vm
    q(4) = rm * wm
    q(5) = em
     
end subroutine p2q

subroutine q2p(nl,q,prim,gama)
    implicit none
    integer :: nl
    real :: prim(nl),q(nl)
    real :: rm,oor,um,vm,wm,pm,em,gama
  
    rm = q(1)
    oor = 1.0/rm
    um = q(2)*oor
    vm = q(3)*oor
    wm = q(4)*oor
    
    em = q(5)
    pm = (gama-1.0) * ( em - 0.5*rm*( um*um + vm*vm + wm*wm ) )

    prim(1) = rm
    prim(2) = um
    prim(3) = vm
    prim(4) = wm
    prim(5) = pm
    
end subroutine q2p

subroutine LxQ(prim,gykb,dq,f,gama)
    use global_variables,only:nl,sml_sss
    implicit none
    real :: prim(nl),gykb(4),dq(nl),f(nl),gama
    real :: nx,ny,nz,nt,oon
    real :: rm,um,vm,wm,pm,c2,cm
    real :: dqn,wp1,wp2
    
    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)
    
    oon = 1.0/max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
    
    nx = nx*oon
    ny = ny*oon
    nz = nz*oon
    nt = nt*oon

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    c2 = gama*pm/rm
    if ( c2 < 0.0 ) then
       write(*,*)'LxQ(prim,gykb,dq,f)'
       write(*,*)'c2<0,gama,pm,rm',gama,pm,rm
    endif

    cm = sqrt(c2)
    
    dqn = nx*dq(2) + ny*dq(3) + nz*dq(4)
    wp1 = dq(1) - dq(5)/c2
    wp2 = dq(5)/(rm*cm)

    f(1) =  nz*dq(3) - ny*dq(4) + nx*wp1
    f(2) = -nz*dq(2) + nx*dq(4) + ny*wp1
    f(3) =  ny*dq(2) - nx*dq(3) + nz*wp1
    f(4) =  dqn + wp2
    f(5) = -dqn + wp2

end subroutine LxQ

subroutine RxQ(prim,gykb,dq,f,gama)
    use global_variables,only:nl,sml_sss
    implicit none
    real :: prim(nl),gykb(4),dq(nl),f(nl),gama
    real :: nx,ny,nz,nt,oon
    real :: rm,um,vm,wm,pm,c2,cm
    real :: dqn,wp1,wp2
    
    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)
    
    oon = 1.0/max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
    
    nx = nx*oon
    ny = ny*oon
    nz = nz*oon
    nt = nt*oon

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    c2 = gama*pm/rm
    if ( c2 < 0.0 ) then
       write(*,*)'RxQ(prim,gykb,dq,f)'
       write(*,*)'c2<0,gama,pm,rm',gama,pm,rm
    endif

    cm = sqrt(c2)
    
    dqn = nx*dq(1) + ny*dq(2) + nz*dq(3)
    wp1 = 0.5*rm*(dq(4) + dq(5))
    wp2 = 0.5*(dq(4) - dq(5))

    f(1) =  dqn                 + wp1/cm
    f(2) = -nz*dq(2) + ny*dq(3) + nx*wp2
    f(3) =  nz*dq(1) - nx*dq(3) + ny*wp2
    f(4) = -ny*dq(1) + nx*dq(2) + nz*wp2
    f(5) =                        wp1*cm

end subroutine RxQ

subroutine VALUE_HALF_NODE_45(n,nmax,ni,q,q_half) !lhy
    use global_variables,only: nijk2nd
    use define_precision_mod
    implicit none

    real(prec) :: A1,B1,A2,B2,C2,D2

    integer :: nmax,n,ni,m,i
    real    :: q(n,-3:nmax+4),q_half(n,-2:nmax+2)

    A1=9.0_prec
    B1=-1.0_prec
    A2=5.0_prec
    B2=15.0_prec
    C2=-5.0_prec
    D2=1.0_prec

!
! deal with the symmetric boundary condiction
!
    if( ni <= nijk2nd )then !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��

        do i=1,ni-1
        do m=1,n
           q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
        enddo
        enddo

        do m=1,n
!           q_half(m, 0) = 1.5*q(m, 1)-0.5*q(m,min(2,ni))  
!           q_half(m,ni) = 1.5*q(m,ni)-0.5*q(m,max(ni-1,1))
           q_half(m, 0) = ( 15._prec*q(m, 1)-10._prec*q(m,min(2,ni))  + 3._prec*q(m,min(3,ni))   )/8._prec
           q_half(m,ni) = ( 15._prec*q(m,ni)-10._prec*q(m,max(ni-1,1))+ 3._prec*q(m,max(ni-2,1)) )/8._prec
           
           !! lhy
           q_half(m,-1) = 0.5*(q(m,-1) + q(m, 0))
           q_half(m,-2) = 0.5*(q(m,-2) + q(m,-1))
           
           q_half(m,ni+1) = 0.5*(q(m,ni+1) + q(m,ni+2))
           q_half(m,ni+2) = 0.5*(q(m,ni+2) + q(m,ni+3))
        enddo

    else

        do i=2,ni-2  !
        do m=1,n
           q_half(m,i) = (A1*(q(m,i) + q(m,i+1)) + B1*(q(m,i+2) + q(m,i-1)))/16.0
        enddo
        enddo

        do m=1,n
           q_half(m,1   ) = (A2*q(m,1 )+B2*q(m,2   )+C2*q(m,3   )+D2*q(m,4   ))/16.0
           q_half(m,ni-1) = (A2*q(m,ni)+B2*q(m,ni-1)+C2*q(m,ni-2)+D2*q(m,ni-3))/16.0

           q_half(m,0   ) = (35.d0*q(m,1 )-35.d0*q(m,2   )+21.d0*q(m,3   )-5.d0*q(m,4   ))/16.0 !cic
           q_half(m,ni  ) = (35.d0*q(m,ni)-35.d0*q(m,ni-1)+21.d0*q(m,ni-2)-5.d0*q(m,ni-3))/16.0 !cic

!           q_half(m,0   ) = (15.d0*q(m,1 )-10.d0*q(m,2   )+3.d0*q(m,3   ))/8.0 !*TGH. NEW WCNS_E_5_BORDER
!           q_half(m,ni  ) = (15.d0*q(m,ni)-10.d0*q(m,ni-1)+3.d0*q(m,ni-2))/8.0 !*TGH. NEW WCNS_E_5_BORDER

           !! lhy
           q_half(m,-1) = (A2* q(m, 0) + B2*q(m,-1)  + C2* q(m,-2) + D2*q(m,-3) )/16.0
           q_half(m,-2) = (A1*(q(m,-2) +    q(m,-1)) + B1*(q(m, 0) +    q(m,-3)))/16.0
           
           q_half(m,ni+1) = (A2* q(m,ni+1) + B2*q(m,ni+2)  + C2* q(m,ni+3) + D2*q(m,ni+4) )/16.0
           q_half(m,ni+2) = (A1*(q(m,ni+2) +    q(m,ni+3)) + B1*(q(m,ni+4) +    q(m,ni+1)))/16.0
        enddo

    endif

    return
end subroutine VALUE_HALF_NODE_45
!=============================================================================!
