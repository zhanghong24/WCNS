!-----------------------------------------------------------------------------!
subroutine differential_model_turb_SA_GCL
	use global_variables,only:method
	implicit none

	call allocate_variables_temp_turbulence

	call set_dq_to_zero_turbulence

! 以下2阶
	call cal_rhs_source_turbulence
	call cal_rhs_viscous_turb_SA_GCL
	call cal_rhs_inviscous_turb_SA_GCL

! 以下4阶或3阶
!	call cal_rhs_source_turbulence_4th
!	call cal_rhs_viscous_turb_SA_GCL_4th
!	call cal_rhs_inviscous_turb_SA_GCL_3th
!	call cal_rhs_inviscous_turb_SA_GCL_4th


	call cal_spec_turbulence

	call lusgs_turbulence

	call update_turbulence

!	call residual_turbulence

	call deallocate_variables_temp_turbulence

	return
end subroutine differential_model_turb_SA_GCL

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_rhs_viscous_turb_SA_GCL
!* 使用几何守恒算法求SA模型粘性量。
!* 程序要用到1~max的viseq，所以viseq必须为已知的
  use define_precision_mod
  use global_variables,only : ni,nj,nk,nmax,method,nlamtur &
               , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz       &
			   , viseq,dmudxyz,dqke,re
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real(prec) :: ev_l,vis,nx,ny,nz
  real(prec) :: dmudx,dmudz,dmudy
  real(prec) :: fv( nlamtur,0:nmax ),dfv( nlamtur,nmax )
  real(prec) :: nxyz(3,nmax),nxyz0(3,0:nmax)

	cmthd = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1

			do i=1,ni
				nxyz(1,i) = kcx(i,j,k)			
				nxyz(2,i) = kcy(i,j,k)			
				nxyz(3,i) = kcz(i,j,k)			
			enddo
			
			call VALUE_HALF_NODE(3,nmax,ni,nxyz,nxyz0)	
					
			do i=cmthd,ni

				im = i
				in = i-1

				nx = nxyz0(1,in)
				ny = nxyz0(2,in)
				nz = nxyz0(3,in)

				do m =1,nlamtur

					dmudx =  dmudxyz(im,j,k,1,m)+dmudxyz(in,j,k,1,m) 
					dmudy =  dmudxyz(im,j,k,2,m)+dmudxyz(in,j,k,2,m)
					dmudz =  dmudxyz(im,j,k,3,m)+dmudxyz(in,j,k,3,m)

					vis  = viseq(im,j,k,m) + viseq(in,j,k,m)

					fv(m,in) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			i=0
			nx = nxyz0(1,i)
			ny = nxyz0(2,i)
			nz = nxyz0(3,i)
			do m =1,nlamtur
			    dmudx = 1.5*dmudxyz(1,j,k,1,m)-0.5*dmudxyz(2,j,k,1,m) 
			    dmudy = 1.5*dmudxyz(1,j,k,2,m)-0.5*dmudxyz(2,j,k,2,m) 
			    dmudz = 1.5*dmudxyz(1,j,k,3,m)-0.5*dmudxyz(2,j,k,3,m) 
				vis   = 1.5*viseq(1,j,k,m)   - 0.5*viseq(2,j,k,m)
				fv(m,i) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis
			enddo

			i=ni
			nx = nxyz0(1,i)
			ny = nxyz0(2,i)
			nz = nxyz0(3,i)
			do m =1,nlamtur
			    dmudx = 1.5*dmudxyz(i,j,k,1,m)-0.5*dmudxyz(i-1,j,k,1,m) 
			    dmudy = 1.5*dmudxyz(i,j,k,2,m)-0.5*dmudxyz(i-1,j,k,2,m) 
			    dmudz = 1.5*dmudxyz(i,j,k,3,m)-0.5*dmudxyz(i-1,j,k,3,m) 
				vis   = 1.5*viseq(i,j,k,m)   - 0.5*viseq(i-1,j,k,m)
				fv(m,i) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis
			enddo
			
			call FLUX_DXYZ(nlamtur,nmax,ni,fv,dfv)

			do m =1,nlamtur
				do i=cmthd,ni1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*dfv(m,i)
				enddo
			enddo
    enddo
    enddo

!------------
! J direction

	do k=cmthd,nk1
	do i=cmthd,ni1

			do j=1,nj
				nxyz(1,j) = etx(i,j,k)			
				nxyz(2,j) = ety(i,j,k)			
				nxyz(3,j) = etz(i,j,k)			
			enddo
			
			call VALUE_HALF_NODE(3,nmax,nj,nxyz,nxyz0)	

			do j=cmthd,nj

				jm = j
				jn = j-1

				nx = nxyz0(1,jn)
				ny = nxyz0(2,jn)
				nz = nxyz0(3,jn)

				do m =1,nlamtur

					dmudx =  dmudxyz(i,jm,k,1,m)+dmudxyz(i,jn,k,1,m) 
					dmudy =  dmudxyz(i,jm,k,2,m)+dmudxyz(i,jn,k,2,m)
					dmudz =  dmudxyz(i,jm,k,3,m)+dmudxyz(i,jn,k,3,m)

					vis  = viseq(i,jm,k,m) + viseq(i,jn,k,m)

					fv(m,jn) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			j=0
			nx = nxyz0(1,j)
			ny = nxyz0(2,j)
			nz = nxyz0(3,j)
			do m =1,nlamtur

				dmudx =  1.5*dmudxyz(i,1,k,1,m)-0.5*dmudxyz(i,2,k,1,m) 
				dmudy =  1.5*dmudxyz(i,1,k,2,m)-0.5*dmudxyz(i,2,k,2,m)
				dmudz =  1.5*dmudxyz(i,1,k,3,m)-0.5*dmudxyz(i,2,k,3,m)
				vis   =  1.5*viseq(i,1,k,m)    -0.5*viseq(i,2,k,m)
				fv(m,j) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis
			enddo

			j=nj
			nx = nxyz0(1,j)
			ny = nxyz0(2,j)
			nz = nxyz0(3,j)
			do m =1,nlamtur
				dmudx =  1.5*dmudxyz(i,j,k,1,m)-0.5*dmudxyz(i,j-1,k,1,m) 
				dmudy =  1.5*dmudxyz(i,j,k,2,m)-0.5*dmudxyz(i,j-1,k,2,m)
				dmudz =  1.5*dmudxyz(i,j,k,3,m)-0.5*dmudxyz(i,j-1,k,3,m)
				vis   =  1.5*viseq(i,j,k,m)    -0.5*viseq(i,j-1,k,m)
				fv(m,j) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis
			enddo

			call FLUX_DXYZ(nlamtur,nmax,nj,fv,dfv)

			do m =1,nlamtur
				do j=cmthd,nj1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*dfv(m,j)
				enddo
			enddo

    enddo
    enddo

!------------
! K direction
    do j=cmthd,nj1
	do i=cmthd,ni1

			do k=1,nk
				nxyz(1,k) = ctx(i,j,k)			
				nxyz(2,k) = cty(i,j,k)			
				nxyz(3,k) = ctz(i,j,k)			
			enddo
			
			call VALUE_HALF_NODE(3,nmax,nk,nxyz,nxyz0)
				
			do k=cmthd,nk

				km = k
				kn = k-1
				nx = nxyz0(1,kn)
				ny = nxyz0(2,kn)
				nz = nxyz0(3,kn)

				do m =1,nlamtur

					dmudx =  dmudxyz(i,j,km,1,m)+dmudxyz(i,j,kn,1,m) 
					dmudy =  dmudxyz(i,j,km,2,m)+dmudxyz(i,j,kn,2,m)
					dmudz =  dmudxyz(i,j,km,3,m)+dmudxyz(i,j,kn,3,m)

					vis  = viseq(i,j,km,m) + viseq(i,j,kn,m)

					fv(m,kn) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			k=0
			nx = nxyz0(1,k)
			ny = nxyz0(2,k)
			nz = nxyz0(3,k)
			do m =1,nlamtur

				dmudx =  1.5*dmudxyz(i,j,1,1,m)-0.5*dmudxyz(i,j,2,1,m) 
				dmudy =  1.5*dmudxyz(i,j,1,2,m)-0.5*dmudxyz(i,j,2,2,m)
				dmudz =  1.5*dmudxyz(i,j,1,3,m)-0.5*dmudxyz(i,j,2,3,m)
				vis   =  1.5*viseq(i,j,1,m)    -0.5*viseq(i,j,2,m)
				fv(m,k) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis
			enddo

			k=nk
			nx = nxyz0(1,k)
			ny = nxyz0(2,k)
			nz = nxyz0(3,k)
			do m =1,nlamtur
				dmudx =  1.5*dmudxyz(i,j,k,1,m)-0.5*dmudxyz(i,j,k-1,1,m) 
				dmudy =  1.5*dmudxyz(i,j,k,2,m)-0.5*dmudxyz(i,j,k-1,2,m)
				dmudz =  1.5*dmudxyz(i,j,k,3,m)-0.5*dmudxyz(i,j,k-1,3,m)
				vis   =  1.5*viseq(i,j,k,m)    -0.5*viseq(i,j,k-1,m)
				fv(m,k) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis
			enddo

			call FLUX_DXYZ(nlamtur,nmax,nk,fv,dfv)

			do m =1,nlamtur
				do k=cmthd,nk1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*dfv(m,k)
				enddo
			enddo

    enddo
    enddo

!------------
  return
end subroutine cal_rhs_viscous_turb_SA_GCL

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine cal_rhs_inviscous_turb_SA_GCL
  use define_precision_mod
  use global_variables,only : nl,ni,nj,nk,nmax,method,nlamtur &
               , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz       &
			   , r,u,v,w,qke,dqke,re,xk,xb,efix

  implicit none

  integer :: i,j,k,m,cmethd
  real(prec) :: fc(1:nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)

  cmethd = 1 + method
	do m=1,nlamtur
    do k= cmethd,nk-1    
      do j= cmethd,nj-1      
        do i= method-1,ni+1    
          q_line(1,i) = r(i,j,k)
          q_line(2,i) = u(i,j,k)
          q_line(3,i) = v(i,j,k)
          q_line(4,i) = w(i,j,k)
          q_line(5,i) = qke(i,j,k,m)
        enddo
        do i=1,ni
           trxyz(1,i) = kcx(i,j,k)
           trxyz(2,i) = kcy(i,j,k)
           trxyz(3,i) = kcz(i,j,k)
        enddo

        call NND_turb_GCL(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

        do i=cmethd,ni-1     
          dqke(i,j,k,m) = dqke(i,j,k,m) +  fc(i)
        enddo

       enddo
    enddo
    do k=cmethd,nk-1        
      do i=cmethd,ni-1     
        do j=method-1,nj+1  
          q_line(1,j) = r(i,j,k)
          q_line(2,j) = u(i,j,k)
          q_line(3,j) = v(i,j,k)
          q_line(4,j) = w(i,j,k)
          q_line(5,j) = qke(i,j,k,m)
        enddo

        do j=1,nj
          trxyz(1,j) = etx(i,j,k)
          trxyz(2,j) = ety(i,j,k)
          trxyz(3,j) = etz(i,j,k)
        enddo

        call NND_turb_GCL(xk,xb,nj,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

        do j=cmethd,nj-1     
          dqke(i,j,k,m) = dqke(i,j,k,m) +  fc(j)
        enddo
      enddo
    enddo

    do j=cmethd,nj-1       
       do i=cmethd,ni-1    
          do k=method-1,nk+1    
            q_line(1,k) = r(i,j,k)
            q_line(2,k) = u(i,j,k)
            q_line(3,k) = v(i,j,k)
            q_line(4,k) = w(i,j,k)
            q_line(5,k) = qke(i,j,k,m)
          enddo
          do k=1,nk
            trxyz(1,k) = ctx(i,j,k)
            trxyz(2,k) = cty(i,j,k)
            trxyz(3,k) = ctz(i,j,k)
          enddo
          call NND_turb_GCl(xk,xb,nk,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)
          do k=cmethd,nk-1      
            dqke(i,j,k,m) = dqke(i,j,k,m) + fc(k)
          enddo
       enddo
    enddo
	enddo
  return
end subroutine cal_rhs_inviscous_turb_SA_GCL

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine NND_turb_GCL(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,dfc)
!* 需要用到虚点上的值
  use define_precision_mod
  implicit none

    integer :: i,ii,m,cmethd,ni,nmax,nl,method
    real(prec),dimension( 1:nl ) :: dql,dqr
    real(prec) :: nx,ny,nz,xk,xb,efix,gamaeq,fl,fr
    real(prec) :: c1,c2,dfc(1:nmax)
    real(prec) :: fc(0:nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)
    real(prec) :: d_q(1:nl,0:nmax+1),d_qcl(1:nl,0:nmax),d_qcr(1:nl,0:nmax)
    real(prec),external :: minmod
	real(prec) :: nxyz0(3,0:nmax)

    c1 = 0.25 * ( 1.0 + xk )
    c2 = 0.25 * ( 1.0 - xk )

    do i=method,ni+1 
       do m=1,nl
           d_q(m,i) = q_line(m,i) - q_line(m,i-1)
       enddo
    enddo

!    do i=1,method
!       do m=1,nl
!          d_q(m,cmethd-1) = d_q(m,cmethd)
!          d_q(m,ni+1  )   = d_q(m,ni)
!       enddo
!    enddo

    do i=method,ni     
       do m=1,nl
          d_qcl(m,i) = minmod( xb*d_q(m,i+1), d_q(m,i) )
          d_qcr(m,i) = minmod( d_q(m,i+1), xb*d_q(m,i) )
       enddo
    enddo

!	do i=1,ni
!	do m=1,3
!       nxyz(m,i) = trxyz(m,i)
!	enddo
!	enddo

	call VALUE_HALF_NODE(3,nmax,ni,trxyz,nxyz0)

    do i=cmethd,ni ! here, cmethod = 2
	
	   ii=i-1     
       nx = nxyz0(1,ii)
       ny = nxyz0(2,ii)
       nz = nxyz0(3,ii)

       do m=1,nl
          dql(m) = q_line(m,ii) + c1 * d_qcr(m,ii) + c2 * d_qcl(m,ii)
          dqr(m) = q_line(m,i ) - c1 * d_qcl(m,i ) - c2 * d_qcr(m,i )
       enddo

       call flux_turbulence(dql,nx,ny,nz,fl, 1)
       call flux_turbulence(dqr,nx,ny,nz,fr,-1)

       fc(ii) = fl + fr
    enddo

	 ii=0    
     nx = nxyz0(1,ii)
     ny = nxyz0(2,ii)
     nz = nxyz0(3,ii)
     do m=1,nl
        dql(m) = 0.5*(q_line(m,1) + q_line(m,0) ) !* 要用虚点值
        dqr(m) = dql(m)
     enddo
     call flux_turbulence(dql,nx,ny,nz,fl, 1)
     call flux_turbulence(dqr,nx,ny,nz,fr,-1)
     fc(ii) = fl + fr

	 ii=ni    
     nx = nxyz0(1,ii)
     ny = nxyz0(2,ii)
     nz = nxyz0(3,ii)
     do m=1,nl
        dql(m) = 0.5*(q_line(m,ni) + q_line(m,ni+1) ) !* 要用虚点值
        dqr(m) = dql(m)
     enddo
     call flux_turbulence(dql,nx,ny,nz,fl, 1)
     call flux_turbulence(dqr,nx,ny,nz,fr,-1)
     fc(ii) = fl + fr

	 call FLUX_DXYZ(1,nmax,ni,fc,dfc)

    return
end subroutine NND_turb_GCL

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_rhs_viscous_turb_SA_GCL_4th
!* 使用几何守恒算法求SA模型粘性量。
!* 程序要用到1~max的viseq，所以viseq必须为已知的
  use define_precision_mod
  use global_variables,only : ni,nj,nk,nmax,method,nlamtur &
               , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz       &
			   , viseq,dmudxyz,dqke,re
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real(prec) :: ev_l,vis,nx,ny,nz
  real(prec) :: dmudx,dmudz,dmudy
  real(prec) :: fv( nlamtur,0:nmax ),dfv( nlamtur,nmax )
  real(prec) :: nxyz(3,nmax),nxyz0(3,0:nmax),tem4(4,nmax),tem4h(4,0:nmax)

	cmthd = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1

			do i=1,ni
				nxyz(1,i) = kcx(i,j,k)			
				nxyz(2,i) = kcy(i,j,k)			
				nxyz(3,i) = kcz(i,j,k)			
			enddo
			
			call VALUE_HALF_NODE(3,nmax,ni,nxyz,nxyz0)	

			do i=1,ni
			do m =1,nlamtur
				tem4(1,i) = dmudxyz(i,j,k,1,m)
				tem4(2,i) = dmudxyz(i,j,k,2,m)
				tem4(3,i) = dmudxyz(i,j,k,3,m)
				tem4(4,i) =   viseq(i,j,k,m)
			enddo
			enddo

			call VALUE_HALF_NODE(4,nmax,ni,tem4,tem4h)	
					
			do i=0,ni
				nx = nxyz0(1,i)
				ny = nxyz0(2,i)
				nz = nxyz0(3,i)

				do m =1,nlamtur
					dmudx =  tem4h(1,i)
					dmudy =  tem4h(2,i)
					dmudz =  tem4h(3,i)
					vis   =  tem4h(4,i)
					fv(m,i) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis
				enddo

			enddo

			call FLUX_DXYZ(nlamtur,nmax,ni,fv,dfv)

			do m =1,nlamtur
				do i=cmthd,ni1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*dfv(m,i)
				enddo
			enddo
    enddo
    enddo

!------------
! J direction

	do k=cmthd,nk1
	do i=cmthd,ni1

			do j=1,nj
				nxyz(1,j) = etx(i,j,k)			
				nxyz(2,j) = ety(i,j,k)			
				nxyz(3,j) = etz(i,j,k)			
			enddo
			
			call VALUE_HALF_NODE(3,nmax,nj,nxyz,nxyz0)	

			do j=1,nj
			do m =1,nlamtur
				tem4(1,j) = dmudxyz(i,j,k,1,m)
				tem4(2,j) = dmudxyz(i,j,k,2,m)
				tem4(3,j) = dmudxyz(i,j,k,3,m)
				tem4(4,j) =   viseq(i,j,k,m)
			enddo
			enddo

			call VALUE_HALF_NODE(4,nmax,nj,tem4,tem4h)

			do j=1,nj

				nx = nxyz0(1,j)
				ny = nxyz0(2,j)
				nz = nxyz0(3,j)

				do m =1,nlamtur

					dmudx =  tem4h(1,j)
					dmudy =  tem4h(2,j)
					dmudz =  tem4h(3,j)
					vis   =  tem4h(4,j)
					fv(m,j) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis

				enddo
			enddo

			call FLUX_DXYZ(nlamtur,nmax,nj,fv,dfv)

			do m =1,nlamtur
				do j=cmthd,nj1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*dfv(m,j)
				enddo
			enddo

    enddo
    enddo

!------------
! K direction
    do j=cmthd,nj1
	do i=cmthd,ni1

			do k=1,nk
				nxyz(1,k) = ctx(i,j,k)			
				nxyz(2,k) = cty(i,j,k)			
				nxyz(3,k) = ctz(i,j,k)			
			enddo
			
			call VALUE_HALF_NODE(3,nmax,nk,nxyz,nxyz0)
				
			do k=1,nk
			do m =1,nlamtur
				tem4(1,k) = dmudxyz(i,j,k,1,m)
				tem4(2,k) = dmudxyz(i,j,k,2,m)
				tem4(3,k) = dmudxyz(i,j,k,3,m)
				tem4(4,k) =   viseq(i,j,k,m)
			enddo
			enddo

			call VALUE_HALF_NODE(4,nmax,nk,tem4,tem4h)

			do k=1,nk

				nx = nxyz0(1,k)
				ny = nxyz0(2,k)
				nz = nxyz0(3,k)

				do m =1,nlamtur

					dmudx =  tem4h(1,k)
					dmudy =  tem4h(2,k)
					dmudz =  tem4h(3,k)
					vis   =  tem4h(4,k)
					fv(m,k) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis

				enddo
			enddo

			call FLUX_DXYZ(nlamtur,nmax,nk,fv,dfv)

			do m =1,nlamtur
				do k=cmthd,nk1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*dfv(m,k)
				enddo
			enddo

    enddo
    enddo

!------------
  return
end subroutine cal_rhs_viscous_turb_SA_GCL_4th

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine cal_rhs_inviscous_turb_SA_GCL_3th
  use define_precision_mod
  use global_variables,only : nl,ni,nj,nk,nmax,method,nlamtur &
               , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz       &
			   , r,u,v,w,qke,dqke,re,xk,xb,efix

  implicit none

  integer :: i,j,k,m,cmethd
  real(prec) :: fc(nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)

  cmethd = 1 + method
	do m=1,nlamtur
    do k= cmethd,nk-1    
      do j= cmethd,nj-1      
        do i= 1,ni    
          q_line(1,i) = r(i,j,k)
          q_line(2,i) = u(i,j,k)
          q_line(3,i) = v(i,j,k)
          q_line(4,i) = w(i,j,k)
          q_line(5,i) = qke(i,j,k,m)
        enddo
        do i=1,ni
           trxyz(1,i) = kcx(i,j,k)
           trxyz(2,i) = kcy(i,j,k)
           trxyz(3,i) = kcz(i,j,k)
        enddo

        call turb_GCL_3th(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

        do i=cmethd,ni-1     
          dqke(i,j,k,m) = dqke(i,j,k,m) +  fc(i)
        enddo

       enddo
    enddo
    do k=cmethd,nk-1        
      do i=cmethd,ni-1     
        do j=1,nj  
          q_line(1,j) = r(i,j,k)
          q_line(2,j) = u(i,j,k)
          q_line(3,j) = v(i,j,k)
          q_line(4,j) = w(i,j,k)
          q_line(5,j) = qke(i,j,k,m)
        enddo

        do j=1,nj
          trxyz(1,j) = etx(i,j,k)
          trxyz(2,j) = ety(i,j,k)
          trxyz(3,j) = etz(i,j,k)
        enddo

        call turb_GCL_3th(xk,xb,nj,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

        do j=cmethd,nj-1     
          dqke(i,j,k,m) = dqke(i,j,k,m) +  fc(j)
        enddo
      enddo
    enddo

    do j=cmethd,nj-1       
       do i=cmethd,ni-1    
          do k=1,nk    
            q_line(1,k) = r(i,j,k)
            q_line(2,k) = u(i,j,k)
            q_line(3,k) = v(i,j,k)
            q_line(4,k) = w(i,j,k)
            q_line(5,k) = qke(i,j,k,m)
          enddo
          do k=1,nk
            trxyz(1,k) = ctx(i,j,k)
            trxyz(2,k) = cty(i,j,k)
            trxyz(3,k) = ctz(i,j,k)
          enddo
          call turb_GCL_3th(xk,xb,nk,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)
          do k=cmethd,nk-1      
            dqke(i,j,k,m) = dqke(i,j,k,m) + fc(k)
          enddo
       enddo
    enddo
	enddo
  return
end subroutine cal_rhs_inviscous_turb_SA_GCL_3th

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine turb_GCL_3th(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,dfc)
!* 需要用到虚点上的值
  use define_precision_mod
  use global_variables,only: nijk2nd
  implicit none

    integer :: i,ii,m,cmethd,ni,nmax,nl,method
    real(prec),dimension( 1:nl ) :: dql,dqr
    real(prec) :: nx,ny,nz,xk,xb,efix,gamaeq,fl,fr
    real(prec) :: c1,c2,dfc(1:nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)
    real(prec) :: fc(0:nmax)
!    real(prec) :: d_q(1:nl,0:nmax+1),d_qcl(1:nl,0:nmax),d_qcr(1:nl,0:nmax)
!    real(prec),external :: minmod
	real(prec) :: nxyz0(3,0:nmax)

	call VALUE_HALF_NODE(3,nmax,ni,trxyz,nxyz0)

	if (ni <= nijk2nd)then  !* 特征注意，可能是2维，将为2阶
        do i=0,ni

           nx = nxyz0(1,i)
           ny = nxyz0(2,i)
           nz = nxyz0(3,i)
		   if(i==0)then
		     do m=1,nl
		      dql(m) = ( 15.*q_line(m, 1)-10.*q_line(m,min(2,ni))  + 3.*q_line(m,min(3,ni))   )/8.
			 enddo

		   elseif(i==ni)then
			 do m=1,nl
			  dql(m) = ( 15.*q_line(m,ni)-10.*q_line(m,max(ni-1,1))+ 3.*q_line(m,max(ni-2,1)) )/8.
			 enddo

		   else

		     do m=1,nl
				dql(m) = 0.5*(q_line(m,i) + q_line(m,i+1))
		     enddo
		   endif

		   do m=1,nl
				dqr(m) = dql(m)
		   enddo

           call flux_turbulence(dql,nx,ny,nz,fl, 1)
		   call flux_turbulence(dqr,nx,ny,nz,fr,-1)

		   fc(i) = fl + fr

		enddo

	    call FLUX_DXYZ(1,nmax,ni,fc,dfc)

		return

	endif
    do i=0,ni

       nx = nxyz0(1,i)
       ny = nxyz0(2,i)
       nz = nxyz0(3,i)

	    do m=1,nl

	        if(i.eq.0 )then

			dql(m) = (35.0*q_line(m,1 )-35.0*q_line(m,2   )+21.0*q_line(m,3   )-5.0*q_line(m,4   ))/16.0 
			dqr(m) = dql(m)

		    elseif(i.eq.ni)then

			dql(m) = (35.0*q_line(m,ni)-35.0*q_line(m,ni-1)+21.0*q_line(m,ni-2)-5.0*q_line(m,ni-3))/16.0 
			dqr(m) = dql(m)

		    elseif(i.eq.1)then

			dql(m) = (5.*q_line(m,1 )+15.*q_line(m,2   )-5.*q_line(m,3   )+q_line(m,4   ))/16.0
			dqr(m) = (3.*q_line(m,i)   + 6.*q_line(m,i+1) - q_line(m,i+2) )/8.0

		    elseif(i.eq.ni-1)then

 			dql(m) = (3.*q_line(m,i+1) + 6.*q_line(m,i  ) - q_line(m,i-1) )/8.0 
			dqr(m) = (5.*q_line(m,ni)+15.*q_line(m,ni-1)-5.*q_line(m,ni-2)+q_line(m,ni-3))/16.0

		    else

 			dql(m) = (3.*q_line(m,i+1) + 6.*q_line(m,i  ) - q_line(m,i-1) )/8.0 
			dqr(m) = (3.*q_line(m,i)   + 6.*q_line(m,i+1) - q_line(m,i+2) )/8.0

		    endif

        enddo

        call flux_turbulence(dql,nx,ny,nz,fl, 1)
		call flux_turbulence(dqr,nx,ny,nz,fr,-1)

		fc(i) = fl + fr
	enddo

	call FLUX_DXYZ(1,nmax,ni,fc,dfc)

    return
end subroutine turb_GCL_3th

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine cal_rhs_inviscous_turb_SA_GCL_4th
  use define_precision_mod
  use global_variables,only : nl,ni,nj,nk,nmax,method,nlamtur &
               , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz       &
			   , r,u,v,w,qke,dqke,re,xk,xb,efix

  implicit none

  integer :: i,j,k,m,cmethd
  real(prec) :: fc(nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)

  cmethd = 1 + method
	do m=1,nlamtur
      do k= cmethd,nk-1    
      do j= cmethd,nj-1      
        do i= 1,ni    
          q_line(1,i) = r(i,j,k)
          q_line(2,i) = u(i,j,k)
          q_line(3,i) = v(i,j,k)
          q_line(4,i) = w(i,j,k)
          q_line(5,i) = qke(i,j,k,m)
        enddo
        do i=1,ni
           trxyz(1,i) = kcx(i,j,k)
           trxyz(2,i) = kcy(i,j,k)
           trxyz(3,i) = kcz(i,j,k)
        enddo

        call turb_GCL_4th(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

        do i=cmethd,ni-1     
          dqke(i,j,k,m) = dqke(i,j,k,m) +  fc(i)
        enddo

      enddo
      enddo

      do k=cmethd,nk-1        
      do i=cmethd,ni-1     
        do j=1,nj  
          q_line(1,j) = r(i,j,k)
          q_line(2,j) = u(i,j,k)
          q_line(3,j) = v(i,j,k)
          q_line(4,j) = w(i,j,k)
          q_line(5,j) = qke(i,j,k,m)
        enddo

        do j=1,nj
          trxyz(1,j) = etx(i,j,k)
          trxyz(2,j) = ety(i,j,k)
          trxyz(3,j) = etz(i,j,k)
        enddo

        call turb_GCL_4th(xk,xb,nj,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

        do j=cmethd,nj-1     
          dqke(i,j,k,m) = dqke(i,j,k,m) +  fc(j)
        enddo

      enddo
      enddo

      do j=cmethd,nj-1       
      do i=cmethd,ni-1    
          do k=1,nk    
            q_line(1,k) = r(i,j,k)
            q_line(2,k) = u(i,j,k)
            q_line(3,k) = v(i,j,k)
            q_line(4,k) = w(i,j,k)
            q_line(5,k) = qke(i,j,k,m)
          enddo
          do k=1,nk
            trxyz(1,k) = ctx(i,j,k)
            trxyz(2,k) = cty(i,j,k)
            trxyz(3,k) = ctz(i,j,k)
          enddo
          call turb_GCL_4th(xk,xb,nk,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)
          do k=cmethd,nk-1      
            dqke(i,j,k,m) = dqke(i,j,k,m) + fc(k)
          enddo
      enddo
      enddo
	enddo

    return
end subroutine cal_rhs_inviscous_turb_SA_GCL_4th

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine turb_GCL_4th(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,dfc)
!* 需要用到虚点上的值
  use define_precision_mod
  use global_variables,only : nijk2nd
  implicit none

    integer :: i,ii,m,cmethd,ni,nmax,nl,method
    real(prec),dimension( 1:nl ) :: dql,dqr
    real(prec) :: nx,ny,nz,xk,xb,efix,gamaeq,fl,fr
    real(prec) :: c1,c2,dfc(1:nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)
    real(prec) :: fc(0:nmax)
!    real(prec) :: d_q(1:nl,0:nmax+1),d_qcl(1:nl,0:nmax),d_qcr(1:nl,0:nmax)
!    real(prec),external :: minmod
	real(prec) :: nxyz0(3,0:nmax)

	call VALUE_HALF_NODE(3,nmax,ni,trxyz,nxyz0)

	if (ni <= nijk2nd)then  !* 特征注意，可能是2维，将为2阶
        do i=0,ni

           nx = nxyz0(1,i)
           ny = nxyz0(2,i)
           nz = nxyz0(3,i)
		   if(i==0)then
		     do m=1,nl
		      dql(m) = ( 15.*q_line(m, 1)-10.*q_line(m,min(2,ni))  + 3.*q_line(m,min(3,ni))   )/8.
			 enddo

		   elseif(i==ni)then
			 do m=1,nl
			  dql(m) = ( 15.*q_line(m,ni)-10.*q_line(m,max(ni-1,1))+ 3.*q_line(m,max(ni-2,1)) )/8.
			 enddo

		   else

		     do m=1,nl
				dql(m) = 0.5*(q_line(m,i) + q_line(m,i+1))
		     enddo
		   endif

		   do m=1,nl
				dqr(m) = dql(m)
		   enddo

           call flux_turbulence(dql,nx,ny,nz,fl, 1)
		   call flux_turbulence(dqr,nx,ny,nz,fr,-1)

		   fc(i) = fl + fr

		enddo

	    call FLUX_DXYZ(1,nmax,ni,fc,dfc)

		return

	endif

    do i=0,ni

       nx = nxyz0(1,i)
       ny = nxyz0(2,i)
       nz = nxyz0(3,i)


	   if(i==0 )then

	      do m=1,nl
			dql(m) = (35.0*q_line(m,1 )-35.0*q_line(m,2   )+21.0*q_line(m,3   )-5.0*q_line(m,4   ))/16.0 
			dqr(m) = dql(m)
		  enddo

		elseif(i==ni)then

	      do m=1,nl
			dql(m) = (35.0*q_line(m,ni)-35.0*q_line(m,ni-1)+21.0*q_line(m,ni-2)-5.0*q_line(m,ni-3))/16.0 
			dqr(m) = dql(m)
		  enddo

	    elseif(i==1)then

	      do m=1,nl
			dql(m) = (5.*q_line(m,1 )+15.*q_line(m,2   )-5.*q_line(m,3   )+q_line(m,4   ))/16.0
			dqr(m) = dql(m)
		  enddo

	    elseif(i==ni-1)then

	      do m=1,nl
			dql(m) = (5.*q_line(m,ni)+15.*q_line(m,ni-1)-5.*q_line(m,ni-2)+q_line(m,ni-3))/16.0
			dqr(m) = dql(m)
		  enddo

	    elseif(i == 2)then

	      do m=1,nl
 			dql(m) = (9.*(q_line(m,i)+q_line(m,i+1)) - (q_line(m,i+2)+q_line(m,i-1)))/16.0
			dqr(m) = (5.*q_line(m,i  )+15.*q_line(m,i+1)-5.*q_line(m,i+2)+q_line(m,i+3) )/16.0
		  enddo

		elseif(i==ni-2)then

	      do m=1,nl
			dql(m) = (5.*q_line(m,i+1)+15.*q_line(m,i  )-5.*q_line(m,i-1)+q_line(m,i-2) )/16.0
 			dqr(m) = (9.*(q_line(m,i)+q_line(m,i+1)) - (q_line(m,i+2)+q_line(m,i-1)))/16.0
		  enddo

		else

	      do m=1,nl
			dql(m) = (5.*q_line(m,i+1)+15.*q_line(m,i  )-5.*q_line(m,i-1)+q_line(m,i-2) )/16.0
			dqr(m) = (5.*q_line(m,i  )+15.*q_line(m,i+1)-5.*q_line(m,i+2)+q_line(m,i+3) )/16.0
		  enddo

		endif

        call flux_turbulence(dql,nx,ny,nz,fl, 1)
		call flux_turbulence(dqr,nx,ny,nz,fr,-1)

		fc(i) = fl + fr
	enddo

	call FLUX_DXYZ(1,nmax,ni,fc,dfc)

    return
end subroutine turb_GCL_4th

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_rhs_source_turbulence_4th
	use global_variables
	implicit none

	if(nameturb=='SA')then
		call cal_rhs_source_SA_4th
	else
		write(*,*)trim(nameturb),'has not been realized in this release version'
		stop
	endif

	return
end subroutine cal_rhs_source_turbulence_4th

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine cal_rhs_source_SA_4th
	use global_variables
	implicit none
	real,parameter :: cb1=0.1355,cb2=0.622,cv1=7.1,eps=1.0e-20 ,    &
	                  ct1=1.0   ,ct2=2.0  ,ct3=0.3,ct4=2.0     ,    &
										cga=2./3.0,ki =0.419,cw2=0.3,cw3=2.0     ,    &
										cw1=cb1/ki/ki+(1.0+cb2)/cga,cv13=cv1**3.0,    &
										ki2=ki*ki,cw36=cw3**6.0,onesix=1.0/6.0   ,    &
										c5=3.5

	integer :: i,j,k,cmethod,n,ni1,nj1,nk1
	real    :: kafan,fv1,fv2,kf3,fw,ft1,ft2,gbar6,rbar5
	real    :: sbar,rbar,gbar
	real    :: dn,dn2,omg,reki2dn2
	real    :: Sp,Sd,Dbar,ro,vsl,vtbar,romu,volume
	real    :: dmudx,dmudy,dmudz,cdmu,cpmu,usound,ux,vy,wz,Sdc
	real    :: dfv2dkafan,dsbardmu,drbardmu,dgbardrbar,dfwdg,dfwdmu
	n       =   1
	cmethod =   method + 1
	ft1     =   0.0
	ft2     =   0.0
 
	call dfidxyz_center_4th(-2,3,ni,nj,nk,u  ,dudxyz ,method,1)
	call dfidxyz_center_4th(-2,3,ni,nj,nk,v  ,dvdxyz ,method,1)
	call dfidxyz_center_4th(-2,3,ni,nj,nk,w  ,dwdxyz ,method,1)
	call dfidxyz_center_4th(-2,3,ni,nj,nk,qke,dmudxyz,method,n)

	ni1 = ni-1
	nj1 = nj-1
	nk1 = nk-1

	do k = cmethod,nk1
		do j = cmethod,nj1
			do i = cmethod,ni1

				ro      =    r(i,j,k)
				volume  =    vol(i,j,k)
				vsl     =		 visl(i,j,k)
				vtbar   =    qke(i,j,k,n)
				romu    =    ro*vtbar

				dn      =    ds_turbulence(i,j,k)
				dn2     =    dn*dn

				reki2dn2=    reynolds*ki2*dn2
				
				kafan   =    romu/vsl  
				kf3     =    kafan**3.0
				fv1     =    kf3/(kf3+cv13)
				fv2     =    1.0-kafan/(1.0+kafan*fv1)

				omg     =    sqrt((dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                  (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
													(dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

				sbar    =    omg+vtbar/reki2dn2*fv2
				sbar    =    amax1(sbar,eps)           ! to avoid sbar=0.0

				rbar    =    vtbar/(sbar*reki2dn2)
				rbar    =    amin1(rbar,10.0)          ! rbar limiter

				rbar5   =    rbar**5.0
				gbar    =    rbar+cw2*rbar*(rbar5-1.0)
				gbar6   =    gbar**6.0

				fw      =    gbar*((1.0+cw36)/(gbar6+cw36))**onesix

				cpmu    =    cb1*(1.0-ft2)*sbar  !*vtbar
				cdmu    =   (cw1*fw-cb1*ft2/ki2)*vtbar/dn2/reynolds

				Sp      =    cpmu*vtbar
				Sd      =    cdmu*vtbar

				dmudx   =    dmudxyz(i,j,k,1,n)
				dmudy   =    dmudxyz(i,j,k,2,n)
				dmudz   =    dmudxyz(i,j,k,3,n)

				Dbar    =    cb2/cga*(dmudx**2.0+dmudy**2.0+dmudz**2.0)*re

				dfv2dkafan    = (3.0*cv13*(fv1/kafan)**2.0-1.0)/(1.0+kafan*fv1)**2.0

				dsbardmu      = (fv2+kafan*dfv2dkafan)/reki2dn2 

				drbardmu			= (rbar-rbar/sbar*dsbardmu)/ro

				dgbardrbar		=  1.0+cw2*(6.0*rbar5-1.0)

				dfwdg         =  fw/gbar*(1.0-gbar6/(gbar6+cw36))

				dfwdmu        =  dfwdg*dgbardrbar*drbardmu
				
				dqke(i,j,k,n) = dqke(i,j,k,n)-ro*(Sp-Sd+Dbar)*volume

				spec(i,j,k,n) = spec(i,j,k,n)-amin1(cpmu-2.0*cdmu,0.0)*volume
				if(ncmpcor==1)then 
					usound =  c(i,j,k)
					ux     =  dudxyz(i,j,k,1)
					vy     =  dvdxyz(i,j,k,2)
					wz     =  dvdxyz(i,j,k,3)
					Sdc    =  re*c5*ro*(vtbar*(ux+vy+wz)/usound)**2
					dqke(i,j,k,n) = dqke(i,j,k,n)-Sdc*volume
				endif
			enddo
		enddo
	enddo

	call viscosity_equivalent_SA

	return
end subroutine cal_rhs_source_SA_4th
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine dfidxyz_center_4th(nsub,nplus,ni,nj,nk,fi,dfidxyz,method,nv)
	use global_variables,only:kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol,nijk2nd
	implicit none

	integer :: nsub,nplus,ni,nj,nk,method,i,j,k,n,nv,m,ni1,nj1,nk1
	real    :: ddx,ddy,ddz,volp
	real    :: dkc,det,dct
	real    :: fim,fjm,fkm
	real    :: fi(nsub:ni+nplus,nsub:nj+nplus,nsub:nk+nplus,nv)
	real    :: dfidxyz(ni,nj,nk,3,nv)

	do m=1,nv
		do k=1,nk
		do j=1,nj
		do i=1,ni
			if(ni <= nijk2nd)then
				dkc = 0.5*(fi(i+1,j,k,m)-fi(i-1,j,k,m))
			elseif(i==1 .or. i==ni)then
				dkc = 0.5*(fi(i+1,j,k,m)-fi(i-1,j,k,m))
			else
				dkc = ( 8.*(fi(i+1,j,k,m)-fi(i-1,j,k,m)) - (fi(i+2,j,k,m)-fi(i-2,j,k,m)) )/12.
			endif

			if(nj <= nijk2nd)then
				det = 0.5*(fi(i,j+1,k,m)-fi(i,j-1,k,m))
			elseif(j==1 .or. j==nj)then
				det = 0.5*(fi(i,j+1,k,m)-fi(i,j-1,k,m))
			else
				det = ( 8.*(fi(i,j+1,k,m)-fi(i,j-1,k,m)) - (fi(i,j+2,k,m)-fi(i,j-2,k,m)) )/12.
			endif

			if(nk <= nijk2nd)then
				dct = 0.5*(fi(i,j,k+1,m)-fi(i,j,k-1,m))
			elseif(k==1 .or. k==nk)then
				dct = 0.5*(fi(i,j,k+1,m)-fi(i,j,k-1,m))
			else
				dct = ( 8.*(fi(i,j,k+1,m)-fi(i,j,k-1,m)) - (fi(i,j,k+2,m)-fi(i,j,k-2,m)) )/12.
			endif

			ddx = dkc*kcx(i,j,k)+det*etx(i,j,k)+dct*ctx(i,j,k)
			ddy = dkc*kcy(i,j,k)+det*ety(i,j,k)+dct*cty(i,j,k)
			ddz = dkc*kcz(i,j,k)+det*etz(i,j,k)+dct*ctz(i,j,k)

			volp =vol(i,j,k)

			dfidxyz(i,j,k,1,m) = ddx/volp
			dfidxyz(i,j,k,2,m) = ddy/volp
			dfidxyz(i,j,k,3,m) = ddz/volp
		enddo
		enddo
		enddo
	enddo

	return
end subroutine dfidxyz_center_4th
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

!=============================================================================!

subroutine WCNSE5_2nd_VIS_virtual_gcl
  use define_precision_mod
  use global_variables,only : u,v,w,t,reynolds,visl,vist,dq,nvis     &
                            , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol &
							, nl,nmax,ni,nj,nk,gama,moo,prl,prt
  use duvwt_all_field
  implicit none
!-----------------------------------------------------------------------------!
!	在原 WCNSE5_VIS 要用到虚一层的值                                          !
!	                                                                          !
!                                                                             !
!		设计：涂国华                                                          !
!		调试：涂国华 2009.03                                                  !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
	real(prec) :: duvwtdxyz (12,0:nmax),nxyz_line(3,0:nmax)
	real(prec) :: duvwt_line(12,nmax) !,uvwt_line(4,nmax)
	real(prec) :: duvwt_half(12,0:nmax),uvwt_half(4,0:nmax)
	real(prec) :: kxyz_line ( 9,nmax),kxyz_half(9,0:nmax)
	real(prec) :: vol_line(nmax),fv(nl,0:nmax),dfv(nl,nmax)
	real(prec) :: vol_half(0:nmax) !,txyz_half(9,0:nmax)
	real(prec) :: vslt1_half(0:nmax),vslt2_half(0:nmax),cp,cp_prl,cp_prt
	real(prec) :: vslt1_line(nmax),vslt2_line(nmax),re
	integer :: mvist
    real(prec) :: uvwt_lvir(4,0:nmax+1) !为了便于使用虚点上的值

	allocate(duvwt(12,ni,nj,nk))   ! overcome stack over problem
	allocate(duvwt_mid(12,0:ni,0:nj,0:nk))   !___ 直接计算半结点1阶导数  2009.2.1

	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)
    cp_prl = cp/prl
    cp_prt = cp/prt

! to calculate the deriative of u,v,w,t in computing	coordinate

	call UVWT_DER_2nd  !*TGH. 直接求节点1阶导数，暂时没有用虚点值
    call UVWT_DER_2nd_half_virtual !*TGH. 直接求半结点1阶导数,要用虚点

!**************** I direction ****************
	do k=1,nk
	do j=1,nj

		do i=1,ni

!			uvwt_line(1,i) = u(i,j,k)
!			uvwt_line(2,i) = v(i,j,k)
!			uvwt_line(3,i) = w(i,j,k)
!			uvwt_line(4,i) = t(i,j,k)

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

!			do m=1,12
			do m=5,12
				duvwt_line(m,i) = duvwt(m,i,j,k)
			enddo

		enddo

		do i=0,ni+1
			uvwt_lvir(1,i) = u(i,j,k)
			uvwt_lvir(2,i) = v(i,j,k)
			uvwt_lvir(3,i) = w(i,j,k)
			uvwt_lvir(4,i) = t(i,j,k)
		enddo

		call VALUE_line_half_2nd_virtual(4,nmax,ni,uvwt_lvir,uvwt_half ) !求原始变量在半结点上的值
		call VALUE_HALF_NODE_2nd_TEM(0,1,1,12,nmax,ni ,duvwt_line,duvwt_half) !把节点上的一阶导数插值到半结点
		call VALUE_HALF_NODE(9 ,nmax,ni,kxyz_line ,kxyz_half )   !把网格导数插值到半结点
		call VALUE_HALF_NODE(1 ,nmax,ni,vol_line  ,vol_half  )
		call VALUE_HALF_2nd(1 ,nmax,ni,vslt1_line,vslt1_half  ) !*坐标：1:ni --> 0:ni
		call VALUE_HALF_2nd(1 ,nmax,ni,vslt2_line,vslt2_half  )

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

!			uvwt_line(1,j) = u(i,j,k)
!			uvwt_line(2,j) = v(i,j,k)
!			uvwt_line(3,j) = w(i,j,k)
!			uvwt_line(4,j) = t(i,j,k)

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

			do m=1,4
				duvwt_line(m,j) = duvwt(m,i,j,k)
			enddo
			do m=9,12
				duvwt_line(m,j) = duvwt(m,i,j,k)
			enddo

		enddo

		do j=0,nj+1
			uvwt_lvir(1,j) = u(i,j,k)
			uvwt_lvir(2,j) = v(i,j,k)
			uvwt_lvir(3,j) = w(i,j,k)
			uvwt_lvir(4,j) = t(i,j,k)
		enddo

		call VALUE_line_half_2nd_virtual(4,nmax,nj,uvwt_lvir,uvwt_half )
		call VALUE_HALF_NODE_2nd_tem(1,0,1,12,nmax,nj,duvwt_line,duvwt_half)
		call VALUE_HALF_NODE(9 ,nmax,nj,kxyz_line ,kxyz_half )
		call VALUE_HALF_NODE(1 ,nmax,nj,vol_line  ,vol_half  )

		call VALUE_HALF_2nd(1 ,nmax,nj,vslt1_line,vslt1_half  ) !*坐标：1:nj --> 0:nj
		call VALUE_HALF_2nd(1 ,nmax,nj,vslt2_line,vslt2_half  )

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

!			uvwt_line(1,k) = u(i,j,k)
!			uvwt_line(2,k) = v(i,j,k)
!			uvwt_line(3,k) = w(i,j,k)
!			uvwt_line(4,k) = t(i,j,k)

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

!			do m=1,12
			do m=1,8
				duvwt_line(m,k) = duvwt(m,i,j,k)
			enddo
		enddo

		do k=0,nk+1
			uvwt_lvir(1,k) = u(i,j,k)
			uvwt_lvir(2,k) = v(i,j,k)
			uvwt_lvir(3,k) = w(i,j,k)
			uvwt_lvir(4,k) = t(i,j,k)
		enddo

		call VALUE_line_half_2nd_virtual(4,nmax,nk,uvwt_lvir,uvwt_half )
		call VALUE_HALF_NODE_2nd_tem(1,1,0,12,nmax,nk,duvwt_line,duvwt_half)
		call VALUE_HALF_NODE(9 ,nmax,nk,kxyz_line ,kxyz_half )
		call VALUE_HALF_NODE(1 ,nmax,nk,vol_line  ,vol_half  )
		call VALUE_HALF_2nd(1 ,nmax,nk,vslt1_line,vslt1_half  ) !*坐标：1:nk --> 0:nk
		call VALUE_HALF_2nd(1 ,nmax,nk,vslt2_line,vslt2_half  )

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
end subroutine WCNSE5_2nd_VIS_virtual_gcl
!=============================================================================!
!=============================================================================!
!=============================================================================!
!=============================================================================!

subroutine UVWT_DER_2nd
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt
  implicit none
	integer :: i,j,k,m
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

			call DUVWT_NODE_LINE_2nd(nmax,ni,4,uvwt,dtem) !*TGH. 暂时没有用虚点值
!			call DUVWT_NODE_LINE_virtual(nmax,ni,4,uvwt,dtem) !*TGH. 要用虚点值

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

			call DUVWT_NODE_LINE_2nd(nmax,nj,4,uvwt,dtem)
!			call DUVWT_NODE_LINE_virtual(nmax,nj,4,uvwt,dtem)

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

			call DUVWT_NODE_LINE_2nd(nmax,nk,4,uvwt,dtem)
!			call DUVWT_NODE_LINE_virtual(nmax,nk,4,uvwt,dtem)

			do k=1,nk
				do m=1,4

					duvwt(m+8,i,j,k)=dtem(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_2nd

!=============================================================================!
!=============================================================================!
subroutine DUVWT_NODE_LINE_2nd(nmax,ni,n1,uvwt,duvwt)  !___DERINODE
!---------------------------------------------------------------!
!* 计算一条线上结点的1阶导数，不用虚点
!---------------------------------------------------------------!
	use global_variables,only : nijk2nd
    use define_precision_mod
    implicit none

!	real(prec) :: AC2,BC2,CC2,DC2,EC2 

	integer :: nmax,i,m,ni,n1,n2
	real(prec) :: uvwt(n1,-1:nmax+1),duvwt(n1,nmax)

		do i=2,ni-1
		do m=1,n1
			duvwt(m,i) = 0.5*(uvwt(m,i+1) - uvwt(m,i-1))
		enddo
		enddo

		do m=1,n1
			duvwt(m,1) = uvwt(m,2 ) - uvwt(m,1)
			duvwt(m,ni)= uvwt(m,ni) - uvwt(m,ni-1)
		enddo
	
end subroutine DUVWT_NODE_LINE_2nd

!=============================================================================!
!=============================================================================!

subroutine UVWT_DER_2nd_half_virtual    !____直接求半结点1阶导数 2009.2.1
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt_mid
  implicit none
	integer :: i,j,k,m
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

			call DUVWT_half_line_2nd_virtual(nmax,ni,4,uvwt,duvwt)     !*TGH. 要用虚点上的值

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

			call DUVWT_half_line_2nd_virtual(nmax,nj,4,uvwt,duvwt)     !*TGH. 要用虚点上的值

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

			call DUVWT_half_line_2nd_virtual(nmax,nk,4,uvwt,duvwt)     !*TGH. 要用虚点上的值

			do k=0,nk
				do m=1,4

					duvwt_mid(m+8,i,j,k)=duvwt(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_2nd_half_virtual

!=============================================================================!
!=============================================================================!

subroutine DUVWT_half_line_2nd_virtual(nmax,ni,n1,uvwt,duvwt)
!---------------------------------------------------------------!
!* 计算一条线上半结点的1阶导数，要用虚点
!---------------------------------------------------------------!
	use global_variables,only : nijk2nd
    use define_precision_mod
    implicit none

	integer :: nmax,i,m,ni,n1,n2
	real(prec) :: uvwt(n1,-1:nmax+1),duvwt(n1,0:nmax)

		do m=1,n1
		do i=0,ni
				duvwt(m,i) = uvwt(m,i+1) - uvwt(m,i)
		enddo
		enddo

  return
end subroutine DUVWT_half_line_2nd_virtual

!=============================================================================!
subroutine VALUE_line_half_2nd_virtual(n,nmax,ni,q,q_half)
!---------------------------------------------------------------!
!* 插值一条线上半结点的值，要用虚点
!---------------------------------------------------------------!
	use global_variables,only : nijk2nd
    use define_precision_mod
	implicit none

	real(prec) A1,B1,dd16
	real(prec) A2,B2,C2,D2
	integer :: nmax,n,ni,m,i
	real    :: q(n,0:nmax+1),q_half(n,0:nmax)


		do i=0,ni
		do m=1,n
				q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
		enddo
		enddo


	return
end subroutine VALUE_line_half_2nd_virtual

!=============================================================================!
!=============================================================================!
!=============================================================================!
subroutine VALUE_HALF_NODE_2nd_TEM(IP,JP,KP,n,nmax,ni,q,q_half) !___INTERV2
!----------------------------------------------------------------------------!
!*TGH 从节点差值半结点 (N,1:NI) --> (N,0:NI)                                 !
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

	return
end subroutine VALUE_HALF_NODE_2nd_tem
!=============================================================================!
!=============================================================================!



