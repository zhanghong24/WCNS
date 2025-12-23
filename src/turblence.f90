!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine turbulence is used to simulate the turbulent effects         !
!	Called by                                                                   !
!			Subroutine inviscd3d in Convect.f90                                     !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			April,2006                                                              !
! Reference                                                                   !
!     Chen Liangzhong's master thesis 
!  Modified by TU Guohua                                                      !
!   Checked by TU Guohua                                                      !
!-----------------------------------------------------------------------------!
subroutine turbulence(nb)
	use global_variables,only:nlamtur,nwallfun,ntrans,nwtmax &
	                         ,nameturb,ndes,mb_dgrid,dgrid
	implicit none

	integer :: nb

	if(nlamtur ==0)then

		call algebric_model_turbulence

	else

	  if(ndes > 0)then
	    if( ndes ==11 .or. (nlamtur == 2 .and. nameturb=='SST') )dgrid => mb_dgrid(nb)%a3d
	  endif

    call differential_model_turbulence
!    call differential_model_turb_SA_GCL


	endif

	if(ntrans>=1)then
		call transition_turbulence
	endif

	if(nwtmax==1)then
		call max_vist
	endif

	call max_control_turbulence

  return
end subroutine turbulence
!-----------------------------------------------------------------------------!
subroutine algebric_model_turbulence
	use global_variables,only:nlamtur
	implicit none
	
	call cal_vist_bl

	return
end subroutine algebric_model_turbulence
!-----------------------------------------------------------------------------!
subroutine differential_model_turbulence
	use global_variables,only:method,ncmpcor,nameturb
	implicit none

	call allocate_variables_temp_turbulence

	call set_dq_to_zero_turbulence
	
!!	call debug_SA

	call cal_rhs_source_turbulence
	

	if( ncmpcor <= 9 .and. nameturb=='SA' )then
	  call cal_rhs_viscous_turbulence_dif_r
	  call cal_spec_turbulence_r
	  ! cal_rhs_viscous_turbulence_dif_r ��cal_rhs_viscous_turbulence_dif�����һ���ܶ�
	  ! cal_cal_spec_turbulence_r �ڿ�������ճ���װ뾶ʱû�г����ܶȣ���Ϊviseq=qke��������r*qke
!	elseif(nameturb =="SST")then
!	  call cal_rhs_viscous_turbulence_dif
!	  call cal_spec_turbulence_sst_tem
	else
	  call cal_rhs_viscous_turbulence_dif
	  call cal_spec_turbulence
	endif

	call cal_rhs_inviscous_turbulence


	if( ncmpcor >=30 .and. ncmpcor<=39 .and. nameturb=='SA')then
	  call viscosity_equiv_lamturb
	  call cal_rhs_viscous_ncmpcor30
	endif

	call lusgs_turbulence

	call update_turbulence

!	call residual_turbulence

	call deallocate_variables_temp_turbulence

	return
end subroutine differential_model_turbulence
subroutine debug_SA
	use global_variables,only:method,ni,nj,nk,qke,dqke,r,visl,vist,nbself,ds_turbulence
	implicit none
	real,parameter :: cv13  = 7.1**3.0
	integer :: i,j,k,nb,cmethod
	real :: fv1,kafan,kf3,romu
  
	cmethod = method - 1
	
	if (nbself == 1) then

	do k = 1,nk+cmethod
	  do j = 1,nj+cmethod
	    do i = 1,ni+cmethod
		   write(52,'(3(1x,i3),4(1x,e12.5))') i,j,k,qke(i,j,k,1),vist(i,j,k),visl(i,j,k),r(i,j,k)

	     enddo
	   enddo
	enddo
	
	write(52,*) "--------------------------------------------------------------------"
	write(52,*) "--------------------------------------------------------------------"
	write(52,*) "--------------------------------------------------------------------"
	!!close(52)
	stop
	
	end if

	return
end subroutine debug_SA
!-----------------------------------------------------------------------------!
subroutine allocate_variables_temp_turbulence
	use global_variables,only:spec,dmudxyz,ni,nj,nk,nlamtur,dqke,viseq, &
	                          drdxyz,dudxyz,dvdxyz,dwdxyz
	implicit none

	integer :: npdi,idimp,jdimp,kdimp

	npdi  = -1
	idimp = ni+1
	jdimp = nj+1
	kdimp = nk+1

  allocate( dqke   (npdi:idimp,npdi:jdimp,npdi:kdimp,nlamtur) )
  allocate( spec   (npdi:idimp,npdi:jdimp,npdi:kdimp,nlamtur) )
  allocate( dmudxyz(ni,nj,nk,3,nlamtur) )
  allocate( drdxyz (ni,nj,nk,3) )
  allocate( dudxyz (ni,nj,nk,3) )
  allocate( dvdxyz (ni,nj,nk,3) )
  allocate( dwdxyz (ni,nj,nk,3) )
  allocate( viseq  (ni,nj,nk,  nlamtur) )

	return
end subroutine allocate_variables_temp_turbulence
!-----------------------------------------------------------------------------!
subroutine deallocate_variables_temp_turbulence
	use global_variables,only:spec,dmudxyz,dqke,viseq,drdxyz,dudxyz,dvdxyz,dwdxyz,timedt_turb
	implicit none

  deallocate( dqke   )
  deallocate( spec   )
  deallocate( dmudxyz)
  deallocate( drdxyz )
  deallocate( dudxyz )
  deallocate( dvdxyz )
  deallocate( dwdxyz )
  deallocate( viseq  )

	return
end subroutine deallocate_variables_temp_turbulence
!-----------------------------------------------------------------------------!
subroutine cal_spec_turbulence
	use global_variables
	implicit none

	integer :: i,j,k,m,cmethod,ni1,nj1,nk1
	real    :: unx,uny,unz,vis,nx,ny,nz,re2

	cmethod = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

	re2     = re*2.0

	do m=1,nlamtur
		do k=cmethod,nk1
			do j=cmethod,nj1
				do i=cmethod,ni1

		      unx = u(i,j,k)*kcx(i,j,k)+v(i,j,k)*kcy(i,j,k)+w(i,j,k)*kcz(i,j,k)
			  uny = u(i,j,k)*etx(i,j,k)+v(i,j,k)*ety(i,j,k)+w(i,j,k)*etz(i,j,k)
		  	  unz = u(i,j,k)*ctx(i,j,k)+v(i,j,k)*cty(i,j,k)+w(i,j,k)*ctz(i,j,k)

!			    spec(i,j,k,m) = spec(i,j,k,m)+abs(unx)+abs(uny)+abs(unz)+ 1.0/dtdt(i,j,k) !old
			    spec(i,j,k,m) = spec(i,j,k,m)+abs(unx)+abs(uny)+abs(unz)+ 1.0/(dtdt(i,j,k)*timedt_turb) !new �Ӵ�ʱ�䲽��
					
					nx = kcx(i,j,k)*kcx(i,j,k)+kcy(i,j,k)*kcy(i,j,k)+kcz(i,j,k)*kcz(i,j,k)
					ny = etx(i,j,k)*etx(i,j,k)+ety(i,j,k)*ety(i,j,k)+etz(i,j,k)*etz(i,j,k)
					nz = ctx(i,j,k)*ctx(i,j,k)+cty(i,j,k)*cty(i,j,k)+ctz(i,j,k)*ctz(i,j,k)

					vis = viseq(i,j,k,m)/vol(i,j,k)/r(i,j,k)
					
					spec(i,j,k,m) = spec(i,j,k,m)+re2*vis*(nx+ny+nz)

				enddo
			enddo
		enddo
	enddo

	return
end subroutine cal_spec_turbulence
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_spec_turbulence_r
! ��cal_spec_turbulence��ȣ���viseq�ٳ����ܶ���
	use global_variables
	implicit none

	integer :: i,j,k,m,cmethod,ni1,nj1,nk1
	real    :: unx,uny,unz,vis,nx,ny,nz,re2

	cmethod = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

	re2     = re*2.0

	do m=1,nlamtur
		do k=cmethod,nk1
			do j=cmethod,nj1
				do i=cmethod,ni1

		      unx = u(i,j,k)*kcx(i,j,k)+v(i,j,k)*kcy(i,j,k)+w(i,j,k)*kcz(i,j,k)
			  uny = u(i,j,k)*etx(i,j,k)+v(i,j,k)*ety(i,j,k)+w(i,j,k)*etz(i,j,k)
		  	  unz = u(i,j,k)*ctx(i,j,k)+v(i,j,k)*cty(i,j,k)+w(i,j,k)*ctz(i,j,k)

!			    spec(i,j,k,m) = spec(i,j,k,m)+abs(unx)+abs(uny)+abs(unz)+ 1.0/dtdt(i,j,k) !old
			    spec(i,j,k,m) = spec(i,j,k,m)+abs(unx)+abs(uny)+abs(unz)+ 1.0/(dtdt(i,j,k)*timedt_turb) !new �Ӵ�ʱ�䲽��
					
					nx = kcx(i,j,k)*kcx(i,j,k)+kcy(i,j,k)*kcy(i,j,k)+kcz(i,j,k)*kcz(i,j,k)
					ny = etx(i,j,k)*etx(i,j,k)+ety(i,j,k)*ety(i,j,k)+etz(i,j,k)*etz(i,j,k)
					nz = ctx(i,j,k)*ctx(i,j,k)+cty(i,j,k)*cty(i,j,k)+ctz(i,j,k)*ctz(i,j,k)

					vis = viseq(i,j,k,m)/vol(i,j,k)
					
					spec(i,j,k,m) = spec(i,j,k,m)+re2*vis*(nx+ny+nz)

				enddo
			enddo
		enddo
	enddo

	return
end subroutine cal_spec_turbulence_r
!-----------------------------------------------------------------------------!
!	Function                                                                    !
!			Subroutine set_dq_to_zero_turbulence is used to set zero for turblence  !
!	Called by                                                                   !
!			Subroutine turbulence                                     !
!	Editor                                                                      !    
!			Chen Liangzhong                                                         !
!	Date                                                                        !
!			April,2006                                                           !                                                          !
!-----------------------------------------------------------------------------!
subroutine set_dq_to_zero_turbulence
	use global_variables,only:nlamtur,dqke,spec,ni,nj,nk
	implicit none
	integer :: i,j,k,n

	do n=1,nlamtur
		do k=-1,nk+1
			do j=-1,nj+1
				do i=-1,ni+1

					dqke(i,j,k,n) = 0.0
					spec(i,j,k,n) = 0.0
					
				enddo
			enddo
		enddo
	enddo

	return
end subroutine set_dq_to_zero_turbulence
!-----------------------------------------------------------------------------!
subroutine cal_rhs_inviscous_turbulence
  use global_variables
  implicit none

  integer :: i,j,k,m,cmethd
  real :: fc(1:nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)

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

        call NND_turbulence(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

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

        call NND_turbulence(xk,xb,nj,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)

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
          call NND_turbulence(xk,xb,nk,nmax,nl,method,cmethd,efix,trxyz,q_line,fc)
          do k=cmethd,nk-1      
            dqke(i,j,k,m) = dqke(i,j,k,m) + fc(k)
          enddo
       enddo
    enddo
  enddo
  return
end subroutine cal_rhs_inviscous_turbulence
!_____________________________________________________________________!
subroutine NND_turbulence(xk,xb,ni,nmax,nl,method,cmethd,efix,trxyz,q_line,dfc)
    implicit none
    integer :: i,m,cmethd,ni,nmax,nl,method
    real,dimension( 1:nl ) :: dql,dqr
    real :: nx,ny,nz,xk,xb,efix,gamaeq,fl,fr
    real :: c1,c2,dfc(1:nmax)
    real :: fc(1:nmax),q_line(1:nl,-1:nmax+1),trxyz(3,nmax)
    real :: d_q(1:nl,0:nmax+1),d_qcl(1:nl,0:nmax),d_qcr(1:nl,0:nmax)
    real,external :: minmod

    c1 = 0.25 * ( 1.0 + xk )
    c2 = 0.25 * ( 1.0 - xk )


    do i=method,ni+1 
       do m=1,nl
           d_q(m,i) = q_line(m,i) - q_line(m,i-1)
       enddo
    enddo

!!----------------------------------------------------
! ������������������ģ����ճ����ɢʱ����������ϵ�ֵ
! ����ֱ�Ӳ������ֵ����Խ���һ���ֵ���������ȶ�
! ȥ���������
! 2009��12��29�գ�Ϳ����
!    do i=1,method
!       do m=1,nl
!          d_q(m,cmethd-1) = d_q(m,cmethd)
!          d_q(m,ni+1  )   = d_q(m,ni)
!       enddo
!    enddo
!!----------------------------------------------------

    do i=method,ni     
       do m=1,nl
          d_qcl(m,i) = minmod( xb*d_q(m,i+1), d_q(m,i) )
          d_qcr(m,i) = minmod( d_q(m,i+1), xb*d_q(m,i) )
       enddo
    enddo

    do i=cmethd,ni      
       nx = trxyz(1,i)
       ny = trxyz(2,i)
       nz = trxyz(3,i)
       do m=1,method
          nx = 0.50 * ( trxyz(1,i-1) + nx )
          ny = 0.50 * ( trxyz(2,i-1) + ny )
          nz = 0.50 * ( trxyz(3,i-1) + nz )
       enddo
       do m=1,nl
          dql(m) = q_line(m,i-1) + c1 * d_qcr(m,i-1) + c2 * d_qcl(m,i-1)
          dqr(m) = q_line(m,i  ) - c1 * d_qcl(m,i  ) - c2 * d_qcr(m,i  )
       enddo

       call flux_turbulence(dql,nx,ny,nz,fl, 1)
       call flux_turbulence(dqr,nx,ny,nz,fr,-1)

       fc(i) = fl + fr
    enddo

    do i=cmethd,ni-1     
       dfc(i) =  fc(i+1) - fc(i)
    enddo
    return
end subroutine NND_turbulence 
!_____________________________________________________________________!
subroutine flux_turbulence(prim,nx,ny,nz,f,npn)
    implicit none
    integer :: npn
    real :: prim(5)
    real :: l1,f
    real :: rm,um,vm,wm,qkem
    real :: nx,ny,nz

    rm   = prim(1)
    um   = prim(2)
    vm   = prim(3)
    wm   = prim(4)
    qkem = prim(5)

    l1  = nx*um + ny*vm + nz*wm

    l1 = 0.5*(l1 + npn*sqrt(l1*l1))


    f = l1 * qkem * rm

    return
end subroutine flux_turbulence
!-----------------------------------------------------------------------------!
subroutine cal_rhs_source_turbulence
	use global_variables
	implicit none

	if(nameturb=='SA')then
		call cal_rhs_source_SA
	elseif(nameturb=='SST')then
		call cal_rhs_source_SST
	else
		write(*,*)trim(nameturb),'has not been realized in this release version'
		stop
	endif

	return
end subroutine cal_rhs_source_turbulence
!-----------------------------------------------------------------------------!
subroutine lusgs_turbulence
  use global_variables , only : ni,nj,nk,dqke,spec,nlamtur,method,dtdt,nl
  implicit none
  integer :: i,j,k,m,cmth,ni1,nj1,nk1,i1,j1,k1,n1
  real    :: prim_i(nl),prim_j(nl),prim_k(nl)
  real    :: gykb_i(4),gykb_j(4),gykb_k(4)
  real    :: dq_i,dq_j,dq_k,de,df,dg,rhs0
	cmth = method+1

	ni1 = ni-1
	nj1 = nj-1
	nk1 = nk-1
	n1  = -1

	do m=1,nlamtur
	  do k=cmth,nk1
       do j=cmth,nj1
          do i=cmth,ni1

						i1 = i-1
						j1 = j-1
						k1 = k-1

            dq_i = dqke(i1,j,k,m)
            dq_j = dqke(i,j1,k,m)
            dq_k = dqke(i,j,k1,m)

            call getprim(i1,j,k,prim_i)
            call getprim(i,j1,k,prim_j)
            call getprim(i,j,k1,prim_k)

            call getrkec_mml(i-method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
            call getrkec_mml(i,j-method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
            call getrkec_mml(i,j,k-method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

            call mxdq_std_turbulence(prim_i,gykb_i,dq_i,de,1)
            call mxdq_std_turbulence(prim_j,gykb_j,dq_j,df,1)
            call mxdq_std_turbulence(prim_k,gykb_k,dq_k,dg,1)

            rhs0 = de + df + dg

            dqke(i,j,k,m) = (rhs0-dqke(i,j,k,m))/spec(i,j,k,m)

          enddo
       enddo
    enddo

    do k=nk1,cmth,-1
       do j=nj1,cmth,-1
          do i=ni1,cmth,-1

						i1 = i+1
						j1 = j+1
						k1 = k+1

            dq_i = dqke(i1,j,k,m)
            dq_j = dqke(i,j1,k,m)
            dq_k = dqke(i,j,k1,m)

            call getprim(i1,j,k,prim_i)
            call getprim(i,j1,k,prim_j)
            call getprim(i,j,k1,prim_k)

            call getrkec_mml(i1,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
            call getrkec_mml(i,j1,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
            call getrkec_mml(i,j,k1,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

            call mxdq_std_turbulence(prim_i,gykb_i,dq_i,de,n1)
            call mxdq_std_turbulence(prim_j,gykb_j,dq_j,df,n1)
            call mxdq_std_turbulence(prim_k,gykb_k,dq_k,dg,n1)

             rhs0 = de + df + dg

             dqke(i,j,k,m) = dqke(i,j,k,m)-rhs0/spec(i,j,k,m)

          enddo
       enddo
	  enddo
	enddo

	return
end subroutine lusgs_turbulence
!-----------------------------------------------------------------------------!
subroutine mxdq_std_turbulence(prim,gykb,dq,f,npn)
	implicit none

  integer :: npn

  real :: f,prim(5),gykb(3),dq
  real :: nx,ny,nz,l1
  real :: um,vm,wm,pm
  
  um = prim(2)
  vm = prim(3)
  wm = prim(4)

  nx = gykb(1)
  ny = gykb(2)
  nz = gykb(3)

  l1 = nx*um + ny*vm + nz*wm
  l1 = 0.5 * ( l1 + npn*(abs(l1)))
 
  f = l1 * dq 

  return
end subroutine mxdq_std_turbulence
!-----------------------------------------------------------------------------!
subroutine update_turbulence
	use global_variables,only:reynolds,r,u,v,w,qke,dqke,nlamtur,ni,nj,nk,method, &
	                          nameturb,kmaxlim,muoo,koo,goo,omgoo
	implicit none
	real,parameter :: minlimiter=1.0e-20
	integer :: nb,mi,i,j,k,cmethod,ni1,nj1,nk1
	real :: parameter_wall(nlamtur),parameter_inflow(nlamtur),up,vp,wp,uvw2

	if(nameturb == 'SA')then
		parameter_inflow(1) = muoo
	elseif(nameturb == 'SST')then
		parameter_inflow(1) = koo
		parameter_inflow(2) = omgoo
	endif

	cmethod = 1  + method
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

	do mi =1,nlamtur
	  do k = cmethod,nk1
	    do j = cmethod,nj1
	      do i = cmethod,ni1
	        qke(i,j,k,mi) = amax1(qke(i,j,k,mi)+dqke(i,j,k,mi)/r(i,j,k) ,& 
					                      minlimiter*parameter_inflow(mi)) 
	      enddo
	    enddo
	  enddo
	enddo

	if (nameturb == 'SA' )then
		call cal_vist_SA
	elseif (nameturb == 'SST' )then
		call cal_vist_SST
	endif

	return
end subroutine update_turbulence
!-----------------------------------------------------------------------------!
subroutine max_control_turbulence
	use global_variables,only:ni,nj,nk,vist,visl,method,kmaxlim
	implicit none

	integer :: i,j,k,cmethod,ni1,nj1,nk1
	real    :: min_t,max_t

	cmethod = method-1
	ni1     = ni+cmethod
	nj1     = nj+cmethod
	nk1     = nk+cmethod

!	min_t = 1.0e-5
!
!	do k=1,nk1
!	  do j=1,nj1
!	    do i=1,ni1
!				
!	      vist(i,j,k) = amax1(vist(i,j,k),min_t*visl(i,j,k))
!	    enddo
!	  enddo
!	enddo

	if(kmaxlim>1.0e-15)then
	  do k = cmethod,nk1
	    do j = cmethod,nj1
	      do i = cmethod,ni1
	        vist(i,j,k) = amin1(vist(i,j,k),kmaxlim*visl(i,j,k)) 
	      enddo
	    enddo
	  enddo
	endif


	return
end subroutine max_control_turbulence
!-----------------------------------------------------------------------------!
subroutine residual_turbulence
	use global_variables
	implicit none
	real                 :: temp1
	integer              :: mi,i,j,k,cmethod,itemp,jtemp,ktemp,ngrid,ni1,nj1,nk1
	real,dimension(2)    :: residual,temprofai,tempdrofai
	integer,dimension(2) :: tempirofai,tempjrofai,tempkrofai, &
	                        tempidrofai,tempjdrofai,tempkdrofai
  
	cmethod = method + 1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

	do mi =1,nlamtur
	   temp1           = 0.0
		 temprofai(mi)   = 0.0
		 tempdrofai(mi)  = 0.0
	   tempirofai(mi)  = -1
		 tempjrofai(mi)  = -1
		 tempkrofai(mi)  = -1
	   tempidrofai(mi) = -1
		 tempjdrofai(mi) = -1
		 tempkdrofai(mi) = -1

     residual(mi)    = 0.0

	   do k =cmethod, nk1
	     do j =cmethod, nj1
	        do i =cmethod, ni1
		        residual(mi) = residual(mi) + (dqke(i,j,k,mi)/dtdt(i,j,k)/vol(i,j,k))**2

		        if(temp1<vist(i,j,k))then
			        temp1 = vist(i,j,k)
			        itemp = i
							jtemp = j
							ktemp = k
		        endif

		        if(temprofai(mi)<qke(i,j,k,mi))then
			        temprofai(mi)  = qke(i,j,k,mi)
			        tempirofai(mi) = i
							tempjrofai(mi) = j
							tempkrofai(mi) = k
		        endif

		        if(tempdrofai(mi)<abs(dqke(i,j,k,mi)))then
			        tempdrofai(mi)  = dqke(i,j,k,mi)
			        tempidrofai(mi) = i
						  tempjdrofai(mi) = j
						  tempkdrofai(mi) = k
		        endif

	         enddo
	      enddo
	   enddo
	enddo

	ngrid = (ni-cmethod) * (nj-cmethod) * (nk-cmethod)
	residual(1:nlamtur) = sqrt(residual(1:nlamtur))/ngrid

	write(*,*)
	write(*,*)'vtmax',temp1,itemp,jtemp,ktemp

	do mi =1,nlamtur
	   write(*,*)'max  ',temprofai(mi),tempirofai(mi),tempjrofai(mi),tempkrofai(mi)
	   write(*,*)'dmax ',tempdrofai(mi),tempidrofai(mi),tempjdrofai(mi),tempkdrofai(mi)
	   write(*,*)'res  ',residual(mi)
	enddo

!  open(1,file='Turb_max_res_dmax.dat',access='append',status='unknown')

!  if(nout == 1 )then
!	   rewind(1)
!	   write(1,*)'��������  ���vist  ���������  �в�   ���ֵ'
!	endif

!	if(nlamtur == 2)then
!	   write(1,10)nout,temp1,tempdrofai(1),residual(1),temprofai(1)
!   else
!	   write(1,20)nout,temp1,tempdrofai(:),residual(:),temprofai(:)
!    endif

!10  format(1x,i6,2x,4(e,2x))
!20  format(1x,i6,2x,7(e,2x))

!    close(1)

	return
end subroutine residual_turbulence
!-----------------------------------------------------------------------------!
subroutine dfidxyz_center(nsub,nplus,ni,nj,nk,fi,dfidxyz,method,nv)
	use global_variables,only:kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol
	implicit none

	integer :: nsub,nplus,ni,nj,nk,method,i,j,k,n,nv,m,ni1,nj1,nk1
	real    :: ddx,ddy,ddz,vol2
	real    :: dkc,det,dct
	real    :: fim,fjm,fkm
	real    :: fi(nsub:ni+nplus,nsub:nj+nplus,nsub:nk+nplus,nv)
	real    :: dfidxyz(ni,nj,nk,3,nv)
!	real,pointer,dimension(:,:,:,:) :: work

!	allocate(work(ni,nj,nk,9))

	do m=1,nv
!		do n=method,0  ! finite volume method
!			do k=1,nk
!				do j=1,nj
!					do i=1,ni
!
!						fim = fi(i,j,k,m)+fi(i-1,j,k,m)
!						fjm = fi(i,j,k,m)+fi(i,j-1,k,m)
!						fkm = fi(i,j,k,m)+fi(i,j,k-1,m)
!
!						work(i,j,k,1) = fim*kcx(i,j,k)
!						work(i,j,k,2) = fjm*etx(i,j,k)
!						work(i,j,k,3) = fkm*ctx(i,j,k)
!
!						work(i,j,k,4) = fim*kcy(i,j,k)
!						work(i,j,k,5) = fjm*ety(i,j,k)
!						work(i,j,k,6) = fkm*cty(i,j,k)
!
!						work(i,j,k,7) = fim*kcz(i,j,k)
!						work(i,j,k,8) = fjm*etz(i,j,k)
!						work(i,j,k,9) = fkm*ctz(i,j,k)
!
!					enddo
!				enddo
!			enddo
!
!			ni1     = ni-1
!			nj1     = nj-1
!			nk1     = nk-1
!
!			do k=1,nk1
!				do j=1,nj1
!					do i=1,ni1
!
!						vol2 = 0.5/vol(i,j,k)
!
!						dfidxyz(i,j,k,1,m) = ((work(i+1,j,k,1)-work(i,j,k,1)) + &
!																  (work(i,j+1,k,2)-work(i,j,k,2)) + &
!																  (work(i,j,k+1,3)-work(i,j,k,3)))*vol2
!
!						dfidxyz(i,j,k,2,m) = ((work(i+1,j,k,4)-work(i,j,k,4)) + &
!																  (work(i,j+1,k,5)-work(i,j,k,5)) + &
!																  (work(i,j,k+1,6)-work(i,j,k,6)))*vol2
!
!						dfidxyz(i,j,k,3,m) = ((work(i+1,j,k,7)-work(i,j,k,7)) + &
!																  (work(i,j+1,k,8)-work(i,j,k,8)) + &
!																  (work(i,j,k+1,9)-work(i,j,k,9)))*vol2
!
!					enddo
!				enddo
!			enddo
!		enddo

		do n=1,method  ! finite difference method
			do k=1,nk
				do j=1,nj
					do i=1,ni

						dkc = fi(i+1,j,k,m)-fi(i-1,j,k,m)
						det = fi(i,j+1,k,m)-fi(i,j-1,k,m)
						dct = fi(i,j,k+1,m)-fi(i,j,k-1,m)

						ddx = dkc*kcx(i,j,k)+det*etx(i,j,k)+dct*ctx(i,j,k)
						ddy = dkc*kcy(i,j,k)+det*ety(i,j,k)+dct*cty(i,j,k)
						ddz = dkc*kcz(i,j,k)+det*etz(i,j,k)+dct*ctz(i,j,k)

						vol2 =0.5/vol(i,j,k)

						dfidxyz(i,j,k,1,m) = ddx*vol2
						dfidxyz(i,j,k,2,m) = ddy*vol2
						dfidxyz(i,j,k,3,m) = ddz*vol2

					enddo
				enddo
			enddo
		enddo

	enddo

!	deallocate(work)

	return
end subroutine dfidxyz_center
!-----------------------------------------------------------------------------!
subroutine cal_rhs_viscous_turbulence_dif
  use global_variables
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real :: ev_l,vis,nx,ny,nz
	real :: dmudx,dmudz,dmudy,Dvis
  real,dimension( nmax+1,nlamtur ) :: fv

	if(ncmpcor>=20 .and. ncmpcor<=29 .and. nameturb=='SST')then
		call cal_rhs_viscous_SST_ncmpcor20
		return
	endif

	cmthd = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1
			do i=cmthd,ni

				do n=method,0   !  FV method

					im = min(i,ni1)
					in = max(1,i-1)

					nx = kcx(i,j,k)
					ny = kcy(i,j,k)
					nz = kcz(i,j,k)

				enddo

				do n=1,method   !  FD method

					im = i
					in = i-1

					nx = 0.5*(kcx(i,j,k)+kcx(in,j,k))
					ny = 0.5*(kcy(i,j,k)+kcy(in,j,k))
					nz = 0.5*(kcz(i,j,k)+kcz(in,j,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(im,j,k,1,m)+dmudxyz(in,j,k,1,m) 
					dmudy =  dmudxyz(im,j,k,2,m)+dmudxyz(in,j,k,2,m)
					dmudz =  dmudxyz(im,j,k,3,m)+dmudxyz(in,j,k,3,m)

					vis  = viseq(im,j,k,m) + viseq(in,j,k,m)

					fv(i,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do i=cmthd,ni1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(i+1,m)-fv(i,m))
				enddo
			enddo
    enddo
  enddo

!------------
! J direction

	do k=cmthd,nk1
		do i=cmthd,ni1
			do j=cmthd,nj

				do n=method,0   !  FV method

					jm = min(j,nj1)
					jn = max(1,j-1)

					nx = etx(i,j,k)
					ny = ety(i,j,k)
					nz = etz(i,j,k)

				enddo

				do n=1,method   !  FD method

					jm = j
					jn = j-1

					nx = 0.5*(etx(i,j,k)+etx(i,jn,k))
					ny = 0.5*(ety(i,j,k)+ety(i,jn,k))
					nz = 0.5*(etz(i,j,k)+etz(i,jn,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,jm,k,1,m)+dmudxyz(i,jn,k,1,m) 
					dmudy =  dmudxyz(i,jm,k,2,m)+dmudxyz(i,jn,k,2,m)
					dmudz =  dmudxyz(i,jm,k,3,m)+dmudxyz(i,jn,k,3,m)

					vis  = viseq(i,jm,k,m) + viseq(i,jn,k,m)

					fv(j,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25

				enddo

			enddo

			do m =1,nlamtur
				do j=cmthd,nj1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(j+1,m)-fv(j,m))
				enddo
			enddo

    enddo
  enddo

!------------
! K direction
  do j=cmthd,nj1
		do i=cmthd,ni1
			do k=cmthd,nk

				do n=method,0   !  FV method

					km = min(k,nk1)
					kn = max(1,k-1)

					nx = ctx(i,j,k)
					ny = cty(i,j,k)
					nz = ctz(i,j,k)

				enddo

				do n=1,method   !  FD method

					km = k
					kn = k-1

					nx = 0.5*(ctx(i,j,k)+ctx(i,j,kn))
					ny = 0.5*(cty(i,j,k)+cty(i,j,kn))
					nz = 0.5*(ctz(i,j,k)+ctz(i,j,kn))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,j,km,1,m)+dmudxyz(i,j,kn,1,m) 
					dmudy =  dmudxyz(i,j,km,2,m)+dmudxyz(i,j,kn,2,m)
					dmudz =  dmudxyz(i,j,km,3,m)+dmudxyz(i,j,kn,3,m)

					vis  = viseq(i,j,km,m) + viseq(i,j,kn,m)

					fv(k,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do k=cmthd,nk1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(k+1,m)-fv(k,m))
				enddo
			enddo

    enddo
  enddo

!------------
  return
end subroutine cal_rhs_viscous_turbulence_dif
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_rhs_viscous_turbulence_dif_r
  use global_variables
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real :: ev_l,vis,nx,ny,nz
	real :: dmudx,dmudz,dmudy
  real,dimension( nmax+1,nlamtur ) :: fv

	cmthd = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1
			do i=cmthd,ni

				do n=method,0   !  FV method

					im = min(i,ni1)
					in = max(1,i-1)

					nx = kcx(i,j,k)
					ny = kcy(i,j,k)
					nz = kcz(i,j,k)

				enddo

				do n=1,method   !  FD method

					im = i
					in = i-1

					nx = 0.5*(kcx(i,j,k)+kcx(in,j,k))
					ny = 0.5*(kcy(i,j,k)+kcy(in,j,k))
					nz = 0.5*(kcz(i,j,k)+kcz(in,j,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(im,j,k,1,m)+dmudxyz(in,j,k,1,m) 
					dmudy =  dmudxyz(im,j,k,2,m)+dmudxyz(in,j,k,2,m)
					dmudz =  dmudxyz(im,j,k,3,m)+dmudxyz(in,j,k,3,m)

					vis  = viseq(im,j,k,m) + viseq(in,j,k,m)

					fv(i,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do i=cmthd,ni1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*r(i,j,k)*(fv(i+1,m)-fv(i,m))
				enddo
			enddo
    enddo
  enddo

!------------
! J direction

	do k=cmthd,nk1
		do i=cmthd,ni1
			do j=cmthd,nj

				do n=method,0   !  FV method

					jm = min(j,nj1)
					jn = max(1,j-1)

					nx = etx(i,j,k)
					ny = ety(i,j,k)
					nz = etz(i,j,k)

				enddo

				do n=1,method   !  FD method

					jm = j
					jn = j-1

					nx = 0.5*(etx(i,j,k)+etx(i,jn,k))
					ny = 0.5*(ety(i,j,k)+ety(i,jn,k))
					nz = 0.5*(etz(i,j,k)+etz(i,jn,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,jm,k,1,m)+dmudxyz(i,jn,k,1,m) 
					dmudy =  dmudxyz(i,jm,k,2,m)+dmudxyz(i,jn,k,2,m)
					dmudz =  dmudxyz(i,jm,k,3,m)+dmudxyz(i,jn,k,3,m)

					vis  = viseq(i,jm,k,m) + viseq(i,jn,k,m)

					fv(j,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25

				enddo

			enddo

			do m =1,nlamtur
				do j=cmthd,nj1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*r(i,j,k)*(fv(j+1,m)-fv(j,m))
				enddo
			enddo

    enddo
  enddo

!------------
! K direction
  do j=cmthd,nj1
		do i=cmthd,ni1
			do k=cmthd,nk

				do n=method,0   !  FV method

					km = min(k,nk1)
					kn = max(1,k-1)

					nx = ctx(i,j,k)
					ny = cty(i,j,k)
					nz = ctz(i,j,k)

				enddo

				do n=1,method   !  FD method

					km = k
					kn = k-1

					nx = 0.5*(ctx(i,j,k)+ctx(i,j,kn))
					ny = 0.5*(cty(i,j,k)+cty(i,j,kn))
					nz = 0.5*(ctz(i,j,k)+ctz(i,j,kn))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,j,km,1,m)+dmudxyz(i,j,kn,1,m) 
					dmudy =  dmudxyz(i,j,km,2,m)+dmudxyz(i,j,kn,2,m)
					dmudz =  dmudxyz(i,j,km,3,m)+dmudxyz(i,j,kn,3,m)

					vis  = viseq(i,j,km,m) + viseq(i,j,kn,m)

					fv(k,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do k=cmthd,nk1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*r(i,j,k)*(fv(k+1,m)-fv(k,m))
				enddo
			enddo

    enddo
  enddo

!------------
  return
end subroutine cal_rhs_viscous_turbulence_dif_r
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_rhs_vis_turb_dif_positive
! �����˷�ֹ��ɢ��Ϊ������

  use global_variables
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real :: ev_l,vis,nx,ny,nz
	real :: dmudx,dmudz,dmudy,Dvis
  real,dimension( nmax+1,nlamtur ) :: fv
  real,pointer :: ftem(:,:,:,:)
  integer :: ncount

	if(ncmpcor>=20 .and. ncmpcor<=29 .and. nameturb=='SST')then
		call cal_rhs_viscous_SST_ncmpcor20 !�˴�Ӧ���޸ģ���ʱû�и�
!		call cal_rhs_viscous_SST_ncmpcor20_positive
		return
	endif

	allocate(ftem(ni,nj,nk,nlamtur))
	ncount = 0

	cmthd = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1
			do i=cmthd,ni

				do n=method,0   !  FV method

					im = min(i,ni1)
					in = max(1,i-1)

					nx = kcx(i,j,k)
					ny = kcy(i,j,k)
					nz = kcz(i,j,k)

				enddo

				do n=1,method   !  FD method

					im = i
					in = i-1

					nx = 0.5*(kcx(i,j,k)+kcx(in,j,k))
					ny = 0.5*(kcy(i,j,k)+kcy(in,j,k))
					nz = 0.5*(kcz(i,j,k)+kcz(in,j,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(im,j,k,1,m)+dmudxyz(in,j,k,1,m) 
					dmudy =  dmudxyz(im,j,k,2,m)+dmudxyz(in,j,k,2,m)
					dmudz =  dmudxyz(im,j,k,3,m)+dmudxyz(in,j,k,3,m)

					vis  = viseq(im,j,k,m) + viseq(in,j,k,m)

					fv(i,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do i=cmthd,ni1
					ftem(i,j,k,m) = fv(i+1,m)-fv(i,m)
				enddo
			enddo
    enddo
  enddo

!------------
! J direction

	do k=cmthd,nk1
		do i=cmthd,ni1
			do j=cmthd,nj

				do n=method,0   !  FV method

					jm = min(j,nj1)
					jn = max(1,j-1)

					nx = etx(i,j,k)
					ny = ety(i,j,k)
					nz = etz(i,j,k)

				enddo

				do n=1,method   !  FD method

					jm = j
					jn = j-1

					nx = 0.5*(etx(i,j,k)+etx(i,jn,k))
					ny = 0.5*(ety(i,j,k)+ety(i,jn,k))
					nz = 0.5*(etz(i,j,k)+etz(i,jn,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,jm,k,1,m)+dmudxyz(i,jn,k,1,m) 
					dmudy =  dmudxyz(i,jm,k,2,m)+dmudxyz(i,jn,k,2,m)
					dmudz =  dmudxyz(i,jm,k,3,m)+dmudxyz(i,jn,k,3,m)

					vis  = viseq(i,jm,k,m) + viseq(i,jn,k,m)

					fv(j,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25

				enddo

			enddo

			do m =1,nlamtur
				do j=cmthd,nj1
					ftem(i,j,k,m)  = ftem(i,j,k,m) + (fv(j+1,m)-fv(j,m))
				enddo
			enddo

    enddo
  enddo

!------------
! K direction
  do j=cmthd,nj1
		do i=cmthd,ni1
			do k=cmthd,nk

				do n=method,0   !  FV method

					km = min(k,nk1)
					kn = max(1,k-1)

					nx = ctx(i,j,k)
					ny = cty(i,j,k)
					nz = ctz(i,j,k)

				enddo

				do n=1,method   !  FD method

					km = k
					kn = k-1

					nx = 0.5*(ctx(i,j,k)+ctx(i,j,kn))
					ny = 0.5*(cty(i,j,k)+cty(i,j,kn))
					nz = 0.5*(ctz(i,j,k)+ctz(i,j,kn))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,j,km,1,m)+dmudxyz(i,j,kn,1,m) 
					dmudy =  dmudxyz(i,j,km,2,m)+dmudxyz(i,j,kn,2,m)
					dmudz =  dmudxyz(i,j,km,3,m)+dmudxyz(i,j,kn,3,m)

					vis  = viseq(i,j,km,m) + viseq(i,j,kn,m)

					fv(k,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do k=cmthd,nk1
				    ftem(i,j,k,m)  = ftem(i,j,k,m) + (fv(k+1,m)-fv(k,m))
!					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(k+1,m)-fv(k,m))
				enddo
			enddo

    enddo
  enddo

  do m=1,nlamtur
  do i=cmthd,ni1
  do j=cmthd,nj1
  do k=cmthd,nk1
	if(ftem(i,j,k,m) > 0.0)then
		dqke(i,j,k,m) = dqke(i,j,k,m)-re*ftem(i,j,k,m)
	else
	  ncount = ncount+1

	  if(ncount<=3)then
	    write(*,*)'����ɢ������˷�������',I,J,K
	  endif
	endif
  enddo
  enddo
  enddo
  enddo

  IF(NCOUNT>0)WRITE(*,*)'����ɢ������˷������Ƶĵ���',ncount

  deallocate(ftem)

!------------
  return
end subroutine cal_rhs_vis_turb_dif_positive
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_rhs_viscous_SST_ncmpcor20
  use global_variables
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real :: ev_l,vis,nx,ny,nz
  real :: dmudx,dmudz,dmudy,drdx,drdy,drdz,roim,roin,rkwim,rkwin
  real :: drdxim,drdxin,drdyim,drdyin,drdzim,drdzin,tem(2)
  real,dimension( nmax+1,nlamtur ) :: fv

	cmthd = 2
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

	tem(1) = 1.0
	tem(2) = 0.5
!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1
		do i=cmthd,ni


			im = i
			in = i-1

			nx = 0.5*(kcx(i,j,k)+kcx(in,j,k))
			ny = 0.5*(kcy(i,j,k)+kcy(in,j,k))
			nz = 0.5*(kcz(i,j,k)+kcz(in,j,k))

			roim = r(im,j,k)
			roin = r(in,j,k)

			drdxim = drdxyz(im,j,k,1)
			drdxin = drdxyz(in,j,k,1)
			drdyim = drdxyz(im,j,k,2)
			drdyin = drdxyz(in,j,k,2)
			drdzim = drdxyz(im,j,k,3)
			drdzin = drdxyz(in,j,k,3)

			do m =1,nlamtur
			    
				rkwim = tem(m)*qke(im,j,k,m)/roim
				rkwin = tem(m)*qke(in,j,k,m)/roin
				
				dmudx =  dmudxyz(im,j,k,1,m) + dmudxyz(in,j,k,1,m) + &
				                rkwim*drdxim + rkwin*drdxin
				dmudy =  dmudxyz(im,j,k,2,m) + dmudxyz(in,j,k,2,m) + &
				                rkwim*drdyim + rkwin*drdyin
				dmudz =  dmudxyz(im,j,k,3,m) + dmudxyz(in,j,k,3,m) + &
				                rkwim*drdzim + rkwin*drdzin

				vis  = viseq(im,j,k,m) + viseq(in,j,k,m)

				fv(i,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
			enddo

		enddo

		do m =1,nlamtur
			do i=cmthd,ni1
				dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(i+1,m)-fv(i,m))
			enddo
		enddo
  enddo
  enddo

!------------
! J direction

	do k=cmthd,nk1
	do i=cmthd,ni1
		do j=cmthd,nj


			jm = j
			jn = j-1

			nx = 0.5*(etx(i,j,k)+etx(i,jn,k))
			ny = 0.5*(ety(i,j,k)+ety(i,jn,k))
			nz = 0.5*(etz(i,j,k)+etz(i,jn,k))

			roim = r(i,jm,k)
			roin = r(i,jn,k)

			drdxim = drdxyz(i,jm,k,1)
			drdxin = drdxyz(i,jn,k,1)
			drdyim = drdxyz(i,jm,k,2)
			drdyin = drdxyz(i,jn,k,2)
			drdzim = drdxyz(i,jm,k,3)
			drdzin = drdxyz(i,jn,k,3)

			do m =1,nlamtur

				rkwim = tem(m)*qke(i,jm,k,m)/roim
				rkwin = tem(m)*qke(i,jn,k,m)/roin

				dmudx =  dmudxyz(i,jm,k,1,m) + dmudxyz(i,jn,k,1,m) + &
				                rkwim*drdxim + rkwin*drdxin
				dmudy =  dmudxyz(i,jm,k,2,m) + dmudxyz(i,jn,k,2,m) + &
				                rkwim*drdyim + rkwin*drdyin
				dmudz =  dmudxyz(i,jm,k,3,m) + dmudxyz(i,jn,k,3,m) + &
				                rkwim*drdzim + rkwin*drdzin

				vis  = viseq(i,jm,k,m) + viseq(i,jn,k,m)

				fv(j,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25

			enddo

		enddo

		do m =1,nlamtur
			do j=cmthd,nj1
				dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(j+1,m)-fv(j,m))
			enddo
		enddo

  enddo
  enddo

!------------
! K direction
  do j=cmthd,nj1
  do i=cmthd,ni1
		do k=cmthd,nk


			km = k
			kn = k-1

			nx = 0.5*(ctx(i,j,k)+ctx(i,j,kn))
			ny = 0.5*(cty(i,j,k)+cty(i,j,kn))
			nz = 0.5*(ctz(i,j,k)+ctz(i,j,kn))

			roim = r(i,j,km)
			roin = r(i,j,kn)

			drdxim = drdxyz(i,j,km,1)
			drdxin = drdxyz(i,j,kn,1)
			drdyim = drdxyz(i,j,km,2)
			drdyin = drdxyz(i,j,kn,2)
			drdzim = drdxyz(i,j,km,3)
			drdzin = drdxyz(i,j,kn,3)

			do m =1,nlamtur

				rkwim = tem(m)*qke(i,j,km,m)/roim
				rkwin = tem(m)*qke(i,j,kn,m)/roin

				dmudx =  dmudxyz(i,j,km,1,m)+dmudxyz(i,j,kn,1,m) + &
				                rkwim*drdxim + rkwin*drdxin
				dmudy =  dmudxyz(i,j,km,2,m)+dmudxyz(i,j,kn,2,m) + &
				                rkwim*drdyim + rkwin*drdyin
				dmudz =  dmudxyz(i,j,km,3,m)+dmudxyz(i,j,kn,3,m) + &
				                rkwim*drdzim + rkwin*drdzin

				vis  = viseq(i,j,km,m) + viseq(i,j,kn,m)

				fv(k,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
			enddo

		enddo

		do m =1,nlamtur
			do k=cmthd,nk1
				dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(k+1,m)-fv(k,m))
			enddo
		enddo

  enddo
  enddo

!------------
  return
end subroutine cal_rhs_viscous_SST_ncmpcor20

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine cal_rhs_viscous_turbulence_old
  use global_variables
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real :: ev_l,vis,nx,ny,nz
	real :: dmudx,dmudz,dmudy
  real,dimension( nmax+1,nlamtur ) :: fv

	cmthd = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1
			do i=cmthd,ni

				do n=method,0   !  FV method

					im = min(i,ni1)
					in = max(1,i-1)

					nx = kcx(i,j,k)
					ny = kcy(i,j,k)
					nz = kcz(i,j,k)

				enddo

				do n=1,method   !  FD method

					im = i
					in = i-1

					nx = 0.5*(kcx(i,j,k)+kcx(in,j,k))
					ny = 0.5*(kcy(i,j,k)+kcy(in,j,k))
					nz = 0.5*(kcz(i,j,k)+kcz(in,j,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(im,j,k,1,m)+dmudxyz(in,j,k,1,m) 
					dmudy =  dmudxyz(im,j,k,2,m)+dmudxyz(in,j,k,2,m)
					dmudz =  dmudxyz(im,j,k,3,m)+dmudxyz(in,j,k,3,m)

					vis  = viseq(im,j,k,m) + viseq(in,j,k,m)

					fv(i,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do i=cmthd,ni1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(i+1,m)-fv(i,m))
				enddo
			enddo
    enddo
  enddo

!------------
! J direction

	do k=cmthd,nk1
		do i=cmthd,ni1
			do j=cmthd,nj

				do n=method,0   !  FV method

					jm = min(j,nj1)
					jn = max(1,j-1)

					nx = etx(i,j,k)
					ny = ety(i,j,k)
					nz = etz(i,j,k)

				enddo

				do n=1,method   !  FD method

					jm = j
					jn = j-1

					nx = 0.5*(etx(i,j,k)+etx(i,jn,k))
					ny = 0.5*(ety(i,j,k)+ety(i,jn,k))
					nz = 0.5*(etz(i,j,k)+etz(i,jn,k))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,jm,k,1,m)+dmudxyz(i,jn,k,1,m) 
					dmudy =  dmudxyz(i,jm,k,2,m)+dmudxyz(i,jn,k,2,m)
					dmudz =  dmudxyz(i,jm,k,3,m)+dmudxyz(i,jn,k,3,m)

					vis  = viseq(i,jm,k,m) + viseq(i,jn,k,m)

					fv(j,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25

				enddo

			enddo

			do m =1,nlamtur
				do j=cmthd,nj1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(j+1,m)-fv(j,m))
				enddo
			enddo

    enddo
  enddo

!------------
! K direction
  do j=cmthd,nj1
		do i=cmthd,ni1
			do k=cmthd,nk

				do n=method,0   !  FV method

					km = min(k,nk1)
					kn = max(1,k-1)

					nx = ctx(i,j,k)
					ny = cty(i,j,k)
					nz = ctz(i,j,k)

				enddo

				do n=1,method   !  FD method

					km = k
					kn = k-1

					nx = 0.5*(ctx(i,j,k)+ctx(i,j,kn))
					ny = 0.5*(cty(i,j,k)+cty(i,j,kn))
					nz = 0.5*(ctz(i,j,k)+ctz(i,j,kn))
				enddo

				do m =1,nlamtur

					dmudx =  dmudxyz(i,j,km,1,m)+dmudxyz(i,j,kn,1,m) 
					dmudy =  dmudxyz(i,j,km,2,m)+dmudxyz(i,j,kn,2,m)
					dmudz =  dmudxyz(i,j,km,3,m)+dmudxyz(i,j,kn,3,m)

					vis  = viseq(i,j,km,m) + viseq(i,j,kn,m)

					fv(k,m) = (dmudx*nx+dmudy*ny+dmudz*nz)*vis*0.25
				enddo

			enddo

			do m =1,nlamtur
				do k=cmthd,nk1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(k+1,m)-fv(k,m))
				enddo
			enddo

    enddo
  enddo

!------------
  return
end subroutine cal_rhs_viscous_turbulence_old
!-----------------------------------------------------------------------------!
subroutine  init_turbulence
	use global_variables,only:nstart,nameturb,ntrans,nreset
#ifdef PARALLEL
    use mod_parallels,only : calc_wall_dists_parallel,read_wall_dists_parallel, &
                             read_turbulence_parallel
#endif
	implicit none

#ifdef PARALLEL
    if (nreset == 0) then
       call calc_wall_dists_parallel
    else
       call read_wall_dists_parallel
    end if
#else
	call init_distance     ! calculating normal distance
#endif

    call modify_des_distance

	if(ntrans==1)then
		call read_trancoe    ! read transition position
	endif

	if((nstart==0) .or. (nstart==1))then
		call init_uniform_turbulence
	else
#ifdef PARALLEL
        call read_turbulence_parallel
#else
		call init_read_turbulence
#endif
	endif

	return
end subroutine  init_turbulence
!-----------------------------------------------------------------------------!
subroutine  init_uniform_turbulence
	use global_variables,only:nameturb
	implicit none

	if(nameturb=='SA') then 
	  call init_SA
	elseif(nameturb=='SST') then 
	  call init_SST
	endif

	return
end subroutine  init_uniform_turbulence
!-----------------------------------------------------------------------------!
subroutine  init_read_turbulence
	use global_variables,only:nlamtur,ni,nj,nk,nblocks,mb_qke,mb_vist,mb_dim,flowname
	implicit none

	integer :: i,j,k,nb,m

	open(1,file=trim(flowname)//".turb",status="old",form='unformatted')

	do nb=1,nblocks
	  do k=1, mb_dim(nb,3)
	    do j=1, mb_dim(nb,2)
	      do i=1, mb_dim(nb,1)    
					do m=1,nlamtur

		        read(1)mb_qke(nb)%a4d(i,j,k,m)

					enddo

					read(1)mb_vist(nb)%a3d(i,j,k)

		    enddo
		  enddo
    enddo
	enddo

  close(1)
	return
end subroutine  init_read_turbulence
!-----------------------------------------------------------------------------!
subroutine write_turbulence
	use global_variables,only:nlamtur,ni,nj,nk,nblocks,mb_qke,mb_vist,mb_dim,flowname
	implicit none

	integer :: i,j,k,nb,m

	open(1,file=trim(flowname)//".turb",status="unknown",form='unformatted')

	do nb=1,nblocks
	  do k=1, mb_dim(nb,3)
	    do j=1, mb_dim(nb,2)
	      do i=1, mb_dim(nb,1)    
					do m=1,nlamtur

		        write(1)mb_qke(nb)%a4d(i,j,k,m)

					enddo

					write(1)mb_vist(nb)%a3d(i,j,k)

		    enddo
		  enddo
    enddo
	enddo

  close(1)
	return
end subroutine write_turbulence
!-----------------------------------------------------------------------------!
subroutine init_distance
	use global_variables,only:nreset,method,nameturb
    implicit none
	integer :: m
	if(nameturb /= 'KG')then
      if( nreset == 0 ) then 
		write(*,*)'Calculating distance'
		call cal_distance
	  else
		write(*,*)'Reading distance'
		call read_distance
	  endif
	
!--------------------
	  do m=method,0
		call cal_distance_center
	  enddo
!--------------------
	endif

	return
end subroutine init_distance
!-----------------------------------------------------------------------------!
subroutine modify_des_distance
    use global_variables,only:mb_dim,ni,nj,nk,nblocks,nameturb, &
                              ndes,x,y,z,mb_distance,mb_dgrid,dgrid
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
    integer :: i,j,k,nb,pnb,m,is,ie,js,je,ks,ke
    real,parameter :: Cdes = 0.65
    real :: temdxyz,temdx,temdy,temdz,dold
    logical :: logi_ddes
    
	if(nameturb /= 'KG')then

      logi_ddes =.false.
      if(nameturb =='SA' .and. ndes==11 )then
        logi_ddes = .true.
      endif
      
      if(nameturb =='SA' .and. ndes==10 )then !DES-97
!        mb_distance = MIN(d, Cdes*Delta)
      
#ifdef PARALLEL
        do pnb=1,pnblocks
           nb = pnbindexs(pnb)
#else
        do nb=1,nblocks
#endif
           call recast_grid(nb)
           do i=1,ni
           do j=1,nj
           do k=1,nk
                is = max(1 ,i-1)
                js = max(1 ,j-1)
                ks = max(1, k-1)
      
                ie = min(ni,i+1)
                je = min(nj,j+1)
                ke = min(nk,k+1)
      
                temdx = (sqrt( (x(i,j,k)-x(is,j,k))**2 + &
                               (y(i,j,k)-y(is,j,k))**2 + &
                               (z(i,j,k)-z(is,j,k))**2 ) + &
                         sqrt( (x(i,j,k)-x(ie,j,k))**2 + &
                               (y(i,j,k)-y(ie,j,k))**2 + &
                               (z(i,j,k)-z(ie,j,k))**2 ) )/dble(ie-is)
      
                temdy = (sqrt( (x(i,j,k)-x(i,js,k))**2 + &
                               (y(i,j,k)-y(i,js,k))**2 + &
                               (z(i,j,k)-z(i,js,k))**2 ) + &
                         sqrt( (x(i,j,k)-x(i,je,k))**2 + &
                               (y(i,j,k)-y(i,je,k))**2 + &
                               (z(i,j,k)-z(i,je,k))**2 ) )/dble(je-js)
      
                temdz = (sqrt( (x(i,j,k)-x(i,j,ks))**2 + &
                               (y(i,j,k)-y(i,j,ks))**2 + &
                               (z(i,j,k)-z(i,j,ks))**2 ) + &
                         sqrt( (x(i,j,k)-x(i,j,ke))**2 + &
                               (y(i,j,k)-y(i,j,ke))**2 + &
                               (z(i,j,k)-z(i,j,ke))**2 ) )/dble(ke-ks)
      
                dold = mb_distance(nb)%a3d(i,j,k)
                mb_distance(nb)%a3d(i,j,k) = min(dold,Cdes*max(temdx,temdy,temdz))
           enddo
           enddo
           enddo
        enddo
      
      elseif( (nameturb =='SA' .or. nameturb == 'SST') .and. ndes > 0 )then !DDES OR DES-SST
!        mb_distance ���䣬��Ϊԭʼ���������
!        dgrid = Delta -->DES-SST ��
!        dgrid = MAX(0, mb_distance - Cdes*Delta) -->DDES-SA-2005
      
        allocate( mb_dgrid(nblocks) )
      
#ifdef PARALLEL
        do pnb=1,pnblocks
           nb = pnbindexs(pnb)
#else
        do nb=1,nblocks
#endif
           call recast_grid(nb)
      
           allocate( mb_dgrid(nb)%a3d(ni,nj,nk) )
      
           dgrid => mb_dgrid(nb)%a3d
           do i=1,ni
           do j=1,nj
           do k=1,nk
                is = max(1 ,i-1)
                js = max(1 ,j-1)
                ks = max(1, k-1)
      
                ie = min(ni,i+1)
                je = min(nj,j+1)
                ke = min(nk,k+1)
      
                temdx = (sqrt( (x(i,j,k)-x(is,j,k))**2 + &
                               (y(i,j,k)-y(is,j,k))**2 + &
                               (z(i,j,k)-z(is,j,k))**2 ) + &
                         sqrt( (x(i,j,k)-x(ie,j,k))**2 + &
                               (y(i,j,k)-y(ie,j,k))**2 + &
                               (z(i,j,k)-z(ie,j,k))**2 ) )/dble(ie-is)
      
                temdy = (sqrt( (x(i,j,k)-x(i,js,k))**2 + &
                               (y(i,j,k)-y(i,js,k))**2 + &
                               (z(i,j,k)-z(i,js,k))**2 ) + &
                         sqrt( (x(i,j,k)-x(i,je,k))**2 + &
                               (y(i,j,k)-y(i,je,k))**2 + &
                               (z(i,j,k)-z(i,je,k))**2 ) )/dble(je-js)
      
                temdz = (sqrt( (x(i,j,k)-x(i,j,ks))**2 + &
                               (y(i,j,k)-y(i,j,ks))**2 + &
                               (z(i,j,k)-z(i,j,ks))**2 ) + &
                         sqrt( (x(i,j,k)-x(i,j,ke))**2 + &
                               (y(i,j,k)-y(i,j,ke))**2 + &
                               (z(i,j,k)-z(i,j,ke))**2 ) )/dble(ke-ks)
      
                temdxyz = max(temdx,temdy,temdz)
      
                dgrid(i,j,k) = temdxyz
      
                if(logi_ddes)then
                  dgrid(i,j,k) = max(0.0,mb_distance(nb)%a3d(i,j,k)-Cdes*temdxyz )
                endif
           enddo
           enddo
           enddo
      
        enddo
      endif
    
    end if

    return
end subroutine modify_des_distance
!-----------------------------------------------------------------------------!
subroutine cal_distance
  use global_variables,only:mb_x,mb_y,mb_z,mb_distance,nblocks,mb_dim,nameturb
  use bl_sub_variables
  implicit none

	integer :: i,j,k,nb,m,iw,jw,kw,nbw
	real    :: xp,yp,zp,ds,xt,yt


	do nb =1,nblocks
	    if( nameturb=="BL"  )then
		  indexwall => mb_indexwall(nb)%a4d
	    endif

	    do k=1, mb_dim(nb,3)
	    do j=1, mb_dim(nb,2)
	    do i=1, mb_dim(nb,1)  
			    xp = mb_x(nb)%a3d(i,j,k)
				yp = mb_y(nb)%a3d(i,j,k)
				zp = mb_z(nb)%a3d(i,j,k)

			    call ds_compute_p(xp,yp,zp,ds,nbw,iw,jw,kw)

			    mb_distance(nb)%a3d(i,j,k) = ds

				if(nameturb=="BL") then
					indexwall(i,j,k,1) = iw
					indexwall(i,j,k,2) = jw
					indexwall(i,j,k,3) = kw
					indexwall(i,j,k,4) = nbw
				endif

	    enddo
	    enddo
	    enddo
    enddo
!--------------------
    call write_distance
!--------------------

	return
end subroutine cal_distance
!-----------------------------------------------------------------------------!
subroutine read_distance
  use global_variables,only:mb_distance,nblocks,mb_dim,nameturb,gridname
  use bl_sub_variables
  implicit none
	integer :: i,j,k,nb,m

	open(1,file=trim(gridname)//".wdst",status='old',form='unformatted')

	do nb =1,nblocks
	  do k=1, mb_dim(nb,3)
	    do j=1, mb_dim(nb,2)
	      do i=1, mb_dim(nb,1)
				 
	        read(1)mb_distance(nb)%a3d(i,j,k)

	      enddo
	    enddo
	  enddo
  enddo

  close(1)

  if( nameturb == "BL") then
	open(2,file='grid/point_to_wall.dat',status='unknown',form='unformatted')

	do nb =1,nblocks
	  do k=1, mb_dim(nb,3)
	    do j=1, mb_dim(nb,2)
	      do i=1, mb_dim(nb,1)
			do m=1,4
				 
	          read(2)mb_indexwall(nb)%a4d(i,j,k,m)

			enddo
	      enddo
	    enddo
	  enddo
    enddo

    close(2)
  endif

	return
end subroutine read_distance
!-----------------------------------------------------------------------------!
subroutine write_distance
  use global_variables,only:mb_distance,nblocks,mb_dim,nameturb,gridname
  use bl_sub_variables
  implicit none
	integer :: i,j,k,nb,m

	open(1,file=trim(gridname)//".wdst",status='unknown',form='unformatted')

	do nb =1,nblocks
	  do k=1, mb_dim(nb,3)
	    do j=1, mb_dim(nb,2)
	      do i=1, mb_dim(nb,1)
				 
	        write(1)mb_distance(nb)%a3d(i,j,k)

	      enddo
	    enddo
	  enddo
  enddo

  close(1)

  if( nameturb == "BL") then
	open(2,file='grid/point_to_wall.dat',status='unknown',form='unformatted')

	do nb =1,nblocks
	  do k=1, mb_dim(nb,3)
	    do j=1, mb_dim(nb,2)
	      do i=1, mb_dim(nb,1)
			do m=1,4
				 
	          write(2)mb_indexwall(nb)%a4d(i,j,k,m)

			enddo
	      enddo
	    enddo
	  enddo
    enddo

    close(2)

  endif

	return
end subroutine write_distance
!-----------------------------------------------------------------------------!
subroutine cal_distance_center
	use global_variables,only:mb_distance,mb_dim,nblocks
	implicit none
	real    :: ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8
	integer :: i,j,k,nb

	do nb =1,nblocks
	  do k=1,mb_dim(nb,3)-1
	    do j=1,mb_dim(nb,2)-1
	       do i=1,mb_dim(nb,1)-1

			    ds1 = mb_distance(nb)%a3d(i  ,j  ,k  )
			    ds2 = mb_distance(nb)%a3d(i+1,j  ,k  )
			    ds3 = mb_distance(nb)%a3d(i+1,j+1,k  )
			    ds4 = mb_distance(nb)%a3d(i  ,j+1,k  )
			    ds5 = mb_distance(nb)%a3d(i  ,j  ,k+1)
			    ds6 = mb_distance(nb)%a3d(i+1,j  ,k+1)
			    ds7 = mb_distance(nb)%a3d(i+1,j+1,k+1)
			    ds8 = mb_distance(nb)%a3d(i  ,j+1,k+1)

			    mb_distance(nb)%a3d(i,j,k) = 0.125 *(ds1+ds2+ds3+ds4+ &
			                                         ds5+ds6+ds7+ds8  )
			  enddo
		  enddo
	  enddo
	enddo
	return
end subroutine cal_distance_center
!-----------------------------------------------------------------------------!
subroutine ds_compute_p(xp,yp,zp,ds,nbw,iw,jw,kw)
  use global_variables,only:nblocks,mb_bc,mb_x,mb_y,mb_z,method
  implicit none
  integer :: nb,nr,bctype,m,nbw,iw,jw,kw
  integer :: s_nd,s_lr
  integer :: i,j,k,cmethod
  integer :: nrmax
  integer :: s_st(3),s_ed(3)
	real :: dxp,dyp,dzp,xp,yp,zp,ds,ds1

	cmethod = 1-method

	ds = 1.0e20

	do nb =1,nblocks
	  nrmax = mb_bc(nb)%nregions                 
	  do nr = 1,nrmax
	    bctype = mb_bc(nb)%bc(nr)%bctype
		  if(bctype == 2 .or.bctype/10==2)then
              do m=1,3
                s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)   
                s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)   
              enddo 
              s_nd  = mb_bc(nb)%bc(nr)%s_nd          
              s_lr  = mb_bc(nb)%bc(nr)%s_lr          
			  if(s_nd==1) then
			    s_ed(3) = s_ed(3) + cmethod
                s_ed(2) = s_ed(2) + cmethod
			    if(s_lr == 1)then
					  s_ed(1) = s_ed(1) + cmethod
						s_st(1) = s_ed(1)
					endif
			  elseif(s_nd==2) then
			    s_ed(3) = s_ed(3) + cmethod
                s_ed(1) = s_ed(1) + cmethod
			    if(s_lr == 1)then
					s_ed(2) = s_ed(2) + cmethod
					s_st(2) = s_ed(2)
				endif
			  else
			     s_ed(1) = s_ed(1) + cmethod
                 s_ed(2) = s_ed(2) + cmethod
			     if(s_lr == 1)then
						s_ed(3) = s_ed(3) + cmethod
						s_st(3) = s_ed(3)
				endif
			  endif 
				        
              do k = s_st(3),s_ed(3) 
              do j = s_st(2),s_ed(2) 
              do i = s_st(1),s_ed(1) 
			        dxp = xp-mb_x(nb)%a3d(i,j,k)
				      dyp = yp-mb_y(nb)%a3d(i,j,k)
				      dzp = zp-mb_z(nb)%a3d(i,j,k)
				      ds1 = sqrt(dxp*dxp + dyp*dyp + dzp*dzp) + 1.0e-20
					  if(ds > ds1)then
			              ds  = ds1
						  iw=i
						  jw=j
						  kw=k
						  nbw = nb
					  endif

	          enddo
			  enddo
			  enddo
		  endif
      enddo
    enddo
	return
end subroutine ds_compute_p
!-----------------------------------------------------------------------------!
subroutine BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k)
  use global_variables,only:mb_qke,nlamtur,method
  implicit none

	integer :: nbs,nbt,i,j,k,is,js,ks,it,jt,kt,m,n

	do m=1,nlamtur
		mb_qke(nbs)%a4d(is,js,ks,m) = mb_qke(nbt)%a4d(it,jt,kt,m)
	enddo


	return
end subroutine BC_connect_turbulence
!-----------------------------------------------------------------------------!
subroutine BC_face_dif_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k)
  use global_variables,only:mb_qke,nlamtur,method
  implicit none

	integer :: nbs,nbt,i,j,k,is,js,ks,it,jt,kt,m,n

	do m=1,nlamtur
		mb_qke(nbs)%a4d(i,j,k,m) =0.5*(mb_qke(nbt)%a4d(it,jt,kt,m) + &
		                               mb_qke(nbs)%a4d(is,js,ks,m) )
	enddo 


	return
end subroutine BC_face_dif_turbulence
!-----------------------------------------------------------------------------!
subroutine BC_wall_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k)
!*TGH. ����ģ�͸����ֵ
  use global_variables,only:mb_qke,nlamtur,method,nameturb,mb_distance, &
	                          reynolds,r,t,visc
  implicit none

	integer :: nbs,is,js,ks,it,jt,kt,m,n,i,j,k,iy1,jy1,ky1
	real    :: vsl,tm
	if(nlamtur>0)then
  
		mb_qke(nbs)%a4d(is,js,ks,1) = -mb_qke(nbs)%a4d(it,jt,kt,1)

		do m=2,nlamtur

			if(nameturb=='SST')then
        tm  = t(i,j,k)
        vsl = tm*sqrt(tm)*(1.0+visc)/(tm+visc)

				iy1 = i
				jy1 = j
				ky1 = k
				do n=1,method
					iy1 = it
					jy1 = jt
					ky1 = kt			
				enddo

				mb_qke(nbs)%a4d(is,js,ks,2) = 1600.0*vsl/r(i,j,k)/reynolds            / &
				                              mb_distance(nbs)%a3d(iy1,jy1,ky1)**2.0  - &
																			mb_qke(nbs)%a4d(it,jt,kt,2)
			elseif(nameturb=='KG')then
		    mb_qke(nbs)%a4d(is,js,ks,2) = -mb_qke(nbs)%a4d(it,jt,kt,2)
			endif

		enddo

	endif

	return
end subroutine BC_wall_turbulence
!-----------------------------------------------------------------------------!
subroutine BC_symmetry_turbulence(nb,is,js,ks,it,jt,kt,i,j,k)
  use global_variables,only:mb_qke,nlamtur,method
  implicit none

	integer :: nb,is,js,ks,it,jt,kt,m,i,j,k,n

	do m=1,nlamtur
		mb_qke(nb)%a4d(is,js,ks,m) = mb_qke(nb)%a4d(it,jt,kt,m)   
	enddo

	return
end subroutine BC_symmetry_turbulence
!-----------------------------------------------------------------------------!
subroutine BC_farfield_turbulence(nb,is,js,ks)
  use global_variables,only:reynolds,mb_qke,nlamtur,nameturb,koo,goo,muoo,omgoo
  implicit none

	integer :: nb,is,js,ks,it,jt,kt,m
	if(nlamtur>0)then
		if(nameturb=='SA')then
			mb_qke(nb)%a4d(is,js,ks,1) = muoo
		elseif(nameturb=='SST')then
			mb_qke(nb)%a4d(is,js,ks,1) = koo
			mb_qke(nb)%a4d(is,js,ks,2) = omgoo
		elseif(nameturb=='KG')then
			mb_qke(nb)%a4d(is,js,ks,1) = koo
			mb_qke(nb)%a4d(is,js,ks,2) = goo
		endif
	endif
	return
end subroutine BC_farfield_turbulence
!-----------------------------------------------------------------------------!
subroutine BC_outflow_turbulence(nb,is,js,ks,it,jt,kt)
  use global_variables,only:mb_qke,nlamtur
  implicit none

	integer :: nb,is,js,ks,it,jt,kt,m

!	if(nlamtur>0)then

	do m=1,nlamtur
		mb_qke(nb)%a4d(is,js,ks,m) = mb_qke(nb)%a4d(it,jt,kt,m)
	enddo

!	endif

	return
end subroutine BC_outflow_turbulence
!-----------------------------------------------------------------------------!
subroutine BC_peroidic_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k)
  use global_variables,only:mb_qke,nlamtur,method
  implicit none

	integer :: nb,is,js,ks,it,jt,kt,m,i,j,k,n,nbs,nbt

	do m=1,nlamtur
		mb_qke(nbs)%a4d(is,js,ks,m) = mb_qke(nbt)%a4d(it,jt,kt,m)
	enddo

	return
end subroutine BC_peroidic_turbulence
!-----------------------------------------------------------------------------!
subroutine BC_jet_turbulence(nb,is,js,ks,it,jt,kt)
  use global_variables,only:mb_qke,nlamtur,nameturb
  implicit none

	integer :: nb,is,js,ks,it,jt,kt,m
	real:: parameter_wall(nlamtur),parameter_inflow(nlamtur)
	real:: cv13,kafan,kf3,fv1,vt

	
	if(nlamtur>0)then
		if(nameturb=='SA')then
		endif

		mb_qke(nb)%a4d(is,js,ks,1) = parameter_inflow(1)

		do m=2,nlamtur
			mb_qke(nb)%a4d(is,js,ks,m) = 100.0
		enddo

	endif

	return
end subroutine BC_jet_turbulence
!-----------------------------------------------------------------------------!
! intruduce transition model for flat flow 
subroutine transition_turbulence !(nb)
	use global_variables
	implicit none
	integer :: i,j,k,cmethod,nb

	cmethod = method - 1
  nb = 1
	if(ntrans==2)then
		if(mb_istran(nb) == 1)then
			do k=1,nk+cmethod
				do j=1,nj+cmethod
				  do i=1,ni+cmethod
						vist(i,j,k) = vist(i,j,k)*mb_trancoe(nb)%a3D(i,j,k)
					enddo
				enddo
			enddo
		endif
	else
		do k=1,nk+cmethod
			do j=1,nj+cmethod
			  do i=1,ni+cmethod
					if(x(i,j,k)<=xtrans)then
						vist(i,j,k) = amin1(vtoo,vist(i,j,k))
					endif
				enddo
			enddo
		enddo
	endif

	return
end subroutine transition_turbulence
!-----------------------------------------------------------------------------!
subroutine max_vist
	use global_variables
	implicit none
	integer :: i,j,k,imax,jmax,kmax,cmethod,ni1,nj1,nk1
	real :: vistmax

	cmethod = method-1
	imax = -1; jmax = -1; kmax = -1
	vistmax = -1.0e20

	ni1 = ni+cmethod
	nj1 = nj+cmethod
	nk1 = nk+cmethod

	do k=1,nk1
		do j=1,nj1
			do i=1,ni1
				if(vistmax < vist(i,j,k))then
					vistmax = vist(i,j,k)
					imax=i
					jmax=j
					kmax=k
				endif
			enddo
		enddo
	enddo
	write(*,*)vistmax,imax,jmax,kmax
	return
end subroutine max_vist
!_____________________________________________________________________!
subroutine  read_trancoe
	use global_variables,only:mb_trancoe,nblocks,mb_istran
	implicit none
	integer :: n_tr,nb_tr,nbb,imax,jmax,kmax,nbq,i,j,k
  logical :: file_exist

	write(*,*)
  write(*,*)' ===================================================================='
	write(*,*)
  write(*,*)'      提示：如果当前计算问题考虑了流动转捩 '
  write(*,*)'            请按帮助文档中的使用手册       '
  write(*,*)'     “4、部分接口说明 （7）转捩描述文件准备文件 trancoe.dat”'
  write(*,*)
  write(*,*)' ===================================================================='
	write(*,*)
!
! check trancoe.dat file exist or not
  inquire( file = 'trancoe.dat' , exist = file_exist )
  if(.not. file_exist) then
		write(*,*)
    write(*,*)' ===================================================================='
		write(*,*)
    write(*,*)'      当前目录不存在描述转捩的数据文件 trancoe.dat '
    write(*,*)'      如果考虑了流动转捩，请按帮助文档中的使用手册'
    write(*,*)'     “4、部分接口说明 （7）转捩描述文件准备文件 trancoe.dat”'
    write(*,*)
    write(*,*)' ===================================================================='
		write(*,*)
  else
	  open(1,file='trancoe.dat')
	  read(1,*)n_tr
	  allocate(mb_istran (nblocks))
	  mb_istran (1:nblocks) = 0

	  if(n_tr==1)then
	 	  allocate(mb_trancoe(nblocks))
		  read(1,*)nb_tr
		  do nbb = 1,nb_tr
			  read(1,*)nbq,imax,jmax,kmax
			  allocate(mb_trancoe(nbq)%a3d(imax,jmax,kmax))
			  mb_istran(nbq) = 1
			  read(1,*)(((mb_trancoe(nbq)%a3d(i,j,k),k=1,kmax),j=1,jmax),i=1,imax)
		  enddo
	  endif
	  close(1)
   endif
	 return
end subroutine read_trancoe
!-----------------------------------------------------------------------------!
subroutine cal_rhs_source_SA
	use global_variables
	implicit none
	real,parameter :: cb1=0.1355,cb2=0.622,cv1=7.1,eps=1.0e-20 ,    &
	                  ct1=1.0   ,ct2=2.0  ,ct3=0.3,ct4=2.0     ,    &
										cga=2./3.0,ki =0.419,cw2=0.3,cw3=2.0     ,    &
										cw1=cb1/ki/ki+(1.0+cb2)/cga,cv13=cv1**3.0,    &
										ki2=ki*ki,cw36=cw3**6.0,onesix=1.0/6.0   ,    &
										c5=3.5
! �ر�ע�⣬����ģ��ϵ��ԭʼ����Ϊ��Ct3=1.1
! Ct3 ֻ�ڿ���ת�������ʹ�ã����� Ft2 ������
	integer :: i,j,k,cmethod,n,ni1,nj1,nk1
	real    :: kafan,fv1,fv2,kf3,fw,ft1,ft2,gbar6,rbar5
	real    :: sbar,rbar,gbar
	real    :: dn,dn2,omg,reki2dn2
	real    :: Sp,Sd,Dbar,ro,vsl,vtbar,romu,volume
	real    :: dmudx,dmudy,dmudz,cdmu,cpmu,usound,Sdc
	real    :: drdx,drdy,drdz,ddmudi,ddmudj,ddmudk,ddmudx,ddmudy,ddmudz
	integer :: iba,ifr,jba,jfr,kba,kfr
	real    :: dfv2dkafan,dsbardmu,drbardmu,dgbardrbar,dfwdg,dfwdmu
	logical :: logi_sound,logi_cmp0,logi_cmp10,logi_cmp20,logi_cmp30,logi_ddes
	real    :: ux,uy,uz,vx,vy,vz,wx,xy,wz,uijuij,dn_ddes,fd_ddes   ! DDES �е�fd����

	n       =   1
	cmethod =   method
	ft1     =   0.0
	ft2     =   0.0
	logi_sound = .false.
	logi_cmp0  = .false.
	logi_cmp10 = .false.
	logi_cmp20 = .false.
	logi_cmp30 = .false.
	logi_ddes  = .false.

	if( (ncmpcor==1) .or. (ncmpcor==21) .or. (ncmpcor==31) )then !�������ٽ��п�ѹ������
	  logi_sound = .true.
	endif
	if( (ncmpcor <= 9) )then !
	  logi_cmp0 = .true.
	endif
	if( (ncmpcor >=10) .and. (ncmpcor<=15) )then !���ֿ�ѹ��������ʽ10 ~ 15
	  logi_cmp10 = .true.
	endif
	if( (ncmpcor >=20) .and. (ncmpcor<=29) )then !��ѹ��������ʽ20 ~ 29
	  logi_cmp20 = .true.
	endif

 	if( (ncmpcor >=30) .and. (ncmpcor<=39) )then !��ѹ��������ʽ30 ~ 39
	  logi_cmp30 = .true.
	endif

	if(ndes == 11)then
	  logi_ddes= .true.
	endif

	call dfidxyz_center(-2,3,ni,nj,nk,u  ,dudxyz ,method,1)
	call dfidxyz_center(-2,3,ni,nj,nk,v  ,dvdxyz ,method,1)
	call dfidxyz_center(-2,3,ni,nj,nk,w  ,dwdxyz ,method,1)
	!!call dfidxyz_center(-1,1,ni,nj,nk,qke,dmudxyz,method,n)
	call dfidxyz_center(-2,3,ni,nj,nk,qke,dmudxyz,method,n)

	if( ncmpcor >=10)then
      call dfidxyz_center(-2,3,ni,nj,nk,r,drdxyz ,method,1)
	endif

	ni1 = ni
	nj1 = nj
	nk1 = nk

!!	do k = cmethod,nk1
!!	do j = cmethod,nj1
!!	do i = cmethod,ni1
	do k = 2,nk-1
	do j = 2,nj-1
	do i = 2,ni-1

		ro      =    r(i,j,k)
		volume  =    vol(i,j,k)
		vsl     =	 visl(i,j,k)
		vtbar   =    qke(i,j,k,n)
		romu    =    ro*vtbar

		dn      =    max(ds_turbulence(i,j,k),re)
		dn2     =    dn*dn

		reki2dn2=    reynolds*ki2*dn2

		if(logi_ddes .or. logi_sound)then
		  ux     =  dudxyz(i,j,k,1)
!		  uy     =  dudxyz(i,j,k,2)
!		  uz     =  dudxyz(i,j,k,1)
!		  vx     =  dvdxyz(i,j,k,1)
		  vy     =  dvdxyz(i,j,k,2)
!		  vz     =  dvdxyz(i,j,k,3)
!		  wx     =  dwdxyz(i,j,k,1)
!		  wy     =  dwdxyz(i,j,k,2)
		  wz     =  dwdxyz(i,j,k,3)
		  uijuij = (ux+vy+wz)*(ux+vy+wz)
		endif
		
		if(logi_ddes)then !DDES ����
			sd = ( vsl + vist(i,j,k) )/reki2dn2/sqrt(uijuij)
			fd_ddes = 1.- tanh(512.0*sd*sd*sd)
			
			dn_ddes = dn - fd_ddes*DGRID(I,J,K)

			reki2dn2=    reynolds*ki2*dn_ddes*dn_ddes
		endif

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
		drdx    =    drdxyz(i,j,k,1)
		drdy    =    drdxyz(i,j,k,2)
		drdz    =    drdxyz(i,j,k,3)

		Dbar    =    cb2/cga*(dmudx**2.0+dmudy**2.0+dmudz**2.0)*re

!		dfv2dkafan    = (3.0*cv13*(fv1/kafan)**2.0-1.0)/(1.0+kafan*fv1)**2.0

!		dsbardmu      = (fv2+kafan*dfv2dkafan)/reki2dn2 

!		drbardmu			= (rbar-rbar/sbar*dsbardmu)/ro

!		dgbardrbar		=  1.0+cw2*(6.0*rbar5-1.0)

!		dfwdg         =  fw/gbar*(1.0-gbar6/(gbar6+cw36))

!		dfwdmu        =  dfwdg*dgbardrbar*drbardmu
				
		dqke(i,j,k,n) = dqke(i,j,k,n)-ro*(Sp-Sd+Dbar)*volume

		spec(i,j,k,n) = spec(i,j,k,n)-amin1(cpmu-2.0*cdmu,0.0)*volume

		if(logi_sound)then !�������ٽ��п�ѹ������
			usound =  c(i,j,k)
			Sdc    =  re*c5*ro*uijuij*(vtbar/usound)**2
!			dqke(i,j,k,n) = dqke(i,j,k,n)-Sdc*volume
			dqke(i,j,k,n) = dqke(i,j,k,n)+Sdc*volume
		endif

		if(logi_cmp20 .or. logi_cmp30 )then
		   ! ��ѹ����
		   ! �ο�����"SA����ģ���Լ������������-��ʱ����.doc"
			sdc =   cb2 * vtbar * ( drdx*dmudx + drdy*dmudy + drdz*dmudz )

		    if(logi_cmp30)then
			  sdc = sdc + cb2 * vtbar *vtbar*(drdx*drdx+drdy*drdy+drdz*drdz) *0.25/ro
		    endif

			dqke(i,j,k,n) = dqke(i,j,k,n)-Sdc/cga*re*volume
		endif

		! ����.��ѹ����
	enddo
	enddo
	enddo


	if( logi_cmp0 )then
	  call viscosity_equivalent_SA_r !���ӳ���� viscosity_equivalent_SA��һ���ܶ���
	else
	  call viscosity_equivalent_SA
	endif

	return
end subroutine cal_rhs_source_SA
!-----------------------------------------------------------------------------!
subroutine init_SA
	use global_variables,only:ni,nj,nk,nblocks,mb_qke,mb_vist,mb_dim,muoo
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
	real,parameter :: cv13  = 7.1**3.0
	integer::i,j,k,nb,pnb
	real:: kafan,kf3,fv1,vt


	kafan = muoo
	kf3   = kafan**3.0
	fv1   = kf3/(kf3+cv13)
	vt    = muoo*fv1

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
	  do k=-1, mb_dim(nb,3)+1
	    do j=-1, mb_dim(nb,2)+1
		    do i=-1, mb_dim(nb,1)+1
			    mb_qke (nb)%a4d(i,j,k,1) = muoo
			    mb_vist(nb)%a3d(i,j,k  ) = vt
				enddo
			enddo
		enddo
	enddo
	return
end subroutine init_SA
!-----------------------------------------------------------------------------!
subroutine cal_vist_SA
	use global_variables,only:method,ni,nj,nk,qke,r,visl,vist
	implicit none
	real,parameter :: cv13  = 7.1**3.0
	integer :: i,j,k,nb,cmethod
	real :: fv1,kafan,kf3,romu
  
	cmethod = method - 1

	do k = 1,nk+cmethod
	  do j = 1,nj+cmethod
	    do i = 1,ni+cmethod

				 romu  = r(i,j,k)*qke(i,j,k,1)

		     kafan = romu/visl(i,j,k)
				 kf3   = kafan**3.0
		     fv1   = kf3/(kf3+cv13)
		     vist(i,j,k) = romu*fv1 

	     enddo
	   enddo
	enddo

	return
end subroutine cal_vist_SA
!-----------------------------------------------------------------------------!
subroutine viscosity_equivalent_SA
	use global_variables,only:method,ni,nj,nk,qke,r,visl,vist,viseq
	implicit none
	real,parameter :: ocga=1.5
	integer :: i,j,k,nb,cmethod,ni1,nj1,nk1
  
	cmethod = method - 1

	ni1 = ni+cmethod
	nj1 = nj+cmethod
	nk1 = nk+cmethod

	do k = 1,nk1
	  do j = 1,nj1
	    do i = 1,ni1
				viseq(i,j,k,1) = ocga*(visl(i,j,k)+qke(i,j,k,1)*r(i,j,k))
	    enddo
	  enddo
	enddo

	return
end subroutine viscosity_equivalent_SA
!-----------------------------------------------------------------------------!
subroutine viscosity_equivalent_SA_r
! viscosity_equivalent_SA_r �� viscosity_equivalent_SA�Ĳ������
! ǰ��û���ܶȣ��������ܶ�
	use global_variables,only:method,ni,nj,nk,qke,r,visl,vist,viseq
	implicit none
	real,parameter :: ocga=1.5
	integer :: i,j,k,nb,cmethod,ni1,nj1,nk1
  
	cmethod = method - 1

	ni1 = ni+cmethod
	nj1 = nj+cmethod
	nk1 = nk+cmethod

	do k = 1,nk1
	  do j = 1,nj1
	    do i = 1,ni1
				viseq(i,j,k,1) = ocga*(visl(i,j,k)/r(i,j,k)+qke(i,j,k,1))
	    enddo
	  enddo
	enddo

	return
end subroutine viscosity_equivalent_SA_r
!-----------------------------------------------------------------------------!
subroutine cal_rhs_source_SST
	use global_variables
	implicit none
	real,parameter :: cgak(2)=[0.85,1.0],cgaw(2)=[0.5,0.856],beta(2)=[0.075,0.0828]
	real,parameter :: CDmin=1.0e-20,a1=0.31,ki=0.41,ki2=ki*ki
	real,parameter :: twothd=0.2/0.3
	real,parameter :: MT02 = 0.25*0.25,alpha1=1.0,alpha2=0.4,alpha3=0.2
	real,parameter :: Cdes1 = 0.78, Cdes2 = 0.61

	integer :: i,j,k,cmethod,n,ni1,nj1,nk1,m
	real    :: ro,dn,dn2,volume,vsl,vst,qk,qw,omg,omg2,re500,cgaw24,cgaw22
	real    :: dkxdw,rowdkdw,CDkw,sqrtkwdn,vlrodn2w,rokdn2
	real    :: Sp(2),Sd(2),Dbar(2),cga(2),F2,arg1,arg2,beita,ga
	real    :: ux,uy,uz,vx,vy,vz,wx,wy,wz,rok,row,cdk,cpk,cdw,cpw,dbw
	real    :: Mat2,usound,FMt,Sdc(2),Spc(2),pdbar,div_u
	real,pointer,dimension(:,:,:) :: F1
	real    :: sck,dck,dcw,ma_turb ! ��ѹ������������
	real    :: dilatation_p !ѹ������ЧӦ
	real    :: tem_kw,tem_kr
	logical :: logi_drdxyz,logi_dila,logi_ncmp20,logi_mt,logi_des
	real    :: rox,roy,roz !Ϊ�����Zemanѹ����������
	real    :: beita0,betax0,betax,obetx,ki2osqrtbetx
	real    :: l_kw,l_kw_des,Cdes_kw
	real    :: x_tem

	betax0 = 0.09
	betax  = betax0
	obetx=1.0/betax
	ki2osqrtbetx=ki2/SQRT(betax)

	allocate(F1(ni,nj,nk))
	
	if(ndes > 0)then
	  logi_des = .true.
	else
	  logi_des = .false.
	endif

	if(ncmpcor==1 .or. ncmpcor-10==1 .or. ncmpcor-20==1)then ! ������������������Beta
	  logi_mt = .true.
	else
	  logi_mt = .false.
	endif
	
	if(ncmpcor >=20)then
	  logi_drdxyz = .true.
	else
	  logi_drdxyz = .false.
	endif

	if(ncmpcor >=10 .and. ncmpcor<=19 )then
	  logi_dila = .true.
	else
      logi_dila = .false.
	endif

	if(ncmpcor >=20 .and. ncmpcor<=29 )then
	  logi_ncmp20 = .true.
	else
      logi_ncmp20 = .false.
	endif

	call dfidxyz_center(-2,3,ni,nj,nk,u  ,dudxyz ,method,1)
	call dfidxyz_center(-2,3,ni,nj,nk,v  ,dvdxyz ,method,1)
	call dfidxyz_center(-2,3,ni,nj,nk,w  ,dwdxyz ,method,1)
	!!call dfidxyz_center(-1,1,ni,nj,nk,qke,dmudxyz,method,2)
	call dfidxyz_center(-2,3,ni,nj,nk,qke,dmudxyz,method,2)

	if( logi_drdxyz )then ! �����ܶ��ݶ�
		call dfidxyz_center(-2,3,ni,nj,nk,r,drdxyz ,method,1)
	endif

	cmethod =   method - 1
	re500  = re*500.0
	cgaw24 = cgaw(2)*4.0
	cgaw22 = cgaw(2)*2.0

	do k = 1,nk
	do j = 1,nj
	do i = 1,ni

		ro       =   r(i,j,k)
		qk       =   qke(i,j,k,1)
		qw       =   qke(i,j,k,2)
		vsl      =	 visl(i,j,k)
		vst      =	 vist(i,j,k)

		rok      =   ro*qk
		row      =   ro*qw

		dn       =   ds_turbulence(i,j,k)
		dn2      =   dn*dn
		dkxdw    =   dmudxyz(i,j,k,1,1)*dmudxyz(i,j,k,1,2) + &
					 dmudxyz(i,j,k,2,1)*dmudxyz(i,j,k,2,2) + &
					 dmudxyz(i,j,k,3,1)*dmudxyz(i,j,k,3,2)
		rowdkdw  =    ro/qw*dkxdw*cgaw22 
	
		CDkw     =    amax1(CDmin,rowdkdw)
				
		sqrtkwdn =    obetx*sqrt(qk)/qw/dn
		vlrodn2w =    re500*vsl/row/dn2 
		rokdn2   =    cgaw24*rok/CDkw/dn2

		arg2     =    amax1(sqrtkwdn,vlrodn2w)
		arg1     =    amin1(arg2,rokdn2)

		if(arg1>=2.0)then
			F1(i,j,k)= 1.0
		else   
			F1(i,j,k)= tanh(arg1**4.0)
		endif

		F2       =    1.0-F1(i,j,k)

		cga(1)   =    F1(i,j,k)*cgak(1)+F2*cgak(2)
		cga(2)   =    F1(i,j,k)*cgaw(1)+F2*cgaw(2)

		do m=1,2
			viseq(i,j,k,m)=vsl+cga(m)*vst
		enddo

	enddo
	enddo
	enddo

	cmethod =   method + 1
	ni1 = ni-1
	nj1 = nj-1
	nk1 = nk-1

	do k = cmethod,nk1
	do j = cmethod,nj1
	do i = cmethod,ni1

		ro       =   r(i,j,k)
		qk       =   qke(i,j,k,1)
		qw       =   qke(i,j,k,2)
		volume   =   vol(i,j,k)
		vst      =	 vist(i,j,k)

		rok      =   ro*qk
		row      =   ro*qw

		F2      =    1.0-F1(i,j,k)
		beita   =    F1(i,j,k)*beta(1)+F2*beta(2)
		cga(2)  =    F1(i,j,k)*cgaw(1)+F2*cgaw(2)

		if(logi_mt)then ! ������������������Beta
		  usound =  c(i,j,k)
		  ma_turb = 2.*qk/usound/usound
		  Fmt = 1.5*max(ma_turb-mt02,0.0)
		  betax = betax0 * (1. + Fmt)
		  beita = beita - betax0*Fmt
		  obetx = 1.0/betax
		  ki2osqrtbetx = ki2/SQRT(betax)
		endif
		
	    ga      =    beita*obetx-cga(2)*ki2osqrtbetx

		ux      =    dudxyz(i,j,k,1)
		uy      =    dudxyz(i,j,k,2)
		uz      =    dudxyz(i,j,k,3)
		vx      =    dvdxyz(i,j,k,1)
		vy      =    dvdxyz(i,j,k,2)
		vz      =    dvdxyz(i,j,k,3)
		wx      =    dwdxyz(i,j,k,1)
		wy      =    dwdxyz(i,j,k,2)
		wz      =    dwdxyz(i,j,k,3)
		div_u   = ux + vy + wz

		dkxdw    =   dmudxyz(i,j,k,1,1)*dmudxyz(i,j,k,1,2) + &
					 dmudxyz(i,j,k,2,1)*dmudxyz(i,j,k,2,2) + &
					 dmudxyz(i,j,k,3,1)*dmudxyz(i,j,k,3,2)

		omg     =    vst*((2.0*(ux*ux+vy*vy+wz*wz+uy*vx+uz*wx+wy*vz) + &
			                    uy*uy+uz*uz+vx*vx+vz*vz+wx*wx+wy*wy) - &
								 twothd*div_u*div_u )*re             - &
									 twothd*div_u*rok

		cdk     =    betax*qw
		cdw     =    beita*qw

		Sd(1)   =    cdk*rok
		Sd(2)   =    cdw*row

!!!x_tem=x(i,j,k)
!!!if(x_tem<3.0)then
!!!	Sp(1)   =    amin1(omg,0.1*Sd(1))  !�ر�ע�⣬Ϊ������פ���������ɽ������޸�
!!!elseif(x_tem >=3.0 .and. x_tem<23.0)then
!!!    sp(1) = min(omg,(19.9*(x_tem-3.0)/20.0 + 0.1)*sd(1))
!!!else
!!!	Sp(1)   =    amin1(omg,20.0*Sd(1))  !�ر�ע�⣬Ϊ������פ���������ɽ������޸�
!!!endif

		Sp(1)   =    amin1(omg,20.0*Sd(1))
		Sp(2)   =    ga*ro*omg/vst*reynolds

		dbw     =    cgaw22*F2*dkxdw/qw/qw
		Dbar(1) =    0.0
		Dbar(2) =    dbw*row

!		spec(i,j,k,1) = spec(i,j,k,1)+(cdk+twothd*amax1(ux+vy+wz,0.0))*volume
!		spec(i,j,k,2) = spec(i,j,k,2)+(2.0*cdw+abs(dbw))*volume
!       �����ǲ���CFL3D���װ뾶�������ǼӴ��װ뾶
		spec(i,j,k,1) = spec(i,j,k,1)+2.0*(cdk+twothd*amax1(ux+vy+wz,0.0))*volume
		spec(i,j,k,2) = spec(i,j,k,2)+2.0*(2.0*cdw+abs(dbw))*volume

		if(logi_dila)then  !ѹ����������
		! ���ο�SST����������ģ���Լ������������-��ʱ����.doc"
		  usound =  c(i,j,k)
		  ma_turb = 2.*qk/usound/usound
		  ma_turb = max(ma_turb,0.0)

		  if(.true.)then !Sarkar��������
		    dilatation_p = (-alpha2*Sp(1) + alpha3*betax*row*qk)*ma_turb
!		    dilatation_p = (-alpha2*omg   + alpha3*betax*row*qk)*ma_turb
		  else !El Baz & Launder ����������
		    dilatation_p = 0.75 * ma_turb * ( 2.666667*rok*(ux+vy+wz) - Sp(1) )
		  endif

		  Sck = 1.0 + alpha1*ma_turb*f2
		  Dck = f2 * dilatation_p
		  Dcw = f2 * ( betax*alpha1*ma_turb*row*qw - reynolds*ro*dilatation_p/vst )

		  Sd(1) = Sd(1) * Sck
		  Dbar(1) = Dbar(1) + Dck
		  Dbar(2) = Dbar(2) + Dcw
		endif

		if(logi_ncmp20)then !�ܶ��ݶ�����
		! ���ο�SST����������ģ���Լ������������-��ʱ����.doc"
		  tem_kw = 2.0*qk/qw
		  tem_kr = qk/ro

		  Dcw = 0.0
		  do m=1,3
			 Dcw = Dcw + drdxyz(i,j,k,m)*(dmudxyz(i,j,k,m,1)+tem_kw*dmudxyz(i,j,k,m,2)+tem_kr*drdxyz(i,j,k,m))
		  enddo
		  Dcw = f2*cgaw(2)*Dcw

		  Dbar(2) = Dbar(2) + Dcw
		ENDIF

		if(logi_des)then
			Cdes_kw = (1.-f2)*Cdes1 + f2*Cdes2
			l_kw = sqrt(qk)/cdk
			l_kw_des = min( l_kw,Cdes_kw*dgrid(i,j,k) )
			
			Sd(1) = Sd(1)*l_kw/l_kw_des
		endif

		do m=1,2
			dqke(i,j,k,m) = dqke(i,j,k,m)-(Sp(m)-Sd(m)+Dbar(m))*volume
		enddo
				
!		if(logi_mt)then !�����е�����
		if(.false.)then !�����е�����
			usound =  c(i,j,k)
			Mat2   = (2.0*qk/usound**2)*F2
			FMt    =  amax1(Mat2-MT02,0.0)
			Sdc(1) =  alpha1*FMt*Sd(1)
			Sdc(2) = -alpha1*FMt*Sd(2)*betax/beita

			pdbar  =  FMt*(alpha3*Sd(1)-alpha2*Sp(1))

			Spc(1) =  pdbar
			Spc(2) = -pdbar*ro/vsl*reynolds

			do m=1,2
				dqke(i,j,k,m) = dqke(i,j,k,m)-(Spc(m)-Sdc(m))*volume
			enddo
		endif
	enddo
	enddo
	enddo

	deallocate(F1)

	return
end subroutine cal_rhs_source_SST
!-----------------------------------------------------------------------------!
subroutine init_SST
	use global_variables,only:ni,nj,nk,nblocks,mb_qke,mb_vist,mb_dim,reynolds, &
	                          koo,omgoo,vtoo
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
	integer::i,j,k,nb,pnb

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
	  do k=-1, mb_dim(nb,3)+1
	    do j=-1, mb_dim(nb,2)+1
		    do i=-1, mb_dim(nb,1)+1   
			    mb_qke (nb)%a4d(i,j,k,1) = koo
			    mb_qke (nb)%a4d(i,j,k,2) = omgoo
			    mb_vist(nb)%a3d(i,j,k  ) = vtoo
				enddo
			enddo
		enddo
	enddo
	return
end subroutine init_SST
!-----------------------------------------------------------------------------!
subroutine cal_vist_SST
	use global_variables,only:method,ni,nj,nk,qke,r,visl,vist,re,ds_turbulence, &
	                          reynolds,dudxyz,dvdxyz,dwdxyz
	implicit none
	real,parameter :: a1=0.31,obetx=2.0/0.09 !
	integer :: cmethod,i,j,k
	real    :: ro,omg,qk,qw,sqrtkwdn,vlrodn2w,vsl,b1,dn,dn2,arg2,F2,a1reynolds
  
	cmethod = method - 1

	a1reynolds = a1*reynolds
	b1=500.0*re

	do k = 1,nk+cmethod
	  do j = 1,nj+cmethod
	    do i = 1,ni+cmethod

				ro      =    r(i,j,k)
				qk      =    qke(i,j,k,1)
				qw      =    qke(i,j,k,2)
				vsl     =    visl(i,j,k)
				dn      =    ds_turbulence(i,j,k)
				dn2     =    dn*dn

				sqrtkwdn=    obetx*sqrt(qk)/qw/dn
				vlrodn2w=    b1*vsl/ro/dn2/qw 

				arg2    =    amax1(sqrtkwdn,vlrodn2w)

				if(arg2>=4.0)then
					F2    =    1.0
				else
					F2    =    tanh(arg2*arg2)
				endif

				omg     =    sqrt((dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                  (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
													(dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

		    vist(i,j,k) = a1reynolds*ro*qk/amax1(a1*qw,omg*F2)
	    enddo
	  enddo
	enddo

	return
end subroutine cal_vist_SST
!_____________________________________________________________________!
subroutine inflow_turbulence
  use global_variables,only:koo,omgoo,goo,muoo,vtoo,reynolds, &
	                     nameturb,lfref,lref
	implicit none
	real :: lref_lfref

	lref_lfref = lref/lfref
	vtoo       = 1.0e-6  ! �ǲ��Ǹ�ȡ" vtoo = 1.e-3 "?
	muoo       = 1.0e-1  !SAģ�͵ĳ�ֵ

	if(nameturb=='SST')then
		omgoo   = 10.0*lref_lfref      !SSTģ�͵ĳ�ֵ
		koo     = omgoo*vtoo/reynolds  !SSTģ�͵ĳ�ֵ
	endif
	return	
end subroutine inflow_turbulence
!_____________________________________________________________________!
!-----------------------------------------------------------------------------!
subroutine viscosity_equiv_lamturb
	use global_variables,only:method,ni,nj,nk,qke,r,visl,vist,viseq
	implicit none
	integer :: i,j,k,nb,cmethod,ni1,nj1,nk1
	real :: ocga

	ocga = 0.75
	do k = 1,nk
	do j = 1,nj
	do i = 1,ni
		viseq(i,j,k,1) = ocga * ( visl(i,j,k)/r(i,j,k) + qke(i,j,k,1) ) * qke(i,j,k,1)
	enddo
	enddo
	enddo

	return
end subroutine viscosity_equiv_lamturb
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_rhs_viscous_ncmpcor30
  use global_variables
  implicit none
  integer :: i,j,k,m,n,cmthd,im,in,jm,jn,km,kn,ni1,nj1,nk1
  real :: ev_l,vis,nx,ny,nz
	real :: drdx,drdz,drdy
  real,dimension( nmax+1,nlamtur ) :: fv

	cmthd = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

!------------
! I direction
	do k=cmthd,nk1
    do j=cmthd,nj1
		do i=cmthd,ni

			im = i
			in = i-1

			nx = 0.5*(kcx(i,j,k)+kcx(in,j,k))
			ny = 0.5*(kcy(i,j,k)+kcy(in,j,k))
			nz = 0.5*(kcz(i,j,k)+kcz(in,j,k))

			do m =1,nlamtur
				drdx =  drdxyz(im,j,k,1)+drdxyz(in,j,k,1) 
				drdy =  drdxyz(im,j,k,2)+drdxyz(in,j,k,2)
				drdz =  drdxyz(im,j,k,3)+drdxyz(in,j,k,3)

				vis  = viseq(im,j,k,m) + viseq(in,j,k,m)

				fv(i,m) = (drdx*nx+drdy*ny+drdz*nz)*vis*0.25
			enddo

			enddo

			do m =1,nlamtur
				do i=cmthd,ni1
					dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(i+1,m)-fv(i,m))
				enddo
			enddo
    enddo
  enddo

!------------
! J direction

	do k=cmthd,nk1
	do i=cmthd,ni1
		do j=cmthd,nj

			jm = j
			jn = j-1

			nx = 0.5*(etx(i,j,k)+etx(i,jn,k))
			ny = 0.5*(ety(i,j,k)+ety(i,jn,k))
			nz = 0.5*(etz(i,j,k)+etz(i,jn,k))

			do m =1,nlamtur

				drdx =  drdxyz(i,jm,k,1)+drdxyz(i,jn,k,1) 
				drdy =  drdxyz(i,jm,k,2)+drdxyz(i,jn,k,2)
				drdz =  drdxyz(i,jm,k,3)+drdxyz(i,jn,k,3)
				vis  = viseq(i,jm,k,m) + viseq(i,jn,k,m)

				fv(j,m) = (drdx*nx+drdy*ny+drdz*nz)*vis*0.25

			enddo

		enddo

		do m =1,nlamtur
			do j=cmthd,nj1
				dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(j+1,m)-fv(j,m))
			enddo
		enddo
  enddo
  enddo

!------------
! K direction
  do j=cmthd,nj1
  do i=cmthd,ni1
		do k=cmthd,nk

			km = k
			kn = k-1

			nx = 0.5*(ctx(i,j,k)+ctx(i,j,kn))
			ny = 0.5*(cty(i,j,k)+cty(i,j,kn))
			nz = 0.5*(ctz(i,j,k)+ctz(i,j,kn))

			do m =1,nlamtur

				drdx =  drdxyz(i,j,km,1)+drdxyz(i,j,kn,1) 
				drdy =  drdxyz(i,j,km,2)+drdxyz(i,j,kn,2)
				drdz =  drdxyz(i,j,km,3)+drdxyz(i,j,kn,3)

				vis  = viseq(i,j,km,m) + viseq(i,j,kn,m)
				fv(k,m) = (drdx*nx+drdy*ny+drdz*nz)*vis*0.25
			enddo

		enddo

		do m =1,nlamtur
			do k=cmthd,nk1
				dqke(i,j,k,m) = dqke(i,j,k,m)-re*(fv(k+1,m)-fv(k,m))
			enddo
		enddo

  enddo
  enddo

!------------
  return
end subroutine cal_rhs_viscous_ncmpcor30
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine cal_spec_turbulence_sst_tem
	use global_variables
	implicit none

	integer :: i,j,k,m,cmethod,ni1,nj1,nk1
	real    :: unx,uny,unz,vis,nx,ny,nz,re2

	cmethod = method+1
	ni1     = ni-1
	nj1     = nj-1
	nk1     = nk-1

	re2     = re*2.0

	do m=1,nlamtur
		do k=cmethod,nk1
			do j=cmethod,nj1
				do i=cmethod,ni1

		      unx = u(i,j,k)*kcx(i,j,k)+v(i,j,k)*kcy(i,j,k)+w(i,j,k)*kcz(i,j,k)
			  uny = u(i,j,k)*etx(i,j,k)+v(i,j,k)*ety(i,j,k)+w(i,j,k)*etz(i,j,k)
		  	  unz = u(i,j,k)*ctx(i,j,k)+v(i,j,k)*cty(i,j,k)+w(i,j,k)*ctz(i,j,k)

!			    spec(i,j,k,m) = spec(i,j,k,m)+abs(unx)+abs(uny)+abs(unz)+ 1.0/dtdt(i,j,k) !old
			    spec(i,j,k,m) = spec(i,j,k,m)+abs(unx)+abs(uny)+abs(unz)+ 1.0/(dtdt(i,j,k)*timedt_turb) !new �Ӵ�ʱ�䲽��
					
					nx = kcx(i,j,k)*kcx(i,j,k)+kcy(i,j,k)*kcy(i,j,k)+kcz(i,j,k)*kcz(i,j,k)
					ny = etx(i,j,k)*etx(i,j,k)+ety(i,j,k)*ety(i,j,k)+etz(i,j,k)*etz(i,j,k)
					nz = ctx(i,j,k)*ctx(i,j,k)+cty(i,j,k)*cty(i,j,k)+ctz(i,j,k)*ctz(i,j,k)

					vis = viseq(i,j,k,m)/vol(i,j,k)/r(i,j,k)
					
					spec(i,j,k,m) = spec(i,j,k,m)+re2*vis*(nx+ny+nz) + m*0.09*qke(i,j,k,1)*vol(i,j,k)

				enddo
			enddo
		enddo
	enddo

	return
end subroutine cal_spec_turbulence_sst_tem
