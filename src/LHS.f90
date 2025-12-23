subroutine l_h_s_tgh
    use define_precision_mod
    use global_variables,only : nblocks,nlhs,nsmooth,nout
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
	integer :: nb,pnb


#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif

        call recast_GRID(nb)
        call recast_field(nb)

        call set_nonconnect_boundary_dq_0(nb) !* �ѷǶԽӱ߽��dq����-rhs������

        call GS_PR_single

        call set_nonconnect_boundary_dq_0(nb) !* ʱ���ƽ������°ѷǶԽӱ߽��dq����

		call update(nb)  !�Ѿ������dq,����q��ԭʼ����


    enddo !*nb = 1,nblocks

    return

end subroutine l_h_s_tgh


subroutine view_all
    use global_variables
#ifdef PARALLEL
    use mod_parallels
#endif
    implicit none
	integer :: nb,pnb
#ifndef PARALLEL
	integer :: myid=0
#endif

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
       call debug_view(nb)  !!(???????????????????)
    end do

    close(20+myid)
    call stop_by_error(203,"writing debug data is finished")

end subroutine view_all


subroutine debug_view(nb)
    use global_variables,only : mb_x,mb_y,mb_z,mb_r,mb_p,mb_dq,mb_dim,nblocks,mb_distance
#ifdef PARALLEL
    use mod_parallels
#endif
    implicit none
    integer :: nb,i,j,k,ni,nj,nk,m
    real    :: x,y,z,tmp1,tmp2,tmp(5)
#ifndef PARALLEL
    integer :: myid=0
#endif

    ni = mb_dim(nb,1)
    nj = mb_dim(nb,2)
    nk = mb_dim(nb,3)
    write(20+myid,*)'variables=x,y,z,tmp1,tmp2'
    write(20+myid,*)'zone i=',ni,'j=',nj,'k=',nk,' f=point'

    do k=1,nk
    do j=1,nj
    do i=1,ni
       x = mb_x(nb)%a3d(i,j,k)
       y = mb_y(nb)%a3d(i,j,k)
       z = mb_z(nb)%a3d(i,j,k)
       tmp1 = mb_r(nb)%a3d(i,j,k)
       tmp2 = mb_distance(nb)%a3d(i,j,k)
       !tmp = mb_dq(nb)%a4d(:,i,j,k)
       write(20+myid,*)x,y,z,tmp1,tmp2  !!, !!(tmp(m),m=1,5)
    end do
    end do
    end do

end subroutine debug_view

subroutine stop_by_error(ierrcode,errmsg)
#ifdef PARALLEL
    use mod_parallels
#endif
    implicit none
    integer :: ierrcode
    character(len=*) :: errmsg

#ifdef PARALLEL
    call error_mpi(ierrcode,errmsg)
#else
    write(*,*) "# "//trim(adjustl(errmsg))//":",ierrcode
    stop
#endif

end subroutine stop_by_error




subroutine spectrum_tgh
    use define_precision_mod
    use global_variables,only : nvis,nlamtur,csrv
    implicit none

    call spectinv

    if ( nvis == 1  ) then
       if ( csrv > 0.0 ) then  !*TGH. csrv  = 2.0-ճ���װ뾶Ȩֵ��Խ��ճ��Խ��*!
	      if(nlamtur < 0) then
		     call spectvisl   !*TGH. �����װ뾶
		  else
		     call spectvist   !*TGH. �����װ뾶
		  endif
	   endif

	else
	     call set_spectvis_to_0
    endif

    return
end subroutine spectrum_tgh

!_____________________________________________________________________!
subroutine spectvist !*TGH.�����װ뾶
    use global_variables,only : csrv,reynolds,r,vol,visl,vist,ni,nj,nk,srva,srvb,srvc,small
    implicit none
    integer :: i,j,k,ii,jj,kk
    real :: rm_vol,ri,rj,rk
    real :: kx,ky,kz,ex,ey,ez,cx,cy,cz,nt
    real :: coef,vis,coef1

    do k=1,nk
    do j=1,nj
    do i=1,ni
       rm_vol = r(i,j,k)* vol(i,j,k)
       call getr0kec_mml(i,j,k,kx,ky,kz,nt,1,0,0)
       call getr0kec_mml(i,j,k,ex,ey,ez,nt,0,1,0)
       call getr0kec_mml(i,j,k,cx,cy,cz,nt,0,0,1)

       ri =  kx*kx + ky*ky + kz*kz
       rj =  ex*ex + ey*ey + ez*ez
       rk =  cx*cx + cy*cy + cz*cz
       vis = visl(i,j,k) + vist(i,j,k)
       coef = 2.0*vis/( reynolds * rm_vol + small )
       coef1 = csrv*coef
       srva(i,j,k) = coef1*ri
       srvb(i,j,k) = coef1*rj
       srvc(i,j,k) = coef1*rk
    enddo
    enddo
    enddo


    do k=1,nk
    do j=1,nj
       do i=0,ni+1,ni+1 !*tgh. I����߽�
		  ii = min(ni,max(1,i))

          rm_vol = r(i,j,k)* vol(ii,j,k)
          call getr0kec_mml(ii,j,k,kx,ky,kz,nt,1,0,0)
          call getr0kec_mml(ii,j,k,ex,ey,ez,nt,0,1,0)
          call getr0kec_mml(ii,j,k,cx,cy,cz,nt,0,0,1)

          ri =  kx*kx + ky*ky + kz*kz
          rj =  ex*ex + ey*ey + ez*ez
          rk =  cx*cx + cy*cy + cz*cz
          vis = visl(ii,j,k) + vist(ii,j,k)
          coef = 2.0*vis/( reynolds * rm_vol  + small )
          coef1 = csrv*coef
          srva(i,j,k) = coef1*ri
          srvb(i,j,k) = coef1*rj
          srvc(i,j,k) = coef1*rk
       enddo
    enddo
    enddo

    do k=1,nk
    do i=1,ni
       do j=0,nj+1,nj+1 !*tgh. j����߽�
		  jj = min(nj,max(1,j))

          rm_vol = r(i,j,k)* vol(i,jj,k)
          call getr0kec_mml(i,jj,k,kx,ky,kz,nt,1,0,0)
          call getr0kec_mml(i,jj,k,ex,ey,ez,nt,0,1,0)
          call getr0kec_mml(i,jj,k,cx,cy,cz,nt,0,0,1)

          ri =  kx*kx + ky*ky + kz*kz
          rj =  ex*ex + ey*ey + ez*ez
          rk =  cx*cx + cy*cy + cz*cz
          vis = visl(i,jj,k) + vist(i,jj,k)
          coef = 2.0*vis/( reynolds * rm_vol  + small )
          coef1 = csrv*coef
          srva(i,j,k) = coef1*ri
          srvb(i,j,k) = coef1*rj
          srvc(i,j,k) = coef1*rk
       enddo
    enddo
    enddo

    do j=1,nj
    do i=1,ni
       do k=0,nk+1,nk+1 !*tgh. k����߽�
		  kk = min(nk,max(1,k))

          rm_vol = r(i,j,k)* vol(i,j,kk)
          call getr0kec_mml(i,j,kk,kx,ky,kz,nt,1,0,0)
          call getr0kec_mml(i,j,kk,ex,ey,ez,nt,0,1,0)
          call getr0kec_mml(i,j,kk,cx,cy,cz,nt,0,0,1)

          ri =  kx*kx + ky*ky + kz*kz
          rj =  ex*ex + ey*ey + ez*ez
          rk =  cx*cx + cy*cy + cz*cz
          vis = visl(i,j,kk) + vist(i,j,kk)
          coef = 2.0*vis/( reynolds * rm_vol  + small )
          coef1 = csrv*coef
          srva(i,j,k) = coef1*ri
          srvb(i,j,k) = coef1*rj
          srvc(i,j,k) = coef1*rk
       enddo
    enddo
    enddo

    return
end subroutine spectvist
!_____________________________________________________________________!
subroutine spectvisl !*TGH.�����װ뾶
    use global_variables,only : csrv,reynolds,r,vol,visl,ni,nj,nk,srva,srvb,srvc,small
    implicit none
    integer :: i,j,k,ii,jj,kk
    real :: rm_vol,ri,rj,rk
    real :: kx,ky,kz,ex,ey,ez,cx,cy,cz,nt
    real :: coef,vis,coef1

    do k=1,nk
    do j=1,nj
    do i=1,ni
       rm_vol = r(i,j,k)* vol(i,j,k)
       call getr0kec_mml(i,j,k,kx,ky,kz,nt,1,0,0)
       call getr0kec_mml(i,j,k,ex,ey,ez,nt,0,1,0)
       call getr0kec_mml(i,j,k,cx,cy,cz,nt,0,0,1)

       ri =  kx*kx + ky*ky + kz*kz
       rj =  ex*ex + ey*ey + ez*ez
       rk =  cx*cx + cy*cy + cz*cz
       vis = visl(i,j,k)
       coef = 2.0*vis/( reynolds * rm_vol + small )
       coef1 = csrv*coef
       srva(i,j,k) = coef1*ri
       srvb(i,j,k) = coef1*rj
       srvc(i,j,k) = coef1*rk
    enddo
    enddo
    enddo


    do k=1,nk
    do j=1,nj
       do i=0,ni+1,ni+1 !*tgh. I����߽�
		  ii = min(ni,max(1,i))

          rm_vol = r(i,j,k)* vol(ii,j,k)
          call getr0kec_mml(ii,j,k,kx,ky,kz,nt,1,0,0)
          call getr0kec_mml(ii,j,k,ex,ey,ez,nt,0,1,0)
          call getr0kec_mml(ii,j,k,cx,cy,cz,nt,0,0,1)

          ri =  kx*kx + ky*ky + kz*kz
          rj =  ex*ex + ey*ey + ez*ez
          rk =  cx*cx + cy*cy + cz*cz
          vis = visl(ii,j,k)
          coef = 2.0*vis/( reynolds * rm_vol  + small )
          coef1 = csrv*coef
          srva(i,j,k) = coef1*ri
          srvb(i,j,k) = coef1*rj
          srvc(i,j,k) = coef1*rk
       enddo
    enddo
    enddo

    do k=1,nk
    do i=1,ni
       do j=0,nj+1,nj+1 !*tgh. j����߽�
		  jj = min(nj,max(1,j))

          rm_vol = r(i,j,k)* vol(i,jj,k)
          call getr0kec_mml(i,jj,k,kx,ky,kz,nt,1,0,0)
          call getr0kec_mml(i,jj,k,ex,ey,ez,nt,0,1,0)
          call getr0kec_mml(i,jj,k,cx,cy,cz,nt,0,0,1)

          ri =  kx*kx + ky*ky + kz*kz
          rj =  ex*ex + ey*ey + ez*ez
          rk =  cx*cx + cy*cy + cz*cz
          vis = visl(i,jj,k)
          coef = 2.0*vis/( reynolds * rm_vol  + small )
          coef1 = csrv*coef
          srva(i,j,k) = coef1*ri
          srvb(i,j,k) = coef1*rj
          srvc(i,j,k) = coef1*rk
       enddo
    enddo
    enddo

    do j=1,nj
    do i=1,ni
       do k=0,nk+1,nk+1 !*tgh. k����߽�
		  kk = min(nk,max(1,k))

          rm_vol = r(i,j,k)* vol(i,j,kk)
          call getr0kec_mml(i,j,kk,kx,ky,kz,nt,1,0,0)
          call getr0kec_mml(i,j,kk,ex,ey,ez,nt,0,1,0)
          call getr0kec_mml(i,j,kk,cx,cy,cz,nt,0,0,1)

          ri =  kx*kx + ky*ky + kz*kz
          rj =  ex*ex + ey*ey + ez*ez
          rk =  cx*cx + cy*cy + cz*cz
          vis = visl(i,j,kk)
          coef = 2.0*vis/( reynolds * rm_vol  + small )
          coef1 = csrv*coef
          srva(i,j,k) = coef1*ri
          srvb(i,j,k) = coef1*rj
          srvc(i,j,k) = coef1*rk
       enddo
    enddo
    enddo

    return
end subroutine spectvisl
!=============================================================================!
!=============================================================================!

!_____________________________________________________________________!
subroutine explicit
    use global_variables
    implicit none
    integer :: i,j,k,m
    real :: coed(nl)


!TGH
    do k=nk,1,-1  !1,nk
       do j=1,nj    !1,nj
          do i=1,ni     !1,ni
!------------------------------------

             call get_rexp(i,j,k,coed)
             do m=1,nl
                dq(m,i,j,k) = - coed(m) * dq(m,i,j,k)
             enddo
          enddo
       enddo
    enddo

    return
end subroutine explicit
!_____________________________________________________________________!
subroutine lusgs
    implicit none
    real :: wmig,beta,sfix
    parameter(wmig=1.0,beta=1.0,sfix=0.0)
    call lusgs_l(wmig,beta,sfix)
    call lusgs_u(wmig,beta,sfix)
    return
end subroutine lusgs
!_____________________________________________________________________!
subroutine lusgs_l(wmig,beta,sfix)
    use define_precision_mod
    use global_variables
    implicit none
    integer :: i,j,k,m,cmethod
    real :: rhs0(nl)
    real :: prim_i(nl),prim_j(nl),prim_k(nl)
    real :: gykb_i(4),gykb_j(4),gykb_k(4)
    real :: dq_i(nl),dq_j(nl),dq_k(nl),de(nl),df(nl),dg(nl)
    real :: wmig,beta,sfix,ra,rb,rc,rva,rvb,rvc
    real :: coed(nl)

    cmethod = 1 + method

if(.false.)then
	i = 0
	do k=1,nk
	do j=1,nj
	do m=1,nl
		dq(m,i,j,k) = 0.0
	enddo
	enddo
	enddo
	j = 0
	do k=1,nk
	do i=1,ni
	do m=1,nl
		dq(m,i,j,k) = 0.0
	enddo
	enddo
	enddo
	k = 0
	do j=1,nj
	do i=1,ni
	do m=1,nl
		dq(m,i,j,k) = 0.0
	enddo
	enddo
	enddo
endif

    do k=1,nk            !___cic 
       do j=1,nj    
          do i=1,ni     
             call getcoed(i,j,k,beta,wmig,coed)
             do m=1,nl
                dq_i(m) = dq(m,i-1,j,k)
                dq_j(m) = dq(m,i,j-1,k)
                dq_k(m) = dq(m,i,j,k-1)


             enddo

!			   if(i==1)dq_i(1:nl)=0.0_prec    !___cic
!			   if(j==1)dq_j(1:nl)=0.0_prec
!			   if(k==1)dq_k(1:nl)=0.0_prec

             call getprim(i-1,j,k,prim_i)
             call getprim(i,j-1,k,prim_j)
             call getprim(i,j,k-1,prim_k)
!             call getprim(MAX(1,i-1),j,k,prim_i)
!             call getprim(i,MAX(1,j-1),k,prim_j)
!             call getprim(i,j,MAX(1,k-1),prim_k)

             call getrkec_mml(i-method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
             call getrkec_mml(i,j-method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k-method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k,-1,ra,rb,rc)

             call mxdq_std(prim_i,gykb_i,dq_i,de,ra, 1)
             call mxdq_std(prim_j,gykb_j,dq_j,df,rb, 1)
             call mxdq_std(prim_k,gykb_k,dq_k,dg,rc, 1)

             !---- solve rhs0
             do m=1,nl
                rhs0(m) = wmig * ( de(m) + df(m) + dg(m) )
             enddo
             
             if (nvis > 0) then
                rva = srva(i-1,j,k)
                rvb = srvb(i,j-1,k)
                rvc = srvc(i,j,k-1)
                do m=1,nl
                   rhs0(m) = rhs0(m) + 0.5*( rva*dq(m,i-1,j,k) + rvb*dq(m,i,j-1,k) + rvc*dq(m,i,j,k-1) )
                end do
             end if
             
             !---- solve dqa
             do m=1,nl
                dq(m,i,j,k) = ( -dq(m,i,j,k) + rhs0(m) ) * coed(m)
             enddo
          enddo
       enddo
    enddo

    return
end subroutine lusgs_l
!_____________________________________________________________________!
subroutine lusgs_u(wmig,beta,sfix)
    use define_precision_mod
    use global_variables
    implicit none
    integer :: i,j,k,m,cmethod
    real :: rhs0(nl)
    real :: prim_i(nl),prim_j(nl),prim_k(nl)
    real :: gykb_i(4),gykb_j(4),gykb_k(4)
    real :: dq_i(nl),dq_j(nl),dq_k(nl),de(nl),df(nl),dg(nl)
    real :: wmig,beta,sfix,ra,rb,rc,rva,rvb,rvc
    real :: coed(nl)

    cmethod = 1 + method

if(.false.)then
	i = ni+1
	do k=1,nk
	do j=1,nj
	do m=1,nl
		dq(m,i,j,k) = 0.0
	enddo
	enddo
	enddo
	j = nj+1
	do k=1,nk
	do i=1,ni
	do m=1,nl
		dq(m,i,j,k) = 0.0
	enddo
	enddo
	enddo
	k = nk+1
	do j=1,nj
	do i=1,ni
	do m=1,nl
		dq(m,i,j,k) = 0.0
	enddo
	enddo
	enddo
endif

    do k=nk,1,-1            !___cic 
       do j=nj,1,-1
          do i=ni,1,-1     
             call getcoed(i,j,k,beta,wmig,coed)
             do m=1,nl
                dq_i(m) = dq(m,i+1,j,k)
                dq_j(m) = dq(m,i,j+1,k)
                dq_k(m) = dq(m,i,j,k+1)

             enddo

!				if(i==ni)dq_i(1:nl)=0.0_prec   !___cic
!				if(j==nj)dq_j(1:nl)=0.0_prec
!				if(k==nk)dq_k(1:nl)=0.0_prec


             call getprim(i+1,j,k,prim_i)
             call getprim(i,j+1,k,prim_j)
             call getprim(i,j,k+1,prim_k)
!            call getprim(MIN(i+1,NI),j,k,prim_i)
!            call getprim(i,MIN(j+1,NJ),k,prim_j)
!            call getprim(i,j,MIN(k+1,NK),prim_k)

             call getrkec_mml(i+1,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
             call getrkec_mml(i,j+1,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k+1,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k,1,ra,rb,rc)

             call mxdq_std(prim_i,gykb_i,dq_i,de,ra,-1)
             call mxdq_std(prim_j,gykb_j,dq_j,df,rb,-1)
             call mxdq_std(prim_k,gykb_k,dq_k,dg,rc,-1)

             !---- solve rhs0
             do m=1,nl
                rhs0(m) = wmig * ( de(m) + df(m) + dg(m) )
             enddo
             
             if (nvis > 0) then
                rva = srva(i+1,j,k)
                rvb = srvb(i,j+1,k)
                rvc = srvc(i,j,k+1)
                do m=1,nl
                   rhs0(m) = rhs0(m) - 0.5*( rva*dq(m,i+1,j,k) + rvb*dq(m,i,j+1,k) + rvc*dq(m,i,j,k+1) )
                end do
             end if
             
             !---- solve dqa
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) - rhs0(m) * coed(m)
             enddo
          enddo
       enddo
    enddo

    return
end subroutine lusgs_u
!_____________________________________________________________________!
subroutine get_rexp(i,j,k,coed)
    use global_variables,ONLY : cfl,dtdt,srva,srvb,srvc,nvis,nl,nm
    implicit none
    integer :: i,j,k,m
    real :: ct
    real :: wmig,beta,ra,rb,rc,rv,rad_ns,rad_chem,ccc
    real :: coed(nl)

!    ct =  1.0/dtdt(i,j,k)
!    rad_ns = srv(i,j,k) * nvis
! !    ccc = 1.0/( ct +  rad_ns ) !* ����ճ���װ뾶��ʱ�䲽�� *
!*tgh. ��ʽ���񲻶ԣ�Ӧ�ø�Ϊ��ccc = 1.0/( ct +  rad_ns/cfl )
!    ccc = 1.0/( ct +  rad_ns/cfl ) !* ����ճ���װ뾶��ʱ�䲽�� *
!*tgh. Modified by TU Guohua, 2009.2

    do m=1,nm
!       coed(m) = ccc
       coed(m) = DTDT(I,J,K)
    enddo

end subroutine get_rexp
!-----------------------------------------------------------------
!_____________________________________________________________________!
subroutine getcoed(i,j,k,beta,wmig,coed)
    use global_variables
    implicit none
    integer :: i,j,k,m
    real :: ct
    real :: wmig,beta,ra,rb,rc,rv,rad_ns,rad_chem,ccc
    real :: coed(nl),odtst

    ra = sra(i,j,k)
    rb = srb(i,j,k)
    rc = src(i,j,k)
    rad_ns = ra + rb + rc
	do m=1,nvis
	    rad_ns = rad_ns + srva(i,j,k) + srvb(i,j,k) + srvc(i,j,k)
	enddo

    ct =  1.0/dtdt(i,j,k)
    
    
    if (ndualtst > 0) then
       odtst = 1.0/dtdts     !!reftime/dtdts
       ct = ct + 1.5*odtst*vol(i,j,k)
    end if


    ccc = 1.0/( ct + beta * wmig * rad_ns )
    do m=1,nm
       coed(m) = ccc
    enddo

    if( nl > 5 ) then
        rad_chem = 0.0
        do m=6,nl
           rad_chem = amax1( srs(m-5,i,j,k), rad_chem )
        enddo
        ra = ccc
       ccc = 1.0/( ct +  (rad_ns + rad_chem) )
       do m=6,nl
          coed(m) = ccc
       enddo
	   srs(1,i,j,k) = ccc/ra
    endif

end subroutine getcoed
!_____________________________________________________________________!
subroutine getrarbrc(i,j,k,dijk,ra,rb,rc)
    use global_variables
    implicit none
    integer :: i,j,k,dijk,i_mml,j_mml,k_mml
    real :: ra,rb,rc
    i_mml = i + dijk
    j_mml = j + dijk
    k_mml = k + dijk
    if ( method == 0 ) then
       i_mml = min(i_mml,ni-1)
       j_mml = min(j_mml,nj-1)
       k_mml = min(k_mml,nk-1)
       i_mml = max(i_mml,   1)
       j_mml = max(j_mml,   1)
       k_mml = max(k_mml,   1)
    endif

    if ( method == 1 ) then !*tgh. Added by TU Guohua
       i_mml = min(i_mml,ni)
       j_mml = min(j_mml,nj)
       k_mml = min(k_mml,nk)
       i_mml = max(i_mml, 1)
       j_mml = max(j_mml, 1)
       k_mml = max(k_mml, 1)
    endif !*tgh. end. Added by TU Guohua


    ra = sra(i_mml,j    ,k    )
    rb = srb(i    ,j_mml,k    )
    rc = src(i    ,j    ,k_mml)
    return
end subroutine getrarbrc
!_____________________________________________________________________!
subroutine mxdq(prim,gykb,dq,f,efix,npn)
    use global_const,only:nl,ns,ms1,beta1,tref,nchem,sml_sss
    implicit none
    integer :: npn,m
    real :: hint,gama,efix,eps
    real :: f(nl),prim(nl),gykb(4),dq(nl),hs(ns),as(ns)
    real :: nx,ny,nz,nt,ct,cgm,cgm1
    real :: l1,l4,l5,x1,x2
    real :: dh,dc,c2dc,ae,af
    real :: rm,um,vm,wm,pm,cm,c2,v2,tm,hm

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    call gettmhmgm(prim,tm,hint,gama)
    c2 = gama*pm/rm
    cm = sqrt(c2)
    v2 = um*um + vm*vm + wm*wm

    hm = hint + 0.5*v2

    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)

    ct = nx*um + ny*vm + nz*wm
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

    l1 = ct
    l4 = ct + cm * cgm
    l5 = ct - cm * cgm
    eps = efix*efix*cgm*cgm

    l1 = 0.5 * ( l1 + npn * sqrt(l1*l1 + eps) )
    l4 = 0.5 * ( l4 + npn * sqrt(l4*l4 + eps) )
    l5 = 0.5 * ( l5 + npn * sqrt(l5*l5 + eps) )

    x1 = ( 2.0*l1 - l4 - l5 )/( 2.0 * c2 )
    x2 = ( l4 - l5 )/( 2.0 * cm )

    cgm1 = 1.0/cgm
    nx = nx * cgm1
    ny = ny * cgm1
    nz = nz * cgm1
    ct = ( ct - nt ) * cgm1
    ae = gama - 1.0

    af = 0.5*ae*v2

    dc = ct * dq(1) - nx * dq(2) - ny * dq(3) - nz * dq(4)
    dh = af * dq(1) - ae * ( um * dq(2) + vm * dq(3) + wm * dq(4) - dq(5) )

    do m=6,nl
       dh = dh + as(m-5) * dq(m)
    enddo

    c2dc = c2 * dc

    f(1) = l1 * dq(1)             -    dh   * x1           -    dc   * x2
    f(2) = l1 * dq(2) + ( nx*c2dc - um*dh ) * x1 + ( nx*dh - um*dc ) * x2
    f(3) = l1 * dq(3) + ( ny*c2dc - vm*dh ) * x1 + ( ny*dh - vm*dc ) * x2
    f(4) = l1 * dq(4) + ( nz*c2dc - wm*dh ) * x1 + ( nz*dh - wm*dc ) * x2
    f(5) = l1 * dq(5) + ( ct*c2dc - hm*dh ) * x1 + ( ct*dh - hm*dc ) * x2

    do m=6,nl
       f(m) = l1 * dq(m) - prim(m) * ( dh * x1 + dc * x2 )
    enddo

    return

end subroutine mxdq
!_____________________________________________________________________!
subroutine mxdq_std(prim,gykb,dq,f,rad,npn)
    use global_const,only:nl,ns,ms1,beta1,tref,nchem,sml_sss
    implicit none
    integer :: npn,m
    real :: hint,gama,rad
    real :: f(nl),prim(nl),gykb(4),dq(nl),hs(ns),as(ns)
    real :: nx,ny,nz,nt,ct,cgm,cgm1
    real :: l1,l4,l5,x1,x2
    real :: dh,dc,c2dc,ae,af
    real :: rm,um,vm,wm,pm,cm,c2,v2,tm,hm

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    call gettmhmgm(prim,tm,hint,gama)
    c2 = gama*pm/rm
    cm = sqrt(c2)
    v2 = um*um + vm*vm + wm*wm

    hm = hint + 0.5*v2

    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)

    ct = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

    l1 = ct
    l4 = ct + cm * cgm
    l5 = ct - cm * cgm

    l1 = 0.5 * ( l1 + npn*rad )
    l4 = 0.5 * ( l4 + npn*rad )
    l5 = 0.5 * ( l5 + npn*rad )

    x1 = ( 2.0*l1 - l4 - l5 )/( 2.0 * c2 )
    x2 = ( l4 - l5 )/( 2.0 * cm )

    cgm1 = 1.0/cgm
    nx = nx * cgm1
    ny = ny * cgm1
    nz = nz * cgm1
    ct = ( ct - nt ) * cgm1
    ae = gama - 1.0

    af = 0.5*ae*v2

    dc = ct * dq(1) - nx * dq(2) - ny * dq(3) - nz * dq(4)
    dh = af * dq(1) - ae * ( um * dq(2) + vm * dq(3) + wm * dq(4) - dq(5) )

    do m=6,nl
       dh = dh + as(m-5) * dq(m)
    enddo

    c2dc = c2 * dc

    f(1) = l1 * dq(1)             -    dh   * x1           -    dc   * x2
    f(2) = l1 * dq(2) + ( nx*c2dc - um*dh ) * x1 + ( nx*dh - um*dc ) * x2
    f(3) = l1 * dq(3) + ( ny*c2dc - vm*dh ) * x1 + ( ny*dh - vm*dc ) * x2
    f(4) = l1 * dq(4) + ( nz*c2dc - wm*dh ) * x1 + ( nz*dh - wm*dc ) * x2
    f(5) = l1 * dq(5) + ( ct*c2dc - hm*dh ) * x1 + ( ct*dh - hm*dc ) * x2

    do m=6,nl
       f(m) = l1 * dq(m) - prim(m) * ( dh * x1 + dc * x2 )
    enddo

    return
end subroutine mxdq_std
!_____________________________________________________________________!
subroutine res_smooth
    use global_variables
    implicit none
    integer :: nb

    do nb=1,nblocks
!       call recast_method(nb)
       call recast_grid(nb)
       call recast_field(nb)

       call res_smooth_i
       call res_smooth_j
       call res_smooth_k
    enddo

    return
end subroutine res_smooth
!_____________________________________________________________________!
subroutine res_smooth_i
    use global_variables
    implicit none
    integer :: i,j,k,m,cmethod
    real :: ex
    real,pointer,dimension(:) :: s_c,s_d,s_f,s_b,s_x
    integer :: nn
    nn = ni
	cmethod = method + 1
    allocate( s_c(nn), s_d(nn), s_f(nn), s_b(nn), s_x(nn) )
    ex = 2.0
    do k=1,nk
       do j=1,nj
          do m=1,nm
             do i=1,ni
                s_b(i) = dq(m,i,j,k)
                s_c(i) = -ex
                s_d(i) = 1.0 + 2.0 * ex
                s_f(i) = -ex
             enddo
             s_c(1 ) = 0.0
             s_f(nn) = 0.0
             call trisys(s_c,s_d,s_f,s_b,s_x,nn)
             do i=cmethod,ni-1
                dq(m,i,j,k) = s_x(i+1)
             enddo
          enddo
       enddo
    enddo

    deallocate( s_c )
    deallocate( s_d )
    deallocate( s_f )
    deallocate( s_b )
    deallocate( s_x )

    return
end subroutine res_smooth_i
!_____________________________________________________________________!
subroutine res_smooth_j
    use global_variables
    implicit none
    integer :: i,j,k,m,cmethod
    real :: ey
    real,pointer,dimension(:) :: s_c,s_d,s_f,s_b,s_x
    integer :: nn
    nn = nj
	cmethod = method + 1
    allocate( s_c(nn), s_d(nn), s_f(nn), s_b(nn), s_x(nn) )
    ey = 2.0
    do k=1,nk
       do i=1,ni
          do m=1,nm
             do j=1,nj
                s_b(j) = dq(m,i,j,k)
                s_c(j) = -ey
                s_d(j) = 1.0 + 2.0 * ey
                s_f(j) = -ey
             enddo
             s_c(1 ) = 0.0
             s_f(nn) = 0.0
             call trisys(s_c,s_d,s_f,s_b,s_x,nn)
             do j=cmethod,nj-1
                dq(m,i,j,k) = s_x(j+1)
             enddo
          enddo
       enddo
    enddo
    deallocate( s_c )
    deallocate( s_d )
    deallocate( s_f )
    deallocate( s_b )
    deallocate( s_x )

    return
end subroutine res_smooth_j
!_____________________________________________________________________!
subroutine res_smooth_k
    use global_variables
    implicit none
    integer :: i,j,k,m,cmethod
    real :: ez
    real,pointer,dimension(:) :: s_c,s_d,s_f,s_b,s_x
    integer :: nn
    nn = nk
	cmethod = method + 1
    allocate( s_c(nn), s_d(nn), s_f(nn), s_b(nn), s_x(nn) )
    ez = 2.0
    do j=1,nj
       do i=1,ni
          do m=1,nm
             do k=1,nk
                s_b(k) = dq(m,i,j,k)
                s_c(k) = -ez
                s_d(k) = 1.0 + 2.0 * ez
                s_f(k) = -ez
             enddo
             s_c(1 ) = 0.0
             s_f(nn) = 0.0
             call trisys(s_c,s_d,s_f,s_b,s_x,nn)
             do k=cmethod,nk-1
                dq(m,i,j,k) = s_x(k+1)
             enddo
          enddo
       enddo
    enddo
    deallocate( s_c )
    deallocate( s_d )
    deallocate( s_f )
    deallocate( s_b )
    deallocate( s_x )

    return
end subroutine res_smooth_k
!_____________________________________________________________________!

subroutine rk_3s_pre
!* �ó�����ʱʹ�ã��ȵ�����Ϻ󣬿ɼ���preprecoss��
!* By TU Guohua
	use global_variables,only : mb_dim,nl,nmax,nblocks
	use rk_3s_global
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
	implicit none
	integer :: m,pnb

!	allocate( q00 (nl,nmax,nmax,nmax) )
	allocate( mb_q00(nblocks) )
#ifdef PARALLEL
    do pnb=1,pnblocks
       m = pnbindexs(pnb)
#else
    do m=1,nblocks
#endif
	    allocate( mb_q00(m)%a4d(nl,-1:mb_dim(m,1)+1,-1:mb_dim(m,2)+1,-1:mb_dim(m,3)+1) )
	enddo

end subroutine rk_3s_pre

!_____________________________________________________________________!
subroutine rk_3s_post
!* �ó�����ʱʹ�ã��ȵ����������������󣬿�ɾ��
	use rk_3s_global
	implicit none
	integer :: m

	 deallocate( mb_q00 )

end subroutine rk_3s_post
!_____________________________________________________________________!
!_____________________________________________________________________!


subroutine rk_3s_tgh
!* �ʺ������Խ�
!* ��ƣ�Ϳ����                             ----------------------------------!
!* ���ԣ�Ϳ����  2009.3                     ----------------------------------!
!* By Tu Guohua  20090402
	use global_variables,only : q,mb_q,ni,nj,nk,nl,nblocks
	use rk_3s_global
	use define_precision_mod
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
	implicit none
	integer :: nb,pnb,i,j,k,m

!	call rk_3s_pre
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
        call recast_field(nb)
		q00 => mb_q00(nb)%a4d

        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
        do m=1,nl
             q00(m,i,j,k) = q(m,i,j,k)
        enddo
        enddo
        enddo
        enddo
	enddo


!   ��һ��  q(1) = q(0) - rhs(q(0))

	call rk_3s1

	call rk_3s2

	call rk_3s3

!	call rk_3s_post

end subroutine rk_3s_tgh

!===========================================================================!
!===========================================================================!

subroutine rk_3s1
	use global_variables,only : q,dq,ni,nj,nk,nl,nblocks
	use rk_3s_global,only : q00
	use define_precision_mod
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs,communicate_pv_npp
#endif
	implicit none
	integer :: i,j,k,m,nb,pnb
	integer :: bg,ied,jed,ked
	real(prec) :: coed(nl)

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif

        call recast_field(nb)

        call set_nonconnect_boundary_dq_0(nb) !* �ѷǶԽӱ߽��dq����-rhs������
	    bg =1
	    ied=ni
	    jed=nj
	    ked=nk
        do k=bg,ked
        do j=bg,jed
        do i=bg,ied
            call get_rexp(i,j,k,coed)
		    do m=1,nl
                dq(m,i,j,k) = - coed(m) * dq(m,i,j,k)
		    enddo
       enddo
       enddo
       enddo

		call update(nb)  !�Ѿ������dq,����q��ԭʼ����

    enddo
    
#ifdef PARALLEL
    call communicate_pv_npp
#else
    call boundary_match_pv
#endif

end subroutine rk_3s1

!===========================================================================!
!===========================================================================!

subroutine rk_3s2
	use define_precision_mod
	use global_variables,only : q,dq,ni,nj,nk,nl,nblocks,nvis,cic1
	use rk_3s_global
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs,exchange_bc,exchange_bc_dq, &
                             communicate_dq_npp,communicate_pv_npp
#endif
	implicit none
	integer :: i,j,k,m,nb,pnb
	integer :: bg,ied,jed,ked
	real(prec) :: coed(nl),rk_mml,rk_q0,rk_q
	
#ifdef PARALLEL
    call exchange_bc
#endif

	call boundary_all

	call initial_q_dq_c_t  !�����غ�����������������ٺ��¶�

	if(nvis == 1)then
	    call r_h_s_vis !* ճ������ɢ
	endif

	if( nvis ==1 )then !*tgh. ��ճʱճ������ɢ���ƽ��
#ifdef PARALLEL
        call communicate_dq_npp
#else
        call boundary_match_dq  !*tgh. �����е�ͬʱƽ��
#endif
	endif

	call r_h_s_invis !* ��ճ����ɢ

	if(cic1 == 1)then !*tgh. �����Խ�ʱ����RHS��������Խ�
#ifdef PARALLEL
        call exchange_bc_dq
#endif
		call bc_connect_cic
	endif
	 
#ifdef PARALLEL
    call communicate_dq_npp
#else
    call boundary_match_dq  !*TGH. ���жԽӱ߽�ռ���ɢȫ��ƽ��
#endif

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif

        call recast_field(nb)
		q00 => mb_q00(nb)%a4d

        call set_nonconnect_boundary_dq_0(nb) !* �ѷǶԽӱ߽��dq����-rhs������
	    bg =1
	    ied=ni
	    jed=nj
	    ked=nk
        do k=bg,ked
        do j=bg,jed
        do i=bg,ied
            call get_rexp(i,j,k,coed)
            do m=1,nl
			    rk_mml = q(m,i,j,k)  - coed(m) * dq(m,i,j,k) ! �õ�q(2)_s
				rk_q0  = q00(m,i,j,k)
                rk_q   = ( 3._prec*rk_q0 + rk_mml )/4._prec              !  �õ�q(2)
                dq(m,i,j,k) = rk_q - rk_q0                       !  �õ�dq(2)
				 q(m,i,j,k) = rk_q0
			enddo
       enddo
       enddo
       enddo

	   call update(nb)  !�Ѿ������dq,����q��ԭʼ����

    enddo

#ifdef PARALLEL
    call communicate_pv_npp
#else
    call boundary_match_pv
#endif


end subroutine rk_3s2

!===========================================================================!
!===========================================================================!

subroutine rk_3s3
	use define_precision_mod
	use global_variables,only : q,dq,ni,nj,nk,nl,nblocks,nvis,cic1
	use rk_3s_global
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs,exchange_bc,exchange_bc_dq, &
                             communicate_dq_npp,communicate_pv_npp
#endif
	implicit none
	integer :: i,j,k,m,nb,pnb
	integer :: bg,ied,jed,ked
	real(prec) :: coed(nl),rk_mml,rk_q0,rk_q
	
#ifdef PARALLEL
    call exchange_bc
#endif

	call boundary_all

	call initial_q_dq_c_t  !�����غ�����������������ٺ��¶�

	if(nvis == 1)then
	    call r_h_s_vis !* ճ������ɢ
	endif

	if( nvis ==1 )then !*tgh. ��ճʱճ������ɢ���ƽ��
#ifdef PARALLEL
        call communicate_dq_npp
#else
        call boundary_match_dq  !*tgh. �����е�ͬʱƽ��
#endif
	endif

	call r_h_s_invis !* ��ճ����ɢ

	if(cic1 == 1)then !*tgh. �����Խ�ʱ����RHS��������Խ�
#ifdef PARALLEL
        call exchange_bc_dq
#endif
		call bc_connect_cic
	endif

#ifdef PARALLEL
    call communicate_dq_npp
#else
    call boundary_match_dq  !*TGH. ���жԽӱ߽�ռ���ɢȫ��ƽ��
#endif

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif

        call recast_field(nb)
		q00 => mb_q00(nb)%a4d

        call set_nonconnect_boundary_dq_0(nb) !* �ѷǶԽӱ߽��dq����-rhs������
	    bg =1
	    ied=ni
	    jed=nj
	    ked=nk
        do k=bg,ked
        do j=bg,jed
        do i=bg,ied
            call get_rexp(i,j,k,coed)
            do m=1,nl
			    rk_mml = q(m,i,j,k)  -  coed(m) * dq(m,i,j,k)        ! �õ�q(3)_s
				rk_q0  = q00(m,i,j,k)
                !!rk_q   = ( rk_q0 + 2.0_prec*q(m,i,j,k) )/3.0_prec  ! �õ� q(3)
                rk_q   = ( rk_q0 + 2.0_prec*rk_mml )/3.0_prec        ! �õ� q(3)
                dq(m,i,j,k) = rk_q - rk_q0                           ! �õ�dq(3)
				q(m,i,j,k) = rk_q0
			enddo
       enddo
       enddo
       enddo

	   call update(nb)  !�Ѿ������dq,����q��ԭʼ����

    enddo

#ifdef PARALLEL
    call communicate_pv_npp
#else
    call boundary_match_pv
#endif

end subroutine rk_3s3

!___________________________________________________________________!
!___________________________________________________________________!
subroutine localdt
    use global_variables
    implicit none
    real :: ra,rb,rc,rv,rabc,ttmin,ttmax
    integer :: i,j,k,cmethod

    cmethod = 1-method
    timedt = large
    
    ttmax = small

    do k=1,nk-cmethod
       do j=1,nj-cmethod
          do i=1,ni-cmethod
             ra = sra(i,j,k)
             rb = srb(i,j,k)
             rc = src(i,j,k)
			 rv = srva(i,j,k) + srvb(i,j,k) + srvc(i,j,k)
             rabc = ra + rb + rc + rv
             dtdt(i,j,k) = cfl/rabc*vol(i,j,k)
             timedt = min(timedt,dtdt(i,j,k))
             
             ttmax = max(ttmax,dtdt(i,j,k))
          enddo
       enddo
    enddo

    ttmin = min(timedt_rate*timedt, ttmax)

    timedt = small
    do k=1,nk-cmethod
       do j=1,nj-cmethod
          do i=1,ni-cmethod
             ra = sra(i,j,k)
             rb = srb(i,j,k)
             rc = src(i,j,k)
			 rv = srva(i,j,k) + srvb(i,j,k) + srvc(i,j,k)
             rabc = ra + rb + rc + rv
             if (timedt_rate > small) then
                dtdt(i,j,k) = min(cfl/rabc,ttmin/vol(i,j,k))
             else
                dtdt(i,j,k) = cfl/rabc
             end if
             timedt = max(timedt,dtdt(i,j,k)*vol(i,j,k))
          enddo
       enddo
    enddo

    return
end subroutine localdt
!___________________________________________________________________!
subroutine localdt0
    use global_variables,only: sra,srb,src,srva,srvb,srvc,method, &
                               timedt,dtdt,ni,nj,nk,cfl,vol,small,timedt_rate
    implicit none
    integer :: i,j,k,cmethod
    real :: ra,rb,rc,rv,rabc

    timedt = small
    cmethod = 1-method
    do k=1,nk-cmethod
       do j=1,nj-cmethod
          do i=1,ni-cmethod
             ra = sra(i,j,k)
             rb = srb(i,j,k)
             rc = src(i,j,k)
			 rv = srva(i,j,k) + srvb(i,j,k) + srvc(i,j,k)
			                     !*TGH. �ر�ע�⣬���ڴ˴�������ճ���װ뾶
			                     !*TGH. ����subroutine get_rexp�в�Ӧ���ڿ���ճ���װ뾶
			                     !*TGH. ��subroutine getcoed�е���ʽ��Ӧ�û���Ҫ����ճ���װ뾶
             rabc = ra + rb + rc + rv
             !!dtdt(i,j,k) = min ( cfl/rabc , 1000.0/vol(i,j,k) )
             !!dtdt(i,j,k) = cfl/rabc
             
             if (timedt_rate > small) then
                dtdt(i,j,k) = min(cfl/rabc,timedt_rate/vol(i,j,k))
             else
                dtdt(i,j,k) = cfl/rabc
             end if
             
             timedt = max(timedt,dtdt(i,j,k)*vol(i,j,k))
          enddo
       enddo
    enddo
    return
end subroutine localdt0
!___________________________________________________________________!
subroutine globaldt
    use global_variables,only: sra,srb,src,srva,srvb,srvc,method,timedt,dtdt,ni,nj,nk,cfl,vol,small
    implicit none
    integer :: i,j,k,cmethod
    real :: ra,rb,rc,rv,rabc,ttmin

    timedt = small
	ttmin  = 100.0
    cmethod = 1-method
    do k=1,nk-cmethod
       do j=1,nj-cmethod
          do i=1,ni-cmethod
             ra = sra(i,j,k)
             rb = srb(i,j,k)
             rc = src(i,j,k)
			 rv = srva(i,j,k) + srvb(i,j,k) + srvc(i,j,k)
             rabc = ra + rb + rc + rv
             ttmin = amin1 ( cfl/rabc , ttmin )
          enddo
       enddo
    enddo

    do k=1,nk-cmethod
       do j=1,nj-cmethod
          do i=1,ni-cmethod
             dtdt(i,j,k) = ttmin / vol(i,j,k)
          enddo
       enddo
    enddo
    return
end subroutine globaldt
!_____________________________________________________________________!

!=============================================================================!
!=============================================================================!
!_____________________________________________________________________!

subroutine GS_pre
!* �ó�����ʱʹ�ã��ȵ�����Ϻ󣬿ɼ���preprecoss��
	use global_variables,only : mb_dim,nl,nmax,nblocks
	use store_rhs_for_GS
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
	implicit none
	integer :: m,pnb

	allocate( mb_rhs0(nblocks) )
#ifdef PARALLEL
    do pnb=1,pnblocks
       m = pnbindexs(pnb)
#else
    do m = 1,nblocks
#endif
	    allocate( mb_rhs0(m)%a4d(nl,-1:mb_dim(m,1)+1,-1:mb_dim(m,2)+1,-1:mb_dim(m,3)+1) )
	enddo

end subroutine GS_pre
!_____________________________________________________________________!

!_____________________________________________________________________!

subroutine GS_post
!* �ó�����ʱʹ�ã��ȵ�����Ϻ󣬿ɼ���preprecoss��
	use global_variables,only : nblocks
	use store_rhs_for_GS
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
	implicit none
	integer :: m,pnb

#ifdef PARALLEL
    do pnb=1,pnblocks
       m = pnbindexs(pnb)
#else
    do m = 1,nblocks
#endif
	    deallocate( mb_rhs0(m)%a4d )
	enddo
	deallocate( mb_rhs0 )

end subroutine GS_post
!_____________________________________________________________________!
!=============================================================================!
!=============================================================================!
subroutine GS_PR_global
!--------------------------------------------------------------!
!* ����LU_SGS����ʾ��ʽԤ��һ��dq �������dq�����£�  ---------!
!* Ȼ����ø�˹-���¶����ɳڽ��м���     ----------------------!
!* ÿ�鵥�����㣬�����ֱ����ʱ���ƽ�ʱ��������Ϣ--------------!
!* ��ƣ�ëö����Ϳ����      ----------------------------------!
!* ���ԣ�Ϳ����  2009.3      ----------------------------------!
!--------------------------------------------------------------!
	use define_precision_mod
	use global_variables,only : nl,ni,nj,nk,nout,nblocks
    use store_rhs_for_GS
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
    real :: wmig,beta
    integer :: nb,pnb

	wmig = 1._prec
	beta = 1._prec

	call GS_pre !* �ó�����ʱʹ�ã��ȵ�����Ϻ󣬿ɼ���preprecoss��

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
 	    call recast_grid(nb)
        call recast_field(nb)
		rhs_nb => mb_rhs0(nb)%a4d
		call store_rhs_nb

        call lusgs_l(wmig,beta,0.0)
        call lusgs_u(wmig,beta,0.0)
	enddo

	call bc_match_connect_dq !* ƽ�����жԽ��棬��ƥ������ϵ�Dq

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
 	    call recast_grid(nb)
	    call recast_field(nb)
		rhs_nb => mb_rhs0(nb)%a4d
		call gs_pr_cic_l(wmig,beta)

		call bc_connect_dq_for_other(nb) !* �����뱾�����ڵ�������������ϵ�Dq

	enddo

	call bc_match_connect_dq !* ƽ�����жԽ��棬��ƥ������ϵ�Dq

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
!!	do nb = nblocks,1,-1
 	    call recast_grid(nb)
	    call recast_field(nb)
		rhs_nb => mb_rhs0(nb)%a4d
		call gs_pr_cic_u(wmig,beta)

		call bc_connect_dq_for_other(nb) !* �����뱾�����ڵ�������������ϵ�Dq
	enddo

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
!!	do nb = nblocks,1,-1
 	    call recast_grid(nb)
	    call recast_field(nb)
        call set_nonconnect_boundary_dq_0(nb) !* ʱ���ƽ������°ѷǶԽӱ߽��dq����
		call update(nb)  !�Ѿ������dq,����q��ԭʼ����

	enddo

	call GS_post !* �ó�����ʱʹ�ã��ȵ�����Ϻ󣬿��Բ�yong

    return
end subroutine GS_PR_global
!=============================================================================!
!=============================================================================!

!=============================================================================!
subroutine GS_PR_Single
!--------------------------------------------------------------!
!* ����LU_SGS����ʾ��ʽԤ��һ��dq �������dq�����£�  ---------!
!* Ȼ����ø�˹-���¶����ɳڽ��м���     ----------------------!
!* ÿ�鵥�����㣬�����ֱ����ʱ���ƽ�ʱ��������Ϣ--------------!
!* ��ƣ�ëö����Ϳ����      ----------------------------------!
!* ���ԣ�Ϳ����  2009.3      ----------------------------------!
!--------------------------------------------------------------!
	use define_precision_mod
	use global_variables,only : nl,ni,nj,nk,nout
    use store_rhs_for_GS
    real :: wmig,beta
	wmig = 1._prec
	beta = 1._prec

	allocate (rhs_nb(nl,ni,nj,nk) )

	call store_rhs_nb  !���浱ǰ���-rhs����û�н���ʱ���ƽ�ǰ��Dq��

    call lusgs_l(wmig,beta,0.0)
    call lusgs_u(wmig,beta,0.0)

    call gs_pr_l(wmig,beta)
	call gs_pr_u(wmig,beta)

	deallocate(rhs_nb)

    return
end subroutine GS_PR_Single

!=============================================================================!
!=============================================================================!

subroutine store_rhs_nb
!--------------------------------------------------------------!
!* ���浱ǰ���-rhs����û�н���ʱ���ƽ�ǰ��Dq��
!* ʹ��ǰ�����Ѿ�recast����
!--------------------------------------------------------------!
	use global_variables,only : nl,ni,nj,nk,dq
    use store_rhs_for_GS
    implicit none
	integer :: i,j,k,m

	do k=1,nk
	do j=1,nj
	do i=1,ni
	do m=1,nl
		rhs_nb(m,i,j,k) = dq(m,i,j,k)
	enddo
	enddo
	enddo
	enddo

	return
end subroutine store_rhs_nb

!=============================================================================!
!=============================================================================!

subroutine bc_match_connect_dq
!* ���öԽ���ͶԽ�����ϵ�dq
	use define_precision_mod
    use global_variables,only : nblocks
    implicit none
    integer :: nb

    do nb = 1,nblocks
		call bc_match_connect_dq_nb(nb)
	enddo

end subroutine bc_match_connect_dq

!=============================================================================!
!=============================================================================!

subroutine bc_match_connect_dq_nb(nb)
	use define_precision_mod
    use global_variables,only : mb_dq,mb_bc,nblocks,nl
    implicit none
    integer :: nb,nr,bctype,m,n
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,i,j,k,nt,ntarg
    integer :: it,jt,kt,it0,jt0,kt0
    integer :: s_st(3),s_ed(3),nrmax

    nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
    !���崦���߽�����
    do nr = 1,nrmax
		do m=1,3
			s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
			s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
		enddo
!		s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
		s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!		s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
		nbs   = mb_bc(nb)%bc(nr)%nbs             !���

		nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
!		t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
		t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!		t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)

        bctype = mb_bc(nb)%bc(nr)%bctype    !* �Խ����� *!
        if( bctype < 0 ) then               !*�Խӱ߽磬ֱ�Ӹ�����ϵ�ֵ,�߽�����ȡƽ��
			do k = s_st(3),s_ed(3)
			do j = s_st(2),s_ed(2)
			do i = s_st(1),s_ed(1)
				is = i + mb_bc(nb)%bc(nr)%s_lr3d(1)
				js = j + mb_bc(nb)%bc(nr)%s_lr3d(2)
				ks = k + mb_bc(nb)%bc(nr)%s_lr3d(3)
				it0 = mb_bc(nb)%bc(nr)%image(i,j,k )
				jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,k )
				kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,k )

				it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)
				jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)
				kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)

				do m=1,nl
					mb_dq(nbs)%a4d(m,is,js,ks) = mb_dq(nbt)%a4d(m,it,jt,kt)
					mb_dq(nbs)%a4d(m,i,j,k) = 0.5_prec*( mb_dq(nbs)%a4d(m,i,j,k) + mb_dq(nbt)%a4d(m,it0,jt0,kt0) )
					mb_dq(nbt)%a4d(m,it0,jt0,kt0) = mb_dq(nbs)%a4d(m,i,j,k)
				enddo

			enddo
			enddo
			enddo

		else !*�ǶԽӱ߽磬�߽�����dq=0��������ϵ�dq=0
			do k = s_st(3),s_ed(3)
			do j = s_st(2),s_ed(2)
			do i = s_st(1),s_ed(1)
				is = i + mb_bc(nb)%bc(nr)%s_lr3d(1)
				js = j + mb_bc(nb)%bc(nr)%s_lr3d(2)
				ks = k + mb_bc(nb)%bc(nr)%s_lr3d(3)
				do m=1,nl
					mb_dq(nbs)%a4d(m,i ,j ,k ) = 0._prec
					mb_dq(nbs)%a4d(m,is,js,ks) = 0._prec
				enddo

			enddo
			enddo
			enddo

		endif !* bctype < 0

	enddo !* nr = 1,nrmax


end subroutine bc_match_connect_dq_nb

!=============================================================================!
!=============================================================================!

subroutine bc_connect_dq_for_other(nb)
!* ������nb�����ڵ�������������ϵ�Dq
	use define_precision_mod
    use global_variables,only : mb_dq,mb_bc,nblocks,nl
    implicit none
    integer :: nb,nr,bctype,m,n
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,i,j,k,nt,ntarg
    integer :: it,jt,kt,it0,jt0,kt0
    integer :: s_st(3),s_ed(3),nrmax

    nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
    !���崦���߽�����
    do nr = 1,nrmax
		do m=1,3
			s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
			s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
		enddo
!		s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
		s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!		s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
		nbs   = mb_bc(nb)%bc(nr)%nbs             !���

		nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
!		t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
		t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!		t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)

        bctype = mb_bc(nb)%bc(nr)%bctype    !* �Խ����� *!
        if( bctype < 0 ) then               !*�Խӱ߽磬ֱ�Ӹ�����ϵ�ֵ,�߽�����ȡƽ��
			do k = s_st(3),s_ed(3)
			do j = s_st(2),s_ed(2)
			do i = s_st(1),s_ed(1)
				is = i - mb_bc(nb)%bc(nr)%s_lr3d(1)
				js = j - mb_bc(nb)%bc(nr)%s_lr3d(2)
				ks = k - mb_bc(nb)%bc(nr)%s_lr3d(3)
				it0 = mb_bc(nb)%bc(nr)%image(i,j,k )
				jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,k )
				kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,k )

				it = it0 + mb_bc(nb)%bc(nr)%t_lr3d(1)
				jt = jt0 + mb_bc(nb)%bc(nr)%t_lr3d(2)
				kt = kt0 + mb_bc(nb)%bc(nr)%t_lr3d(3)

				do m=1,nl
					mb_dq(nbt)%a4d(m,it,jt,kt) = mb_dq(nbs)%a4d(m,is,js,ks)
				enddo

			enddo
			enddo
			enddo

		endif !* bctype < 0

	enddo !* nr = 1,nrmax


end subroutine bc_connect_dq_for_other

!=============================================================================!
!=============================================================================!

subroutine gs_pr_l(wmig,beta)
!-----------------------------------------------------------------------------!
!* �򻯰�ĸ�˹-���¶����ɳڵ��������ϡ�ɨ��     -----------------------------!
!* �򻯲��֣������D�ǲ����������ֵ���ѣ��õ���     -------------------------!
!*                     D = diag(d); d = 1/J + Dt * (Ra +Rb +Rc +Rv)   --------!
!-----------------------------------------------------------------------------!
    use define_precision_mod
    use global_variables,only : nl,ni,nj,nk,dq,method
    use store_rhs_for_GS
    implicit none
    integer :: i,j,k,m
    real :: rhs0(nl)
    real :: prim_i(nl),prim_j(nl),prim_k(nl)
    real :: gykb_i(4),gykb_j(4),gykb_k(4)
    real :: dq_i(nl),dq_j(nl),dq_k(nl)
	real :: del(nl),dfl(nl),dgl(nl),der(nl),dfr(nl),dgr(nl)
    real :: wmig,beta,sfix,ra,rb,rc
    real :: coed(nl)


    do k=1,nk
       do j=1,nj
          do i=1,ni
             call getcoed(i,j,k,beta,wmig,coed)
!*TGH. ��õ�coed = Dt * (1/J + Dt*BETA*WMIG*(Ra + Rb + Rc + Rv ) !

			!-- ������߲��֣��Ѿ��Ǹ��º��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i-1,j,k)
                dq_j(m) = dq(m,i,j-1,k)
                dq_k(m) = dq(m,i,j,k-1)
             enddo

			   if(i==1)dq_i(1:nl)=0.0_prec
			   if(j==1)dq_j(1:nl)=0.0_prec
			   if(k==1)dq_k(1:nl)=0.0_prec

             call getprim(i-1,j,k,prim_i)
             call getprim(i,j-1,k,prim_j)
             call getprim(i,j,k-1,prim_k)

             call getrkec_mml(i-method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0) !*�õ�������
             call getrkec_mml(i,j-method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k-method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k,-1,ra,rb,rc) !�õ��װ뾶,����-1��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,del,ra, 1) !* (A_i~+) * (dQ_i-1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfl,rb, 1) !* (A_j~+) * (dQ_j-1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgl,rc, 1) !* (A_k~+) * (dQ_k-1)

			!-- �����ұ߲��֣�����ǰ��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i+1,j,k)
                dq_j(m) = dq(m,i,j+1,k)
                dq_k(m) = dq(m,i,j,k+1)
             enddo

				if(i==ni)dq_i(1:nl)=0.0_prec   !___cic
				if(j==nj)dq_j(1:nl)=0.0_prec
				if(k==nk)dq_k(1:nl)=0.0_prec

             call getprim(i+1,j,k,prim_i)
             call getprim(i,j+1,k,prim_j)
             call getprim(i,j,k+1,prim_k)

             call getrkec_mml(i+method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
             call getrkec_mml(i,j+method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k+method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k, 1,ra,rb,rc) !�õ��װ뾶������ 1 ��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,der,ra,-1) !* (A_i~-) * (dQ_i+1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfr,rb,-1) !* (A_j~-) * (dQ_j+1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgr,rc,-1) !* (A_k~-) * (dQ_k+1)

             !---- solve rhs0
             do m=1,nl
                rhs0(m) = wmig * ( del(m) + dfl(m) + dgl(m) - der(m) - dfr(m) - dgr(m) )
             enddo
             !---- solve dqa
             do m=1,nl
                dq(m,i,j,k) = ( rhs0(m) - rhs_nb(m,i,j,k) ) * coed(m)
             enddo
          enddo
       enddo
    enddo

    return
end subroutine gs_pr_l
!=============================================================================!
!=============================================================================!
subroutine gs_pr_cic_l(wmig,beta)
!-----------------------------------------------------------------------------!
!* �򻯰�ĸ�˹-���¶����ɳڵ��������ϡ�ɨ��     -----------------------------!
!* �ʺ϶Խӱ߽磬�ԽӴ����ڵ㴦��                -----------------------------!
!* �򻯲��֣������D�ǲ����������ֵ���ѣ��õ���     -------------------------!
!*                     D = diag(d); d = 1/J + Dt * (Ra +Rb +Rc +Rv)   --------!
!-----------------------------------------------------------------------------!
    use define_precision_mod
    use global_variables,only : nl,ni,nj,nk,dq,method
    use store_rhs_for_GS
    implicit none
    integer :: i,j,k,m
    real :: rhs0(nl)
    real :: prim_i(nl),prim_j(nl),prim_k(nl)
    real :: gykb_i(4),gykb_j(4),gykb_k(4)
    real :: dq_i(nl),dq_j(nl),dq_k(nl)
	real :: del(nl),dfl(nl),dgl(nl),der(nl),dfr(nl),dgr(nl)
    real :: wmig,beta,sfix,ra,rb,rc
    real :: coed(nl)


    do k=1,nk
       do j=1,nj
          do i=1,ni
             call getcoed(i,j,k,beta,wmig,coed)
!*TGH. ��õ�coed = Dt * (1/J + Dt*BETA*WMIG*(Ra + Rb + Rc + Rv ) !

			!-- ������߲��֣��Ѿ��Ǹ��º��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i-1,j,k)
                dq_j(m) = dq(m,i,j-1,k)
                dq_k(m) = dq(m,i,j,k-1)
             enddo

             call getprim(i-1,j,k,prim_i)
             call getprim(i,j-1,k,prim_j)
             call getprim(i,j,k-1,prim_k)

             call getrkec_mml(i-method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0) !*�õ�������
             call getrkec_mml(i,j-method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k-method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k,-1,ra,rb,rc) !�õ��װ뾶,����-1��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,del,ra, 1) !* (A_i~+) * (dQ_i-1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfl,rb, 1) !* (A_j~+) * (dQ_j-1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgl,rc, 1) !* (A_k~+) * (dQ_k-1)

			!-- �����ұ߲��֣�����ǰ��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i+1,j,k)
                dq_j(m) = dq(m,i,j+1,k)
                dq_k(m) = dq(m,i,j,k+1)
             enddo

             call getprim(i+1,j,k,prim_i)
             call getprim(i,j+1,k,prim_j)
             call getprim(i,j,k+1,prim_k)

             call getrkec_mml(i+method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
             call getrkec_mml(i,j+method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k+method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k, 1,ra,rb,rc) !�õ��װ뾶������ 1 ��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,der,ra,-1) !* (A_i~-) * (dQ_i+1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfr,rb,-1) !* (A_j~-) * (dQ_j+1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgr,rc,-1) !* (A_k~-) * (dQ_k+1)

             !---- solve rhs0
             do m=1,nl
                rhs0(m) = wmig * ( del(m) + dfl(m) + dgl(m) - der(m) - dfr(m) - dgr(m) )
             enddo
             !---- solve dqa
             do m=1,nl
                dq(m,i,j,k) = ( rhs0(m) - rhs_nb(m,i,j,k) ) * coed(m)
             enddo
          enddo
       enddo
    enddo

    return
end subroutine gs_pr_cic_l

!=============================================================================!
!=============================================================================!

subroutine gs_pr_u(wmig,beta)
!-----------------------------------------------------------------------------!
!* �򻯰�ĸ�˹-���¶����ɳڵ��������¡�ɨ��    ------------------------------!
!* �ʺ϶Խӱ߽磬�ԽӴ����ڵ㴦��                -----------------------------!
!* �򻯲��֣������D�ǲ����������ֵ���ѣ��õ���     -------------------------!
!*                     D = diag(d); d = 1/J + Dt * (Ra +Rb +Rc +Rv)   --------!
!-----------------------------------------------------------------------------!
    use define_precision_mod
    use global_variables,only : nl,ni,nj,nk,dq,method
    use store_rhs_for_GS
    implicit none
    integer :: i,j,k,m
    real :: rhs0(nl)
    real :: prim_i(nl),prim_j(nl),prim_k(nl)
    real :: gykb_i(4),gykb_j(4),gykb_k(4)
    real :: dq_i(nl),dq_j(nl),dq_k(nl)
	real :: del(nl),dfl(nl),dgl(nl),der(nl),dfr(nl),dgr(nl)
    real :: wmig,beta,sfix,ra,rb,rc
    real :: coed(nl)


    do k=nk,1,-1
    do j=nj,1,-1
    do i=ni,1,-1
             call getcoed(i,j,k,beta,wmig,coed)
!*TGH. ��õ�coed = Dt * (1/J + Dt*BETA*WMIG*(Ra + Rb + Rc + Rv ) !

			!-- ������߲��֣�����ǰ��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i-1,j,k)
                dq_j(m) = dq(m,i,j-1,k)
                dq_k(m) = dq(m,i,j,k-1)
             enddo

			   if(i==1)dq_i(1:nl)=0.0_prec
			   if(j==1)dq_j(1:nl)=0.0_prec
			   if(k==1)dq_k(1:nl)=0.0_prec

             call getprim(i-1,j,k,prim_i)
             call getprim(i,j-1,k,prim_j)
             call getprim(i,j,k-1,prim_k)

             call getrkec_mml(i-method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0) !*�õ�������
             call getrkec_mml(i,j-method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k-method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k,-1,ra,rb,rc) !�õ��װ뾶,����-1��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,del,ra, 1) !* (A_i~+) * (dQ_i-1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfl,rb, 1) !* (A_j~+) * (dQ_j-1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgl,rc, 1) !* (A_k~+) * (dQ_k-1)

			!-- �����ұ߲��֣����º��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i+1,j,k)
                dq_j(m) = dq(m,i,j+1,k)
                dq_k(m) = dq(m,i,j,k+1)
             enddo

				if(i==ni)dq_i(1:nl)=0.0_prec   !___cic
				if(j==nj)dq_j(1:nl)=0.0_prec
				if(k==nk)dq_k(1:nl)=0.0_prec

             call getprim(i+1,j,k,prim_i)
             call getprim(i,j+1,k,prim_j)
             call getprim(i,j,k+1,prim_k)

             call getrkec_mml(i+method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
             call getrkec_mml(i,j+method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k+method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k, 1,ra,rb,rc) !�õ��װ뾶������ 1 ��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,der,ra,-1) !* (A_i~-) * (dQ_i+1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfr,rb,-1) !* (A_j~-) * (dQ_j+1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgr,rc,-1) !* (A_k~-) * (dQ_k+1)

             !---- solve rhs0
             do m=1,nl
                rhs0(m) = wmig * ( del(m) + dfl(m) + dgl(m) - der(m) - dfr(m) - dgr(m) )
             enddo
             !---- solve dqa
             do m=1,nl
                dq(m,i,j,k) = ( rhs0(m) - rhs_nb(m,i,j,k) ) * coed(m)
             enddo
    enddo
    enddo
    enddo

    return
end subroutine gs_pr_u

!=============================================================================!
!=============================================================================!

subroutine gs_pr_cic_u(wmig,beta)
!-----------------------------------------------------------------------------!
!* �򻯰�ĸ�˹-���¶����ɳڵ��������¡�ɨ��    ------------------------------!
!* �򻯲��֣������D�ǲ����������ֵ���ѣ��õ���     -------------------------!
!*                     D = diag(d); d = 1/J + Dt * (Ra +Rb +Rc +Rv)   --------!
!-----------------------------------------------------------------------------!
    use define_precision_mod
    use global_variables,only : nl,ni,nj,nk,dq,method
    use store_rhs_for_GS
    implicit none
    integer :: i,j,k,m
    real :: rhs0(nl)
    real :: prim_i(nl),prim_j(nl),prim_k(nl)
    real :: gykb_i(4),gykb_j(4),gykb_k(4)
    real :: dq_i(nl),dq_j(nl),dq_k(nl)
	real :: del(nl),dfl(nl),dgl(nl),der(nl),dfr(nl),dgr(nl)
    real :: wmig,beta,sfix,ra,rb,rc
    real :: coed(nl)


    do k=nk,1,-1
    do j=nj,1,-1
    do i=ni,1,-1
             call getcoed(i,j,k,beta,wmig,coed)
!*TGH. ��õ�coed = Dt * (1/J + Dt*BETA*WMIG*(Ra + Rb + Rc + Rv ) !

			!-- ������߲��֣�����ǰ��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i-1,j,k)
                dq_j(m) = dq(m,i,j-1,k)
                dq_k(m) = dq(m,i,j,k-1)
             enddo

             call getprim(i-1,j,k,prim_i)
             call getprim(i,j-1,k,prim_j)
             call getprim(i,j,k-1,prim_k)

             call getrkec_mml(i-method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0) !*�õ�������
             call getrkec_mml(i,j-method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k-method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k,-1,ra,rb,rc) !�õ��װ뾶,����-1��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,del,ra, 1) !* (A_i~+) * (dQ_i-1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfl,rb, 1) !* (A_j~+) * (dQ_j-1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgl,rc, 1) !* (A_k~+) * (dQ_k-1)

			!-- �����ұ߲��֣����º��dq --!
             do m=1,nl
                dq_i(m) = dq(m,i+1,j,k)
                dq_j(m) = dq(m,i,j+1,k)
                dq_k(m) = dq(m,i,j,k+1)
             enddo

             call getprim(i+1,j,k,prim_i)
             call getprim(i,j+1,k,prim_j)
             call getprim(i,j,k+1,prim_k)

             call getrkec_mml(i+method,j,k,gykb_i(1),gykb_i(2),gykb_i(3),gykb_i(4),1,0,0)
             call getrkec_mml(i,j+method,k,gykb_j(1),gykb_j(2),gykb_j(3),gykb_j(4),0,1,0)
             call getrkec_mml(i,j,k+method,gykb_k(1),gykb_k(2),gykb_k(3),gykb_k(4),0,0,1)

             call getrarbrc(i,j,k, 1,ra,rb,rc) !�õ��װ뾶������ 1 ��ʾ��ƫһ����

             call mxdq_std(prim_i,gykb_i,dq_i,der,ra,-1) !* (A_i~-) * (dQ_i+1)
             call mxdq_std(prim_j,gykb_j,dq_j,dfr,rb,-1) !* (A_j~-) * (dQ_j+1)
             call mxdq_std(prim_k,gykb_k,dq_k,dgr,rc,-1) !* (A_k~-) * (dQ_k+1)

             !---- solve rhs0
             do m=1,nl
                rhs0(m) = wmig * ( del(m) + dfl(m) + dgl(m) - der(m) - dfr(m) - dgr(m) )
             enddo
             !---- solve dqa
             do m=1,nl
                dq(m,i,j,k) = ( rhs0(m) - rhs_nb(m,i,j,k) ) * coed(m)
             enddo
    enddo
    enddo
    enddo

    return
end subroutine gs_pr_cic_u

!=============================================================================!
!=============================================================================!
